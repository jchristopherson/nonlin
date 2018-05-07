! nonlin_solve_newton.f90

submodule (nonlin_solve) nonlin_solve_newton
contains
! ------------------------------------------------------------------------------
    module subroutine ns_solve(this, fcn, x, fvec, ib, err)
        ! Arguments
        class(newton_solver), intent(inout) :: this
        class(vecfcn_helper), intent(in) :: fcn
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: fvec
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: half = 0.5d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: mintol = 1.0d-12
        real(real64), parameter :: factor = 1.0d2

        ! Local Variables
        logical :: check, xcnvrg, fcnvrg, gcnvrg
        integer(int32) :: i, neqn, nvar, lwork, flag, neval, iter, maxeval, njac
        integer(int32), allocatable, dimension(:) :: ipvt
        real(real64), allocatable, dimension(:) :: dir, grad, xold, work
        real(real64), allocatable, dimension(:,:) :: jac
        real(real64) :: ftol, xtol, gtol, f, fold, stpmax, xnorm, fnorm, temp, test
        type(iteration_behavior) :: lib
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        class(line_search), allocatable :: ls

        ! Initialization
        xcnvrg = .false.
        fcnvrg = .false.
        gcnvrg = .false.
        neqn = fcn%get_equation_count()
        nvar = fcn%get_variable_count()
        neval = 0
        iter = 0
        njac = 0
        ftol = this%get_fcn_tolerance()
        xtol = this%get_var_tolerance()
        gtol = this%get_gradient_tolerance()
        maxeval = this%get_max_fcn_evals()
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = njac
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = gcnvrg
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        if (this%get_use_line_search()) then
            if (.not.this%is_line_search_defined()) &
                call this%set_default_line_search()
            call this%get_line_search(ls)
        end if

        ! Input Checking
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("ns_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (nvar /= neqn) then
            ! ERROR: # of equations doesn't match # of variables
            write(errmsg, '(AI0AI0A)') "The number of equations (", neqn, &
                ") does not match the number of unknowns (", nvar, ")."
            call errmgr%report_error("ns_solve", trim(errmsg), &
                NL_INVALID_INPUT_ERROR)
            return
        end if
        flag = 0
        if (size(x) /= nvar) then
            flag = 3
        else if (size(fvec) /= neqn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("ns_solve", trim(errmsg), &
                NL_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(ipvt(nvar), stat = flag)
        if (flag == 0) allocate(dir(nvar), stat = flag)
        if (flag == 0) allocate(grad(nvar), stat = flag)
        if (flag == 0) allocate(xold(nvar), stat = flag)
        if (flag == 0) allocate(jac(nvar, neqn), stat = flag)
        if (flag == 0) then
            call fcn%jacobian(x, jac, fv = fvec, olwork = lwork)
            allocate(work(lwork), stat = flag)
        end if
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("ns_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Test to see if the initial guess is a root
        call fcn%fcn(x, fvec)
        f = half * dot_product(fvec, fvec)
        neval = neval + 1
        test = zero
        do i = 1, neqn
            test = max(abs(fvec(i)), test)
        end do
        if (test < ftol) then
            fcnvrg = .true.
        end if

        ! Process
        flag = 0 ! Used to check for convergence errors
        if (.not.fcnvrg) then
            ! Compute the maximum step size for the line search process
            stpmax = factor * max(norm2(x), real(nvar, real64))

            ! Main Iteration Loop
            do
                ! Increment the iteration counter
                iter = iter + 1

                ! Compute the Jacobian
                call fcn%jacobian(x, jac, fvec, work)
                njac = njac + 1

                ! Compute the gradient
                do i = 1, nvar
                    grad(i) = dot_product(jac(:,i), fvec)
                end do

                ! Compute the LU factorization of the Jacobian
                call lu_factor(jac, ipvt, errmgr)
                if (errmgr%has_warning_occurred()) then
                    ! The Jacobian is singular - warning was issued already, so
                    ! simply exit the routine.  Do not return as a return at
                    ! this point would not allow for proper updating of the
                    ! iteration tracking parameters
                    exit
                end if

                ! Store previous iteration values
                xold = x
                fold = f

                ! Define the right-hand-side for the linear system
                dir = -fvec

                ! Solve the linear system of equations
                call solve_lu(jac, ipvt, dir)

                ! Apply the line search if needed
                if (this%get_use_line_search()) then
                    ! Define the step length for the line search
                    temp = dot_product(dir, dir)
                    if (temp > stpmax) dir = dir * (stpmax / temp)

                    ! Apply the line search
                    call limit_search_vector(dir, stpmax)
                    call ls%search(fcn, xold, grad, dir, x, fvec, &
                        fold, f, lib, errmgr)
                    neval = neval + lib%fcn_count
                else
                    ! No line search - just update the solution estimate
                    x = x + dir
                    call fcn%fcn(x, fvec)
                    f = half * dot_product(fvec, fvec)
                    neval = neval + 1
                end if

                ! Check for convergence
                call test_convergence(x, xold, fvec, grad, .true., xtol, &
                    ftol, gtol, check, xcnvrg, fcnvrg, gcnvrg, xnorm, fnorm)
                if (check) then
                    ! The solution has converged
                    exit
                else if (gcnvrg) then
                    ! The solution appears to have settled at a point where
                    ! the gradient has a zero slope
                    write(errmsg, '(AI0AE8.3AE8.3)') &
                        "It appears the solution has settled to " // &
                        "a point where the slope of the gradient " // &
                        "is effectively zero.  " // new_line('c') // &
                        "Function evaluations performed: ", neval, &
                        "." // new_line('c') // &
                        "Change in Variable: ", xnorm, &
                        new_line('c') // "Residual: ", fnorm
                    call errmgr%report_warning("ns_solve", trim(errmsg), &
                        NL_SPURIOUS_CONVERGENCE_ERROR)
                end if

                ! Print status
                if (this%get_print_status()) then
                    call print_status(iter, neval, njac, xnorm, fnorm)
                end if

                ! Ensure we haven't made too many function evaluations
                if (neval >= maxeval) then
                    flag = 1
                    exit
                end if
            end do
        end if

        ! Report out iteration statistics
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = njac
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = gcnvrg
        end if

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, '(AI0AE8.3AE8.3)') "The algorithm failed to " // &
                "converge.  Function evaluations performed: ", neval, &
                "." // new_line('c') // "Change in Variable: ", xnorm, &
                new_line('c') // "Residual: ", fnorm
            call errmgr%report_error("ns_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if
    end subroutine
end submodule
