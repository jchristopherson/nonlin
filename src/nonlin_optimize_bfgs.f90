! nonlin_optimize_bfgs.f90

submodule (nonlin_optimize) nonlin_optimize_bfgs
contains
    module subroutine bfgs_solve(this, fcn, x, fout, ib, err)
        ! Arguments
        class(bfgs), intent(inout) :: this
        class(fcnnvar_helper), intent(in) :: fcn
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(out), optional :: fout
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: negone = -1.0d0
        real(real64), parameter :: factor = 1.0d2
        real(real64), parameter :: small = 1.0d-10

        ! Local Variables
        logical :: xcnvrg, gcnvrg
        integer(int32) :: i, n, maxeval, neval, ngrad, flag, iter
        real(real64) :: xtol, gtol, fp, stpmax, fret, xtest, gtest, temp, ydx
        real(real64), allocatable, dimension(:) :: g, dx, u, v, y, gold, xnew, bdx
        real(real64), allocatable, dimension(:,:) :: b, r
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        type(iteration_behavior) :: lib
        class(line_search), allocatable :: ls

        ! Initialization
        n = fcn%get_variable_count()
        maxeval = this%get_max_fcn_evals()
        gtol = this%get_tolerance()
        xtol = this%get_var_tolerance()
        iter = 0
        neval = 0
        ngrad = 0
        xcnvrg = .false.
        gcnvrg = .false.
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = 0
            ib%gradient_count = ngrad
            ib%converge_on_fcn = .false.
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

        ! Input Check
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("bfgs_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (size(x) /= n) then
            write(errmsg, '(AI0AI0A)') &
                "It was expected to receive a coordinate vector of length ", &
                n, " , but a vector of length ", size(x), " was received."
            call errmgr%report_error("bfgs_solve", trim(errmsg), &
                NL_INVALID_INPUT_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(g(n), stat = flag)
        if (flag == 0) allocate(dx(n), stat = flag)
        if (flag == 0) allocate(u(n), stat = flag)
        if (flag == 0) allocate(v(n), stat = flag)
        if (flag == 0) allocate(y(n), stat = flag)
        if (flag == 0) allocate(bdx(n), stat = flag)
        if (flag == 0) allocate(gold(n), stat = flag)
        if (flag == 0) allocate(xnew(n), stat = flag)
        if (flag == 0) allocate(b(n,n), stat = flag)
        if (flag == 0) allocate(r(n,n), stat = flag)
        if (flag /= 0) then
            ! ERROR: Memory Error
            call errmgr%report_error("bfgs_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Process
        fp = fcn%fcn(x)
        call fcn%gradient(x, g, fp)
        neval = 1
        ngrad = 1

        ! Check for a "zero" gradient at the initial point
        gtest = norm2(g)
        if (gtest < gtol) then
            gcnvrg = .true.
        end if

        ! Main Loop
        flag = 0
        if (.not.gcnvrg) then
            do
                ! Update the iteration counter
                iter = iter + 1

                ! Define the initial direction, and a limit on the line search
                ! step
                if (iter == 1) then
                    dx = -g
                    stpmax = factor * max(norm2(x), real(n, real64))
                end if

                ! Perform the line search
                if (this%get_use_line_search()) then
                    call limit_search_vector(dx, stpmax)
                    call ls%search(fcn, x, g, dx, xnew, fp, fret, lib, errmgr)
                    neval = neval + lib%fcn_count
                    fp = fret
                else
                    xnew = x + dx
                    fp = fcn%fcn(xnew)
                    neval = neval + 1
                end if

                ! Update the gradient and line direction
                do i = 1, n
                    dx(i) = xnew(i) - x(i)
                    x(i) = xnew(i)
                    gold(i) = g(i)
                end do
                call fcn%gradient(x, g, fp)
                ngrad = ngrad + 1

                ! Test for convergence on the change in X
                xtest = zero
                do i = 1, n
                    temp = abs(dx(i)) / max(abs(x(i)), one)
                    xtest = max(temp, xtest)
                end do
                if (xtest < xtol) then
                    xcnvrg = .true.
                    exit
                end if

                ! Test for convergence on the gradient
                gtest = norm2(g)
                if (gtest < gtol) then
                    gcnvrg = .true.
                    exit
                end if

                ! Perform the BFGS update
                y = g - gold
                ydx = dot_product(y, dx)

                ! Establish an initial approximation to the Hessian matrix
                if (iter == 1) then
                    temp = sqrt(dot_product(y, y) / ydx)
                    call dlaset('A', n, n, zero, temp, r, n)
                end if

                ! Compute: B = R**T * R
                call tri_mtx_mult(.true., one, r, zero, b)

                ! Compute bdx = B * dX (B is symmetric)
                call dsymv('u', n, one, b, n, dx, 1, zero, bdx, 1)

                ! Perform the actual update
                if (ydx > small) then
                    ! Compute the rank 1 update and downdate
                    u = y / sqrt(ydx)
                    v = bdx / sqrt(dot_product(dx, bdx))
                    call cholesky_rank1_update(r, u)
                    call cholesky_rank1_downdate(r, v)
                end if ! Else just skip the update

                ! Compute the solution to: B * dx = -g = (R**T * R) * dx
                dx = -g
                call solve_cholesky(.true., r, dx)

                ! Print iteration status
                if (this%get_print_status()) then
                    print *, ""
                    print '(AI0)', "Iteration: ", iter
                    print '(AI0)', "Function Evaluations: ", neval
                    print '(AE8.3)', "Function Value: ", fp
                    print '(AE8.3)', "Change in Variable: ", xtest
                    print '(AE8.3)', "Gradient: ", gtest
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
            ib%jacobian_count = 0
            ib%gradient_count = ngrad
            ib%converge_on_fcn = .false.
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = gcnvrg
        end if

        ! Get the function value at the computed minimum
        if (present(fout)) fout = fp

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, '(AI0AE8.3AE8.3AE8.3)') &
                "The algorithm failed to converge." // new_line('c') // &
                "Function evaluations performed: ", neval, new_line('c') // &
                "Function Value: ", fp, new_line('c') // &
                "Change in Variable: ", xtest, new_line('c') // &
                "Gradient: ", gtest
            call errmgr%report_error("bfgs_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if
    end subroutine
end submodule
