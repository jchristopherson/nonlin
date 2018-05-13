! nonlin_solve_quasi_newton.f90

submodule (nonlin_solve) nonlin_solve_quasi_newton
contains
! ------------------------------------------------------------------------------
    module subroutine qns_solve(this, fcn, x, fvec, ib, err)
        ! Arguments
        class(quasi_newton_solver), intent(inout) :: this
        class(vecfcn_helper), intent(in) :: fcn
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: fvec
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: half = 0.5d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: factor = 1.0d2

        ! Local Variables
        logical :: restart, xcnvrg, fcnvrg, gcnvrg, check
        integer(int32) :: i, neqn, nvar, flag, lw1, lw2, lw3, neval, iter, &
            maxeval, jcount, njac
        real(real64), allocatable, dimension(:) :: work, tau, dx, df, fvold, &
            xold, s
        real(real64), allocatable, dimension(:,:) :: q, r, b
        real(real64) :: test, f, fold, temp, ftol, xtol, gtol, &
            stpmax, x2, xnorm, fnorm
        type(iteration_behavior) :: lib
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        class(line_search), allocatable :: ls

        ! Initialization
        restart = .true.
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

        ! Input Check
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("qns_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (nvar /= neqn) then
            ! ERROR: # of equations doesn't match # of variables
            write(errmsg, '(AI0AI0A)') "The number of equations (", neqn, &
                ") does not match the number of unknowns (", nvar, ")."
            call errmgr%report_error("qns_solve", trim(errmsg), &
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
            call errmgr%report_error("qns_solve", trim(errmsg), &
                NL_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(q(neqn, neqn), stat = flag)
        if (flag == 0) allocate(r(neqn, nvar), stat = flag)
        if (flag == 0) allocate(tau(min(neqn, nvar)), stat = flag)
        if (flag == 0) allocate(b(neqn, nvar), stat = flag)
        if (flag == 0) allocate(df(neqn), stat = flag)
        if (flag == 0) allocate(fvold(neqn), stat = flag)
        if (flag == 0) allocate(xold(nvar), stat = flag)
        if (flag == 0) allocate(dx(nvar), stat = flag)
        if (flag == 0) allocate(s(neqn), stat = flag)
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("qns_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if
        call qr_factor(r, tau, work, lw1)
        call form_qr(r, tau, q, work, lw2)
        call fcn%jacobian(x, b, fv = fvec, olwork = lw3)
        allocate(work(max(lw1, lw2, lw3)), stat = flag)
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("qns_solve", &
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
            ! Determine the maximum line search step
            stpmax = factor * max(norm2(x), real(nvar, real64))

            ! Main Iteration Loop
            do
                ! Update the iteration counter
                iter = iter + 1

                ! Compute or update the Jacobian
                if (restart) then
                    ! Compute the Jacobian
                    call fcn%jacobian(x, b, fvec, work)
                    njac = njac + 1

                    ! Compute the QR factorization, and form Q & R
                    r = b ! Copy the Jacobian - we'll need it later
                    call qr_factor(r, tau, work)
                    call form_qr(r, tau, q, work)

                    ! Reset the Jacobian iteration counter
                    jcount = 0
                else
                    ! Apply the rank 1 update to Q and R
                    df = fvec - fvold
                    dx = x - xold
                    x2 = dot_product(dx, dx)

                    ! Compute S = ALPHA * (DF - B * DX)
                    s = (df - matmul(b, dx))
                    call recip_mult_array(x2, s)

                    ! Compute the new Q and R matrices for the rank1 update:
                    ! B' = B + ALPHA * S * DX**T
                    call rank1_update(one, s, dx, b)
                    call qr_rank1_update(q, r, s, dx, work) ! S & DX overwritten

                    ! Increment the counter tracking how many iterations have
                    ! passed since the last Jacobian recalculation
                    jcount = jcount + 1
                end if

                ! Compute GRAD = B**T * F, store in DX
                call mtx_mult(.true., one, b, fvec, zero, dx)

                ! Store FVEC and X
                xold = x
                fvold = fvec
                fold = f

                ! Solve the linear system: B * DX = -F for DX noting that
                ! B = Q * R.  As such, form -Q**T * F, and store in DF
                call mtx_mult(.true., -one, q, fvec, zero, df)

                ! Now we have R * DX = -Q**T * F, and since R is upper
                ! triangular, the solution is readily computed.  The solution
                ! will be stored in the first NVAR elements of DF
                call solve_triangular_system(.true., .false., .true., r, &
                    df(1:nvar))

                ! Ensure the new solution estimate is heading in a sensible
                ! direction.  If not, it is likely time to update the Jacobian
                temp = dot_product(dx, df(1:nvar))
                if (temp >= zero) then
                    restart = .true.
                    if (this%get_print_status()) then
                        call print_status(iter, neval, njac, xnorm, fnorm)
                    end if
                    cycle
                end if

                ! Apply the line search if needed
                if (this%get_use_line_search()) then
                    ! Define the step length for the line search
                    temp = dot_product(df(1:nvar), df(1:nvar))
                    if (temp > stpmax) df(1:nvar) = df(1:nvar) * (stpmax / temp)

                    ! Apply the line search
                    call limit_search_vector(df(1:nvar), stpmax)
                    call ls%search(fcn, xold, dx, df(1:nvar), x, fvec, fold, &
                        f, lib, errmgr)
                    neval = neval + lib%fcn_count
                else
                    ! No line search - just update the solution estimate
                    x = x + df(1:nvar)
                    call fcn%fcn(x, fvec)
                    f = half * dot_product(fvec, fvec)
                    neval = neval + 1
                end if

                ! Test for convergence
                if (lib%converge_on_zero_diff .and. &
                        this%get_use_line_search()) then
                    call test_convergence(x, xold, fvec, dx, .true., xtol, &
                        ftol, gtol, check, xcnvrg, fcnvrg, gcnvrg, xnorm, fnorm)
                else
                    call test_convergence(x, xold, fvec, dx, .false., xtol, &
                        ftol, gtol, check, xcnvrg, fcnvrg, gcnvrg, xnorm, fnorm)
                end if
                if (.not.check) then
                    ! The solution did not converge, figure out why
                    if (gcnvrg) then
                        ! The slope of the gradient is sufficiently close to
                        ! zero to cause issue.
                        if (restart) then
                            ! We've already tried recalculating a new Jacobian,
                            ! issue a warning
                            write(errmsg, '(AI0AE8.3AE8.3)') &
                                "It appears the solution has settled to " // &
                                "a point where the slope of the gradient " // &
                                "is effectively zero.  " // new_line('c') // &
                                "Function evaluations performed: ", neval, &
                                "." // new_line('c') // &
                                "Change in Variable: ", xnorm, &
                                new_line('c') // "Residual: ", fnorm
                            call errmgr%report_warning("nqs_solve", &
                                trim(errmsg), NL_SPURIOUS_CONVERGENCE_ERROR)
                            exit
                        else
                            ! Try computing a new Jacobian
                            restart = .true.
                        end if
                    else
                        ! We have not converged, but we're not stuck with a
                        ! zero slope gradient vector either.  Go ahead and
                        ! continue the iteration process without recomputing
                        ! the Jacobian - unless the user dictates a
                        ! recaclulation.
                        if (jcount >= this%m_jDelta) then
                            restart = .true.
                        else
                            restart = .false.
                        end if
                    end if
                else
                  ! The solution has converged.  It's OK to exit
                  exit
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
            call errmgr%report_error("qns_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    pure module function qns_get_jac_interval(this) result(n)
        class(quasi_newton_solver), intent(in) :: this
        integer(int32) :: n
        n = this%m_jDelta
    end function

! --------------------
    module subroutine qns_set_jac_interval(this, n)
        class(quasi_newton_solver), intent(inout) :: this
        integer(int32), intent(in) :: n
        this%m_jDelta = n
    end subroutine
end submodule
