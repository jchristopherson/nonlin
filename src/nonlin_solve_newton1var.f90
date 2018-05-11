! nonlin_solve_newton1var.f90

submodule (nonlin_solve) nonlin_solve_newton1var
contains
! ------------------------------------------------------------------------------
    module subroutine newt1var_solve(this, fcn, x, lim, f, ib, err)
        ! Arguments
        class(newton_1var_solver), intent(inout) :: this
        class(fcn1var_helper), intent(in) :: fcn
        real(real64), intent(inout) :: x
        type(value_pair), intent(in) :: lim
        real(real64), intent(out), optional :: f
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: p5 = 0.5d0
        real(real64), parameter :: two = 2.0d0

        ! Local Variables
        logical :: fcnvrg, xcnvrg, dcnvrg
        integer(int32) :: neval, ndiff, maxeval, flag, iter
        real(real64) :: ftol, xtol, dtol, xh, xl, fh, fl, x1, x2, eps, dxold, &
            dx, df, temp, ff
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg

        ! Initialization
        fcnvrg = .false.
        xcnvrg = .false.
        dcnvrg = .false.
        neval = 0
        ndiff = 0
        iter = 0
        ftol = this%get_fcn_tolerance()
        xtol = this%get_var_tolerance()
        dtol = this%get_diff_tolerance()
        maxeval = this%get_max_fcn_evals()
        if (present(f)) f = zero
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = ndiff
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = dcnvrg
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        x1 = min(lim%x1, lim%x2)
        x2 = max(lim%x1, lim%x2)
        eps = epsilon(eps)

        ! Input Check
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("brent_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (abs(x1 - x2) < eps) then
            ! ERROR: Search limits are too tight
            write(errmsg, '(AE8.3AE8.3)') "Search limits have no " // &
                "appreciable difference between them.  Lower Limit: ", x1, &
                ", Upper Limit: ", x2
            call errmgr%report_error("brent_solve", trim(errmsg), &
                NL_INVALID_OPERATION_ERROR)
            return
        end if

        ! See if the root is one of the end points
        flag = 0
        fl = fcn%fcn(x1)
        fh = fcn%fcn(x2)
        neval = 2
        if (abs(fl) < ftol) then
            x = x1
            if (present(f)) f = fl
            if (present(ib)) then
                ib%converge_on_fcn = .true.
                ib%fcn_count = 2
            end if
            return
        end if
        if (abs(fh) < ftol) then
            x = x2
            if (present(f)) f = fh
            if (present(ib)) then
                ib%converge_on_fcn = .true.
                ib%fcn_count = 2
            end if
            return
        end if

        ! Process
        if (fl < zero) then
            xl = x1
            xh = x2
        else
            xl = x2
            xh = x1
        end if
        x = p5 * (x1 + x2)
        dxold = abs(x2 - x1)
        dx = dxold
        ff = fcn%fcn(x)
        df = fcn%diff(x, ff)
        neval = neval + 1
        ndiff = ndiff + 1
        do
            ! Increment the iteration counter
            iter = iter + 1

            ! Bisect if the Newton step went out of range, or if the rate
            ! of change was too slow
            if ((((x - xh) * df - ff) * ((x - xl) * df - ff) > zero) .or. &
                (abs(two * ff) > abs(dxold * df))) &
            then
                ! Bisection
                dxold = dx
                dx = p5 * (xh - xl)
                x = xl + dx
                if (abs(xl - x) < xtol) then
                    ! Convergence as the change in root is within tolerance
                    xcnvrg = .true.
                    exit
                end if
            else
                ! Newton's Method
                dxold = dx
                dx = ff / df
                temp = x
                x = x - dx
                if (abs(temp - x) < xtol) then
                    ! Convergence as the change in root is within tolerance
                    xcnvrg = .true.
                    exit
                end if
            end if

            ! Update function values
            ff = fcn%fcn(x)
            df = fcn%diff(x, ff)
            neval = neval + 1
            ndiff = ndiff + 1

            ! Check for convergence
            if (abs(ff) < ftol) then
                fcnvrg = .true.
                exit
            end if
            if (abs(dx) < xtol) then
                xcnvrg = .true.
                exit
            end if
            if (abs(df) < dtol) then
                dcnvrg = .true.
                exit
            end if

            ! Update the bracket on the root
            if (ff < zero) then
                xl = x
            else
                xh = x
            end if

            ! Print status
            if (this%get_print_status()) then
                call print_status(iter, neval, ndiff, dx, ff)
            end if

            ! Ensure we haven't made too many function evaluations
            if (neval >= maxeval) then
                flag = 1
                exit
            end if
        end do

        ! Ensure the function value is current with the estimate of the root
        if (present(f)) then
            f = fcn%fcn(x)
            neval = neval + 1
        end if

        ! Report out iteration statistics and other optional outputs
        if (present(f)) f = ff
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = ndiff
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = dcnvrg
        end if

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, '(AI0AE8.3AE8.3)') "The algorithm failed to " // &
                "converge.  Function evaluations performed: ", neval, &
                "." // new_line('c') // "Root estimate: ", x, &
                new_line('c') // "Residual: ", ff
            call errmgr%report_error("newt1var_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if
    end subroutine
end submodule
