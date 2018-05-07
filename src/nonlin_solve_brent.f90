! nonlin_solve_brent.f90

submodule (nonlin_solve) nonlin_solve_brent
contains
! ------------------------------------------------------------------------------
    module subroutine brent_solve(this, fcn, x, lim, f, ib, err)
        ! Arguments
        class(brent_solver), intent(inout) :: this
        class(fcn1var_helper), intent(in) :: fcn
        real(real64), intent(inout) :: x
        type(value_pair), intent(in) :: lim
        real(real64), intent(out), optional :: f
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: half = 0.5d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: two = 2.0d0
        real(real64), parameter :: three = 3.0d0

        ! Local Variables
        logical :: fcnvrg, xcnvrg
        integer(int32) :: neval, maxeval, flag, iter
        real(real64) :: ftol, xtol, a, b, c, fa, fb, fc, p, q, r, s, xm, e, d, &
            mn1, mn2, eps, tol1, temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg

        ! Initialization
        fcnvrg = .false.
        xcnvrg = .false.
        x = zero
        a = min(lim%x1, lim%x2)
        b = max(lim%x1, lim%x2)
        neval = 0
        iter = 0
        eps = epsilon(eps)
        ftol = this%get_fcn_tolerance()
        xtol = this%get_var_tolerance()
        maxeval = this%get_max_fcn_evals()
        if (present(f)) f = zero
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = 0
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = .false.
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("brent_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (abs(a - b) < eps) then
            ! ERROR: Search limits are too tight
            write(errmsg, '(AE8.3AE8.3)') "Search limits have no " // &
                "appreciable difference between them.  Lower Limit: ", a, &
                ", Upper Limit: ", b
            call errmgr%report_error("brent_solve", trim(errmsg), &
                NL_INVALID_OPERATION_ERROR)
            return
        end if

        ! Process
        flag = 0
        fa = fcn%fcn(a)
        fb = fcn%fcn(b)
        neval = 2
        fc = fb
        do
            ! Increment the iteration counter
            iter = iter + 1

            ! Adjust the bounding interval
            if ((fb > zero .and. fc >= zero) .or. &
                    (fb < zero .and. fc < zero)) then
                c = a
                fc = fa
                d = b - a
                e = d
            end if
            if (abs(fc) < abs(fb)) then
                a = b
                b = c
                c = a
                fa = fb
                fb = fc
                fc = fa
            end if

            ! Convergence Check
            tol1 = two * eps * abs(b) + half * xtol
            xm = half * (c - b)
            if (abs(fb) < ftol) then
                x = b
                fcnvrg = .true.
                exit
            end if
            if (abs(xm) <= tol1) then
                x = b
                xcnvrg = .true.
                exit
            end if

            ! Actual Method
            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
                ! Attempt the inverse quadratic interpolation to determine
                ! the root
                s = fb / fa
                if (abs(a - c) < eps) then ! a == c
                    p = two * xm * s
                    q = one - s
                else
                    q = fa / fc
                    r = fb / fc
                    p = s * (two * xm * q * (q - r) - (b - a) * (r - one))
                    q = (q - one) * (r - one) * (s - one)
                end if

                ! Ensure we're within bounds
                if (p > zero) q = -q
                p = abs(p)
                mn1 = three * xm * q - abs(tol1 * q)
                mn2 = abs(e * q)
                if (mn1 < mn2) then
                    temp = mn1
                else
                    temp = mn2
                end if
                if (two * p < temp) then
                    ! Accept the interpolation
                    e = d
                    d = p / q
                else
                    ! The interpolation failed, use bisection
                    d = xm
                    e = d
                end if
            else
                ! The bounds are decreasing too slowly, use bisection
                d = xm
                e = d
            end if

            ! Move the last best guess to the lower limit parameter (A)
            a = b
            fa = fb
            if (abs(d) > tol1) then
                b = b + d
            else
                b = b + sign(tol1, xm)
            end if
            fb = fcn%fcn(b)
            neval = neval + 1

            ! Print iteration status
            if (this%get_print_status()) then
                call print_status(iter, neval, 0, xm, fb)
            end if

            ! Ensure we haven't made too many function evaluations
            if (neval >= maxeval) then
                flag = 1
                exit
            end if
        end do

        ! Report out iteration statistics and other optional outputs
        if (present(f)) f = fb
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = 0
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = .false.
        end if

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, '(AI0AE8.3AE8.3)') "The algorithm failed to " // &
                "converge.  Function evaluations performed: ", neval, &
                "." // new_line('c') // "Change in Variable: ", xm, &
                new_line('c') // "Residual: ", fb
            call errmgr%report_error("brent_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if
    end subroutine
end submodule
