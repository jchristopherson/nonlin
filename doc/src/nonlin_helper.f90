module nonlin_helper
    use iso_fortran_env
    implicit none
    private
    public :: print_status
    public :: test_convergence

contains
! ------------------------------------------------------------------------------
    !> @brief Prints the iteration status.
    !!
    !! @param[in] iter The iteration number.
    !! @param[in] nfeval The number of function evaluations.
    !! @param[in] njaceval The number of Jacobian evaluations.
    !! @param[in] xnorm The change in variable value.
    !! @param[in] fnorm The residual.
    subroutine print_status(iter, nfeval, njaceval, xnorm, fnorm)
        ! Arguments
        integer(int32), intent(in) :: iter, nfeval, njaceval
        real(real64), intent(in) :: xnorm, fnorm

        ! Process
        print *, ""
        print 100, "Iteration: ", iter
        print 100, "Function Evaluations: ", nfeval
        if (njaceval > 0) print 100, "Jacobian Evaluations: ", njaceval
        print 101, "Change in Variable: ", xnorm
        print 101, "Residual: ", fnorm

        ! Formatting
100     format(A, I0)
101     format(A, E10.3)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine test_convergence(x, xo, f, g, lg, xtol, ftol, gtol, c, cx, cf, &
            cg, xnorm, fnorm)
        !! Tests for convergence.
        real(real64), intent(in), dimension(:) :: x
            !! The current solution estimate.
        real(real64), intent(in), dimension(:) :: xo
            !! The previous solution estimate.
        real(real64), intent(in), dimension(:) :: f
            !! The current residual based upon x.
        real(real64), intent(in), dimension(:) :: g
            !! The current estimate of the gradient vector at x.
        logical, intent(in) :: lg
            !! Set to true if the gradient slope check should be performed; 
            !! else, false.
        real(real64), intent(in) :: xtol
            !! The tolerance on the change in variable.
        real(real64), intent(in) :: ftol
            !! The tolerance on the residual.
        real(real64), intent(in) :: gtol
            !! The tolerance on the slope of the gradient.
        logical, intent(out) :: c
            !! True if the solution converged on either the residual or
            !! change in variable.
        logical, intent(out) :: cx
            !! True if convergence occurred due to change in variable.
        logical, intent(out) :: cf
            !! True if convergence occurred due to residual.
        logical, intent(out) :: cg
            !! True if convergence occured due to slope of the gradient.
        real(real64), intent(out) :: xnorm
            !! The largest magnitude component of the scaled change in variable 
            !! vector.
        real(real64), intent(out) :: fnorm
            !! The largest magnitude residual component.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: half = 0.5d0

        ! Local Variables
        integer(int32) :: i, nvar, neqn
        real(real64) :: test, dxmax, fc, den

        ! Initialization
        nvar = size(x)
        neqn = size(f)
        cx = .false.
        cf = .false.
        cg = .false.
        c = .false.
        fc = half * dot_product(f, f)
        fnorm = zero
        xnorm = zero

        ! Test for convergence on residual
        do i = 1, neqn
            fnorm = max(abs(f(i)), fnorm)
        end do
        if (fnorm < ftol) then
            cf = .true.
            c = .true.
            return
        end if

        ! Test the change in solution
        do i = 1, nvar
            test = abs(x(i) - xo(i)) / max(abs(x(i)), one)
            xnorm = max(test, xnorm)
        end do
        if (xnorm < xtol) then
            cx = .true.
            c = .true.
            return
        end if

        ! Test for zero gradient slope - do not set convergence flag
        if (lg) then
            test = zero
            den = max(fc, half * nvar)
            do i = 1, nvar
                dxmax = abs(g(i)) * max(abs(x(i)), one) / den
                test = max(test, dxmax)
            end do
            if (test < gtol) then
                cg = .true.
            end if
        end if
    end subroutine

! ------------------------------------------------------------------------------
end module