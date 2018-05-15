! nonlin_bfgs_example.f90

program example
    use iso_fortran_env
    use nonlin_core
    use nonlin_optimize
    implicit none

    ! Local Variables
    type(bfgs) :: solver
    type(fcnnvar_helper) :: obj
    procedure(fcnnvar), pointer :: fcn
    procedure(gradientfcn), pointer :: grad
    real(real64) :: x(2), f

    ! Tell the solver where to find the function
    fcn => beale
    grad => bealegrad
    call obj%set_fcn(fcn, 2)
    call obj%set_gradient_fcn(grad)

    ! Define an initial guess
    x = 1.0d0

    ! Compute the solution
    call solver%solve(obj, x, f)

    ! Display the output
    print '(AF7.5AF7.5A)', "Minimum: (", x(1), ", ", x(2), ")"
    print '(AE9.3)', "Function Value: ", f
contains
    ! The Beale function:
    ! f(x) = (1.5 - x + xy)**2 + (2.25 - x + xy**2)**2 + (2.625 - x + xy**3)**2
    ! The minimum is at x = 3, y = 0.5, and f(3, 0.5) = 0
    function beale(x) result(f)
        real(real64), intent(in), dimension(:) :: x
        real(real64) :: f
        f = (1.5d0 - x(1) + x(1) * x(2))**2 + &
            (2.25d0 - x(1) + x(1) * x(2)**2)**2 + &
            (2.625d0 - x(1) + x(1) * x(2)**3)**2
    end function

    ! The gradient
    subroutine bealegrad(x, g)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: g

        g(1) = 2.0d0 * (x(2)**3 - 1.0d0) * (x(1) * x(2)**3 - x(1) + 2.625d0) + &
            2.0d0 * (x(2)**2 - 1.0d0) * (x(1) * x(2)**2 - x(1) + 2.25d0) + &
            2.0d0 * (x(2) - 1.0d0) * (x(1) * x(2) - x(1) + 1.5d0)

        g(2) = 6.0d0 * x(1) * x(2)**2 * (x(1) * x(2)**3 - x(1) + 2.625d0) + &
            4.0d0 * x(1) * x(2) * (x(1) * x(2)**2 - x(1) + 2.25d0) + &
            2.0d0 * x(1) * (x(1) * x(2) - x(1) + 1.5d0)
    end subroutine
end program
