! nonlin_newton_solve_jacobiaon.f90

program example
    use iso_fortran_env
    use nonlin_core
    use nonlin_solve
    implicit none

    ! Local Variables
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    procedure(jacobianfcn), pointer :: jac
    type(newton_solver) :: solver
    real(real64) :: x(2), f(2)

    ! Assign the function and Jacobian routines
    fcn => fcns
    jac => fcnjac
    call obj%set_fcn(fcn, 2, 2)
    call obj%set_jacobian(jac)

    ! Define an initial guess
    x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]

    ! Solve the system of equations
    call solver%solve(obj, x, f)

    ! Display the output
    print '(AF7.5AF7.5A)', "Solution: (", x(1), ", ", x(2), ")"
    print '(AE9.3AE9.3A)', "Residual: (", f(1), ", ", f(2), ")"
contains
    ! The system of equations (source: https://www.mathworks.com/help/optim/ug/fsolve.html)
    ! 2 * x1 - x2 = exp(-x1)
    ! -x1 + 2 * x2 = exp(-x2)
    subroutine fcns(x, f)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        f(1) = 2.0d0 * x(1) - x(2) - exp(-x(1))
        f(2) = -x(1) + 2.0d0 * x(2) - exp(-x(2))
    end subroutine

    ! The Jacobian matrix:
    !     | exp(-x1) + 2          -1     |
    ! J = |                              |
    !     |     -1          exp(-x2) + 2 |
    subroutine fcnjac(x, jac)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: jac
        jac(1,1) = exp(-x(1)) + 2.0d0
        jac(2,1) = -1.0d0
        jac(1,2) = -1.0d0
        jac(2,2) = exp(-x(2)) + 2.0d0
    end subroutine
end program
