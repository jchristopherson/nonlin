! nonlin_newton_example.f90

program example
    use iso_fortran_env
    use nonlin_core
    use nonlin_solve
    implicit none

    ! Local Variables
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    type(newton_solver) :: solver
    real(real64) :: x(2), f(2)

    ! Assign a pointer to the subroutine containing the equations to solve
    fcn => fcns
    call obj%set_fcn(fcn, 2, 2) ! There are 2 equations with 2 unknowns

    ! Define an initial guess
    x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]

    ! Solve the system of equations
    call solver%solve(obj, x, f)

    ! Display the output
    print '(AF7.5AF7.5A)', "Solution: (", x(1), ", ", x(2), ")"
    print '(AE9.3AE9.3A)', "Residual: (", f(1), ", ", f(2), ")"
contains
    ! Define the routine containing the equations to solve.  The equations are:
    ! x**2 + y**2 = 34
    ! x**2 - 2 * y**2 = 7
    subroutine fcns(x, f)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    end subroutine
end program
