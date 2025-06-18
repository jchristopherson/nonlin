! nonlin_newton_solve_jacobiaon.f90

program example
    use iso_fortran_env
    use nonlin
    use example_problems
    implicit none

    ! Local Variables
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    procedure(jacobianfcn), pointer :: jac
    type(newton_solver) :: solver
    real(real64) :: x(2), f(2)

    ! Assign the function and Jacobian routines
    fcn => misc_2fcn_01
    jac => misc_2fcn_01_jac
    call obj%set_fcn(fcn, 2, 2)
    call obj%set_jacobian(jac)

    ! Define an initial guess
    x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]

    ! Solve the system of equations
    call solver%solve(obj, x, f)

    ! Display the output
    print 100, "Solution: (", x(1), ", ", x(2), ")"
    print 101, "Residual: (", f(1), ", ", f(2), ")"

    ! Formatting
100 format(A, F7.5, A, F7.5, A)
101 format(A, E9.3, A, E9.3, A)
end program
