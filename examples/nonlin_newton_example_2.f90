! nonlin_newton_example_2.f90

program example
    use iso_fortran_env
    use nonlin_core
    use nonlin_solve
    use example_problems
    implicit none

    ! Variables
    type(vecfcn_helper) :: obj
    type(newton_solver) :: solver
    procedure(vecfcn), pointer :: fcn
    real(real64) :: x(2), f(2)

    ! Initialization
    fcn => powell_bad
    call obj%set_fcn(fcn, 2, 2)
    x = [0.0d0, 1.0d0]
    call solver%set_print_status(.true.)

    ! Solve the equations
    call solver%solve(obj, x, f)

    ! Display the output
    print *, ""
    print 100, "Solution: (", x(1), ", ", x(2), ")"
    print 100, "Residual: (", f(1), ", ", f(2), ")"

    ! Formatting
100 format(A, E12.6, A, E12.6, A)
end program
