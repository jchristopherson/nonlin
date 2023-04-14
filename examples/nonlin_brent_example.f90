! nonlin_brent_example.f90

program example
    use iso_fortran_env
    use nonlin_core
    use nonlin_solve
    use example_problems
    implicit none

    ! Local variables
    type(fcn1var_helper) :: obj
    procedure(fcn1var), pointer :: fcn
    type(brent_solver) :: solver
    real(real64) :: x, f
    type(value_pair) :: limits

    ! Define the search limits
    limits%x1 = 1.5d0
    limits%x2 = 5.0d0

    ! Establish the function
    fcn => sinx_div_x
    call obj%set_fcn(fcn)

    ! Solve the equation
    call solver%solve(obj, x, limits, f)

    ! Print the output and the residual
    print 100, "The solution: ", x
    print 101, "The residual: ", f

    ! Formatting
100 format(A, F7.5)
101 format(A, E10.3)
end program
