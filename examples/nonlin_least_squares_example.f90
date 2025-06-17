! nonlin_least_squares_example.f90

program example
    use iso_fortran_env
    use nonlin
    use example_problems
    implicit none

    ! Local Variables
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    type(least_squares_solver) :: solver
    real(real64) :: x(4), f(21) ! There are 4 coefficients and 21 data points

    ! Locate the routine containing the equations to solve
    fcn => lsq_poly_fit_fcn
    call obj%set_fcn(fcn, 21, 4)

    ! Define an initial guess
    x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0, 1.0d0, 1.0d0]

    ! Solve
    call solver%solve(obj, x, f)

    ! Display the output
    print 100, "c0: ", x(4)
    print 100, "c1: ", x(3)
    print 100, "c2: ", x(2)
    print 100, "c3: ", x(1)
    print 101, "Max Residual: ", maxval(abs(f))

    ! Formatting
100 format(A, F12.10)
101 format(A, F7.5)
end program
