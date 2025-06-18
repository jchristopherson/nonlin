! nonlin_bfgs_example.f90

program example
    use iso_fortran_env
    use nonlin
    use example_problems
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
    print 100, "Minimum: (", x(1), ", ", x(2), ")"
    print 101, "Function Value: ", f

    ! Formatting
100 format(A, F7.5, A, F7.5, A)
101 format(A, E9.3)
end program
