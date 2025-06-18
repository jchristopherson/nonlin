! nonlin_quasi_newton_example.f90

program  example
    use iso_fortran_env
    use nonlin
    use example_problems
    implicit none

    ! Local Variables
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    type(iteration_behavior) :: ib
    type(quasi_newton_solver) :: solver
    real(real64) :: x(2), f(2)

    ! Locate the routine containing the equations to solve
    fcn => misc_2fcn
    call obj%set_fcn(fcn, 2, 2)

    ! Define an initial guess
    x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]

    ! Defining solver parameters.  This step is optional as the defaults are
    ! typically sufficient; however, this is being done for illustration
    ! purposes.
    !
    ! Establish how many iterations are allowed to pass before the solver
    ! forces a re-evaluation of the Jacobian matrix.  Notice, the solver may
    ! choose to re-evaluate the Jacobian sooner than this, but that is
    ! dependent upon the behavior of the problem.
    call solver%set_jacobian_interval(20)

    ! Establish convergence criteria.  Again, this step is optional as the
    ! defaults are typically sufficient; however, this is being done for
    ! illustration purposes.
    call solver%set_fcn_tolerance(1.0d-8)
    call solver%set_var_tolerance(1.0d-12)
    call solver%set_gradient_tolerance(1.0d-12)

    ! Solve
    call solver%solve(obj, x, f, ib)

    ! Display the output
    print 100, "Solution: (", x(1), ", ", x(2), ")"
    print 101, "Residual: (", f(1), ", ", f(2), ")"
    print 102, "Iterations: ", ib%iter_count
    print 102, "Function Evaluations: ", ib%fcn_count
    print 102, "Jacobian Evaluations: ", ib%jacobian_count

    ! Formatting
100 format(A, F7.5, A, F7.5, A)
101 format(A, E9.3, A, E9.3, A)
102 format(A, I0)
end program
