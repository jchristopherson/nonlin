! nonlin_quasi_newton_example.f90

program  example
    use nonlin_types, only : dp, i32, vecfcn_helper, vecfcn, iteration_behavior
    use nonlin_solve, only : quasi_newton_solver
    implicit none

    ! Local Variables
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    type(iteration_behavior) :: ib
    type(quasi_newton_solver) :: solver
    real(dp) :: x(2), f(2)

    ! Locate the routine containing the equations to solve
    fcn => fcns
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
    print '(AF7.5AF7.5A)', "Solution: (", x(1), ", ", x(2), ")"
    print '(AE9.3AE9.3A)', "Residual: (", f(1), ", ", f(2), ")"
    print '(AI0)', "Iterations: ", ib%iter_count
    print '(AI0)', "Function Evaluations: ", ib%fcn_count
    print '(AI0)', "Jacobian Evaluations: ", ib%jacobian_count

contains
    ! Define the routine containing the equations to solve.  The equations are:
    ! x**2 + y**2 = 34
    ! x**2 - 2 * y**2 = 7
    subroutine fcns(x, f)
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: f
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    end subroutine
end program