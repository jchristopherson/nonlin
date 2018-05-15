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
    type(iteration_behavior) :: tracking

    ! Assign a pointer to the subroutine containing the equations to solve
    fcn => fcns
    call obj%set_fcn(fcn, 2, 2) ! There are 2 equations with 2 unknowns

    ! Set up solver parameters
    call solver%set_max_fcn_evals(1000) ! Specify the maximum number of function evaluations before iteration termination
    call solver%set_fcn_tolerance(1.0d-10)
    call solver%set_var_tolerance(1.0d-10)
    call solver%set_gradient_tolerance(1.0d-10)

    ! Tell the solver to print out a status update at each iteration
    call solver%set_print_status(.true.)

    ! Define an initial guess
    x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]

    ! Solve the system of equations, but include solver statistics tracking
    call solver%solve(obj, x, f, tracking)

    ! Display the output
    print *, ""
    print '(AF7.5AF7.5A)', "Solution: (", x(1), ", ", x(2), ")"
    print '(AE9.3AE9.3A)', "Residual: (", f(1), ", ", f(2), ")"
    print '(AI0)', "Iteration Count: ", tracking%iter_count
    print '(AI0)', "Function Evaluations: ", tracking%fcn_count
    print '(AI0)', "Jacobian Evaluations: ", tracking%jacobian_count
    print '(AL1)', "Converge on Function Value: ", tracking%converge_on_fcn
    print '(AL1)', "Converge on Change in Variable: ", tracking%converge_on_chng
    print '(AL1)', "Converge on Zero Slope Gradient Vector: ", tracking%converge_on_zero_diff
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
