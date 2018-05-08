! nonlin_newton1d_example.f90

program example
    use iso_fortran_env
    use nonlin_core
    use nonlin_solve
    implicit none

    ! Local variables
    type(fcn1var_helper) :: obj
    procedure(fcn1var), pointer :: fcn
    type(newton_1var_solver) :: solver
    real(real64) :: x, f
    type(value_pair) :: limits
    type(iteration_behavior) :: tracking

    ! Define the search limits
    limits%x1 = 1.5d0
    limits%x2 = 5.0d0

    ! Establish the function
    fcn => fcn1
    call obj%set_fcn(fcn)

    ! Solve the equation
    call solver%solve(obj, x, limits, f, ib = tracking)

    ! Print the output and the residual
    print '(AF7.5)', "The solution: ", x
    print '(AE10.3)', "The residual: ", f
    print '(AI0)', "Iterations: ", tracking%iter_count
    print '(AI0)', "Function Evaluations: ", tracking%fcn_count
    print '(AI0)', "Derivative Evaluations: ", tracking%jacobian_count
    print '(AL1)', "Converge on Function Value: ", tracking%converge_on_fcn
    print '(AL1)', "Converge on Change in Variable: ", tracking%converge_on_chng
    print '(AL1)', "Converge on Derivative: ", tracking%converge_on_zero_diff

contains
    ! The function:
    ! f(x) = sin(x) / x, solution: x = x * pi for n = 1, 2, 3, ...
    function fcn1(x) result(f)
        real(real64), intent(in) :: x
        real(real64) :: f
        f = sin(x) / x
    end function
end program
