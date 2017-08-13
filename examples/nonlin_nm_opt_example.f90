! nonlin_nm_opt_example.f90

program example
    use nonlin_optimize, only : nelder_mead
    use nonlin_types, only : dp, i32, fcnnvar, fcnnvar_helper, &
        iteration_behavior
    implicit none

    ! Local Variables
    type(nelder_mead) :: solver
    type(fcnnvar_helper) :: obj
    procedure(fcnnvar), pointer :: fcn
    real(dp) :: x(2), fout
    type(iteration_behavior) :: ib

    ! Initialization
    fcn => rosenbrock
    call obj%set_fcn(fcn, 2)

    ! Define an initial guess - the solution is (1, 1)
    call random_number(x)

    ! Call the solver
    call solver%solve(obj, x, fout, ib)

     ! Display the output
     print '(AF7.5AF7.5A)', "Minimum: (", x(1), ", ", x(2), ")"
     print '(AE9.3)', "Function Value: ", fout
     print '(AI0)', "Iterations: ", ib%iter_count
     print '(AI0)', "Function Evaluations: ", ib%fcn_count
contains
    ! Rosenbrock's Function
    function rosenbrock(x) result(f)
        real(dp), intent(in), dimension(:) :: x
        real(dp) :: f
        f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    end function
end program