! nonlin_test_optimize.f90

module nonlin_test_optimize
    use linalg_constants, only : dp, i32
    use nonlin_optimize
    use nonlin_types
    implicit none
    private
    public :: test_nelder_mead_1

contains
! ******************************************************************************
! TEST FUNCTIONS
! ------------------------------------------------------------------------------
! REF: https://en.wikipedia.org/wiki/Test_functions_for_optimization
! ------------------------------------------------------------------------------
    function rosenbrock(x) result(f)
        real(dp), intent(in), dimension(:) :: x
        real(dp) :: f
        f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    end function

! ------------------------------------------------------------------------------

! ******************************************************************************
! NELDER-MEAD METHOD
! ------------------------------------------------------------------------------
    subroutine test_nelder_mead_1()
        ! Local Variables
        type(nelder_mead) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(dp) :: x(2), fout
        type(iteration_behavior) :: ib

        ! Initialization
        fcn => rosenbrock
        call obj%set_fcn(fcn, 2)

        ! Establish convergence properties
        call solver%set_max_fcn_evals(1000)

        ! Define an initial guess
        x = -1.0d0 ! Equivalent to [-1, -1]

        ! Call the solver
        call solver%solve(obj, x, fout, ib)

        ! Display the output
        print '(AF8.3AF8.3A)', "Rosenbrock Minimum: (", x(1), ", ", x(2), ")"
        print '(AE9.3)', "Residual: ", fout
        print '(AI0)', "Iterations: ", ib%iter_count
        print '(AI0)', "Function Evaluations: ", ib%fcn_count
    end subroutine

! ------------------------------------------------------------------------------    
end module
