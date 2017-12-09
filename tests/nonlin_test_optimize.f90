! nonlin_test_optimize.f90

module nonlin_test_optimize
    use nonlin_optimize
    use nonlin_types
    use test_core
    implicit none
    private
    public :: test_nelder_mead_1
    public :: test_nelder_mead_2
    public :: test_bfgs_1
    public :: test_bfgs_2

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
    function beale(x) result(f)
        real(dp), intent(in), dimension(:) :: x
        real(dp) :: f
        f = (1.5d0 - x(1) + x(1) * x(2))**2 + &
            (2.25d0 - x(1) + x(1) * x(2)**2)**2 + &
            (2.625d0 - x(1) + x(1) * x(2)**3)**2
    end function

! ******************************************************************************
! NELDER-MEAD METHOD
! ------------------------------------------------------------------------------
    function test_nelder_mead_1() result(check)
        ! Parameters
        real(dp), parameter :: tol = 1.0d-6

        ! Local Variables
        type(nelder_mead) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(dp) :: x(2), fout, xans(2)
        type(iteration_behavior) :: ib
        logical :: check

        ! Initialization
        fcn => rosenbrock
        call obj%set_fcn(fcn, 2)

        ! Define an initial guess - the solution is (1, 1)
        call random_number(x)

        ! Call the solver
        call solver%solve(obj, x, fout, ib)

        ! Test
        xans = 1.0d0
        if (is_mtx_equal(x, xans, tol)) then
            check = .true.
        else
            check = .false.
            print '(A)', "Test Failed: Nelder-Mead Test - Rosenbrock Function"
            print '(AF8.5AF8.5A)', "Expected: (", xans(1), ", ", xans(2), ")"
            print '(AF8.5AF8.5A)', "Computed: (", x(1), ", ", x(2), ")"
        end if

        ! ! Display the output
        ! print '(AF8.5AF8.5A)', "Rosenbrock Minimum: (", x(1), ", ", x(2), ")"
        ! print '(AE9.3)', "Function Value: ", fout
        ! print '(AI0)', "Iterations: ", ib%iter_count
        ! print '(AI0)', "Function Evaluations: ", ib%fcn_count
    end function

! ------------------------------------------------------------------------------
    function test_nelder_mead_2() result(check)
        ! Parameters
        real(dp), parameter :: tol = 1.0d-5

        ! Local Variables
        type(nelder_mead) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(dp) :: x(2), fout, xans(2)
        type(iteration_behavior) :: ib
        logical :: check

        ! Initialization
        fcn => beale
        call obj%set_fcn(fcn, 2)

        ! Define an initial guess
        call random_number(x)

        ! Call the solver
        call solver%solve(obj, x, fout, ib)

        ! Test
        xans = [3.0d0, 0.5d0]
        if (is_mtx_equal(x, xans, tol)) then
            check = .true.
        else
            check = .false.
            print '(A)', "Test Failed: Nelder-Mead Test - Beale's Function"
            print '(AF8.5AF8.5A)', "Expected: (", xans(1), ", ", xans(2), ")"
            print '(AF8.5AF8.5A)', "Computed: (", x(1), ", ", x(2), ")"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_bfgs_1() result(check)
        ! Parameters
        real(dp), parameter :: tol = 1.0d-5

        ! Local Variables
        type(bfgs) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(dp) :: x(2), fout, xans(2)
        type(iteration_behavior) :: ib
        logical :: check

        ! Initialization
        fcn => rosenbrock
        call obj%set_fcn(fcn, 2)

        ! Define an initial guess - the solution is (1, 1)
        call random_number(x)

        ! Call the solver
        call solver%solve(obj, x, fout, ib)

        ! Test
        xans = 1.0d0
        if (is_mtx_equal(x, xans, tol)) then
            check = .true.
        else
            check = .false.
            print '(A)', "Test Failed: BFGS Test - Rosenbrock Function"
            print '(AF8.5AF8.5A)', "Expected: (", xans(1), ", ", xans(2), ")"
            print '(AF8.5AF8.5A)', "Computed: (", x(1), ", ", x(2), ")"
        end if

        ! ! Display the output
        ! print '(AF8.5AF8.5A)', "Rosenbrock Minimum: (", x(1), ", ", x(2), ")"
        ! print '(AE9.3)', "Function Value: ", fout
        ! print '(AI0)', "Iterations: ", ib%iter_count
        ! print '(AI0)', "Function Evaluations: ", ib%fcn_count
    end function

! ------------------------------------------------------------------------------
    function test_bfgs_2() result(check)
        ! Parameters
        real(dp), parameter :: tol = 1.0d-5

        ! Local Variables
        type(bfgs) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(dp) :: x(2), fout, xans(2)
        type(iteration_behavior) :: ib
        logical :: check

        ! Initialization
        fcn => beale
        call obj%set_fcn(fcn, 2)

        ! Define an initial guess
        call random_number(x)

        ! Call the solver
        call solver%solve(obj, x, fout, ib)

        ! Test
        xans = [3.0d0, 0.5d0]
        if (is_mtx_equal(x, xans, tol)) then
            check = .true.
        else
            check = .false.
            print '(A)', "Test Failed: BFGS Test - Beale's Function"
            print '(AF8.5AF8.5A)', "Expected: (", xans(1), ", ", xans(2), ")"
            print '(AF8.5AF8.5A)', "Computed: (", x(1), ", ", x(2), ")"
        end if
    end function

! ------------------------------------------------------------------------------
end module
