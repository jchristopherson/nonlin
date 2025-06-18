! nonlin_test_optimize.f90

module nonlin_test_optimize
    use iso_fortran_env
    use nonlin
    use fortran_test_helper
    implicit none
    private
    public :: test_nelder_mead_1
    public :: test_nelder_mead_2
    public :: test_nelder_mead_3
    public :: test_bfgs_1
    public :: test_bfgs_2
    public :: test_bfgs_3

contains
! ******************************************************************************
! TEST FUNCTIONS
! ------------------------------------------------------------------------------
! REF: https://en.wikipedia.org/wiki/Test_functions_for_optimization
! ------------------------------------------------------------------------------
    function rosenbrock(x, args) result(f)
        real(real64), intent(in), dimension(:) :: x
        class(*), intent(inout), optional :: args
        real(real64) :: f
        f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    end function

! ------------------------------------------------------------------------------
    function rosenbrock2(x, args) result(f)
        real(real64), intent(in), dimension(:) :: x
        class(*), intent(inout), optional :: args
        real(real64) :: f, a
        select type (args)
        type is (real(real64))
            a = args
        end select
        f = a * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    end function

! ------------------------------------------------------------------------------
    function beale(x, args) result(f)
        real(real64), intent(in), dimension(:) :: x
        class(*), intent(inout), optional :: args
        real(real64) :: f
        f = (1.5d0 - x(1) + x(1) * x(2))**2 + &
            (2.25d0 - x(1) + x(1) * x(2)**2)**2 + &
            (2.625d0 - x(1) + x(1) * x(2)**3)**2
    end function

! ******************************************************************************
! NELDER-MEAD METHOD
! ------------------------------------------------------------------------------
    function test_nelder_mead_1() result(check)
        ! Parameters
        real(real64), parameter :: tol = 1.0d-5

        ! Local Variables
        type(nelder_mead) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(real64) :: x(2), fout, xans(2)
        type(iteration_behavior) :: ib
        logical :: check

        ! Initialization
        fcn => rosenbrock
        call obj%set_fcn(fcn, 2)

        ! Define an initial guess - the solution is (1, 1)
        x = 0.0d0

        ! Call the solver
        call solver%solve(obj, x, fout, ib)

        ! Test
        xans = 1.0d0
        if (assert(x, xans, tol)) then
            check = .true.
        else
            check = .false.
            print 100, "Test Failed: Nelder-Mead Test - Rosenbrock Function"
            print 101, "Expected: (", xans(1), ", ", xans(2), ")"
            print 101, "Computed: (", x(1), ", ", x(2), ")"
        end if

        ! ! Display the output
        ! print 101, "Rosenbrock Minimum: (", x(1), ", ", x(2), ")"
        ! print '(AE9.3)', "Function Value: ", fout
        ! print '(AI0)', "Iterations: ", ib%iter_count
        ! print '(AI0)', "Function Evaluations: ", ib%fcn_count

        ! Formatting
100     format(A)
101     format(A, F8.5, A, F8.5, A)
    end function

! ------------------------------------------------------------------------------
    function test_nelder_mead_2() result(check)
        ! Parameters
        real(real64), parameter :: tol = 1.0d-5

        ! Local Variables
        type(nelder_mead) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(real64) :: x(2), fout, xans(2)
        type(iteration_behavior) :: ib
        logical :: check

        ! Initialization
        fcn => beale
        call obj%set_fcn(fcn, 2)

        ! Define an initial guess
        x = 1.0d0

        ! Call the solver
        call solver%solve(obj, x, fout, ib)

        ! Test
        xans = [3.0d0, 0.5d0]
        if (assert(x, xans, tol)) then
            check = .true.
        else
            check = .false.
            print 100, "Test Failed: Nelder-Mead Test - Beale's Function"
            print 101, "Expected: (", xans(1), ", ", xans(2), ")"
            print 101, "Computed: (", x(1), ", ", x(2), ")"
        end if

        ! Formatting
100     format(A)
101     format(A, F8.5, A, F8.5, A)
    end function

! ------------------------------------------------------------------------------
    function test_nelder_mead_3() result(check)
        ! Parameters
        real(real64), parameter :: tol = 1.0d-5

        ! Local Variables
        type(nelder_mead) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(real64) :: x(2), fout, xans(2), a
        type(iteration_behavior) :: ib
        logical :: check

        ! Initialization
        a = 1.0d2
        fcn => rosenbrock2
        call obj%set_fcn(fcn, 2)

        ! Define an initial guess - the solution is (1, 1)
        x = 0.0d0

        ! Call the solver
        call solver%solve(obj, x, fout, ib, args = a)

        ! Test
        xans = 1.0d0
        if (assert(x, xans, tol)) then
            check = .true.
        else
            check = .false.
            print 100, "Test Failed: Nelder-Mead Test - Rosenbrock Function -2"
            print 101, "Expected: (", xans(1), ", ", xans(2), ")"
            print 101, "Computed: (", x(1), ", ", x(2), ")"
        end if

        ! ! Display the output
        ! print 101, "Rosenbrock Minimum: (", x(1), ", ", x(2), ")"
        ! print '(AE9.3)', "Function Value: ", fout
        ! print '(AI0)', "Iterations: ", ib%iter_count
        ! print '(AI0)', "Function Evaluations: ", ib%fcn_count

        ! Formatting
100     format(A)
101     format(A, F8.5, A, F8.5, A)
    end function

! ------------------------------------------------------------------------------
    function test_bfgs_1() result(check)
        ! Parameters
        real(real64), parameter :: tol = 1.0d-5

        ! Local Variables
        type(bfgs) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(real64) :: x(2), fout, xans(2)
        type(iteration_behavior) :: ib
        logical :: check

        ! Initialization
        fcn => rosenbrock
        call obj%set_fcn(fcn, 2)

        ! Define an initial guess - the solution is (1, 1)
        x = 0.0d0

        ! Call the solver
        call solver%solve(obj, x, fout, ib)

        ! Test
        xans = 1.0d0
        if (assert(x, xans, tol)) then
            check = .true.
        else
            check = .false.
            print 100, "Test Failed: BFGS Test - Rosenbrock Function"
            print 101, "Expected: (", xans(1), ", ", xans(2), ")"
            print 101, "Computed: (", x(1), ", ", x(2), ")"
        end if

        ! ! Display the output
        ! print 101, "Rosenbrock Minimum: (", x(1), ", ", x(2), ")"
        ! print '(AE9.3)', "Function Value: ", fout
        ! print '(AI0)', "Iterations: ", ib%iter_count
        ! print '(AI0)', "Function Evaluations: ", ib%fcn_count

        ! Formatting
100     format(A)
101     format(A, F8.5, A, F8.5, A)
    end function

! ------------------------------------------------------------------------------
    function test_bfgs_2() result(check)
        ! Parameters
        real(real64), parameter :: tol = 1.0d-5

        ! Local Variables
        type(bfgs) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(real64) :: x(2), fout, xans(2)
        type(iteration_behavior) :: ib
        logical :: check

        ! Initialization
        fcn => beale
        call obj%set_fcn(fcn, 2)

        ! Define an initial guess
        x = 1.0d0

        ! Call the solver
        call solver%solve(obj, x, fout, ib)

        ! Test
        xans = [3.0d0, 0.5d0]
        if (assert(x, xans, tol)) then
            check = .true.
        else
            check = .false.
            print 100, "Test Failed: BFGS Test - Beale's Function"
            print 101, "Expected: (", xans(1), ", ", xans(2), ")"
            print 101, "Computed: (", x(1), ", ", x(2), ")"
        end if

        ! Formatting
100     format(A)
101     format(A, F8.5, A, F8.5, A)
    end function

! ------------------------------------------------------------------------------
    function test_bfgs_3() result(check)
        ! Parameters
        real(real64), parameter :: tol = 1.0d-5

        ! Local Variables
        type(bfgs) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        real(real64) :: x(2), fout, xans(2), a
        type(iteration_behavior) :: ib
        logical :: check

        ! Initialization
        a = 1.0d2
        fcn => rosenbrock2
        call obj%set_fcn(fcn, 2)

        ! Define an initial guess - the solution is (1, 1)
        x = 0.0d0

        ! Call the solver
        call solver%solve(obj, x, fout, ib, args = a)

        ! Test
        xans = 1.0d0
        if (assert(x, xans, tol)) then
            check = .true.
        else
            check = .false.
            print 100, "Test Failed: BFGS Test - Rosenbrock Function -2"
            print 101, "Expected: (", xans(1), ", ", xans(2), ")"
            print 101, "Computed: (", x(1), ", ", x(2), ")"
        end if

        ! Formatting
100     format(A)
101     format(A, F8.5, A, F8.5, A)
    end function

! ------------------------------------------------------------------------------
end module
