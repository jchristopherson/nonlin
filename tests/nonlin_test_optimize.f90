! MIT License
!
! Copyright (c) 2017-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
! nonlin_test_optimize.f90

module nonlin_test_optimize
    use iso_fortran_env
    use nonlin
    use autodiff, only : dual, scalar_function, operator(+), operator(-), &
        operator(*), operator(**)
    use fortran_test_helper
    implicit none
    private
    public :: test_nelder_mead_1
    public :: test_nelder_mead_2
    public :: test_nelder_mead_3
    public :: test_bfgs_1
    public :: test_bfgs_2
    public :: test_bfgs_3
    public :: test_nelder_mead_ad
    public :: test_bfgs_ad

contains
! ******************************************************************************
! TEST FUNCTIONS
! ------------------------------------------------------------------------------
! REF: https://en.wikipedia.org/wiki/Test_functions_for_optimization
! ------------------------------------------------------------------------------
    function rosenbrock(x, args) result(f)
        real(real64), intent(in), dimension(:) :: x
        class(*), intent(in), optional :: args
        real(real64) :: f
        f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    end function

! ------------------------------------------------------------------------------
    function rosenbrock2(x, args) result(f)
        real(real64), intent(in), dimension(:) :: x
        class(*), intent(in), optional :: args
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
        class(*), intent(in), optional :: args
        real(real64) :: f
        f = (1.5d0 - x(1) + x(1) * x(2))**2 + &
            (2.25d0 - x(1) + x(1) * x(2)**2)**2 + &
            (2.625d0 - x(1) + x(1) * x(2)**3)**2
    end function

! ------------------------------------------------------------------------------
    function rosenbrock_ad(x, args) result(f)
        type(dual), intent(in), dimension(:) :: x
        class(*), intent(in), optional :: args
        type(dual) :: f
        f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
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
    function test_nelder_mead_ad() result(check)
        real(real64), parameter :: tol = 1.0d-5
        type(nelder_mead) :: solver
        type(fcnnvar_helper) :: obj
        procedure(scalar_function), pointer :: fcn
        real(real64) :: x(2), fout, xans(2)
        logical :: check

        fcn => rosenbrock_ad
        call obj%set_auto_fcn(fcn, 2)
        x = 0.0d0
        call solver%solve(obj, x, fout)

        xans = 1.0d0
        check = assert(x, xans, tol) .and. abs(fout) <= tol
        if (.not.check) then
            print 100, "Test Failed: Nelder-Mead AD Test - Rosenbrock Function"
            print 101, "Expected: (", xans(1), ", ", xans(2), ")"
            print 101, "Computed: (", x(1), ", ", x(2), ")"
            print 102, "Function Value: ", fout
        end if

100     format(A)
101     format(A, F8.5, A, F8.5, A)
102     format(A, E12.4)
    end function

! ------------------------------------------------------------------------------
    function test_bfgs_ad() result(check)
        real(real64), parameter :: tol = 1.0d-10
        type(bfgs) :: solver
        type(fcnnvar_helper) :: obj
        procedure(scalar_function), pointer :: fcn
        real(real64) :: x(2), g(2), xans(2), gans(2), fout
        logical :: check

        fcn => rosenbrock_ad
        call obj%set_auto_fcn(fcn, 2)

        x = 0.0d0
        call obj%gradient(x, g)
        gans = [-2.0d0, 0.0d0]
        check = assert(g, gans, tol)
        if (.not.check) then
            print 100, "Test Failed: BFGS AD Gradient Test"
            print 101, "Expected: (", gans(1), ", ", gans(2), ")"
            print 101, "Computed: (", g(1), ", ", g(2), ")"
            return
        end if

        x = 0.0d0
        call solver%set_use_line_search(.true.)
        call solver%solve(obj, x, fout)
        xans = 1.0d0
        check = assert(x, xans, 1.0d-5) .and. abs(fout) <= 1.0d-10
        if (.not.check) then
            print 100, "Test Failed: BFGS AD Test - Rosenbrock Function"
            print 101, "Expected: (", xans(1), ", ", xans(2), ")"
            print 101, "Computed: (", x(1), ", ", x(2), ")"
            print 102, "Function Value: ", fout
        end if

100     format(A)
101     format(A, F8.5, A, F8.5, A)
102     format(A, E12.4)
    end function

! ------------------------------------------------------------------------------
end module
