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
! nonlin_test_solve.f90

module nonlin_test_solve
    use iso_fortran_env
    use nonlin
    use, intrinsic :: ieee_arithmetic
    use autodiff, only : dual, vector_function, assignment(=), &
        operator(+), operator(-), operator(*), operator(**)
    implicit none
    private
    public :: test_quasinewton_1
    public :: test_quasinewton_2
    public :: test_quasinewton_3
    public :: test_quasinewton_4
    public :: test_newton_1
    public :: test_newton_2
    public :: test_newton_3
    public :: test_newton_4
    public :: test_least_squares_1
    public :: test_least_squares_2
    public :: test_least_squares_3
    public :: test_least_squares_4
    public :: test_brent_1
    public :: test_brent_2
    public :: test_newton_1var_1
    public :: test_newton_1var_2
    public :: test_scalar_solver_edge_cases
    public :: test_constrained_least_squares_1
    public :: test_constrained_least_squares_2
    public :: test_constrained_least_squares_3
    public :: test_constrained_least_squares_4
    public :: test_constrained_least_squares_bounds
    public :: test_vector_solvers_ad
    public :: test_vector_solver_robustness
    integer :: vector_ad_calls = 0, vector_jac_calls = 0
contains
! ******************************************************************************
! TEST FUNCTIONS
! ------------------------------------------------------------------------------
    ! System of Equations #1:
    !
    ! x**2 + y**2 = 34
    ! x**2 - 2 * y**2 = 7
    !
    ! Solution:
    ! x = +/-5
    ! y = +/-3
    subroutine fcn1(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    end subroutine

    subroutine fcn1a(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        real(real64) :: a
        select type (args)
        type is (real(real64))
            a = args
        end select
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - a * x(2)**2 - 7.0d0
    end subroutine

    ! Jacobian:
    !
    !     | 2x  2y |
    ! J = |        |
    !     | 2x  -4y|
    subroutine jac1(x, j, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: j
        class(*), intent(in), optional :: args
        j = 2.0d0 * reshape([x(1), x(1), x(2), -2.0d0 * x(2)], [2, 2])
    end subroutine

    subroutine jac1a(x, j, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: j
        class(*), intent(in), optional :: args
        real(real64) :: a
        select type (args)
        type is (real(real64))
            a = args
        end select
        j = 2.0d0 * reshape([x(1), x(1), x(2), -a * x(2)], [2, 2])
    end subroutine

    pure function is_ans_1(x, tol) result(c)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(in) :: tol
        logical :: c
        real(real64), parameter :: x1 = 5.0d0
        real(real64), parameter :: x2 = 3.0d0
        real(real64) :: ax1, ax2
        c = .true.
        ax1 = abs(x(1)) - x1
        ax2 = abs(x(2)) - x2
        if (abs(ax1) > tol .or. abs(ax2) > tol) c = .false.
    end function

! ------------------------------------------------------------------------------
    ! System of Equations #2 (Poorly Scaled Problem)
    ! REF: http://folk.uib.no/ssu029/Pdf_file/Hiebert82.pdf
    !
    ! x2 - 10 = 0
    ! x1 * x2 - 5e4 = 0
    !
    ! Solution:
    ! x1 = 5e3
    ! x2 = 10
    subroutine fcn2(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        f(1) = x(2) - 10.0d0
        f(2) = x(1) * x(2) - 5e4
    end subroutine

    pure function is_ans_2(x, tol) result(c)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(in) :: tol
        real(real64) :: ax1, ax2
        logical :: c
        real(real64), parameter :: x1 = 5.0d3
        real(real64), parameter :: x2 = 1.0d1
        c = .true.
        ax1 = abs(x(1)) - x1
        ax2 = abs(x(2)) - x2
        if (abs(ax1) > tol .or. abs(ax2) > tol) c = .false.
    end function

! ******************************************************************************
! LEAST SQUARES FUNCTIONS
! ------------------------------------------------------------------------------
    subroutine lsfcn1(x, f, args)
        ! Arguments
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args

        ! Local Variables
        real(real64), dimension(21) :: xp, yp

        ! Data to fit
        xp = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
            0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
            1.8d0, 1.9d0, 2.0d0]
        yp = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
            1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
            1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
            3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
            5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
            10.30252404d0]

        ! We'll apply a cubic polynomial model to this data:
        ! y = c1 * x**3 + c2 * x**2 + c3 * x + c4
        f = x(1) * xp**3 + x(2) * xp**2 + x(3) * xp + x(4) - yp

        ! For reference, the data was generated by adding random errors to
        ! the following polynomial: y = x**3 - 0.3 * x**2 + 1.2 * x + 0.3
    end subroutine

! ******************************************************************************
! 1 VARIABLE FUNCTIONS
! ------------------------------------------------------------------------------
    ! f(x) = sin(x) / x, SOLUTION: x = n * pi for n = 0, 1, 2, 3, ...
    function f1var_1(x, args) result(f)
        real(real64), intent(in) :: x
        class(*), intent(in), optional :: args
        real(real64) :: f
        f = sin(x) / x
    end function

    function f1var_1a(x, args) result(f)
        real(real64), intent(in) :: x
        class(*), intent(in), optional :: args
        real(real64) :: f
        real(real64) :: a
        select type (args)
        type is (real(real64))
            a = args
        end select
        f = a * sin(x) / x
    end function

    function f1var_shifted(x, args) result(f)
        real(real64), intent(in) :: x
        class(*), intent(in), optional :: args
        real(real64) :: f
        f = x - 1.0d0
    end function

! ******************************************************************************
! SOLVER TEST ROUTINES
! ------------------------------------------------------------------------------
    function test_quasinewton_1() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(quasi_newton_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2)
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        fcn => fcn1
        jac => jac1
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Quasi-Newton Solver Failed: Test 1-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_quasinewton_2() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        type(quasi_newton_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2)
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        fcn => fcn2
        call obj%set_fcn(fcn, 2, 2)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Turn off the line search - this set of functions is too poorly scaled
        ! for the current implementation of the line search algorithm to offer
        ! much help.  This seems to indicate a need for improvement in the
        ! line search code - perhaps variable scaling?
        call solver%set_use_line_search(.false.)

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib)
            if (.not.is_ans_2(x, 1.0d-6)) then
                check = .false.
                print 100, "Quasi-Newton Solver Failed: Test 2-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_quasinewton_3() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(quasi_newton_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2), a
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        a = 2.0d0
        fcn => fcn1a
        jac => jac1a
        call obj%set_fcn(fcn, 2, 2)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib, args = a)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Quasi-Newton Solver Failed: Test 3a-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Now try with a Jacobian matrix
        call obj%set_jacobian(jac)
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib, args = a)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Quasi-Newton Solver Failed: Test 3b-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_newton_1() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(newton_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2)
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        fcn => fcn1
        jac => jac1
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Newton Solver Failed: Test 1-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_newton_2() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        type(newton_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2)
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        fcn => fcn2
        call obj%set_fcn(fcn, 2, 2)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Turn off the line search - this set of functions is too poorly scaled
        ! for the current implementation of the line search algorithm to offer
        ! much help.  This seems to indicate a need for improvement in the
        ! line search code - perhaps variable scaling?
        call solver%set_use_line_search(.false.)

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib)
            if (.not.is_ans_2(x, 1.0d-6)) then
                check = .false.
                print 100, "Newton Solver Failed: Test 2-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_newton_3() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(newton_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2), a
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        a = 2.0d0
        fcn => fcn1a
        jac => jac1a
        call obj%set_fcn(fcn, 2, 2)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib, args = a)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Newton Solver Failed: Test 3a-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Now try with a user-defined Jacobian
        call obj%set_jacobian(jac)
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib, args = a)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Newton Solver Failed: Test 3b-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_least_squares_1() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(least_squares_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2)
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        fcn => fcn1
        jac => jac1
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Least Squares Solver Failed: Test 1-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_least_squares_2() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        type(least_squares_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2)
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        fcn => fcn2
        call obj%set_fcn(fcn, 2, 2)

        ! Increase the number of iterations allowed
        call solver%set_max_fcn_evals(1000)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib)
            if (.not.is_ans_2(x, 1.0d-6)) then
                check = .false.
                print 100, "Least Squares Solver Failed: Test 2-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    subroutine test_least_squares_3()
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        type(least_squares_solver) :: solver
        real(real64) :: x(4), f(21) ! There are 4 coefficients and 21 data points

        ! Initialization
        fcn => lsfcn1
        x = 1.0d0   ! Set X to an initial guess of [1, 1, 1, 1]
        call obj%set_fcn(fcn, 21, 4)

        ! Compute the solution, and store the polynomial coefficients in X
        call solver%solve(obj, x, f)

        ! Print out the coefficients
        !print *, x
    end subroutine

! ------------------------------------------------------------------------------
    function test_least_squares_4() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(least_squares_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2), a
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        a = 2.0d0
        fcn => fcn1
        jac => jac1
        call obj%set_fcn(fcn, 2, 2)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib, args = a)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Least Squares Solver Failed: Test 4a-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Try with a user-defined Jacobian
        call obj%set_jacobian(jac)
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib, args = a)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Least Squares Solver Failed: Test 4b-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_brent_1() result(check)
        ! Local Variables
        type(brent_solver) :: solver
        type(fcn1var_helper) :: obj
        procedure(fcn1var), pointer :: fcn
        real(real64) :: x, f
        type(value_pair) :: limits
        logical :: check

        ! Parameters
        real(real64), parameter :: pi = 3.141592653589793d0
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        fcn => f1var_1
        call obj%set_fcn(fcn)

        ! Define the search limits
        limits%x1 = 1.5d0
        limits%x2 = 5.0d0

        ! Compute the solution
        call solver%solve(obj, x, limits, f)

        ! The solution on this interval should be: pi
        if (abs(x - pi) > tol) then
            check = .false.
            print 100, &
                "Test Failed: Brent's Method Test 1.  Expected: ", pi, &
                ", Found: ", x
        end if

        ! Formatting
100     format(A, F8.5, A, F8.5)
    end function

! ------------------------------------------------------------------------------
    function test_brent_2() result(check)
        ! Local Variables
        type(brent_solver) :: solver
        type(fcn1var_helper) :: obj
        procedure(fcn1var), pointer :: fcn
        real(real64) :: x, f, a
        class(*), allocatable :: args
        type(value_pair) :: limits
        logical :: check

        ! Parameters
        real(real64), parameter :: pi = 3.141592653589793d0
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        a = 2.0d0
        allocate(args, source = a)
        fcn => f1var_1a
        call obj%set_fcn(fcn)

        ! Define the search limits
        limits%x1 = 1.5d0
        limits%x2 = 5.0d0

        ! Compute the solution
        call solver%solve(obj, x, limits, f, args = args)

        ! The solution on this interval should be: pi
        if (abs(x - pi) > tol) then
            check = .false.
            print 100, &
                "Test Failed: Brent's Method Test 2.  Expected: ", pi, &
                ", Found: ", x
        end if

        ! Formatting
100     format(A, F8.5, A, F8.5)
    end function

! ------------------------------------------------------------------------------
    function test_newton_4() result(check)
        use powell_badly_scaled_module

        ! Local Variables
        logical :: check
        type(newton_solver) :: solver
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        real(real64), dimension(2) :: x, f, sol
        integer :: i
        real(real64), parameter :: tol = 1.0d-5

        ! Initialization
        check = .true.

        ! Define the initial conditions
        x = powell_badly_scaled_start()

        ! Set up the solver
        fcn => powell_badly_scaled
        jac => powell_badly_scaled_jacobian
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)

        ! Solve
        call solver%solve(obj, x, f)

        ! Test the solution
        sol = powell_badly_scaled_solution()
        do i = 1, size(sol)
            if (abs(x(i) - sol(i)) > tol) then
                check = .false.
                print 100, "Test Failed: Newton's Method, Test 4."
                print 101, "Expected: ", sol(i), " for root ", &
                    i, ", but found: ", x(i)
            end if
        end do

        ! Formatting
100     format(A)
101     format(A, E12.5, A, I0, A, E12.5)
    end function

! ------------------------------------------------------------------------------
    function test_quasinewton_4() result(check)
        use powell_badly_scaled_module

        ! Local Variables
        logical :: check
        type(quasi_newton_solver) :: solver
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        real(real64), dimension(2) :: x, f, sol
        integer :: i
        real(real64), parameter :: tol = 1.0d-5

        ! Initialization
        check = .true.

        ! Define the initial conditions
        x = powell_badly_scaled_start()

        ! Set up the solver
        fcn => powell_badly_scaled
        jac => powell_badly_scaled_jacobian
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)

        ! Solve - Do not use line search as the line search will fail on this
        ! problem.
        call solver%set_use_line_search(.false.)
        call solver%solve(obj, x, f)

        ! Test the solution
        sol = powell_badly_scaled_solution()
        do i = 1, size(sol)
            if (abs(x(i) - sol(i)) > tol) then
                check = .false.
                print 100, "Test Failed: Quasi-Newton's Method, Test 4."
                print 101, "Expected: ", sol(i), " for root ", &
                    i, ", but found: ", x(i)
            end if
        end do

        ! Formatting
100     format(A)
101     format(A, E12.5, A, I0, A, E12.5)
    end function

! ------------------------------------------------------------------------------
    function test_newton_1var_1() result(check)
        ! Local Variables
        type(newton_1var_solver) :: solver
        type(fcn1var_helper) :: obj
        procedure(fcn1var), pointer :: fcn
        real(real64) :: x, f
        type(value_pair) :: limits
        logical :: check

        ! Parameters
        real(real64), parameter :: pi = 3.141592653589793d0
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        fcn => f1var_1
        call obj%set_fcn(fcn)

        ! Define the search limits
        limits%x1 = 1.5d0
        limits%x2 = 5.0d0

        ! Compute the solution
        call solver%solve(obj, x, limits, f)

        ! The solution on this interval should be: pi
        if (abs(x - pi) > tol) then
            check = .false.
            print 100, &
                "Test Failed: Newton's (1 Variable) Method Test 1.  Expected: ", pi, &
                ", Found: ", x
        end if

        ! Formatting
100     format(A, F8.5, A, F8.5)
    end function

! ------------------------------------------------------------------------------
    function test_newton_1var_2() result(check)
        ! Local Variables
        type(newton_1var_solver) :: solver
        type(fcn1var_helper) :: obj
        procedure(fcn1var), pointer :: fcn
        real(real64) :: x, f, a
        class(*), allocatable :: args
        type(value_pair) :: limits
        logical :: check

        ! Parameters
        real(real64), parameter :: pi = 3.141592653589793d0
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        a = 2.0d0
        allocate(args, source = a)
        fcn => f1var_1a
        call obj%set_fcn(fcn)

        ! Define the search limits
        limits%x1 = 1.5d0
        limits%x2 = 5.0d0

        ! Compute the solution
        call solver%solve(obj, x, limits, f, args = args)

        ! The solution on this interval should be: pi
        if (abs(x - pi) > tol) then
            check = .false.
            print 100, &
                "Test Failed: Newton's (1 Variable) Method Test 2.  Expected: ", pi, &
                ", Found: ", x
        end if

        ! Formatting
100     format(A, F8.5, A, F8.5)
    end function

! ------------------------------------------------------------------------------
    function test_scalar_solver_edge_cases() result(check)
        type(brent_solver) :: brent
        type(newton_1var_solver) :: newton
        type(fcn1var_helper) :: obj
        procedure(fcn1var), pointer :: fcn
        type(value_pair) :: limits
        type(iteration_behavior) :: ib
        real(real64) :: x, f, scale
        class(*), allocatable :: args
        logical :: check

        check = .true.
        fcn => f1var_shifted
        call obj%set_fcn(fcn)

        ! Both interval endpoints should be accepted, even with reversed bounds.
        limits%x1 = 2.0d0
        limits%x2 = 1.0d0
        call brent%solve(obj, x, limits, f, ib)
        if (abs(x - 1.0d0) > 1.0d-12 .or. abs(f) > 1.0d-12 .or. &
                .not.ib%converge_on_fcn) then
            check = .false.
            print *, "Test Failed: Brent solver did not accept a reversed endpoint root."
        end if
        call newton%solve(obj, x, limits, f, ib)
        if (abs(x - 1.0d0) > 1.0d-12 .or. abs(f) > 1.0d-12 .or. &
                .not.ib%converge_on_fcn) then
            check = .false.
            print *, "Test Failed: Newton solver did not accept a reversed endpoint root."
        end if

        ! Small residual scales must not make an interior root numerically invisible.
        fcn => f1var_1a
        call obj%set_fcn(fcn)
        scale = 1.0d-20
        allocate(args, source = scale)
        limits%x1 = 1.5d0
        limits%x2 = 5.0d0
        call brent%set_fcn_tolerance(1.0d-30)
        call newton%set_fcn_tolerance(1.0d-30)
        call newton%set_diff_tolerance(1.0d-30)

        call brent%solve(obj, x, limits, f, ib, args = args)
        if (abs(x - 3.141592653589793d0) > 1.0d-6 .or. &
                abs(f) > 1.0d-30 .or. &
                .not.(ib%converge_on_fcn .or. ib%converge_on_chng)) then
            check = .false.
            print *, "Test Failed: Brent solver was not robust to a small function scale."
        end if

        call newton%solve(obj, x, limits, f, ib, args = args)
        if (abs(x - 3.141592653589793d0) > 1.0d-6 .or. &
                abs(f) > 1.0d-30 .or. &
                .not.(ib%converge_on_fcn .or. ib%converge_on_chng)) then
            check = .false.
            print *, "Test Failed: Newton solver was not robust to a small function scale."
        end if
    end function

! ------------------------------------------------------------------------------
    function test_constrained_least_squares_1() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(constrained_least_squares_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2), big
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        fcn => fcn1
        jac => jac1
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)
        big = huge(big)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Set limits
        call solver%set_upper_limits([big, big])
        call solver%set_lower_limits([-big, -big])

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Constrained Least Squares Solver Failed: Test 1-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_constrained_least_squares_2() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        type(constrained_least_squares_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2)
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        fcn => fcn2
        call obj%set_fcn(fcn, 2, 2)

        ! Increase the number of iterations allowed
        call solver%set_max_fcn_evals(5000)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib)
            if (.not.is_ans_2(x, 1.0d-6)) then
                check = .false.
                print 100, "Constrained Least Squares Solver Failed: Test 2-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_constrained_least_squares_3() result(rst)
        ! Local Variables
        logical :: rst
        real(real64), parameter :: tol = 1.0d-5
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        type(least_squares_solver) :: solver
        type(constrained_least_squares_solver) :: csolver
        real(real64) :: x(4), f(21), xc(4) ! There are 4 coefficients and 21 data points

        ! Initialization
        rst = .true.
        fcn => lsfcn1
        x = 1.0d0   ! Set X to an initial guess of [1, 1, 1, 1]
        xc = 1.0d0
        call obj%set_fcn(fcn, 21, 4)

        ! Compute the solution, and store the polynomial coefficients in X
        call solver%solve(obj, x, f)
        call csolver%solve(obj, xc, f)

        ! Test
        if ( &
            abs(x(1) - xc(1)) > tol .or. &
            abs(x(2) - xc(2)) > tol .or. &
            abs(x(3) - xc(3)) > tol .or. &
            abs(x(4) - xc(4)) > tol &
        ) &
        then
            rst = .false.
            print '(A)', "Constrained Least Squares Test 3 Failed."
        end if
    end function

! ------------------------------------------------------------------------------
    function test_constrained_least_squares_4() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(constrained_least_squares_solver) :: solver
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2), ic(2, 2), a
        integer :: i
        logical :: check

        ! Initialization
        check = .true.
        a = 2.0d0
        fcn => fcn1
        jac => jac1
        call obj%set_fcn(fcn, 2, 2)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib, args = a)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Constrained Least Squares Solver Failed: Test 4a-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Try with a user-defined Jacobian
        call obj%set_jacobian(jac)
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib, args = a)
            if (.not.is_ans_1(x, 1.0d-6)) then
                check = .false.
                print 100, "Constrained Least Squares Solver Failed: Test 4b-", i
                print 101, "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print 101, "Solution:", x(1), ", ", x(2)
                print 101, "Residual:", f(1), ", ", f(2)
                print 102, "Converged on residual: ", ib%converge_on_fcn
                print 102, "Converged on solution change: ", &
                    ib%converge_on_chng
                print 102, "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print 100, "Iterations: ", ib%iter_count
                print 100, "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Formatting
100     format(A, I0)
101     format(A, F9.5, A, F9.5)
102     format(A, L)
    end function

! ------------------------------------------------------------------------------
    function test_constrained_least_squares_bounds() result(check)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        type(constrained_least_squares_solver) :: solver
        real(real64) :: x(2), f(2), low(2), high(2)
        logical :: check

        ! Initialization
        check = .true.
        fcn => fcn1
        call obj%set_fcn(fcn, 2, 2)

        ! Set tight, active bounds around the known solution (5,3)
        low = [4.0d0, 2.0d0]
        high = [5.6d0, 3.6d0]
        call solver%set_lower_limits(low)
        call solver%set_upper_limits(high)

        ! Starting point outside bounds
        x = [1.0d0, 1.0d0]

        ! Solve
        call solver%solve(obj, x, f)

        ! Bound enforcement should hold
        if (any(x < low - 1.0d-10) .or. any(x > high + 1.0d-10)) then
            check = .false.
            print 100, "Constrained Least Squares Bounds Solver Failed"
            print 101, "Result:", x(1), x(2)
            print 101, "Expected range:", low(1), low(2), high(1), high(2)
        end if

        ! The solution should be in the feasible region
        if (x(1) < low(1) .or. x(1) > high(1) .or. x(2) < low(2) .or. x(2) > high(2)) then
            check = .false.
        end if

        ! Formatting
100     format(A)
101     format(A, 4F10.6)
    end function

! ------------------------------------------------------------------------------
    subroutine vector_ad_system(x, f, args)
        type(dual), intent(in), dimension(:) :: x
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        real(real64) :: coefficient
        coefficient = 2.0d0
        if (present(args)) then
            select type (args)
            type is (real(real64))
                coefficient = args
            end select
        end if
        vector_ad_calls = vector_ad_calls + 1
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - coefficient * x(2)**2 - 7.0d0
    end subroutine

    subroutine vector_user_jacobian(x, jac, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: jac
        class(*), intent(in), optional :: args
        real(real64) :: coefficient
        coefficient = 2.0d0
        if (present(args)) then
            select type (args)
            type is (real(real64))
                coefficient = args
            end select
        end if
        vector_jac_calls = vector_jac_calls + 1
        jac = 2.0d0 * reshape([x(1), x(1), x(2), -coefficient * x(2)], [2,2])
    end subroutine

    subroutine make_vector_solver(which, solver)
        integer, intent(in) :: which
        class(equation_solver), allocatable, intent(out) :: solver
        select case (which)
        case (1, 2)
            allocate(newton_solver :: solver)
        case (3, 4)
            allocate(quasi_newton_solver :: solver)
        case (5)
            allocate(least_squares_solver :: solver)
        case (6)
            allocate(constrained_least_squares_solver :: solver)
        end select
        select type (solver)
        class is (line_search_solver)
            call solver%set_use_line_search(which == 1 .or. which == 3)
        end select
        call solver%set_max_fcn_evals(5000)
    end subroutine

    subroutine vector_ad_nonfinite(x, f, args)
        type(dual), intent(in), dimension(:) :: x
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        integer :: mode
        mode = 1
        if (present(args)) then
            select type (args)
            type is (integer)
                mode = args
            end select
        end if
        vector_ad_calls = vector_ad_calls + 1
        f%x = 1.0d0
        f%dx = 0.0d0
        select case (mode)
        case (1)
            f%x = ieee_value(0.0d0, ieee_quiet_nan)
        case (2)
            f%x = ieee_value(0.0d0, ieee_positive_inf)
        case (3)
            if (any(x%dx /= 0.0d0)) f%dx = ieee_value(0.0d0, ieee_quiet_nan)
        case (4)
            if (any(x%dx /= 0.0d0)) f%dx = ieee_value(0.0d0, ieee_positive_inf)
        end select
    end subroutine

    subroutine vector_ad_fit(x, f, args)
        type(dual), intent(in), dimension(:) :: x
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        real(real64), parameter, dimension(3) :: observations = [1.0d0, 2.0d0, 4.0d0]
        integer :: row
        do row = 1, size(f)
            f(row) = real(row - 1, real64) * x(1) + x(2) - observations(row)
        end do
    end subroutine

    subroutine vector_ad_rank_deficient(x, f, args)
        type(dual), intent(in), dimension(:) :: x
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        f(1) = x(1) + x(2) - 3.0d0
        f(2) = 2.0d0 * f(1)
    end subroutine

    subroutine vector_ad_scaled(x, f, args)
        type(dual), intent(in), dimension(:) :: x
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        f(1) = x(2) - 10.0d0
        f(2) = x(1) * x(2) - 5.0d4
    end subroutine

    subroutine vector_ad_single(x, f, args)
        type(dual), intent(in), dimension(:) :: x
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        f(1) = x(1)**2 - 4.0d0
    end subroutine

    function test_vector_solver_robustness() result(rst)
        class(equation_solver), allocatable :: solver
        type(vecfcn_helper) :: obj
        type(iteration_behavior) :: ib
        procedure(vector_function), pointer :: auto_fcn
        real(real64), dimension(2) :: x, values, low, high
        real(real64), dimension(1) :: single_x, single_f
        real(real64), dimension(3) :: fit_values, expected_fit
        integer :: which, mode
        logical :: rst, good

        rst = .true.
        do which = 1, 6
            call make_vector_solver(which, solver)
            call solver%set_max_fcn_evals(2)
            auto_fcn => vector_ad_system
            call obj%set_auto_fcn(auto_fcn, 2, 2)
            x = [0.5d0, 1.0d0]
            call solver%solve(obj, x, values, ib)
            good = all(ieee_is_nan(values)) .and. .not.ib%converge_on_fcn &
                .and. .not.ib%converge_on_chng .and. .not.ib%converge_on_zero_diff
            good = good .and. ib%fcn_count >= 2 .and. ib%fcn_count <= 30
            if (.not.good) print *, "AD evaluation-budget test failed:", which
            rst = good .and. rst

            auto_fcn => vector_ad_nonfinite
            call obj%set_auto_fcn(auto_fcn, 2, 2)
            do mode = 1, 4
                vector_ad_calls = 0
                x = 1.0d0
                call solver%solve(obj, x, values, ib, args = mode)
                good = all(ieee_is_nan(values)) .and. .not.ib%converge_on_fcn &
                    .and. .not.ib%converge_on_chng .and. .not.ib%converge_on_zero_diff
                good = good .and. ib%fcn_count == 1
                if (mode <= 2) then
                    good = good .and. ib%jacobian_count == 0 .and. vector_ad_calls == 1
                else
                    good = good .and. ib%jacobian_count == 1 .and. vector_ad_calls == 3
                end if
                if (.not.good) print *, "AD nonfinite-value test failed:", which, mode
                rst = good .and. rst
            end do

            x = 1.0d0
            call solver%solve(obj, x, values, ib, args = 5)
            if (which <= 4) then
                good = all(ieee_is_nan(values)) .and. .not.ib%converge_on_fcn &
                    .and. .not.ib%converge_on_chng
                good = good .and. ib%fcn_count <= 30 .and. ib%jacobian_count <= 2
            else
                good = all(ieee_is_finite(values)) .and. maxval(abs(values - 1.0d0)) < 1.0d-12
                good = good .and. (ib%converge_on_zero_diff .or. ib%converge_on_chng)
            end if
            if (.not.good) print *, "AD zero-Jacobian test failed:", which
            rst = good .and. rst

            call solver%set_max_fcn_evals(5000)
            auto_fcn => vector_ad_scaled
            call obj%set_auto_fcn(auto_fcn, 2, 2)
            x = 0.5d0
            call solver%solve(obj, x, values, ib)
            good = all(ieee_is_finite(x)) .and. all(ieee_is_finite(values))
            good = good .and. maxval(abs(x - [5.0d3, 10.0d0])) < 1.0d-5
            good = good .and. maxval(abs(values)) < 1.0d-5
            if (.not.good) print *, "AD poorly-scaled test failed:", which, x, values
            rst = good .and. rst

            auto_fcn => vector_ad_single
            call obj%set_auto_fcn(auto_fcn, 1, 1)
            single_x = 0.5d0
            call solver%solve(obj, single_x, single_f, ib)
            good = all(ieee_is_finite(single_x)) .and. all(ieee_is_finite(single_f))
            good = good .and. abs(single_x(1) - 2.0d0) < 1.0d-6 .and. abs(single_f(1)) < 1.0d-6
            if (.not.good) print *, "AD one-variable test failed:", which, single_x, single_f
            rst = good .and. rst

            if (which <= 4) then
                auto_fcn => vector_ad_rank_deficient
                call obj%set_auto_fcn(auto_fcn, 2, 2)
                call solver%set_max_fcn_evals(5)
                x = 0.0d0
                call solver%solve(obj, x, values, ib)
                good = all(ieee_is_nan(values)) .or. &
                    (all(ieee_is_finite(values)) .and. maxval(abs(values)) < 1.0d-6)
                if (.not.good) print *, "AD singular-system termination test failed:", which
                rst = good .and. rst
            end if
        end do

        auto_fcn => vector_ad_fit
        call obj%set_auto_fcn(auto_fcn, 3, 2)
        expected_fit = [-1.0d0 / 6.0d0, 1.0d0 / 3.0d0, -1.0d0 / 6.0d0]
        do which = 5, 6
            call make_vector_solver(which, solver)
            x = 0.0d0
            call solver%solve(obj, x, fit_values, ib)
            good = all(ieee_is_finite(x)) .and. all(ieee_is_finite(fit_values))
            good = good .and. maxval(abs(x - [1.5d0, 5.0d0 / 6.0d0])) < 1.0d-6
            good = good .and. maxval(abs(fit_values - expected_fit)) < 1.0d-6
            if (.not.good) print *, "AD rectangular least-squares test failed:", which, x, fit_values
            rst = good .and. rst
        end do

        call make_vector_solver(5, solver)
        auto_fcn => vector_ad_rank_deficient
        call obj%set_auto_fcn(auto_fcn, 2, 2)
        x = 0.0d0
        call solver%solve(obj, x, values, ib)
        good = all(ieee_is_finite(x)) .and. all(ieee_is_finite(values))
        good = good .and. abs(sum(x) - 3.0d0) < 1.0d-6
        good = good .and. maxval(abs(values)) < 1.0d-6
        if (.not.good) print *, "AD rank-deficient least-squares test failed:", x, values
        rst = good .and. rst

        call make_vector_solver(6, solver)
        auto_fcn => vector_ad_fit
        call obj%set_auto_fcn(auto_fcn, 3, 2)
        low = [-2.0d0, -2.0d0]
        high = [1.0d0, 3.0d0]
        select type (solver)
        type is (constrained_least_squares_solver)
            call solver%set_lower_limits(low)
            call solver%set_upper_limits(high)
        end select
        x = [0.0d0, 0.0d0]
        call solver%solve(obj, x, fit_values, ib)
        good = all(ieee_is_finite(x)) .and. all(ieee_is_finite(fit_values))
        good = good .and. all(x >= low) .and. all(x <= high)
        good = good .and. maxval(abs(x - [1.0d0, 4.0d0 / 3.0d0])) < 1.0d-5
        if (.not.good) print *, "AD active-bound least-squares test failed:", x, fit_values
        rst = good .and. rst

        low(1) = 1.0d0
        select type (solver)
        type is (constrained_least_squares_solver)
            call solver%set_lower_limits(low)
        end select
        x = [10.0d0, -10.0d0]
        call solver%solve(obj, x, fit_values, ib)
        good = all(ieee_is_finite(x)) .and. all(ieee_is_finite(fit_values))
        good = good .and. all(x >= low) .and. all(x <= high)
        good = good .and. maxval(abs(x - [1.0d0, 4.0d0 / 3.0d0])) < 1.0d-5
        if (.not.good) print *, "AD fixed-bound least-squares test failed:", x, fit_values
        rst = good .and. rst
    end function

    function test_vector_solvers_ad() result(rst)
        class(equation_solver), allocatable :: solver
        type(vecfcn_helper) :: obj
        type(iteration_behavior) :: ib
        procedure(vector_function), pointer :: auto_fcn
        procedure(jacobianfcn), pointer :: user_jac
        real(real64), dimension(2) :: x, values, exact, root
        real(real64), dimension(2,3) :: starts
        integer :: which, start, jac_mode
        logical :: rst, good

        rst = .true.
        auto_fcn => vector_ad_system
        starts(:,1) = [0.5d0, 1.0d0]
        starts(:,2) = [8.0d0, 6.0d0]
        starts(:,3) = [-4.0d0, -2.0d0]
        root = [sqrt(20.5d0), sqrt(13.5d0)]
        call obj%set_auto_fcn(auto_fcn, 2, 2)
        do which = 1, 6
            call make_vector_solver(which, solver)
            do jac_mode = 1, 2
                nullify(user_jac)
                if (jac_mode == 2) user_jac => vector_user_jacobian
                call obj%set_jacobian(user_jac)
                do start = 1, size(starts, 2)
                    x = starts(:,start)
                    vector_ad_calls = 0
                    vector_jac_calls = 0
                    call solver%solve(obj, x, values, ib, args = 1.0d0)
                    call fcn1a(x, exact, args = 1.0d0)
                    good = all(ieee_is_finite(x)) .and. all(ieee_is_finite(values))
                    good = good .and. maxval(abs(abs(x) - root)) < 1.0d-5
                    good = good .and. maxval(abs(exact)) < 1.0d-5
                    good = good .and. maxval(abs(values - exact)) < 1.0d-10
                    good = good .and. (ib%converge_on_fcn .or. ib%converge_on_chng &
                        .or. ib%converge_on_zero_diff)
                    good = good .and. ib%jacobian_count > 0 .and. ib%fcn_count > 1
                    if (jac_mode == 1) then
                        good = good .and. vector_jac_calls == 0
                        good = good .and. vector_ad_calls == ib%fcn_count + 2 * ib%jacobian_count
                    else
                        good = good .and. vector_jac_calls == ib%jacobian_count
                        good = good .and. vector_ad_calls == ib%fcn_count
                    end if
                    if (.not.good) then
                        print *, "AD solver convergence failed:", which, jac_mode, start
                        print *, "Solution/residual:", x, values
                        print *, "Counts:", vector_ad_calls, vector_jac_calls, &
                            ib%fcn_count, ib%jacobian_count
                    end if
                    rst = good .and. rst
                end do
            end do
            nullify(user_jac)
            call obj%set_jacobian(user_jac)
            x = [5.0d0, 3.0d0]
            vector_ad_calls = 0
            call solver%solve(obj, x, values, ib)
            good = all(ieee_is_finite(values)) .and. maxval(abs(values)) < 1.0d-12
            good = good .and. maxval(abs(x - [5.0d0, 3.0d0])) < 1.0d-12
            good = good .and. vector_ad_calls == ib%fcn_count + 2 * ib%jacobian_count
            if (which <= 4) good = good .and. ib%jacobian_count == 0 .and. ib%converge_on_fcn
            if (.not.good) print *, "AD exact-root test failed:", which
            rst = good .and. rst
        end do
    end function
end module
