! nonlin_test_solve.f90

module nonlin_test_solve
    use iso_fortran_env
    use nonlin
    use ferror, only : errors
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
        class(*), intent(inout), optional :: args
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    end subroutine

    subroutine fcn1a(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(inout), optional :: args
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
        class(*), intent(inout), optional :: args
        j = 2.0d0 * reshape([x(1), x(1), x(2), -2.0d0 * x(2)], [2, 2])
    end subroutine

    subroutine jac1a(x, j, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: j
        class(*), intent(inout), optional :: args
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
        class(*), intent(inout), optional :: args
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
        class(*), intent(inout), optional :: args

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
        class(*), intent(inout), optional :: args
        real(real64) :: f
        f = sin(x) / x
    end function

    function f1var_1a(x, args) result(f)
        real(real64), intent(in) :: x
        class(*), intent(inout), optional :: args
        real(real64) :: f
        real(real64) :: a
        select type (args)
        type is (real(real64))
            a = args
        end select
        f = a * sin(x) / x
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
        integer(int32) :: i
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
        integer(int32) :: i
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
        integer(int32) :: i
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
        integer(int32) :: i
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
        integer(int32) :: i
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
        integer(int32) :: i
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
        integer(int32) :: i
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
        integer(int32) :: i
        logical :: check
        type(errors) :: errmgr

        ! Initialization
        check = .true.
        fcn => fcn2
        call obj%set_fcn(fcn, 2, 2)

        ! Do not terminate testing if the solution does not converge.  This
        ! routine may have a bit of issue with this problem.  It can take
        ! many iterations to converge.
        call errmgr%set_exit_on_error(.false.)

        ! Increase the number of iterations allowed
        call solver%set_max_fcn_evals(1000)

        ! Generate a set of initial conditions
        ic(1,:) = 0.5d0
        ic(2,:) = 1.0d0

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib, err = errmgr)
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
        integer(int32) :: i
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
        type(value_pair) :: limits
        logical :: check

        ! Parameters
        real(real64), parameter :: pi = 3.141592653589793d0
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        a = 2.0d0
        fcn => f1var_1a
        call obj%set_fcn(fcn)

        ! Define the search limits
        limits%x1 = 1.5d0
        limits%x2 = 5.0d0

        ! Compute the solution
        call solver%solve(obj, x, limits, f, args = a)

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
        integer(int32) :: i
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
        integer(int32) :: i
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
        type(value_pair) :: limits
        logical :: check

        ! Parameters
        real(real64), parameter :: pi = 3.141592653589793d0
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        a = 2.0d0
        fcn => f1var_1a
        call obj%set_fcn(fcn)

        ! Define the search limits
        limits%x1 = 1.5d0
        limits%x2 = 5.0d0

        ! Compute the solution
        call solver%solve(obj, x, limits, f, args = a)

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
end module
