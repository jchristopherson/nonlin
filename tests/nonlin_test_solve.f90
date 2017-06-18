! nonlin_test_solve.f90

module nonlin_test_solve
    use linalg_constants, only : dp, i32
    use nonlin_types
    use nonlin_solve
    implicit none
    private
    public :: test_quasinewton_1
contains
! ------------------------------------------------------------------------------
    ! System of Equations #1:
    !
    ! x**2 + y**2 = 34
    ! x**2 - 2 * y**2 = 7
    !
    ! Solution:
    ! x = +/-5
    ! y = +/-3
    subroutine fcn1(x, f)
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: f
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    end subroutine

    ! Jacobian:
    !
    !     | 2x  2y |
    ! J = |        |
    !     | 2x  -4y|
    subroutine jac1(x, j)
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(out), dimension(:,:) :: j
        j = 2.0d0 * reshape([x(1), x(1), x(2), -2.0d0 * x(2)], [2, 2])
    end subroutine

    pure function is_ans_1(x, tol) result(c)
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(in) :: tol
        logical :: c
        real(dp), parameter :: x1 = 5.0d0
        real(dp), parameter :: x2 = 3.0d0
        real(dp) :: ax1, ax2
        c = .true.
        ax1 = abs(x(1)) - x1
        ax2 = abs(x(2)) - x2
        if (abs(ax1) > tol .or. abs(ax2) > tol) c = .false.
    end function

! ------------------------------------------------------------------------------
    subroutine test_quasinewton_1()
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(quasi_newton_solver) :: solver
        type(iteration_behavior) :: ib
        real(dp) :: x(2), f(2), ic(10, 2)
        integer(i32) :: i
        logical :: check

        ! Initialization
        check = .true.
        fcn => fcn1
        jac => jac1
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)

        !call solver%set_use_line_search(.false.)
        call solver%set_print_status(.true.)

        ! Generate a set of random initial conditions
        call random_number(ic)

        ! Process - Cycle over each different initial condition set
        do i = 1, size(ic, 1)
            x = ic(i,:)
            call solver%solve(obj, x, f, ib)
            if (.not.is_ans_1(x, 1.0d-5)) then
                check = .false.
                print '(AI0)', "Quasi-Newton Solver Failed: Test 1-", i
                print '(AF9.5AF9.5)', "Initial Condition: ", ic(i,1), ", ", &
                    ic(i,2)
                print '(AF9.5AF9.5)', "Solution:", x(1), ", ", x(2)
                print '(AF9.5AF9.5)', "Residual:", f(1), ", ", f(2)
                print '(AL)', "Converged on residual: ", ib%converge_on_fcn
                print '(AL)', "Converged on solution change: ", &
                    ib%converge_on_chng
                print '(AL)', "Converge on zero gradient: ", &
                    ib%converge_on_zero_diff
                print '(AI0)', "Iterations: ", ib%iter_count
                print '(AI0)', "Function Evaluations: ", ib%fcn_count
            end if
        end do

        ! Inform user of a succussful test
        if (check) then
            print '(A)', "Test Passed: Quasi-Newton Test 1"
        end if
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
