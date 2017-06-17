! nonlin_test_quasinewton.f90

module nonlin_test_quasinewton
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

! ------------------------------------------------------------------------------
    subroutine test_quasinewton_1()
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(quasi_newton_solver) :: solver
        type(iteration_behavior) :: ib
        real(dp) :: x(2), f(2)

        ! Initialization
        fcn => fcn1
        jac => jac1
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)
        x = [1.0d0, 1.0d0]

        ! Process
        call solver%solve(obj, x, f, ib)
        print '(AF7.3AF7.3)', "Solution 1:", x(1), ", ", x(2)
        print '(AF7.3AF7.3)', "Residual 1:", f(1), ", ", f(2)
        print '(AL)', "Converged on residual: ", ib%converge_on_fcn
        print '(AL)', "Converged on function change: ", ib%converge_on_chng
        print '(AL)', "Converge on zero gradient: ", ib%converge_on_zero_diff
        print '(AI0)', "Iterations: ", ib%iter_count
        print '(AI0)', "Function Evaluations: ", ib%fcn_count
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
