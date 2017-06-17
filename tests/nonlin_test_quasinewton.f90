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

! ------------------------------------------------------------------------------
    subroutine test_quasinewton_1()
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        type(quasi_newton_solver) :: solver
        type(iteration_behavior) :: ib
        real(dp) :: x(2), f(2)

        ! Initialization
        fcn => fcn1
        call obj%set_fcn(fcn, 2, 2)
        x = 1.0d0

        ! Process
        call solver%solve(obj, x, f, ib)
        print '(AF6.3AF6.3)', "Solution 1:", x(1), ", ", x(2)
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
