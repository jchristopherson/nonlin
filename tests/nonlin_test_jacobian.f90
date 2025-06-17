! nonlin_test_jacobian.f90

! Tests the Jacobian calculation utilities.
module nonlin_test_jacobian
    use iso_fortran_env
    use nonlin
    use fortran_test_helper
    implicit none
    private
    public :: test_jacobian_1

contains
! ******************************************************************************
! FUNCTIONS:
! ------------------------------------------------------------------------------
    ! FUNCTIONS:
    ! x = r * cos(theta)
    ! y = r * sin(theta)
    subroutine fcn1(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(inout), optional :: args
        real(real64) :: r, theta
        r = x(1)
        theta = x(2)
        f(1) = r * cos(theta)
        f(2) = r * sin(theta)
    end subroutine

    ! JACOBIAN:
    !     |cos  -r*sin|
    ! J = |           |
    !     |sin   r*cos|
    subroutine jac1(x, jac, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: jac
        class(*), intent(inout), optional :: args
        real(real64) :: r, theta
        r = x(1)
        theta = x(2)
        jac = reshape(&
            [cos(theta), sin(theta), &
            -r * sin(theta), r * cos(theta)], &
            [2, 2])
    end subroutine

! ------------------------------------------------------------------------------
    function test_jacobian_1() result(rst)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        real(real64) :: numjac(2, 2), exact(2, 2), x(2)
        logical :: rst
        integer(int32) :: i

        ! Parameters
        real(real64), parameter :: tol = 1.0d-4

        ! Initialization
        rst = .true.
        fcn => fcn1
        jac => jac1
        call obj%set_fcn(fcn, 2, 2)

        ! Test point 1 (0, 0)
        x = 0.0d0
        call obj%jacobian(x, numjac)
        call jac(x, exact)
        if (.not.assert(numjac, exact, tol)) then
            rst = .false.
            print 100, "Test Failed: Jacobian Test 1A"
            print 100, "Expected:"
            do i = 1, size(exact, 1)
                print *, exact(i,:)
            end do
            print 100, "Found:"
            do i = 1, size(exact, 1)
                print *, numjac(i,:)
            end do
        end if

        ! Test point 2 (1, 0)
        x = [1.0d0, 0.0d0]
        call obj%jacobian(x, numjac)
        call jac(x, exact)
        if (.not.assert(numjac, exact, tol)) then
            rst = .false.
            print 100, "Test Failed: Jacobian Test 1B"
            print 100, "Expected:"
            do i = 1, size(exact, 1)
                print *, exact(i,:)
            end do
            print 100, "Found:"
            do i = 1, size(exact, 1)
                print *, numjac(i,:)
            end do
        end if

        ! Test point 3 (0, 1)
        x = [0.0d0, 1.0d0]
        call obj%jacobian(x, numjac)
        call jac(x, exact)
        if (.not.assert(numjac, exact, tol)) then
            rst = .false.
            print 100, "Test Failed: Jacobian Test 1C"
            print 100, "Expected:"
            do i = 1, size(exact, 1)
                print *, exact(i,:)
            end do
            print 100, "Found:"
            do i = 1, size(exact, 1)
                print *, numjac(i,:)
            end do
        end if

        ! Test point 4 (0.5, -0.5)
        x = [0.5d0, -0.5d0]
        call obj%jacobian(x, numjac)
        call jac(x, exact)
        if (.not.assert(numjac, exact, tol)) then
            rst = .false.
            print 100, "Test Failed: Jacobian Test 1D"
            print 100, "Expected:"
            do i = 1, size(exact, 1)
                print *, exact(i,:)
            end do
            print 100, "Found:"
            do i = 1, size(exact, 1)
                print *, numjac(i,:)
            end do
        end if

        ! Formatting
100     format(A)
    end function

! ------------------------------------------------------------------------------
end module
