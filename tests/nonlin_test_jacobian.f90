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
! nonlin_test_jacobian.f90

! Tests the Jacobian calculation utilities.
module nonlin_test_jacobian
    use iso_fortran_env
    use nonlin
    use fortran_test_helper
    use autodiff, only : dual, vector_function, assignment(=), &
        operator(+), operator(-), operator(*), operator(**)
    implicit none
    private
    public :: test_jacobian_1
    public :: test_jacobian_2
    public :: test_vecfcn_ad
    integer :: ad_calls = 0, user_jac_calls = 0

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
        class(*), intent(in), optional :: args
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
        class(*), intent(in), optional :: args
        real(real64) :: r, theta
        r = x(1)
        theta = x(2)
        jac = reshape(&
            [cos(theta), sin(theta), &
            -r * sin(theta), r * cos(theta)], &
            [2, 2])
    end subroutine

! ------------------------------------------------------------------------------
    ! FUNCTIONS:
    ! x = Y * r * cos(theta)
    ! y = Y * r * sin(theta)
    subroutine fcn2(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        real(real64) :: r, theta, y
        r = x(1)
        theta = x(2)
        select type (args)
        type is (real(real64))
            y = args
        end select
        f(1) = y * r * cos(theta) 
        f(2) = y * r * sin(theta)
    end subroutine

    ! JACOBIAN:
    !     |Y * cos  -Y * r * sin|
    ! J = |                     |
    !     |Y * sin   Y * r * cos|
    subroutine jac2(x, jac, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: jac
        class(*), intent(in), optional :: args
        real(real64) :: r, theta, y
        r = x(1)
        theta = x(2)
        select type (args)
        type is (real(real64))
            y = args
        end select
        jac = reshape(&
            [y * cos(theta), y * sin(theta), &
            -y * r * sin(theta), y * r * cos(theta)], &
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
    function test_jacobian_2() result(rst)
        ! Local Variables
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        real(real64) :: numjac(2, 2), exact(2, 2), x(2), y
        logical :: rst
        integer(int32) :: i

        ! Parameters
        real(real64), parameter :: tol = 1.0d-4

        ! Initialization
        rst = .true.
        fcn => fcn2
        jac => jac2
        call obj%set_fcn(fcn, 2, 2)
        call random_number(y)

        ! Test point 1 (0, 0)
        x = 0.0d0
        call obj%jacobian(x, numjac, args = y)
        call jac(x, exact, args = y)
        if (.not.assert(numjac, exact, tol)) then
            rst = .false.
            print 100, "Test Failed: Jacobian Test 2A"
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
        call obj%jacobian(x, numjac, args = y)
        call jac(x, exact, args = y)
        if (.not.assert(numjac, exact, tol)) then
            rst = .false.
            print 100, "Test Failed: Jacobian Test 2B"
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
        call obj%jacobian(x, numjac, args = y)
        call jac(x, exact, args = y)
        if (.not.assert(numjac, exact, tol)) then
            rst = .false.
            print 100, "Test Failed: Jacobian Test 2C"
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
        call obj%jacobian(x, numjac, args = y)
        call jac(x, exact, args = y)
        if (.not.assert(numjac, exact, tol)) then
            rst = .false.
            print 100, "Test Failed: Jacobian Test 2D"
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
    subroutine ad_vector(x, f, args)
        type(dual), intent(in), dimension(:) :: x
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        real(real64) :: offset
        integer :: row
        offset = 0.0d0
        if (present(args)) then
            select type (args)
            type is (real(real64))
                offset = args
            end select
        end if
        ad_calls = ad_calls + 1
        do row = 1, size(f)
            f(row) = x(1)**2 + real(row, real64) * x(size(x)) - offset
        end do
    end subroutine

    subroutine counted_jacobian(x, jac, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: jac
        class(*), intent(in), optional :: args
        integer :: row
        user_jac_calls = user_jac_calls + 1
        jac = 0.0d0
        jac(:,1) = 2.0d0 * x(1)
        do row = 1, size(jac, 1)
            jac(row,size(x)) = jac(row,size(x)) + real(row, real64)
        end do
    end subroutine

    function test_vecfcn_ad() result(rst)
        type(vecfcn_helper) :: obj
        procedure(vector_function), pointer :: auto_fcn
        procedure(vecfcn), pointer :: real_fcn
        procedure(jacobianfcn), pointer :: user_jac
        real(real64), dimension(2) :: x, values, expected
        real(real64), dimension(2,2) :: actual_jac, exact_jac
        real(real64), dimension(1) :: small_x, small_f
        real(real64), dimension(1,1) :: small_jac
        real(real64), dimension(3) :: tall_f
        real(real64), dimension(3,2) :: tall_jac
        integer :: previous_calls
        logical :: rst

        rst = .true.
        auto_fcn => ad_vector
        real_fcn => fcn1
        user_jac => counted_jacobian
        call obj%set_auto_fcn(auto_fcn, 2, 2)
        rst = obj%is_fcn_defined() .and. .not.obj%is_jacobian_defined()
        x = [2.0d0, 3.0d0]
        expected = [7.0d0, 10.0d0]
        call obj%fcn(x, values)
        rst = assert(values, expected, 1.0d-12) .and. rst
        call obj%jacobian(x, actual_jac, args = 5.0d0)
        call counted_jacobian(x, exact_jac)
        rst = assert(actual_jac, exact_jac, 1.0d-12) .and. rst
        rst = assert(x, [2.0d0, 3.0d0], 1.0d-12) .and. rst
        call obj%fcn(x, values, args = 5.0d0)
        rst = assert(values, expected - 5.0d0, 1.0d-12) .and. rst

        call obj%set_fcn(real_fcn, 2, 2)
        call obj%fcn(x, values)
        rst = assert(values, expected, 1.0d-12) .and. rst
        call obj%set_jacobian(user_jac)
        previous_calls = ad_calls
        user_jac_calls = 0
        call obj%jacobian(x, actual_jac, fv = expected)
        rst = assert(actual_jac, exact_jac, 1.0d-12) .and. rst
        rst = ad_calls == previous_calls .and. user_jac_calls == 1 .and. rst
        nullify(user_jac)
        call obj%set_jacobian(user_jac)

        call obj%set_auto_fcn(auto_fcn, 1, 1)
        small_x = 2.0d0
        call obj%fcn(small_x, small_f)
        rst = assert(small_f, [6.0d0], 1.0d-12) .and. rst
        call obj%jacobian(small_x, small_jac)
        rst = abs(small_jac(1,1) - 5.0d0) < 1.0d-12 .and. rst
        call obj%set_auto_fcn(auto_fcn, 3, 2)
        call obj%fcn(x, tall_f)
        rst = assert(tall_f, [7.0d0, 10.0d0, 13.0d0], 1.0d-12) .and. rst
        call obj%jacobian(x, tall_jac)
        rst = assert(tall_jac(:,1), [4.0d0, 4.0d0, 4.0d0], 1.0d-12) .and. rst
        rst = assert(tall_jac(:,2), [1.0d0, 2.0d0, 3.0d0], 1.0d-12) .and. rst

        nullify(auto_fcn)
        call obj%set_auto_fcn(auto_fcn, 2, 2)
        call obj%fcn(x, values)
        call fcn1(x, expected)
        rst = assert(values, expected, 1.0d-12) .and. rst
        call obj%jacobian(x, actual_jac, fv = values)
        call jac1(x, exact_jac)
        rst = assert(actual_jac, exact_jac, 1.0d-4) .and. rst
        nullify(real_fcn)
        call obj%set_fcn(real_fcn, 2, 2)
        rst = .not.obj%is_fcn_defined() .and. rst
        if (.not.rst) print *, "Test Failed: vecfcn_helper AD handling"
    end function
end module
