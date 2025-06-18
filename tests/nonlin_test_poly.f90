! nonlin_test_poly.f90

module nonlin_test_poly
    use iso_fortran_env
    use nonlin
    use fortran_test_helper
    implicit none
contains
! ******************************************************************************
! TEST ROUTINES
! ------------------------------------------------------------------------------
    ! Test fit a polynomial to the same data set that we used for testing the
    ! Levenberg-Marquardt curve fitting.
    subroutine test_poly_fit()
        ! Local Variables
        ! integer(int32) :: i
        real(real64), dimension(21) :: xp, yp, yf, yc, err
        real(real64) :: res
        type(polynomial) :: p

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

        ! Create a copy of yp as it will be overwritten in the fit command
        yc = yp

        ! Fit the polynomial
        call p%fit(xp, yp, 3)

        ! Evaluate the polynomial at xp, and then determine the residual
        yf = p%evaluate(xp)
        err = abs(yf - yc)
        res = maxval(err)

        ! Print out the coefficients
        ! print 100, "Polynomial Coefficients (c0 + c1*x + c2*x**2 + c3*x**3):"
        ! do i = 1, 4
        !     print '(AI0AF12.9)', "c", i - 1, " = ", p%get(i)
        ! end do
        ! print '(AE9.4)', "Residual: ", res
    end subroutine

! ------------------------------------------------------------------------------
    ! Tests the polynomial root finding capabilities.
    function test_poly_roots() result(check)
        ! Parameters
        real(real64), parameter :: tol = 1.0d-6
        integer(int32), parameter :: order = 10

        ! Local Variables
        integer(int32) :: i
        type(polynomial) :: p
        real(real64), dimension(order+1) :: coeff
        complex(real64), allocatable, dimension(:) :: rts, sol
        logical :: check

        ! Define the polynomial
        call random_number(coeff)
        call p%initialize(order)
        do i = 1, size(coeff)
            call p%set(i, coeff(i))
        end do

        ! Compute the roots via the polynomial routine
        rts = p%roots()

        ! Compute the value of the polynomial at each root and ensure it
        ! is sufficiently close to zero.
        sol = p%evaluate(rts)
        check = .true.
        do i = 1, size(sol)
            if (abs(sol(i)) > tol) then
                check = .false.
                print 100, "Test Failed: Polynomial Roots Test 1"
                exit
            end if
        end do

        ! Formatting
100     format(A)
    end function

! ------------------------------------------------------------------------------
    ! Tests the polynomial addition routine
    function test_poly_add() result(check)
        ! Parameters
        integer(int32), parameter :: order1 = 10
        integer(int32), parameter :: order2 = 20
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(int32) :: i
        type(polynomial) :: p1, p2, p3
        real(real64), dimension(order1+1) :: c1
        real(real64), dimension(order2+1) :: c2
        real(real64), allocatable, dimension(:) :: p
        logical :: check

        ! Define the polynomials
        call random_number(c1)
        call random_number(c2)
        call p1%initialize(order1)
        call p2%initialize(order2)
        do i = 1, size(c1)
            call p1%set(i, c1(i))
        end do
        do i = 1, size(c2)
            call p2%set(i, c2(i))
        end do

        ! Add the two polynomials
        p3 = p1 + p2
        p = p3%get_all()

        ! Compute the actual solution, and compare
        check = .true.
        if (size(c1) > size(c2)) then
            ! Use C1 to store the solution
            do i = 1, min(size(c1), size(c2))
                c1(i) = c1(i) + c2(i)
            end do
            if (.not.assert(p, c1, tol)) then
                check = .false.
                print 100, "Test Failed: Polynomial Addition Test 1"
            end if
        else
            ! Use C2 to store the solution
            do i = 1, min(size(c1), size(c2))
                c2(i) = c2(i) + c1(i)
            end do
            if (.not.assert(p, c2, tol)) then
                check = .false.
                print 100, "Test Failed: Polynomial Addition Test 1"
            end if
        end if

        ! Formatting
100     format(A)
    end function

! ------------------------------------------------------------------------------
    ! Tests the polynomial subtraction routine
    function test_poly_subtract() result(check)
        ! Parameters
        integer(int32), parameter :: order1 = 10
        integer(int32), parameter :: order2 = 20
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(int32) :: i
        type(polynomial) :: p1, p2, p3
        real(real64), dimension(order1+1) :: c1
        real(real64), dimension(order2+1) :: c2
        real(real64), allocatable, dimension(:) :: p
        logical :: check

        ! Define the polynomials
        call random_number(c1)
        call random_number(c2)
        call p1%initialize(order1)
        call p2%initialize(order2)
        do i = 1, size(c1)
            call p1%set(i, c1(i))
        end do
        do i = 1, size(c2)
            call p2%set(i, c2(i))
        end do

        ! Subtract the two polynomials
        p3 = p1 - p2
        p = p3%get_all()

        ! Compute the actual solution, and compare
        check = .true.
        if (size(c1) > size(c2)) then
            ! Use C1 to store the solution
            do i = 1, min(size(c1), size(c2))
                c1(i) = c1(i) - c2(i)
            end do
            if (.not.assert(p, c1, tol)) then
                check = .false.
                print 100, "Test Failed: Polynomial Subtraction Test 1"
            end if
        else
            ! Use C2 to store the solution
            do i = 1, min(size(c1), size(c2))
                c2(i) = c1(i) - c2(i)
            end do
            do i = size(c1) + 1, size(c2)
                c2(i) = -c2(i)
            end do
            if (.not.assert(p, c2, tol)) then
                check = .false.
                print 100, "Test Failed: Polynomial Subtraction Test 1"
            end if
        end if

        ! Formatting
100     format(A)
    end function

! ------------------------------------------------------------------------------
    ! Tests the polynomial multiplication routine
    function test_poly_multiply() result(check)
        ! Local Variables
        logical :: check
        type(polynomial) :: p1, p2, p3, ans
        real(real64), allocatable, dimension(:) :: a, b
        real(real64), parameter :: tol = 1.0d-8

        ! Initialization
        call p1%initialize(3)
        call p2%initialize(2)
        call ans%initialize(5)

        ! Set p1 = 5 + 10x**2 + 6x**3
        call p1%set(1, 5.0d0)
        call p1%set(2, 0.0d0)
        call p1%set(3, 10.0d0)
        call p1%set(4, 6.0d0)

        ! Set p2 = 1 + 2x + 4x**2
        call p2%set(1, 1.0d0)
        call p2%set(2, 2.0d0)
        call p2%set(3, 4.0d0)

        ! Answer: ans = 5 + 10x + 30x**2 + 26x**3 + 52x**4 + 24x**5
        call ans%set(1, 5.0d0)
        call ans%set(2, 10.0d0)
        call ans%set(3, 30.0d0)
        call ans%set(4, 26.0d0)
        call ans%set(5, 52.0d0)
        call ans%set(6, 24.0d0)

        ! Compute p1 * p2 = p3
        p3 = p1 * p2

        ! Test
        a = p3%get_all()
        b = ans%get_all()
        if (.not.assert(a, b, tol)) then
            check = .false.
            print 100, "Test Failed: Polynomial Multiplication"
            print 100, "Expected:"
            print *, b
            print 100, "Computed:"
            print *, a
        else
            check = .true.
        end if

        ! Formatting
100     format(A)
    end function

! ------------------------------------------------------------------------------
end module
