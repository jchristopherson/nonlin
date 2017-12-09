! nonlin_test_poly.f90

module nonlin_test_poly
    use nonlin_types, only : dp, i32
    use nonlin_polynomials
    use test_core
    implicit none
contains
! ******************************************************************************
! TEST ROUTINES
! ------------------------------------------------------------------------------
    ! Test fit a polynomial to the same data set that we used for testing the
    ! Levenberg-Marquardt curve fitting.
    subroutine test_poly_fit()
        ! Local Variables
        ! integer(i32) :: i
        real(dp), dimension(21) :: xp, yp, yf, yc, err
        real(dp) :: res
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
        ! print '(A)', "Polynomial Coefficients (c0 + c1*x + c2*x**2 + c3*x**3):"
        ! do i = 1, 4
        !     print '(AI0AF12.9)', "c", i - 1, " = ", p%get(i)
        ! end do
        ! print '(AE9.4)', "Residual: ", res
    end subroutine

! ------------------------------------------------------------------------------
    ! Tests the polynomial root finding capabilities.
    function test_poly_roots() result(check)
        ! Parameters
        real(dp), parameter :: tol = 1.0d-8
        integer(i32), parameter :: order = 10

        ! Local Variables
        integer(i32) :: i
        type(polynomial) :: p
        real(dp), dimension(order+1) :: coeff
        complex(dp), allocatable, dimension(:) :: rts, sol
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
                print '(A)', "Test Failed: Polynomial Roots Test 1"
                exit
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    ! Tests the polynomial addition routine
    function test_poly_add() result(check)
        ! Parameters
        integer(i32), parameter :: order1 = 10
        integer(i32), parameter :: order2 = 20
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(i32) :: i
        type(polynomial) :: p1, p2, p3
        real(dp), dimension(order1+1) :: c1
        real(dp), dimension(order2+1) :: c2
        real(dp), allocatable, dimension(:) :: p
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
            do i = 1, size(c2)
                c1(i) = c1(i) + c2(i)
            end do
            if (.not.is_mtx_equal(p, c1, tol)) then
                check = .false.
                print '(A)', "Test Failed: Polynomial Addition Test 1"
            end if
        else
            ! Use C2 to store the solution
            do i = 1, size(c1)
                c2(i) = c2(i) + c1(i)
            end do
            if (.not.is_mtx_equal(p, c2, tol)) then
                check = .false.
                print '(A)', "Test Failed: Polynomial Addition Test 1"
            end if
        end if
    end function

! ------------------------------------------------------------------------------
    ! Tests the polynomial subtraction routine
    function test_poly_subtract() result(check)
        ! Parameters
        integer(i32), parameter :: order1 = 10
        integer(i32), parameter :: order2 = 20
        real(dp), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(i32) :: i
        type(polynomial) :: p1, p2, p3
        real(dp), dimension(order1+1) :: c1
        real(dp), dimension(order2+1) :: c2
        real(dp), allocatable, dimension(:) :: p
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
            do i = 1, size(c2)
                c1(i) = c1(i) - c2(i)
            end do
            if (.not.is_mtx_equal(p, c1, tol)) then
                check = .false.
                print '(A)', "Test Failed: Polynomial Subtraction Test 1"
            end if
        else
            ! Use C2 to store the solution
            do i = 1, size(c1)
                c2(i) = c1(i) - c2(i)
            end do
            do i = size(c1) + 1, size(c2)
                c2(i) = -c2(i)
            end do
            if (.not.is_mtx_equal(p, c2, tol)) then
                check = .false.
                print '(A)', "Test Failed: Polynomial Subtraction Test 1"
            end if
        end if
    end function

! ------------------------------------------------------------------------------
    ! Tests the polynomial multiplication routine
    function test_poly_multiply() result(check)
        ! Local Variables
        logical :: check
        type(polynomial) :: p1, p2, p3, ans
        real(dp), allocatable, dimension(:) :: a, b
        real(dp), parameter :: tol = 1.0d-8

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
        if (.not.is_mtx_equal(a, b, tol)) then
            check = .false.
            print '(A)', "Test Failed: Polynomial Multiplication"
            print '(A)', "Expected:"
            print *, b
            print '(A)', "Computed:"
            print *, a
        else
            check = .true.
        end if
    end function

! ------------------------------------------------------------------------------
end module
