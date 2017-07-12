! nonlin_test_poly.f90

module nonlin_test_poly
    use linalg_constants, only : dp, i32
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
    subroutine test_poly_roots()
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

        if (check) print '(A)', "Test Passed: Polynomial Roots Test 1"
    end subroutine

! ------------------------------------------------------------------------------
    ! Tests the polynomial addition routine
    subroutine test_poly_add()
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

        if (check) print '(A)', "Test Passed: Polynomial Addition Test 1"
    end subroutine

! ------------------------------------------------------------------------------
    ! Tests the polynomial subtraction routine
    subroutine test_poly_subtract()
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

        if (check) print '(A)', "Test Passed: Polynomial Subtraction Test 1"
    end subroutine

! ------------------------------------------------------------------------------
end module