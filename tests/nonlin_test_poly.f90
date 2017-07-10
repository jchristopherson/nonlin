! nonlin_test_poly.f90

module nonlin_test_poly
    use linalg_constants, only : dp, i32
    use nonlin_polynomials
    implicit none
contains
! ******************************************************************************
! TEST ROUTINES
! ------------------------------------------------------------------------------
    ! Test fit a polynomial to the same data set that we used for testing the
    ! Levenberg-Marquardt curve fitting.
    subroutine test_poly_fit()
        ! Local Variables
        integer(i32) :: i
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
        print '(A)', "Polynomial Coefficients (c0 + c1*x + c2*x**2 + c3*x**3):"
        do i = 1, 4
            print '(AI0AF12.9)', "c", i - 1, " = ", p%get(i)
        end do
        print '(AE9.4)', "Residual: ", res
    end subroutine

! ------------------------------------------------------------------------------
    ! Tests the polynomial root finding capabilities.
    subroutine test_poly_roots()
        ! Parameters
        complex(dp), parameter :: j = complex(0.0d0, 1.0d0)

        ! Local Variables
        type(polynomial) :: poly
        real(dp) :: p, q, r, a, b
        complex(dp) :: y1, y2, y3, x1, x2, x3, aa, bb
        complex(dp), allocatable, dimension(:) :: rts

        ! Define the polynomial
        call random_number(p)
        call random_number(q)
        call random_number(r)
        call poly%initialize(3)
        call poly%set(1, r)
        call poly%set(2, q)
        call poly%set(3, p)
        call poly%set(4, 1.0d0)

        ! Roots of a cubic equation: x**3 + p*x**2 + q*x + r = 0 are given
        ! as follows:
        !
        ! a = (1/3) * (3*q - p**2)
        ! b = (1/27) * (2*p**3 - 9*p*q + 27*r)
        !
        ! A = (-b/2 + (b**2/4 + a**3/27)**(1/2))**(1/3)
        ! B = (-b/2 - (b**2/4 + a**3/27)**(1/2))**(1/3)
        !
        ! y1 = A + B
        ! y2 = -(1/2) * (A + B) + i * sqrt(3)/2 * (A - B)
        ! y3 = -(1/2) * (A + B) - i * sqrt(3)/2 * (A - B)
        !
        ! x = y - p / 3
        !
        ! REF: http://web.cs.iastate.edu/~cs577/handouts/polyroots.pdf

        ! Compute the solution
        a = (1.0d0 / 3.0d0) * (3.0d0 * q - p**2)
        b = (1.0d0 / 27.0d0) * (2.0d0 * p**3 - 9.0d0 * p * q + 27.0d0 * r)
        aa = (-b / 2.0d0 + &
            sqrt(complex(b**2 / 4.0d0 + a**3 / 27.0d0, 0.0d0)))**(1.0d0 / 3.0d0)
        bb = (-b / 2.0d0 - &
            sqrt(complex(b**2 / 4.0d0 + a**3 / 27.0d0, 0.0d0)))**(1.0d0 / 3.0d0)
        y1 = aa + bb
        y2 = -(1.0d0 / 2.0d0) * (aa + bb) + &
            j * (sqrt(3.0d0) / 2.0d0) * (aa - bb)
        y3 = -(1.0d0 / 2.0d0) * (aa + bb) - &
            j * (sqrt(3.0d0) / 2.0d0) * (aa - bb)
        x1 = y1 - p / 3.0d0
        x2 = y2 - p / 3.0d0
        x3 = y3 - p / 3.0d0

        ! Print the polynomial
        print *, p
        print *, q
        print *, r

        ! Compute the roots via the polynomial routine
        rts = poly%roots()

        ! Display the solution
        print '(A)', "ROOTS (SOLUTION):"
        print *, x1
        print *, x2
        print *, x3
        print *, ""
        print '(A)', "ROOTS COMPUTED:"
        print *, rts(1)
        print *, rts(2)
        print *, rts(3)
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
