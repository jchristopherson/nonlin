! nonlin_least_squares_example.f90

program example
    use iso_fortran_env
    use nonlin_core
    use nonlin_least_squares, only : least_squares_solver
    implicit none

    ! Local Variables
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    type(least_squares_solver) :: solver
    real(real64) :: x(4), f(21) ! There are 4 coefficients and 21 data points

    ! Locate the routine containing the equations to solve
    fcn => fcns
    call obj%set_fcn(fcn, 21, 4)

    ! Define an initial guess
    x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0, 1.0d0, 1.0d0]

    ! Solve
    call solver%solve(obj, x, f)

    ! Display the output
    print '(AF12.10)', "c0: ", x(4)
    print '(AF12.10)', "c1: ", x(3)
    print '(AF12.10)', "c2: ", x(2)
    print '(AF12.10)', "c3: ", x(1)
    print '(AF7.5)', "Max Residual: ", maxval(abs(f))

contains
    ! The function containing the data to fit
    subroutine fcns(x, f)
        ! Arguments
        real(real64), intent(in), dimension(:) :: x  ! Contains the coefficients
        real(real64), intent(out), dimension(:) :: f

        ! Local Variables
        real(real64), dimension(21) :: xp, yp

        ! Data to fit (21 data points)
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
        ! y = c3 * x**3 + c2 * x**2 + c1 * x + c0
        f = x(1) * xp**3 + x(2) * xp**2 + x(3) * xp + x(4) - yp

        ! For reference, the data was generated by adding random errors to
        ! the following polynomial: y = x**3 - 0.3 * x**2 + 1.2 * x + 0.3
    end subroutine
end program
