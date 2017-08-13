! nonlin_poly_fit_example.f90

program example
    use nonlin_types, only : dp, i32
    use nonlin_polynomials
    implicit none

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
    print '(AI0AF12.10)', ("c", i - 1, " = ", p%get(i), i = 1, 4)
    print '(AF7.5)', "Max Residual: ", res
end program