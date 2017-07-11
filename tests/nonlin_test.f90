! nonlin_test.f90

! The testing application for the NONLIN library.
program main
    ! Imported Modules
    use nonlin_test_jacobian
    use nonlin_test_solve
    use nonlin_test_poly

    ! Introduce the testing application
    print '(A)', "Hello from the NONLIN test application."

    ! Tests
    call test_jacobian_1()
    call test_quasinewton_1()
    call test_quasinewton_2()
    call test_newton_1()
    call test_newton_2()
    call test_least_squares_1()
    call test_least_squares_2()
    call test_least_squares_3()
    call test_brent_1()
    call test_poly_fit()
    call test_poly_roots()
    call test_poly_add()
end program
