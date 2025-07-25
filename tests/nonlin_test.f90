! nonlin_test.f90

! The testing application for the NONLIN library.
program main
    ! Imported Modules
    use nonlin_test_jacobian
    use nonlin_test_solve
    use nonlin_test_poly
    use nonlin_test_optimize

    ! Local Variables
    logical :: rst, overall

    ! Initialization
    overall = .true.

    ! Tests
    rst = test_jacobian_1()
    if (.not.rst) overall = .false.

    rst = test_jacobian_2()
    if (.not.rst) overall = .false.

    rst = test_quasinewton_1()
    if (.not.rst) overall = .false.

    rst = test_quasinewton_2()
    if (.not.rst) overall = .false.

    rst = test_newton_1()
    if (.not.rst) overall = .false.

    rst = test_newton_2()
    if (.not.rst) overall = .false.

    rst = test_least_squares_1()
    if (.not.rst) overall = .false.

    rst = test_least_squares_2()
    if (.not.rst) overall = .false.

    call test_least_squares_3()

    rst = test_brent_1()
    if (.not.rst) overall = .false.

    call test_poly_fit()
    rst = test_poly_roots()
    if (.not.rst) overall = .false.

    rst = test_poly_add()
    if (.not.rst) overall = .false.

    rst = test_poly_subtract()
    if (.not.rst) overall = .false.

    rst = test_poly_multiply()
    if (.not.rst) overall = .false.

    rst = test_nelder_mead_1()
    if (.not.rst) overall = .false.

    rst = test_nelder_mead_2()
    if (.not.rst) overall = .false.

    rst = test_nelder_mead_3()
    if (.not.rst) overall = .false.

    rst = test_bfgs_1()
    if (.not.rst) overall = .false.

    rst = test_bfgs_2()
    if (.not.rst) overall = .false.

    rst = test_bfgs_3()
    if (.not.rst) overall = .false.

    rst = test_newton_3()
    if (.not.rst) overall = .false.

    rst = test_quasinewton_3()
    if (.not.rst) overall = .false.

    rst = test_quasinewton_4()
    if (.not.rst) overall = .false.

    rst = test_newton_4()
    if (.not.rst) overall = .false.

    rst = test_least_squares_4()
    if (.not.rst) overall = .false.

    rst = test_brent_2()
    if (.not.rst) overall = .false.

    rst = test_newton_1var_1()
    if (.not.rst) overall = .false.

    rst = test_newton_1var_2()
    if (.not.rst) overall = .false.


    ! End
    if (.not.overall) stop -1
end program
