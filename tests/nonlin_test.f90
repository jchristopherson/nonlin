! nonlin_test.f90

! The testing application for the NONLIN library.
program main
    ! Imported Modules
    use test_jacobian

    ! Introduce the testing application
    print '(A)', "Hello from the NONLIN test application."

    ! Tests
    call test_jacobian_1()
end program
