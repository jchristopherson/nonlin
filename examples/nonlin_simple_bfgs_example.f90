program example
    use iso_fortran_env
    use autodiff
    use nonlin
    implicit none

    real(real64) :: x(2), xo(2)
    type(dual) :: f, xv(2)

    ! Define an initial guess at the solution
    call random_number(xo)

    ! Compute the solution via BFGS
    x = solve_bfgs(rosenbrock, xo)

    ! Evaluate the function at x, and then display the results
    xv = x
    f = rosenbrock(xv)

    print "(A, F0.3, A, F0.3, A)", "x = [", x(1), ", ", x(2), "]"
    print "(A, E0.3)", "f(x) = ", f%x

contains
    function rosenbrock(x, args) result(f)
        type(dual), intent(in), dimension(:) :: x
        class(*), intent(in), optional :: args
        type(dual) :: f
        f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    end function
end program