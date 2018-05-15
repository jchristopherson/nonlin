! powell_badly_scaled.f90
! REF: http://people.sc.fsu.edu/~jburkardt/f_src/test_nonlin/test_nonlin.f90
! Problem #3, Powell's badly scaled function

! Computes the function
module powell_badly_scaled_module
    use iso_fortran_env
contains
subroutine powell_badly_scaled(x, f)
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: f
    f(1) = 1.0d4 * x(1) * x(2) - 1.0d0
    f(2) = exp(-x(1)) + exp(-x(2)) - 1.0001d0
end subroutine

! Computes the Jacobian matrix
subroutine powell_badly_scaled_jacobian(x, j)
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:,:) :: j
    j(1,1) = 1.0d4 * x(2)
    j(2,1) = -exp(-x(1))
    j(1,2) = 1.0d4 * x(1)
    j(2,2) = -exp(-x(2))
end subroutine

! Returns suitable initial conditions
function powell_badly_scaled_start() result(ic)
    real(real64) :: ic(2)
    ic = [0.0d0, 1.0d0]
end function

! Returns the solution
function powell_badly_scaled_solution() result(x)
    real(real64) :: x(2)
    x = [1.098159d-5, 9.106146d0]
end function
end module
