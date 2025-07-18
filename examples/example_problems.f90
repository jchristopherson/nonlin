module example_problems
    use iso_fortran_env
contains
! ------------------------------------------------------------------------------
! f(x) = sin(x) / x, solution: x = x * pi for n = 1, 2, 3, ...
function sinx_div_x(x, args) result(f)
    real(real64), intent(in) :: x
    class(*), intent(inout), optional :: args
    real(real64) :: f
    f = sin(x) / x
end function

! ------------------------------------------------------------------------------
! The Beale function:
! f(x) = (1.5 - x + xy)**2 + (2.25 - x + xy**2)**2 + (2.625 - x + xy**3)**2
! The minimum is at x = 3, y = 0.5, and f(3, 0.5) = 0
function beale(x, args) result(f)
    real(real64), intent(in), dimension(:) :: x
    class(*), intent(inout), optional :: args
    real(real64) :: f
    f = (1.5d0 - x(1) + x(1) * x(2))**2 + &
        (2.25d0 - x(1) + x(1) * x(2)**2)**2 + &
        (2.625d0 - x(1) + x(1) * x(2)**3)**2
end function

! The gradient
subroutine bealegrad(x, g, args)
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: g
    class(*), intent(inout), optional :: args

    g(1) = 2.0d0 * (x(2)**3 - 1.0d0) * (x(1) * x(2)**3 - x(1) + 2.625d0) + &
        2.0d0 * (x(2)**2 - 1.0d0) * (x(1) * x(2)**2 - x(1) + 2.25d0) + &
        2.0d0 * (x(2) - 1.0d0) * (x(1) * x(2) - x(1) + 1.5d0)

    g(2) = 6.0d0 * x(1) * x(2)**2 * (x(1) * x(2)**3 - x(1) + 2.625d0) + &
        4.0d0 * x(1) * x(2) * (x(1) * x(2)**2 - x(1) + 2.25d0) + &
        2.0d0 * x(1) * (x(1) * x(2) - x(1) + 1.5d0)
end subroutine

! ------------------------------------------------------------------------------
! REF: http://people.sc.fsu.edu/~jburkardt/f_src/test_nonlin/test_nonlin.f90
! Problem #3, Powell's badly scaled function
!
! Solution:
! x(1) = 1.098159e-5
! x(2) = 9.106146
subroutine powell_bad(x, f, args)
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: f
    class(*), intent(inout), optional :: args
    f(1) = 1.0d4 * x(1) * x(2) - 1.0d0
    f(2) = exp(-x(1)) + exp(-x(2)) - 1.0001d0
end subroutine

! ------------------------------------------------------------------------------
! x**2 + y**2 = 34
! x**2 - 2 * y**2 = 7
subroutine misc_2fcn(x, f, args)
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: f
    class(*), intent(inout), optional :: args
    f(1) = x(1)**2 + x(2)**2 - 34.0d0
    f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
end subroutine

! ------------------------------------------------------------------------------
! The system of equations (source: https://www.mathworks.com/help/optim/ug/fsolve.html)
! 2 * x1 - x2 = exp(-x1)
! -x1 + 2 * x2 = exp(-x2)
subroutine misc_2fcn_01(x, f, args)
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: f
    class(*), intent(inout), optional :: args
    f(1) = 2.0d0 * x(1) - x(2) - exp(-x(1))
    f(2) = -x(1) + 2.0d0 * x(2) - exp(-x(2))
end subroutine

! The Jacobian matrix:
!     | exp(-x1) + 2          -1     |
! J = |                              |
!     |     -1          exp(-x2) + 2 |
subroutine misc_2fcn_01_jac(x, jac, args)
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:,:) :: jac
    class(*), intent(inout), optional :: args
    jac(1,1) = exp(-x(1)) + 2.0d0
    jac(2,1) = -1.0d0
    jac(1,2) = -1.0d0
    jac(2,2) = exp(-x(2)) + 2.0d0
end subroutine

! ------------------------------------------------------------------------------
! f(x) = x**3 - 2 * x - 1
function fcn_1var(x, args) result(f)
    real(real64), intent(in) :: x
    class(*), intent(inout), optional :: args
    real(real64) :: f
    f = x**3 - 2.0d0 * x - 1.0d0
end function

! ------------------------------------------------------------------------------
! Rosenbrock's Function
function rosenbrock(x, args) result(f)
    real(real64), intent(in), dimension(:) :: x
    class(*), intent(inout), optional :: args
    real(real64) :: f
    f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
end function

! ------------------------------------------------------------------------------
! Least-Squares Polynomial Fitting Function
subroutine lsq_poly_fit_fcn(x, f, args)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x  ! Contains the coefficients
    real(real64), intent(out), dimension(:) :: f
    class(*), intent(inout), optional :: args

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

! ------------------------------------------------------------------------------
end module