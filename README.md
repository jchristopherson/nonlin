# nonlin
A library that provides routines to compute the solutions to systems of nonlinear equations.

## Status
[![CMake](https://github.com/jchristopherson/nonlin/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/nonlin/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/nonlin/workflows/fpm/badge.svg)](https://github.com/jchristopherson/nonlin/actions)

## Documentation
Documentation can be found [here](https://jchristopherson.github.io/nonlin/)

## Building NONLIN
[CMake](https://cmake.org/)This library can be built using CMake.  For instructions see [Running CMake](https://cmake.org/runningcmake/).

[FPM](https://github.com/fortran-lang/fpm) can also be used to build this library using the provided fpm.toml.
```txt
fpm build
```
The NONLIN library can be used within your FPM project by adding the following to your fpm.toml file.
```toml
[dependencies]
nonlin = { git = "https://github.com/jchristopherson/nonlin" }
```
## External Libraries
Here is a list of external code libraries utilized by this library.
- [BLAS](http://www.netlib.org/blas/)
- [LAPACK](http://www.netlib.org/lapack/)
- [FERROR](https://github.com/jchristopherson/ferror)
- [LINALG](https://github.com/jchristopherson/linalg)

## Example 1
This example solves a set of two equations of two unknowns using a Quasi-Newton type solver.  In this example, the solver is left to compute the derivatives numerically.

```fortran
program  example
    use iso_fortran_env
    use nonlin_core
    use nonlin_solve, only : quasi_newton_solver
    implicit none

    ! Local Variables
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    type(iteration_behavior) :: ib
    type(quasi_newton_solver) :: solver
    real(real64) :: x(2), f(2)

    ! Locate the routine containing the equations to solve
    fcn => fcns
    call obj%set_fcn(fcn, 2, 2)

    ! Define an initial guess
    x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]

    ! Defining solver parameters.  This step is optional as the defaults are
    ! typically sufficient; however, this is being done for illustration 
    ! purposes.
    !
    ! Establish how many iterations are allowed to pass before the solver
    ! forces a re-evaluation of the Jacobian matrix.  Notice, the solver may
    ! choose to re-evaluate the Jacobian sooner than this, but that is 
    ! dependent upon the behavior of the problem.
    call solver%set_jacobian_interval(20)

    ! Establish convergence criteria.  Again, this step is optional as the
    ! defaults are typically sufficient; however, this is being done for
    ! illustration purposes.
    call solver%set_fcn_tolerance(1.0d-8)
    call solver%set_var_tolerance(1.0d-12)
    call solver%set_gradient_tolerance(1.0d-12)

    ! Solve
    call solver%solve(obj, x, f, ib)

    ! Display the output
    print '(AF7.5AF7.5A)', "Solution: (", x(1), ", ", x(2), ")"
    print '(AE9.3AE9.3A)', "Residual: (", f(1), ", ", f(2), ")"
    print '(AI0)', "Iterations: ", ib%iter_count
    print '(AI0)', "Function Evaluations: ", ib%fcn_count
    print '(AI0)', "Jacobian Evaluations: ", ib%jacobian_count

contains
    ! Define the routine containing the equations to solve.  The equations are:
    ! x**2 + y**2 = 34
    ! x**2 - 2 * y**2 = 7
    subroutine fcns(x, f)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    end subroutine
end program
```
The above program produces the following output.
```text
Solution: (5.00000, 3.00000)
Residual: (0.323E-11, 0.705E-11)
Iterations: 11
Function Evaluations: 15
Jacobian Evaluations: 1
```

## Example 2
This example uses a least-squares approach to determine the coefficients of a polynomial that best fits a set of data.

```fortran
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
```
The above program produces the following output.
```text
c0: 1.1866142244
c1: 0.4466134462
c2: -.1223202909
c3: 1.0647627571
Max Residual: 0.50636
```

The following graph illustrates the fit.
![](images/Curve_Fit_Example_1.png?raw=true)

## Example 3
This example utilizes the polynomial type to fit a polynomial to the data set utilized in Example 2.
```fortran
program example
    use iso_fortran_env
    use nonlin_polynomials
    implicit none

    ! Local Variables
    integer(int32) :: i
    real(real64), dimension(21) :: xp, yp, yf, yc, err
    real(real64) :: res
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
```
The above program yields the following coefficients.
```text
c0 = 1.1866141861
c1 = 0.4466136311
c2 = -.1223204989
c3 = 1.0647628218
Max Residual: 0.50636
```
Notice, as expected, the results are very similar to the output of Example 2.

## Example 4
This example uses the Nelder-Mead simplex method to find the minimum of the Rosenbrock function.
```fortran
program example
    use iso_fortran_env
    use nonlin_optimize, only : nelder_mead
    use nonlin_core
    implicit none

    ! Local Variables
    type(nelder_mead) :: solver
    type(fcnnvar_helper) :: obj
    procedure(fcnnvar), pointer :: fcn
    real(real64) :: x(2), fout
    type(iteration_behavior) :: ib

    ! Initialization
    fcn => rosenbrock
    call obj%set_fcn(fcn, 2)

    ! Define an initial guess - the solution is (1, 1)
    call random_number(x)

    ! Call the solver
    call solver%solve(obj, x, fout, ib)

     ! Display the output
     print '(AF7.5AF7.5A)', "Minimum: (", x(1), ", ", x(2), ")"
     print '(AE9.3)', "Function Value: ", fout
     print '(AI0)', "Iterations: ", ib%iter_count
     print '(AI0)', "Function Evaluations: ", ib%fcn_count
contains
    ! Rosenbrock's Function
    function rosenbrock(x) result(f)
        real(real64), intent(in), dimension(:) :: x
        real(real64) :: f
        f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    end function
end program
```
The above program produces the following output:
```text
Minimum: (1.00000, 1.00000)
Function Value: 0.121E-12
Iterations: 52
Function Evaluations: 101
```
Notice, the convergence tolerance was set to its default value (1e-12).
