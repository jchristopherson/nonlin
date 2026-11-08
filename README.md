# nonlin

A Fortran library for nonlinear equations, nonlinear least squares, local optimization, and polynomial operations.

## Status

[![CMake](https://github.com/jchristopherson/nonlin/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/nonlin/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/nonlin/workflows/fpm/badge.svg)](https://github.com/jchristopherson/nonlin/actions)

## Documentation

See the [NONLIN API documentation](https://jchristopherson.github.io/nonlin/).

## Capabilities Matrix

Here, M is the number of residual equations and N is the number of unknowns. AD means forward-mode automatic differentiation through [autodiff](https://github.com/jchristopherson/autodiff).

| Solver | Problem | Finite-difference derivatives | User derivatives | AD | Bounds or bracket | Parallel Jacobian |
| --- | --- | --- | --- | --- | --- | --- |
| `brent_solver` | One scalar root | Not needed | Not needed | No | Root bracket | No |
| `newton_1var_solver` | One scalar root | Yes | First derivative | No | Root bracket | No |
| `newton_solver` | Nonlinear system, M = N | Yes | Jacobian | Yes | No | AD and finite differences |
| `quasi_newton_solver` | Nonlinear system, M = N | Yes | Jacobian | Yes | No | AD and finite differences |
| `least_squares_solver` | Minimize sum of squared residuals, M >= N | Yes | Jacobian | Yes | No | AD and finite differences |
| `constrained_least_squares_solver` | Minimize sum of squared residuals, M >= N | Yes | Jacobian | Yes | Variable box bounds | AD and finite differences |
| `nelder_mead` | Scalar objective, N variables | Not needed | Not needed | Function values only | No | No |
| `bfgs` | Scalar objective, N variables | Yes | Gradient | Yes | No | No |

The library also provides polynomial evaluation, fitting, arithmetic, and root calculations through `polynomial`. Underdetermined residual systems (M < N), general nonlinear constraints, and global optimization are not provided by the solvers above.

### Derivatives and Threading

- Vector solvers use `vecfcn_helper`: register a real-valued callback with `set_fcn`, an AD-compatible callback with `set_auto_fcn`, and an optional analytic Jacobian with `set_jacobian`.
- Scalar-objective optimizers use `fcnnvar_helper`: register callbacks with `set_fcn` or `set_auto_fcn`, and an optional analytic gradient with `set_gradient_fcn`.
- Scalar root solvers use `fcn1var_helper`; `set_diff` supplies an optional analytic first derivative for Newton's method.
- Analytic Jacobians/gradients take priority over AD; finite differences are the fallback. For function values, `vecfcn_helper` prefers the AD callback, while `fcnnvar_helper` prefers the ordinary callback if both are registered. Register only the AD callback for AD-only scalar-objective evaluation; when supplying both callbacks, ensure they describe the same function.
- Vector Jacobian evaluation is serial by default. Set `obj%set_thread_count(n)` to request OpenMP workers for AD or finite-difference columns. This does not parallelize entire solver iterations or automatically parallelize user-supplied Jacobians.
- Multi-threaded callbacks must be thread-safe. Do not evaluate the same helper concurrently from separate callers; its workspaces are mutable. More threads are not necessarily faster for small or inexpensive problems.

## Available Solvers

These are local methods: starting points, smoothness, scaling, and tolerances affect convergence. For least squares, convergence can occur at a nonzero residual; inspect both the returned residual/objective and `iteration_behavior` rather than treating every stopping flag as an exact root.

| Solver | Method and Description | Recommendation |
| --- | --- | --- |
| `brent_solver` | Bracket-preserving combination of bisection, secant steps, and inverse quadratic interpolation; no derivative required. | First choice for a continuous scalar function with a known sign-changing bracket. |
| `newton_1var_solver` | Safeguarded scalar Newton iteration within a root bracket; uses an analytic derivative or finite differences. | Use for smooth scalar equations when derivative information is reliable; prefer Brent when derivatives are unreliable or expensive. |
| `newton_solver` | Solves the linearized square system using a freshly evaluated Jacobian each iteration; optional backtracking line search, enabled by default. | First choice for smooth square systems, especially with AD or an analytic Jacobian and a reasonable starting point. |
| `quasi_newton_solver` | Broyden rank-one Jacobian updates with periodic/restarted Jacobian evaluations and optional backtracking, enabled by default. | Use when repeated Jacobian evaluations are costly. Tune `set_jacobian_interval`; prefer Newton when the approximation becomes unreliable. |
| `least_squares_solver` | Levenberg-Marquardt method with damping and variable scaling. | First choice for unconstrained data fitting or overdetermined systems; also useful when a square system needs a least-squares approach. |
| `constrained_least_squares_solver` | Bound-aware trust-region dogleg method with backtracking and projected-gradient steps at active bounds. | Use for residual minimization with lower/upper variable limits. Set bounds with `set_lower_limits` and `set_upper_limits`. |
| `nelder_mead` | Derivative-free simplex search for a local minimum of a scalar objective. | Use for modest-dimensional problems when useful derivatives are unavailable; generally less efficient than gradient methods for smooth, larger problems. |
| `bfgs` | Quasi-Newton scalar minimization with gradient evaluations, curvature updates, and line search. | First choice for smooth unconstrained scalar objectives with AD or analytic gradients. Prefer a least-squares solver when the objective naturally consists of residuals. |

## Building NONLIN

### Prerequisites

- A Fortran 2018-capable compiler with OpenMP support, such as a recent GNU Fortran (`gfortran`). Use a consistent compiler/toolchain for NONLIN and its dependencies because Fortran module files are compiler-specific.
- Git and internet access for dependencies that are not already installed or cached.
- BLAS and LAPACK libraries compatible with the compiler. FPM also needs ARPACK for the current `linalg` dependency. For example, on Debian/Ubuntu, install `gfortran`, `git`, `libblas-dev`, `liblapack-dev`, and `libarpack2-dev`, plus the chosen build tool. CMake can fetch/build ARPACK when required by `linalg`.
- For CMake: [CMake](https://cmake.org/) 3.24 or newer and a build backend such as Make or Ninja. The current CMake configuration requires OpenMP.
- For FPM: a recent [Fortran Package Manager](https://fpm.fortran-lang.org/) with [OpenMP metapackage support](https://fpm.fortran-lang.org/spec/metapackages.html#openmp).

On Windows, run commands from a shell where the compiler and build backend are on PATH (for example, a configured MinGW-w64 or compiler developer shell). Native Windows PowerShell and MSYS2 shells may require different CMake generators; use the one matching your toolchain.

Start from the repository root:

```sh
git clone https://github.com/jchristopherson/nonlin.git
cd nonlin
```

### CMake

Configure an out-of-source Release build with tests and examples explicitly enabled:

```sh
cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DBUILD_NONLIN_EXAMPLES=ON
cmake --build build-cmake --parallel
ctest --test-dir build-cmake --output-on-failure
```

`build-cmake` keeps CMake output separate from FPM's `build` directory. CMake uses installed `linalg` (version 3.0.0 or newer) and `autodiff` (version 1.0.0 or newer) packages when found; otherwise it downloads and builds them. Tests also require `fortran_test_helper`, which is fetched when necessary. Transitive dependencies may also be downloaded on the first configuration.

| CMake Option | Purpose |
| --- | --- |
| `CMAKE_Fortran_COMPILER` | Select the compiler, e.g. `-DCMAKE_Fortran_COMPILER=gfortran`; set this when creating a new build directory. |
| `CMAKE_BUILD_TYPE` | `Release` for optimized builds or `Debug` for debugging with single-configuration generators. |
| `BUILD_SHARED_LIBS` | `OFF` (default) builds static libraries; `ON` builds shared libraries. |
| `BUILD_TESTING` | Enable/disable the test executables and CTest registration. |
| `BUILD_NONLIN_EXAMPLES` | Enable/disable example executables (default `OFF`). |
| `CMAKE_PREFIX_PATH` | Search prefix for preinstalled dependency packages. |
| `BLA_VENDOR` | Select a BLAS/LAPACK vendor through CMake's discovery, e.g. `OpenBLAS`. |

For a library-only build, configure with `-DBUILD_TESTING=OFF -DBUILD_NONLIN_EXAMPLES=OFF`. To select a backend explicitly, add `-G Ninja` or, for MinGW Make, `-G "MinGW Makefiles"` to the configure command. Use a fresh build directory when changing compiler or generator.

For multi-configuration generators, select the configuration at build and test time instead:

```sh
cmake --build build-cmake --config Release --parallel
ctest --test-dir build-cmake -C Release --output-on-failure
```

Install to a writable prefix after a successful build:

```sh
cmake --install build-cmake --prefix ./install --config Release
```

The installation includes libraries, Fortran module files, and CMake/pkg-config metadata. Shared-library users must also make the installed runtime libraries discoverable by their operating system.

### FPM

The supplied manifest declares the Fortran dependencies and OpenMP support. FPM downloads the Git dependencies and enables compiler-appropriate OpenMP flags automatically. NONLIN links against `blas` and `lapack`, and the current `linalg` dependency also links against `arpack`; install these first and make their library directories visible to the linker. A library-only build can succeed without these libraries, but linking tests or applications requires them.

```sh
fpm build --profile release --build-dir build-fpm
fpm test --profile debug --build-dir build-fpm
fpm test --profile debug --build-dir build-fpm -- jacobian_threads
```

Use a recent FPM version supporting `--build-dir`; `build-fpm` avoids collisions with existing CMake output or old dependency caches. The last command runs the threaded AD/finite-difference Jacobian regression through the existing test executable. The manifest disables automatic example and executable discovery; use CMake with `BUILD_NONLIN_EXAMPLES=ON` to build the supplied examples.

Select another compiler with `--compiler`, for example `fpm build --compiler gfortran --profile release --build-dir build-fpm`. If BLAS/LAPACK/ARPACK are installed outside the default search path, pass their directory with `--link-flag "-L/path/to/libraries"` when building or testing with GNU Fortran. Linker errors such as `cannot find -lblas`, `-llapack`, or `-larpack` mean these libraries or their search paths are missing. On Windows, an OpenBLAS installation may expose only an `openblas` library name, not the `blas`/`lapack` names the manifest requests; use compatible linker aliases or choose CMake's vendor discovery instead. CMake is usually easier with a vendor-specific BLAS/LAPACK package such as MKL.

To install the library and module files to a local prefix:

```sh
fpm install --profile release --build-dir build-fpm --prefix ./install-fpm
```

To use NONLIN in another FPM project, add this to that project's manifest:

```toml
[dependencies]
nonlin = { git = "https://github.com/jchristopherson/nonlin" }
```

For reproducible builds, add a `rev` specifying the NONLIN commit you have validated. If your own source directly imports `autodiff`, declare it as a direct dependency too. The application still needs compatible BLAS/LAPACK/ARPACK libraries available to its linker.

## External Libraries

| Dependency | Role |
| --- | --- |
| [BLAS](https://www.netlib.org/blas/) and [LAPACK](https://www.netlib.org/lapack/) | Dense numerical linear algebra. |
| [ARPACK](https://github.com/opencollab/arpack-ng) | Eigenvalue routines required by the current `linalg` dependency. |
| [linalg](https://github.com/jchristopherson/linalg) | Linear-algebra interfaces and algorithms, with additional transitive dependencies. |
| [autodiff](https://github.com/jchristopherson/autodiff) | Forward-mode automatic differentiation. |
| [fortran_test_helper](https://github.com/jchristopherson/fortran_test_helper) | Test assertions; required only for testing. |

## Example 1

This example solves a set of two equations of two unknowns using a Quasi-Newton type solver.  In this example, the solver is left to compute the derivatives numerically.

```fortran
program  example
    use iso_fortran_env
    use nonlin
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
    print '(A,F7.5,A,F7.5,A)', "Solution: (", x(1), ", ", x(2), ")"
    print '(A,E9.3,A,E9.3,A)', "Residual: (", f(1), ", ", f(2), ")"
    print '(A,I0)', "Iterations: ", ib%iter_count
    print '(A,I0)', "Function Evaluations: ", ib%fcn_count
    print '(A,I0)', "Jacobian Evaluations: ", ib%jacobian_count

contains
    ! Define the routine containing the equations to solve.  The equations are:
    ! x**2 + y**2 = 34
    ! x**2 - 2 * y**2 = 7
    subroutine fcns(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
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
    use nonlin
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
    print '(A,F12.10)', "c0: ", x(4)
    print '(A,F12.10)', "c1: ", x(3)
    print '(A,F12.10)', "c2: ", x(2)
    print '(A,F12.10)', "c3: ", x(1)
    print '(A,F7.5)', "Max Residual: ", maxval(abs(f))

contains
    ! The function containing the data to fit
    subroutine fcns(x, f, args)
        ! Arguments
        real(real64), intent(in), dimension(:) :: x  ! Contains the coefficients
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args

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
![Cubic polynomial fit to the example data](images/Curve_Fit_Example_1.png?raw=true)

## Example 3

This example utilizes the polynomial type to fit a polynomial to the data set utilized in Example 2.

```fortran
program example
    use iso_fortran_env
    use nonlin
    implicit none

    ! Local Variables
    integer :: i
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
    print '(A,I0,A,F12.10)', ("c", i - 1, " = ", p%get(i), i = 1, 4)
    print '(A,F7.5)', "Max Residual: ", res
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
    use nonlin
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
     print '(A,F7.5,A,F7.5,A)', "Minimum: (", x(1), ", ", x(2), ")"
     print '(A,E9.3)', "Function Value: ", fout
     print '(A,I0)', "Iterations: ", ib%iter_count
     print '(A,I0)', "Function Evaluations: ", ib%fcn_count
contains
    ! Rosenbrock's Function
    function rosenbrock(x, args) result(f)
        real(real64), intent(in), dimension(:) :: x
        class(*), intent(in), optional :: args
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
