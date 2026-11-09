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
| `continuation_solver` | Nonlinear system, M = N | Yes | Jacobian | Yes | No | AD and finite differences |
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
| `continuation_solver` | Selectable fixed-parameter or pseudo-arclength continuation with a damped Newton corrector. The default homotopy is $H(x,\lambda)=(1-\lambda)(x-x_0)+\lambda F(x)$; its derivatives use `vecfcn_helper`. | Use fixed-parameter mode for well-behaved paths. Enable pseudo-arclength with `set_use_pseudo_arclength(.true.)` to follow regular folds where lambda reverses direction. It does not guarantee passage through singular branch points or selection of a desired branch. |
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

## Choosing an API

NONLIN offers simple wrapper functions and configurable solver objects:

| API | Use it when | Trade-off |
| --- | --- | --- |
| Simple functions such as `solve_brent`, `solve_newton`, `solve_quasi_newton`, `solve_nelder_mead`, `solve_continuation`, `solve_continuation_pseudo_arclength`, and `solve_bfgs` | You want a short call with default solver settings and an initial guess. Callback procedures are passed directly. | Less setup; solver controls, custom derivatives, and iteration statistics are not exposed by the wrapper. |
| Solver objects such as `brent_solver`, `newton_solver`, `quasi_newton_solver`, `continuation_solver`, `nelder_mead`, and `bfgs` | You need to configure tolerances, bounds/brackets, derivative callbacks, line search behavior, continuation steps, or inspect iteration statistics. | More setup; you provide helper objects and configure the solver explicitly. |

For example, the simple BFGS wrapper takes an AD-compatible objective and an initial guess:

```fortran
x = solve_bfgs(rosenbrock, initial_guess)
```

The solver-object API makes setup and controls explicit:

```fortran
call objective%set_fcn(fcn, 2)
call solver%set_fcn_tolerance(1.0d-10)
call solver%solve(objective, x, value)
```

The simple wrappers cover scalar root-finding (`solve_brent`, `solve_newton_1var`), square nonlinear systems (`solve_newton`, `solve_quasi_newton`), and unconstrained scalar-objective minimization (`solve_nelder_mead`, `solve_bfgs`).  Use solver objects for least-squares problems or whenever you need more control or reporting.

## Example 1: Quasi-Newton System Solve

Solve two equations in two unknowns. This version estimates the Jacobian numerically; configure the solver or provide derivatives with `vecfcn_helper` when needed.

```fortran
type(vecfcn_helper) :: equations
type(quasi_newton_solver) :: solver
procedure(vecfcn), pointer :: fcn
real(real64) :: x(2), residual(2)

fcn => misc_2fcn
call equations%set_fcn(fcn, 2, 2)
x = 1.0d0
call solver%solve(equations, x, residual)
```

The unknowns converge to `(5, 3)`. See [examples/nonlin_quasi_newton_example.f90](examples/nonlin_quasi_newton_example.f90) for the full program.

## Example 2: Polynomial Fit with Least Squares

Fit a cubic model to 21 measured data points by minimizing the residual vector over four coefficients. Least squares is available through the solver-object API.

```fortran
type(vecfcn_helper) :: residuals
type(least_squares_solver) :: solver
procedure(vecfcn), pointer :: fcn
real(real64) :: coefficients(4), residual(21)

fcn => lsq_poly_fit_fcn
call residuals%set_fcn(fcn, 21, 4)
coefficients = 1.0d0
call solver%solve(residuals, coefficients, residual)
```

Source: [examples/nonlin_least_squares_example.f90](examples/nonlin_least_squares_example.f90).

## Example 3: Polynomial Utilities

Use the `polynomial` type for direct polynomial fitting without configuring a nonlinear solver. `x_data` and `y_data` are the sampled points.

```fortran
type(polynomial) :: p

call p%fit(x_data, y_data, 3)
y_fit = p%evaluate(x_data)
```

Source: [examples/nonlin_poly_fit_example.f90](examples/nonlin_poly_fit_example.f90).

## Example 4: Nelder-Mead Minimization

Minimize the Rosenbrock objective without derivatives:

```fortran
type(fcnnvar_helper) :: objective
type(nelder_mead) :: solver
procedure(fcnnvar), pointer :: fcn
real(real64) :: x(2), value

fcn => rosenbrock
call objective%set_fcn(fcn, 2)
x = [0.0d0, 0.0d0]
call solver%solve(objective, x, value)
```

Source: [examples/nonlin_nm_opt_example.f90](examples/nonlin_nm_opt_example.f90).

## Example 5: Simple API with BFGS

The convenience wrapper minimizes a smooth objective with default solver settings. Its callback accepts dual numbers, allowing NONLIN to obtain gradients by AD without a procedure pointer or a helper object.

```fortran
use nonlin
use autodiff, only : dual

real(real64) :: x(2), initial_guess(2)
initial_guess = [0.0d0, 0.0d0]
x = solve_bfgs(rosenbrock, initial_guess)

contains
    function rosenbrock(x, args) result(f)
        type(dual), intent(in), dimension(:) :: x
        class(*), intent(in), optional :: args
        type(dual) :: f
        f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    end function
```

See [examples/nonlin_simple_bfgs_example.f90](examples/nonlin_simple_bfgs_example.f90) for the complete program and [examples/nonlin_bfgs_example.f90](examples/nonlin_bfgs_example.f90) for the configurable solver-object approach.

## Example 6: Powell's Singular Function

This four-equation system has the root `x = [0, 0, 0, 0]`. Its last two equations are squares, so their first derivatives vanish at the root. The root Jacobian has rank 2 rather than 4, making this a singular endpoint rather than a regular fold.

```fortran
type(vecfcn_helper) :: equations
type(continuation_solver) :: solver
procedure(vector_function), pointer :: fcn
real(real64) :: x(4), residual(4)

fcn => fcn_
call equations%set_auto_fcn(fcn, 4, 4)
x = 1.0d0
call solver%solve(equations, x, residual)
```

The default fixed-parameter continuation succeeds on this example. Pseudo-arclength is not a suitable solver here: at the endpoint, the augmented derivative `[J_F, H_lambda]` also loses row rank, so the solution path has no unique tangent for the pseudo-arclength corrector to follow. This differs from a regular fold, where `J_F` may be singular while the augmented derivative remains full row rank. See [examples/nonlin_continuation_example.f90](examples/nonlin_continuation_example.f90) for the full Powell system and convergence-history output.
