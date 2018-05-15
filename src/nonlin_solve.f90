! nonlin_solve.f90

!> @brief \b nonlin_solve
!!
!! @par Purpose
!! To provide various routines capapble of solving systems of nonlinear
!! equations.
module nonlin_solve
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use nonlin_constants
    use nonlin_core
    use nonlin_linesearch, only : line_search, limit_search_vector
    use ferror, only : errors
    use linalg_core, only : qr_factor, form_qr, qr_rank1_update, lu_factor, &
        rank1_update, mtx_mult, recip_mult_array, solve_triangular_system, &
        solve_lu
    implicit none
    private
    public :: line_search_solver
    public :: quasi_newton_solver
    public :: newton_solver
    public :: brent_solver
    public :: newton_1var_solver
    public :: test_convergence

! ******************************************************************************
! NONLIN_SOLVE_LINE_SEARCH.F90
! ------------------------------------------------------------------------------
    !> @brief A class describing nonlinear solvers that use a line search
    !! algorithm to improve convergence behavior.
    type, abstract, extends(equation_solver) :: line_search_solver
        private
        !> The line search module.
        class(line_search), allocatable :: m_lineSearch
        !> Set to true if a line search should be used regardless of the status
        !! of m_lineSearch
        logical :: m_useLineSearch = .true.
    contains
        !> @brief Gets the line search module.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine get_line_search(class(line_search_solver) this, class(line_search) ls)
        !! @endcode
        !!
        !! @param[in] this The line_search_solver object.
        !! @param[out] ls The line_search object.
        procedure, public :: get_line_search => lss_get_line_search
        !> @brief Sets the line search module.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_line_search(class(line_search_solver) this, class(line_search) ls)
        !! @endcode
        !!
        !! @param[in,out] this The line_search_solver object.
        !! @param[in] ls The line_search object.
        procedure, public :: set_line_search => lss_set_line_search
        !> @brief Establishes a default line_search object for the line search
        !! module.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_default_line_search(class(line_search_solver) this)
        !! @endcode
        !!
        !! @param[in,out] this The line_search_solver object.
        procedure, public :: set_default_line_search => lss_set_default
        !> @brief Tests to see if a line search module is defined.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function is_line_search_defined(class(line_search_solver) this)
        !! @endcode
        !!
        !! @param[in] this The line_search_solver object.
        !! @return Returns true if a module is defined; else, false.
        procedure, public :: is_line_search_defined => &
            lss_is_line_search_defined
        !> @brief Gets a value determining if a line-search should be employed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function get_use_line_search(class(line_search_solver) this)
        !! @endcode
        !!
        !! @param[in] this The line_search_solver object.
        !! @return Returns true if a line search should be used; else, false.
        procedure, public :: get_use_line_search => lss_get_use_search
        !> @brief Sets a value determining if a line-search should be employed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_use_line_search(class(line_search_solver) this, logical x)
        !! @endcode
        !!
        !! @param[in,out] this The line_search_solver object.
        !! @param[in] x Set to true if a line search should be used; else, false.
        procedure, public :: set_use_line_search => lss_set_use_search
    end type

! ------------------------------------------------------------------------------
    interface
        module subroutine lss_get_line_search(this, ls)
            class(line_search_solver), intent(in) :: this
            class(line_search), intent(out), allocatable :: ls
        end subroutine

        module subroutine lss_set_line_search(this, ls)
            class(line_search_solver), intent(inout) :: this
            class(line_search), intent(in) :: ls
        end subroutine

        module subroutine lss_set_default(this)
            class(line_search_solver), intent(inout) :: this
            type(line_search) :: ls
        end subroutine

        pure module function lss_is_line_search_defined(this) result(x)
            class(line_search_solver), intent(in) :: this
            logical :: x
        end function

        pure module function lss_get_use_search(this) result(x)
            class(line_search_solver), intent(in) :: this
            logical :: x
        end function

        module subroutine lss_set_use_search(this, x)
            class(line_search_solver), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine
    end interface

! ******************************************************************************
! NONLIN_SOLVE_QUASI_NEWTON.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a quasi-Newton type solver based upon Broyden's method.
    type, extends(line_search_solver) :: quasi_newton_solver
        private
        !> The number of iterations that may pass between Jacobian calculation.
        integer(int32) :: m_jDelta = 5
    contains
        !> @brief Applies the quasi-Newton's method developed by Broyden in
        !! conjunction with a backtracking type line search to solve N equations
        !! of N unknowns.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine solve(class(quasi_newton_solver) this, class(vecfcn_helper) fcn, real(real64) x(:), real(real64) fvec(:), optional type(iteration_behavior) ib, optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver-based object.
        !! @param[in] fcn The vecfcn_helper object containing the equations to
        !!  solve.
        !! @param[in,out] x On input, an N-element array containing an initial
        !!  estimate to the solution.  On output, the updated solution estimate.
        !!  N is the number of variables.
        !! @param[out] fvec An N-element array that, on output, will contain
        !!  the values of each equation as evaluated at the variable values
        !!  given in @p x.
        !! @param[out] ib An optional output, that if provided, allows the
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can be
        !!  used to retrieve information relating to any errors encountered during
        !!  execution.  If not provided, a default implementation of the errors
        !!  class is used internally to provide error handling.  Possible errors and
        !!  warning messages that may be encountered are as follows.
        !!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
        !!  - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is different
        !!      than the number of variables.
        !!  - NL_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
        !!      correctly.
        !!  - NL_DIVERGENT_BEHAVIOR_ERROR: Occurs if the direction vector is
        !!      pointing in an apparent uphill direction.
        !!  - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
        !!      the allowed number of iterations.
        !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - NL_SPURIOUS_CONVERGENCE_ERROR: Occurs as a warning if the slope of the
        !!      gradient vector becomes sufficiently close to zero.
        !!
        !! @par Example
        !! The following code provides an example of how to solve a system of N
        !! equations of N unknonwns using this Quasi-Newton method.
        !! @code{.f90}
        !! program main
        !!     use iso_fortran_env
        !!     use nonlin_core, only : vecfcn, vecfcn_helper
        !!     use nonlin_solve, only : quasi_newton_solver
        !!
        !!     type(vecfcn_helper) :: obj
        !!     procedure(vecfcn), pointer :: fcn
        !!     type(quasi_newton_solver) :: solver
        !!     real(real64) :: x(2), f(2)
        !!
        !!     ! Set the initial conditions to [1, 1]
        !!     x = 1.0d0
        !!
        !!     ! Define the function
        !!     fcn => fcn1
        !!     call obj%set_fcn(fcn, 2, 2)
        !!
        !!     ! Solve the system of equations.  The solution overwrites X
        !!     call solver%solve(obj, x, f)
        !!
        !!     ! Print the output and the residual:
        !!     print '(AF5.3AF5.3A)', "The solution: (", x(1), ", ", x(2), ")"
        !!     print '(AE8.3AE8.3A)', "The residual: (", f(1), ", ", f(2), ")"
        !! contains
        !!     ! System of Equations:
        !!     !
        !!     ! x**2 + y**2 = 34
        !!     ! x**2 - 2 * y**2 = 7
        !!     !
        !!     ! Solution:
        !!     ! x = +/-5
        !!     ! y = +/-3
        !!     subroutine fcn1(x, f)
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64), intent(out), dimension(:) :: f
        !!         f(1) = x(1)**2 + x(2)**2 - 34.0d0
        !!         f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
        !!     end subroutine
        !! end program
        !! @endcode
        !! The above program returns the following results.
        !! @code{.txt}
        !! The solution: (5.000, 3.000)
        !! The residual: (.604E-10, .121E-09)
        !! @endcode
        !!
        !! @par See Also
        !! - [Broyden's Paper](http://www.ams.org/journals/mcom/1965-19-092/S0025-5718-1965-0198670-6/S0025-5718-1965-0198670-6.pdf)
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Broyden%27s_method)
        !! - [Numerical Recipes](http://numerical.recipes/)
        procedure, public :: solve => qns_solve
        !> @brief Gets the number of iterations that may pass before forcing a
        !! recalculation of the Jacobian matrix.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_jacobian_interval(class(quasi_newton_solver) this)
        !! @endcode
        !!
        !! @param[in] this The quasi_newton_solver object.
        !! @return The number of iterations.
        procedure, public :: get_jacobian_interval => qns_get_jac_interval
        !> @brief Sets the number of iterations that may pass before forcing a
        !! recalculation of the Jacobian matrix.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_jacobian_interval(class(quasi_newton_solver) this, integer(int32) n)
        !! @endcode
        !!
        !! @param[in,out] this The quasi_newton_solver object.
        !! @param[in] n The number of iterations.
        procedure, public :: set_jacobian_interval => qns_set_jac_interval
    end type

! ------------------------------------------------------------------------------
    interface
        module subroutine qns_solve(this, fcn, x, fvec, ib, err)
            class(quasi_newton_solver), intent(inout) :: this
            class(vecfcn_helper), intent(in) :: fcn
            real(real64), intent(inout), dimension(:) :: x
            real(real64), intent(out), dimension(:) :: fvec
            type(iteration_behavior), optional :: ib
            class(errors), intent(inout), optional, target :: err
        end subroutine

        pure module function qns_get_jac_interval(this) result(n)
            class(quasi_newton_solver), intent(in) :: this
            integer(int32) :: n
        end function

        module subroutine qns_set_jac_interval(this, n)
            class(quasi_newton_solver), intent(inout) :: this
            integer(int32), intent(in) :: n
        end subroutine
    end interface

! ******************************************************************************
! ------------------------------------------------------------------------------
    !> @brief Defines a Newton solver.
    type, extends(line_search_solver) :: newton_solver
    contains
        !> @brief Applies Newton's method in conjunction with a backtracking type
        !! line search to solve N equations of N unknowns.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine solve(class(newton_solver) this, class(vecfcn_helper) fcn, real(real64) x(:), real(real64) fvec(:), optional type(iteration_behavior) ib, optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver-based object.
        !! @param[in] fcn The vecfcn_helper object containing the equations to
        !!  solve.
        !! @param[in,out] x On input, an N-element array containing an initial
        !!  estimate to the solution.  On output, the updated solution estimate.
        !!  N is the number of variables.
        !! @param[out] fvec An N-element array that, on output, will contain
        !!  the values of each equation as evaluated at the variable values
        !!  given in @p x.
        !! @param[out] ib An optional output, that if provided, allows the
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can be
        !!  used to retrieve information relating to any errors encountered during
        !!  execution.  If not provided, a default implementation of the errors
        !!  class is used internally to provide error handling.  Possible errors and
        !!  warning messages that may be encountered are as follows.
        !!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
        !!  - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is different
        !!      than the number of variables.
        !!  - NL_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
        !!      correctly.
        !!  - NL_DIVERGENT_BEHAVIOR_ERROR: Occurs if the direction vector is
        !!      pointing in an apparent uphill direction.
        !!  - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
        !!      the allowed number of iterations.
        !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - NL_SPURIOUS_CONVERGENCE_ERROR: Occurs as a warning if the slope of the
        !!      gradient vector becomes sufficiently close to zero.
        !!
        !! @par Example
        !! The following code provides an example of how to solve a system of N
        !! equations of N unknonwns using Newton's method.
        !! @code{.f90}
        !! program main
        !!     use nonlin_core, only : vecfcn, vecfcn_helper
        !!     use nonlin_solve, only : newton_solver
        !!
        !!     type(vecfcn_helper) :: obj
        !!     procedure(vecfcn), pointer :: fcn
        !!     type(newton_solver) :: solver
        !!     real(real64) :: x(2), f(2)
        !!
        !!     ! Set the initial conditions to [1, 1]
        !!     x = 1.0d0
        !!
        !!     ! Define the function
        !!     fcn => fcn1
        !!     call obj%set_fcn(fcn, 2, 2)
        !!
        !!     ! Solve the system of equations.  The solution overwrites X
        !!     call solver%solve(obj, x, f)
        !!
        !!     ! Print the output and the residual:
        !!     print '(AF5.3AF5.3A)', "The solution: (", x(1), ", ", x(2), ")"
        !!     print '(AE8.3AE8.3A)', "The residual: (", f(1), ", ", f(2), ")"
        !! contains
        !!     ! System of Equations:
        !!     !
        !!     ! x**2 + y**2 = 34
        !!     ! x**2 - 2 * y**2 = 7
        !!     !
        !!     ! Solution:
        !!     ! x = +/-5
        !!     ! y = +/-3
        !!     subroutine fcn1(x, f)
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64), intent(out), dimension(:) :: f
        !!         f(1) = x(1)**2 + x(2)**2 - 34.0d0
        !!         f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
        !!     end subroutine
        !! end program
        !! @endcode
        !! The above program returns the following results.
        !! @code{.txt}
        !! The solution: (5.000, 3.000)
        !! The residual: (.000E+00, .000E+00)
        !! @endcode
        !!
        !! @par See Also
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Newton%27s_method)
        procedure, public :: solve => ns_solve
    end type

! ------------------------------------------------------------------------------
    interface
        module subroutine ns_solve(this, fcn, x, fvec, ib, err)
            class(newton_solver), intent(inout) :: this
            class(vecfcn_helper), intent(in) :: fcn
            real(real64), intent(inout), dimension(:) :: x
            real(real64), intent(out), dimension(:) :: fvec
            type(iteration_behavior), optional :: ib
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ******************************************************************************
! NONLIN_SOLVE_BRENT.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a solver based upon Brent's method for solving an equation
    !! of one variable without using derivatives.
    type, extends(equation_solver_1var) :: brent_solver
    contains
        !> @brief Solves the equation.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine solve(class(brent_solver) this, class(fcn1var_helper) fcn, real(real64) x, type(value_pair) lim, optional real(real64) f, optional type(iteration_behavior) ib, optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The brent_solver object.
        !! @param[in] fcn The fcn1var_helper object containing the equation
        !!  to solve.
        !! @param[in,out] x A parameter used to return the solution.  Notice, any
        !!  input value will be ignored as this routine relies upon the
        !!  search limits in @p lim to provide a starting point.
        !! @param[in] lim A value_pair object defining the search limits.
        !! @param[out] f An optional parameter used to return the function
        !!  residual as computed at @p x.
        !! @param[out] ib An optional output, that if provided, allows the
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can be
        !!  used to retrieve information relating to any errors encountered during
        !!  execution.  If not provided, a default implementation of the errors
        !!  class is used internally to provide error handling.  Possible errors and
        !!  warning messages that may be encountered are as follows.
        !!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
        !!  - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is different
        !!      than the number of variables.
        !!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
        !!      the allowed number of iterations.
        !!
        !! @par Example
        !! The following code provides an example of how to solve an equation of
        !! one variable using Brent's method.
        !! @code{.f90}
        !! program main
        !!     use nonlin_core, only : fcn1var, fcn1var_helper, value_pair
        !!     use nonlin_solve, only : brent_solver
        !!
        !!     type(fcn1var_helper) :: obj
        !!     procedure(fcn1var), pointer :: fcn
        !!     type(brent_solver) :: solver
        !!     real(real64) :: x, f
        !!     type(value_pair) :: limits
        !!
        !!     ! Define the solution limits
        !!     limits%x1 = 1.5d0
        !!     limits%x2 = 5.0d0
        !!
        !!     ! Define the function
        !!     fcn => fcn1
        !!     call obj%set_fcn(fcn)
        !!
        !!     ! Solve the equation
        !!     call solver%solve(obj, x, limits, f)
        !!
        !!     ! Print the output and the residual:
        !!     print '(AF5.3)', "The solution: ", x
        !!     print '(AE9.3)', "The residual: ", f
        !! contains
        !!     ! f(x) = sin(x) / x, SOLUTION: x = n * pi for n = 1, 2, 3, ...
        !!     function fcn1(x) result(f)
        !!         real(real64), intent(in) :: x
        !!         real(real64) :: f
        !!         f = sin(x) / x
        !!     end function
        !! end program
        !! @endcode
        !! The above program returns the following results.
        !! @code{.txt}
        !! The solution: 3.142
        !! The residual: -.751E-11
        !! @endcode
        !!
        !! @par See Also
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Brent%27s_method)
        !! - [Numerical Recipes](http://numerical.recipes/)
        !! - R.P. Brent, "Algorithms for Minimization without Derivatives,"
        !!      Dover Publications, January 2002. ISBN 0-486-41998-3.
        !!      Further information available
        !!      [here](https://maths-people.anu.edu.au/~brent/pub/pub011.html).
        procedure, public :: solve => brent_solve
    end type

! ------------------------------------------------------------------------------
    interface
        module subroutine brent_solve(this, fcn, x, lim, f, ib, err)
            class(brent_solver), intent(inout) :: this
            class(fcn1var_helper), intent(in) :: fcn
            real(real64), intent(inout) :: x
            type(value_pair), intent(in) :: lim
            real(real64), intent(out), optional :: f
            type(iteration_behavior), optional :: ib
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ******************************************************************************
! NONLIN_SOLVE_NEWTON1VAR.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a solver based upon Newtons's method for solving an
    !! equation of one variable.  The algorithm uses a bisection method in
    !! conjunction with Newton's method in order to keep bounds upon the
    !! Newton iterations.
    type, extends(equation_solver_1var) :: newton_1var_solver
    contains
        !> @brief Solves the equation.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine solve(class(newton_1var_solver) this, class(fcn1var_helper) fcn, real(real64) x, type(value_pair) lim, optional real(real64) f, optional type(iteration_behavior) ib, optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The brent_solver object.
        !! @param[in] fcn The fcn1var_helper object containing the equation
        !!  to solve.
        !! @param[in,out] x A parameter used to return the solution.  Notice, any
        !!  input value will be ignored as this routine relies upon the
        !!  search limits in @p lim to provide a starting point.
        !! @param[in] lim A value_pair object defining the search limits.
        !! @param[out] f An optional parameter used to return the function
        !!  residual as computed at @p x.
        !! @param[out] ib An optional output, that if provided, allows the
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can be
        !!  used to retrieve information relating to any errors encountered during
        !!  execution.  If not provided, a default implementation of the errors
        !!  class is used internally to provide error handling.  Possible errors and
        !!  warning messages that may be encountered are as follows.
        !!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
        !!  - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is different
        !!      than the number of variables.
        !!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
        !!      the allowed number of iterations.
        !!
        !! @par Example
        !! @code{.f90}
        !! program example
        !!     use iso_fortran_env
        !!     use nonlin_core
        !!     use nonlin_solve
        !!     implicit none
        !!
        !!     ! Local variables
        !!     type(fcn1var_helper) :: obj
        !!     procedure(fcn1var), pointer :: fcn
        !!     type(newton_1var_solver) :: solver
        !!     real(real64) :: x, f
        !!     type(value_pair) :: limits
        !!     type(iteration_behavior) :: tracking
        !!
        !!     ! Define the search limits
        !!     limits%x1 = 1.5d0
        !!     limits%x2 = 5.0d0
        !!
        !!     ! Establish the function
        !!     fcn => fcn1
        !!     call obj%set_fcn(fcn)
        !!
        !!     ! Solve the equation
        !!     call solver%solve(obj, x, limits, f, ib = tracking)
        !!
        !!     ! Print the output and the residual
        !!     print '(AF7.5)', "The solution: ", x
        !!     print '(AE10.3)', "The residual: ", f
        !!     print '(AI0)', "Iterations: ", tracking%iter_count
        !!     print '(AI0)', "Function Evaluations: ", tracking%fcn_count
        !!     print '(AI0)', "Derivative Evaluations: ", tracking%jacobian_count
        !!     print '(AL1)', "Converge on Function Value: ", tracking%converge_on_fcn
        !!     print '(AL1)', "Converge on Change in Variable: ", tracking%converge_on_chng
        !!     print '(AL1)', "Converge on Derivative: ", tracking%converge_on_zero_diff
        !!
        !! contains
        !!     ! The function:
        !!     ! f(x) = sin(x) / x, solution: x = n * pi for n = 1, 2, 3, ...
        !!     function fcn1(x) result(f)
        !!         real(real64), intent(in) :: x
        !!         real(real64) :: f
        !!         f = sin(x) / x
        !!     end function
        !! end program
        !! @endcode
        !! The above program produces the following output.
        !! @code{.txt}
        !! The solution: 3.14159
        !! The residual:  0.314E-11
        !! Iterations: 3
        !! Function Evaluations: 7
        !! Derivative Evaluations: 4
        !! Converge on Function Value: T
        !! Converge on Change in Variable: F
        !! Converge on Derivative: F
        !! @endcode
        procedure, public :: solve => newt1var_solve
    end type
! ------------------------------------------------------------------------------
    interface
        module subroutine newt1var_solve(this, fcn, x, lim, f, ib, err)
            class(newton_1var_solver), intent(inout) :: this
            class(fcn1var_helper), intent(in) :: fcn
            real(real64), intent(inout) :: x
            type(value_pair), intent(in) :: lim
            real(real64), intent(out), optional :: f
            type(iteration_behavior), optional :: ib
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

contains
! ******************************************************************************
! GENERAL ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Tests for convergence.
    !!
    !! @param[in] x The current solution estimate.
    !! @param[in] xo The previous solution estimate.
    !! @param[in] f The current residual based upon @p x.
    !! @param[in] g The current estimate of the gradient vector at @p x.
    !! @param[in] lg Set to true if the gradient slope check should be
    !!  performed; else, false.
    !! @param[in] xtol The tolerance on the change in variable.
    !! @param[in] ftol The tolerance on the residual.
    !! @param[in] gtol The tolerance on the slope of the gradient.
    !! @param[out] c True if the solution converged on either the residual or
    !!  change in variable.
    !! @param[out] cx True if convergence occurred due to change in variable.
    !! @param[out] cf True if convergence occurred due to residual.
    !! @param[out] cg True if convergence occured due to slope of the gradient.
    !! @param[out] xnorm The largest magnitude component of the scaled change
    !!  in variable vector.
    !! @param[out] fnorm The largest magnitude residual component
    subroutine test_convergence(x, xo, f, g, lg, xtol, ftol, gtol, c, cx, cf, &
            cg, xnorm, fnorm)
        ! Arguments
        real(real64), intent(in), dimension(:) :: x, xo, f, g
        real(real64), intent(in) :: xtol, ftol, gtol
        logical, intent(in) :: lg
        logical, intent(out) :: c, cx, cf, cg
        real(real64), intent(out) :: xnorm, fnorm

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: half = 0.5d0

        ! Local Variables
        integer(int32) :: i, nvar, neqn
        real(real64) :: test, dxmax, fc, den

        ! Initialization
        nvar = size(x)
        neqn = size(f)
        cx = .false.
        cf = .false.
        cg = .false.
        c = .false.
        fc = half * dot_product(f, f)
        fnorm = zero
        xnorm = zero

        ! Test for convergence on residual
        do i = 1, neqn
            fnorm = max(abs(f(i)), fnorm)
        end do
        if (fnorm < ftol) then
            cf = .true.
            c = .true.
            return
        end if

        ! Test the change in solution
        do i = 1, nvar
            test = abs(x(i) - xo(i)) / max(abs(x(i)), one)
            xnorm = max(test, xnorm)
        end do
        if (xnorm < xtol) then
            cx = .true.
            c = .true.
            return
        end if

        ! Test for zero gradient slope - do not set convergence flag
        if (lg) then
            test = zero
            den = max(fc, half * nvar)
            do i = 1, nvar
                dxmax = abs(g(i)) * max(abs(x(i)), one) / den
                test = max(test, dxmax)
            end do
            if (test < gtol) then
                cg = .true.
            end if
        end if
    end subroutine

! ------------------------------------------------------------------------------
end module
