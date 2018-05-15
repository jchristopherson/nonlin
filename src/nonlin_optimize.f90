! nonlin_optimize.f90

! REF:
! http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#Nelder-Mead_Simplex
! http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
! https://scicomp.stackexchange.com/questions/14787/fortran-library-for-minimization-or-maximization-of-functions-optimization-prob
! http://ab-initio.mit.edu/wiki/index.php/NLopt

!> @brief \b nonlin_optimize
!!
!! @par Purpose
!! To provide various optimization routines.
module nonlin_optimize
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use linalg_constants, only : LA_MATRIX_FORMAT_ERROR
    use ferror, only : errors
    use nonlin_linesearch, only : line_search, limit_search_vector
    use nonlin_constants
    use nonlin_core
    use linalg_core, only : rank1_update, tri_mtx_mult, cholesky_rank1_update, &
        cholesky_rank1_downdate, solve_cholesky
    implicit none
    private
    public :: nelder_mead
    public :: line_search_optimizer
    public :: bfgs

! ******************************************************************************
! NONLIN_OPTIMIZE_NELDER_MEAD.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a solver based upon Nelder and Mead's simplex algorithm
    !! for minimization of functions of multiple variables.
    type, extends(equation_optimizer) :: nelder_mead
        private
        !> The simplex vertices.
        real(real64), allocatable, dimension(:,:) :: m_simplex
        !> A scaling parameter used to define the size of the simplex in each
        !! coordinate direction.
        real(real64) :: m_initSize = 1.0d0
    contains
        !> @brief Utilizes the Nelder-Mead simplex method for finding a minimum
        !! value of the specified function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine solve(class(nelder_mead) this, class(fcnnvar_helper) fcn, real(real64) x(:), real(real64) fout, optional type(iteration_behavior) ib, optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The nelder_mead object.
        !! @param[in] fcn The fcnnvar_helper object containing the equation to
        !!  optimize.
        !! @param[in,out] x On input, the initial guess at the optimal point.
        !!  On output, the updated optimal point estimate.
        !! @param[out] fout An optional output, that if provided, returns the
        !!  value of the function at @p x.
        !! @param[out] ib An optional output, that if provided, allows the
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can be
        !!  used to retrieve information relating to any errors encountered during
        !!  execution.  If not provided, a default implementation of the errors
        !!  class is used internally to provide error handling.  Possible errors and
        !!  warning messages that may be encountered are as follows.
        !!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
        !!  - NL_INVALID_INPUT_ERROR: Occurs if @p x is not appropriately sized for
        !!      the problem as defined in @p fcn.
        !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
        !!      the allowed number of iterations.
        !!
        !! @par Example
        !! The following example illustrates how to find the minimum of Rosenbrock's
        !! function using this Nelder-Mead solver.
        !! @code{.f90}
        !! program example
        !!     use nonlin_optimize, only : nelder_mead
        !!     use nonlin_core, only : fcnnvar, fcnnvar_helper, iteration_behavior
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     type(nelder_mead) :: solver
        !!     type(fcnnvar_helper) :: obj
        !!     procedure(fcnnvar), pointer :: fcn
        !!     real(real64) :: x(2), fout
        !!     type(iteration_behavior) :: ib
        !!
        !!     ! Initialization
        !!     fcn => rosenbrock
        !!     call obj%set_fcn(fcn, 2)
        !!
        !!     ! Define an initial guess - the solution is (1, 1)
        !!     call random_number(x)
        !!
        !!     ! Call the solver
        !!     call solver%solve(obj, x, fout, ib)
        !!
        !!     ! Display the output
        !!     print '(AF8.5AF8.5A)', "Rosenbrock Minimum: (", x(1), ", ", x(2), ")"
        !!     print '(AE9.3)', "Function Value: ", fout
        !!     print '(AI0)', "Iterations: ", ib%iter_count
        !!     print '(AI0)', "Function Evaluations: ", ib%fcn_count
        !! contains
        !!     ! Rosenbrock's Function
        !!     function rosenbrock(x) result(f)
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64) :: f
        !!         f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
        !!     end function
        !! end
        !! @endcode
        !! The above program yields the following output.
        !! @code{.txt}
        !! Rosenbrock Minimum: ( 1.00000,  1.00000)
        !! Function Value: 0.264E-12
        !! Iterations: 59
        !! Function Evaluations: 112
        !! @endcode
        !!
        !! @par Remarks
        !! The implementation of the Nelder-Mead algorithm presented here is a
        !! slight modification of the original work of Nelder and Mead.  The
        !! Numerical Recipes implementation is also quite similar.  In fact, the
        !! Numerical Recipes section relating to reflection, contraction, etc.
        !! is leveraged for this implemetation.
        !!
        !! @par See Also
        !!  - Nelder, John A.; R. Mead (1965). "A simplex method for function
        !!      minimization". Computer Journal. 7: 308–313.
        !!  - [Gao, Fuchang, Han, Lixing (2010). "Implementing the Nelder-Mead
        !!      simplex algorithm with adaptive parameters."]
        !!      (http://www.webpages.uidaho.edu/~fuchang/res/ANMS.pdf)
        !!  - [Wikipedia](https://en.wikipedia.org/wiki/Nelder–Mead_method)
        !!  - [Numerical Recipes](http://numerical.recipes/)
        procedure, public :: solve => nm_solve
        !> @brief Gets an N-by-(N+1) matrix containing the current simplex.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64)(:,:) function get_simplex(class(nelder_mead) this)
        !! @endcode
        !!
        !! @param[in] this The nelder_mead object.
        !! @return The N-by-(N+1) matrix containing the simplex.  Each vertex of the
        !!  simplex is stored as its own column of this matrix.
        procedure, public :: get_simplex => nm_get_simplex
        !> @brief Sets an N-by-(N+1) matrix as the current simplex.  Notice, if
        !! this matrix is different in size from the problem dimensionallity,
        !! the Nelder-Mead routine will replace it with an appropriately sized
        !! matrix.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_simplex(class(nelder_mead) this, real(real64) x(:,:))
        !! @endcode
        !!
        !! @param[in,out] this The nelder_mead object.
        !! @param[in] x The simplex matrix.  Each column of the matrix must contain
        !!  the coordinates of each vertex of the simplex.
        !!
        !! @par example
        !! @code{.f90}
        !! @endcode
        procedure, public :: set_simplex => nm_set_simplex
        !> @brief Gets the size of the initial simplex that will be utilized by
        !! the Nelder-Mead algorithm in the event that the user does not supply
        !! a simplex geometry, or if the user supplies an invalid simplex
        !! geometry.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) get_initial_size(class(nelder_mead) this)
        !! @endcode
        !!
        !! @param[in] this The nelder_mead object.
        !! @return The size of the simplex (length of an edge).
        procedure, public :: get_initial_size => nm_get_size
        !> @brief Sets the size of the initial simplex that will be utilized by
        !! the Nelder-Mead algorithm in the event that the user does not supply
        !! a simplex geometry, or if the user supplies an invalid simplex
        !! geometry.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_initial_size(class(nelder_mead) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The nelder_mead object.
        !! @param[in] x The size of the simplex (length of an edge).
        !!
        !! @par example
        !! @code{.f90}
        !! @endcode
        procedure, public :: set_initial_size => nm_set_size

        !> @brief Extrapolates by the specified factor through the simplex across
        !! from the largest point.  If the extrapolation results in a better
        !! estimate, the current high point is replaced with the new estimate.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function extrapolate(class(nelder_mead) this, class(fcnnvar_helper) fcn, real(real64) y(:), real(real64) pcent(:), integer(int32) ihi, real(real64) fac, integer(int32) neval, real(real64) work(:))
        !! @endcode
        !!
        !! @param[in,out] this The nelder_mead object.
        !! @param[in] fcn The function to evaluate.
        !! @param[in,out] y An array containing the function values at each vertex.
        !! @param[in,out] pcent An array containing the centroid of vertex position
        !!  information.
        !! @param[in] ihi The index of the largest magnitude vertex.
        !! @param[in,out] neval The number of function evaluations.
        !! @param[out] work An N-element workspace array where N is the number of
        !!  dimensions of the problem.
        !! @return The new function estimate.
        procedure, private :: extrapolate => nm_extrapolate
    end type

! ------------------------------------------------------------------------------
    interface
        module subroutine nm_solve(this, fcn, x, fout, ib, err)
            class(nelder_mead), intent(inout) :: this
            class(fcnnvar_helper), intent(in) :: fcn
            real(real64), intent(inout), dimension(:) :: x
            real(real64), intent(out), optional :: fout
            type(iteration_behavior), optional :: ib
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module function nm_extrapolate(this, fcn, y, pcent, ihi, fac, neval, &
                work) result(ytry)
            class(nelder_mead), intent(inout) :: this
            class(fcnnvar_helper), intent(in) :: fcn
            real(real64), intent(inout), dimension(:) :: y, pcent
            integer(int32), intent(in) :: ihi
            real(real64), intent(in) :: fac
            integer(int32), intent(inout) :: neval
            real(real64), intent(out), dimension(:) :: work
            real(real64) :: ytry
        end function

        pure module function nm_get_simplex(this) result(p)
            class(nelder_mead), intent(in) :: this
            real(real64), allocatable, dimension(:,:) :: p
        end function

        module subroutine nm_set_simplex(this, x)
            class(nelder_mead), intent(inout) :: this
            real(real64), dimension(:,:) :: x
        end subroutine

        pure module function nm_get_size(this) result(x)
            class(nelder_mead), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine nm_set_size(this, x)
            class(nelder_mead), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine
    end interface

! ******************************************************************************
! NONLIN_OPTIMIZE_LINE_SEARCH.F90
! ------------------------------------------------------------------------------
    !> @brief A class describing equation optimizers that use a line search
    !! algorithm to improve convergence behavior.
    type, abstract, extends(equation_optimizer) :: line_search_optimizer
        private
        !> The line search object.
        class(line_search), allocatable :: m_lineSearch
        !> Set to true if a line search should be used regardless of the status
        !! of m_lineSearch
        logical :: m_useLineSearch = .true.
        !> The convergence criteria on change in variable
        real(real64) :: m_xtol = 1.0d-12
    contains
        !> @brief Gets the line search module.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine get_line_search(class(line_search_optimizer) this, allocatable class(line_search) ls)
        !! @endcode
        !!
        !! @param[in] this The line_search_optimizer object.
        !! @param[out] ls The line_search object.
        procedure, public :: get_line_search => lso_get_line_search
        !> @brief Sets the line search module.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_line_search(class(line_search_optimizer) this, class(line_search) ls)
        !! @endcode
        !!
        !! @param[in,out] this The line_search_optimizer object.
        !! @param[in] ls The line_search object.
        procedure, public :: set_line_search => lso_set_line_search
        !> @brief Establishes a default line_search object for the line search
        !! module.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_default_line_search(class(line_search_optimizer) this)
        !! @endcode
        !!
        !! @param[in,out] this The line_search_optimizer object.
        procedure, public :: set_default_line_search => lso_set_default
        !> @brief Tests to see if a line search module is defined.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function is_line_search_defined(class(line_search_optimizer) this)
        !! @endcode
        !!
        !! @param[in] this The line_search_optimizer object.
        !! @return Returns true if a module is defined; else, false.
        procedure, public :: is_line_search_defined => &
            lso_is_line_search_defined
        !> @brief Gets a value determining if a line-search should be employed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical get_use_line_search(class(line_search_optimizer) this)
        !! @endcode
        !!
        !! @param[in] this The line_search_optimizer object.
        !! @return Returns true if a line search should be used; else, false.
        procedure, public :: get_use_line_search => lso_get_use_search
        !> @brief Sets a value determining if a line-search should be employed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_use_line_search(class(line_search_optimizer) this, logical x)
        !! @endcode
        !!
        !! @param[in,out] this The line_search_optimizer object.
        !! @param[in] x Set to true if a line search should be used; else, false.
        procedure, public :: set_use_line_search => lso_set_use_search
        !> @brief Gets the convergence on change in variable tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_var_tolerance(class(line_search_optimizer) this)
        !! @endcode
        !!
        !! @param[in] this The line_search_optimizer object.
        !! @return The tolerance value.
        procedure, public :: get_var_tolerance => lso_get_var_tol
        !> @brief Sets the convergence on change in variable tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_var_tolerance(class(line_search_optimizer) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The line_search_optimizer object.
        !! @param[in] x The tolerance value.
        procedure, public :: set_var_tolerance => lso_set_var_tol
    end type

! ------------------------------------------------------------------------------
    interface
        module subroutine lso_get_line_search(this, ls)
            class(line_search_optimizer), intent(in) :: this
            class(line_search), intent(out), allocatable :: ls
        end subroutine

        module subroutine lso_set_line_search(this, ls)
            class(line_search_optimizer), intent(inout) :: this
            class(line_search), intent(in) :: ls
        end subroutine

        module subroutine lso_set_default(this)
            class(line_search_optimizer), intent(inout) :: this
            type(line_search) :: ls
        end subroutine

        pure module function lso_is_line_search_defined(this) result(x)
            class(line_search_optimizer), intent(in) :: this
            logical :: x
        end function

        pure module function lso_get_use_search(this) result(x)
            class(line_search_optimizer), intent(in) :: this
            logical :: x
        end function

        module subroutine lso_set_use_search(this, x)
            class(line_search_optimizer), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine

        pure module function lso_get_var_tol(this) result(x)
            class(line_search_optimizer), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine lso_set_var_tol(this, x)
            class(line_search_optimizer), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines a Broyden–Fletcher–Goldfarb–Shanno (BFGS) solver for
    !! minimization of functions of multiple variables.
    type, extends(line_search_optimizer) :: bfgs
    contains
        !> @brief Utilizes the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
        !! for finding a minimum value of the specified function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine solve(class(bfgs) this, class(fcnnvar_helper) fcn, real(real64) x(:), real(real64) fout, optional type(iteration_behavior) ib, optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The bfgs_mead object.
        !! @param[in] fcn The fcnnvar_helper object containing the equation to
        !!  optimize.
        !! @param[in,out] x On input, the initial guess at the optimal point.
        !!  On output, the updated optimal point estimate.
        !! @param[out] fout An optional output, that if provided, returns the
        !!  value of the function at @p x.
        !! @param[out] ib An optional output, that if provided, allows the
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can be
        !!  used to retrieve information relating to any errors encountered during
        !!  execution.  If not provided, a default implementation of the errors
        !!  class is used internally to provide error handling.  Possible errors and
        !!  warning messages that may be encountered are as follows.
        !!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
        !!  - NL_INVALID_INPUT_ERROR: Occurs if @p x is not appropriately sized for
        !!      the problem as defined in @p fcn.
        !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
        !!      available.
        !!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
        !!      the allowed number of iterations.
        !!
        !! @par Example
        !! The following example illustrates how to find the minimum of Rosenbrock's
        !! function using this BFGS solver.
        !! @code{.f90}
        !! program example
        !!     use nonlin_optimize, only : bfgs
        !!     use nonlin_core, only : fcnnvar, fcnnvar_helper, iteration_behavior
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     type(bfgs) :: solver
        !!     type(fcnnvar_helper) :: obj
        !!     procedure(fcnnvar), pointer :: fcn
        !!     real(real64) :: x(2), fout
        !!     type(iteration_behavior) :: ib
        !!
        !!     ! Initialization
        !!     fcn => rosenbrock
        !!     call obj%set_fcn(fcn, 2)
        !!
        !!     ! Define an initial guess - the solution is (1, 1)
        !!     call random_number(x)
        !!
        !!     ! Call the solver
        !!     call solver%solve(obj, x, fout, ib)
        !!
        !!     ! Display the output
        !!     print '(AF8.5AF8.5A)', "Rosenbrock Minimum: (", x(1), ", ", x(2), ")"
        !!     print '(AE9.3)', "Function Value: ", fout
        !!     print '(AI0)', "Iterations: ", ib%iter_count
        !!     print '(AI0)', "Function Evaluations: ", ib%fcn_count
        !! contains
        !!     ! Rosenbrock's Function
        !!     function rosenbrock(x) result(f)
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64) :: f
        !!         f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
        !!     end function
        !! end
        !! @endcode
        !! The above program yields the following output:
        !! @code{.txt}
        !! Rosenbrock Minimum: ( 1.00000,  0.99999)
        !! Function Value: 0.200E-10
        !! Iterations: 47
        !! Function Evaluations: 70
        !! @endcode
        !!
        !! @par See Also
        !! - [Wikipedia - BFGS Methods](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm)
        !! - [Wikipedia - Quasi-Newton Methods](https://en.wikipedia.org/wiki/Quasi-Newton_method)
        !! - [minFunc](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html)
        procedure, public :: solve => bfgs_solve
    end type

! ------------------------------------------------------------------------------
    interface
        module subroutine bfgs_solve(this, fcn, x, fout, ib, err)
            class(bfgs), intent(inout) :: this
            class(fcnnvar_helper), intent(in) :: fcn
            real(real64), intent(inout), dimension(:) :: x
            real(real64), intent(out), optional :: fout
            type(iteration_behavior), optional :: ib
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

end module
