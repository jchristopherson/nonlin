! nonlin_linesearch.f90

module nonlin_linesearch
    use linalg_constants, only : dp, i32
    use nonlin_types, only : vecfcn_helper
    use ferror, only : errors
    implicit none
    private
    public :: line_search

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief Defines a type capable of performing an inexact, backtracking line
    !! search based on the Armijo-Goldstein condition to find a point as far 
    !! along the specified direction vector that is usable for unconstrained
    !! minimization problems.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Backtracking_line_search)
    !! - [Oxfford Lecture Notes]
    !!      (https://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf)
    !! - [Northwestern University - Line Search]
    !!      (https://optimization.mccormick.northwestern.edu/index.php/Line_search_methods)
    !! - [Northwestern University - Trust Region Methods]
    !!      (https://optimization.mccormick.northwestern.edu/index.php/Trust-region_methods)
    !! - [Wolfram](https://reference.wolfram.com/language/tutorial/UnconstrainedOptimizationLineSearchMethods.html)
    !! - [Numerical Recipes](http://numerical.recipes/)
    type line_search
        !> The maximum number of function evaluations allowed during a single
        !! line search.
        integer(i32) :: m_maxEval = 100
        !> Defines the scaling of the product of the gradient and direction 
        !! vectors as part of the Armijo-Goldstein condition such that
        !! F(X + LAMBDA * P) <= F(X) + LAMBDA * ALPHA * P**T * G, where P is the
        !! search direction vector, G is the gradient vector, and LAMBDA is the
        !! scaling factor.  The parameter must exist on the set (0, 1).  A value
        !! of 1e-4 is typically a good starting point.
        real(dp) :: m_alpha = 1.0d-4
        !> Defines a minimum factor X used to determine a minimum value LAMBDA
        !! such that MIN(LAMBDA) = X * LAMBDA, where LAMBDA defines the distance
        !! along the line search direction assuming a value of one means the
        !! full length of the direction vector is traversed.  As such, the value
        !! must exist on the set (0, 1); however, for practical considerations,
        !! the minimum value should be limited to 0.1 such that the value must
        !! exist on the set [0.1, 1).
        real(dp) :: m_factor = 0.1d0
    contains
        !> @brief Gets the maximum number of function evaluations allowed during
        !! a single line search.
        procedure, public :: get_max_fcn_evals => ls_get_max_eval
        !> @brief Sets the maximum number of function evaluations allowed during
        !! a single line search.
        procedure, public :: set_max_fcn_evals => ls_set_max_eval
    end type

contains
! ******************************************************************************
! LINE SEARCH MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Gets the maximum number of function evaluations allowed during a
    !! single line search.
    !!
    !! @param[in] this The line_search object.
    !! @return The maximum number of function evaluations.
    pure function ls_get_max_eval(this) result(n)
        class(line_search), intent(in) :: this
        integer(i32) :: n
        n = this%m_maxEval
    end function

! --------------------
    !> @brief Sets the maximum number of function evaluations allowed during a
    !! single line search.
    !!
    !! @param[in,out] this The line_search object.
    !! @param[in] x The maximum number of function evaluations.
    subroutine ls_set_max_eval(this, x)
        class(line_search), intent(inout) :: this
        integer(i32), intent(in) :: x
        this%m_maxEval = x
    end subroutine

! ------------------------------------------------------------------------------

! --------------------

! ------------------------------------------------------------------------------

! --------------------

! ------------------------------------------------------------------------------
    !
    subroutine ls_search(this, fcn, xold, grad, dir, x, fvec, fold, fx, err)
        ! Arguments
        class(line_search), intent(in) :: this
        class(vecfcn_helper), intent(in) :: fcn
        real(dp), intent(in), dimension(:) :: xold, dir, grad
        real(dp), intent(out), dimension(:) :: x, fvec
        real(dp), intent(in), optional :: fold
        real(dp), intent(out), optional :: fx
        class(errors), intent(in), optional, target :: err

        ! Parameters

        ! Local Variables
        logical :: xcnvrg, fcnvrg
        integer(i32) :: i, m, n, neval, niter, flag, maxeval
        real(dp) :: a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2, &
            slope, temp, test, tmplam, alpha, tolx, lambdamin, f

        ! Initialization
        xcnvrg = .false.
        fcnvrg = .false.
        neval = 0
        niter = 0
        
    end subroutine

! ------------------------------------------------------------------------------
end module
