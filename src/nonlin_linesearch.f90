! nonlin_linesearch.f90

module nonlin_linesearch
    use linalg_constants, only : dp, i32
    use nonlin_types
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
        !> @brief Gets the scaling of the product of the gradient and direction 
        !! vectors.
        procedure, public :: get_scaling_factor => ls_get_scale
        !> @brief Sets the scaling of the product of the gradient and direction 
        !! vectors.
        procedure, public :: set_scaling_factor => ls_set_scale
        !> @brief Gets a distance factor defining the minimum distance along the 
        !! search direction vector is practical.
        procedure, public :: get_distance_factor => ls_get_dist
        !> @brief Sets a distance factor defining the minimum distance along the 
        !! search direction vector is practical.
        procedure, public :: set_distance_factor => ls_set_dist

        procedure, public :: search => ls_search
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
    !> @brief Gets the scaling of the product of the gradient and direction 
    !! vectors as part of the Armijo-Goldstein condition such that
    !! F(X + LAMBDA * P) <= F(X) + LAMBDA * ALPHA * P**T * G, where P is the
    !! search direction vector, G is the gradient vector, and LAMBDA is the
    !! scaling factor.
    !!
    !! @param[in] this The line_search object.
    !! @return The scaling factor.
    pure function ls_get_scale(this) result(x)
        class(line_search), intent(in) :: this
        real(dp) :: x
        x = this%m_alpha
    end function

! --------------------
    !> @brief Sets the scaling of the product of the gradient and direction 
    !! vectors as part of the Armijo-Goldstein condition such that
    !! F(X + LAMBDA * P) <= F(X) + LAMBDA * ALPHA * P**T * G, where P is the
    !! search direction vector, G is the gradient vector, and LAMBDA is the
    !! scaling factor.
    !!
    !! @param[in,out] this The line_search object.
    !! @param[in] x The scaling factor.
    subroutine ls_set_scale(this, x)
        class(line_search), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_alpha = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets a distance factor defining the minimum distance along the 
    !! search direction vector is practical.
    !!
    !! @param[in] this The line_search object.
    !! @return The distance factor.  A value of 1 indicates the full length
    !!  of the vector.
    pure function ls_get_dist(this) result(x)
        class(line_search), intent(in) :: this
        real(dp) :: x
        x = this%m_factor
    end function

! --------------------
    !> @brief Sets a distance factor defining the minimum distance along the 
    !! search direction vector is practical.
    !!
    !! @param[in,out] this The line_search object.
    !! @param[in] x The distance factor.  A value of 1 indicates the full length
    !!  of the vector.  Notice, this value is restricted to lie in the set
    !!  [0.1, 1.0)
    subroutine ls_set_dist(this, x)
        class(line_search), intent(inout) :: this
        real(dp), intent(in) :: x
        if (x <= 0.0d0) then
            this%m_factor = 0.1d0
        else if (x >= 1.0d0) then
            this%m_factor = 0.99d0
        else
            this%m_factor = x
        end if
    end subroutine

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
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: p5 = 0.5d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: two = 2.0d0
        real(dp), parameter :: three = 3.0d0

        ! Local Variables
        logical :: xcnvrg, fcnvrg
        integer(i32) :: i, m, n, neval, niter, flag, maxeval
        real(dp) :: a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2, &
            slope, temp, test, tmplam, alpha, tolx, lambdamin, f, fo
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        xcnvrg = .false.
        fcnvrg = .false.
        neval = 0
        niter = 0
        m = fcn%get_function_count()
        n = fcn%get_variable_count()
        tolx = two * epsilon(tolx)
        alpha = this%m_alpha
        lambdamin = this%m_factor
        maxeval = this%m_maxEval
        if (present(fx)) fx = zero
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        flag = 0
        if (size(xold) /= n) then
            flag = 3
        else if (size(grad) /= n) then
            flag = 4
        else if (size(dir) /= n) then
            flag = 5
        else if (size(x) /= n) then
            flag = 6
        else if (size(fvec) /= m) then
            flag = 7
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("ls_search", trim(errmsg), &
                NL_ARRAY_SIZE_ERROR)
            return
        end if

        ! Compute 1/2 F * F (* = dot product) if not provided
        if (present(fold)) then
            fo = fold
        else
            ! Evaluate the function, and compute the dot product
            call fcn%fcn(xold, fvec)
            fo = p5 * dot_product(fvec, fvec)
            neval = neval + 1
        end if

        ! Compute the slope parameter
        slope = dot_product(grad, dir)
        if (slope >= zero) then
            ! ERROR: The slope should not be pointing uphill - invalid direction
            call errmgr%report_error("ls_search", "The search direction " // &
                "vector appears to be pointing in an uphill direction - " // &
                "away from a minimum.", NL_INVALID_OPERATION_ERROR)
            return
        end if

        ! Compute the minimum lambda value (length along the search direction)
        test = zero
        do i = 1, n
            temp = abs(dir(i)) / max(abs(xold(i)), one)
            if (temp > test) test = temp
        end do
        alamin = tolx / test
        alam = one

        ! Iteration Loop
        do
            ! Step along the specified direction by the amount ALAM
            x = xold + alam * dir
            call fcn%fcn(x, fvec)
            f = p5 * dot_product(fvec, fvec)
            neval = neval + 1
            niter = niter + 1

            ! Check the step
            if (alam < alamin) then
                ! The change in X between steps is sufficiently small to 
                ! consider the iteration converged
                x = xold
                xcnvrg = .true.
                exit
            else if (f <= fo + alpha * alam * slope) then
                ! The function has converged
                fcnvrg = .true.
                exit
            else
                ! Convergence has not yet occurred, continue backtracking
                if (niter == 1) then
                    ! Use the quadratic function approximation
                    tmplam = -slope / (two * (f - fo - slope))
                else
                    ! Use the cubic function approximation
                    rhs1 = f - fo - alam * slope
                    rhs2 = f2 - fo - alam2 * slope
                    a = (rhs1 / alam**2 - rhs2 / alam2**2) / (alam - alam2)
                    b = (-alam2 * rhs1 / alam**2 + alam * rhs2 / alam2**2) / &
                        (alam - alam2)
                    if (a == zero) then
                        tmplam = -slope / (two * b)
                    else
                        disc = b**2 - three * a * slope
                        if (disc < zero) then
                            tmplam = p5 * alam
                        else if (b <= zero) then
                            tmplam = (-b + sqrt(disc)) / (three * a)
                        else
                            tmplam = -slope / (b + sqrt(disc))
                        end if
                    end if
                    if (tmplam > p5 * alam) tmplam = p5 * alam
                end if
            end if
            
            ! Set up parameters for the cubic model as we've already been 
            ! through once with the quadratic model without success.
            alam2 = alam
            f2 = f2
            alam = max(tmplam, lambdamin * alam)

            ! Ensure we haven't performed too many function evaluations
            if (neval >= maxeval) then
                ! ERROR: Too many function evaluations
                write(errmsg, '(AI0A)') "The line search failed to " // &
                    "converge.  Function evaluations performed: ", neval, "."
                call errmgr%report_error("ls_search", errmsg, &
                    NL_CONVERGENCE_ERROR)
                return
            end if
        end do
    end subroutine

! ------------------------------------------------------------------------------
end module
