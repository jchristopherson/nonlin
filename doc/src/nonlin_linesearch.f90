! REFERENCES
! - https://scicomp.stackexchange.com/questions/26330/backtracking-armijo-line-search-algorithm
! - https://ctk.math.ncsu.edu/
module nonlin_linesearch
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use nonlin_error_handling
    use nonlin_types
    use nonlin_multi_eqn_mult_var
    use nonlin_multi_var
    use ferror, only : errors
    implicit none
    private
    public :: line_search
    public :: limit_search_vector

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    type line_search
        !! Defines a type capable of performing an inexact, backtracking line
        !! search to find a point as far along the specified direction vector 
        !! that is usable for unconstrained minimization problems.
        !!
        !! See Also:
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Backtracking_line_search)
        !!
        !! - [Oxfford Lecture Notes](https://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf)
        !!
        !! - [Wolfram](https://reference.wolfram.com/language/tutorial/UnconstrainedOptimizationLineSearchMethods.html)
        !!
        !! - [Numerical Recipes](http://numerical.recipes/)
        integer(int32), private :: m_maxEval = 100
            !! The maximum number of function evaluations allowed during a 
            !! single line search.
        real(real64), private :: m_alpha = 1.0d-4
            !! Defines the scaling of the product of the gradient and direction
            !! vectors such that \( f(\vec{x} + \lambda \vec{p}) <= 
            !! f(\vec{x}) + \lambda \alpha \vec{p}^{T} \vec{g} \) where 
            !! \(\vec{p}\) is the search direction vector, \(\vec{g}\) is the 
            !! gradient vector, and \(\lambda\) is the scaling factor.  The 
            !! parameter must exist on the set (0, 1).  A value of 1e-4 is 
            !! typically a good starting point.
        real(real64), private :: m_factor = 0.1d0
            !! Defines a minimum factor X used to determine a minimum value 
            !! \(\lambda\) where \(\lambda\) defines the distance along the 
            !! line search direction assuming a value of one means the full 
            !! length of the direction vector is traversed.  As such, the value 
            !! must exist on the set (0, 1); however, for practical 
            !! considerations, the minimum value should be limited to 0.1 such 
            !! that the value must exist on the set [0.1, 1).
    contains
        procedure, public :: get_max_fcn_evals => ls_get_max_eval
        procedure, public :: set_max_fcn_evals => ls_set_max_eval
        procedure, public :: get_scaling_factor => ls_get_scale
        procedure, public :: set_scaling_factor => ls_set_scale
        procedure, public :: get_distance_factor => ls_get_dist
        procedure, public :: set_distance_factor => ls_set_dist
        generic, public :: search => ls_search_mimo, ls_search_miso

        procedure :: ls_search_mimo
        procedure :: ls_search_miso
    end type

contains
! ******************************************************************************
! LINE_SEARCH MEMBERS
! ------------------------------------------------------------------------------
    pure function ls_get_max_eval(this) result(n)
        !! Gets the maximum number of function evaluations allowed during a
        !! single line search.
        class(line_search), intent(in) :: this
            !! The [[line_search]] object.
        integer(int32) :: n
            !! The maximum number of function evaluations.
        n = this%m_maxEval
    end function

! --------------------
    subroutine ls_set_max_eval(this, x)
        !! Sets the maximum number of function evaluations allowed during a
        !! single line search.
        class(line_search), intent(inout) :: this
            !! The [[line_search]] object.
        integer(int32), intent(in) :: x
            !! The maximum number of function evaluations.
        this%m_maxEval = x
    end subroutine

! ------------------------------------------------------------------------------
    pure function ls_get_scale(this) result(x)
        !! Gets the scaling of the product of the gradient and direction
        !! vectors \(\alpha\) such that \( f(\vec{x} + \lambda \vec{p}) \le 
        !! f(\vec{x}) + \lambda \alpha \vec{p}^{T} \vec{g} \), where \(\vec{p}\) 
        !! is the search direction vector, \(\vec{g}\) is the gradient vector, 
        !! and \(\lambda\) is the scaling factor.
        class(line_search), intent(in) :: this
            !! The [[line_search]] object.
        real(real64) :: x
            !! The scaling factor.
        x = this%m_alpha
    end function

! --------------------
    subroutine ls_set_scale(this, x)
        !! Sets the scaling of the product of the gradient and direction
        !! vectors \(\alpha\) such that \( f(\vec{x} + \lambda \vec{p}) \le 
        !! f(\vec{x}) + \lambda \alpha \vec{p}^{T} \vec{g} \), where \(\vec{p}\) 
        !! is the search direction vector, \(\vec{g}\) is the gradient vector, 
        !! and \(\lambda\) is the scaling factor.
        class(line_search), intent(inout) :: this
            !! The [[line_search]] object.
        real(real64), intent(in) :: x
            !! The scaling factor.
        this%m_alpha = x
    end subroutine

! ------------------------------------------------------------------------------
    pure function ls_get_dist(this) result(x)
        !! Gets a distance factor defining the minimum distance along the
        !! search direction vector is practical.
        class(line_search), intent(in) :: this
            !! The [[line_search]] object.
        real(real64) :: x
            !! The distance factor.  A value of 1 indicates the full length
            !! of the vector.
        x = this%m_factor
    end function

! --------------------
    subroutine ls_set_dist(this, x)
        !! Sets a distance factor defining the minimum distance along the
        !! search direction vector is practical.
        class(line_search), intent(inout) :: this
            !! The [[line_search]] object.
        real(real64), intent(in) :: x
            !! The distance factor.  A value of 1 indicates the full length
            !! of the vector.  Notice, this value is restricted to lie in the 
            !! set [0.1, 1.0).
        if (x <= 0.0d0) then
            this%m_factor = 0.1d0
        else if (x >= 1.0d0) then
            this%m_factor = 0.99d0
        else
            this%m_factor = x
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine ls_search_mimo(this, fcn, xold, grad, dir, x, fvec, fold, fx, &
            ib, err)
        !! Utilizes an inexact, backtracking line search to find a point as
        !! far along the specified direction vector that is usable for 
        !! unconstrained minimization problems.
        class(line_search), intent(in) :: this
            !! The [[line_search]] object.
        class(vecfcn_helper), intent(in) :: fcn
            !! A [[vecfcn_helper]] object containing the system of equations.
        real(real64), intent(in), dimension(:) :: xold
            !! An N-element array defining the initial point, where N is the 
            !! number of variables.
        real(real64), intent(in), dimension(:) :: grad
            !! An N-element array defining the gradient of fcn evaluated at 
            !! xold.
        real(real64), intent(in), dimension(:) :: dir
            !! An N-element array defining the search direction.
        real(real64), intent(out), dimension(:) :: x
            !! An N-element array where the updated solution point will be 
            !! written.
        real(real64), intent(out), dimension(:) :: fvec
            !! An M-element array containing the M equation values evaluated at
            !! x, where M is the number of equations.
        real(real64), intent(in), optional :: fold
            !! An optional input that provides the value resulting from:
            !! \( \frac{1}{2} \vec{f_{old}} \cdot \vec{f_{old}} \).  If not 
            !! provided, fcn is evalauted at xold, and the aforementioned 
            !! relationship is computed. 
        real(real64), intent(out), optional :: fx
            !! The result of the operation: \( \frac{1}{2} \vec{f} \cdot 
            !! \vec{f} \).
        type(iteration_behavior), optional :: ib
            !! An optional output, that if provided, allows the caller to
            !! obtain iteration performance statistics.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: p5 = 0.5d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: two = 2.0d0
        real(real64), parameter :: three = 3.0d0

        ! Local Variables
        logical :: xcnvrg, fcnvrg
        integer(int32) :: i, m, n, neval, niter, flag, maxeval
        real(real64) :: alam, alam1, alamin, f1, slope, temp, test, tmplam, &
            alpha, tolx, lambdamin, f, fo
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        xcnvrg = .false.
        fcnvrg = .false.
        neval = 0
        niter = 0
        m = fcn%get_equation_count()
        n = fcn%get_variable_count()
        tolx = two * epsilon(tolx)
        alpha = this%m_alpha
        lambdamin = this%m_factor
        maxeval = this%m_maxEval
        if (present(fx)) fx = zero
        if (present(ib)) then
            ib%iter_count = niter
            ib%fcn_count = neval
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = .false.
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("ls_search_mimo", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
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
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("ls_search_mimo", trim(errmsg), &
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
            call errmgr%report_error("ls_search_mimo", &
                "The search direction vector appears to be pointing in " // &
                "an uphill direction -  away from a minimum.", &
                NL_DIVERGENT_BEHAVIOR_ERROR)
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
        flag = 0 ! Used to check for convergence errors
        do
            ! Step along the specified direction by the amount ALAM
            x = xold + alam * dir
            call fcn%fcn(x, fvec)
            f = p5 * dot_product(fvec, fvec)
            neval = neval + 1
            niter = niter + 1

            ! Check the step
            if (alam < alamin) then
                ! Either the solution has converged due to negligible change in
                ! the root values, or the line search may have fully
                ! backtracked.  In the event the solution fully backtracked,
                ! we'll issue a warning to inform the user of the potential
                ! issue.
                if (norm2(x - xold) == zero) then
                    ! The line search fully backtracked
                    call errmgr%report_warning("ls_search_mimo", &
                        "The line search appears to have fully " // &
                        "backtracked.  As such, check results carefully, " // &
                        "and/or consider attempting the solve without " // &
                        "the line search.", &
                        NL_CONVERGENCE_ERROR)
                end if
                x = xold
                xcnvrg = .true.
                exit
            else if (f <= fo + alpha * alam * slope) then
                ! The function has converged
                fcnvrg = .true.
                exit
            else
                ! Convergence has not yet occurred, continue backtracking
                tmplam = min_backtrack_search(niter, fo, f, f1, alam, &
                    alam1, slope)
            end if

            ! Set up parameters for the cubic model as we've already been
            ! through once with the quadratic model without success.
            alam1 = alam
            f1 = f
            alam = max(tmplam, lambdamin * alam)

            ! Ensure we haven't performed too many function evaluations
            if (neval >= maxeval) then
                ! ERROR: Too many function evaluations
                flag = 1
                exit
            end if
        end do
        if (present(fx)) fx = f

        ! Report out iteration statistics
        if (present(ib)) then
            ib%iter_count = niter
            ib%fcn_count = neval
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = .false.
        end if

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, 100) "The line search failed to " // &
                "converge.  Function evaluations performed: ", neval, "."
            call errmgr%report_error("ls_search_mimo", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine ls_search_miso(this, fcn, xold, grad, dir, x, fold, fx, &
            ib, err)
        !! Utilizes an inexact, backtracking line search to find a point as far 
        !! along the specified direction vector that is usable for unconstrained
        !! minimization problems.
        class(line_search), intent(in) :: this
            !! The [[line_search]] object.
        class(fcnnvar_helper), intent(in) :: fcn
            !! A [[fcnnvar_helper]] object containing the system of equations.
        real(real64), intent(in), dimension(:) :: xold
            !! An N-element array defining the initial point, where N is the 
            !! number of variables.
        real(real64), intent(in), dimension(:) :: grad
            !! An N-element array defining the gradient of fcn evaluated at 
            !! xold.
        real(real64), intent(in), dimension(:) :: dir
            !! An N-element array defining the search direction.
        real(real64), intent(out), dimension(:) :: x
            !! An N-element array where the updated solution point will be 
            !! written.
        real(real64), intent(in), optional :: fold
            !! An optional input that provides the function value at xold.  If 
            !! not provided, fcn is evalauted at xold.
        real(real64), intent(out), optional :: fx
            !! The value of the function as evaluated at x.
        type(iteration_behavior), optional :: ib
            !! An optional output, that if provided, allows the caller to
            !! obtain iteration performance statistics.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: p5 = 0.5d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: two = 2.0d0
        real(real64), parameter :: three = 3.0d0

        ! Local Variables
        logical :: xcnvrg, fcnvrg
        integer(int32) :: i, n, neval, niter, flag, maxeval
        real(real64) :: alam, alam1, alamin, f1, slope, temp, test, tmplam, &
            alpha, tolx, lambdamin, f, fo
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 128) :: errmsg

        ! Initialization
        xcnvrg = .false.
        fcnvrg = .false.
        neval = 0
        niter = 0
        n = fcn%get_variable_count()
        tolx = two * epsilon(tolx)
        alpha = this%m_alpha
        lambdamin = this%m_factor
        maxeval = this%m_maxEval
        if (present(fx)) fx = zero
        if (present(ib)) then
            ib%iter_count = niter
            ib%fcn_count = neval
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = .false.
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("ls_search_miso", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        flag = 0
        if (size(xold) /= n) then
            flag = 3
        else if (size(grad) /= n) then
            flag = 4
        else if (size(dir) /= n) then
            flag = 5
        else if (size(x) /= n) then
            flag = 6
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, 100) "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("ls_search_miso", trim(errmsg), &
                NL_ARRAY_SIZE_ERROR)
            return
        end if

        ! Establish the "old" function value
        if (present(fold)) then
            fo = fold
        else
            ! Evaluate the function
            fo = fcn%fcn(xold)
            neval = neval + 1
        end if

        ! Compute the slope parameter
        slope = dot_product(grad, dir)
        if (slope >= zero) then
            ! ERROR: The slope should not be pointing uphill - invalid direction
            call errmgr%report_error("ls_search_miso", &
                "The search direction vector appears to be pointing in " // &
                "an uphill direction -  away from a minimum.", &
                NL_DIVERGENT_BEHAVIOR_ERROR)
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
        flag = 0 ! Used to check for convergence errors
        do
            ! Step along the specified direction by the amount ALAM
            x = xold + alam * dir
            f = fcn%fcn(x)
            neval = neval + 1
            niter = niter + 1

            ! Check the step
            if (alam < alamin) then
                ! Either the solution has converged due to negligible change in
                ! the root values, or the line search may have fully
                ! backtracked.  In the event the solution fully backtracked,
                ! we'll issue a warning to inform the user of the potential
                ! issue.
                if (norm2(x - xold) == zero) then
                    ! The line search fully backtracked
                    call errmgr%report_warning("ls_search_miso", &
                        "The line search appears to have fully " // &
                        "backtracked.  As such, check results carefully, " // &
                        "and/or consider attempting the solve without " // &
                        "the line search.", &
                        NL_CONVERGENCE_ERROR)
                end if
                x = xold
                xcnvrg = .true.
                exit
            else if (f <= fo + alpha * alam * slope) then
                ! The function has converged
                fcnvrg = .true.
                exit
            else
                ! Convergence has not yet occurred, continue backtracking
                tmplam = min_backtrack_search(niter, fo, f, f1, alam, &
                    alam1, slope)
            end if

            ! Set up parameters for the cubic model as we've already been
            ! through once with the quadratic model without success.
            alam1 = alam
            f1 = f
            alam = max(tmplam, lambdamin * alam)

            ! Ensure we haven't performed too many function evaluations
            if (neval >= maxeval) then
                ! ERROR: Too many function evaluations
                flag = 1
                exit
            end if
        end do
        if (present(fx)) fx = f

        ! Report out iteration statistics
        if (present(ib)) then
            ib%iter_count = niter
            ib%fcn_count = neval
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = .false.
        end if

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, 100) "The line search failed to " // &
                "converge.  Function evaluations performed: ", neval, "."
            call errmgr%report_error("ls_search_miso", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if

        ! Formatting
100     format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    pure function min_backtrack_search(mode, f0, f, f1, alam, alam1, &
            slope) result(lam)
        !! Minimizes either the quadratic or cubic representation for a
        !! backtracking-type line search.
        integer(int32), intent(in) :: mode
            !! Set to 1 to apply the quadratic model; else, any other value
            !! will apply the cubic model.
        real(real64), intent(in) :: f0
            !! The previous function value.
        real(real64), intent(in) :: f
            !! The current function value.
        real(real64), intent(in) :: f1
            !! The predicted function value.
        real(real64), intent(in) :: alam
            !! The step length scaling factor at f.
        real(real64), intent(in) :: alam1
            !! The step length scaling factor at f1.
        real(real64), intent(in) :: slope
            !! The slope of the direction vector.
        real(real64) :: lam
            !! The new step length scaling factor.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: p5 = 0.5d0
        real(real64), parameter :: two = 2.0d0
        real(real64), parameter :: three = 3.0d0

        ! Local Variables
        real(real64) :: rhs1, rhs2, a, b, disc

        ! Process
        if (mode == 1) then
            ! Use a quadratic function approximation
            lam = -slope / (two * (f - f0 - slope))
        else
            ! Use a cubic function
            rhs1 = f - f0 - alam * slope
            rhs2 = f1 - f0 - alam1 * slope
            a = (rhs1 / alam**2 - rhs2 / alam1**2) / (alam - alam1)
            b = (-alam1 * rhs1 / alam**2 + alam * rhs2 / alam1**2) / &
                (alam - alam1)
            if (a == zero) then
                lam = -slope / (two * b)
            else
                disc = b**2 - three * a * slope
                if (disc < zero) then
                    lam = p5 * alam
                else if (b <= zero) then
                    lam = (-b + sqrt(disc)) / (three * a)
                else
                    lam = -slope / (b + sqrt(disc))
                end if
            end if
            if (lam > p5 * alam) lam = p5 * alam
        end if
    end function

! ------------------------------------------------------------------------------
    subroutine limit_search_vector(x, lim)
        !! Provides a means of scaling the length of the search direction
        !! vector.
        real(real64), intent(inout), dimension(:) :: x
            !! On input, the search direction vector.  On output, the search 
            !! direction vector limited in length to that specified by lim.
            !!  If the vector is originally shorter than the limit length, no 
            !! change is made.
        real(real64), intent(in) :: lim
            !! The length limit value.

        ! Local Variables
        real(real64) :: mag

        ! Process
        mag = norm2(x)
        if (mag == 0.0d0) return
        if (mag > lim) x = (lim / mag) * x
    end subroutine

! ------------------------------------------------------------------------------
end module
