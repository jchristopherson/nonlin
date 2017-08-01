! nonlin_c_binding.f90

!> @brief \b nonlin_c_binding
!!
!! @par Purpose
!! Provides C bindings to the nonlin library.
module nonlin_c_binding
    use, intrinsic :: iso_c_binding
    use linalg_constants, only : dp, i32
    use nonlin_types
    use nonlin_linesearch
    use nonlin_solve
    use nonlin_least_squares
    use nonlin_polynomials
    use ferror, only : errors
    use ferror_c_binding, only : errorhandler, get_errorhandler, &
        update_errorhandler
    implicit none


! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
interface
    !> @brief The C-friendly interface to fcn1var.
    !!
    !! @param[in] x The independent variable.
    !!
    !! @return The value of the function @p x.
    function cfcn1var(x) result(f)
        ! This is required as opposed to fcn1var in order to allow the C
        ! input to be passed by value, not as a pointer.
        use linalg_constants, only : dp
        real(dp), intent(in), value :: x
        real(dp) :: f
    end function

    !> @brief The C-friendly interface to vecfcn.
    !!
    !! @param[in] neqn The number of equations.
    !! @param[in] nvar The number of variables.
    !! @param[in] x The NVAR-element array containing the independent variables.
    !! @param[out] f The NEQN-element array containing the function values.
    subroutine cvecfcn(neqn, nvar, x, f)
        use linalg_constants, only : dp, i32
        integer(i32), intent(in), value :: neqn, nvar
        real(dp), intent(in) :: x(nvar)
        real(dp), intent(out) :: f(neqn)
    end subroutine

    !>  @brief The C-friendly interface to jacobianfcn.
    !!
    !! @param[in] neqn The number of equations.
    !! @param[in] nvar The number of variables.
    !! @param[in] x The NVAR-element array containing the independent variables.
    !! @param[out] jac An NEQN-byNVAR matrix where the Jacobian will be written.
    subroutine cjacobianfcn(neqn, nvar, x, jac)
        use linalg_constants, only : dp, i32
        integer(i32), intent(in), value :: neqn, nvar
        real(dp), intent(in) :: x(nvar)
        real(dp), intent(out) :: jac(neqn, nvar)
    end subroutine
end interface

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief Defines a set of solver control information.
    type, bind(c) :: solver_control
        !> The maximum number of function evaluations allowed.
        integer(i32) :: max_evals
        !> The convergence criteria on function values.
        real(dp) :: fcn_tolerance
        !> The convergence criteria on change in variable values.
        real(dp) :: var_tolerance
        !> The convergence criteria for the slope of the gradient vector.
        real(dp) :: grad_tolerance
        !> Controls whether iteration status is printed.
        logical(c_bool) :: print_status
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a set of line search controls.
    type, bind(c) :: line_search_control
        !> The maximum number of function evaluations allowed per search.
        integer(i32) :: max_evals
        !> Defines the scaling of the product of the gradient and direction
        !! vectors such that F(X + LAMBDA * P) <=
        !! F(X) + LAMBDA * ALPHA * P**T * G, where P is the search direction
        !! vector, G is the gradient vector, and LAMBDA is the scaling factor.
        !! The parameter must exist on the set (0, 1).  A value of 1e-4 is
        !! typically a good starting point.
        real(dp) :: alpha
        !> Defines a minimum factor X used to determine a minimum value LAMBDA
        !! such that MIN(LAMBDA) = X * LAMBDA, where LAMBDA defines the distance
        !! along the line search direction assuming a value of one means the
        !! full length of the direction vector is traversed.  As such, the value
        !! must exist on the set (0, 1); however, for practical considerations,
        !! the minimum value should be limited to 0.1 such that the value must
        !! exist on the set [0.1, 1).
        real(dp) :: factor
    end type

! ------------------------------------------------------------------------------
    !> @brief A C compatible type encapsulating a polynomial object.
    type, bind(C) :: c_polynomial
        !> @brief The size of the polynomial object, in bytes.
        integer(i32) :: n
        !> @brief A pointer to the polynomial object.
        type(c_ptr) :: ptr
    end type

! ------------------------------------------------------------------------------
    !> @brief A container allowing the use of cfcn1var in the solver codes.
    type, extends(fcn1var_helper) :: cfcn1var_helper
        private
        !> A pointer to the target cfcn1var routine.
        procedure(cfcn1var), pointer, nopass :: m_cfcn => null()
    contains
        !> @brief Executes the routine containing the function to evaluate.
        procedure, public :: fcn => cf1h_fcn
        !> @brief Tests if the pointer to the function containing the equation
        !! to solve has been assigned.
        procedure, public :: is_fcn_defined => cf1h_is_fcn_defined
        !> @brief Establishes a pointer to the routine containing the equations
        !! to solve.
        procedure, public :: set_cfcn => cf1h_set_fcn
    end type

! ------------------------------------------------------------------------------
    !> @brief A container allowing the use of cvecfcn in the solver codes.
    type, extends(vecfcn_helper) :: cvecfcn_helper
        private
        !> A pointer to the target cvecfcn routine.
        procedure(cvecfcn), pointer, nopass :: m_cfcn => null()
        !> A pointer to the Jacobian routine.
        procedure(cjacobianfcn), pointer, nopass :: m_cjac => null()
    contains
        !> @brief Establishes a pointer to the routine containing the system of
        !!  equations to solve.
        procedure, public :: set_cfcn => cvfh_set_fcn
        !> @brief Establishes a pointer to the routine for computing the
        !! Jacobian matrix of the system of equations.  If no routine is
        !! defined, the Jacobian matrix will be computed numerically (this is
        !! the default state).
        procedure, public :: set_cjacobian => cvfh_set_jac
        !> @brief Tests if the pointer to the subroutine containing the system
        !! of equations to solve has been assigned.
        procedure, public :: is_fcn_defined => cvfh_is_fcn_defined
        !> @brief Tests if the pointer to the subroutine containing the system
        !! of equations to solve has been assigned.
        procedure, public :: is_jacobian_defined => cvfh_is_jac_defined
        !> @brief Executes the routine containing the system of equations to
        !! solve.  No action is taken if the pointer to the subroutine has not
        !! been defined.
        procedure, public :: fcn => cvfh_fcn
        !> @brief Executes the routine containing the Jacobian matrix if
        !! supplied.  If not supplied, the Jacobian is computed via finite
        !! differences.
        procedure, public :: jacobian => cvfh_jac_fcn
    end type

contains
! ******************************************************************************
! CFCN1VAR_HELPER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Executes the routine containing the function to evaluate.
    !!
    !! @param[in] this The cfcn1var_helper object.
    !! @param[in] x The value of the independent variable at which the function
    !!  should be evaluated.
    !! @return The value of the function at @p x.
    function cf1h_fcn(this, x) result(f)
        class(cfcn1var_helper), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: f
        if (associated(this%m_cfcn)) then
            f = this%m_cfcn(x)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the function containing the equation to
    !! solve has been assigned.
    !!
    !! @param[in] this The cfcn1var_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function cf1h_is_fcn_defined(this) result(x)
        class(cfcn1var_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_cfcn)
    end function

! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine containing the equations to
    !! solve.
    !!
    !! @param[in,out] this The cfcn1var_helper object.
    !! @param[in] fcn The function pointer.
    subroutine cf1h_set_fcn(this, fcn)
        class(cfcn1var_helper), intent(inout) :: this
        procedure(cfcn1var), intent(in), pointer :: fcn
        this%m_cfcn => fcn
    end subroutine

! ******************************************************************************
! CVECFCN_HELPER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine containing the system of
    !!  equations to solve.
    !!
    !! @param[in,out] this The cvecfcn_helper object.
    !! @param[in] fcn The function pointer.
    !! @param[in] nfcn The number of functions.
    !! @param[in] nvar The number of variables.
    subroutine cvfh_set_fcn(this, fcn, nfcn, nvar)
        class(cvecfcn_helper), intent(inout) :: this
        procedure(cvecfcn), intent(in), pointer :: fcn
        integer(i32), intent(in) :: nfcn, nvar
        procedure(vecfcn), pointer :: nptr
        nptr => null()
        call this%set_fcn(nptr, nfcn, nvar)
        this%m_cfcn => fcn
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the subroutine containing the system of
    !! equations to solve has been assigned.
    !!
    !! @param[in] this The cvecfcn_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function cvfh_is_fcn_defined(this) result(x)
        class(cvecfcn_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_cfcn)
    end function

! ------------------------------------------------------------------------------
    !> @brief Executes the routine containing the system of equations to solve.
    !! No action is taken if the pointer to the subroutine has not been defined.
    !!
    !! @param[in] this The cvecfcn_helper object.
    !! @param[in] x An N-element array containing the independent variables.
    !! @param[out] f An M-element array that, on output, contains the values
    !!  of the M functions.
    subroutine cvfh_fcn(this, x, f)
        class(cvecfcn_helper), intent(in) :: this
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: f
        integer(i32) :: neqn, nvar
        neqn = this%get_equation_count()
        nvar = this%get_variable_count()
        if (this%is_fcn_defined()) then
            call this%m_cfcn(neqn, nvar, x, f)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine for computing the Jacobian
    !! matrix of the system of equations.  If no routine is defined, the
    !! Jacobian matrix will be computed numerically (this is the default state).
    !!
    !! @param[in,out] this The cvecfcn_helper object.
    !! @param[in] jac The function pointer.
    subroutine cvfh_set_jac(this, jac)
        class(cvecfcn_helper), intent(inout) :: this
        procedure(cjacobianfcn), intent(in), pointer :: jac
        this%m_cjac => jac
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the subroutine containing the system of
    !! equations to solve has been assigned.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function cvfh_is_jac_defined(this) result(x)
        class(cvecfcn_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_cjac)
    end function

! ------------------------------------------------------------------------------
    !> @brief Executes the routine containing the Jacobian matrix if supplied.
    !! If not supplied, the Jacobian is computed via finite differences.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @param[in] x An N-element array containing the independent variabls
    !!  defining the point about which the derivatives will be calculated.
    !! @param[out] jac An M-by-N matrix where, on output, the Jacobian will
    !!  be written.
    !! @param[in] fv An optional M-element array containing the function values
    !!  at @p x.  If not supplied, the function values are computed at @p x.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.  Notice, a workspace array is only utilized if the user does
    !!  not provide a routine for computing the Jacobian.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional integer output that can be used to determine
    !!  error status.  If not used, and an error is encountered, the routine
    !!  simply returns silently.  If used, the following error codes identify
    !!  error status:
    !!  - 0: No error has occurred.
    !!  - n: A positive integer denoting the index of an invalid input.
    !!  - -1: Indicates internal memory allocation failed.
    subroutine cvfh_jac_fcn(this, x, jac, fv, work, olwork, err)
        ! Arguments
        class(cvecfcn_helper), intent(in) :: this
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), dimension(:,:) :: jac
        real(dp), intent(in), dimension(:), optional, target :: fv
        real(dp), intent(out), dimension(:), optional, target :: work
        integer(i32), intent(out), optional :: olwork, err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: j, m, n, lwork, flag
        real(dp) :: eps, epsmch, h, temp
        real(dp), pointer, dimension(:) :: fptr, f1ptr
        real(dp), allocatable, target, dimension(:) :: wrk

        ! Initialization
        if (present(err)) err = 0
        m = this%get_equation_count()
        n = this%get_variable_count()

        ! Input Checking
        flag = 0
        if (size(x) /= n) then
            flag = 2
        else if (size(jac, 1) /= m .or. size(jac, 2) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: Incorrectly sized input arrays
            if (present(err)) err = flag
            return
        end if

        ! Process
        if (.not.this%is_fcn_defined()) return
        if (associated(this%m_cjac)) then
            ! Workspace Query
            if (present(olwork)) then
                olwork = 0
                return
            end if

            ! Call the user-defined Jacobian routine
            call this%m_cjac(m, n, x, jac)
        else
            ! Compute the Jacobian via finite differences
            if (present(fv)) then
                lwork = m
            else
                lwork = 2 * m
            end if

            if (present(olwork)) then
                ! The user is just making a workspace query.  Simply return the
                ! workspace length, and exit the routine.
                olwork = lwork
                return
            end if

            ! Local Memory Allocation
            if (present(work)) then
                if (size(work) < lwork) then
                    ! ERROR: Workspace is too small
                    if (present(err)) err = 5
                    return
                end if
                f1ptr => work(1:m)
                if (present(fv)) then
                    if (size(fv) < m) then
                        ! ERROR: Function vector too small
                        if (present(err)) err = 4
                        return
                    end if
                    fptr => fv(1:m)
                else
                    fptr => work(m+1:2*m)
                    call this%fcn(x, fptr)
                end if
            else
                allocate(wrk(lwork), stat = flag)
                if (flag /= 0) then
                    ! ERROR: Memory issues
                    if (present(err)) err = -1
                    return
                end if
                f1ptr => wrk(1:m)
                if (present(fv)) then
                    fptr => fv(1:m)
                else
                    fptr => wrk(m+1:2*m)
                    call this%fcn(x, fptr)
                end if
            end if

            ! Establish step size factors
            epsmch = epsilon(epsmch)
            eps = sqrt(epsmch)

            ! Compute the derivatives via finite differences
            do j = 1, n
                temp = x(j)
                h = eps * abs(temp)
                if (h == zero) h = eps
                x(j) = temp + h
                call this%fcn(x, f1ptr)
                x(j) = temp
                jac(:,j) = (f1ptr - fptr) / h
            end do
        end if
    end subroutine

! ******************************************************************************
! SOLVER ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Solves an equation of one variable using Brent's method.
    !!
    !! @param[in] fcn A pointer to the routine containing the function to solve.
    !! @param[in] lim A value_pair object defining the search limits.
    !! @param[out] x On output, the solution.
    !! @param[out] f On output, the residual as computed at @p x.
    !! @param[in] tol A solver_control object defining the solver control
    !!  parameters.
    !! @param[out] ib On output, an iteration_behavior object containing the
    !!  iteration performance statistics.
    !! @param[in] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is different
    !!      than the number of variables.
    !!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
    !!      the allowed number of iterations.
    subroutine brent_solver_c(fcn, lim, x, f, tol, ib, err) &
            bind(C, name = "solve_brent")
        ! Arguments
        type(c_funptr), intent(in), value :: fcn
        type(value_pair), intent(in), value :: lim
        real(dp), intent(out) :: x, f
        type(solver_control), intent(in) :: tol
        type(iteration_behavior), intent(out) :: ib
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        procedure(cfcn1var), pointer :: fptr
        class(errors), allocatable :: eptr
        type(brent_solver) :: solver
        type(cfcn1var_helper) :: obj

        ! Initialization
        call c_f_procpointer(fcn, fptr)
        call solver%set_max_fcn_evals(tol%max_evals)
        call solver%set_fcn_tolerance(tol%fcn_tolerance)
        call solver%set_var_tolerance(tol%var_tolerance)
        call solver%set_print_status(logical(tol%print_status))
        call obj%set_cfcn(fptr)

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call solver%solve(obj, x, lim, f, ib, eptr)
            call update_errorhandler(eptr, err)
        else
            call solver%solve(obj, x, lim, f, ib)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Applies the quasi-Newton's method developed by Broyden in
    !! conjunction with a backtracking type line search to solve N equations
    !! of N unknowns.
    !!
    !! @param[in] fcn A pointer to the routine containing the system of
    !!  equations to solve.
    !! @param[in] jac A pointer to a routine used to compute the Jacobian of
    !!  the system of equations.  To let the program compute the Jacobian
    !!  numerically, simply pass NULL.
    !! @param[in] n The number of equations, and the number of unknowns.
    !! @param[in,out] x On input, an N-element array containing an initial
    !!  estimate to the solution.  On output, the updated solution estimate.
    !!  N is the number of variables.
    !! @param[out] fvec An N-element array that, on output, will contain
    !!  the values of each equation as evaluated at the variable values
    !!  given in @p x.
    !! @param[in] tol A solver_control object defining the solver control
    !!  parameters.
    !! @param[in] lsearch A pointer to a line_search_control object defining
    !!  the line search control parameters.  If no line search is desired,
    !!  simply pass NULL.
    !! @param[out] ib On output, an iteration_behavior object containing the
    !!  iteration performance statistics.
    !! @param[in] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
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
    subroutine quasi_newton_c(fcn, jac, n, x, fvec, tol, lsearch, ib, err) &
            bind(C, name = "solve_quasi_newton")
        ! Arguments
        type(c_funptr), intent(in), value :: fcn, jac
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: x(n)
        real(dp), intent(out) :: fvec(n)
        type(solver_control), intent(in) :: tol
        type(c_ptr), intent(in), value :: lsearch
        type(iteration_behavior), intent(out) :: ib
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        procedure(cvecfcn), pointer :: fptr
        procedure(cjacobianfcn), pointer :: jptr
        class(errors), allocatable :: eptr
        type(quasi_newton_solver) :: solver
        type(cvecfcn_helper) :: obj
        type(line_search) :: ls
        type(line_search_control), pointer :: lsc

        ! Initialization
        call c_f_procpointer(fcn, fptr)
        call obj%set_cfcn(fptr, n, n)
        if (c_associated(jac)) then
            call c_f_procpointer(jac, jptr)
            call obj%set_cjacobian(jptr)
        end if
        call solver%set_max_fcn_evals(tol%max_evals)
        call solver%set_fcn_tolerance(tol%fcn_tolerance)
        call solver%set_var_tolerance(tol%var_tolerance)
        call solver%set_gradient_tolerance(tol%grad_tolerance)
        call solver%set_print_status(logical(tol%print_status))
        if (c_associated(lsearch)) then
            ! Use a line search
            call c_f_pointer(lsearch, lsc)
            call ls%set_max_fcn_evals(lsc%max_evals)
            call ls%set_scaling_factor(lsc%alpha)
            call ls%set_distance_factor(lsc%factor)
            call solver%set_use_line_search(.true.)
            call solver%set_line_search(ls)
        else
            ! Do not use a line search
            call solver%set_use_line_search(.false.)
        end if

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call solver%solve(obj, x, fvec, ib, eptr)
            call update_errorhandler(eptr, err)
        else
            call solver%solve(obj, x, fvec, ib)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Applies Newton's method in conjunction with a backtracking
    !! type line search to solve N equations of N unknowns.
    !!
    !! @param[in] fcn A pointer to the routine containing the system of
    !!  equations to solve.
    !! @param[in] jac A pointer to a routine used to compute the Jacobian of
    !!  the system of equations.  To let the program compute the Jacobian
    !!  numerically, simply pass NULL.
    !! @param[in] n The number of equations, and the number of unknowns.
    !! @param[in,out] x On input, an N-element array containing an initial
    !!  estimate to the solution.  On output, the updated solution estimate.
    !!  N is the number of variables.
    !! @param[out] fvec An N-element array that, on output, will contain
    !!  the values of each equation as evaluated at the variable values
    !!  given in @p x.
    !! @param[in] tol A solver_control object defining the solver control
    !!  parameters.
    !! @param[in] lsearch A pointer to a line_search_control object defining
    !!  the line search control parameters.  If no line search is desired,
    !!  simply pass NULL.
    !! @param[out] ib On output, an iteration_behavior object containing the
    !!  iteration performance statistics.
    !! @param[in] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
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
    subroutine newton_c(fcn, jac, n, x, fvec, tol, lsearch, ib, err) &
            bind(C, name = "solve_newton")
        ! Arguments
        type(c_funptr), intent(in), value :: fcn, jac
        integer(i32), intent(in), value :: n
        real(dp), intent(inout) :: x(n)
        real(dp), intent(out) :: fvec(n)
        type(solver_control), intent(in) :: tol
        type(c_ptr), intent(in), value :: lsearch
        type(iteration_behavior), intent(out) :: ib
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        procedure(cvecfcn), pointer :: fptr
        procedure(cjacobianfcn), pointer :: jptr
        class(errors), allocatable :: eptr
        type(newton_solver) :: solver
        type(cvecfcn_helper) :: obj
        type(line_search) :: ls
        type(line_search_control), pointer :: lsc

        ! Initialization
        call c_f_procpointer(fcn, fptr)
        call obj%set_cfcn(fptr, n, n)
        if (c_associated(jac)) then
            call c_f_procpointer(jac, jptr)
            call obj%set_cjacobian(jptr)
        end if
        call solver%set_max_fcn_evals(tol%max_evals)
        call solver%set_fcn_tolerance(tol%fcn_tolerance)
        call solver%set_var_tolerance(tol%var_tolerance)
        call solver%set_gradient_tolerance(tol%grad_tolerance)
        call solver%set_print_status(logical(tol%print_status))
        if (c_associated(lsearch)) then
            ! Use a line search
            call c_f_pointer(lsearch, lsc)
            call ls%set_max_fcn_evals(lsc%max_evals)
            call ls%set_scaling_factor(lsc%alpha)
            call ls%set_distance_factor(lsc%factor)
            call solver%set_use_line_search(.true.)
            call solver%set_line_search(ls)
        else
            ! Do not use a line search
            call solver%set_use_line_search(.false.)
        end if

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call solver%solve(obj, x, fvec, ib, eptr)
            call update_errorhandler(eptr, err)
        else
            call solver%solve(obj, x, fvec, ib)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Applies the Levenberg-Marquardt method to solve the nonlinear
    !! least-squares problem.
    !!
    !! @param[in] fcn A pointer to the routine containing the system of
    !!  equations to solve.
    !! @param[in] jac A pointer to a routine used to compute the Jacobian of
    !!  the system of equations.  To let the program compute the Jacobian
    !!  numerically, simply pass NULL.
    !! @param[in] neqn The number of equations.
    !! @param[in] nvar The number of unknowns.  This must be less than or equal
    !!  to @p neqn.
    !! @param[in,out] x On input, an N-element array containing an initial
    !!  estimate to the solution.  On output, the updated solution estimate.
    !!  N is the number of variables.
    !! @param[out] fvec An N-element array that, on output, will contain
    !!  the values of each equation as evaluated at the variable values
    !!  given in @p x.
    !! @param[in] tol A solver_control object defining the solver control
    !!  parameters.
    !! @param[out] ib On output, an iteration_behavior object containing the
    !!  iteration performance statistics.
    !! @param[in] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is less than
    !!      than the number of variables.
    !!  - NL_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      correctly.
    !!  - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
    !!      the allowed number of iterations.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - NL_TOLERANCE_TOO_SMALL_ERROR: Occurs if the requested tolerance is
    !!      to small to be practical for the problem at hand.
    subroutine levmarq_c(fcn, jac, neqn, nvar, x, fvec, tol, ib, err) &
            bind(C, name = "solve_nl_least_squares")
        ! Arguments
        type(c_funptr), intent(in), value :: fcn, jac
        integer(i32), intent(in), value :: neqn, nvar
        real(dp), intent(inout) :: x(nvar)
        real(dp), intent(out) :: fvec(neqn)
        type(solver_control), intent(in) :: tol
        type(iteration_behavior), intent(out) :: ib
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        procedure(cvecfcn), pointer :: fptr
        procedure(cjacobianfcn), pointer :: jptr
        class(errors), allocatable :: eptr
        type(least_squares_solver) :: solver
        type(cvecfcn_helper) :: obj

        ! Initialization
        call c_f_procpointer(fcn, fptr)
        call obj%set_cfcn(fptr, neqn, nvar)
        if (c_associated(jac)) then
            call c_f_procpointer(jac, jptr)
            call obj%set_cjacobian(jptr)
        end if
        call solver%set_max_fcn_evals(tol%max_evals)
        call solver%set_fcn_tolerance(tol%fcn_tolerance)
        call solver%set_var_tolerance(tol%var_tolerance)
        call solver%set_gradient_tolerance(tol%grad_tolerance)
        call solver%set_print_status(logical(tol%print_status))

        ! Process
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call solver%solve(obj, x, fvec, ib, eptr)
            call update_errorhandler(eptr, err)
        else
            call solver%solve(obj, x, fvec, ib)
        end if
    end subroutine

! ******************************************************************************
! OPTIMIZATION ROUTINES
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------


! ******************************************************************************
! POLYNOMIAL ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Initializes a new polynomial object.
    !!
    !! @param[out] obj The c_polynomial object to initialize.
    !! @param[in] order The order of the polynomial.  This value must be > 0.
    subroutine alloc_polynomial(obj, order) bind(C, name = "alloc_polynomial")
        ! Arguments
        type(c_polynomial), intent(out) :: obj
        integer(i32), intent(in), value :: order

        ! Local Variables
        integer(c_short), allocatable, target, dimension(:) :: temp
        type(polynomial) :: poly

        ! Process
        call poly%initialize(order)
        temp = transfer(poly, temp)
        obj%n = size(temp)
        obj%ptr = c_loc(temp(1))
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Retrieves the polynomial object from the C compatible 
    !! c_polynomial data structure.
    !!
    !! @param[in] obj The C compatible c_polynomial data structure.
    !! @param[out] poly The resulting polynomials object.
    subroutine get_polynomial(obj, poly)
        ! Arguments
        type(c_polynomial), intent(in), target :: obj
        class(polynomial), intent(out), allocatable :: poly

        ! Local Variables
        integer(c_short), pointer, dimension(:) :: temp
        type(polynomial) :: item
        type(c_ptr) :: testptr

        ! Process
        testptr = c_loc(obj) ! Ensures that obj wasn't passed as NULL from C
        if (.not.c_associated(testptr)) return
        if (.not.c_associated(obj%ptr)) return
        call c_f_pointer(obj%ptr, temp, shape = [obj%n])
        item = transfer(temp, item)
        allocate(poly, source = item)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Updates the c_polynomial object.
    !!
    !! @param[in] poly The polynomial object.
    !! @param[out] cobj The c_polynomial object to update.
    subroutine update_polynomial(poly, cobj)
        ! Arguments
        class(polynomial), intent(in) :: poly
        type(c_polynomial), intent(out) :: cobj

        ! Local Variables
        integer(c_short), allocatable, target, dimension(:) :: temp

        ! Process
        temp = transfer(poly, temp)
        cobj%n = size(temp)
        cobj%ptr = c_loc(temp(1))
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the order of the polynomial.
    !!
    !! @param[in] poly The c_polynomial object.
    !! @return The order of the polynomial object.
    function get_polynomial_order_c(poly) result(n) &
            bind(C, name = "get_polynomial_order")
        ! Arguments
        type(c_polynomial), intent(in) :: poly
        integer(i32) :: n

        ! Local Variables
        class(polynomial), allocatable :: pptr

        ! Process
        n = 0
        call get_polynomial(poly, pptr)
        if (.not.allocated(pptr)) return
        n = pptr%order()
    end function

! ------------------------------------------------------------------------------
    !> @brief Fits a polynomial of the specified order to a data set.
    !!
    !! @param[out] poly The c_polynomial object to initialize.
    !! @param[in] n The size of the arrays.
    !! @param[in] x An N-element array containing the independent variable data
    !!  points.  Notice, must be N > @p order.
    !! @param[in,out] y On input, an N-element array containing the dependent
    !!  variable data points.  On output, the contents are overwritten.
    !! @param[in] order The order of the polynomial (must be >= 1).
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if a zero or negative polynomial order
    !!      was specified, or if order is too large for the data set.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
    !!  - NL_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are different sizes.
    subroutine fit_polynomial(poly, n, x, y, order, err) &
            bind(C, name = "fit_polynomial")
        ! Arguments
        type(c_polynomial), intent(out) :: poly
        integer(i32), intent(in), value :: n, order
        real(dp), intent(in) :: x(n)
        real(dp), intent(inout) :: y(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr
        class(polynomial), allocatable :: pptr

        ! Process
        call get_polynomial(poly, pptr)
        if (.not.allocated(pptr)) return
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call pptr%fit(x, y, order, eptr)
            call update_errorhandler(eptr, err)
        else
            call pptr%fit(x, y, order)
        end if
        call update_polynomial(pptr, poly)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Fits a polynomial of the specified order that passes through zero
    !! to a data set.
    !!
    !! @param[out] poly The c_polynomial object to initialize.
    !! @param[in] n The size of the arrays.
    !! @param[in] x An N-element array containing the independent variable data
    !!  points.  Notice, must be N > @p order.
    !! @param[in,out] y On input, an N-element array containing the dependent
    !!  variable data points.  On output, the contents are overwritten.
    !! @param[in] order The order of the polynomial (must be >= 1).
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if a zero or negative polynomial order
    !!      was specified, or if order is too large for the data set.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
    !!  - NL_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are different sizes.
    subroutine fit_polynomial_thru_zero(poly, n, x, y, order, err) &
            bind(C, name = "fit_polynomial_thru_zero")
        ! Arguments
        type(c_polynomial), intent(out) :: poly
        integer(i32), intent(in), value :: n, order
        real(dp), intent(in) :: x(n)
        real(dp), intent(inout) :: y(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(errors), allocatable :: eptr
        class(polynomial), allocatable :: pptr

        ! Process
        call get_polynomial(poly, pptr)
        if (.not.allocated(pptr)) return
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call pptr%fit_thru_zero(x, y, order, eptr)
            call update_errorhandler(eptr, err)
        else
            call pptr%fit_thru_zero(x, y, order)
        end if
        call update_polynomial(pptr, poly)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Evaluates a polynomial at the specified points.
    !!
    !! @param[in] poly The c_polynomial object.
    !! @param[in] n The number of points to evaluate.
    !! @param[in] x An N-element array containing the points at which to
    !!  evaluate the polynomial.
    !! @param[out] y An N-element array where the resulting polynomial outputs
    !!  will be written.
    subroutine evaluate_polynomial(poly, n, x, y) &
            bind(C, name = "evaluate_polynomial")
        ! Arguments
        type(c_polynomial), intent(in) :: poly
        integer(i32), intent(in), value :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: y(n)

        ! Local Variables
        class(polynomial), allocatable :: pptr

        ! Process
        call get_polynomial(poly, pptr)
        if (.not.allocated(pptr)) return
        y = pptr%evaluate(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Evaluates a polynomial at the specified points.
    !!
    !! @param[in] poly The c_polynomial object.
    !! @param[in] n The number of points to evaluate.
    !! @param[in] x An N-element array containing the points at which to
    !!  evaluate the polynomial.
    !! @param[out] y An N-element array where the resulting polynomial outputs
    !!  will be written.
    subroutine evaluate_polynomial_cmplx(poly, n, x, y) &
            bind(C, name = "evaluate_polynomial_cmplx")
        ! Arguments
        type(c_polynomial), intent(in) :: poly
        integer(i32), intent(in), value :: n
        complex(dp), intent(in) :: x(n)
        complex(dp), intent(out) :: y(n)

        ! Local Variables
        class(polynomial), allocatable :: pptr

        ! Process
        call get_polynomial(poly, pptr)
        if (.not.allocated(pptr)) return
        y = pptr%evaluate(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes all the roots of a polynomial by computing the
    !! eigenvalues of the polynomial companion matrix.
    !!
    !! @param[in] poly The c_polynomial object.
    !! @param[in] n The size of @p rts.  This value should be the same as the
    !!  order of the polynomial.
    !! @param[out] rts An N-element array where the roots of the polynomial
    !!  will be written.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !!  - LA_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - LA_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    subroutine polynomial_roots_c(poly, n, rts, err) &
            bind(C, name = "polynomial_roots")
        ! Arguments
        type(c_polynomial), intent(in) :: poly
        integer(i32), intent(in), value :: n
        complex(dp), intent(out) :: rts(n)
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(polynomial), allocatable :: pptr
        class(errors), allocatable :: eptr
        complex(dp), allocatable, dimension(:) :: roots
        integer(i32) :: m

        ! Process
        call get_polynomial(poly, pptr)
        if (.not.allocated(pptr)) return
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            roots = pptr%roots(eptr)
            call update_errorhandler(eptr, err)
        else
            roots = pptr%roots()
        end if
        m = min(n, size(roots))
        rts(1:m) = roots(1:m)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the requested polynomial coefficient by index.  The
    !! coefficient index is established as follows: c(1) + c(2) * x +
    !! c(3) * x**2 + ... c(n) * x**n-1.
    !!
    !! @param[in] poly The c_polynomial object.
    !! @param[in] ind The polynomial coefficient index (0 < ind <= order + 1).
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !! - NL_INVALID_INPUT_ERROR: Occurs if the requested index is less than or
    !!      equal to zero, or if the requested index exceeds the number of
    !!      polynomial coefficients.
    function get_polynomial_coefficient(poly, ind, err) result(x) &
            bind(C, name = "get_polynomial_coefficient")
        ! Arguments
        type(c_polynomial), intent(in) :: poly
        integer(i32), intent(in), value :: ind
        type(errorhandler), intent(inout) :: err
        real(dp) :: x

        ! Local Variables
        class(polynomial), allocatable :: pptr
        class(errors), allocatable :: eptr

        ! Process
        x = 0.0d0
        call get_polynomial(poly, pptr)
        if (.not.allocated(pptr)) return
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            x = pptr%get(ind, eptr)
            call update_errorhandler(eptr, err)
        else
            x = pptr%get(ind)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Sets the requested polynomial coefficient by index.  The
    !! coefficient index is established as follows: c(1) + c(2) * x +
    !! c(3) * x**2 + ... c(n) * x**n-1.
    !!
    !! @param[in,out] poly The c_polynomial object.
    !! @param[in] ind The polynomial coefficient index (0 < ind <= order + 1).
    !! @param[in] x The polynomial coefficient.
    !! @param[in,out] err The errorhandler object.  If no error handling is
    !!  desired, simply pass NULL, and errors will be dealt with by the default
    !!  internal error handler.  Possible errors that may be encountered are as
    !!  follows.
    !! - NL_INVALID_INPUT_ERROR: Occurs if the requested index is less than or
    !!      equal to zero, or if the requested index exceeds the number of
    !!      polynomial coefficients.
    subroutine set_polynomial_coefficient(poly, ind, x, err) &
            bind(C, name = "set_polynomial_coefficient")
        ! Arguments
        type(c_polynomial), intent(inout) :: poly
        integer(i32), intent(in), value :: ind
        real(dp), intent(in), value :: x
        type(errorhandler), intent(inout) :: err

        ! Local Variables
        class(polynomial), allocatable :: pptr
        class(errors), allocatable :: eptr

        ! Process
        call get_polynomial(poly, pptr)
        if (.not.allocated(pptr)) return
        call get_errorhandler(err, eptr)
        if (allocated(eptr)) then
            call pptr%set(ind, x, eptr)
            call update_errorhandler(eptr, err)
        else
            call pptr%set(ind, x)
        end if
        call update_polynomial(pptr, poly)
    end subroutine

! ------------------------------------------------------------------------------
end module
