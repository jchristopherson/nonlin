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
    use ferror, only : errors
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
    !! @param[in] x The NVAR-element array containing the indendent variables.
    !! @param[out] f The NEQN-element array containing the function values.
    subroutine cvecfcn(neqn, nvar, x, f)
        use linalg_constants, only : dp, i32
        integer(i32), intent(in) :: neqn, nvar
        real(dp), intent(in) :: x(nvar)
        real(dp), intent(out) :: f(neqn)
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
    contains
        !> @brief Establishes a pointer to the routine containing the system of
        !!  equations to solve.
        procedure, public :: set_cfcn => cvfh_set_fcn
        ! !> @brief Establishes a pointer to the routine for computing the 
        ! !! Jacobian matrix of the system of equations.  If no routine is 
        ! !! defined, the Jacobian matrix will be computed numerically (this is 
        ! !! the default state).
        ! procedure, public :: set_jacobian => vfh_set_jac

        !> @brief Tests if the pointer to the subroutine containing the system
        !! of equations to solve has been assigned.
        procedure, public :: is_fcn_defined => cvfh_is_fcn_defined
        ! !> @brief Tests if the pointer to the subroutine containing the system
        ! !! of equations to solve has been assigned.
        ! procedure, public :: is_jacobian_defined => vfh_is_jac_defined

        !> @brief Executes the routine containing the system of equations to
        !! solve.  No action is taken if the pointer to the subroutine has not 
        !! been defined.
        procedure, public :: fcn => cvfh_fcn
        ! !> @brief Executes the routine containing the Jacobian matrix if 
        ! !! supplied.  If not supplied, the Jacobian is computed via finite 
        ! !! differences.
        ! procedure, public :: jacobian => vfh_jac_fcn
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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

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
    !! @param[in] err A pointer to the C error handler object.  If no error
    !!  handling is desired, simply pass NULL, and errors will be dealt with
    !!  by the default internal error handler.  Possible errors that may be
    !!  encountered are as follows.
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
        type(c_ptr), intent(in), value :: err

        ! Local Variables
        procedure(cfcn1var), pointer :: fptr
        type(errors), pointer :: eptr
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
        if (c_associated(err)) then
            call c_f_pointer(err, eptr)
            call solver%solve(obj, x, lim, f, ib, eptr)
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
    !! @param[in] err A pointer to the C error handler object.  If no error
    !!  handling is desired, simply pass NULL, and errors will be dealt with
    !!  by the default internal error handler.  Possible errors that may be
    !!  encountered are as follows.
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
        type(c_ptr), intent(in), value :: err

        ! Local Variables
        procedure(cvecfcn), pointer :: fptr
        procedure(jacobianfcn), pointer :: jptr
        type(errors), pointer :: eptr
        type(quasi_newton_solver) :: solver
        type(cvecfcn_helper) :: obj
        type(line_search) :: ls
        type(line_search_control), pointer :: lsc

        ! Initialization
        call c_f_procpointer(fcn, fptr)
        call obj%set_cfcn(fptr, n, n)
        if (c_associated(jac)) then
            call c_f_procpointer(jac, jptr)
            call obj%set_jacobian(jptr)
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
        if (c_associated(err)) then
            call c_f_pointer(err, eptr)
            call solver%solve(obj, x, fvec, ib, eptr)
        else
            call solver%solve(obj, x, fvec, ib)
        end if
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
