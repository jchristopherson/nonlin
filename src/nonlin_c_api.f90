! nonlin_c_api.f90

!> @brief Defines a C-friendly API to the NONLIN library.
module nonlin_c_api
    use iso_c_binding
    use ferror
    use nonlin_constants
    use nonlin_core
    use nonlin_solve
    use nonlin_least_squares
    implicit none

    !> @brief A type for providing a set of iteration control parameters.
    type, bind(C) :: iteration_controls
        !> @brief Defines the maximum number of allowable function evaluations.
        integer(c_int) :: max_function_evals
        !> @brief Defines convergence criteria based upon function value.
        real(c_double) :: function_tolerance
        !> @brief Defines convergence criteria based upon the solution value.
        real(c_double) :: solution_tolerance
        !> @brief Defines convergence criteria based upon the slope of the 
        !! gradient vector.
        real(c_double) :: gradient_tolerance
        !> @brief Gets a logical value determining if iteration status should
        !! be printed.
        logical(c_bool) :: print_status
    end type

    !> @brief A type providing information on the iteration process.
    type, bind(C) :: iteration_process
        !> @brief The number of iteration performed.
        integer(c_int) :: iteration_count
        !> @brief The number of function evaluations performed.  This typically
        !! does not include derivative evaulations.
        integer(c_int) :: function_eval_count
        !> @brief The number of Jacobian evaluations performed.
        integer(c_int) :: jacobian_eval_count
        !> @brief The number of gradient vector evaluations performed.
        integer(c_int) :: gradient_eval_count
        !> @brief True if the solution converged as a result of a zero-valued
        !! function; else, false.
        logical(c_bool) :: converge_on_function
        !> @brief True if the solution converged as a result of no appreciable
        !! change in the solution between iterations; else, false.
        logical(c_bool) :: converge_on_solution_change
        !> @brief True if the solution appears to have settled on a stationary
        !! point such that the gradient of the function is zero-valued; else,
        !! false.
        logical(c_bool) :: converge_on_gradient
    end type
contains
! ------------------------------------------------------------------------------
    !> @brief Defines standard solver settings.
    !!
    !! @param[out] x The iteration_controls item to populate.
    subroutine c_set_default_solver_settings(x) &
            bind(C, name = "c_set_default_solver_settings")
        ! Arguments
        type(iteration_controls), intent(out) :: x

        ! Process
        x%max_function_evals = 100
        x%function_tolerance = 1.0d-8
        x%solution_tolerance = 1.0d-12
        x%gradient_tolerance = 1.0d-12
        x%print_status = logical(.false., c_bool)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Utilizes Newton's method to sovle an equation of one variable.  
    !! The derivative calculation will be numerically estimated.
    !!
    !! @param[in] fcn The function to solve.
    !! @param[in,out] xi On input, the initial guess at the solution.  On
    !!  output, the solution.
    !! @param[out] fo The value of the function at the solution.
    !! @param[in] limits A set of limits on the solution bounding the solution
    !!  range thereby preventing the solver from wandering too far off course.
    !! @param[in] cntrls The iteration controls.
    !! @param[out] stats The iteration status.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !! - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !! - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within 
    !!      the allowed number of iterations.
    function c_solver_newton_1var(fcn, xi, fo, limits, cntrls, stats) &
            bind(C, name = "c_solver_newton_1var") result(flag)
        ! Arguments
        type(c_funptr), intent(in), value :: fcn
        real(c_double), intent(inout) :: xi
        real(c_double), intent(out) :: fo
        type(value_pair), intent(in) :: limits
        type(iteration_controls), intent(in) :: cntrls
        type(iteration_process), intent(out) :: stats
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        type(fcn1var_helper) :: obj
        procedure(fcn1var), pointer :: fptr
        type(iteration_behavior) :: tracking
        type(newton_1var_solver) :: solver

        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)
        call c_f_procpointer(fcn, fptr)
        call obj%set_fcn(fptr)

        ! Set the solver parameters
        call solver%set_max_fcn_evals(cntrls%max_function_evals)
        call solver%set_fcn_tolerance(cntrls%function_tolerance)
        call solver%set_var_tolerance(cntrls%solution_tolerance)
        call solver%set_diff_tolerance(cntrls%gradient_tolerance)
        call solver%set_print_status(logical(cntrls%print_status))

        ! Solve
        call solver%solve(obj, xi, limits, fo, ib = tracking, err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if

        ! Retrieve iteration status
        stats%iteration_count = tracking%iter_count
        stats%function_eval_count = tracking%fcn_count
        stats%jacobian_eval_count = tracking%jacobian_count
        stats%gradient_eval_count = tracking%gradient_count
        stats%converge_on_function = logical(tracking%converge_on_fcn, c_bool)
        stats%converge_on_solution_change = logical(tracking%converge_on_chng, c_bool)
        stats%converge_on_gradient = logical(tracking%converge_on_zero_diff, c_bool)
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
