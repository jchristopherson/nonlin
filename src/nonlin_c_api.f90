! nonlin_c_api.f90

!> @brief Defines a C-friendly API to the NONLIN library.
module nonlin_c_api
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env
    use ferror
    use nonlin_constants
    use nonlin_core
    use nonlin_solve
    use nonlin_least_squares
    use nonlin_linesearch
    implicit none

    interface
        !> @brief Describes a function of one variable.
        !!
        !! @param[in] x The independent variable.
        !! @return The value of the function at @p x.
        function c_fcn1var(x) result(f)
            use, intrinsic :: iso_c_binding
            real(c_double), intent(in), value :: x
            real(c_double) :: f
        end function

        !> @brief Describes a function of N variables.
        !!
        !! @param[in] n The number of independent variables.
        !! @param[in] x An N-element array containing the independent variables.
        !! @return The value of the function at @p x.
        function c_fcnnvar(n, x) result(f)
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in) :: x(n)
            real(c_double) :: f
        end function

        !> @brief Describes a vector-valued function.
        !!
        !! @param[in] neqn The number of equations.
        !! @param[in] nvar The number of independent variables.
        !! @param[in] x An NVAR-element array containing the independent
        !!  variables.
        !! @param[out] f An NEQN-element array containing the values of the
        !!  function at @p x.
        subroutine c_vecfcn(neqn, nvar, x, f)
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(in), value :: neqn, nvar
            real(c_double), intent(in) :: x(nvar)
            real(c_double), intent(out) :: f(neqn)
        end subroutine

        !> @brief Describes a routin capable of computing the gradient vector
        !! of an equation of N variables.
        !!
        !! @param[in] n The number of independent variables.
        !! @param[in] x An N-element array containing the independent variables.
        !! @param[out] g An N-element array where the gradient vector will be
        !!  written as output.
        subroutine c_gradientfcn(n, x, g)
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out) :: g(n)
        end subroutine

        !> @brief Describes a routine capable of computing the Jacobian matrix
        !! of a system of equations.
        !!
        !! @param[in] neqn The number of equations.
        !! @param[in] nvar The number of independent variables.
        !! @param[in] x An NVAR-element array containing the independent 
        !!  variables.
        !! @param[out] jac An NEQN-by-NVAR matrix where the Jacobian will be
        !!  written.
        subroutine c_jacobianfcn(neqn, nvar, x, jac)
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(in), value :: neqn, nvar
            real(c_double), intent(in) :: x(nvar)
            real(c_double), intent(out) :: jac(neqn, nvar)
        end subroutine
    end interface

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

    !> @brief A type for providing controls on line search parameters.
    type, bind(C) :: line_search_controls
        !> @brief Defines whether line-searching should be employed.
        logical(c_bool) :: enable
        !> @brief Defines the maximum number of allowable function evaluations
        !! during a single line search.
        integer(c_int) :: max_function_evals
        !> @brief Defines the scaling of the product of the gradient and 
        !! direction vectors such that F(X + LAMBDA * P) <=
        !! F(X) + LAMBDA * ALPHA * P**T * G, where P is the search direction
        !! vector, G is the gradient vector, and LAMBDA is the scaling factor.
        !! The parameter must exist on the set (0, 1).  A value of 1e-4 is
        !! typically a good starting point.
        real(c_double) :: alpha
        !> @brief Defines a minimum factor X used to determine a minimum value 
        !! LAMBDA such that MIN(LAMBDA) = X * LAMBDA, where LAMBDA defines the 
        !! distance along the line search direction assuming a value of one 
        !! means the full length of the direction vector is traversed.  As such,
        !! the value must exist on the set (0, 1); however, for practical 
        !! considerations, the minimum value should be limited to 0.1 such that 
        !! the value must exist on the set [0.1, 1).
        real(c_double) :: factor
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
    !> @brief Defines standard line search settings.
    !!
    !! @param[out] x The line_search_controls item to populate.
    subroutine c_set_default_line_search_settings(x) &
            bind(C, name = "c_set_default_line_search_settings")
        ! Arguments
        type(line_search_controls), intent(out) :: x

        ! Process
        x%enable = logical(.false., c_bool)
        x%alpha = 1.0d-4
        x%factor = 1.0d-1
        x%max_function_evals = 100
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Utilizes Newton's method to sovle an equation of one variable.  
    !! The derivative calculation will be numerically estimated.
    !!
    !! @param[in] fcn The function to solve.
    !! @param[in] limits A set of limits on the solution bounding the solution
    !!  range thereby preventing the solver from wandering too far off course.
    !! @param[out] x The solution.
    !! @param[out] f The value of the function at the solution.
    !! @param[in] cntrls The iteration controls.
    !! @param[out] stats The iteration status.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !! - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !! - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within 
    !!      the allowed number of iterations.
    function c_solver_newton_1var(fcn, limits, x, f, cntrls, stats) &
            bind(C, name = "c_solver_newton_1var") result(flag)
        ! Arguments
        type(c_funptr), intent(in), value :: fcn
        type(value_pair), intent(in), value :: limits
        real(c_double), intent(out) :: x, f
        type(iteration_controls), intent(in) :: cntrls
        type(iteration_process), intent(out) :: stats
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        type(fcn1var_helper) :: obj
        procedure(c_fcn1var), pointer :: cfptr
        procedure(fcn1var), pointer :: fptr
        type(iteration_behavior) :: tracking
        type(newton_1var_solver) :: solver

        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)
        call c_f_procpointer(fcn, cfptr)
        fptr => fun
        call obj%set_fcn(fptr)

        ! Set the solver parameters
        call solver%set_max_fcn_evals(cntrls%max_function_evals)
        call solver%set_fcn_tolerance(cntrls%function_tolerance)
        call solver%set_var_tolerance(cntrls%solution_tolerance)
        call solver%set_diff_tolerance(cntrls%gradient_tolerance)
        call solver%set_print_status(logical(cntrls%print_status))

        ! Solve
        call solver%solve(obj, x, limits, f, ib = tracking, err = err)

        ! Retrieve iteration status
        stats%iteration_count = tracking%iter_count
        stats%function_eval_count = tracking%fcn_count
        stats%jacobian_eval_count = tracking%jacobian_count
        stats%gradient_eval_count = tracking%gradient_count
        stats%converge_on_function = logical(tracking%converge_on_fcn, c_bool)
        stats%converge_on_solution_change = logical(tracking%converge_on_chng, c_bool)
        stats%converge_on_gradient = logical(tracking%converge_on_zero_diff, c_bool)

        ! Check for errors
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    contains
        function fun(xi) result(fx)
            real(real64), intent(in) :: xi
            real(real64) :: fx
            fx = cfptr(xi)
        end function
    end function

! ------------------------------------------------------------------------------
    !> @brief Utilizes Brent's method to solve an equation of one variable.
    !!
    !! @param[in] fcn The function to solve.
    !! @param[in,out] x On input, the initial guess.  On output, the solution.
    !! @param[out] f The value of the function at the solution.
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
    function c_solver_brent_1var(fcn, x, f, limits, cntrls, stats) &
            bind(C, name = "c_solver_brent_1var") result(flag)
        ! Arguments
        type(c_funptr), intent(in), value :: fcn
        real(c_double), intent(inout) :: x
        real(c_double), intent(out) :: f
        type(value_pair), intent(in), value :: limits
        type(iteration_controls), intent(in) :: cntrls
        type(iteration_process), intent(out) :: stats
        integer(c_int) :: flag

        type(errors) :: err
        type(fcn1var_helper) :: obj
        procedure(c_fcn1var), pointer :: cfptr
        procedure(fcn1var), pointer :: fptr
        type(iteration_behavior) :: tracking
        type(brent_solver) :: solver

        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)
        call c_f_procpointer(fcn, cfptr)
        fptr => fun
        call obj%set_fcn(fptr)

        ! Set the solver parameters
        call solver%set_max_fcn_evals(cntrls%max_function_evals)
        call solver%set_fcn_tolerance(cntrls%function_tolerance)
        call solver%set_var_tolerance(cntrls%solution_tolerance)
        call solver%set_diff_tolerance(cntrls%gradient_tolerance)
        call solver%set_print_status(logical(cntrls%print_status))

        ! Solve
        call solver%solve(obj, x, limits, f, ib = tracking, err = err)

        ! Retrieve iteration status
        stats%iteration_count = tracking%iter_count
        stats%function_eval_count = tracking%fcn_count
        stats%jacobian_eval_count = tracking%jacobian_count
        stats%gradient_eval_count = tracking%gradient_count
        stats%converge_on_function = logical(tracking%converge_on_fcn, c_bool)
        stats%converge_on_solution_change = logical(tracking%converge_on_chng, c_bool)
        stats%converge_on_gradient = logical(tracking%converge_on_zero_diff, c_bool)

        ! Check for errors
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if

    contains
        function fun(xi) result(fx)
            real(real64), intent(in) :: xi
            real(real64) :: fx
            fx = cfptr(xi)
        end function
    end function

! ------------------------------------------------------------------------------
    !> @brief Utilizes Broyden's Quasi-Newton method to solve a system of N
    !! equations of N unknowns.  A backtracking type line search is also 
    !! employed.
    !!
    !! @param[in] fcn The function to solve.
    !! @param[in] jac A function for evaluating the Jacobian.  If null, the
    !!  Jacobian is estimated numerically.
    !! @param[in] n The number of equations.
    !! @param[in,out] x On input, an N-element array containing an initial
    !!  estimate to the solution.  On output, the updated solution estimate.
    !! @param[out] f An N-element array that, on output, will contain
    !!  the values of each equation as evaluated at the variable values
    !!  given in @p x.
    !! @param[in] cntrls The iteration controls.
    !! @param[in] ls The line search controls.
    !! @param[out] stats The iteration status.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !! - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !!  - NL_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      correctly.
    !! - NL_DIVERGENT_BEHAVIOR_ERROR: Occurs if the direction vector is
    !!      pointing in an apparent uphill direction.
    !! - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
    !!      the allowed number of iterations.
    !! - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !! - NL_SPURIOUS_CONVERGENCE_ERROR: Occurs as a warning if the slope of the
    !!      gradient vector becomes sufficiently close to zero.
    function c_solver_quasi_newton(fcn, jac, n, x, f, cntrls, ls, stats) &
            bind(C, name = "c_solver_quasi_newton") result(flag)
        ! Arguments
        type(c_funptr), intent(in), value :: fcn, jac
        integer(c_int), intent(in), value :: n
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: f(n)
        type(iteration_controls), intent(in) :: cntrls
        type(line_search_controls), intent(in) :: ls
        type(iteration_process), intent(out) :: stats
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        type(vecfcn_helper) :: obj
        procedure(c_vecfcn), pointer :: cfptr
        procedure(vecfcn), pointer :: fptr
        procedure(c_jacobianfcn), pointer :: cjptr
        procedure(jacobianfcn), pointer :: jptr
        type(iteration_behavior) :: tracking
        type(quasi_newton_solver) :: solver
        type(line_search) :: search
        
        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)
        call c_f_procpointer(fcn, cfptr)
        fptr => fun
        call obj%set_fcn(fptr, 2, 2)
        if (c_associated(jac)) then
            call c_f_procpointer(jac, cjptr)
            jptr => jacfun
            call obj%set_jacobian(jptr)
        end if

        ! Set the solver parameters
        call solver%set_max_fcn_evals(cntrls%max_function_evals)
        call solver%set_fcn_tolerance(cntrls%function_tolerance)
        call solver%set_var_tolerance(cntrls%solution_tolerance)
        call solver%set_gradient_tolerance(cntrls%gradient_tolerance)
        call solver%set_print_status(logical(cntrls%print_status))

        ! Set line search parameters
        call solver%set_use_line_search(logical(ls%enable))
        if (ls%enable) then
            call search%set_max_fcn_evals(ls%max_function_evals)
            call search%set_scaling_factor(ls%alpha)
            call search%set_distance_factor(ls%factor)
            call solver%set_line_search(search)
        end if

        ! Solve
        call solver%solve(obj, x, f, ib = tracking, err = err)

        ! Retrieve the ieration status
        stats%iteration_count = tracking%iter_count
        stats%function_eval_count = tracking%fcn_count
        stats%jacobian_eval_count = tracking%jacobian_count
        stats%gradient_eval_count = tracking%gradient_count
        stats%converge_on_function = logical(tracking%converge_on_fcn, c_bool)
        stats%converge_on_solution_change = logical(tracking%converge_on_chng, c_bool)
        stats%converge_on_gradient = logical(tracking%converge_on_zero_diff, c_bool)

        ! Check for errors
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    contains
        subroutine fun(xx, fx)
            real(real64), intent(in), dimension(:) :: xx
            real(real64), intent(out), dimension(:) :: fx
            call cfptr(n, n, xx, fx)
        end subroutine

        subroutine jacfun(xx, jx)
            real(real64), intent(in), dimension(:) :: xx
            real(real64), intent(out), dimension(:,:) :: jx
            call cjptr(n, n, xx, jx)
        end subroutine
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
