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
    use nonlin_optimize
    use nonlin_polynomials
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

    !> @brief A type for representing a polynomial in C that is compatible with
    !! the polynomial type in this library.
    type, bind(C) :: c_polynomial
        !> @brief The size of the polynomial object, in bytes.
        integer(c_int) :: size_in_bytes
        !> @brief A pointer to the underlying Fortran polynomial object.
        type(c_ptr) :: ptr
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
        x%max_function_evals = 500
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
        call obj%set_fcn(fptr, n, n)
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
    !> @brief Utilizes Newton's method to solve a system of N equations of N 
    !! unknowns in conjuction with a backtracking type line search.
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
    function c_solver_newton(fcn, jac, n, x, f, cntrls, ls, stats) &
            bind(C, name = "c_solver_newton") result(flag)
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
        type(newton_solver) :: solver
        type(line_search) :: search
        
        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)
        call c_f_procpointer(fcn, cfptr)
        fptr => fun
        call obj%set_fcn(fptr, n, n)
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
    !> @brief Utilizes the Levenberg-Marquardt method to solve a least-squares
    !! problem of M equations of N unknowns.  There must be at least as many
    !! equations as unknowns for this solver.
    !!
    !! @param[in] fcn The function to solve.
    !! @param[in] jac A function for evaluating the Jacobian.  If null, the
    !!  Jacobian is estimated numerically.
    !! @param[in] neqn The number of equations.
    !! @param[in] nvar The number of variables.
    !! @param[in,out] x On input, an NVAR element array containing the initial 
    !!  estimate to the solution.  On output, the solution.
    !! @param[out] f An NEQN-element array that, on output, will contain the
    !!  values of each equation as evaluated at the output solution given in
    !!  @p x.
    !! @param[in] cntrls The iteration controls.
    !! @param[out] stats The iteration status.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !! - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !! - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is less than
    !!      than the number of variables.
    !! - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
    !!      the allowed number of iterations.
    !! - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !! - NL_TOLERANCE_TOO_SMALL_ERROR: Occurs if the requested tolerance is
    !!      to small to be practical for the problem at hand.
    function c_solver_least_squares(fcn, jac, neqn, nvar, x, f, cntrls, stats) &
            bind(C, name = "c_solver_least_squares") result(flag)
        ! Arguments
        type(c_funptr), intent(in), value :: fcn, jac
        integer(c_int), intent(in), value :: neqn, nvar
        real(c_double), intent(inout) :: x(nvar)
        real(c_double), intent(out) :: f(neqn)
        type(iteration_controls), intent(in) :: cntrls
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
        type(least_squares_solver) :: solver
        
        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)
        call c_f_procpointer(fcn, cfptr)
        fptr => fun
        call obj%set_fcn(fptr, neqn, nvar)
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
            call cfptr(neqn, nvar, xx, fx)
        end subroutine

        subroutine jacfun(xx, jx)
            real(real64), intent(in), dimension(:) :: xx
            real(real64), intent(out), dimension(:,:) :: jx
            call cjptr(neqn, nvar, xx, jx)
        end subroutine
    end function

! ------------------------------------------------------------------------------
    !> @brief Utilizes Nelder and Mead's simplex algorithm to minimize a
    !! function of N variables.
    !!
    !! @param[in] fcn The function to minimize.
    !! @param[in] n The number of variables.
    !! @param[in,out] x On input, an N-element array containing an initial
    !!  estimate to the solution.  On output, the updated solution estimate.
    !! @param[out] f On output, the value of the function at @p x.
    !! @param[in] cntrls The iteration controls.
    !! @param[out] stats The iteration status.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !! - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !! - NL_INVALID_INPUT_ERROR: Occurs if @p x is not appropriately sized for
    !!      the problem as defined in @p fcn.
    !! - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !! - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
    !!      the allowed number of iterations.
    function c_solver_nelder_mead(fcn, n, x, f, cntrls, stats) &
            bind(C, name = "c_solver_nelder_mead") result(flag)
        ! Arguments
        type(c_funptr), intent(in), value :: fcn
        integer(c_int), intent(in), value :: n
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: f
        type(iteration_controls), intent(in) :: cntrls
        type(iteration_process), intent(out) :: stats
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        type(fcnnvar_helper) :: obj
        procedure(c_fcnnvar), pointer :: cfptr
        procedure(fcnnvar), pointer :: fptr
        type(iteration_behavior) :: tracking
        type(nelder_mead) :: solver

        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)
        call c_f_procpointer(fcn, cfptr)
        fptr => fun
        call obj%set_fcn(fptr, n)

        ! Set solver parameters
        call solver%set_max_fcn_evals(cntrls%max_function_evals)
        call solver%set_tolerance(cntrls%function_tolerance)
        call solver%set_print_status(logical(cntrls%print_status))

        ! Solve
        call solver%solve(obj, x, f, tracking, err)

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
        function fun(xx) result(fx)
            real(real64), intent(in), dimension(:) :: xx
            real(real64) :: fx
            fx = cfptr(n, xx)
        end function
    end function

! ------------------------------------------------------------------------------
    !> @brief Utilizes a Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm to 
    !! minimize a function of N variables.
    !!
    !! @param[in] fcn The function to minimize.
    !! @param[in] grad A function for evaluating the gradiant.  If null, the
    !!  gradient is estimated numerically.
    !! @param[in] n The number of variables.
    !! @param[in,out] x On input, an N-element array containing an initial
    !!  estimate to the solution.  On output, the updated solution estimate.
    !! @param[out] f On output, the value of the function at @p x.
    !! @param[in] cntrls The iteration controls.
    !! @param[in] ls The line search controls.
    !! @param[out] stats The iteration status.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !! - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !! - NL_INVALID_INPUT_ERROR: Occurs if @p x is not appropriately sized for
    !!      the problem as defined in @p fcn.
    !! - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !! - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
    !!      the allowed number of iterations.
    function c_solver_bfgs(fcn, grad, n, x, f, cntrls, ls, stats) &
            bind(C, name = "c_solver_bfgs") result(flag)
        ! Arguments
        type(c_funptr), intent(in), value :: fcn, grad
        integer(c_int), intent(in), value :: n
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(out) :: f
        type(iteration_controls), intent(in) :: cntrls
        type(line_search_controls), intent(in) :: ls
        type(iteration_process), intent(out) :: stats
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        type(fcnnvar_helper) :: obj
        procedure(c_fcnnvar), pointer :: cfptr
        procedure(fcnnvar), pointer :: fptr
        procedure(c_gradientfcn), pointer :: cgptr
        procedure(gradientfcn), pointer :: gptr
        type(iteration_behavior) :: tracking
        type(bfgs) :: solver
        type(line_search) :: search

        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)
        call c_f_procpointer(fcn, cfptr)
        fptr => fun
        call obj%set_fcn(fptr, n)
        if (c_associated(grad)) then
            call c_f_procpointer(grad, cgptr)
            gptr => grd
            call obj%set_gradient_fcn(gptr)
        end if

        ! Set the solver parameters
        call solver%set_max_fcn_evals(cntrls%max_function_evals)
        call solver%set_tolerance(cntrls%function_tolerance)
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
        call solver%solve(obj, x, f, tracking, err)

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
        function fun(xx) result(fx)
            real(real64), intent(in), dimension(:) :: xx
            real(real64) :: fx
            fx = cfptr(n, xx)
        end function

        subroutine grd(xx, gx)
            real(real64), intent(in), dimension(:) :: xx
            real(real64), intent(out), dimension(:) :: gx
            call cgptr(n, xx, gx)
        end subroutine
    end function

! ******************************************************************************
! POLYNOMIAL SUPPORT
! ------------------------------------------------------------------------------
    !> @brief Initializes a new C-compatible polynomial object.
    !!
    !! @param[in] order The order of the polynomial.  This must be at least 1.
    !! @param[out] poly A pointer to the polynomial object.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if a zero or negative polynomial order
    !!      was specified.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
    function c_init_polynomial(order, poly) &
            bind(C, name = "c_init_polynomial") result(flag)
        ! Arguments
        integer(c_int), intent(in), value :: order
        type(c_polynomial), intent(out) :: poly
        integer(c_int) :: flag

        ! Local Variables
        integer(int32) :: i
        type(errors) :: err
        type(polynomial) :: fpoly
        integer(int8), pointer, dimension(:) :: map
        integer(int8), allocatable, dimension(:) :: buffer

        ! Initialization
        call err%set_exit_on_error(.false.)
        flag = NL_NO_ERROR

        ! Initialize the Fortran polynomial
        call fpoly%initialize(order, err = err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if

        ! Initialize the pointer
        buffer = transfer(fpoly, buffer)
        allocate(map(size(buffer)))
        do i = 1, size(buffer)
            map(i) = buffer(i)
        end do

        ! Store the results
        poly%size_in_bytes = size(map)
        poly%ptr = c_loc(map)
    end function

! ------------------------------------------------------------------------------
    !> @brief Frees memory allocated for the polynomial object.
    !!
    !! @param[in,out] poly The polynomial object to free.
    subroutine c_free_polynomial(poly) bind(C, name = "c_free_polynomial")
        ! Arguments
        type(c_polynomial), intent(inout) :: poly

        ! Local Variables
        integer(int8), pointer, dimension(:) :: map

        ! Ensure there's something to work with
        if (.not.c_associated(poly%ptr) .or. poly%size_in_bytes == 0) return

        ! Obtain the pointer
        call c_f_pointer(poly%ptr, map, [poly%size_in_bytes])
        if (.not.associated(map)) return

        ! Free memory
        deallocate(map)

        ! Make the C pointer NULL
        poly%ptr = c_null_ptr
        poly%size_in_bytes = 0
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the order of the polynomial.
    !!
    !! @param[in] poly The polynomial object.
    !! @return The order of the polynomial.
    function c_get_polynomial_order(poly) &
            bind(C, name = "c_get_polynomial_order") result(n)
        ! Arguments
        type(c_polynomial), intent(in) :: poly
        integer(c_int) :: n

        ! Local Variables
        integer(int8), pointer, dimension(:) :: map
        type(polynomial) :: fpoly

        ! Initialization
        n = 0
        
        ! Ensure there's something to work with
        if (.not.c_associated(poly%ptr) .or. poly%size_in_bytes == 0) return

        ! Obtain the pointer
        call c_f_pointer(poly%ptr, map, [poly%size_in_bytes])
        if (.not.associated(map)) return

        ! Reconstruct the Fortran polynomial object
        fpoly = transfer(map, fpoly)

        ! Obtain the polynomial order
        n = fpoly%order()
    end function

! ------------------------------------------------------------------------------
    !> @brief Fits a data set to the polynomial.
    !!
    !! @param[in,out] poly The polynomial object.
    !! @param[in] n The number of data points to fit.
    !! @param[in] x An N-element array of the independent variable data points.
    !! @param[in] y An N-element array of the dependent variable data points.
    !! @param[in] zero Set to true to force the fit thru zero; else, set to
    !!  false.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if a zero or negative polynomial order
    !!      was specified, or if order is too large for the data set.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
    function c_fit_polynomial(poly, n, x, y, zero) &
            bind(C, name = "c_fit_polynomial") result(flag)
        ! Arguments
        type(c_polynomial), intent(inout) :: poly
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n), y(n)
        logical(c_bool), intent(in), value :: zero
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(int32) :: i, order
        integer(int8), pointer, dimension(:) :: map
        integer(int8), allocatable, dimension(:) :: buffer
        type(polynomial) :: fpoly
        real(real64), allocatable, dimension(:) :: ycopy

        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)
        ycopy = y   ! Prevents overwritting of y

        ! Ensure there's something to work with
        if (.not.c_associated(poly%ptr) .or. poly%size_in_bytes == 0) return

        ! Obtain the pointer
        call c_f_pointer(poly%ptr, map, [poly%size_in_bytes])
        if (.not.associated(map)) return

        ! Reconstruct the Fortran polynomial object
        fpoly = transfer(map, fpoly)

        ! Fit the data
        order = fpoly%order()
        if (zero) then
            call fpoly%fit_thru_zero(x, ycopy, order, err)
        else
            call fpoly%fit(x, ycopy, order, err)
        end if
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if

        ! Reconstruct the C polynomial type
        deallocate(map)
        buffer = transfer(fpoly, buffer)
        allocate(map(size(buffer)))
        do i = 1, size(buffer)
            map(i) = buffer(i)
        end do

        ! Store the results
        poly%size_in_bytes = size(map)
        poly%ptr = c_loc(map)
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets an array containing the polynomial coefficients in ascending
    !! order such that f(x) = c0 + c1 * x + c2 * x**2 .... + cN * x**N.
    !!
    !! @param[in] poly The polynomial object.
    !! @param[in] nc The number of elements in @p c.  Ideally, this value is
    !!  one greater than the order of the polynomial.
    !! @param[out] c An NC-element array where the coefficients will be written.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !! - NL_INVALID_INPUT_ERROR: Occurs if @p nc is smaller than one greater 
    !!      than the order of the polynomial.
    function c_get_polynomial_coefficients(poly, nc, c) &
            bind(C, name = "c_get_polynomial_coefficients") result(flag)
        ! Arguments
        type(c_polynomial), intent(in) :: poly
        integer(c_int), intent(in), value :: nc
        real(c_double), intent(out) :: c(nc)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(int32) :: i, order, n
        integer(int8), pointer, dimension(:) :: map
        type(polynomial) :: fpoly

        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Ensure there's something to work with
        if (.not.c_associated(poly%ptr) .or. poly%size_in_bytes == 0) return

        ! Obtain the pointer
        call c_f_pointer(poly%ptr, map, [poly%size_in_bytes])
        if (.not.associated(map)) return

        ! Reconstruct the Fortran polynomial object
        fpoly = transfer(map, fpoly)

        ! Get the order of the polynomial, and then copy over the coefficients
        order = fpoly%order()
        n = order + 1
        if (nc < n) then
            flag = NL_INVALID_INPUT_ERROR
            return
        end if
        do i = 1, n
            c(i) = fpoly%get(i)
        end do
    end function

! ------------------------------------------------------------------------------
    !> @brief Sets the coefficients of the polynomial by using an array 
    !! containing the polynomial coefficients in ascending order such that
    !! f(x) = c0 + c1 * x + c2 * x**2 .... + cN * x**N.
    !!
    !! @param[in,out] poly The polynomial object.
    !! @param[in] nc The number of elements in @P c.  This value must be 
    !!  one greater than the order of the polynomial.
    !! @param[in] c The NC-element array containing the new polynomial 
    !!  coefficients in ascending order.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !! - NL_INVALID_INPUT_ERROR: Occurs if @p nc is not equal to than one
    !!      greater than the order of the polynomial.
    function c_set_polynomial_coefficients(poly, nc, c) &
            bind(C, name = "c_set_polynomial_coefficients") result(flag)
        ! Arguments
        type(c_polynomial), intent(inout) :: poly
        integer(c_int), intent(in), value :: nc
        real(c_double) :: c(nc)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(int32) :: i, order, n
        integer(int8), pointer, dimension(:) :: map
        integer(int8), allocatable, dimension(:) :: buffer
        type(polynomial) :: fpoly

        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Ensure there's something to work with
        if (.not.c_associated(poly%ptr) .or. poly%size_in_bytes == 0) return

        ! Obtain the pointer
        call c_f_pointer(poly%ptr, map, [poly%size_in_bytes])
        if (.not.associated(map)) return

        ! Reconstruct the Fortran polynomial object
        fpoly = transfer(map, fpoly)

        ! Get the order of the polynomial, and then copy over the coefficients
        order = fpoly%order()
        n = order + 1
        if (nc /= n) then
            flag = NL_INVALID_INPUT_ERROR
            return
        end if
        do i = 1, n
            call fpoly%set(i, c(i))
        end do

        ! Reconstruct the C polynomial type
        deallocate(map)
        buffer = transfer(fpoly, buffer)
        allocate(map(size(buffer)))
        do i = 1, size(buffer)
            map(i) = buffer(i)
        end do

        ! Store the results
        poly%size_in_bytes = size(map)
        poly%ptr = c_loc(map)
    end function

! ------------------------------------------------------------------------------
    !> @brief Evaluates the polynomial at the specified values.
    !!
    !! @param[in] poly The polynomial object.
    !! @param[in] n The number of points at which to evaluate the polynomial.
    !! @param[in] x An N-element array containing the values at which to
    !!  evaluate the polynomial.
    !! @param[out] y An N-element array where the results of the polynomial
    !!  evaluation will be written.
    subroutine c_evaluate_polynomial_real(poly, n, x, y) &
            bind(C, name = "c_evaluate_polynomial_real")
        ! Arguments
        type(c_polynomial), intent(in) :: poly
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out) :: y(n)

        ! Local Variables
        integer(int8), pointer, dimension(:) :: map
        type(polynomial) :: fpoly

        ! Ensure there's something to work with
        if (.not.c_associated(poly%ptr) .or. poly%size_in_bytes == 0) return

        ! Obtain the pointer
        call c_f_pointer(poly%ptr, map, [poly%size_in_bytes])
        if (.not.associated(map)) return

        ! Reconstruct the Fortran polynomial object
        fpoly = transfer(map, fpoly)

        ! Evaluate the polynomial at X
        y = fpoly%evaluate(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Evaluates the polynomial at the specified values.
    !!
    !! @param[in] poly The polynomial object.
    !! @param[in] n The number of points at which to evaluate the polynomial.
    !! @param[in] x An N-element array containing the values at which to
    !!  evaluate the polynomial.
    !! @param[out] y An N-element array where the results of the polynomial
    !!  evaluation will be written.
    subroutine c_evaluate_polynomial_complex(poly, n, x, y) &
            bind(C, name = "c_evaluate_polynomial_complex")
        ! Arguments
        type(c_polynomial), intent(in) :: poly
        integer(c_int), intent(in), value :: n
        complex(c_double), intent(in) :: x(n)
        complex(c_double), intent(out) :: y(n)

        ! Local Variables
        integer(int8), pointer, dimension(:) :: map
        type(polynomial) :: fpoly

        ! Ensure there's something to work with
        if (.not.c_associated(poly%ptr) .or. poly%size_in_bytes == 0) return

        ! Obtain the pointer
        call c_f_pointer(poly%ptr, map, [poly%size_in_bytes])
        if (.not.associated(map)) return

        ! Reconstruct the Fortran polynomial object
        fpoly = transfer(map, fpoly)

        ! Evaluate the polynomial at X
        y = fpoly%evaluate(x)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the roots of the polynomial.
    !!
    !! @param[in] poly The polynomial object.
    !! @param[in] n THe size of @p rts.  This must be equal to the order of the
    !!  polynomial.
    !! @param[out] rts An N-element array where the roots of the polynomial will
    !!  be written.
    !!
    !! @return An error flag with the following possible values.
    !! - NL_NO_ERROR: No error has occurred - successful execution.
    !! - NL_INVALID_INPUT_ERROR: Occurs if @p n is not equal to the order of
    !!      the polynomial.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    function c_polynomial_roots(poly, n, rts) &
            bind(C, name = "c_polynomial_roots") result(flag)
        ! Arguments
        type(c_polynomial), intent(in) :: poly
        integer(c_int), intent(in), value :: n
        complex(c_double), intent(out) :: rts(n)
        integer(c_int) :: flag

        ! Local Variables
        type(errors) :: err
        integer(int8), pointer, dimension(:) :: map
        type(polynomial) :: fpoly
        integer(int32) :: order

        ! Initialization
        flag = NL_NO_ERROR
        call err%set_exit_on_error(.false.)

        ! Ensure there's something to work with
        if (.not.c_associated(poly%ptr) .or. poly%size_in_bytes == 0) return

        ! Obtain the pointer
        call c_f_pointer(poly%ptr, map, [poly%size_in_bytes])
        if (.not.associated(map)) return

        ! Reconstruct the Fortran polynomial object
        fpoly = transfer(map, fpoly)

        ! Ensure the output array is properly sized
        order = fpoly%order()
        if (n /= order) then
            flag = NL_INVALID_INPUT_ERROR
            return
        end if

        ! Compute the polynomial roots
        rts = fpoly%roots(err)
        if (err%has_error_occurred()) then
            flag = err%get_error_flag()
            return
        end if
    end function

! ------------------------------------------------------------------------------
end module
