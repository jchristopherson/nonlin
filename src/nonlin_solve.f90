! nonlin_solve.f90

!> @brief \b nonlin_solve
!!
!! @par Purpose
!! To provide various routines capapble of solving systems of nonlinear 
!! equations.
module nonlin_solve
    use linalg_constants, only : dp, i32
    use nonlin_types
    use nonlin_linesearch, only : line_search
    use ferror, only : errors
    use linalg_factor, only : qr_factor, form_qr, qr_rank1_update
    use linalg_core, only : rank1_update, mtx_mult
    use linalg_solve, only : solve_triangular_system
    implicit none
    private
    public :: equation_solver
    public :: line_search_solver
    public :: quasi_newton_solver
    public :: nonlin_solver

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief A base class for various nonlinear equation solvers.
    type, abstract :: equation_solver
        private
        !> The maximum number of function evaluations allowed per solve.
        integer(i32) :: m_maxEval = 100
        !> The convergence criteria on function values.
        real(dp) :: m_fcnTol = 1.0d-8
        !> The convergence criteria on change in variable values.
        real(dp) :: m_xtol = 1.0d-12
        !> The convergence criteria for the slope of the gradient vector.
        real(dp) :: m_gtol = 1.0d-12
        !> Set to true to print iteration status; else, false.
        logical :: m_printStatus = .false.
    contains
        !> @brief Gets the maximum number of function evaluations allowed during
        !! a single solve.
        procedure, public :: get_max_fcn_evals => es_get_max_eval
        !> @brief Sets the maximum number of function evaluations allowed during
        !! a single solve.
        procedure, public :: set_max_fcn_evals => es_set_max_eval
        !> @brief Gets the convergence on function value tolerance.
        procedure, public :: get_fcn_tolerance => es_get_fcn_tol
        !> @brief Sets the convergence on function value tolerance.
        procedure, public :: set_fcn_tolerance => es_set_fcn_tol
        !> @brief Gets the convergence on change in variable tolerance.
        procedure, public :: get_var_tolerance => es_get_var_tol
        !> @brief Sets the convergence on change in variable tolerance.
        procedure, public :: set_var_tolerance => es_set_var_tol
        !> @brief Gets the convergence on slope of the gradient vector 
        !! tolerance.
        procedure, public :: get_gradient_tolerance => es_get_grad_tol
        !> @brief Sets the convergence on slope of the gradient vector 
        !! tolerance.
        procedure, public :: set_gradient_tolerance => es_set_grad_tol
        !> @brief Gets a logical value determining if iteration status should be
        !! printed.
        procedure, public :: get_print_status => es_get_print_status
        !> @brief Sets a logical value determining if iteration status should be
        !! printed.
        procedure, public :: set_print_status => es_set_print_status
        !> @brief Solves the system of equations.
        procedure(nonlin_solver), deferred, public, pass :: solve
    end type

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
        procedure, public :: get_line_search => lss_get_line_search
        !> @brief Sets the line search module.
        procedure, public :: set_line_search => lss_set_line_search
        !> @brief Establishes a default line_search object for the line search 
        !! module.
        procedure, public :: set_default_line_search => lss_set_default
        !> @brief Tests to see if a line search module is defined.
        procedure, public :: is_line_search_defined => &
            lss_is_line_search_defined
        !> @brief Gets a value determining if a line-search should be employed.
        procedure, public :: get_use_line_search => lss_get_use_search
        !> @brief Sets a value determining if a line-search should be employed.
        procedure, public :: set_use_line_search => lss_set_use_search
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a quasi-Newton type solver based upon Broyden's method.
    type, extends(line_search_solver) :: quasi_newton_solver
        private
        !> The number of iterations that may pass between Jacobian calculation.
        integer(i32) :: m_jDelta = 5
    contains
        !> @brief Solves the system of equations.
        procedure, public :: solve => qns_solve
        !> @brief Gets the number of iterations that may pass before forcing a 
        !! recalculation of the Jacobian matrix.
        procedure, public :: get_jacobian_interval => qns_get_jac_interval
        !> @brief Sets the number of iterations that may pass before forcing a 
        !! recalculation of the Jacobian matrix.
        procedure, public :: set_jacobian_interval => qns_set_jac_interval
    end type

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        !> @brief Describes the interface of a nonlinear equation solver.
        !!
        !! @param[in,out] this The equation_solver-based object.
        !! @param[in] fcn The vecfcn_helper object containing the equations to
        !!  solve.
        !! @param[in,out] x On input, an N-element array containing an initial
        !!  estimate to the solution.  On output, the updated solution estimate.
        !!  N is the number of variables.
        !! @param[out] fvec An M-element array that, on output, will contain
        !!  the values of each equation as evaluated at the variable values
        !!  given in @p x.
        !! @param[out] ib An optional output, that if provided, allows the 
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can
        !!  be used to retrieve information relating to any errors encountered 
        !!  during execution.  If not provided, a default implementation of the
        !!  errors class is used internally to provide error handling.  The
        !!  possible error codes returned will likely vary from solver to 
        !!  solver.
        subroutine nonlin_solver(this, fcn, x, fvec, ib, err)
            use linalg_constants, only : dp, i32
            use nonlin_types, only : vecfcn_helper, iteration_behavior
            use ferror, only : errors
            import equation_solver
            class(equation_solver), intent(inout) :: this
            class(vecfcn_helper), intent(in) :: fcn
            real(dp), intent(inout), dimension(:) :: x
            real(dp), intent(out), dimension(:) :: fvec
            type(iteration_behavior), optional :: ib
            class(errors), intent(in), optional, target :: err
        end subroutine
    end interface

contains
! ******************************************************************************
! EQUATION_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Gets the maximum number of function evaluations allowed during a
    !! single solve.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The maximum number of function evaluations.
    pure function es_get_max_eval(this) result(n)
        class(equation_solver), intent(in) :: this
        integer(i32) :: n
        n = this%m_maxEval
    end function

! --------------------
    !> @brief Sets the maximum number of function evaluations allowed during a
    !! single solve.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] n The maximum number of function evaluations.
    subroutine es_set_max_eval(this, n)
        class(equation_solver), intent(inout) :: this
        integer(i32), intent(in) :: n
        this%m_maxEval = n
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on function value tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    pure function es_get_fcn_tol(this) result(x)
        class(equation_solver), intent(in) :: this
        real(dp) :: x
        x = this%m_fcnTol
    end function

! --------------------
    !> @brief Sets the convergence on function value tolerance.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] x The tolerance value.
    subroutine es_set_fcn_tol(this, x)
        class(equation_solver), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_fcnTol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on change in variable tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    pure function es_get_var_tol(this) result(x)
        class(equation_solver), intent(in) :: this
        real(dp) :: x
        x = this%m_xtol
    end function

! --------------------
    !> @brief Sets the convergence on change in variable tolerance.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] x The tolerance value.
    subroutine es_set_var_tol(this, x)
        class(equation_solver), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_xtol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on slope of the gradient vector tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    pure function es_get_grad_tol(this) result(x)
        class(equation_solver), intent(in) :: this
        real(dp) :: x
        x = this%m_gtol
    end function

! --------------------
    !> @brief Sets the convergence on slope of the gradient vector tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    subroutine es_set_grad_tol(this, x)
        class(equation_solver), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_gtol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets a logical value determining if iteration status should be
    !! printed.
    !!
    !! @param[in] this The equation_solver object.
    !! @return True if the iteration status should be printed; else, false.
    pure function es_get_print_status(this) result(x)
        class(equation_solver), intent(in) :: this
        logical :: x
        x = this%m_printStatus
    end function

! --------------------
    !> @brief Sets a logical value determining if iteration status should be
    !! printed.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] x True if the iteration status should be printed; else, false.
    subroutine es_set_print_status(this, x)
        class(equation_solver), intent(inout) :: this
        logical, intent(in) :: x
        this%m_printStatus = x
    end subroutine

! ******************************************************************************
! LINE_SEARCH_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Gets the line search module.
    !!
    !! @param[in] this The line_search_solver object.
    !! @param[out] The line_search object.
    subroutine lss_get_line_search(this, ls)
        class(line_search_solver), intent(in) :: this
        class(line_search), intent(out), allocatable :: ls
        if (allocated(this%m_lineSearch)) &
            allocate(ls, source = this%m_lineSearch)
    end subroutine

! ----------------------
    !> @brief Sets the line search module.
    !!
    !! @param[in,out] this The line_search_solver object.
    !! @param[in] ls The line_search object.
    subroutine lss_set_line_search(this, ls)
        class(line_search_solver), intent(inout) :: this
        class(line_search), intent(in) :: ls
        if (allocated(this%m_lineSearch)) deallocate(this%m_lineSearch)
        allocate(this%m_lineSearch, source = ls)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Establishes a default line_search object for the line search 
    !! module.
    !!
    !! @param[in,out] this The line_search_solver object.
    subroutine lss_set_default(this)
        class(line_search_solver), intent(inout) :: this
        type(line_search) :: ls
        call this%set_line_search(ls)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Tests to see if a line search module is defined.
    !!
    !! @param[in] this The line_search_solver object.
    !! @return Returns true if a module is defined; else, false.
    pure function lss_is_line_search_defined(this) result(x)
        class(line_search_solver), intent(in) :: this
        logical :: x
        x = allocated(this%m_lineSearch)
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets a value determining if a line-search should be employed.
    !!
    !! @param[in] this The line_search_solver object.
    !! @return Returns true if a line search should be used; else, false.
    pure function lss_get_use_search(this) result(x)
        class(line_search_solver), intent(in) :: this
        logical :: x
        x = this%m_useLineSearch
    end function

! --------------------
    !> @brief Sets a value determining if a line-search should be employed.
    !!
    !! @param[in,out] this The line_search_solver object.
    !! @param[in] x Set to true if a line search should be used; else, false.
    subroutine lss_set_use_search(this, x)
        class(line_search_solver), intent(inout) :: this
        logical, intent(in) :: x
        this%m_useLineSearch = x
    end subroutine

! ******************************************************************************
! QUASI_NEWTON_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Applies the quasi-Newton's method developed by Broyden in
    !! conjunction with a backtracking type line search to solve N equations
    !! of N unknowns.
    !!
    !! @param[in,out] this The equation_solver-based object.
    !! @param[in] fcn The vecfcn_helper object containing the equations to
    !!  solve.
    !! @param[in,out] x On input, an N-element array containing an initial
    !!  estimate to the solution.  On output, the updated solution estimate.
    !!  N is the number of variables.
    !! @param[out] fvec An M-element array that, on output, will contain
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
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Broyden%27s_method)
    !! - [Numerical Recipes](http://numerical.recipes/)
    subroutine qns_solve(this, fcn, x, fvec, ib, err)
        ! Arguments
        class(quasi_newton_solver), intent(inout) :: this
        class(vecfcn_helper), intent(in) :: fcn
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: fvec
        type(iteration_behavior), optional :: ib
        class(errors), intent(in), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: half = 0.5d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: factor = 1.0d2

        ! Local Variables
        logical :: restart, xcnvrg, fcnvrg, gcnvrg, check
        integer(i32) :: i, neqn, nvar, flag, lw1, lw2, lw3, neval, iter, &
            maxeval, jcount, njac
        real(dp), allocatable, dimension(:) :: work, tau, dx, df, fvold, &
            xold, s
        real(dp), allocatable, dimension(:,:) :: q, r, b
        real(dp) :: test, f, fold, alpha, temp, ftol, xtol, gtol, eps, &
            stpmax, x2, xnorm, fnorm
        type(iteration_behavior) :: lib
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        class(line_search), allocatable :: ls

        ! Initialization
        restart = .true.
        xcnvrg = .false.
        fcnvrg = .false.
        gcnvrg = .false.
        neqn = fcn%get_equation_count()
        nvar = fcn%get_variable_count()
        neval = 0
        iter = 0
        njac = 0
        ftol = this%get_fcn_tolerance()
        xtol = this%get_var_tolerance()
        gtol = this%get_gradient_tolerance()
        maxeval = this%get_max_fcn_evals()
        eps = epsilon(eps)
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = njac
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = gcnvrg
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        if (this%get_use_line_search()) then
            if (.not.this%is_line_search_defined()) &
                call this%set_default_line_search()
            call this%get_line_search(ls)
        end if

        ! Input Check
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("qns_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (nvar /= neqn) then
            ! ERROR: # of equations doesn't match # of variables
            write(errmsg, '(AI0AI0A)') "The number of equations (", neqn, &
                ") does not match the number of unknowns (", nvar, ")."
            call errmgr%report_error("qns_solve", trim(errmsg), &
                NL_INVALID_INPUT_ERROR)
            return
        end if
        flag = 0
        if (size(x) /= nvar) then
            flag = 3
        else if (size(fvec) /= neqn) then
            flag = 4
        end if
        if (flag /= 0) then
            ! One of the input arrays is not sized correctly
            write(errmsg, '(AI0A)') "Input number ", flag, &
                " is not sized correctly."
            call errmgr%report_error("qns_solve", trim(errmsg), &
                NL_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(q(neqn, neqn), stat = flag)
        if (flag == 0) allocate(r(neqn, nvar), stat = flag)
        if (flag == 0) allocate(tau(min(neqn, nvar)), stat = flag)
        if (flag == 0) allocate(b(neqn, nvar), stat = flag)
        if (flag == 0) allocate(df(neqn), stat = flag)
        if (flag == 0) allocate(fvold(neqn), stat = flag)
        if (flag == 0) allocate(xold(nvar), stat = flag)
        if (flag == 0) allocate(dx(nvar), stat = flag)
        if (flag == 0) allocate(s(neqn), stat = flag)
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("qns_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if
        call qr_factor(r, tau, work, lw1)
        call form_qr(r, tau, q, work, lw2)
        call fcn%jacobian(x, b, fv = fvec, olwork = lw3)
        allocate(work(max(lw1, lw2, lw3)), stat = flag)
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("qns_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Test to see if the initial guess is a root
        call fcn%fcn(x, fvec)
        f = half * dot_product(fvec, fvec)
        neval = neval + 1
        test = zero
        do i = 1, neqn
            test = max(abs(fvec(i)), test)
        end do
        if (test < ftol) then
            fcnvrg = .true.
        end if

        ! Process
        flag = 0 ! Used to check for convergence errors
        if (.not.fcnvrg) then
            ! Determine the maximum line search step
            stpmax = factor * max(norm2(x), real(nvar, dp))

            ! Main Iteration Loop
            do
                ! Update the iteration counter
                iter = iter + 1

                ! Compute or update the Jacobian
                if (restart) then
                    ! Compute the Jacobian
                    call fcn%jacobian(x, b, fvec, work)
                    njac = njac + 1

                    ! Compute the QR factorization, and form Q & R
                    r = b ! Copy the Jacobian - we'll need it later
                    call qr_factor(r, tau, work)
                    call form_qr(r, tau, q, work)

                    ! Reset the Jacobian iteration counter
                    jcount = 0
                else
                    ! Apply the rank 1 update to Q and R
                    df = fvec - fvold
                    dx = x - xold
                    x2 = dot_product(dx, dx)
                    alpha = one / x2

                    ! Compute S = ALPHA * (DF - B * DX)
                    s = alpha * (df - matmul(b, dx))

                    ! Compute the new Q and R matrices for the rank1 update:
                    ! B' = B + ALPHA * S * DX**T
                    call rank1_update(alpha, s, dx, b)
                    call qr_rank1_update(q, r, s, dx, work) ! S & DX overwritten

                    ! Increment the counter tracking how many iterations have
                    ! passed since the last Jacobian recalculation
                    jcount = jcount + 1
                end if

                ! Compute GRAD = B**T * F, store in DX
                call mtx_mult(.true., one, b, fvec, zero, dx)

                ! Store FVEC and X
                xold = x
                fvold = fvec
                fold = f

                ! Solve the linear system: B * DX = -F for DX noting that
                ! B = Q * R.  As such, form -Q**T * F, and store in DF
                call mtx_mult(.true., -one, q, fvec, zero, df)

                ! Now we have R * DX = -Q**T * F, and since R is upper 
                ! triangular, the solution is readily computed.  The solution
                ! will be stored in the first NVAR elements of DF
                call solve_triangular_system(.true., .false., .true., r, &
                    df(1:nvar))
                
                ! Ensure the new solution estimate is heading in a sensible
                ! direction.  If not, it is likely time to update the Jacobian
                temp = dot_product(dx, df(1:nvar))
                if (temp >= zero) then
                    restart = .true.
                    if (this%get_print_status()) then
                        call print_status(iter, neval, njac, xnorm, fnorm)
                    end if
                    cycle
                end if
                
                ! Apply the line search if needed
                if (this%get_use_line_search()) then
                    ! Define the step length for the line search
                    temp = dot_product(df(1:nvar), df(1:nvar))
                    if (temp > stpmax) df(1:nvar) = df(1:nvar) * (stpmax / temp)

                    ! Apply the line search
                    call ls%search(fcn, xold, dx, df(1:nvar), x, fvec, fold, &
                        f, lib, errmgr)
                    neval = neval + lib%fcn_count
                else
                    ! No line search - just update the solution estimate
                    x = x + df(1:nvar)
                    call fcn%fcn(x, fvec)
                    f = half * dot_product(fvec, fvec)
                    neval = neval + 1
                end if

                ! Test for convergence
                if (lib%converge_on_zero_diff .and. &
                        this%get_use_line_search()) then
                    call test_convergence(x, xold, fvec, dx, .true., xtol, &
                        ftol, gtol, check, xcnvrg, fcnvrg, gcnvrg, xnorm, fnorm)
                else
                    call test_convergence(x, xold, fvec, dx, .false., xtol, &
                        ftol, gtol, check, xcnvrg, fcnvrg, gcnvrg, xnorm, fnorm)
                end if
                if (.not.check) then
                    ! The solution did not converge, figure out why
                    if (gcnvrg) then
                        ! The slope of the gradient is sufficiently close to 
                        ! zero to cause issue.
                        if (restart) then
                            ! We've already tried recalculating a new Jacobian
                            exit
                        else
                            ! Try computing a new Jacobian
                            restart = .true.
                        end if
                    else
                        ! We have not converged, but we're not stuck with a 
                        ! zero slope gradient vector either.  Go ahead and
                        ! continue the iteration process without recomputing
                        ! the Jacobian - unless the user dictates a 
                        ! recaclulation.
                        if (jcount >= this%m_jDelta) then
                            restart = .true.
                        else
                            restart = .false.
                        end if
                    end if
                else
                    ! The solution has converged.  It's OK to exit
                    exit
                end if

                ! Print status
                if (this%get_print_status()) then
                    call print_status(iter, neval, njac, xnorm, fnorm)
                end if

                ! Ensure we haven't made too many function evaluations
                if (neval >= maxeval) then
                    flag = 1
                    exit
                end if
            end do
        end if

        ! Report out iteration statistics
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = njac
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = gcnvrg
        end if

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, '(AI0AE8.3AE8.3)') "The algorithm failed to " // &
                "converge.  Function evaluations performed: ", neval, &
                "." // new_line('c') // "Change in Variable: ", xnorm, &
                new_line('c') // "Residual: ", fnorm
            call errmgr%report_error("qns_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the number of iterations that may pass before forcing a 
    !! recalculation of the Jacobian matrix.
    !!
    !! @param[in] this The quasi_newton_solver object.
    !! @return The number of iterations.
    pure function qns_get_jac_interval(this) result(n)
        class(quasi_newton_solver), intent(in) :: this
        integer(i32) :: n
        n = this%m_jDelta
    end function

! --------------------
    !> @brief Sets the number of iterations that may pass before forcing a 
    !! recalculation of the Jacobian matrix.
    !!
    !! @param[in,out] this The quasi_newton_solver object.
    !! @param[in] n The number of iterations.
    subroutine qns_set_jac_interval(this, n)
        class(quasi_newton_solver), intent(inout) :: this
        integer(i32), intent(in) :: n
        this%m_jDelta = n
    end subroutine

! ------------------------------------------------------------------------------




! ******************************************************************************
! GENERAL ROUTINES
! ------------------------------------------------------------------------------
    !> @breif Tests for convergence.
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
    !! @param[out] cs True if convergence occurred due to change in variable.
    !! @param[out] cf True if convergence occurred due to residual.
    !! @param[out] cg True if convergence occured due to slope of the gradient.
    !! @param[out] xnorm The largest magnitude component of the scaled change 
    !!  in variable vector.
    !! @param[out] fnorm The largest magnitude residual component
    subroutine test_convergence(x, xo, f, g, lg, xtol, ftol, gtol, c, cx, cf, &
            cg, xnorm, fnorm)
        ! Arguments
        real(dp), intent(in), dimension(:) :: x, xo, f, g
        real(dp), intent(in) :: xtol, ftol, gtol
        logical, intent(in) :: lg
        logical, intent(out) :: c, cx, cf, cg
        real(dp), intent(out) :: xnorm, fnorm

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: half = 0.5d0

        ! Local Variables
        integer(i32) :: i, nvar, neqn
        real(dp) :: test, dxmax, fc, den

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
        print '(AE8.3AE8.3)', "FNORM: ", fnorm, ", FTOL: ", ftol
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
        print '(AE8.3AE8.3)', "XNORM: ", xnorm, ", XTOL: ", xtol
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
    !> @brief Prints the iteration status.
    !!
    !! @param[in] iter The iteration number.
    !! @param[in] nfeval The number of function evaluations.
    !! @param[in] njaceval The number of Jacobian evaluations.
    !! @param[in] xnorm The change in variable value.
    !! @param[in] fnorm The residual.
    subroutine print_status(iter, nfeval, njaceval, xnorm, fnorm)
        ! Arguments
        integer(i32), intent(in) :: iter, nfeval, njaceval
        real(dp), intent(in) :: xnorm, fnorm

        ! Process
        print *, ""
        print '(AI0)', "Iteration: ", iter
        print '(AI0)', "Function Evaluations: ", nfeval
        if (njaceval > 0) print '(AI0)', "Jacobian Evaluations: ", njaceval
        print '(AE8.3)', "Change in Variable: ", xnorm
        print '(AE8.3)', "Residual: ", fnorm
    end subroutine

! ------------------------------------------------------------------------------
end module
