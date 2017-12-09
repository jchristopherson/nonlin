! nonlin_solve.f90

!> @brief \b nonlin_solve
!!
!! @par Purpose
!! To provide various routines capapble of solving systems of nonlinear
!! equations.
module nonlin_solve
    use nonlin_types
    use nonlin_linesearch, only : line_search, limit_search_vector
    use ferror, only : errors
    use linalg_factor, only : qr_factor, form_qr, qr_rank1_update, lu_factor
    use linalg_core, only : rank1_update, mtx_mult, recip_mult_array
    use linalg_solve, only : solve_triangular_system, solve_lu
    implicit none
    private
    public :: line_search_solver
    public :: quasi_newton_solver
    public :: newton_solver
    public :: brent_solver

! ******************************************************************************
! TYPES
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

! ------------------------------------------------------------------------------
    !> @brief Defines a Newton solver.
    type, extends(line_search_solver) :: newton_solver
    contains
        !> @brief Solves the system of equations.
        procedure, public :: solve => ns_solve
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a solver based upon Brent's method for solving an equation
    !! of one variable without using derivatives.
    type, extends(equation_solver_1var) :: brent_solver
    contains
        !> @brief Solves the equation.
        procedure, public :: solve => brent_solve
    end type
 

contains
! ******************************************************************************
! LINE_SEARCH_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Gets the line search module.
    !!
    !! @param[in] this The line_search_solver object.
    !! @param[out] ls The line_search object.
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
    !! @param[out] fvec An N-element array that, on output, will contain
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
    !!  - NL_SPURIOUS_CONVERGENCE_ERROR: Occurs as a warning if the slope of the
    !!      gradient vector becomes sufficiently close to zero.
    !!
    !! @par Usage
    !! The following code provides an example of how to solve a system of N
    !! equations of N unknonwns using this Quasi-Newton method.
    !! @code{.f90}
    !! program main
    !!     use linalg_constants, only : dp
    !!     use nonlin_types, only : vecfcn, vecfcn_helper
    !!     use nonlin_solve, only : quasi_newton_solver
    !!
    !!     type(vecfcn_helper) :: obj
    !!     procedure(vecfcn), pointer :: fcn
    !!     type(quasi_newton_solver) :: solver
    !!     real(dp) :: x(2), f(2)
    !!
    !!     ! Set the initial conditions to [1, 1]
    !!     x = 1.0d0
    !!
    !!     ! Define the function
    !!     fcn => fcn1
    !!     call obj%set_fcn(fcn, 2, 2)
    !!
    !!     ! Solve the system of equations.  The solution overwrites X
    !!     call solver%solve(obj, x, f)
    !!
    !!     ! Print the output and the residual:
    !!     print '(AF5.3AF5.3A)', "The solution: (", x(1), ", ", x(2), ")"
    !!     print '(AE8.3AE8.3A)', "The residual: (", f(1), ", ", f(2), ")"
    !! contains
    !!     ! System of Equations:
    !!     !
    !!     ! x**2 + y**2 = 34
    !!     ! x**2 - 2 * y**2 = 7
    !!     !
    !!     ! Solution:
    !!     ! x = +/-5
    !!     ! y = +/-3
    !!     subroutine fcn1(x, f)
    !!         real(dp), intent(in), dimension(:) :: x
    !!         real(dp), intent(out), dimension(:) :: f
    !!         f(1) = x(1)**2 + x(2)**2 - 34.0d0
    !!         f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    !!     end subroutine
    !! end program
    !! @endcode
    !! The above program returns the following results.
    !! @code{.txt}
    !! The solution: (5.000, 3.000)
    !! The residual: (.604E-10, .121E-09)
    !! @endcode
    !!
    !! @par See Also
    !! - [Broyden's Paper](http://www.ams.org/journals/mcom/1965-19-092/S0025-5718-1965-0198670-6/S0025-5718-1965-0198670-6.pdf)
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Broyden%27s_method)
    !! - [Numerical Recipes](http://numerical.recipes/)
    subroutine qns_solve(this, fcn, x, fvec, ib, err)
        ! Arguments
        class(quasi_newton_solver), intent(inout) :: this
        class(vecfcn_helper), intent(in) :: fcn
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: fvec
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

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
        real(dp) :: test, f, fold, temp, ftol, xtol, gtol, &
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
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = njac
            ib%gradient_count = 0
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

                    ! Compute S = ALPHA * (DF - B * DX)
                    s = (df - matmul(b, dx))
                    call recip_mult_array(x2, s)

                    ! Compute the new Q and R matrices for the rank1 update:
                    ! B' = B + ALPHA * S * DX**T
                    call rank1_update(one, s, dx, b)
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
                    call limit_search_vector(df(1:nvar), stpmax)
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
                            ! We've already tried recalculating a new Jacobian,
                            ! issue a warning
                            write(errmsg, '(AI0AE8.3AE8.3)') &
                                "It appears the solution has settled to " // &
                                "a point where the slope of the gradient " // &
                                "is effectively zero.  " // new_line('c') // &
                                "Function evaluations performed: ", neval, &
                                "." // new_line('c') // &
                                "Change in Variable: ", xnorm, &
                                new_line('c') // "Residual: ", fnorm
                            call errmgr%report_warning("nqs_solve", &
                                trim(errmsg), NL_SPURIOUS_CONVERGENCE_ERROR)
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
            ib%gradient_count = 0
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

! ******************************************************************************
! NEWTON_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Applies Newton's method in conjunction with a backtracking type 
    !! line search to solve N equations of N unknowns.
    !!
    !! @param[in,out] this The equation_solver-based object.
    !! @param[in] fcn The vecfcn_helper object containing the equations to
    !!  solve.
    !! @param[in,out] x On input, an N-element array containing an initial
    !!  estimate to the solution.  On output, the updated solution estimate.
    !!  N is the number of variables.
    !! @param[out] fvec An N-element array that, on output, will contain
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
    !!  - NL_SPURIOUS_CONVERGENCE_ERROR: Occurs as a warning if the slope of the
    !!      gradient vector becomes sufficiently close to zero.
    !!
    !! @par Usage
    !! The following code provides an example of how to solve a system of N
    !! equations of N unknonwns using Newton's method.
    !! @code{.f90}
    !! program main
    !!     use linalg_constants, only : dp
    !!     use nonlin_types, only : vecfcn, vecfcn_helper
    !!     use nonlin_solve, only : newton_solver
    !!
    !!     type(vecfcn_helper) :: obj
    !!     procedure(vecfcn), pointer :: fcn
    !!     type(newton_solver) :: solver
    !!     real(dp) :: x(2), f(2)
    !!
    !!     ! Set the initial conditions to [1, 1]
    !!     x = 1.0d0
    !!
    !!     ! Define the function
    !!     fcn => fcn1
    !!     call obj%set_fcn(fcn, 2, 2)
    !!
    !!     ! Solve the system of equations.  The solution overwrites X
    !!     call solver%solve(obj, x, f)
    !!
    !!     ! Print the output and the residual:
    !!     print '(AF5.3AF5.3A)', "The solution: (", x(1), ", ", x(2), ")"
    !!     print '(AE8.3AE8.3A)', "The residual: (", f(1), ", ", f(2), ")"
    !! contains
    !!     ! System of Equations:
    !!     !
    !!     ! x**2 + y**2 = 34
    !!     ! x**2 - 2 * y**2 = 7
    !!     !
    !!     ! Solution:
    !!     ! x = +/-5
    !!     ! y = +/-3
    !!     subroutine fcn1(x, f)
    !!         real(dp), intent(in), dimension(:) :: x
    !!         real(dp), intent(out), dimension(:) :: f
    !!         f(1) = x(1)**2 + x(2)**2 - 34.0d0
    !!         f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    !!     end subroutine
    !! end program
    !! @endcode
    !! The above program returns the following results.
    !! @code{.txt}
    !! The solution: (5.000, 3.000)
    !! The residual: (.000E+00, .000E+00)
    !! @endcode
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Newton%27s_method)
    subroutine ns_solve(this, fcn, x, fvec, ib, err)
        ! Arguments
        class(newton_solver), intent(inout) :: this
        class(vecfcn_helper), intent(in) :: fcn
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: fvec
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: half = 0.5d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: mintol = 1.0d-12
        real(dp), parameter :: factor = 1.0d2

        ! Local Variables
        logical :: check, xcnvrg, fcnvrg, gcnvrg
        integer(i32) :: i, neqn, nvar, lwork, flag, neval, iter, maxeval, njac
        integer(i32), allocatable, dimension(:) :: ipvt
        real(dp), allocatable, dimension(:) :: dir, grad, xold, work
        real(dp), allocatable, dimension(:,:) :: jac
        real(dp) :: ftol, xtol, gtol, f, fold, stpmax, xnorm, fnorm, temp, test
        type(iteration_behavior) :: lib
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        class(line_search), allocatable :: ls

        ! Initialization
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
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = njac
            ib%gradient_count = 0
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

        ! Input Checking
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("ns_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (nvar /= neqn) then
            ! ERROR: # of equations doesn't match # of variables
            write(errmsg, '(AI0AI0A)') "The number of equations (", neqn, &
                ") does not match the number of unknowns (", nvar, ")."
            call errmgr%report_error("ns_solve", trim(errmsg), &
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
            call errmgr%report_error("ns_solve", trim(errmsg), &
                NL_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(ipvt(nvar), stat = flag)
        if (flag == 0) allocate(dir(nvar), stat = flag)
        if (flag == 0) allocate(grad(nvar), stat = flag)
        if (flag == 0) allocate(xold(nvar), stat = flag)
        if (flag == 0) allocate(jac(nvar, neqn), stat = flag)
        if (flag == 0) then
            call fcn%jacobian(x, jac, fv = fvec, olwork = lwork)
            allocate(work(lwork), stat = flag)
        end if
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("ns_solve", &
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
            ! Compute the maximum step size for the line search process
            stpmax = factor * max(norm2(x), real(nvar, dp))

            ! Main Iteration Loop
            do
                ! Increment the iteration counter
                iter = iter + 1

                ! Compute the Jacobian
                call fcn%jacobian(x, jac, fvec, work)
                njac = njac + 1

                ! Compute the gradient
                do i = 1, nvar
                    grad(i) = dot_product(jac(:,i), fvec)
                end do

                ! Compute the LU factorization of the Jacobian
                call lu_factor(jac, ipvt, errmgr)
                if (errmgr%has_warning_occurred()) then
                    ! The Jacobian is singular - warning was issued already, so
                    ! simply exit the routine.  Do not return as a return at 
                    ! this point would not allow for proper updating of the
                    ! iteration tracking parameters
                    exit
                end if

                ! Store previous iteration values
                xold = x
                fold = f

                ! Define the right-hand-side for the linear system
                dir = -fvec

                ! Solve the linear system of equations
                call solve_lu(jac, ipvt, dir)

                ! Apply the line search if needed
                if (this%get_use_line_search()) then
                    ! Define the step length for the line search
                    temp = dot_product(dir, dir)
                    if (temp > stpmax) dir = dir * (stpmax / temp)

                    ! Apply the line search
                    call limit_search_vector(dir, stpmax)
                    call ls%search(fcn, xold, grad, dir, x, fvec, &
                        fold, f, lib, errmgr)
                    neval = neval + lib%fcn_count
                else
                    ! No line search - just update the solution estimate
                    x = x + dir
                    call fcn%fcn(x, fvec)
                    f = half * dot_product(fvec, fvec)
                    neval = neval + 1
                end if

                ! Check for convergence
                call test_convergence(x, xold, fvec, grad, .true., xtol, &
                    ftol, gtol, check, xcnvrg, fcnvrg, gcnvrg, xnorm, fnorm)
                if (check) then
                    ! The solution has converged
                    exit
                else if (gcnvrg) then
                    ! The solution appears to have settled at a point where
                    ! the gradient has a zero slope
                    write(errmsg, '(AI0AE8.3AE8.3)') &
                        "It appears the solution has settled to " // &
                        "a point where the slope of the gradient " // &
                        "is effectively zero.  " // new_line('c') // &
                        "Function evaluations performed: ", neval, &
                        "." // new_line('c') // &
                        "Change in Variable: ", xnorm, &
                        new_line('c') // "Residual: ", fnorm
                    call errmgr%report_warning("ns_solve", trim(errmsg), &
                        NL_SPURIOUS_CONVERGENCE_ERROR)
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
            ib%gradient_count = 0
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
            call errmgr%report_error("ns_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if
    end subroutine

! ******************************************************************************
! BRENT_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Solves the equation.
    !!
    !! @param[in,out] this The brent_solver object.
    !! @param[in] fcn The fcn1var_helper object containing the equation
    !!  to solve.
    !! @param[in,out] x A parameter used to return the solution.  Notice, any
    !!  input value will be ignored as this routine relies upon the
    !!  search limits in @p lim to provide a starting point.
    !! @param[in] lim A value_pair object defining the search limits.
    !! @param[out] f An optional parameter used to return the function
    !!  residual as computed at @p x.
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
    !!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
    !!      the allowed number of iterations.
    !!
    !! @par Usage
    !! The following code provides an example of how to solve an equation of
    !! one variable using Brent's method.
    !! @code{.f90}
    !! program main
    !!     use linalg_constants, only : dp
    !!     use nonlin_types, only : fcn1var, fcn1var_helper, value_pair
    !!     use nonlin_solve, only : brent_solver
    !!
    !!     type(fcn1var_helper) :: obj
    !!     procedure(fcn1var), pointer :: fcn
    !!     type(brent_solver) :: solver
    !!     real(dp) :: x, f
    !!     type(value_pair) :: limits
    !!
    !!     ! Define the solution limits
    !!     lmiits%x1 = 1.5d0
    !!     limits%x2 = 5.0d0
    !!
    !!     ! Define the function
    !!     fcn => fcn1
    !!     call obj%set_fcn(fcn)
    !!
    !!     ! Solve the equation
    !!     call solver%solve(obj, x, limits, f)
    !!
    !!     ! Print the output and the residual:
    !!     print '(AF5.3)', "The solution: ", x
    !!     print '(AE9.3)', "The residual: ", f
    !! contains
    !!     ! f(x) = sin(x) / x, SOLUTION: x = n * pi for n = 1, 2, 3, ...
    !!     function fcn1(x) result(f)
    !!         real(dp), intent(in) :: x
    !!         real(dp) :: f
    !!         f = sin(x) / x
    !!     end function
    !! end program
    !! @endcode
    !! The above program returns the following results.
    !! @code{.txt}
    !! The solution: 3.142
    !! The residual: -.751E-11
    !! @endcode
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Brent%27s_method)
    !! - [Numerical Recipes](http://numerical.recipes/)
    !! - R.P. Brent, "Algorithms for Minimization without Derivatives,"
    !!      Dover Publications, January 2002. ISBN 0-486-41998-3.
    !!      Further information available 
    !!      [here](https://maths-people.anu.edu.au/~brent/pub/pub011.html).
    subroutine brent_solve(this, fcn, x, lim, f, ib, err)
        ! Arguments
        class(brent_solver), intent(inout) :: this
        class(fcn1var_helper), intent(in) :: fcn
        real(dp), intent(inout) :: x
        type(value_pair), intent(in) :: lim
        real(dp), intent(out), optional :: f
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: half = 0.5d0
        real(dp), parameter :: one = 1.0d0
        real(dp), parameter :: two = 2.0d0
        real(dp), parameter :: three = 3.0d0

        ! Local Variables
        logical :: fcnvrg, xcnvrg
        integer(i32) :: neval, maxeval, flag, iter
        real(dp) :: ftol, xtol, a, b, c, fa, fb, fc, p, q, r, s, xm, e, d, &
            mn1, mn2, eps, tol1, temp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg

        ! Initialization
        fcnvrg = .false.
        xcnvrg = .false.
        x = zero
        a = min(lim%x1, lim%x2)
        b = max(lim%x1, lim%x2)
        neval = 0
        iter = 0
        eps = epsilon(eps)
        ftol = this%get_fcn_tolerance()
        xtol = this%get_var_tolerance()
        maxeval = this%get_max_fcn_evals()
        if (present(f)) f = zero
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = 0
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = .false.
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        
        ! Input Check
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("brent_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (abs(a - b) < eps) then
            ! ERROR: Search limits are too tight
            write(errmsg, '(AE8.3AE8.3)') "Search limits have no " // &
                "appreciable difference between them.  Lower Limit: ", a, &
                ", Upper Limit: ", b
            call errmgr%report_error("brent_solve", trim(errmsg), &
                NL_INVALID_OPERATION_ERROR)
            return
        end if

        ! Process
        flag = 0
        fa = fcn%fcn(a)
        fb = fcn%fcn(b)
        neval = 2
        fc = fb
        do
            ! Increment the iteration counter
            iter = iter + 1

            ! Adjust the bounding interval
            if ((fb > zero .and. fc >= zero) .or. &
                    (fb < zero .and. fc < zero)) then
                c = a
                fc = fa
                d = b - a
                e = d
            end if
            if (abs(fc) < abs(fb)) then
                a = b
                b = c
                c = a
                fa = fb
                fb = fc
                fc = fa
            end if

            ! Convergence Check
            tol1 = two * eps * abs(b) + half * xtol
            xm = half * (c - b)
            if (abs(fb) < ftol) then
                x = b
                fcnvrg = .true.
                exit
            end if
            if (abs(xm) <= tol1) then
                x = b
                xcnvrg = .true.
                exit
            end if

            ! Actual Method
            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
                ! Attempt the inverse quadratic interpolation to determine
                ! the root
                s = fb / fa
                if (abs(a - c) < eps) then ! a == c
                    p = two * xm * s
                    q = one - s
                else
                    q = fa / fc
                    r = fb / fc
                    p = s * (two * xm * q * (q - r) - (b - a) * (r - one))
                    q = (q - one) * (r - one) * (s - one)
                end if

                ! Ensure we're within bounds
                if (p > zero) q = -q
                p = abs(p)
                mn1 = three * xm * q - abs(tol1 * q)
                mn2 = abs(e * q)
                if (mn1 < mn2) then
                    temp = mn1
                else
                    temp = mn2
                end if
                if (two * p < temp) then
                    ! Accept the interpolation
                    e = d
                    d = p / q
                else
                    ! The interpolation failed, use bisection
                    d = xm
                    e = d
                end if
            else
                ! The bounds are decreasing too slowly, use bisection
                d = xm
                e = d
            end if

            ! Move the last best guess to the lower limit parameter (A)
            a = b
            fa = fb
            if (abs(d) > tol1) then
                b = b + d
            else
                b = b + sign(tol1, xm)
            end if
            fb = fcn%fcn(b)
            neval = neval + 1

            ! Print iteration status
            if (this%get_print_status()) then
                call print_status(iter, neval, 0, xm, fb)
            end if

            ! Ensure we haven't made too many function evaluations
            if (neval >= maxeval) then
                flag = 1
                exit
            end if
        end do

        ! Report out iteration statistics and other optional outputs
        if (present(f)) f = fb
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = 0
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = .false.
        end if

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, '(AI0AE8.3AE8.3)') "The algorithm failed to " // &
                "converge.  Function evaluations performed: ", neval, &
                "." // new_line('c') // "Change in Variable: ", xm, &
                new_line('c') // "Residual: ", fb
            call errmgr%report_error("brent_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if
    end subroutine


! ******************************************************************************
! GENERAL ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Tests for convergence.
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
    !! @param[out] cx True if convergence occurred due to change in variable.
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
end module
