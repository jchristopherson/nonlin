module nonlin_solve
    use iso_fortran_env
    use nonlin_error_handling
    use nonlin_multi_eqn_mult_var
    use nonlin_single_var
    use nonlin_linesearch
    use nonlin_helper
    use nonlin_types
    use ferror
    use linalg, only : qr_factor, form_qr, qr_rank1_update, lu_factor, &
        rank1_update, mtx_mult, recip_mult_array, solve_triangular_system, &
        solve_lu
    implicit none
    private
    public :: line_search_solver
    public :: quasi_newton_solver
    public :: newton_solver
    public :: brent_solver
    public :: newton_1var_solver

    type, abstract, extends(equation_solver) :: line_search_solver
        !! A class describing nonlinear solvers that use a line search
        !! algorithm to improve convergence behavior.
        class(line_search), private, allocatable :: m_lineSearch
            !! The line search module.
        logical, private :: m_useLineSearch = .true.
            !! Set to true if a line search should be used regardless of the 
            !! status of m_lineSearch
    contains
        procedure, public :: get_line_search => lss_get_line_search
        procedure, public :: set_line_search => lss_set_line_search
        procedure, public :: set_default_line_search => lss_set_default
        procedure, public :: is_line_search_defined => &
            lss_is_line_search_defined
        procedure, public :: get_use_line_search => lss_get_use_search
        procedure, public :: set_use_line_search => lss_set_use_search
    end type

    type, extends(line_search_solver) :: quasi_newton_solver
        !! Defines a quasi-Newton type solver based upon Broyden's method.
        integer(int32), private :: m_jDelta = 5
            !! The number of iterations that may pass between Jacobian
            !! calculation.
    contains
        procedure, public :: solve => qns_solve
        procedure, public :: get_jacobian_interval => qns_get_jac_interval
        procedure, public :: set_jacobian_interval => qns_set_jac_interval
    end type
    
    type, extends(line_search_solver) :: newton_solver
        !! Defines a Newton solver.
    contains
        procedure, public :: solve => ns_solve
    end type

    type, extends(equation_solver_1var) :: brent_solver
        !! Defines a solver based upon Brent's method for solving an equation
        !! of one variable without using derivatives.
    contains
        procedure, public :: solve => brent_solve
    end type

    type, extends(equation_solver_1var) :: newton_1var_solver
        !! Defines a solver based upon Newtons's method for solving an
        !! equation of one variable.  The algorithm uses a bisection method in
        !! conjunction with Newton's method in order to keep bounds upon the
        !! Newton iterations.
    contains
        procedure, public :: solve => newt1var_solve
    end type
    
contains
! ******************************************************************************
! LINE_SEARCH_SOLVER
! ------------------------------------------------------------------------------
    subroutine lss_get_line_search(this, ls)
        !! Gets the line search module.
        class(line_search_solver), intent(in) :: this
            !! The [[line_search_solver]] object.
        class(line_search), intent(out), allocatable :: ls
            !! The [[line_search]] object.
        if (allocated(this%m_lineSearch)) &
            allocate(ls, source = this%m_lineSearch)
    end subroutine

! ----------------------
    subroutine lss_set_line_search(this, ls)
        !! Sets the line search module.
        class(line_search_solver), intent(inout) :: this
            !! The [[line_search_solver]] object.
        class(line_search), intent(in) :: ls
            !! The [[line_search]] object.
        if (allocated(this%m_lineSearch)) deallocate(this%m_lineSearch)
        allocate(this%m_lineSearch, source = ls)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine lss_set_default(this)
        !! Establishes a default line_search object for the line search
        !! module.
        class(line_search_solver), intent(inout) :: this
            !! The [[line_search_solver]] object.
        type(line_search) :: ls
        call this%set_line_search(ls)
    end subroutine

! ------------------------------------------------------------------------------
    pure function lss_is_line_search_defined(this) result(x)
        !! Tests to see if a line search module is defined.
        class(line_search_solver), intent(in) :: this
            !! The [[line_search_solver]] object.
        logical :: x
            !! Returns true if a module is defined; else, false.
        x = allocated(this%m_lineSearch)
    end function

! ------------------------------------------------------------------------------
    pure function lss_get_use_search(this) result(x)
        !! Gets a value determining if a line-search should be employed.
        class(line_search_solver), intent(in) :: this
            !! The [[line_search_solver]] object.
        logical :: x
            !! Returns true if a line search should be used; else, false.
        x = this%m_useLineSearch
    end function

! --------------------
    subroutine lss_set_use_search(this, x)
        !! Sets a value determining if a line-search should be employed.
        class(line_search_solver), intent(inout) :: this
            !! The [[line_search_solver]] object.
        logical, intent(in) :: x
            !! Set to true if a line search should be used; else, false.
        this%m_useLineSearch = x
    end subroutine

! ******************************************************************************
! QUASI_NEWTON_SOLVER
! ------------------------------------------------------------------------------
    subroutine qns_solve(this, fcn, x, fvec, ib, args, err)
        !! Applies the quasi-Newton's method developed by Broyden in 
        !! conjunction with a backtracking type line search to solve N equations
        !! of N unknowns.
        !!
        !! See Also:
        !!
        !! - [Broyden's Paper](http://www.ams.org/journals/mcom/1965-19-092/S0025-5718-1965-0198670-6/S0025-5718-1965-0198670-6.pdf)
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Broyden%27s_method)
        !!
        !! - [Numerical Recipes](http://numerical.recipes/)
        class(quasi_newton_solver), intent(inout) :: this
            !! The [[quasi_newton_solver]] object.
        class(vecfcn_helper), intent(in) :: fcn
            !! The [[vecfcn_helper]] object containing the equations to solve.
        real(real64), intent(inout), dimension(:) :: x
            !! On input, an N-element array containing an initial estimate to 
            !! the solution.  On output, the updated solution estimate.  N is 
            !! the number of variables.
        real(real64), intent(out), dimension(:) :: fvec
            !! An N-element array that, on output, will contain the values of 
            !! each equation as evaluated at the variable values given in x.
        type(iteration_behavior), optional :: ib
            !! An optional output, that if provided, allows the caller to 
            !! obtain iteration performance statistics.
        class(*), intent(inout), optional :: args
            !! An optional argument to allow the user to communicate with
            !! fcn.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: half = 0.5d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: factor = 1.0d2

        ! Local Variables
        logical :: restart, xcnvrg, fcnvrg, gcnvrg, check
        integer(int32) :: i, neqn, nvar, flag, lw1, lw2, lw3, neval, iter, &
            maxeval, jcount, njac
        real(real64), allocatable, dimension(:) :: work, tau, dx, df, fvold, &
            xold, s
        real(real64), allocatable, dimension(:,:) :: q, r, b
        real(real64) :: test, f, fold, temp, ftol, xtol, gtol, &
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
            write(errmsg, 100) "The number of equations (", neqn, &
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
            write(errmsg, 101) "Input number ", flag, &
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
        call fcn%jacobian(x, b, fv = fvec, olwork = lw3, args = args)
        allocate(work(max(lw1, lw2, lw3)), stat = flag)
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("qns_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Test to see if the initial guess is a root
        call fcn%fcn(x, fvec, args)
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
            stpmax = factor * max(norm2(x), real(nvar, real64))

            ! Main Iteration Loop
            do
                ! Update the iteration counter
                iter = iter + 1

                ! Compute or update the Jacobian
                if (restart) then
                    ! Compute the Jacobian
                    call fcn%jacobian(x, b, fvec, work, args = args)
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
                        f, lib, args = args, err = errmgr)
                    neval = neval + lib%fcn_count
                else
                    ! No line search - just update the solution estimate
                    x = x + df(1:nvar)
                    call fcn%fcn(x, fvec, args)
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
                            write(errmsg, 102) &
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
            write(errmsg, 102) "The algorithm failed to " // &
                "converge.  Function evaluations performed: ", neval, &
                "." // new_line('c') // "Change in Variable: ", xnorm, &
                new_line('c') // "Residual: ", fnorm
            call errmgr%report_error("qns_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if

        ! Format
100     format(A, I0, A, I0, A)
101     format(A, I0, A)
102     format(A, I0, A, E10.3, A, E10.3)
    end subroutine

! ------------------------------------------------------------------------------
    pure function qns_get_jac_interval(this) result(n)
        !! Gets the number of iterations that may pass before forcing a
        !! recalculation of the Jacobian matrix.
        class(quasi_newton_solver), intent(in) :: this
            !! The [[quasi_newton_solver]] object.
        integer(int32) :: n
            !! The number of iterations.
        n = this%m_jDelta
    end function

! --------------------
    subroutine qns_set_jac_interval(this, n)
        !! Sets the number of iterations that may pass before forcing a
        !! recalculation of the Jacobian matrix.
        class(quasi_newton_solver), intent(inout) :: this
            !! The [[quasi_newton_solver]] object.
        integer(int32), intent(in) :: n
            !! The number of iterations.
        this%m_jDelta = n
    end subroutine

! ******************************************************************************
! NEWTON_SOLVER
! ------------------------------------------------------------------------------
    subroutine ns_solve(this, fcn, x, fvec, ib, args, err)
        !! Applies Newton's method in conjunction with a backtracking type
        !! line search to solve N equations of N unknowns.
        class(newton_solver), intent(inout) :: this
            !! The [[newton_solver]] object.
        class(vecfcn_helper), intent(in) :: fcn
            !! The [[vecfcn_helper]] object containing the equations to solve.
        real(real64), intent(inout), dimension(:) :: x
            !! On input, an N-element array containing an initial estimate to 
            !! the solution.  On output, the updated solution estimate.  N is 
            !! the number of variables.
        real(real64), intent(out), dimension(:) :: fvec
            !! An N-element array that, on output, will contain the values of 
            !! each equation as evaluated at the variable values given in x.
        type(iteration_behavior), optional :: ib
            !! An optional output, that if provided, allows the caller to 
            !! obtain iteration performance statistics.
        class(*), intent(inout), optional :: args
            !! An optional argument to allow the user to communicate with fcn.
        class(errors), intent(inout), optional, target :: err
            !! An error-handling object.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: half = 0.5d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: mintol = 1.0d-12
        real(real64), parameter :: factor = 1.0d2

        ! Local Variables
        logical :: check, xcnvrg, fcnvrg, gcnvrg
        integer(int32) :: i, neqn, nvar, lwork, flag, neval, iter, maxeval, njac
        integer(int32), allocatable, dimension(:) :: ipvt
        real(real64), allocatable, dimension(:) :: dir, grad, xold, work
        real(real64), allocatable, dimension(:,:) :: jac
        real(real64) :: ftol, xtol, gtol, f, fold, stpmax, xnorm, fnorm, temp, test
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
            write(errmsg, 100) "The number of equations (", neqn, &
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
            write(errmsg, 101) "Input number ", flag, &
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
            call fcn%jacobian(x, jac, fv = fvec, olwork = lwork, args = args)
            allocate(work(lwork), stat = flag)
        end if
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("ns_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Test to see if the initial guess is a root
        call fcn%fcn(x, fvec, args)
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
            stpmax = factor * max(norm2(x), real(nvar, real64))

            ! Main Iteration Loop
            do
                ! Increment the iteration counter
                iter = iter + 1

                ! Compute the Jacobian
                call fcn%jacobian(x, jac, fvec, work, args = args)
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
                        fold, f, lib, args = args, err = errmgr)
                    neval = neval + lib%fcn_count
                else
                    ! No line search - just update the solution estimate
                    x = x + dir
                    call fcn%fcn(x, fvec, args)
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
                    write(errmsg, 102) &
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
            write(errmsg, 102) "The algorithm failed to " // &
                "converge.  Function evaluations performed: ", neval, &
                "." // new_line('c') // "Change in Variable: ", xnorm, &
                new_line('c') // "Residual: ", fnorm
            call errmgr%report_error("ns_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if

        ! Formatting
100     format(A, I0, A, I0, A)
101     format(A, I0, A)
102     format(A, I0, A, E10.3, A, E10.3)
    end subroutine

! ******************************************************************************
! BRENT_SOLVER
! ------------------------------------------------------------------------------
    subroutine brent_solve(this, fcn, x, lim, f, ib, args, err)
        !! Solves an equation of one variable using Brent's method.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Brent%27s_method)
        !!
        !! - [Numerical Recipes](http://numerical.recipes/)
        !!
        !! - R.P. Brent, "Algorithms for Minimization without Derivatives,"
        !!      Dover Publications, January 2002. ISBN 0-486-41998-3.
        !!      Further information available
        !!      [here](https://maths-people.anu.edu.au/~brent/pub/pub011.html).
        class(brent_solver), intent(inout) :: this
            !! The [[brent_solver]] object.
        class(fcn1var_helper), intent(in) :: fcn
            !! The [[fcn1var_helper]] object containing the equation to solve.
        real(real64), intent(inout) :: x
            !! A parameter used to return the solution.  Notice, any input 
            !! value will be ignored as this routine relies upon the search 
            !! limits in lim to provide a starting point.
        type(value_pair), intent(in) :: lim
            !! A [[value_pair]] object defining the search limits.
        real(real64), intent(out), optional :: f
            !! An optional parameter used to return the function residual as 
            !! computed at x.
        type(iteration_behavior), optional :: ib
            !! An optional output, that if provided, allows the caller to 
            !! obtain iteration performance statistics.
        class(*), intent(inout), optional :: args
            !! An optional argument to allow the user to communicate with fcn.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: half = 0.5d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: two = 2.0d0
        real(real64), parameter :: three = 3.0d0

        ! Local Variables
        logical :: fcnvrg, xcnvrg
        integer(int32) :: neval, maxeval, flag, iter
        real(real64) :: ftol, xtol, a, b, c, fa, fb, fc, p, q, r, s, xm, e, d, &
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
            write(errmsg, 100) "Search limits have no " // &
                "appreciable difference between them.  Lower Limit: ", a, &
                ", Upper Limit: ", b
            call errmgr%report_error("brent_solve", trim(errmsg), &
                NL_INVALID_OPERATION_ERROR)
            return
        end if

        ! Process
        flag = 0
        fa = fcn%fcn(a, args)
        fb = fcn%fcn(b, args)
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
            fb = fcn%fcn(b, args)
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
            write(errmsg, 101) "The algorithm failed to " // &
                "converge.  Function evaluations performed: ", neval, &
                "." // new_line('c') // "Change in Variable: ", xm, &
                new_line('c') // "Residual: ", fb
            call errmgr%report_error("brent_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if

        ! Formatting
100     format(A, E10.3, A, E10.3)
101     format(A, I0, A, E10.3, A, E10.3)
    end subroutine

! ******************************************************************************
! NEWTON_1VAR_SOLVER
! ------------------------------------------------------------------------------
    subroutine newt1var_solve(this, fcn, x, lim, f, ib, args, err)
        !! Solves an equation of one variable using Newton's method.
        class(newton_1var_solver), intent(inout) :: this
            !! The [[newton_1var_solver]] object.
        class(fcn1var_helper), intent(in) :: fcn
            !! The [[fcn1var_helper]] object containing the equation to solve.
        real(real64), intent(inout) :: x
            !! A parameter used to return the solution.  Notice, any input 
            !! value will be ignored as this routine relies upon the search 
            !! limits in lim to provide a starting point.
        type(value_pair), intent(in) :: lim
            !! A value_pair object defining the search limits.
        real(real64), intent(out), optional :: f
            !! An optional parameter used to return the function residual as 
            !! computed at x.
        type(iteration_behavior), optional :: ib
            !! An optional output, that if provided, allows the caller to 
            !! obtain iteration performance statistics.
        class(*), intent(inout), optional :: args
            !! An optional argument to allow the user to communicate with fcn.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: p5 = 0.5d0
        real(real64), parameter :: two = 2.0d0

        ! Local Variables
        logical :: fcnvrg, xcnvrg, dcnvrg
        integer(int32) :: neval, ndiff, maxeval, flag, iter
        real(real64) :: ftol, xtol, dtol, xh, xl, fh, fl, x1, x2, eps, dxold, &
            dx, df, temp, ff
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg

        ! Initialization
        fcnvrg = .false.
        xcnvrg = .false.
        dcnvrg = .false.
        neval = 0
        ndiff = 0
        iter = 0
        ftol = this%get_fcn_tolerance()
        xtol = this%get_var_tolerance()
        dtol = this%get_diff_tolerance()
        maxeval = this%get_max_fcn_evals()
        if (present(f)) f = zero
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = ndiff
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = dcnvrg
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        x1 = min(lim%x1, lim%x2)
        x2 = max(lim%x1, lim%x2)
        eps = epsilon(eps)

        ! Input Check
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("brent_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (abs(x1 - x2) < eps) then
            ! ERROR: Search limits are too tight
            write(errmsg, 100) "Search limits have no " // &
                "appreciable difference between them.  Lower Limit: ", x1, &
                ", Upper Limit: ", x2
            call errmgr%report_error("brent_solve", trim(errmsg), &
                NL_INVALID_OPERATION_ERROR)
            return
        end if

        ! See if the root is one of the end points
        flag = 0
        fl = fcn%fcn(x1, args)
        fh = fcn%fcn(x2, args)
        neval = 2
        if (abs(fl) < ftol) then
            x = x1
            if (present(f)) f = fl
            if (present(ib)) then
                ib%converge_on_fcn = .true.
                ib%fcn_count = 2
            end if
            return
        end if
        if (abs(fh) < ftol) then
            x = x2
            if (present(f)) f = fh
            if (present(ib)) then
                ib%converge_on_fcn = .true.
                ib%fcn_count = 2
            end if
            return
        end if

        ! Process
        if (fl < zero) then
            xl = x1
            xh = x2
        else
            xl = x2
            xh = x1
        end if
        x = p5 * (x1 + x2)
        dxold = abs(x2 - x1)
        dx = dxold
        ff = fcn%fcn(x, args)
        df = fcn%diff(x, f = ff, args = args)
        neval = neval + 1
        ndiff = ndiff + 1
        do
            ! Increment the iteration counter
            iter = iter + 1

            ! Bisect if the Newton step went out of range, or if the rate
            ! of change was too slow
            if ((((x - xh) * df - ff) * ((x - xl) * df - ff) > zero) .or. &
                (abs(two * ff) > abs(dxold * df))) &
            then
                ! Bisection
                dxold = dx
                dx = p5 * (xh - xl)
                x = xl + dx
                if (abs(xl - x) < xtol) then
                    ! Convergence as the change in root is within tolerance
                    xcnvrg = .true.
                    exit
                end if
            else
                ! Newton's Method
                dxold = dx
                dx = ff / df
                temp = x
                x = x - dx
                if (abs(temp - x) < xtol) then
                    ! Convergence as the change in root is within tolerance
                    xcnvrg = .true.
                    exit
                end if
            end if

            ! Update function values
            ff = fcn%fcn(x, args)
            df = fcn%diff(x, f = ff, args = args)
            neval = neval + 1
            ndiff = ndiff + 1

            ! Check for convergence
            if (abs(ff) < ftol) then
                fcnvrg = .true.
                exit
            end if
            if (abs(dx) < xtol) then
                xcnvrg = .true.
                exit
            end if
            if (abs(df) < dtol) then
                dcnvrg = .true.
                exit
            end if

            ! Update the bracket on the root
            if (ff < zero) then
                xl = x
            else
                xh = x
            end if

            ! Print status
            if (this%get_print_status()) then
                call print_status(iter, neval, ndiff, dx, ff)
            end if

            ! Ensure we haven't made too many function evaluations
            if (neval >= maxeval) then
                flag = 1
                exit
            end if
        end do

        ! Ensure the function value is current with the estimate of the root
        if (present(f)) then
            f = fcn%fcn(x, args)
            neval = neval + 1
        end if

        ! Report out iteration statistics and other optional outputs
        if (present(f)) f = ff
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = ndiff
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = dcnvrg
        end if

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, 101) "The algorithm failed to " // &
                "converge.  Function evaluations performed: ", neval, &
                "." // new_line('c') // "Root estimate: ", x, &
                new_line('c') // "Residual: ", ff
            call errmgr%report_error("newt1var_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if

        ! Format
100     format(A, E10.3, A, E10.3)
101     format(A, I0, A, E10.3, A, E10.3)
    end subroutine

end module
