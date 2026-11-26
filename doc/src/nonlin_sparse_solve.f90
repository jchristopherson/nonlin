module nonlin_sparse_solve
    use iso_fortran_env, only : real64
    use ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
    use nonlin_error_handling
    use nonlin_linesearch
    use nonlin_helper
    use nonlin_types
    use nonlin_multi_eqn_mult_var
    use nonlin_solve
    use linalg
    implicit none
    private
    public :: sparse_newton_solver

    type, extends(line_search_solver) :: sparse_newton_solver
        !! Defines a Newton solver for sparse systems of nonlinear equations. At
        !! each iteration the correction is obtained from the linearized system
        !! $$ J(x_k) \Delta x_k = -F(x_k), \qquad x_{k+1} = x_k + \Delta x_k $$
        !! where \(J(x_k)\) is the Jacobian matrix of \(F\). The linear solve
        !! uses preconditioned GMRES, and the optional line search uses the
        !! actual slope \(F^T J \Delta x\) for the inexact correction.
        type(csr_matrix), private :: jac
            !! The Jacobian matrix
        real(real64), private, allocatable, dimension(:) :: xold
            !! A workspace array used to store the previous solution values
        type(msr_matrix), private :: precond
            !! The ILU preconditioner
        integer, private, allocatable, dimension(:) :: ilupvt
            !! The pivot array from the incomplete LU factorization
        real(real64), private, allocatable, dimension(:) :: dir
            !! The solution direction vector workspace
    contains
        procedure, public :: solve => sns_solve
        procedure, private :: allocate_workspace => sns_alloc_workspaces
    end type

contains
! ******************************************************************************
! SPARSE_NEWTON_SOLVER
! ------------------------------------------------------------------------------
    subroutine sns_solve(this, fcn, x, fvec, ib, args)
        !! Applies Newton's method in conjunction with a backtracking type
        !! line search to solve N equations of N unknowns.
        !! The sparse Jacobian must be supplied by FCN. Linear-factorization,
        !! Krylov, nonfinite-value, and evaluation-budget failures return NaN
        !! residuals rather than terminating the program. As the merit function
        !! uses a sum of squares, finite residuals whose squared norm overflows
        !! are also reported as a failure.
        class(sparse_newton_solver), intent(inout) :: this
            !! The [[sparse_newton_solver]] object.
        class(vecfcn_helper), intent(inout) :: fcn
            !! The [[vecfcn_helper]] object containing the equations to solve.
        real(real64), intent(inout), dimension(:) :: x
            !! On input, an N-element array containing an initial estimate to
            !! the solution.  On output, the updated solution estimate.  N is
            !! the number of variables.
        real(real64), intent(out), dimension(:) :: fvec
            !! An N-element array that, on output, will contain the values of
            !! each equation as evaluated at the variable values given in x.
            !! If the iteration process does not converge, each element is
            !! returned as NaN, and, if supplied, ib will indicate a
            !! failure to converge.
        type(iteration_behavior), intent(out), optional :: ib
            !! An optional output, that if provided, allows the caller to
            !! obtain iteration performance statistics.
        class(*), intent(in), optional :: args
            !! An optional argument to allow the user to communicate with fcn.

        ! Parameters
        real(real64), parameter :: mintol = 1.0d-12
        real(real64), parameter :: factor = 1.0d2

        ! Local Variables
        logical :: xcnvrg, fcnvrg, check, cg
        integer :: i, nvar, neqn, neval, iter, njac, maxeval, flag, lflag, lsmax
        real(real64) :: ftol, xtol, rtol, f, fold, test, rnorm, xnorm, fnorm, nan, slope
        real(real64), allocatable :: jdir(:)
        class(line_search), allocatable :: ls
        type(iteration_behavior) :: lib

        ! Initialization
        xcnvrg = .false.
        fcnvrg = .false.
        neqn = fcn%get_equation_count()
        nvar = fcn%get_variable_count()
        neval = 0
        iter = 0
        njac = 0
        ftol = this%get_fcn_tolerance()
        xtol = this%get_var_tolerance()
        rtol = sqrt(epsilon(rtol))
        maxeval = this%get_max_fcn_evals()
        if (maxeval < 1) error stop NL_INVALID_INPUT_ERROR
        nan = ieee_value(nan, ieee_quiet_nan)
        if (present(ib)) then
            call ib%initialize(maxeval)
        end if
        if (this%get_use_line_search()) then
            if (.not.this%is_line_search_defined()) then
                call this%set_default_line_search()
            end if
            call this%get_line_search(ls)
            lsmax = ls%get_max_fcn_evals()
            if (lsmax < 1) error stop NL_INVALID_INPUT_ERROR
        end if

        ! Input Checking
        if (.not.fcn%is_fcn_defined()) error stop NL_UNDEFINED_FUNCTION_ERROR
        if (.not.fcn%get_use_sparse()) error stop NL_INVALID_OPERATION_ERROR
        if (.not.fcn%is_sparse_jacobian_defined()) error stop NL_UNDEFINED_FUNCTION_ERROR
        if (nvar < 1 .or. nvar /= neqn) error stop NL_INVALID_INPUT_ERROR
        flag = 0
        if (size(x) /= nvar) error stop NL_ARRAY_SIZE_ERROR
        if (size(fvec) /= neqn) error stop NL_ARRAY_SIZE_ERROR

        ! Allocate necessary workspaces
        call this%allocate_workspace(nvar)

        ! Test to see if the initial guess is a root
        call fcn%fcn(x, fvec, args = args)
        neval = 1
        flag = 0
        if (.not.all(ieee_is_finite(x)) .or. .not.all(ieee_is_finite(fvec))) then
            flag = 1
            go to 900
        end if
        f = 0.5d0 * dot_product(fvec, fvec)
        if (.not.ieee_is_finite(f)) then
            flag = 1
            go to 900
        end if
        test = 0.0d0
        do i = 1, neqn
            test = max(abs(fvec(i)), test)
        end do
        if (test < ftol) then
            fcnvrg = .true.
        end if

        ! Process
        flag = 0 ! Used to check for convergence errors
        if (.not.fcnvrg) then
            ! Main Iteration Loop
            do
                if (neval >= maxeval) then
                    flag = 1
                    exit
                end if
                ! Increment the iteration counter
                iter = iter + 1

                ! Compute the Jacobian
                call fcn%sparse_jacobian(x, this%jac, args = args)
                njac = njac + 1
                if (.not.allocated(this%jac%values)) then
                    flag = 1
                    exit
                end if
                if (.not.all(ieee_is_finite(this%jac%values)) .or. &
                    all(this%jac%values == 0.0d0)) &
                then
                    flag = 1
                    exit
                end if

                ! Solve the linear system
                this%xold = x
                fold = f

                call lu_factor(this%jac, this%precond, this%ilupvt, &
                    droptol = epsilon(rnorm), flag = lflag)
                if (lflag /= NL_NO_ERROR) then
                    flag = 1
                    exit
                end if
                call pgmres_solve(this%jac, this%precond, this%ilupvt, &
                    -fvec, this%dir, tol = rtol, maxits = 200, flag = lflag)
                if (lflag /= NL_NO_ERROR) then
                    flag = 1
                    exit
                end if
                jdir = matmul(this%jac, this%dir)
                if (.not.all(ieee_is_finite(jdir))) then
                    flag = 1
                    exit
                end if
                rnorm = norm2(fvec + jdir) / norm2(fvec)
                if (.not.ieee_is_finite(rnorm) .or. rnorm > rtol) then
                    flag = 1
                    exit
                end if

                ! Apply the line search, if needed
                if (this%get_use_line_search()) then
                    slope = dot_product(fvec, jdir)
                    if (.not.ieee_is_finite(slope) .or. slope >= 0.0d0) then
                        flag = 1
                        exit
                    end if
                    call ls%set_max_fcn_evals(min(lsmax, maxeval - neval))
                    call ls%search_directional(fcn, this%xold, this%dir, slope, &
                        x, fvec, fold = fold, fx = f, ib = lib, args = args)
                    neval = neval + lib%fcn_count
                    if (.not.ieee_is_finite(f)) then
                        ! The line search failed to locate an acceptable point
                        flag = 1
                        exit
                    end if
                else
                    ! No line search - just update the solution estimate
                    x = x + this%dir
                    call fcn%fcn(x, fvec, args = args)
                    f = 0.5d0 * dot_product(fvec, fvec)
                    neval = neval + 1
                end if

                if (.not.all(ieee_is_finite(x)) .or. .not.all(ieee_is_finite(fvec)) .or. &
                    .not.ieee_is_finite(f)) then
                    flag = 1
                    exit
                end if

                ! Check for convergence
                call test_convergence(x, this%xold, fvec, fvec, .false., xtol, &
                    ftol, rtol, check, xcnvrg, fcnvrg, cg, xnorm, fnorm)
                if (present(ib)) then
                    call ib%log_iteration(fnorm, xnorm, nan)
                end if
                if (check) then
                    ! The solution has converged
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


    900     continue
        ! Report out iteration statistics
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = njac
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = .false.
        end if

        ! Check for convergence issues
        if (flag /= 0) then
            ! The solver failed to converge - signal via NaN.  The
            ! iteration_behavior flags above already denote non-convergence.
            fvec = nan
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine sns_alloc_workspaces(this, n)
        !! Allocates internal workspaces.
        class(sparse_newton_solver), intent(inout) :: this
            !! The [[sparse_newton_solver]] object.
        integer, intent(in) :: n
            !! The number of equations.

        if (allocated(this%xold)) then
            if (size(this%xold) /= n) then
                deallocate(this%xold)
                allocate(this%xold(n))
            end if
        else
            allocate(this%xold(n))
        end if

        if (allocated(this%ilupvt)) then
            if (size(this%ilupvt) /= n) then
                deallocate(this%ilupvt)
                allocate(this%ilupvt(n))
            end if
        else
            allocate(this%ilupvt(n))
        end if

        if (allocated(this%dir)) then
            if (size(this%dir) /= n) then
                deallocate(this%dir)
                allocate(this%dir(n))
            end if
        else
            allocate(this%dir(n))
        end if
    end subroutine

! ------------------------------------------------------------------------------
end module