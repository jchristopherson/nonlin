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
    use sparskit, only : ilut, pgmres
    implicit none
    private
    public :: sparse_newton_solver

    type, extends(line_search_solver) :: sparse_newton_solver
        !! Defines a Newton solver for sparse systems of nonlinear equations. At
        !! each iteration the correction is obtained from the linearized system
        !! $$ J(x_k) \Delta x_k = -F(x_k), \qquad x_{k+1} = x_k + \Delta x_k $$
        !! where \(J(x_k)\) is the Jacobian matrix of \(F\). The linear solve
        !! uses ILUT-preconditioned GMRES as implemented in SPARSKIT, and the
        !! optional line search uses the actual slope \(F^T J \Delta x\) for the
        !! inexact correction. All workspaces are owned by this type.
        type(csr_matrix), private :: jac
            !! The Jacobian matrix
        real(real64), private, allocatable, dimension(:) :: xold
            !! A workspace array used to store the previous solution values
        integer, private, allocatable, dimension(:) :: ilupvt
            !! The ILUT diagonal pointer array (JU)
        real(real64), private, allocatable, dimension(:) :: alu
            !! The ILUT factor values
        integer, private, allocatable, dimension(:) :: jlu
            !! The ILUT factor column indices
        real(real64), private, allocatable, dimension(:) :: w
            !! ILUT workspace of length N+1
        integer, private, allocatable, dimension(:) :: jw
            !! ILUT integer workspace of length 2N
        real(real64), private, allocatable, dimension(:) :: dir
            !! The solution direction vector workspace
        real(real64), private, allocatable, dimension(:) :: jdir
            !! The product of the Jacobian and the solution direction
        real(real64), private, allocatable, dimension(:) :: rhs
            !! The right-hand side workspace for the GMRES solve
        real(real64), private, allocatable, dimension(:,:) :: vv
            !! The GMRES Krylov basis workspace of size N x (IM+1)
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
        integer :: i, k, nvar, neqn, neval, iter, njac, maxeval, flag, lflag, lsmax, lfil, im
        real(real64) :: ftol, xtol, rtol, f, fold, test, rnorm, xnorm, fnorm, nan, slope
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
                if (.not.sns_is_valid_csr(this%jac, nvar)) then
                    flag = 1
                    exit
                end if
                if (.not.all(ieee_is_finite(this%jac%values)) .or. &
                    all(this%jac%values == 0.0d0)) &
                then
                    flag = 1
                    exit
                end if

                ! Size the workspaces from the Jacobian's row fill
                lfil = 1
                do i = 1, nvar
                    lfil = max(lfil, this%jac%row_indices(i+1) - this%jac%row_indices(i))
                end do
                call this%allocate_workspace(nvar, lfil)

                ! Solve the linear system
                this%xold = x
                fold = f

                call ilut(nvar, this%jac%values, this%jac%column_indices, &
                    this%jac%row_indices, lfil, epsilon(rnorm), this%alu, this%jlu, &
                    this%ilupvt, size(this%alu), this%w, this%jw, lflag)
                if (lflag /= 0) then
                    flag = 1
                    exit
                end if
                im = size(this%vv, 2) - 1
                this%dir = 0.0d0
                this%rhs = -fvec
                call pgmres(nvar, im, this%rhs, this%dir, this%vv, rtol, 200, 0, &
                    this%jac%values, this%jac%column_indices, this%jac%row_indices, &
                    this%alu, this%jlu, this%ilupvt, lflag)
                if (lflag /= 0 .and. lflag /= -1) then
                    flag = 1
                    exit
                end if

                ! Compute the Jacobian-direction product
                this%jdir = 0.0d0
                do i = 1, nvar
                    do k = this%jac%row_indices(i), this%jac%row_indices(i+1) - 1
                        this%jdir(i) = this%jdir(i) + this%jac%values(k) * &
                            this%dir(this%jac%column_indices(k))
                    end do
                end do
                if (.not.all(ieee_is_finite(this%jdir))) then
                    flag = 1
                    exit
                end if
                rnorm = norm2(fvec + this%jdir) / norm2(fvec)
                if (.not.ieee_is_finite(rnorm) .or. rnorm > rtol) then
                    flag = 1
                    exit
                end if

                ! Apply the line search, if needed
                if (this%get_use_line_search()) then
                    slope = dot_product(fvec, this%jdir)
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
    subroutine sns_alloc_workspaces(this, n, lfil)
        !! Allocates internal workspaces. The ILUT factor arrays are sized
        !! for at most LFIL fill entries in each of the L and U rows, and
        !! they are only grown, never shrunk, so they can be reused.
        class(sparse_newton_solver), intent(inout) :: this
            !! The [[sparse_newton_solver]] object.
        integer, intent(in) :: n
            !! The number of equations.
        integer, intent(in) :: lfil
            !! The maximum number of fill entries per ILUT row.

        ! Local Variables
        integer :: im, iwk

        im = min(n, 50)
        iwk = n * (2 * lfil + 1) + 1

        call resize_real(this%xold, n)
        call resize_real(this%dir, n)
        call resize_real(this%jdir, n)
        call resize_real(this%rhs, n)
        call resize_real(this%w, n + 1)
        call resize_int(this%ilupvt, n)
        call resize_int(this%jw, 2 * n)

        if (.not.allocated(this%vv)) then
            allocate(this%vv(n, im + 1))
        else if (size(this%vv, 1) /= n .or. size(this%vv, 2) /= im + 1) then
            deallocate(this%vv)
            allocate(this%vv(n, im + 1))
        end if

        if (.not.allocated(this%alu)) then
            allocate(this%alu(iwk))
        else if (size(this%alu) < iwk) then
            deallocate(this%alu)
            allocate(this%alu(iwk))
        end if

        if (.not.allocated(this%jlu)) then
            allocate(this%jlu(iwk))
        else if (size(this%jlu) < iwk) then
            deallocate(this%jlu)
            allocate(this%jlu(iwk))
        end if
    end subroutine

! ------------------------------------------------------------------------------
    pure logical function sns_is_valid_csr(a, n) result(valid)
        !! Checks that a square CSR matrix has consistent 1-based index arrays
        !! before they are passed to SPARSKIT.
        type(csr_matrix), intent(in) :: a
            !! The matrix to check.
        integer, intent(in) :: n
            !! The expected number of rows and columns.

        ! Local Variables
        integer :: i, k, nnz

        valid = .false.
        if (.not.allocated(a%row_indices) .or. .not.allocated(a%column_indices)) return
        if (.not.allocated(a%values)) return
        if (size(a%row_indices) /= n + 1) return
        if (a%row_indices(1) /= 1) return
        do i = 1, n
            if (a%row_indices(i+1) < a%row_indices(i)) return
        end do
        nnz = a%row_indices(n+1) - 1
        if (nnz > size(a%column_indices) .or. nnz > size(a%values)) return
        do k = 1, nnz
            if (a%column_indices(k) < 1 .or. a%column_indices(k) > n) return
        end do
        valid = .true.
    end function

! ------------------------------------------------------------------------------
    subroutine resize_real(x, n)
        !! Ensures a real workspace array is allocated with N elements.
        real(real64), allocatable, intent(inout) :: x(:)
        integer, intent(in) :: n

        if (allocated(x)) then
            if (size(x) == n) return
            deallocate(x)
        end if
        allocate(x(n))
    end subroutine

! ------------------------------------------------------------------------------
    subroutine resize_int(x, n)
        !! Ensures an integer workspace array is allocated with N elements.
        integer, allocatable, intent(inout) :: x(:)
        integer, intent(in) :: n

        if (allocated(x)) then
            if (size(x) == n) return
            deallocate(x)
        end if
        allocate(x(n))
    end subroutine

! ------------------------------------------------------------------------------
end module