program nonlin_test_sparse_solver
    use iso_fortran_env, only : real64
    use ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
    use nonlin_sparse_solve, only : sparse_newton_solver
    use nonlin_multi_eqn_mult_var, only : vecfcn_helper, vecfcn
    use nonlin_types, only : iteration_behavior
    use nonlin_linesearch, only : line_search
    use linalg, only : csr_matrix, create_empty_csr_matrix
    implicit none

    type problem_data
        integer :: kind
        real(real64) :: target = 1.0d0
    end type

    type(sparse_newton_solver) :: solver
    type(vecfcn_helper) :: fcn
    type(iteration_behavior) :: ib
    type(line_search) :: ls
    type(problem_data) :: data
    procedure(vecfcn), pointer :: residual
    procedure(sparse_jacobian), pointer :: jac
    character(32) :: test_name
    real(real64), allocatable :: x(:), f(:)
    integer :: n, nfcn, njac

    call get_command_argument(1, test_name)
    n = 2
    data%kind = 1
    data%target = 4.0d0
    nfcn = 0
    njac = 0
    if (test_name == "linear" .or. test_name == "large_sparse") then
        n = 64
        if (test_name == "large_sparse") n = 10000
        data%kind = 2
    else if (test_name == "singular") then
        data%kind = 3
    else if (test_name == "zero_jacobian") then
        data%kind = 4
    else if (test_name == "nonfinite_jacobian") then
        data%kind = 5
    else if (test_name == "nonfinite_residual") then
        data%kind = 6
    else if (test_name == "scaled") then
        data%kind = 7
    else if (test_name == "ill_conditioned") then
        data%kind = 8
    else if (test_name == "large_residual") then
        data%kind = 9
    else if (test_name == "nonfinite_trial") then
        data%kind = 10
    else if (test_name == "unallocated_jacobian") then
        data%kind = 11
    else if (test_name == "wrong_jacobian_shape") then
        data%kind = 12
    else if (test_name == "malformed_jacobian") then
        data%kind = 13
    end if
    allocate(x(n), f(n))
    x = 1.0d0
    residual => residual_fcn
    jac => sparse_jacobian
    call fcn%set_fcn(residual, n, n)
    call fcn%set_sparse_jacobian(jac)

    select case (test_name)
    case ("root")
        x = 2.0d0
        call solver%solve(fcn, x, f, ib, data)
        call require(ib%converge_on_fcn .and. ib%iter_count == 0, "initial root not recognized")
        call require(nfcn == 1 .and. njac == 0, "initial root performed extra evaluations")
        call require(all(f == 0.0d0), "initial root residual incorrect")
    case ("nonlinear", "no_line_search")
        if (test_name == "no_line_search") call solver%set_use_line_search(.false.)
        call solver%solve(fcn, x, f, ib, data)
        call require(ib%converge_on_fcn, "nonlinear solve did not converge on residual")
        call require(maxval(abs(x - 2.0d0)) < 1.0d-7, "incorrect nonlinear solution")
        call require(maxval(abs(f)) < 1.0d-8, "incorrect nonlinear residual")
        call require(ib%fcn_count == nfcn .and. ib%jacobian_count == njac, "incorrect evaluation counts")
    case ("backtrack")
        x = 0.1d0
        call solver%solve(fcn, x, f, ib, data)
        call require(ib%converge_on_fcn, "backtracking solve did not converge")
        call require(maxval(abs(x - 2.0d0)) < 1.0d-7, "backtracking solution incorrect")
        call require(nfcn > ib%iter_count + 1, "expected backtracking did not occur")
        call require(ib%fcn_count == nfcn, "backtracking function count incorrect")
    case ("linear", "large_sparse")
        x = 0.0d0
        call solver%solve(fcn, x, f, ib, data)
        call require(ib%converge_on_fcn, "sparse linear solve did not converge")
        call require(maxval(abs(x - 1.0d0)) < 1.0d-8, "sparse linear solution incorrect")
        call require(maxval(abs(f)) < 1.0d-8, "sparse linear residual incorrect")
        call require(ib%jacobian_count == njac .and. ib%fcn_count == nfcn, "incorrect evaluation counts")
    case ("optional_arguments")
        residual => noargs_residual
        jac => noargs_jacobian
        call fcn%set_fcn(residual, n, n)
        call fcn%set_sparse_jacobian(jac)
        call solver%solve(fcn, x, f)
        call require(maxval(abs(x - 2.0d0)) < 1.0d-7 .and. maxval(abs(f)) < 1.0d-8, &
            "solve without optional arguments failed")
    case ("scaled")
        x = 0.0d0
        call solver%solve(fcn, x, f, ib, data)
        call require(ib%converge_on_fcn, "scaled system did not converge on residual")
        call require(maxval(abs(x - 1.0d0)) < 1.0d-6, "scaled system solution inaccurate")
        call require(maxval(abs(f)) < 1.0d-8, "scaled system true residual too large")
    case ("ill_conditioned")
        x = 0.0d0
        call solver%solve(fcn, x, f, ib, data)
        call require(ib%converge_on_fcn, "ill-conditioned system did not converge on residual")
        call require(maxval(abs(x - 1.0d0)) < 1.0d-5, "ill-conditioned system solution inaccurate")
        call require(maxval(abs(f)) < 1.0d-8, "ill-conditioned system true residual too large")
    case ("large_residual")
        x = 0.0d0
        call solver%solve(fcn, x, f, ib, data)
        if (ib%converge_on_fcn) then
            call require(all(ieee_is_finite(x)) .and. all(ieee_is_finite(f)), &
                "large-residual solve returned nonfinite solution")
            call require(maxval(abs(x - 1.0d0)) < 1.0d-8, "large-residual solution inaccurate")
        else
            call require(all(ieee_is_nan(f)), "large-residual failure not reported")
        end if
    case ("reuse")
        call solver%solve(fcn, x, f, ib, data)
        call require(ib%converge_on_fcn, "first solve failed")
        deallocate(x, f)
        n = 7
        allocate(x(n), f(n))
        x = 0.0d0
        data%kind = 2
        call fcn%set_fcn(residual, n, n)
        call solver%solve(fcn, x, f, ib, data)
        call require(ib%converge_on_fcn .and. maxval(abs(x - 1.0d0)) < 1.0d-8, &
            "reused solver failed after resize")
    case ("zero_jacobian", "nonfinite_jacobian", "nonfinite_residual", "singular", &
        "unallocated_jacobian", "wrong_jacobian_shape", "malformed_jacobian")
        call solver%solve(fcn, x, f, ib, data)
        call require(all(ieee_is_nan(f)), "failed solve did not return NaN residual")
        call require(.not.ib%converge_on_fcn .and. .not.ib%converge_on_chng, &
            "failed solve reported convergence")
    case ("nonfinite_initial")
        x(1) = ieee_value(0.0d0, ieee_quiet_nan)
        call solver%solve(fcn, x, f, ib, data)
        call require(all(ieee_is_nan(f)), "nonfinite initial point did not fail")
        call require(.not.ib%converge_on_fcn .and. njac == 0, "nonfinite initial point continued")
    case ("nonfinite_trial")
        x = 0.1d0
        call solver%solve(fcn, x, f, ib, data)
        call require(all(ieee_is_nan(f)), "nonfinite trial did not produce failure")
        call require(.not.ib%converge_on_fcn .and. .not.ib%converge_on_chng, &
            "nonfinite trial reported convergence")
        call require(ib%fcn_count == nfcn .and. nfcn <= solver%get_max_fcn_evals(), &
            "nonfinite trial exceeded the evaluation budget")
        call require(nfcn == 2, "nonfinite trial was repeatedly evaluated")
    case ("max_evaluations", "backtrack_budget")
        if (test_name == "backtrack_budget") then
            x = 0.1d0
            call solver%set_max_fcn_evals(2)
        else
            call solver%set_max_fcn_evals(1)
        end if
        call solver%solve(fcn, x, f, ib, data)
        call require(nfcn <= solver%get_max_fcn_evals() .and. ib%fcn_count == nfcn, &
            "function evaluation budget exceeded")
        call require(all(ieee_is_nan(f)), "budget exhaustion did not return NaNs")
    case ("missing_fcn")
        block
            type(vecfcn_helper) :: empty
            call solver%solve(empty, x, f)
        end block
    case ("sparse_disabled")
        call fcn%set_use_sparse(.false.)
        call solver%solve(fcn, x, f)
    case ("missing_jacobian")
        nullify(jac)
        call fcn%set_sparse_jacobian(jac)
        call fcn%set_use_sparse(.true.)
        call solver%solve(fcn, x, f)
    case ("nonsquare")
        call fcn%set_fcn(residual, n, n + 1)
        call solver%solve(fcn, x, f)
    case ("empty")
        call fcn%set_fcn(residual, 0, 0)
        call solver%solve(fcn, x(:0), f(:0))
    case ("bad_x")
        call solver%solve(fcn, x(:1), f)
    case ("bad_f")
        call solver%solve(fcn, x, f(:1))
    case ("bad_budget")
        call solver%set_max_fcn_evals(0)
        call solver%solve(fcn, x, f)
    case ("bad_search_budget")
        call ls%set_max_fcn_evals(0)
        call solver%set_line_search(ls)
        call solver%solve(fcn, x, f)
    case default
        print *, "Unknown case: ", trim(test_name)
        stop 1
    end select

    if (test_name == "missing_fcn" .or. test_name == "sparse_disabled" .or. &
        test_name == "missing_jacobian" .or. test_name == "nonsquare" .or. &
        test_name == "empty" .or. &
        test_name == "bad_x" .or. test_name == "bad_f" .or. &
        test_name == "bad_budget" .or. test_name == "bad_search_budget") then
        print *, "Invalid input unexpectedly returned"
        stop 1
    end if

contains
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not.condition) then
            print *, trim(message)
            stop 1
        end if
    end subroutine

    subroutine residual_fcn(x, f, args)
        real(real64), intent(in) :: x(:)
        real(real64), intent(out) :: f(:)
        class(*), intent(in), optional :: args
        integer :: i
        nfcn = nfcn + 1
        select type (args)
        type is (problem_data)
            select case (args%kind)
            case (1, 4, 5, 10, 11, 12, 13)
                f = x**2 - args%target
                if (args%kind == 10 .and. any(x > 1.0d0)) &
                    f = ieee_value(0.0d0, ieee_quiet_nan)
            case (2)
                do i = 1, size(x)
                    f(i) = 2.0d0 * (x(i) - 1.0d0)
                    if (i > 1) f(i) = f(i) - 0.5d0 * (x(i-1) - 1.0d0)
                    if (i < size(x)) f(i) = f(i) - 0.5d0 * (x(i+1) - 1.0d0)
                end do
            case (3)
                f = [x(1) + x(2) - 1.0d0, x(1) + x(2) - 2.0d0]
            case (6)
                f = ieee_value(0.0d0, ieee_quiet_nan)
            case (7)
                f = [1.0d-4 * (x(1) - 1.0d0), 1.0d4 * (x(2) - 1.0d0)]
            case (8)
                f = [x(1) + x(2) - 2.0d0, &
                    x(1) + (1.0d0 + 1.0d-6) * x(2) - (2.0d0 + 1.0d-6)]
            case (9)
                f = 1.0d155 * (x - 1.0d0)
            end select
        class default
            error stop 1
        end select
    end subroutine

    subroutine sparse_jacobian(x, matrix, args)
        real(real64), intent(in) :: x(:)
        type(csr_matrix), intent(out) :: matrix
        class(*), intent(in), optional :: args
        integer :: i, k, n
        njac = njac + 1
        n = size(x)
        select type (args)
        type is (problem_data)
            select case (args%kind)
            case (1, 4, 5, 10)
                matrix = create_empty_csr_matrix(n, n, n)
                do i = 1, n
                    matrix%row_indices(i) = i
                    matrix%column_indices(i) = i
                end do
                matrix%row_indices(n+1) = n+1
                matrix%values = 2.0d0 * x
                if (args%kind == 4) matrix%values = 0.0d0
                if (args%kind == 5) matrix%values = ieee_value(0.0d0, ieee_quiet_nan)
            case (2)
                matrix = create_empty_csr_matrix(n, n, 3*n-2)
                k = 1
                do i = 1, n
                    matrix%row_indices(i) = k
                    if (i > 1) then
                        matrix%column_indices(k) = i-1
                        matrix%values(k) = -0.5d0
                        k = k + 1
                    end if
                    matrix%column_indices(k) = i
                    matrix%values(k) = 2.0d0
                    k = k + 1
                    if (i < n) then
                        matrix%column_indices(k) = i+1
                        matrix%values(k) = -0.5d0
                        k = k + 1
                    end if
                end do
                matrix%row_indices(n+1) = k
            case (3)
                matrix = create_empty_csr_matrix(2, 2, 4)
                matrix%row_indices = [1, 3, 5]
                matrix%column_indices = [1, 2, 1, 2]
                matrix%values = 1.0d0
            case (7)
                matrix = create_empty_csr_matrix(2, 2, 2)
                matrix%row_indices = [1, 2, 3]
                matrix%column_indices = [1, 2]
                matrix%values = [1.0d-4, 1.0d4]
            case (8)
                matrix = create_empty_csr_matrix(2, 2, 4)
                matrix%row_indices = [1, 3, 5]
                matrix%column_indices = [1, 2, 1, 2]
                matrix%values = [1.0d0, 1.0d0, 1.0d0, 1.0d0 + 1.0d-6]
            case (9)
                matrix = create_empty_csr_matrix(n, n, n)
                do i = 1, n
                    matrix%row_indices(i) = i
                    matrix%column_indices(i) = i
                end do
                matrix%row_indices(n+1) = n+1
                matrix%values = 1.0d155
            case (11)
                return
            case (12)
                matrix = create_empty_csr_matrix(2, 3, 2)
                matrix%row_indices = [1, 2, 3]
                matrix%column_indices = [1, 2]
                matrix%values = 2.0d0
            case (13)
                matrix = create_empty_csr_matrix(2, 2, 2)
                matrix%row_indices = [1, 4, 3]
                matrix%column_indices = [1, 2]
                matrix%values = 2.0d0
            end select
        class default
            error stop 1
        end select
    end subroutine

    subroutine noargs_residual(x, f, args)
        real(real64), intent(in) :: x(:)
        real(real64), intent(out) :: f(:)
        class(*), intent(in), optional :: args
        f = x**2 - 4.0d0
    end subroutine

    subroutine noargs_jacobian(x, matrix, args)
        real(real64), intent(in) :: x(:)
        type(csr_matrix), intent(out) :: matrix
        class(*), intent(in), optional :: args
        integer :: i
        matrix = create_empty_csr_matrix(size(x), size(x), size(x))
        do i = 1, size(x)
            matrix%row_indices(i) = i
            matrix%column_indices(i) = i
        end do
        matrix%row_indices(size(x)+1) = size(x)+1
        matrix%values = 2.0d0 * x
    end subroutine
end program
