! MIT License
!
! Copyright (c) 2017-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
! nonlin_test_continuation.f90
module nonlin_test_continuation
    use iso_fortran_env
    use ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    use nonlin
    use autodiff, only : dual, vector_function, assignment(=), operator(-), &
        operator(*), operator(**)
    implicit none
    private
    public :: test_continuation_default
    public :: test_continuation_custom_homotopy
    public :: test_continuation_evaluation_limit
    public :: test_continuation_fd_fallback
    public :: test_continuation_pseudo_arclength_fold
    public :: test_continuation_rank_deficient_qr

    type, extends(continuation_solver) :: shifted_continuation_solver
        integer :: homotopy_calls = 0
        integer :: jacobian_calls = 0
    contains
        procedure :: homotopy => shifted_homotopy
        procedure :: homotopy_jacobian => shifted_homotopy_jacobian
    end type

    type, extends(continuation_solver) :: fold_continuation_solver
        real(real64) :: minimum_lambda = huge(1.0d0)
        real(real64) :: maximum_lambda = -huge(1.0d0)
    contains
        procedure :: homotopy => fold_homotopy
        procedure :: homotopy_jacobian => fold_homotopy_jacobian
        procedure :: homotopy_lambda_derivative => fold_lambda_derivative
    end type

contains
    subroutine ad_quadratic(x, f, args)
        type(dual), intent(in), dimension(:) :: x
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        real(real64) :: root_squared
        root_squared = 2.0d0
        if (present(args)) then
            select type (args)
            type is (real(real64))
                root_squared = args
            end select
        end if
        f(1) = x(1)**2 - root_squared
    end subroutine

    subroutine real_quadratic(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        real(real64) :: root_squared
        root_squared = 2.0d0
        if (present(args)) then
            select type (args)
            type is (real(real64))
                root_squared = args
            end select
        end if
        f(1) = x(1)**2 - root_squared
    end subroutine

    subroutine rank_deficient_system(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        f(1) = x(1) + x(2) - 1.0d0
        f(2) = 2.0d0 * x(1) + 2.0d0 * x(2) - 2.0d0
    end subroutine

    subroutine fold_target(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        real(real64) :: p
        p = x(1) * (x(1) - 1.0d0) * (x(1) - 2.0d0)
        f(1) = 1.0d0 - p
    end subroutine

    subroutine fold_homotopy(this, fcn, x, x0, lambda, h, args)
        class(fold_continuation_solver), intent(inout) :: this
        class(vecfcn_helper), intent(inout) :: fcn
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(in), dimension(:) :: x0
        real(real64), intent(in) :: lambda
        real(real64), intent(out), dimension(:) :: h
        class(*), intent(in), optional :: args
        real(real64) :: p

        this%minimum_lambda = min(this%minimum_lambda, lambda)
        this%maximum_lambda = max(this%maximum_lambda, lambda)
        p = x(1) * (x(1) - 1.0d0) * (x(1) - 2.0d0)
        h(1) = lambda - p
    end subroutine

    subroutine fold_homotopy_jacobian(this, fcn, x, x0, lambda, jac, args)
        class(fold_continuation_solver), intent(inout) :: this
        class(vecfcn_helper), intent(inout) :: fcn
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(in), dimension(:) :: x0
        real(real64), intent(in) :: lambda
        real(real64), intent(out), dimension(:,:) :: jac
        class(*), intent(in), optional :: args
        jac(1,1) = -(3.0d0 * x(1)**2 - 6.0d0 * x(1) + 2.0d0)
    end subroutine

    subroutine fold_lambda_derivative(this, fcn, x, x0, lambda, deriv, args)
        class(fold_continuation_solver), intent(inout) :: this
        class(vecfcn_helper), intent(inout) :: fcn
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(in), dimension(:) :: x0
        real(real64), intent(in) :: lambda
        real(real64), intent(out), dimension(:) :: deriv
        class(*), intent(in), optional :: args
        deriv = 1.0d0
    end subroutine

    subroutine shifted_homotopy(this, fcn, x, x0, lambda, h, args)
        class(shifted_continuation_solver), intent(inout) :: this
        class(vecfcn_helper), intent(inout) :: fcn
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(in), dimension(:) :: x0
        real(real64), intent(in) :: lambda
        real(real64), intent(out), dimension(:) :: h
        class(*), intent(in), optional :: args
        real(real64), allocatable :: fx(:)

        this%homotopy_calls = this%homotopy_calls + 1
        allocate(fx(size(h)))
        call fcn%fcn(x, fx, args)
        h = (1.0d0 - lambda) * (x - x0) + lambda * fx
        h = h + 0.25d0 * lambda * (1.0d0 - lambda)
    end subroutine

    subroutine shifted_homotopy_jacobian(this, fcn, x, x0, lambda, jac, args)
        class(shifted_continuation_solver), intent(inout) :: this
        class(vecfcn_helper), intent(inout) :: fcn
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(in), dimension(:) :: x0
        real(real64), intent(in) :: lambda
        real(real64), intent(out), dimension(:,:) :: jac
        class(*), intent(in), optional :: args
        real(real64), allocatable :: jac_f(:,:)
        integer :: i

        this%jacobian_calls = this%jacobian_calls + 1
        allocate(jac_f(size(jac,1), size(jac,2)))
        call fcn%jacobian(x, jac_f, args = args)
        jac = lambda * jac_f
        do i = 1, size(jac,1)
            jac(i,i) = jac(i,i) + 1.0d0 - lambda
        end do
    end subroutine

    function test_continuation_default() result(check)
        type(vecfcn_helper) :: equations
        type(continuation_solver) :: solver
        type(iteration_behavior) :: ib
        type(convergence_history) :: history
        procedure(vector_function), pointer :: fcn
        real(real64) :: x(1), residual(1), target
        logical :: check

        fcn => ad_quadratic
        call equations%set_auto_fcn(fcn, 1, 1)
        target = 2.0d0
        x = 1.0d0
        call solver%solve(equations, x, residual, ib, args = target)
        history = ib%get_convergence_history()

        check = abs(x(1) - sqrt(target)) < 1.0d-7
        check = check .and. abs(residual(1)) < 1.0d-8
        check = check .and. ib%converge_on_fcn
        check = check .and. ib%jacobian_count > 0 .and. ib%fcn_count > 0
        check = check .and. size(history%log_residual) == ib%iter_count
        check = check .and. size(history%log_step_change) == ib%iter_count
        check = check .and. size(history%log_gradient) == ib%iter_count
        if (.not.check) print *, "Continuation default path failed:", x, residual
    end function

    function test_continuation_custom_homotopy() result(check)
        type(vecfcn_helper) :: equations
        type(shifted_continuation_solver) :: solver
        procedure(vector_function), pointer :: fcn
        real(real64) :: x(1), residual(1), target
        logical :: check

        fcn => ad_quadratic
        call equations%set_auto_fcn(fcn, 1, 1)
        target = 2.0d0
        x = 1.0d0
        call solver%solve(equations, x, residual, args = target)

        check = abs(x(1) - sqrt(target)) < 1.0d-7
        check = check .and. abs(residual(1)) < 1.0d-8
        check = check .and. solver%homotopy_calls > 0
        check = check .and. solver%jacobian_calls > 0
        if (.not.check) print *, "Continuation custom homotopy failed:", x, residual
    end function

    function test_continuation_fd_fallback() result(check)
        type(vecfcn_helper) :: equations
        type(continuation_solver) :: solver
        procedure(vecfcn), pointer :: fcn
        real(real64) :: x(1), residual(1), target
        logical :: check

        fcn => real_quadratic
        call equations%set_fcn(fcn, 1, 1)
        target = 2.0d0
        x = 1.0d0
        call solver%solve(equations, x, residual, args = target)

        check = abs(x(1) - sqrt(target)) < 1.0d-7
        check = check .and. abs(residual(1)) < 1.0d-8
        if (.not.check) print *, "Continuation finite-difference fallback failed:", x, residual
    end function

    function test_continuation_evaluation_limit() result(check)
        type(vecfcn_helper) :: equations
        type(continuation_solver) :: solver
        type(iteration_behavior) :: ib
        procedure(vector_function), pointer :: fcn
        real(real64) :: x(1), residual(1)
        logical :: check

        fcn => ad_quadratic
        call equations%set_auto_fcn(fcn, 1, 1)
        call solver%set_max_fcn_evals(1)
        x = 1.0d0
        call solver%solve(equations, x, residual, ib)

        check = all(ieee_is_nan(residual)) .and. .not.ib%converge_on_fcn
        check = check .and. ib%fcn_count == 1
        if (.not.check) print *, "Continuation evaluation limit handling failed."
    end function

    function test_continuation_pseudo_arclength_fold() result(check)
        type(vecfcn_helper) :: equations
        type(fold_continuation_solver) :: solver
        type(iteration_behavior) :: ib
        type(convergence_history) :: history
        procedure(vecfcn), pointer :: fcn
        real(real64) :: x(1), residual(1)
        logical :: check

        fcn => fold_target
        call equations%set_fcn(fcn, 1, 1)
        call solver%set_use_pseudo_arclength(.true.)
        call solver%set_initial_step(0.06d0)
        call solver%set_max_fcn_evals(2000)
        call solver%set_max_corrector_iterations(16)
        x = 0.0d0
        call solver%solve(equations, x, residual, ib)
        history = ib%get_convergence_history()

        check = ib%converge_on_fcn .and. abs(residual(1)) < 1.0d-8
        check = check .and. x(1) > 2.0d0 .and. x(1) < 2.5d0
        check = check .and. solver%minimum_lambda < 0.0d0
        check = check .and. solver%maximum_lambda > 0.38d0
        check = check .and. ib%fcn_count > 0 .and. ib%jacobian_count > 0
        check = check .and. size(history%log_residual) == ib%iter_count
        check = check .and. size(history%log_step_change) == ib%iter_count
        check = check .and. size(history%log_gradient) == ib%iter_count
        check = check .and. .not.any(ieee_is_nan(history%log_residual))
        if (.not.check) print *, &
            "Pseudo-arclength did not traverse the fold:", x, residual, &
            solver%minimum_lambda, solver%maximum_lambda, ib%fcn_count, &
            ib%jacobian_count, ib%iter_count
    end function

    function test_continuation_rank_deficient_qr() result(check)
        type(vecfcn_helper) :: equations
        type(continuation_solver) :: solver
        type(iteration_behavior) :: ib
        procedure(vecfcn), pointer :: fcn
        real(real64) :: x(2), residual(2)
        logical :: check

        fcn => rank_deficient_system
        call equations%set_fcn(fcn, 2, 2)
        call solver%set_initial_step(0.25d0)
        x = 0.0d0
        call solver%solve(equations, x, residual, ib)

        check = ib%converge_on_fcn .and. all(ieee_is_finite(x))
        check = check .and. all(ieee_is_finite(residual))
        check = check .and. abs(sum(x) - 1.0d0) < 1.0d-7
        check = check .and. maxval(abs(residual)) < 1.0d-8
        if (.not.check) print *, &
            "Rank-deficient pivoted QR continuation failed:", x, residual
    end function
end module