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
! nonlin_continuation.f90
module nonlin_continuation
    use iso_fortran_env
    use ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use linalg, only : lu_factor, solve_lu
    use nonlin_multi_eqn_mult_var
    use nonlin_error_handling
    use nonlin_types
    implicit none
    private
    public :: continuation_solver

    type, extends(equation_solver) :: continuation_solver
        !! Adaptive parameter continuation for square nonlinear systems.
        !! The default homotopy connects the initial guess x0 to the target
        !! system F(x)=0. Its Jacobian uses vecfcn_helper, which selects an
        !! analytic Jacobian, AD, or finite differences for the target system.
        !! A custom homotopy must satisfy H(x0,0)=0 and H(x,1)=F(x).
        !!
        !! This is fixed-parameter, not pseudo-arclength, continuation. It may
        !! fail when the selected path folds or its Jacobian becomes singular.
        !! Extend this type and override homotopy and homotopy_jacobian to
        !! change the path without replacing the continuation algorithm.
        real(real64), private :: m_initialStep = 0.2d0
            !! Initial and maximum continuation step in lambda.
        real(real64), private :: m_minimumStep = 1.0d-5
            !! Minimum continuation step attempted before failure.
        integer, private :: m_maxCorrectorIterations = 12
            !! Maximum Newton corrections attempted at each lambda value.
    contains
        procedure, public :: solve => cs_solve
        procedure, public :: homotopy => cs_homotopy
        procedure, public :: homotopy_jacobian => cs_homotopy_jacobian
        procedure, public :: get_initial_step => cs_get_initial_step
        procedure, public :: set_initial_step => cs_set_initial_step
        procedure, public :: get_minimum_step => cs_get_minimum_step
        procedure, public :: set_minimum_step => cs_set_minimum_step
        procedure, public :: get_max_corrector_iterations => cs_get_max_corrector_iterations
        procedure, public :: set_max_corrector_iterations => cs_set_max_corrector_iterations
    end type

contains
! ------------------------------------------------------------------------------
    subroutine cs_solve(this, fcn, x, fvec, ib, args)
        !! Follow the homotopy from \(\lambda=0\) to \(\lambda=1\) using damped 
        !! Newton corrections. Failed corrections are retried with smaller 
        !! \(\lambda\) steps; successful steps allow the step size to grow 
        !! again.
        class(continuation_solver), intent(inout) :: this
            !! The [[continuation_solver]] object.
        class(vecfcn_helper), intent(inout) :: fcn
            !! The [[vecfcn_helper]] object containing the equations to solve.
        real(real64), intent(inout), dimension(:) :: x
            !! On input, an N-element array containing an initial estimate to
            !! the solution.  On output, the updated solution estimate.
        real(real64), intent(out), dimension(:) :: fvec
            !! On success, the target residual F(x). On failure, NaNs.
        type(iteration_behavior), intent(out), optional :: ib
            !! An optional output, that if provided, allows the caller to obtain
            !! iteration performance statistics.
        class(*), intent(in), optional :: args
            !! An optional argument to allow the user to communicate with the
            !! routine.

        integer :: n, neval, njac, niter, maxeval
        integer :: corrector_iteration, backtrack_iteration
        real(real64) :: lambda, target_lambda, step, alpha, residual_norm
        real(real64) :: trial_norm, ftol
        logical :: corrected, accepted, reached_end
        real(real64), allocatable :: x0(:), xbase(:), trial_x(:), residual(:)
        real(real64), allocatable :: trial_residual(:), jac(:,:), direction(:)
        real(real64), allocatable :: lu(:,:)
        integer, allocatable :: pivots(:)

        ! Initialize statistics before validating inputs.
        n = fcn%get_variable_count()
        neval = 0
        njac = 0
        niter = 0
        reached_end = .false.
        maxeval = this%get_max_fcn_evals()
        ftol = this%get_fcn_tolerance()
        if (present(ib)) then
            ib%iter_count = 0
            ib%fcn_count = 0
            ib%jacobian_count = 0
            ib%gradient_count = 0
            ib%converge_on_fcn = .false.
            ib%converge_on_chng = .false.
            ib%converge_on_zero_diff = .false.
        end if

        ! This continuation method follows a square system from the supplied
        ! initial point, which also anchors the default homotopy.
        if (.not.fcn%is_fcn_defined()) error stop NL_UNDEFINED_FUNCTION_ERROR
        if (n /= fcn%get_equation_count()) error stop NL_INVALID_INPUT_ERROR
        if (size(x) /= n) error stop 3
        if (size(fvec) /= n) error stop 4

        allocate(x0(n), xbase(n), trial_x(n), residual(n), trial_residual(n), &
            jac(n,n), direction(n))
        x0 = x
        if (.not.all(ieee_is_finite(x0)) .or. maxeval <= 0) then
            fvec = ieee_value(0.0d0, ieee_quiet_nan)
            go to 900
        end if

        lambda = 0.0d0
        step = this%m_initialStep

        ! Correct the solution at each target lambda. On failure, restore the
        ! last accepted point and retry with a smaller continuation step.
        do while (lambda < 1.0d0)
            target_lambda = min(1.0d0, lambda + step)
            xbase = x
            corrected = .false.

            do corrector_iteration = 1, this%m_maxCorrectorIterations
                if (neval >= maxeval) exit
                call this%homotopy(fcn, x, x0, target_lambda, residual, args)
                neval = neval + 1
                if (.not.all(ieee_is_finite(residual))) exit
                residual_norm = maxval(abs(residual))
                if (residual_norm <= ftol) then
                    corrected = .true.
                    exit
                end if

                ! The default binding gets J_F from vecfcn_helper, which uses
                ! AD whenever an AD-compatible function is registered.
                call this%homotopy_jacobian(fcn, x, x0, target_lambda, jac, args)
                njac = njac + 1
                if (.not.all(ieee_is_finite(jac))) exit

                ! Solve J_H*dx=-H. A singular factorization fails this lambda
                ! attempt and triggers a smaller continuation step.
                call lu_factor(jac, ipvt = pivots, lu = lu)
                if (.not.all(ieee_is_finite(lu))) exit
                if (any(pivots < 1) .or. any(pivots > n)) exit
                direction = solve_lu(lu, pivots, -residual)
                if (.not.all(ieee_is_finite(direction))) exit

                ! Backtrack until the corrector reduces the homotopy residual.
                accepted = .false.
                alpha = 1.0d0
                do backtrack_iteration = 1, 12
                    if (neval >= maxeval) exit
                    trial_x = x + alpha * direction
                    if (.not.all(ieee_is_finite(trial_x))) then
                        alpha = 0.5d0 * alpha
                        cycle
                    end if
                    call this%homotopy(fcn, trial_x, x0, target_lambda, &
                        trial_residual, args)
                    neval = neval + 1
                    if (all(ieee_is_finite(trial_residual))) then
                        trial_norm = maxval(abs(trial_residual))
                        if (trial_norm < residual_norm) then
                            x = trial_x
                            residual = trial_residual
                            accepted = .true.
                            exit
                        end if
                    end if
                    alpha = 0.5d0 * alpha
                end do
                niter = niter + 1
                if (.not.accepted) exit
                if (maxval(abs(residual)) <= ftol) then
                    corrected = .true.
                    exit
                end if
            end do

            if (corrected) then
                lambda = target_lambda
                step = min(this%m_initialStep, 1.5d0 * step)
            else
                x = xbase
                if (neval >= maxeval) exit
                step = 0.5d0 * step
                if (step < this%m_minimumStep) exit
            end if
        end do

        if (lambda >= 1.0d0) then
            reached_end = .true.
            fvec = residual
        else
            fvec = ieee_value(0.0d0, ieee_quiet_nan)
        end if

900     continue
        if (present(ib)) then
            ib%iter_count = niter
            ib%fcn_count = neval
            ib%jacobian_count = njac
            ib%gradient_count = 0
            ib%converge_on_fcn = reached_end
            ib%converge_on_chng = .false.
            ib%converge_on_zero_diff = .false.
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine cs_homotopy(this, fcn, x, x0, lambda, h, args)
        !! Evaluate \(H(x,\lambda)=(1-\lambda)*(x-x0)+\lambda*F(x)\).
        !! Extend [[continuation_solver]] and override this binding to change
        !! the path while reusing the continuation and corrector logic. Keep
        !! \(H(x0,0)=0\) and \(H(x,1)=F(x)\), and override homotopy_jacobian 
        !! with its derivative with respect to \(x\).
        class(continuation_solver), intent(inout) :: this
            !! The continuation solver instance.
        class(vecfcn_helper), intent(inout) :: fcn
            !! Helper containing the target system F.
        real(real64), intent(in), dimension(:) :: x
            !! Point at which the homotopy is evaluated.
        real(real64), intent(in), dimension(:) :: x0
            !! Initial point that anchors the default path.
        real(real64), intent(in) :: lambda
            !! Continuation parameter in [0,1].
        real(real64), intent(out), dimension(:) :: h
            !! Homotopy residual.
        class(*), intent(in), optional :: args
            !! User arguments forwarded to F.
        real(real64), allocatable :: f(:)

        allocate(f(size(h)))
        call fcn%fcn(x, f, args)
        h = (1.0d0 - lambda) * (x - x0) + lambda * f
    end subroutine

! ------------------------------------------------------------------------------
    subroutine cs_homotopy_jacobian(this, fcn, x, x0, lambda, jac, args)
        !! Evaluate \(J_H=(1-\lambda)I+\lambda*J_F\). [[vecfcn_helper]] selects 
        !! an analytic Jacobian, AD, or finite differences for \(J_F\). A custom
        !! homotopy override should also override this binding accordingly.
        class(continuation_solver), intent(inout) :: this
            !! The continuation solver instance.
        class(vecfcn_helper), intent(inout) :: fcn
            !! Helper containing the target system and its Jacobian behavior.
        real(real64), intent(inout), dimension(:) :: x
            !! Point at which the Jacobian is evaluated; unchanged on return.
        real(real64), intent(in), dimension(:) :: x0
            !! Initial point that anchors the default path.
        real(real64), intent(in) :: lambda
            !! Continuation parameter in [0,1].
        real(real64), intent(out), dimension(:,:) :: jac
            !! Jacobian of the homotopy with respect to x.
        class(*), intent(in), optional :: args
            !! User arguments forwarded to F.
        integer :: i
        real(real64), allocatable :: jac_f(:,:)

        allocate(jac_f(size(jac,1), size(jac,2)))
        call fcn%jacobian(x, jac_f, args = args)
        jac = lambda * jac_f
        do i = 1, size(jac,1)
            jac(i,i) = jac(i,i) + 1.0d0 - lambda
        end do
    end subroutine

! ------------------------------------------------------------------------------
    pure function cs_get_initial_step(this) result(step)
        !! Get the initial and maximum continuation step.
        class(continuation_solver), intent(in) :: this
            !! Solver whose step size is queried.
        real(real64) :: step
            !! Current initial/maximum lambda increment.
        step = this%m_initialStep
    end function

! ------------------------------------------------------------------------------
    subroutine cs_set_initial_step(this, step)
        !! Set the initial and maximum continuation step, clamped to (0,1].
        class(continuation_solver), intent(inout) :: this
            !! Solver whose step size is changed.
        real(real64), intent(in) :: step
            !! Requested lambda increment. Values outside (0,1] are clamped.
        this%m_initialStep = max(this%m_minimumStep, min(1.0d0, step))
    end subroutine

! ------------------------------------------------------------------------------
    pure function cs_get_minimum_step(this) result(step)
        !! Get the minimum continuation step before failure is reported.
        class(continuation_solver), intent(in) :: this
            !! Solver whose minimum step is queried.
        real(real64) :: step
            !! Current minimum lambda increment.
        step = this%m_minimumStep
    end function

! ------------------------------------------------------------------------------
    subroutine cs_set_minimum_step(this, step)
        !! Set the smallest continuation step that the solver will attempt.
        class(continuation_solver), intent(inout) :: this
            !! Solver whose minimum step is changed.
        real(real64), intent(in) :: step
            !! Requested minimum lambda increment, clamped to [epsilon,1].
        this%m_minimumStep = max(epsilon(1.0d0), min(1.0d0, step))
        this%m_initialStep = max(this%m_initialStep, this%m_minimumStep)
    end subroutine

! ------------------------------------------------------------------------------
    pure function cs_get_max_corrector_iterations(this) result(n)
        !! Get the maximum Newton corrections allowed per continuation step.
        class(continuation_solver), intent(in) :: this
            !! Solver whose corrector limit is queried.
        integer :: n
            !! Current corrector iteration limit.
        n = this%m_maxCorrectorIterations
    end function

! ------------------------------------------------------------------------------
    subroutine cs_set_max_corrector_iterations(this, n)
        !! Set the maximum Newton corrections allowed per continuation step.
        class(continuation_solver), intent(inout) :: this
            !! Solver whose corrector limit is changed.
        integer, intent(in) :: n
            !! Requested limit. Values below one are clamped to one.
        this%m_maxCorrectorIterations = max(1, n)
    end subroutine

! ------------------------------------------------------------------------------
end module