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
        !! Fixed-parameter continuation is the default. Pseudo-arclength mode
        !! augments the corrector and can follow regular folds in the path.
        !! Extend this type and override homotopy and homotopy_jacobian to
        !! change the path without replacing the continuation algorithm.
        real(real64), private :: m_initialStep = 0.2d0
            !! Initial and maximum continuation step in lambda.
        real(real64), private :: m_minimumStep = 1.0d-5
            !! Minimum continuation step attempted before failure.
        integer, private :: m_maxCorrectorIterations = 12
            !! Maximum Newton corrections attempted at each lambda value.
        logical, private :: m_usePseudoArclength = .false.
            !! Select pseudo-arclength instead of fixed-lambda continuation.
    contains
        procedure, public :: solve => cs_solve
        procedure, private :: solve_fixed_lambda => cs_solve_fixed
        procedure, private :: solve_pseudo_arclength => cs_solve_pseudo
        procedure, public :: homotopy => cs_homotopy
        procedure, public :: homotopy_jacobian => cs_homotopy_jacobian
        procedure, public :: homotopy_lambda_derivative => cs_homotopy_lambda_derivative
        procedure, public :: get_use_pseudo_arclength => cs_get_use_pseudo_arclength
        procedure, public :: set_use_pseudo_arclength => cs_set_use_pseudo_arclength
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
        !! Dispatch to the selected continuation strategy.
        class(continuation_solver), intent(inout) :: this
        class(vecfcn_helper), intent(inout) :: fcn
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: fvec
        type(iteration_behavior), intent(out), optional :: ib
        class(*), intent(in), optional :: args

        if (this%m_usePseudoArclength) then
            call this%solve_pseudo_arclength(fcn, x, fvec, ib, args)
        else
            call this%solve_fixed_lambda(fcn, x, fvec, ib, args)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine cs_solve_fixed(this, fcn, x, fvec, ib, args)
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
        real(real64) :: trial_norm, ftol, nan
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
        nan = ieee_value(nan, ieee_quiet_nan)
        if (present(ib)) then
            ib = iteration_behavior(maxeval)
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
            return
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
                trial_norm = residual_norm
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
                if (present(ib)) then
                    call ib%log_iteration(trial_norm, nan, nan)
                end if
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
    subroutine cs_solve_pseudo(this, fcn, x, fvec, ib, args)
        !! Trace the solution curve using a tangent predictor and an augmented
        !! arclength corrector. Lambda is allowed to reverse direction at folds.
        class(continuation_solver), intent(inout) :: this
        class(vecfcn_helper), intent(inout) :: fcn
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: fvec
        type(iteration_behavior), intent(out), optional :: ib
        class(*), intent(in), optional :: args

        integer :: n, neval, njac, niter, maxeval
        integer :: corrector_iteration, backtrack_iteration
        real(real64) :: ftol, xtol, ds, lambda_tol, lambda, ds_used
        real(real64) :: arc_error, merit, trial_merit, alpha, tangent_norm
        real(real64) :: lambda_new, nan
        logical :: success, corrected, accepted, endpoint_step, endpoint_ok
        real(real64), allocatable :: x0(:), z(:), tangent(:), predictor(:)
        real(real64), allocatable :: candidate(:), residual(:), trial_residual(:)
        real(real64), allocatable :: jac(:,:), lambda_jac(:), augmented(:,:)
        real(real64), allocatable :: rhs(:), delta(:), lu(:,:)
        integer, allocatable :: pivots(:)

        n = fcn%get_variable_count()
        neval = 0
        njac = 0
        niter = 0
        success = .false.
        fvec = ieee_value(0.0d0, ieee_quiet_nan)
        maxeval = this%get_max_fcn_evals()
        ftol = this%get_fcn_tolerance()
        xtol = this%get_var_tolerance()
        nan = ieee_value(nan, ieee_quiet_nan)
        if (present(ib)) then
            ib = iteration_behavior(maxeval)
        end if

        if (.not.fcn%is_fcn_defined()) error stop NL_UNDEFINED_FUNCTION_ERROR
        if (n /= fcn%get_equation_count()) error stop NL_INVALID_INPUT_ERROR
        if (size(x) /= n) error stop 3
        if (size(fvec) /= n) error stop 4
        if (.not.all(ieee_is_finite(x)) .or. maxeval <= 0) return

        allocate(x0(n), z(n+1), tangent(n+1), predictor(n+1), candidate(n+1))
        allocate(residual(n), trial_residual(n), jac(n,n), lambda_jac(n))
        allocate(augmented(n+1,n+1), rhs(n+1), delta(n+1))
        x0 = x
        lambda = 0.0d0
        ds = this%m_initialStep
        lambda_tol = max(ftol, 10.0d0 * epsilon(1.0d0))

        ! The initial tangent is oriented toward increasing lambda. For the
        ! default homotopy, J_H is the identity at lambda=0.
        if (.not.get_derivatives(x, lambda, jac, lambda_jac)) go to 900
        call lu_factor(jac, ipvt = pivots, lu = lu)
        if (.not.all(ieee_is_finite(lu))) go to 900
        if (any(pivots < 1) .or. any(pivots > n)) go to 900
        delta(1:n) = solve_lu(lu, pivots, -lambda_jac)
        delta(n+1) = 1.0d0
        if (.not.all(ieee_is_finite(delta))) go to 900
        tangent_norm = norm2(delta)
        if (tangent_norm <= tiny(tangent_norm)) go to 900
        tangent = delta / tangent_norm

        do while (neval < maxeval)
            z(1:n) = x
            z(n+1) = lambda
            ds_used = ds
            endpoint_step = .false.

            ! Aim the predictor at lambda=1 whenever the tangent points toward
            ! the endpoint. At folds, its lambda component naturally changes sign.
            if (tangent(n+1) > 0.0d0 .and. lambda < 1.0d0) then
                if (lambda + ds_used * tangent(n+1) >= 1.0d0) then
                    ds_used = (1.0d0 - lambda) / tangent(n+1)
                    endpoint_step = .true.
                end if
            end if
            if (ds_used <= 0.0d0) exit

            predictor = z + ds_used * tangent
            candidate = predictor
            corrected = .false.

            do corrector_iteration = 1, this%m_maxCorrectorIterations
                if (.not.evaluate_homotopy(candidate(1:n), candidate(n+1), &
                    residual)) exit
                arc_error = dot_product(tangent, candidate - predictor)
                merit = sqrt(dot_product(residual, residual) + arc_error**2)
                if (maxval(abs(residual)) <= ftol .and. abs(arc_error) <= xtol) then
                    corrected = .true.
                    exit
                end if

                if (.not.get_derivatives(candidate(1:n), candidate(n+1), &
                    jac, lambda_jac)) exit
                augmented(1:n,1:n) = jac
                augmented(1:n,n+1) = lambda_jac
                augmented(n+1,:) = tangent
                rhs(1:n) = -residual
                rhs(n+1) = -arc_error
                call lu_factor(augmented, ipvt = pivots, lu = lu)
                if (.not.all(ieee_is_finite(lu))) exit
                if (any(pivots < 1) .or. any(pivots > n+1)) exit
                delta = solve_lu(lu, pivots, rhs)
                if (.not.all(ieee_is_finite(delta))) exit

                accepted = .false.
                alpha = 1.0d0
                trial_merit = merit
                do backtrack_iteration = 1, 12
                    candidate = candidate + alpha * delta
                    if (.not.all(ieee_is_finite(candidate))) then
                        alpha = 0.5d0 * alpha
                        cycle
                    end if
                    if (.not.evaluate_homotopy(candidate(1:n), candidate(n+1), &
                        trial_residual)) exit
                    arc_error = dot_product(tangent, candidate - predictor)
                    trial_merit = sqrt(dot_product(trial_residual, trial_residual) + &
                        arc_error**2)
                    if (ieee_is_finite(trial_merit) .and. trial_merit < merit) then
                        accepted = .true.
                        exit
                    end if
                    alpha = 0.5d0 * alpha
                end do
                niter = niter + 1
                if (present(ib)) then
                    call ib%log_iteration(trial_merit, nan, nan)
                end if
                if (.not.accepted) exit
            end do

            lambda_new = candidate(n+1)
            if (corrected .and. (endpoint_step .or. lambda_new >= 1.0d0)) then
                ! Interpolate the last arclength segment to lambda=1 and accept
                ! only after checking the actual target-system residual.
                if (lambda_new > lambda) then
                    alpha = min(1.0d0, max(0.0d0, (1.0d0-lambda) / &
                        (lambda_new-lambda)))
                    x = x + alpha * (candidate(1:n)-x)
                    if (neval >= maxeval) exit
                    call fcn%fcn(x, fvec, args)
                    neval = neval + 1
                    success = all(ieee_is_finite(fvec)) .and. maxval(abs(fvec)) <= ftol
                    if (success) exit
                end if
                x = z(1:n)
                lambda = z(n+1)
                ds = 0.5d0 * ds
                if (ds < this%m_minimumStep) exit
                cycle
            end if

            if (corrected) then
                x = candidate(1:n)
                lambda = lambda_new
                if (abs(1.0d0-lambda) <= lambda_tol) then
                    if (neval >= maxeval) exit
                    call fcn%fcn(x, fvec, args)
                    neval = neval + 1
                    success = all(ieee_is_finite(fvec)) .and. maxval(abs(fvec)) <= ftol
                    if (success) exit
                end if

                ! The bordered system selects the tangent consistent with the
                ! previous direction, preventing arbitrary sign reversals.
                if (.not.get_derivatives(x, lambda, jac, lambda_jac)) exit
                augmented = 0.0d0
                augmented(1:n,1:n) = jac
                augmented(1:n,n+1) = lambda_jac
                augmented(n+1,:) = tangent
                rhs = 0.0d0
                rhs(n+1) = 1.0d0
                call lu_factor(augmented, ipvt = pivots, lu = lu)
                if (.not.all(ieee_is_finite(lu))) exit
                if (any(pivots < 1) .or. any(pivots > n+1)) exit
                tangent = solve_lu(lu, pivots, rhs)
                if (.not.all(ieee_is_finite(tangent))) exit
                tangent_norm = norm2(tangent)
                if (tangent_norm <= tiny(tangent_norm)) exit
                tangent = tangent / tangent_norm
                ds = min(this%m_initialStep, 1.25d0 * ds)
            else
                x = z(1:n)
                lambda = z(n+1)
                ds = 0.5d0 * ds
                if (ds < this%m_minimumStep) exit
            end if
        end do

900     continue
        if (present(ib)) then
            ib%iter_count = niter
            ib%fcn_count = neval
            ib%jacobian_count = njac
            ib%gradient_count = 0
            ib%converge_on_fcn = success
            ib%converge_on_chng = .false.
            ib%converge_on_zero_diff = .false.
        end if

    contains
        logical function evaluate_homotopy(point, parameter, value)
            real(real64), intent(in), dimension(:) :: point
            real(real64), intent(in) :: parameter
            real(real64), intent(out), dimension(:) :: value
            if (neval >= maxeval) then
                evaluate_homotopy = .false.
                return
            end if
            call this%homotopy(fcn, point, x0, parameter, value, args)
            neval = neval + 1
            evaluate_homotopy = all(ieee_is_finite(value))
        end function

        logical function get_derivatives(point, parameter, point_jac, lambda_jacobian)
            real(real64), intent(inout), dimension(:) :: point
            real(real64), intent(in) :: parameter
            real(real64), intent(out), dimension(:,:) :: point_jac
            real(real64), intent(out), dimension(:) :: lambda_jacobian
            call this%homotopy_jacobian(fcn, point, x0, parameter, point_jac, args)
            call this%homotopy_lambda_derivative(fcn, point, x0, parameter, &
                lambda_jacobian, args)
            njac = njac + 1
            get_derivatives = all(ieee_is_finite(point_jac)) .and. &
                all(ieee_is_finite(lambda_jacobian))
        end function
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
    subroutine cs_homotopy_lambda_derivative(this, fcn, x, x0, lambda, deriv, args)
        !! Evaluate the default path derivative H_lambda=F(x)-(x-x0).
        !! Override this binding along with homotopy and homotopy_jacobian when
        !! defining a custom homotopy.
        class(continuation_solver), intent(inout) :: this
        class(vecfcn_helper), intent(inout) :: fcn
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(in), dimension(:) :: x0
        real(real64), intent(in) :: lambda
        real(real64), intent(out), dimension(:) :: deriv
        class(*), intent(in), optional :: args
        real(real64), allocatable :: fx(:)

        allocate(fx(size(deriv)))
        call fcn%fcn(x, fx, args)
        deriv = fx - (x - x0)
    end subroutine

! ------------------------------------------------------------------------------
    pure function cs_get_use_pseudo_arclength(this) result(use_pseudo)
        !! Query whether pseudo-arclength continuation is enabled.
        class(continuation_solver), intent(in) :: this
        logical :: use_pseudo
        use_pseudo = this%m_usePseudoArclength
    end function

! ------------------------------------------------------------------------------
    subroutine cs_set_use_pseudo_arclength(this, use_pseudo)
        !! Select pseudo-arclength continuation instead of fixed-lambda steps.
        class(continuation_solver), intent(inout) :: this
        logical, intent(in) :: use_pseudo
        this%m_usePseudoArclength = use_pseudo
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