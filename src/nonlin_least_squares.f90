! nonlin_least_squares.f90

!> @brief \b nonlin_least_squares
!!
!! @par Purpose
!! To provide routines capable of solving the nonlinear least squares problem.
module nonlin_least_squares
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use nonlin_constants
    use nonlin_core
    use ferror, only : errors
    implicit none
    private
    public :: least_squares_solver

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief Defines a Levenberg-Marquardt based solver for unconstrained
    !! least-squares problems.
    type, extends(equation_solver) :: least_squares_solver
        private
        !> Initial step bounding factor
        real(real64) :: m_factor = 100.0d0
    contains
        !> @brief Gets a factor used to scale the bounds on the initial step.
        procedure, public :: get_step_scaling_factor => lss_get_factor
        !> @brief Sets a factor used to scale the bounds on the initial step.
        procedure, public :: set_step_scaling_factor => lss_set_factor
        !> @brief Solves the system of equations.
        procedure, public :: solve => lss_solve
    end type

contains
! ******************************************************************************
! LEAST_SQUARES_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Gets a factor used to scale the bounds on the initial step.
    !!
    !! @param[in] this The least_squares_solver object.
    !! @return The factor.
    !!
    !! @par Remarks
    !! This factor is used to set the bounds on the initial step such that the
    !! initial step is bounded as the product of the factor with the Euclidean
    !! norm of the vector resulting from multiplication of the diagonal
    !! scaling matrix and the solution estimate.  If zero, the factor itself
    !! is used.
    pure function lss_get_factor(this) result(x)
        class(least_squares_solver), intent(in) :: this
        real(real64) :: x
        x = this%m_factor
    end function

! --------------------
    !> @brief Sets a factor used to scale the bounds on the initial step.
    !!
    !! @param[in] this The least_squares_solver object.
    !! @param[in] x The factor.  Notice, the factor is limited to the interval
    !!  [0.1, 100].
    !!
    !! @par Remarks
    !! This factor is used to set the bounds on the initial step such that the
    !! initial step is bounded as the product of the factor with the Euclidean
    !! norm of the vector resulting from multiplication of the diagonal
    !! scaling matrix and the solution estimate.  If zero, the factor itself
    !! is used.
    subroutine lss_set_factor(this, x)
        class(least_squares_solver), intent(inout) :: this
        real(real64), intent(in) :: x
        if (x < 0.1d0) then
            this%m_factor = 0.1d0
        else if (x > 1.0d2) then
            this%m_factor = 1.0d2
        else
            this%m_factor = x
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Applies the Levenberg-Marquardt method to solve the nonlinear
    !! least-squares problem.
    !!
    !! @param[in,out] this The least_squares_solver object.
    !! @param[in] fcn The vecfcn_helper object containing the equations to
    !!  solve.
    !! @param[in,out] x On input, an M-element array containing an initial
    !!  estimate to the solution.  On output, the updated solution estimate.
    !!  M is the number of variables.
    !! @param[out] fvec An N-element array that, on output, will contain
    !!  the values of each equation as evaluated at the variable values
    !!  given in @p x.  N is the number of equations.
    !! @param[out] ib An optional output, that if provided, allows the
    !!  caller to obtain iteration performance statistics.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - NL_INVALID_OPERATION_ERROR: Occurs if no equations have been defined.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if the number of equations is less than
    !!      than the number of variables.
    !!  - NL_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not sized
    !!      correctly.
    !!  - NL_CONVERGENCE_ERROR: Occurs if the line search cannot converge within
    !!      the allowed number of iterations.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - NL_TOLERANCE_TOO_SMALL_ERROR: Occurs if the requested tolerance is
    !!      to small to be practical for the problem at hand.
    !!
    !! @par Remarks
    !! This routines is based upon the MINPACK routine LMDIF.
    !!
    !! @par Example 1
    !! The following code provides an example of how to solve a system of N
    !! equations of N unknonwns using the Levenberg-Marquardt method.
    !! @code{.f90}
    !! program main
    !!     use nonlin_core, only : vecfcn, vecfcn_helper
    !!     use nonlin_least_squares, only : least_squares_solver
    !!
    !!     type(vecfcn_helper) :: obj
    !!     procedure(vecfcn), pointer :: fcn
    !!     type(least_squares_solver) :: solver
    !!     real(real64) :: x(2), f(2)
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
    !!         real(real64), intent(in), dimension(:) :: x
    !!         real(real64), intent(out), dimension(:) :: f
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
    !! @par Example 2
    !! @code{.f90}
    !! program example
    !!     use nonlin_core, only : vecfcn_helper, vecfcn
    !!     use nonlin_least_squares, only : least_squares_solver
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     type(vecfcn_helper) :: obj
    !!     procedure(vecfcn), pointer :: fcn
    !!     type(least_squares_solver) :: solver
    !!     real(real64) :: x(4), f(21) ! There are 4 coefficients and 21 data points
    !!
    !!     ! Locate the routine containing the equations to solve
    !!     fcn => fcns
    !!     call obj%set_fcn(fcn, 21, 4)
    !!
    !!     ! Define an initial guess
    !!     x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0, 1.0d0, 1.0d0]
    !!
    !!     ! Solve
    !!     call solver%solve(obj, x, f)
    !!
    !!     ! Display the output
    !!     print "(AF12.8)", "c1: ", x(1)
    !!     print "(AF12.8)", "c2: ", x(2)
    !!     print "(AF12.8)", "c3: ", x(3)
    !!     print "(AF12.8)", "c4: ", x(4)
    !!     print "(AF9.5)", "Max Residual: ", maxval(abs(f))
    !!
    !! contains
    !!     ! The function containing the data to fit
    !!     subroutine fcns(x, f)
    !!         ! Arguments
    !!         real(real64), intent(in), dimension(:) :: x  ! Contains the coefficients
    !!         real(real64), intent(out), dimension(:) :: f
    !!
    !!         ! Local Variables
    !!         real(real64), dimension(21) :: xp, yp
    !!
    !!         ! Data to fit (21 data points)
    !!         xp = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
    !!             0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
    !!             1.8d0, 1.9d0, 2.0d0]
    !!         yp = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
    !!             1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
    !!             1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
    !!             3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
    !!             5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
    !!             10.30252404d0]
    !!
    !!         ! We'll apply a cubic polynomial model to this data:
    !!         ! y = c1 * x**3 + c2 * x**2 + c3 * x + c4
    !!         f = x(1) * xp**3 + x(2) * xp**2 + x(3) * xp + x(4) - yp
    !!
    !!         ! For reference, the data was generated by adding random errors to
    !!         ! the following polynomial: y = x**3 - 0.3 * x**2 + 1.2 * x + 0.3
    !!     end subroutine
    !! end program
    !! @endcode
    !! The above program returns the following results.
    !! @code{.txt}
    !! c1:   1.06476276
    !! c2:  -0.12232029
    !! c3:   0.44661345
    !! c4:   1.18661422
    !! Max Residual:   0.50636
    !! @endcode
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm)
    !! - [MINPACK (Wikipedia)](https://en.wikipedia.org/wiki/MINPACK)
    subroutine lss_solve(this, fcn, x, fvec, ib, err)
        ! Arguments
        class(least_squares_solver), intent(inout) :: this
        class(vecfcn_helper), intent(in) :: fcn
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: fvec
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: p0001 = 1.0d-4
        real(real64), parameter :: p1 = 0.1d0
        real(real64), parameter :: qtr = 0.25d0
        real(real64), parameter :: half = 0.5d0
        real(real64), parameter :: p75 = 0.75d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: hndrd = 1.0d2

        ! Local Variables
        logical :: xcnvrg, fcnvrg, gcnvrg
        integer(int32) :: i, neqn, nvar, flag, neval, iter, j, l, maxeval, &
            njac, lwork
        integer(int32), allocatable, dimension(:) :: jpvt
        real(real64) :: ftol, xtol, gtol, fac, eps, fnorm, par, xnorm, delta, &
            sm, temp, gnorm, pnorm, fnorm1, actred, temp1, temp2, prered, &
            dirder, ratio
        real(real64), allocatable, dimension(:,:) :: jac
        real(real64), allocatable, dimension(:) :: diag, qtf, wa1, wa2, wa3, wa4, w
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg

        ! Initialization
        xcnvrg = .false.
        fcnvrg = .false.
        gcnvrg = .false.
        neqn = fcn%get_equation_count()
        nvar = fcn%get_variable_count()
        neval = 0
        iter = 0
        njac = 0
        fac = this%m_factor
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

        ! Input Check
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("lss_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (nvar > neqn) then
            ! ERROR: System is underdetermined
            call errmgr%report_error("lss_solve", "The solver cannot " // &
                "solve the underdetermined problem.  The number of " // &
                "unknowns must not exceed the number of equations.", &
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
            call errmgr%report_error("lss_solve", trim(errmsg), &
                NL_ARRAY_SIZE_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(jpvt(nvar), stat = flag)
        if (flag == 0) allocate(jac(neqn, nvar), stat = flag)
        if (flag == 0) allocate(diag(nvar), stat = flag)
        if (flag == 0) allocate(qtf(nvar), stat = flag)
        if (flag == 0) allocate(wa1(nvar), stat = flag)
        if (flag == 0) allocate(wa2(nvar), stat = flag)
        if (flag == 0) allocate(wa3(nvar), stat = flag)
        if (flag == 0) allocate(wa4(neqn), stat = flag)
        if (flag == 0) then
            call fcn%jacobian(x, jac, fv = fvec, olwork = lwork)
            allocate(w(lwork), stat = flag)
        end if
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("qns_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Evaluate the function at the starting point, and calculate its norm
        call fcn%fcn(x, fvec)
        neval = 1
        fnorm = norm2(fvec)

        ! Process
        par = zero
        iter = 1
        flag = 0
        do
            ! Evaluate the Jacobian
            call fcn%jacobian(x, jac, fvec, w)
            njac = njac + 1

            ! Compute the QR factorization of the Jacobian
            call lmfactor(jac, .true., jpvt, wa1, wa2, wa3)

            ! On the first iteration, scale the problem according to the norms
            ! of each of the columns of the initial Jacobian
            if (iter == 1) then
                do j = 1, nvar
                    diag(j) = wa2(j)
                    if (wa2(j) == zero) diag(j) = one
                end do
                wa3 = diag * x
                xnorm = norm2(wa3)
                delta = fac * xnorm
                if (delta == zero) delta = fac
            end if

            ! Form Q**T * FVEC, and store the first N components in QTF
            wa4 = fvec
            do j = 1, nvar
                if (jac(j,j) /= zero) then
                    sm = zero
                    do i = j, neqn
                        sm = sm + jac(i,j) * wa4(i)
                    end do
                    temp = -sm / jac(j,j)
                    wa4(j:neqn) = wa4(j:neqn) + jac(j:neqn,j) * temp
                end if ! LINE 120
                jac(j,j) = wa1(j)
                qtf(j) = wa4(j)
            end do

            ! Compute the norm of the scaled gradient
            gnorm = zero
            if (fnorm /= zero) then
                do j = 1, nvar
                    l = jpvt(j)
                    if (wa2(l) == zero) cycle
                    sm = zero
                    do i = 1, j
                        sm = sm + jac(i,j) * (qtf(i) / fnorm)
                    end do
                    gnorm = max(gnorm, abs(sm / wa2(l)))
                end do
            end if ! LINE 170

            ! Test for convergence of the gradient norm
            if (gnorm <= gtol) then
                gcnvrg = .true.
                exit
            end if

            ! Rescale if necessary
            do j = 1, nvar
                diag(j) = max(diag(j), wa2(j))
            end do

            ! Inner Loop
            do
                ! Determine the Levenberg-Marquardt parameter
                call lmpar(jac, jpvt, diag, qtf, delta, par, wa1, wa2, wa3, wa4)

                ! Store the direction P, and X + P.  Calculate the norm of P
                do j = 1, nvar
                    wa1(j) = -wa1(j)
                    wa2(j) = x(j) + wa1(j)
                    wa3(j) = diag(j) * wa1(j)
                end do
                pnorm = norm2(wa3)

                ! On the first iteration, adjust the initial step bounds
                if (iter == 1) delta = min(delta, pnorm)

                ! Evaluate the function at X + P, and calculate its norm
                call fcn%fcn(wa2, wa4)
                neval = neval + 1
                fnorm1 = norm2(wa4)

                ! Compute the scaled actual reduction
                actred = -one
                if (p1 * fnorm1 < fnorm) actred = one - (fnorm1 / fnorm)**2

                ! Compute the scaled predicted reduction and the scaled
                ! directional derivative
                do j = 1, nvar
                    wa3(j) = zero
                    l = jpvt(j)
                    temp = wa1(l)
                    wa3(1:j) = wa3(1:j) + jac(1:j,j) * temp
                end do
                temp1 = norm2(wa3) / fnorm
                temp2 = (sqrt(par) * pnorm) / fnorm
                prered = temp1**2 + temp2**2 / half
                dirder = -(temp1**2 + temp2**2)

                ! Compute the ratio of the actual to the predicted reduction
                ratio = zero
                if (prered /= zero) ratio = actred / prered

                ! Update the step bounds
                if (ratio <= qtr) then
                    if (actred >= zero) temp = half
                    if (actred < zero) temp = half * dirder / &
                        (dirder + half * actred)
                    if (p1 * fnorm1 >= fnorm .or. temp < p1) temp = p1
                    delta = temp * min(delta, pnorm / p1)
                    par = par / temp
                else
                    if (par /= zero .and. ratio < p75) then
                        ! NO ACTION REQUIRED
                    else
                        delta = pnorm / half
                        par = half * par
                    end if
                end if ! LINE 240

                ! Test for successful iteration
                if (ratio >= p0001) then
                    do j = 1, nvar
                        x(j) = wa2(j)
                        wa2(j) = diag(j) * x(j)
                    end do
                    fvec = wa4
                    xnorm = norm2(wa2)
                    fnorm = fnorm1
                    iter = iter + 1
                end if ! LINE 290

                ! Tests for convergence
                if (abs(actred) <= ftol .and. prered <= ftol .and. &
                    half * ratio <= one) fcnvrg = .true.
                if (delta <= xtol * xnorm) xcnvrg = .true.
                if (fcnvrg .or. xcnvrg) exit

                ! Tests for termination and stringent tolerances
                if (neval >= maxeval) flag = NL_CONVERGENCE_ERROR
                if (abs(actred) <= eps .and. prered <= eps .and. &
                    half * ratio <= one) flag = NL_TOLERANCE_TOO_SMALL_ERROR
                if (delta <= eps * xnorm) flag = NL_TOLERANCE_TOO_SMALL_ERROR
                if (gnorm <= eps) flag = NL_TOLERANCE_TOO_SMALL_ERROR
                if (flag /= 0) exit

                if (ratio >= p0001) exit
            end do ! End of the inner loop

            ! Check for termination criteria
            if (fcnvrg .or. xcnvrg .or. gcnvrg .or. flag /= 0) exit

            ! Print the iteration status
            if (this%get_print_status()) then
                call print_status(iter, neval, njac, xnorm, fnorm)
            end if
        end do

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
            call errmgr%report_error("lss_solve", trim(errmsg), &
                flag)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Completes the solution of the Levenberg-Marquardt problem when
    !! provided with a QR factored form of the system Jacobian matrix.  The form
    !! of the problem at this stage is J*X = B (J = Jacobian), and D*X = 0,
    !! where D is a diagonal matrix.
    !!
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R1 of the
    !! QR factorization.  On output, the upper triangular portion is unaltered,
    !! but the strict lower triangle contains the strict upper triangle
    !! (transposed) of the matrix S.
    !! @param[in] ipvt An N-element array tracking the pivoting operations from
    !! the original QR factorization.
    !! @param[in] diag An N-element array containing the diagonal components of
    !! the matrix D.
    !! @param[in] qtb An N-element array containing the first N elements of
    !! Q1**T * B.
    !! @param[in] delta A positive input variable that specifies an upper bounds
    !! on the Euclidean norm of D*X.
    !! @param[in,out] par On input, the initial estimate of the
    !! Levenberg-Marquardt parameter.  On output, the final estimate.
    !! @param[out] x The N-element array that is the solution of A*X = B, and of
    !! D*X = 0.
    !! @param[out] sdiag An N-element array containing the diagonal elements of
    !! the matrix S.
    !! @param[out] wa1 An N-element workspace array.
    !! @param[out] wa2 An N-element workspace array.
    !!
    !! @par Remarks
    !! This routines is based upon the MINPACK routine LMPAR.
    subroutine lmpar(r, ipvt, diag, qtb, delta, par, x, sdiag, wa1, wa2)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: r
        integer(int32), intent(in), dimension(:) :: ipvt
        real(real64), intent(in), dimension(:) :: diag, qtb
        real(real64), intent(in) :: delta
        real(real64), intent(inout) :: par
        real(real64), intent(out), dimension(:) :: x, sdiag, wa1, wa2

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: p001 = 1.0d-3
        real(real64), parameter :: p1 = 0.1d0

        ! Local Variables
        integer :: iter, j, jm1, jp1, k, l, nsing, n
        real(real64) :: dxnorm, dwarf, fp, gnorm, parc, parl, paru, sm, temp

        ! Initialization
        n = size(r, 2)  ! NOTE: R is M-by-N
        dwarf = tiny(dwarf)
        nsing = n

        ! Compute and store in X, the Gauss-Newton direction.  If the Jacobian
        ! is rank deficient, obtain a least-squares solution.
        do j = 1, n
            wa1(j) = qtb(j)
            if (r(j,j) == zero .and. nsing == n) nsing = j - 1
            if (nsing < n) wa1(j) = zero
        end do ! LINE 10

        if (nsing >= 1) then
            do k = 1, nsing
                j = nsing - k + 1
                wa1(j) = wa1(j) / r(j,j)
                temp = wa1(j)
                jm1 = j - 1
                if (jm1 >= 1) then
                    wa1(1:jm1) = wa1(1:jm1) - r(1:jm1,j) * temp
                end if
            end do ! LINE 40
        end if

        ! LINE 50
        do j = 1, n
            l = ipvt(j)
            x(l) = wa1(j)
        end do

        ! Initialize the iteration counter, evaluate the function at the origin,
        ! and test for acceptance of the Gauss-Newton direction.
        iter = 0
        wa2(1:n) = diag * x
        dxnorm = norm2(wa2(1:n))
        fp = dxnorm - delta
        if (fp <= p1 * delta) then
            ! LINE 220
            if (iter == 0) par = zero
            return
        end if

        ! If the Jacobian is not rank deficient, the Newton step provides a
        ! lower bound, PARL, for the zero of the function; else, set this bound
        ! to zero
        parl = zero
        if (nsing == n) then
            do j = 1, n
                l = ipvt(j)
                wa1(j) = diag(l) * (wa2(l) / dxnorm)
            end do ! LINE 80

            do j = 1, n
                sm = zero
                jm1 = j - 1
                if (jm1 >= 1) then
                    sm = dot_product(r(1:jm1,j), wa1(1:jm1))
                end if
                wa1(j) = (wa1(j) - sm) / r(j,j)
            end do ! LINE 110
            temp = norm2(wa1)
            parl = ((fp / delta) / temp) / temp
        end if

        ! Calculate an upper bound, PARU, for the zero of the function
        do j = 1, n
            sm = dot_product(r(1:j,j), qtb(1:j))
            l = ipvt(j)
            wa1(j) = sm / diag(l)
        end do ! LINE 140
        gnorm = norm2(wa1)
        paru = gnorm / delta
        if (paru == zero) paru = dwarf / min(delta, p1)

        ! If the input PAR lies outside of the interval (PARL,PARU), set
        ! PAR to the closer end point.
        par = max(par, parl)
        par = min(par, paru)
        if (par == zero) par = gnorm / dxnorm

        ! Iteration
        do
            iter = iter + 1

            ! Evaluate the function at the current value of PAR
            if (par == zero) par = max(dwarf, p001 * paru)
            temp = sqrt(par)
            wa1 = temp * diag
            call lmsolve(r(1:n,1:n), ipvt, wa1, qtb, x, sdiag, wa2)
            wa2 = diag * x
            dxnorm = norm2(wa2)
            temp = fp
            fp = dxnorm - delta

            ! If the function is small enough, accept the current value of PAR.
            ! Also test for the exceptional cases where PARL is zero, or the
            ! number of iterations has reached 10
            if (abs(fp) <= p1 * delta &
                .or. parl == zero .and. fp <= temp &
                .and. temp < zero .or. iter == 10) exit

            ! Compute the Newton correction
            do j = 1, n
                l = ipvt(j)
                wa1(j) = diag(l) * (wa2(l) / dxnorm)
            end do ! LINE 180
            do j = 1, n
                wa1(j) = wa1(J) / sdiag(j)
                temp = wa1(j)
                jp1 = j + 1
                if (n < jp1) cycle
                wa1 = wa1 - r(1:n,j) * temp
            end do ! LINE 210
            temp = norm2(wa1)
            parc = ((fp / delta) / temp) / temp

            ! Depending on the sign of the function update PARL or PARU
            if (fp > zero) parl = max(parl, par)
            if (fp < zero) paru = min(paru, par)

            ! Compute an improved estimate for PAR
            par = max(parl, par + parc)
        end do ! LINE 220 (End of iteration)
        if (iter == zero) par = zero
        return
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Computes the QR factorization of an M-by-N matrix.
    !!
    !! @param[in,out] a On input, the M-by-N matrix to factor.  On output, the
    !! strict upper triangular portion contains matrix R1 of the factorization,
    !! the lower trapezoidal portion contains the factored form of Q1, and the
    !! diagonal contains the corresponding elementary reflector.
    !! @param[in] pivot Set to true to utilize column pivoting; else, set to
    !! false for no pivoting.
    !! @param[out] ipvt An N-element array that is used to contain the pivot
    !! indices unless @p pivot is set to false.  In such event, this array is
    !! unused.
    !! @param[out] rdiag An N-element array used to store the diagonal elements
    !! of the R1 matrix.
    !! @param[out] acnorm An N-element array used to contain the norms of each
    !! column in the event column pivoting is used.  If pivoting is not used,
    !! this array is unused.
    !! @param[out] wa An N-element workspace array.
    !!
    !! @par Remarks
    !! This routines is based upon the MINPACK routine QRFAC.
    subroutine lmfactor(a, pivot, ipvt, rdiag, acnorm, wa)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: a
        logical, intent(in) :: pivot
        integer(int32), intent(out), dimension(:) :: ipvt
        real(real64), intent(out), dimension(:) :: rdiag, acnorm, wa

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: p05 = 5.0d-2
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: m, n, i, j, jp1, k, kmax, minmn
        real(real64) :: ajnorm, epsmch, sm, temp

        ! Initialization
        m = size(a, 1)
        n = size(a, 2)
        minmn = min(m, n)
        epsmch = epsilon(epsmch)

        ! Compute the initial column norms, and initialize several arrays
        do j = 1, n
            acnorm(j) = norm2(a(:,j))
            rdiag(j) = acnorm(j)
            wa(j) = rdiag(j)
            if (pivot) ipvt(j) = j
        end do ! LINE 10

        ! Reduce A to R with Householder transformations
        do j = 1, minmn
            if (pivot) then
                ! Bring the column of largest norm into the pivot position
                kmax = j
                do k = j, n
                    if (rdiag(k) > rdiag(kmax)) kmax = k
                end do ! LINE 20
                if (kmax /= j) then
                    do i = 1, m
                        temp = a(i,j)
                        a(i,j) = a(i,kmax)
                        a(i,kmax) = temp
                    end do ! LINE 30
                    rdiag(kmax) = rdiag(j)
                    wa(kmax) = wa(j)
                    k = ipvt(j)
                    ipvt(j) = ipvt(kmax)
                    ipvt(kmax) = k
                end if
            end if ! LINE 40

            ! Compute the Householder transformation to reduce the J-th column
            ! of A to a multiple of the J-th unit vector
            ajnorm = norm2(a(j:m,j))
            if (ajnorm /= zero) then
                if (a(j,j) < zero) ajnorm = -ajnorm
                a(j:m,j) = a(j:m,j) / ajnorm
                a(j,j) = a(j,j) + one

                ! Apply the transformation to the remaining columns and update
                ! the norms
                jp1 = j + 1
                if (n >= jp1) then
                    do k = jp1, n
                        sm = dot_product(a(j:m,j), a(j:m,k))
                        temp = sm / a(j,j)
                        a(j:m,k) = a(j:m,k) - temp * a(j:m,j)
                        if (.not.pivot .or. rdiag(k) == zero) cycle
                        temp = a(j,k) / rdiag(k)
                        rdiag(k) = rdiag(k) * sqrt(max(zero, one - temp**2))
                        if (p05 * (rdiag(k) / wa(k))**2 > epsmch) cycle
                        rdiag(k) = norm2(a(jp1:m,k))
                        wa(k) = rdiag(k)
                    end do ! LINE 90
                end if
            end if ! LINE 100
            rdiag(j) = -ajnorm
        end do ! LINE 110
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Solves the QR factored system A*X = B, coupled with the diagonal
    !! system D*X = 0 in the least-squares sense.
    !!
    !! @param[in,out] r On input, the N-by-N upper triangular matrix R1 of the
    !! QR factorization.  On output, the upper triangular portion is unaltered,
    !! but the strict lower triangle contains the strict upper triangle
    !! (transposed) of the matrix S.
    !! @param[in] ipvt An N-element array tracking the pivoting operations from
    !! the original QR factorization.
    !! @param[in] diag An N-element array containing the diagonal components of
    !! the matrix D.
    !! @param[in] qtb An N-element array containing the first N elements of
    !! Q1**T * B.
    !! @param[out] x The N-element array that is the solution of A*X = B, and of
    !! D*X = 0.
    !! @param[out] sdiag An N-element array containing the diagonal elements of
    !! the matrix S.
    !! @param[out] wa An N-element workspace array.
    !!
    !! @par Remarks
    !! This routines is based upon the MINPACK routine QRSOLV.
    subroutine lmsolve(r, ipvt, diag, qtb, x, sdiag, wa)
        ! Arguments
        real(real64), intent(inout), dimension(:,:) :: r
        integer(int32), intent(in), dimension(:) :: ipvt
        real(real64), intent(in), dimension(:) :: diag, qtb
        real(real64), intent(out), dimension(:) :: x, sdiag, wa

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: qtr = 0.25d0
        real(real64), parameter :: half = 0.5d0

        ! Local Variables
        integer(int32) :: n, i, j, jp1, k, kp1, l, nsing
        real(real64) :: cs, ctan, qtbpj, sn, sm, tn, temp

        ! Initialization
        n = size(r, 1)

        ! Copy R and Q**T*B to preserve inputs and initialize S
        do j = 1, n
            r(j:n,j) = r(j,j:n)
            x(j) = r(j,j)
            wa(j) = qtb(j)
        end do ! LINE 20

        ! Eliminate the diagonal matrix D using a Givens rotation
        do j = 1, n
            ! Prepare the row of D to be eliminated, locating the diagonal
            ! element using P from the QR factorization
            l = ipvt(j)
            if (diag(l) /= zero) then
                sdiag(j:n) = zero
                sdiag(j) = diag(l)

                ! The transformations to eliminate the row of D modify only a
                ! single element of Q**T * B beyond the first N, which is
                ! initially zero.
                qtbpj = zero
                do k = j, n
                    ! Determine a Givens rotation which eliminates the
                    ! appropriate element in the current row of D
                    if (sdiag(k) == zero) cycle
                    if (abs(r(k,k)) < abs(sdiag(k))) then
                        ctan = r(k,k) / sdiag(k)
                        sn = half / sqrt(qtr + qtr * ctan**2)
                        cs = sn * ctan
                    else
                        tn = sdiag(k) / r(k,k)
                        cs = half / sqrt(qtr + qtr * tn**2)
                        sn = cs * tn
                    end if

                    ! Compute the modified diagonal element of R and the
                    ! modified element of Q**T * B
                    r(k,k) = cs * r(k,k) + sn * sdiag(k)
                    temp = cs * wa(k) + sn * qtbpj
                    qtbpj = -sn * wa(k) + cs * qtbpj
                    wa(k) = temp

                    ! Accumulate the transformation in the row of S
                    kp1 = k + 1
                    if (n < kp1) cycle
                    do i = kp1, n
                        temp = cs * r(i,k) + sn * sdiag(i)
                        sdiag(i) = -sn * r(i,k) + cs * sdiag(i)
                        r(i,k) = temp
                    end do ! LINE 60
                end do ! LINE 80
            end if

            ! Store the diagonal element of S and restore the corresponding
            ! diagonal element of R
            sdiag(j) = r(j,j)
            r(j,j) = x(j)
        end do ! LINE 100

        ! Solve the triangular system.  If the system is singular, then obtain a
        ! least-squares solution
        nsing = n
        do j = 1, n
            if (sdiag(j) == zero .and. nsing == n) nsing = j - 1
            if (nsing < n) wa(j) = zero
        end do ! LINE 110
        if (nsing >= 1) then
            do k = 1, nsing
                j = nsing - k + 1
                sm = zero
                jp1 = j + 1
                if (nsing >= jp1) then
                    sm = dot_product(r(jp1:nsing,j), wa(jp1:nsing))
                end if
                wa(j) = (wa(j) - sm) / sdiag(j)
            end do ! LINE 140
        end if

        ! Permute the components of Z back to components of X
        do j = 1, n
            l = ipvt(j)
            x(l) = wa(j)
        end do ! LINE 160
    end subroutine

! ------------------------------------------------------------------------------
end module
