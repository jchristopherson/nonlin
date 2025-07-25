module nonlin_least_squares
    use iso_fortran_env
    use nonlin_multi_eqn_mult_var
    use nonlin_error_handling
    use nonlin_types
    use nonlin_helper
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
        !! Defines a Levenberg-Marquardt based solver for unconstrained
        !! least-squares problems.
        real(real64), private :: m_factor = 100.0d0
            !! Initial step bounding factor
    contains
        procedure, public :: get_step_scaling_factor => lss_get_factor
        procedure, public :: set_step_scaling_factor => lss_set_factor
        procedure, public :: solve => lss_solve
    end type

contains
! ******************************************************************************
! LEAST_SQUARES_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    pure function lss_get_factor(this) result(x)
        !! Gets a factor used to scale the bounds on the initial step.
        !!
        !! This factor is used to set the bounds on the initial step such that 
        !! the initial step is bounded as the product of the factor with the 
        !! Euclidean norm of the vector resulting from multiplication of the 
        !! diagonal scaling matrix and the solution estimate.  If zero, the 
        !! factor itself is used.
        class(least_squares_solver), intent(in) :: this
            !! The [[least_squares_solver]] object.
        real(real64) :: x
            !! The factor.
        x = this%m_factor
    end function

! --------------------
    subroutine lss_set_factor(this, x)
        !! Sets a factor used to scale the bounds on the initial step.
        !!
        !! This factor is used to set the bounds on the initial step such that 
        !! the initial step is bounded as the product of the factor with the 
        !! Euclidean norm of the vector resulting from multiplication of the 
        !! diagonal scaling matrix and the solution estimate.  If zero, the 
        !! factor itself is used.
        class(least_squares_solver), intent(inout) :: this
            !! The [[least_squares_solver]] object.
        real(real64), intent(in) :: x
            !! The factor.
        if (x < 0.1d0) then
            this%m_factor = 0.1d0
        else if (x > 1.0d2) then
            this%m_factor = 1.0d2
        else
            this%m_factor = x
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine lss_solve(this, fcn, x, fvec, ib, args, err)
        !! Applies the Levenberg-Marquardt method to solve the nonlinear
        !! least-squares problem.  This routines is based upon the MINPACK 
        !! routine LMDIF.
        class(least_squares_solver), intent(inout) :: this
            !! The [[least_squares_solver]] object.
        class(vecfcn_helper), intent(in) :: fcn
            !! The [[vecfcn_helper]] object containing the equations to solve.
        real(real64), intent(inout), dimension(:) :: x
            !! On input, an N-element array containing an initial estimate 
            !! to the solution.  On output, the updated solution estimate.
            !! N is the number of variables.
        real(real64), intent(out), dimension(:) :: fvec
            !! An M-element array that, on output, will contain the values 
            !! of each equation as evaluated at the variable values given 
            !! in x.
        type(iteration_behavior), optional :: ib
            !! An optional output, that if provided, allows the caller to 
            !! obtain iteration performance statistics.
        class(*), intent(inout), optional :: args
                !! An optional argument to allow the user to communicate with
                !! the routine.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

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
            write(errmsg, 100) "Input number ", flag, &
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
        call fcn%fcn(x, fvec, args)
        neval = 1
        fnorm = norm2(fvec)

        ! Process
        par = zero
        iter = 1
        flag = 0
        do
            ! Evaluate the Jacobian
            call fcn%jacobian(x, jac, fvec, w, args = args)
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
                call fcn%fcn(wa2, wa4, args)
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
            write(errmsg, 101) "The algorithm failed to " // &
                "converge.  Function evaluations performed: ", neval, &
                "." // new_line('c') // "Change in Variable: ", xnorm, &
                new_line('c') // "Residual: ", fnorm
            call errmgr%report_error("lss_solve", trim(errmsg), &
                flag)
        end if

        ! Formatting
100     format(A, I0, A)
101     format(A, I0, A, E10.3, A, E10.3)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine lmpar(r, ipvt, diag, qtb, delta, par, x, sdiag, wa1, wa2)
        !! Completes the solution of the Levenberg-Marquardt problem when
        !! provided with a QR factored form of the system Jacobian matrix.  The
        !! form of the problem at this stage is J*X = B (J = Jacobian), and 
        !! D*X = 0, where D is a diagonal matrix.
        !!
        !! This routines is based upon the MINPACK routine LMPAR.
        real(real64), intent(inout), dimension(:,:) :: r
            !! On input, the N-by-N upper triangular matrix R1 of the
            !! QR factorization.  On output, the upper triangular portion is 
            !! unaltered, but the strict lower triangle contains the strict 
            !! upper triangle (transposed) of the matrix S.
        integer(int32), intent(in), dimension(:) :: ipvt
            !! An N-element array tracking the pivoting operations from the 
            !! original QR factorization.
        real(real64), intent(in), dimension(:) :: diag
            !! An N-element array containing the diagonal components of the 
            !! matrix D.
        real(real64), intent(in), dimension(:) :: qtb
            !! An N-element array containing the first N elements of Q1**T * B.
        real(real64), intent(in) :: delta
            !! A positive input variable that specifies an upper bounds on the 
            !! Euclidean norm of D*X.
        real(real64), intent(inout) :: par
            !! On input, the initial estimate of the Levenberg-Marquardt 
            !! parameter.  On output, the final estimate.
        real(real64), intent(out), dimension(:) :: x
            !! The N-element array that is the solution of A*X = B, and of
            !! D*X = 0.
        real(real64), intent(out), dimension(:) :: sdiag
            !! An N-element array containing the diagonal elements of the 
            !! matrix S.
        real(real64), intent(out), dimension(:) :: wa1
            !! An N-element workspace array.
        real(real64), intent(out), dimension(:) :: wa2
            !! An N-element workspace array.

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
            wa2(1:n) = diag * x
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
    subroutine lmfactor(a, pivot, ipvt, rdiag, acnorm, wa)
        !! Computes the QR factorization of an M-by-N matrix.
        !!
        !! This routines is based upon the MINPACK routine QRFAC.
        real(real64), intent(inout), dimension(:,:) :: a
            !! On input, the M-by-N matrix to factor.  On output, the strict 
            !! upper triangular portion contains matrix R1 of the factorization,
            !! the lower trapezoidal portion contains the factored form of Q1, 
            !! and the diagonal contains the corresponding elementary reflector.
        logical, intent(in) :: pivot
            !! Set to true to utilize column pivoting; else, set to false for 
            !! no pivoting.
        integer(int32), intent(out), dimension(:) :: ipvt
            !! An N-element array that is used to contain the pivot indices 
            !! unless pivot is set to false.  In such event, this array is
            !! unused.
        real(real64), intent(out), dimension(:) :: rdiag
            !! An N-element array used to store the diagonal elements of the 
            !! R1 matrix.
        real(real64), intent(out), dimension(:) :: acnorm
            !! An N-element array used to contain the norms of each column in 
            !! the event column pivoting is used.  If pivoting is not used,
            !! this array is unused.
        real(real64), intent(out), dimension(:) :: wa
            !! An N-element workspace array.

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
    subroutine lmsolve(r, ipvt, diag, qtb, x, sdiag, wa)
        !! Solves the QR factored system A*X = B, coupled with the diagonal 
        !! system D*X = 0 in the least-squares sense.
        !!
        !! This routines is based upon the MINPACK routine QRSOLV.
        real(real64), intent(inout), dimension(:,:) :: r
            !! On input, the N-by-N upper triangular matrix R1 of the QR 
            !! factorization.  On output, the upper triangular portion is 
            !! unaltered, but the strict lower triangle contains the strict 
            !! upper triangle (transposed) of the matrix S.
        integer(int32), intent(in), dimension(:) :: ipvt
            !! An N-element array tracking the pivoting operations from the 
            !! original QR factorization.
        real(real64), intent(in), dimension(:) :: diag
            !! An N-element array containing the diagonal components of the 
            !! matrix D.
        real(real64), intent(in), dimension(:) :: qtb
            !! An N-element array containing the first N elements of Q1**T * B.
        real(real64), intent(out), dimension(:) :: x
            !! The N-element array that is the solution of A*X = B, and of
            !! D*X = 0.
        real(real64), intent(out), dimension(:) :: sdiag
            !! An N-element array containing the diagonal elements of the 
            !! matrix S.
        real(real64), intent(out), dimension(:) :: wa
            !! An N-element workspace array.

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
