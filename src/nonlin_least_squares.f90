! nonlin_least_squares.f90

!> @brief \b nonlin_least_squares
!!
!! @par Purpose
!! To provide routines capable of solving the nonlinear least squares problem.
module nonlin_least_squares
    use linalg_constants, only : dp, i32
    use nonlin_types
    use ferror, only : errors
    implicit none
    private
    public :: least_squares_solver

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !
    type, extends(equation_solver) :: least_squares_solver
    contains
        !> @brief Solves the system of equations.
        procedure, public :: solve => lss_solve
    end type

contains
! ******************************************************************************
! LEAST_SQUARES_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    !
    subroutine lss_solve(this, fcn, x, fvec, ib, err)
        ! Arguments
        class(least_squares_solver), intent(inout) :: this
        class(vecfcn_helper), intent(in) :: fcn
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: fvec
        type(iteration_behavior), optional :: ib
        class(errors), intent(in), optional, target :: err

        ! Local Variables
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
        real(dp), intent(inout), dimension(:,:) :: r
        integer(i32), intent(in), dimension(:) :: ipvt
        real(dp), intent(in), dimension(:) :: diag, qtb
        real(dp), intent(in) :: delta
        real(dp), intent(inout) :: par
        real(dp), intent(out), dimension(:) :: x, sdiag, wa1, wa2

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: p001 = 1.0d-3
        real(dp), parameter :: p1 = 0.1d0

        ! Local Variables
        integer :: iter, j, jm1, jp1, k, l, nsing, n
        real(dp) :: dxnorm, dwarf, fp, gnorm, parc, parl, paru, sm, temp

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
        wa2 = diag * x
        dxnorm = norm2(wa2)
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
            call lmsolve(r, ipvt, wa1, qtb, x, sdiag, wa2)
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
        real(dp), intent(inout), dimension(:,:) :: a
        logical, intent(in) :: pivot
        integer(i32), intent(out), dimension(:) :: ipvt
        real(dp), intent(out), dimension(:) :: rdiag, acnorm, wa

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: p05 = 5.0d-2
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: m, n, i, j, jp1, k, kmax, minmn
        real(dp) :: ajnorm, epsmch, sm, temp

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
        real(dp), intent(inout), dimension(:,:) :: r
        integer(i32), intent(in), dimension(:) :: ipvt
        real(dp), intent(in), dimension(:) :: diag, qtb
        real(dp), intent(out), dimension(:) :: x, sdiag, wa

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: qtr = 0.25d0
        real(dp), parameter :: half = 0.5d0

        ! Local Variables
        integer(i32) :: n, i, j, jp1, k, kp1, l, nsing
        real(dp) :: cs, ctan, qtbpj, sn, sm, tn, temp

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
