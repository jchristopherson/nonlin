! nonlin_optimize.f90

! REF:
! http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#Nelder-Mead_Simplex
! http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
! https://scicomp.stackexchange.com/questions/14787/fortran-library-for-minimization-or-maximization-of-functions-optimization-prob
! http://ab-initio.mit.edu/wiki/index.php/NLopt

module nonlin_optimize
    use linalg_constants, only : dp, i32
    use ferror, only : errors
    use nonlin_types, only : fcnnvar_helper, optimize_equation, &
        iteration_behavior, NL_OUT_OF_MEMORY_ERROR, NL_CONVERGENCE_ERROR, &
        NL_INVALID_INPUT_ERROR
    implicit none

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !
    type, extends(optimize_equation) :: nelder_mead
        private
        !> The simplex vertices.
        real(dp), allocatable, dimension(:,:) :: m_simplex
        !> A scaling parameter used to define the size of the simplex in each
        !! coordinate direction.
        real(dp) :: m_initSize = 1.0d0
    contains
        !> @brief Optimizes the equation.
        procedure, public :: solve => nm_solve
        !> @brief Extrapolates by the specified factor through the simplex
        !! across from the largest point.  If the extrapolation results in a
        !! better estimate, the current high point is replaced with the new
        !! estimate.
        procedure, private :: extrapolate => nm_extrapolate
    end type

contains
! ******************************************************************************
! NELDER_MEAD MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Utilizes the Nelder-Mead simplex method for finding a minimum
    !! value of the specified function.
    !!
    !! @param[in,out] this The nelder_mead object.
    !! @param[in] fcn The fcnnvar_helper object containing the equation to
    !!  optimize.
    !! @param[in,out] x On input, the initial guess at the optimal point.
    !!  On output, the updated optimal point estimate.
    !! @param[out] fout An optional output, that if provided, returns the
    !!  value of the function at @p x.
    !! @param[out] ib An optional output, that if provided, allows the
    !!  caller to obtain iteration performance statistics.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if @p x is not appropriately sized for
    !!      the problem as defined in @p fcn.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
    !!      available.
    !!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm cannot converge within
    !!      the allowed number of iterations.
    subroutine nm_solve(this, fcn, x, fout, ib, err)
        ! Arguments
        class(nelder_mead), intent(inout) :: this
        class(fcnnvar_helper), intent(in) :: fcn
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), optional :: fout
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: negone = -1.0d0
        real(dp), parameter :: half = 0.5d0
        real(dp), parameter :: two = 2.0d0

        ! Local Variables
        logical :: buildSimplex, fcnvrg
        integer(i32) :: i, ihi, ilo, inhi, ndim, npts, flag, neval, iter, &
            maxeval
        real(dp) :: ftol, rtol, ytry, ysave, swp
        real(dp), allocatable, dimension(:) :: y, psum, pmin, work
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg

        ! Initialization
        ndim = fcn%get_variable_count()
        npts = ndim + 1
        buildSimplex = .true.
        maxeval = this%get_max_fcn_evals()
        ftol = this%get_tolerance()
        iter = 0
        neval = 0
        fcnvrg = .false.
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = .false.
            ib%converge_on_zero_diff = .false.
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(x) /= ndim) then
            write(errmsg, '(AI0AI0A)') &
                "It was expected to receive a coordinate vector of length ", &
                ndim, " , but a vector of length ", size(x), " was received."
            call errmgr%report_error("nm_solve", trim(errmsg), &
                NL_INVALID_INPUT_ERROR)
            return
        end if

        ! Ensure that if an initial simplex was defined, that it is
        ! appropriately sized.  If not, simply create a new simplex of the
        ! appropriate size.
        if (allocated(this%m_simplex)) then
            ! This matrix must be NPTS-BY-NDIM
            if (size(this%m_simplex, 1) /= npts .or. &
                size(this%m_simplex, 2) /= ndim) then
                    deallocate(this%m_simplex)
                buildSimplex = .true.
            else
                ! The simplex is appropriately sized
                buildSimplex = .false.
            end if
        end if

        ! Local Memory Allocation
        allocate(y(npts), stat = flag)
        if (flag == 0) allocate(psum(ndim), stat = flag)
        if (flag == 0) allocate(pmin(ndim), stat = flag)
        if (flag == 0) allocate(work(ndim), stat = flag)
        if (buildSimplex .and. flag == 0) allocate(this%m_simplex(npts, ndim))
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("nm_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Define the initial simplex, if needed
        if (buildSimplex) then
            this%m_simplex(1,:) = x
            do i = 2, npts
                this%m_simplex(i,:) = x + this%m_initSize
            end do
        end if

        ! Evaluate the function at each vertex of the simplex
        do i = 1, npts
            y(i) = fcn%fcn(this%m_simplex(i,:))
        end do
        if (present(fout)) fout = y(1)
        neval = npts

        do i = 1, ndim
            psum(i) = sum(this%m_simplex(:,i))
        end do

        ! Main Loop
        flag = 0 ! Used to check for convergence errors
        do
            ! Update the iteration counter
            iter = iter + 1

            ! Determine the characteristics of each vertex
            ilo = 1
            if (y(1) > y(2)) then
                ihi = 1
                inhi = 2
            else
                ihi = 2
                inhi = 1
            end if
            do i = 1, npts
                if (y(i) <= y(ilo)) ilo = i
                if (y(i) > y(ihi)) then
                    inhi = ihi
                    ihi = i
                else if (y(i) > y(inhi)) then
                    if (i /= ihi) inhi = i
                end if
            end do

            ! Compute the fractional range from highest to lowest, and return
            ! if the result is within tolerance
            rtol = y(ihi) - y(ilo)
            if (rtol < ftol) then
                swp = y(1)
                y(1) = y(ilo)
                y(ilo) = swp
                do i = 1, ndim
                    swp = this%m_simplex(1,i)
                    this%m_simplex(1,i) = this%m_simplex(ilo,i)
                    this%m_simplex(ilo,i) = swp
                    x(i) = this%m_simplex(1,i)
                end do
                if (present(fout)) fout = y(1)
                fcnvrg = .true.
                exit
            end if

            ! Start of a new iteration.  Start by attempting to extrapolate by
            ! a factor of -1 (reflect the simplex at its highest point)
            ytry = this%extrapolate(fcn, y, psum, ihi, negone, neval, work)
            if (ytry <= y(ilo)) then
                ! The result of the extrapolation is better than the current
                ! best point.  As a result, try a factor of 2
                ytry = this%extrapolate(fcn, y, psum, ihi, two, neval, work)
            else if (ytry >= y(inhi)) then
                ! The reflected point is worse than the second highest, so look
                ! for an intermediate lower point (contract the simplex)
                ysave = y(ihi)
                ytry = this%extrapolate(fcn, y, psum, ihi, half, neval, work)
                if (ytry >= ysave) then
                    ! Cannot improve on the high point.  Try to contract around
                    ! the low point.
                    do i = 1, npts
                        if (i /= ilo) then
                            psum = half * (this%m_simplex(i,:) + &
                                this%m_simplex(ilo,:))
                            this%m_simplex(i,:) = psum
                            y(i) = fcn%fcn(psum)
                        end if
                    end do
                    neval = neval + npts
                    do i = 1, ndim
                        psum(i) = sum(this%m_simplex(:,i))
                    end do
                end if
            else
                ! Correct the evaluation count
                neval = neval - 1
            end if

            ! Print iteration status
            if (this%get_print_status()) then
                print *, ""
                print '(AI0)', "Iteration: ", iter
                print '(AI0)', "Function Evaluations: ", neval
                print '(AE8.3)', "Convergence Parameter: ", rtol
            end if

            ! Ensure we haven't made too many function evaluations
            if (neval >= maxeval) then
                flag = 1
                exit
            end if
        end do

        ! Report out iteration statistics
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = .false.
            ib%converge_on_zero_diff = .false.
        end if

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, '(AI0AE8.3AE8.3)') &
                "The algorithm failed to converge." // new_line('c') // &
                "Function evaluations performed: ", neval, new_line('c') // &
                "Convergence Parameter: ", rtol, new_line('c') // &
                "Convergence Criteria: ", ftol
            call errmgr%report_error("nm_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Extrapolates by the specified factor through the simplex across
    !! from the largest point.  If the extrapolation results in a better
    !! estimate, the current high point is replaced with the new estimate.
    !!
    !! @param[in,out] this The nelder_mead object.
    !! @param[in] fcn The function to evaluate.
    !! @param[in,out] y An array containing the function values at each vertex.
    !! @param[in,out] psum An array containing the summation of vertex position
    !!  information.
    !! @param[in] ihi The index of the largest magnitude vertex.
    !! @param[in,out] neval The number of function evaluations.
    !! @param[out] work An N-element workspace array where N is the number of
    !!  dimensions of the problem.
    !! @return The new function estimate.
    function nm_extrapolate(this, fcn, y, psum, ihi, fac, neval, work) &
            result(ytry)
        ! Arguments
        class(nelder_mead), intent(inout) :: this
        class(fcnnvar_helper), intent(in) :: fcn
        real(dp), intent(inout), dimension(:) :: y, psum
        integer(i32), intent(in) :: ihi
        real(dp), intent(in) :: fac
        integer(i32), intent(inout) :: neval
        real(dp), intent(out), dimension(:) :: work
        real(dp) :: ytry

        ! Parameters
        real(dp), parameter :: one = 1.0d0

        ! Local Variables
        integer(i32) :: i, ndim
        real(dp) :: fac1, fac2

        ! Initialization
        ndim = size(this%m_simplex, 2)

         ! Define a trial point
         fac1 = (one - fac) / ndim
         fac2 = fac1 - fac
         do i = 1, ndim
            work(i) = psum(i) * fac1 - this%m_simplex(ihi,i) * fac2
         end do

         ! Evaluate the function at the trial point, and then replace if the
         ! trial provides an improvement
         ytry = fcn%fcn(work)
         neval = neval + 1
         if (ytry < y(ihi)) then
            y(ihi) = ytry
            do i = 1, ndim
                psum(i) = psum(i) + work(i) - this%m_simplex(ihi,i)
                this%m_simplex(ihi,i) = work(i)
            end do
         end if
    end function

! ------------------------------------------------------------------------------
! EXAMPLE NELDER-MEAD (FROSEN)
! http://www.mikehutt.com/neldermead.html

! Michael F. Hutt
! http://www.mikehutt.com
! Mar. 31, 1998
! $Id: frosen.f90,v 1.4 2007/07/10 12:45:32 mike Exp $
!
! This program will attempt to minimize Rosenbrock's function using the 
! Nelder-Mead simplex method. The program was originally developed in C. 
! To be consistent with the way arrays are handled in C, all arrays will 
! start from 0.
!
! to compile this program with g77 use:
! g77 -ffree-form -o frosen.exe frosen.f

! * Copyright (c) 1998-2004 <Michael F. Hutt>
! *
! * Permission is hereby granted, free of charge, to any person obtaining
! * a copy of this software and associated documentation files (the
! * "Software"), to deal in the Software without restriction, including
! * without limitation the rights to use, copy, modify, merge, publish,
! * distribute, sublicense, and/or sell copies of the Software, and to
! * permit persons to whom the Software is furnished to do so, subject to
! * the following conditions:
! *
! * The above copyright notice and this permission notice shall be
! * included in all copies or substantial portions of the Software.
! *
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
! * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
! * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!
! Converted program from Fortran 77 to Fortran 90
! This program will attempt to minimize Rosenbrock's function using the 
! Nelder-Mead simplex method. The program was originally developed in C. 
! To be consistent with the way arrays are handled in C, all arrays will
! start from 0.
! compiles with ELF90
! ======================================================================
! Start of main program
! program frosen
!   implicit none
!   real :: sp(0:1) 
!   integer :: d=2, iprint
!   real :: e=1.0e-4, scale=1.0
!   sp(0) = -1.2
!   sp(1) = 1.0
!   iprint = 0 ! set to 0 to get step by step print out

!   CALL simplex(sp,d,e,scale,iprint)
! stop "OK."

! ! End of main program

! CONTAINS

! ! ======================================================================
! ! This is the function to be minimized

! real function func(x) result(rosen) 
!   implicit none
!   real, intent (in) :: x(:)
!   rosen = (100*((x(2)-x(1)**2)**2)+(1.0-x(1))**2)

!   return 
! end function func

! ! ======================================================================
! ! This is the simplex routine

! subroutine simplex(start, n, EPSILON, scale, iprint)
!   implicit none

!   integer, intent (in) :: n, iprint
!   real, intent (in), dimension(0:n-1) :: start
!   real, intent (in) :: EPSILON, scale

! ! Define Constants
!   integer, parameter :: MAX_IT = 1000
!   real, parameter :: ALPHA=1.0
!   real, parameter :: BETA=0.5
!   real, parameter :: GAMMA=2.0

! ! ======================================================================
! ! Variable Definitions
! ! 
! ! Integer vs = vertex with the smallest value
! ! Integer vh = vertex with next smallest value 
! ! Integer vg = vertex with largest value 
! ! Integer i,j,m,row
! ! Integer k = track the number of function evaluations 
! ! Integer itr = track the number of iterations
! ! real v = holds vertices of simplex 
! ! real pn,qn = values used to create initial simplex 
! ! real f = value of function at each vertex 
! ! real fr = value of function at reflection point 
! ! real fe = value of function at expansion point 
! ! real fc = value of function at contraction point 
! ! real vr = reflection - coordinates 
! ! real ve = expansion - coordinates 
! ! real vc = contraction - coordinates 
! ! real vm = centroid - coordinates 
! ! real min
! ! real fsum,favg,s,cent
! ! real vtmp = temporary array passed to FUNC
! ! ======================================================================

!   Integer :: vs,vh,vg
!   Integer :: i,j,k,itr,m,row
!   real, dimension(:,:), allocatable :: v
!   real, dimension(:), allocatable  :: f
!   real, dimension(:), allocatable :: vr
!   real, dimension(:), allocatable :: ve
!   real, dimension(:), allocatable :: vc
!   real, dimension(:), allocatable :: vm
!   real, dimension(:), allocatable :: vtmp
!   real :: pn,qn
!   real :: fr,fe,fc
!   real :: min,fsum,favg,cent,s

!   allocate (v(0:n,0:n-1))
!   allocate (f(0:n))
!   allocate (vr(0:n-1))
!   allocate (ve(0:n-1))
!   allocate (vc(0:n-1))
!   allocate (vm(0:n-1))
!   allocate (vtmp(0:n-1))

! ! create the initial simplex
! ! assume one of the vertices is 0.0

!   pn = scale*(sqrt(n+1.)-1.+n)/(n*sqrt(2.))
!   qn = scale*(sqrt(n+1.)-1.)/(n*sqrt(2.))

!   DO i=0,n-1
!     v(0,i) = start(i)
!   END DO

!   DO i=1,n
!     DO j=0,n-1
!       IF (i-1 == j) THEN
!         v(i,j) = pn + start(j)
!       ELSE
!         v(i,j) = qn + start(j)
!       END IF
!     END DO
!   END DO


! ! find the initial function values

!   DO j=0,n
! ! put coordinates into single dimension array
! ! to pass it to FUNC
!     DO m=0,n-1
!       vtmp(m) = v(j,m)
!     END DO
!     f(j) = FUNC(vtmp)
!   END DO

! ! Print out the initial simplex
! ! Print out the initial function values

!   IF (iprint == 0) THEN
!     Write(*,*) "Initial Values"
!     Write(*,300) ((v(i,j),j=0,n-1),f(i),i=0,n)
!   END IF

!   k = n+1

! ! begin main loop of the minimization

! DO itr=1,MAX_IT
! ! find the index of the largest value
!   vg = 0
!   DO j=0,n
!     IF (f(j) .GT. f(vg)) THEN
!       vg = j
!     END IF
!   END DO

! ! find the index of the smallest value
!   vs = 0
!   DO j=0,n
!     If (f(j) .LT. f(vs)) Then
!       vs = j
!     END IF
!   END DO

! ! find the index of the second largest value
!   vh = vs
!   Do j=0,n
!     If ((f(j) .GT. f(vh)) .AND. (f(j) .LT. f(vg))) Then
!       vh = j
!     END IF
!   END DO

! ! calculate the centroid
!   DO j=0,n-1
!   cent = 0.0
!     DO m=0,n
!       If (m .NE. vg) Then
!         cent = cent + v(m,j)
!       END IF
!     END DO
!     vm(j) = cent/n
!   END DO

! ! reflect vg to new vertex vr
!   DO j=0,n-1
!     vr(j) = (1+ALPHA)*vm(j) - ALPHA*v(vg,j)
!   END DO
!   fr = FUNC(vr)
!   k = k+1

!   If ((fr .LE. f(vh)) .AND. (fr .GT. f(vs))) Then
!     DO j=0,n-1
!       v(vg,j) = vr(j)
!     END DO
!     f(vg) = fr
!   END IF

! ! investigate a step further in this direction
!   If (fr .LE. f(vs)) Then
!     DO j=0,n-1
!       ve(j) = GAMMA*vr(j) + (1-GAMMA)*vm(j)
!     END DO
!     fe = FUNC(ve)
!     k = k+1

! ! by making fe < fr as opposed to fe < f(vs), Rosenbrocks function
! ! takes 62 iterations as opposed to 64. 

!     If (fe .LT. fr) Then
!       DO j=0,n-1
!         v(vg,j) = ve(j)
!       END DO
!       f(vg) = fe
!     Else
!       DO j=0,n-1
!         v(vg,j) = vr(j)
!       END DO
!       f(vg) = fr
!     END IF
!   END IF

! ! check to see if a contraction is necessary
!   If (fr .GT. f(vh)) Then
!     DO j=0,n-1
!       vc(j) = BETA*v(vg,j) + (1-BETA)*vm(j)
!     END DO
!     fc = FUNC(vc)
!     k = k+1
!     If (fc .LT. f(vg)) Then
!       DO j=0,n-1
!         v(vg,j) = vc(j)
!       END DO
!     f(vg) = fc

! ! at this point the contraction is not successful,
! ! we must halve the distance from vs to all the
! ! vertices of the simplex and then continue.
! ! 10/31/97 - modified C program to account for 
! ! all vertices.

!   Else
!     DO row=0,n
!       If (row .NE. vs) Then
!         DO j=0,n-1
!           v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0
!         END DO
!       END IF
!     END DO
!     DO m=0,n-1
!       vtmp(m) = v(vg,m)
!     END DO
!     f(vg) = FUNC(vtmp)
!     k = k+1

!     DO m=0,n-1
!       vtmp(m) = v(vh,m)
!     END DO
!     f(vh) = FUNC(vtmp)
!     k = k+1
!     END IF
!   END IF

! ! print out the value at each iteration 
!   IF (iprint == 0) THEN
!     Write(*,*) "Iteration ",itr
!     Write(*,300) ((v(i,j),j=0,n-1),f(i),i=0,n)
!   END IF

! ! test for convergence
!   fsum = 0.0
!   DO j=0,n
!     fsum = fsum + f(j)
!   END DO
!   favg = fsum/(n+1.)
!   s = 0.0
!   DO j=0,n
!     s = s + ((f(j)-favg)**2.)/n
!   END DO
!   s = sqrt(s)
!   If (s .LT. EPSILON) Then
!     EXIT ! Nelder Mead has converged - exit main loop
!   END IF
! END DO
! ! end main loop of the minimization
! ! find the index of the smallest value

!   vs = 0
!   DO j=0,n
!     If (f(j) .LT. f(vs)) Then
!       vs = j
!     END IF
!   END DO

! !  print out the minimum

!   DO m=0,n-1
!     vtmp(m) = v(vs,m)
!   END DO

!   min = FUNC(vtmp)
!   k = k+1
!   write(*,*)'The minimum was found at ',v(vs,0),v(vs,1)
!   write(*,250)'The value at the minimum is ',min
!   write(*,*)'The number of function evaluations was',k
!   write(*,*)'The number of iterations was',itr
! 250  FORMAT(A29,F7.4)
! 300  FORMAT(F11.6,F11.6,F11.6)

!   return
!   end subroutine simplex
! ! ======================================================================
! end program frosen


! ------------------------------------------------------------------------------
! Other References:
! - http://people.sc.fsu.edu/~jburkardt/cpp_src/asa047/asa047.html
! - https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
! - http://www.brnt.eu/phd/node10.html#SECTION00622200000000000000
! - http://www.webpages.uidaho.edu/~fuchang/res/ANMS.pdf

! ------------------------------------------------------------------------------
end module
