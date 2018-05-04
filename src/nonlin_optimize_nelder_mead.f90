! nonlin_optimize_nelder_mead.f90

submodule (nonlin_optimize) nonlin_optimize_nelder_mead
contains
! ------------------------------------------------------------------------------
    module subroutine nm_solve(this, fcn, x, fout, ib, err)
        ! Arguments
        class(nelder_mead), intent(inout) :: this
        class(fcnnvar_helper), intent(in) :: fcn
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(out), optional :: fout
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: negone = -1.0d0
        real(real64), parameter :: half = 0.5d0
        real(real64), parameter :: two = 2.0d0

        ! Local Variables
        logical :: buildSimplex, fcnvrg
        integer(int32) :: i, ihi, ilo, ihi2, ndim, npts, flag, neval, iter, &
            maxeval
        real(real64) :: ftol, rtol, ftry, fsave, fval, swp
        real(real64), allocatable, dimension(:) :: f, pcent, pmin, work
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
            ib%gradient_count = 0
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
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("nm_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
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
            ! This matrix must be NDIM-by-NPTS
            if (size(this%m_simplex, 1) /= ndim .or. &
                size(this%m_simplex, 2) /= npts) then
                    deallocate(this%m_simplex)
                buildSimplex = .true.
            else
                ! The simplex is appropriately sized
                buildSimplex = .false.
            end if
        end if

        ! Local Memory Allocation
        allocate(f(npts), stat = flag)
        if (flag == 0) allocate(pcent(ndim), stat = flag)
        if (flag == 0) allocate(pmin(ndim), stat = flag)
        if (flag == 0) allocate(work(ndim), stat = flag)
        if (buildSimplex .and. flag == 0) allocate(this%m_simplex(ndim, npts))
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("nm_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Define the initial simplex, if needed
        if (buildSimplex) then
            this%m_simplex(:,1) = x
            do i = 2, npts
                this%m_simplex(:,i) = x
            end do
            do i = 1, ndim
                this%m_simplex(i,i+1) = this%m_simplex(i,i+1) + this%m_initSize
            end do
        end if

        ! Evaluate the function at each vertex of the simplex
        do i = 1, npts
            f(i) = fcn%fcn(this%m_simplex(:,i))
        end do
        neval = npts
        fval = f(1)

        do i = 1, ndim
            pcent(i) = sum(this%m_simplex(i,:))
        end do

        ! Main Loop
        flag = 0 ! Used to check for convergence errors
        do
            ! Update the iteration counter
            iter = iter + 1

            ! Determine the characteristics of each vertex
            ilo = 1
            if (f(1) > f(2)) then
                ihi = 1
                ihi2 = 2
            else
                ihi = 2
                ihi2 = 1
            end if
            do i = 1, npts
                if (f(i) <= f(ilo)) ilo = i
                if (f(i) > f(ihi)) then
                    ihi2 = ihi
                    ihi = i
                else if (f(i) > f(ihi2)) then
                    if (i /= ihi) ihi2 = i
                end if
            end do

            ! Check for convergence.  Nelder and Mead recommend using the
            ! following convergence test: sqrt(sum(f - favg)**2 / n); however,
            ! it seems that a sufficient check may be made using only the
            ! extreme function values of the simplex (highest and lowest valued
            ! points).
            rtol = abs(f(ihi) - f(ilo))
            if (rtol < ftol) then
                swp = f(1)
                f(1) = f(ilo)
                f(ilo) = swp
                do i = 1, ndim
                    swp = this%m_simplex(i,1)
                    this%m_simplex(i,1) = this%m_simplex(i,ilo)
                    this%m_simplex(i,ilo) = swp
                    x(i) = this%m_simplex(i,1)
                end do
                fval = f(1)
                fcnvrg = .true.
                exit
            end if

            ! Start of a new iteration by reflecting the simplex at its largest
            ! point.
            ftry = this%extrapolate(fcn, f, pcent, ihi, negone, neval, work)
            if (ftry <= f(ilo)) then
                ! The result of the reflection is better than the current
                ! best point.  As a result, try a factor of 2 in the reflected
                ! direction.  Again, the highest point is of interest.
                ftry = this%extrapolate(fcn, f, pcent, ihi, two, neval, work)
            else if (ftry >= f(ihi2)) then
                ! The reflected point is worse than the second highest, so look
                ! for an intermediate lower point (contract the simplex)
                fsave = f(ihi)
                ftry = this%extrapolate(fcn, f, pcent, ihi, half, neval, work)
                if (ftry >= fsave) then
                    ! Cannot improve on the high point.  Try to contract around
                    ! the low point.
                    do i = 1, npts
                        if (i /= ilo) then
                            pcent = half * (this%m_simplex(:,i) + &
                                this%m_simplex(:,ilo))
                            this%m_simplex(:,i) = pcent
                            f(i) = fcn%fcn(pcent)
                        end if
                    end do
                    neval = neval + npts
                    do i = 1, ndim
                        pcent(i) = sum(this%m_simplex(i,:))
                    end do
                end if
            end if

            ! Print iteration status
            if (this%get_print_status()) then
                print *, ""
                print '(AI0)', "Iteration: ", iter
                print '(AI0)', "Function Evaluations: ", neval
                print '(AE8.3)', "Function Value: ", fval
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
            ib%gradient_count = 0
            ib%converge_on_fcn = fcnvrg
            ib%converge_on_chng = .false.
            ib%converge_on_zero_diff = .false.
        end if

        ! Get the function value at the computed minimum
        if (present(fout)) fout = fval

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
    module function nm_extrapolate(this, fcn, y, pcent, ihi, fac, neval, &
            work) result(ytry)
        ! Arguments
        class(nelder_mead), intent(inout) :: this
        class(fcnnvar_helper), intent(in) :: fcn
        real(real64), intent(inout), dimension(:) :: y, pcent
        integer(int32), intent(in) :: ihi
        real(real64), intent(in) :: fac
        integer(int32), intent(inout) :: neval
        real(real64), intent(out), dimension(:) :: work
        real(real64) :: ytry

        ! Parameters
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: i, ndim
        real(real64) :: fac1, fac2

        ! Initialization
        ndim = size(this%m_simplex, 1)

         ! Define a trial point
         fac1 = (one - fac) / ndim
         fac2 = fac1 - fac
         do i = 1, ndim
            work(i) = pcent(i) * fac1 - this%m_simplex(i,ihi) * fac2
         end do

         ! Evaluate the function at the trial point, and then replace if the
         ! trial provides an improvement
         ytry = fcn%fcn(work)
         neval = neval + 1
         if (ytry < y(ihi)) then
            y(ihi) = ytry
            do i = 1, ndim
                pcent(i) = pcent(i) + work(i) - this%m_simplex(i,ihi)
                this%m_simplex(i,ihi) = work(i)
            end do
         end if
    end function

! ------------------------------------------------------------------------------
    pure module function nm_get_simplex(this) result(p)
        class(nelder_mead), intent(in) :: this
        real(real64), allocatable, dimension(:,:) :: p
        integer(int32) :: m, n
        if (allocated(this%m_simplex)) then
            m = size(this%m_simplex, 1)
            n = size(this%m_simplex, 2)
            allocate(p(m,n))
            p = this%m_simplex
        end if
    end function

! --------------------
    module subroutine nm_set_simplex(this, x)
        class(nelder_mead), intent(inout) :: this
        real(real64), dimension(:,:) :: x
        integer(int32) :: m, n
        m = size(x, 1)
        n = size(x, 2)
        if (allocated(this%m_simplex)) then
            if (size(this%m_simplex, 1) /= m .or. &
                size(this%m_simplex, 2) /= n) then
                deallocate(this%m_simplex)
                allocate(this%m_simplex(m, n))
            end if
        else
            allocate(this%m_simplex(m, n))
        end if
        this%m_simplex = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function nm_get_size(this) result(x)
        class(nelder_mead), intent(in) :: this
        real(real64) :: x
        x = this%m_initSize
    end function

! --------------------
    module subroutine nm_set_size(this, x)
        class(nelder_mead), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_initSize = x
    end subroutine
end submodule
