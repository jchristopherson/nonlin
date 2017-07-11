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
        iteration_behavior, NL_OUT_OF_MEMORY_ERROR
    implicit none

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !
    type, extends(optimize_equation) :: nelder_mead
        private
        !> The simplex vertices.
        real(dp), allocatable, dimension(:,:) :: m_simplex
    contains
        !> @brief Optimizes the equation.
        procedure, public :: optimize => nm_optimize
    end type

contains
! ******************************************************************************
! NELDER_MEAD MEMBERS
! ------------------------------------------------------------------------------
    !
    subroutine nm_optimize(this, fcn, x, fout, ib, err)
        ! Arguments
        class(nelder_mead), intent(inout) :: this
        class(fcnnvar_helper), intent(in) :: fcn
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), optional :: fout
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        logical :: buildSimplex, fcnvrg
        integer(i32) :: i, ihi, ilo, inhi, ndim, npts, flag, neval, iter, &
            maxeval
        real(dp) :: ftol
        real(dp), allocatable, dimension(:) :: y, psum, pmin
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg

        ! Initialization
        ndim = size(x)
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
        if (buildSimplex .and. flag == 0) allocate(this%m_simplex(npts, ndim))
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("nm_optimize", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Define the initial simplex, if needed
        if (buildSimplex) then
            ! TO DO: Build the initial estimate of the simplex geometry
        end if

        ! Evaluate the function at each vertex of the simplex
        do i = 1, npts
            y(i) = fcn%fcn(this%m_simplex(i,:))
        end do
        neval = npts
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
