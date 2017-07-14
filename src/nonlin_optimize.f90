! nonlin_optimize.f90

! REF:
! http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#Nelder-Mead_Simplex
! http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
! https://scicomp.stackexchange.com/questions/14787/fortran-library-for-minimization-or-maximization-of-functions-optimization-prob
! http://ab-initio.mit.edu/wiki/index.php/NLopt

!> @brief \b nonlin_optimize
!!
!! @par Purpose
!! To provide various optimization routines.
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
    !> @brief Defines a solver based upon Nelder and Mead's simplex algorithm
    !! for minimization of functions of multiple variables.
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
        !> @brief Gets an (N+1)-by-N matrix containing the current simplex.
        procedure, public :: get_simplex => nm_get_simplex
        !> @brief Sets an (N+1)-by-N matrix containing the current simplex.
        procedure, public :: set_simplex => nm_set_simplex
        !> @brief Gets the size of the initial simplex.
        procedure, public :: get_initial_size => nm_get_size
        !> @brief Sets the size of the initial simplex.
        procedure, public :: set_initial_size => nm_set_size

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
    !!
    !! @par Usage
    !! The following example illustrates how to find the minimum of Rosenbrock's
    !! function using this Nelder-Mead solver.
    !! @code{.f90}
    !! program example
    !!     use linalg_constants, only : dp, i32
    !!     use nonlin_optimize, only : nelder_mead
    !!     use nonlin_types, only : fcnnvar, fcnnvar_helper, iteration_behavior
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     type(nelder_mead) :: solver
    !!     type(fcnnvar_helper) :: obj
    !!     procedure(fcnnvar), pointer :: fcn
    !!     real(dp) :: x(2), fout
    !!     type(iteration_behavior) :: ib
    !!
    !!     ! Initialization
    !!     fcn => rosenbrock
    !!     call obj%set_fcn(fcn, 2)
    !!
    !!     ! Define an initial guess - the solution is (1, 1)
    !!     call random_number(x)
    !!
    !!     ! Call the solver
    !!     call solver%solve(obj, x, fout, ib)
    !!
    !!     ! Display the output
    !!     print '(AF8.5AF8.5A)', "Rosenbrock Minimum: (", x(1), ", ", x(2), ")"
    !!     print '(AE9.3)', "Function Value: ", fout
    !!     print '(AI0)', "Iterations: ", ib%iter_count
    !!     print '(AI0)', "Function Evaluations: ", ib%fcn_count
    !! contains
    !!     ! Rosenbrock's Function
    !!     function rosenbrock(x) result(f)
    !!         real(dp), intent(in), dimension(:) :: x
    !!         real(dp) :: f
    !!         f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    !!     end function
    !! end
    !! @endcode
    !! The above program yields the following output:
    !! @code{.txt}
    !! Rosenbrock Minimum: ( 1.00000,  1.00000)
    !! Function Value: 0.264E-12
    !! Iterations: 59
    !! Function Evaluations: 105
    !! @endcode
    !!
    !! @par Remarks
    !! The implementation of the Nelder-Mead algorithm presented here is a 
    !! slight modification of the original work of Nelder and Mead.  The 
    !! Numerical Recipes implementation is also quite similar.  In fact, the
    !! Numerical Recipes section relating to reflection, contraction, etc.
    !! is leveraged for this implemetation.
    !!
    !! @par See Also
    !!  - Nelder, John A.; R. Mead (1965). "A simplex method for function 
    !!      minimization". Computer Journal. 7: 308–313.
    !!  - [Gao, Fuchang, Han, Lixing (2010). "Implementing the Nelder-Mead 
    !!      simplex algorithm with adaptive parameters."]
    !!      (http://www.webpages.uidaho.edu/~fuchang/res/ANMS.pdf)
    !!  - [Wikipedia](https://en.wikipedia.org/wiki/Nelder–Mead_method)
    !!  - [Numerical Recipes](http://numerical.recipes/)
    subroutine nm_solve(this, fcn, x, fout, ib, err)
        ! Arguments
        class(nelder_mead), intent(inout) :: this
        class(fcnnvar_helper), intent(in) :: fcn
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), optional :: fout
        type(iteration_behavior), optional :: ib
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0
        real(dp), parameter :: negone = -1.0d0
        real(dp), parameter :: half = 0.5d0
        real(dp), parameter :: two = 2.0d0

        ! Local Variables
        logical :: buildSimplex, fcnvrg
        integer(i32) :: i, ihi, ilo, ihi2, ndim, npts, flag, neval, iter, &
            maxeval
        real(dp) :: ftol, rtol, ftry, fsave, fval, swp
        real(dp), allocatable, dimension(:) :: f, pcent, pmin, work
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
            else
                ! Correct the evaluation count
                neval = neval - 1
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
    !> @brief Extrapolates by the specified factor through the simplex across
    !! from the largest point.  If the extrapolation results in a better
    !! estimate, the current high point is replaced with the new estimate.
    !!
    !! @param[in,out] this The nelder_mead object.
    !! @param[in] fcn The function to evaluate.
    !! @param[in,out] y An array containing the function values at each vertex.
    !! @param[in,out] pcent An array containing the centroid of vertex position
    !!  information.
    !! @param[in] ihi The index of the largest magnitude vertex.
    !! @param[in,out] neval The number of function evaluations.
    !! @param[out] work An N-element workspace array where N is the number of
    !!  dimensions of the problem.
    !! @return The new function estimate.
    function nm_extrapolate(this, fcn, y, pcent, ihi, fac, neval, work) &
            result(ytry)
        ! Arguments
        class(nelder_mead), intent(inout) :: this
        class(fcnnvar_helper), intent(in) :: fcn
        real(dp), intent(inout), dimension(:) :: y, pcent
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
    !> @brief Gets an N-by-(N+1) matrix containing the current simplex.
    !!
    !! @param[in] this The nelder_mead object.
    !! @return The N-by-(N+1) matrix containing the simplex.  Each vertex of the
    !!  simplex is stored as its own column of this matrix.
    pure function nm_get_simplex(this) result(p)
        class(nelder_mead), intent(in) :: this
        real(dp), allocatable, dimension(:,:) :: p
        integer(i32) :: m, n
        if (allocated(this%m_simplex)) then
            m = size(this%m_simplex, 1)
            n = size(this%m_simplex, 2)
            allocate(p(m,n))
            p = this%m_simplex
        end if
    end function

! --------------------
    !> @brief Sets an N-by-(N+1) matrix as the current simplex.  Notice, if this
    !! matrix is different in size from the problem dimensionallity, the 
    !! Nelder-Mead routine will replace it with an appropriately sized matrix.
    !!
    !! @param[in,out] this The nelder_mead object.
    !! @param[in] x The simplex matrix.  Each column of the matrix must contain
    !!  the coordinates of each vertex of the simplex.
    subroutine nm_set_simplex(this, x)
        class(nelder_mead), intent(inout) :: this
        real(dp), dimension(:,:) :: x
        integer(i32) :: m, n
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
    !> @brief Gets the size of the initial simplex that will be utilized by
    !! the Nelder-Mead algorithm in the event that the user does not supply
    !! a simplex geometry, or if the user supplies an invalid simplex geometry.
    !!
    !! @param[in] this The nelder_mead object.
    !! @return The size of the simplex (length of an edge).
    pure function nm_get_size(this) result(x)
        class(nelder_mead), intent(in) :: this
        real(dp) :: x
        x = this%m_initSize
    end function

! --------------------
    !> @brief Sets the size of the initial simplex that will be utilized by
    !! the Nelder-Mead algorithm in the event that the user does not supply
    !! a simplex geometry, or if the user supplies an invalid simplex geometry.
    !!
    !! @param[in,out] this The nelder_mead object.
    !! @param[in] x The size of the simplex (length of an edge).
    subroutine nm_set_size(this, x)
        class(nelder_mead), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_initSize = x
    end subroutine

! ------------------------------------------------------------------------------
end module
