! nonlin_optimize.f90

! REF:
! http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#Nelder-Mead_Simplex
! http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
! https://scicomp.stackexchange.com/questions/14787/fortran-library-for-minimization-or-maximization-of-functions-optimization-prob
! http://ab-initio.mit.edu/wiki/index.php/NLopt

module nonlin_optimize
    use iso_fortran_env
    use ferror
    use nonlin_linesearch
    use nonlin_error_handling
    use nonlin_multi_var
    use nonlin_types
    use lapack
    use linalg, only : rank1_update, tri_mtx_mult, cholesky_rank1_update, &
        cholesky_rank1_downdate, solve_cholesky
    use linalg_errors, only : LA_MATRIX_FORMAT_ERROR
    implicit none
    private
    public :: nelder_mead
    public :: line_search_optimizer
    public :: bfgs
    
    type, extends(equation_optimizer) :: nelder_mead
        !! Defines a solver based upon Nelder and Mead's simplex algorithm
        !! for minimization of functions of multiple variables.
        real(real64), private, allocatable, dimension(:,:) :: m_simplex
            !! The simplex vertices.
        real(real64), private :: m_initSize = 1.0d0
            !! A scaling parameter used to define the size of the simplex in 
            !! each coordinate direction.
    contains
        procedure, public :: solve => nm_solve
        procedure, public :: get_simplex => nm_get_simplex
        procedure, public :: set_simplex => nm_set_simplex
        procedure, public :: get_initial_size => nm_get_size
        procedure, public :: set_initial_size => nm_set_size
        procedure, private :: extrapolate => nm_extrapolate
    end type

    type, abstract, extends(equation_optimizer) :: line_search_optimizer
        !! A class describing equation optimizers that use a line search
        !! algorithm to improve convergence behavior.
        class(line_search), private, allocatable :: m_lineSearch
            !! The line search object.
        logical, private :: m_useLineSearch = .true.
            !! Set to true if a line search should be used regardless of the 
            !! status of m_lineSearch
        real(real64), private :: m_xtol = 1.0d-12
            !! The convergence criteria on change in variable.
    contains
        procedure, public :: get_line_search => lso_get_line_search
        procedure, public :: set_line_search => lso_set_line_search
        procedure, public :: set_default_line_search => lso_set_default
        procedure, public :: is_line_search_defined => &
            lso_is_line_search_defined
        procedure, public :: get_use_line_search => lso_get_use_search
        procedure, public :: set_use_line_search => lso_set_use_search
        procedure, public :: get_var_tolerance => lso_get_var_tol
        procedure, public :: set_var_tolerance => lso_set_var_tol
    end type

    type, extends(line_search_optimizer) :: bfgs
        !! Defines a Broyden–Fletcher–Goldfarb–Shanno (BFGS) solver for
        !! minimization of functions of multiple variables.
        !!
        !! See Also:
        !!
        !! - [Wikipedia - BFGS Methods](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm)
        !!
        !! - [Wikipedia - Quasi-Newton Methods](https://en.wikipedia.org/wiki/Quasi-Newton_method)
        !!
        !! - [minFunc](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html)
    contains
        procedure, public :: solve => bfgs_solve
    end type

    interface
        subroutine DSYMV(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
            use iso_fortran_env, only : int32, real64
            character, intent(in) :: uplo
            integer(int32), intent(in) :: n, lda, incx, incy
            real(real64), intent(in) :: alpha, beta, a(lda,*), x(*)
            real(real64), intent(inout) :: y(*)
        end subroutine
    end interface

contains
! ******************************************************************************
! NELDER_MEAD
! ------------------------------------------------------------------------------
    subroutine nm_solve(this, fcn, x, fout, ib, args, err)
        !! Utilizes the Nelder-Mead simplex method for finding a minimum
        !! value of the specified function.
        !!
        !! The implementation of the Nelder-Mead algorithm presented here is a
        !! slight modification of the original work of Nelder and Mead.  The
        !! Numerical Recipes implementation is also quite similar.  In fact, the
        !! Numerical Recipes section relating to reflection, contraction, etc.
        !! is leveraged for this implemetation.
        !!
        !! See Also:
        !!
        !!  - Nelder, John A.; R. Mead (1965). "A simplex method for function
        !!      minimization". Computer Journal. 7: 308–313.
        !!
        !!  - [Gao, Fuchang, Han, Lixing (2010). "Implementing the Nelder-Mead
        !!      simplex algorithm with adaptive parameters."]
        !!      (http://www.webpages.uidaho.edu/~fuchang/res/ANMS.pdf)
        !!
        !!  - [Wikipedia](https://en.wikipedia.org/wiki/Nelder–Mead_method)
        !!
        !!  - [Numerical Recipes](http://numerical.recipes/)
        class(nelder_mead), intent(inout) :: this
            !! The [[nelder_mead]] object.
        class(fcnnvar_helper), intent(in) :: fcn
            !! The [[fcnnvar_helper]] object containing the equation to 
            !! optimize.
        real(real64), intent(inout), dimension(:) :: x
            !! On input, the initial guess at the optimal point.  On output, 
            !! the updated optimal point estimate.
        real(real64), intent(out), optional :: fout
            !! An optional output, that if provided, returns the value of the 
            !! function at x.
        type(iteration_behavior), optional :: ib
            !! An optional output, that if provided, allows the caller to 
            !! obtain iteration performance statistics.
        class(*), intent(inout), optional :: args
            !! An optional argument to allow the user to communicate with fcn.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

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
            write(errmsg, 100) &
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
            f(i) = fcn%fcn(this%m_simplex(:,i), args)
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
            ftry = this%extrapolate(fcn, f, pcent, ihi, negone, neval, work, &
                args = args)
            if (ftry <= f(ilo)) then
                ! The result of the reflection is better than the current
                ! best point.  As a result, try a factor of 2 in the reflected
                ! direction.  Again, the highest point is of interest.
                ftry = this%extrapolate(fcn, f, pcent, ihi, two, neval, work, &
                    args = args)
            else if (ftry >= f(ihi2)) then
                ! The reflected point is worse than the second highest, so look
                ! for an intermediate lower point (contract the simplex)
                fsave = f(ihi)
                ftry = this%extrapolate(fcn, f, pcent, ihi, half, neval, work, &
                    args = args)
                if (ftry >= fsave) then
                    ! Cannot improve on the high point.  Try to contract around
                    ! the low point.
                    do i = 1, npts
                        if (i /= ilo) then
                            pcent = half * (this%m_simplex(:,i) + &
                                this%m_simplex(:,ilo))
                            this%m_simplex(:,i) = pcent
                            f(i) = fcn%fcn(pcent, args)
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
                print 101, "Iteration: ", iter
                print 101, "Function Evaluations: ", neval
                print 102, "Function Value: ", fval
                print 102, "Convergence Parameter: ", rtol
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
            write(errmsg, 103) &
                "The algorithm failed to converge." // new_line('c') // &
                "Function evaluations performed: ", neval, new_line('c') // &
                "Convergence Parameter: ", rtol, new_line('c') // &
                "Convergence Criteria: ", ftol
            call errmgr%report_error("nm_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if

        ! Formatting
100     format(A, I0, A, I0, A)
101     format(A, I0)
102     format(A, E10.3)
103     format(A, I0, A, E10.3, A, E10.3)
    end subroutine

! ------------------------------------------------------------------------------
    function nm_extrapolate(this, fcn, y, pcent, ihi, fac, neval, &
            work, args) result(ytry)
        !! Extrapolates by the specified factor through the simplex across
        !! from the largest point.  If the extrapolation results in a better
        !! estimate, the current high point is replaced with the new estimate.
        class(nelder_mead), intent(inout) :: this
            !! The [[nelder_mead]] object.
        class(fcnnvar_helper), intent(in) :: fcn
            !! The function to evaluate.
        real(real64), intent(inout), dimension(:) :: y
            !! An array containing the function values at each vertex.
        real(real64), intent(inout), dimension(:) :: pcent
            !! An array containing the centroid of vertex position information.
        integer(int32), intent(in) :: ihi
            !! The index of the largest magnitude vertex.
        real(real64), intent(in) :: fac
            !! A scaling factor.
        integer(int32), intent(inout) :: neval
            !! The number of function evaluations.
        real(real64), intent(out), dimension(:) :: work
            !! An N-element workspace array where N is the number of dimensions
            !! of the problem.
        class(*), intent(inout), optional :: args
            !! An optional argument to allow the user to communicate with
            !! the routine.
        real(real64) :: ytry
            !! The new function estimate.

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
         ytry = fcn%fcn(work, args)
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
    pure function nm_get_simplex(this) result(p)
        !! Gets an N-by-(N+1) matrix containing the current simplex.
        class(nelder_mead), intent(in) :: this
            !! The [[nelder_mead]] object.
        real(real64), allocatable, dimension(:,:) :: p
            !! The N-by-(N+1) matrix containing the simplex.  Each vertex of 
            !! the simplex is stored as its own column of this matrix.
        integer(int32) :: m, n
        if (allocated(this%m_simplex)) then
            m = size(this%m_simplex, 1)
            n = size(this%m_simplex, 2)
            allocate(p(m,n))
            p = this%m_simplex
        end if
    end function

! --------------------
    subroutine nm_set_simplex(this, x)
        !! Sets an N-by-(N+1) matrix as the current simplex.  Notice, if
        !! this matrix is different in size from the problem dimensionallity,
        !! the Nelder-Mead routine will replace it with an appropriately sized
        !! matrix.
        class(nelder_mead), intent(inout) :: this
            !! The [[nelder_mead]] object.
        real(real64), dimension(:,:) :: x
            !! The simplex matrix.  Each column of the matrix must contain the
            !! coordinates of each vertex of the simplex.
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
    pure function nm_get_size(this) result(x)
        !! Gets the size of the initial simplex that will be utilized by
        !! the Nelder-Mead algorithm in the event that the user does not supply
        !! a simplex geometry, or if the user supplies an invalid simplex
        !! geometry.
        class(nelder_mead), intent(in) :: this
            !! The [[nelder_mead]] object.
        real(real64) :: x
            !! The size of the simplex (length of an edge).
        x = this%m_initSize
    end function

! --------------------
    subroutine nm_set_size(this, x)
        !! Sets the size of the initial simplex that will be utilized by
        !! the Nelder-Mead algorithm in the event that the user does not supply
        !! a simplex geometry, or if the user supplies an invalid simplex
        !! geometry.
        class(nelder_mead), intent(inout) :: this
            !! The [[nelder_mead]] object.
        real(real64), intent(in) :: x
            !! The size of the simplex (length of an edge).
        this%m_initSize = x
    end subroutine

! ******************************************************************************
! LINE_SEARCH_OPTIMIZER
! ------------------------------------------------------------------------------
    subroutine lso_get_line_search(this, ls)
        !! Gets the line search module.
        class(line_search_optimizer), intent(in) :: this
            !! The [[line_search_optimizer]] object.
        class(line_search), intent(out), allocatable :: ls
            !! The [[line_search]] object.
        if (allocated(this%m_lineSearch)) &
            allocate(ls, source = this%m_lineSearch)
    end subroutine

! ----------------------
    subroutine lso_set_line_search(this, ls)
        !! Sets the line search module.
        class(line_search_optimizer), intent(inout) :: this
            !! The [[line_search_optimizer]] object.
        class(line_search), intent(in) :: ls
            !! The [[line_search]] object.
        if (allocated(this%m_lineSearch)) deallocate(this%m_lineSearch)
        allocate(this%m_lineSearch, source = ls)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine lso_set_default(this)
        !! Establishes a default line_search object for the line search
        !! module.
        class(line_search_optimizer), intent(inout) :: this
            !! The [[line_search_optimizer]] object.
        type(line_search) :: ls
        call this%set_line_search(ls)
    end subroutine

! ------------------------------------------------------------------------------
    pure function lso_is_line_search_defined(this) result(x)
        !! Tests to see if a line search module is defined.
        class(line_search_optimizer), intent(in) :: this
            !! The [[line_search_optimizer]] object.
        logical :: x
            !! Returns true if a module is defined; else, false.
        x = allocated(this%m_lineSearch)
    end function

! ------------------------------------------------------------------------------
    pure function lso_get_use_search(this) result(x)
        !! Gets a value determining if a line-search should be employed.
        class(line_search_optimizer), intent(in) :: this
            !! The [[line_search_optimizer]] object.
        logical :: x
            !! Returns true if a line search should be used; else, false.
        x = this%m_useLineSearch
    end function

! --------------------
    subroutine lso_set_use_search(this, x)
        !! Sets a value determining if a line-search should be employed.
        class(line_search_optimizer), intent(inout) :: this
            !! The [[line_search_optimizer]] object.
        logical, intent(in) :: x
            !! Set to true if a line search should be used; else, false.
        this%m_useLineSearch = x
    end subroutine

! ------------------------------------------------------------------------------
    pure function lso_get_var_tol(this) result(x)
        !! Gets the convergence on change in variable tolerance.
        class(line_search_optimizer), intent(in) :: this
            !! The [[line_search_optimizer]] object.
        real(real64) :: x
            !! The tolerance value.
        x = this%m_xtol
    end function

! --------------------
    subroutine lso_set_var_tol(this, x)
        !! Sets the convergence on change in variable tolerance.
        class(line_search_optimizer), intent(inout) :: this
            !! The [[line_search_optimizer]] object.
        real(real64), intent(in) :: x
            !! The tolerance value.
        this%m_xtol = x
    end subroutine

! ******************************************************************************
! BFGS
! ------------------------------------------------------------------------------
    subroutine bfgs_solve(this, fcn, x, fout, ib, args, err)
        !! Utilizes the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
        !! for finding a minimum value of the specified function.
        class(bfgs), intent(inout) :: this
            !! The [[bfgs]] object.
        class(fcnnvar_helper), intent(in) :: fcn
            !! The [[fcnnvar_helper]] object containing the equation to 
            !! optimize.
        real(real64), intent(inout), dimension(:) :: x
            !! On input, the initial guess at the optimal point.  On output, 
            !! the updated optimal point estimate.
        real(real64), intent(out), optional :: fout
            !! An optional output, that if provided, returns the value of the 
            !! function at x.
        type(iteration_behavior), optional :: ib
            !! An optional output, that if provided, allows the caller to 
            !! obtain iteration performance statistics.
        class(*), intent(inout), optional :: args
            !! An optional argument to allow the user to communicate with fcn.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0
        real(real64), parameter :: negone = -1.0d0
        real(real64), parameter :: factor = 1.0d2
        real(real64), parameter :: small = 1.0d-10

        ! Local Variables
        logical :: xcnvrg, gcnvrg
        integer(int32) :: i, n, maxeval, neval, ngrad, flag, iter
        real(real64) :: xtol, gtol, fp, stpmax, fret, xtest, gtest, temp, ydx
        real(real64), allocatable, dimension(:) :: g, dx, u, v, y, gold, xnew, bdx
        real(real64), allocatable, dimension(:,:) :: b, r
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        type(iteration_behavior) :: lib
        class(line_search), allocatable :: ls

        ! Initialization
        n = fcn%get_variable_count()
        maxeval = this%get_max_fcn_evals()
        gtol = this%get_tolerance()
        xtol = this%get_var_tolerance()
        iter = 0
        neval = 0
        ngrad = 0
        xcnvrg = .false.
        gcnvrg = .false.
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = 0
            ib%gradient_count = ngrad
            ib%converge_on_fcn = .false.
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = gcnvrg
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        if (this%get_use_line_search()) then
            if (.not.this%is_line_search_defined()) &
                call this%set_default_line_search()
            call this%get_line_search(ls)
        end if

        ! Input Check
        if (.not.fcn%is_fcn_defined()) then
            ! ERROR: No function is defined
            call errmgr%report_error("bfgs_solve", &
                "No function has been defined.", &
                NL_INVALID_OPERATION_ERROR)
            return
        end if
        if (size(x) /= n) then
            write(errmsg, 100) &
                "It was expected to receive a coordinate vector of length ", &
                n, " , but a vector of length ", size(x), " was received."
            call errmgr%report_error("bfgs_solve", trim(errmsg), &
                NL_INVALID_INPUT_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(g(n), stat = flag)
        if (flag == 0) allocate(dx(n), stat = flag)
        if (flag == 0) allocate(u(n), stat = flag)
        if (flag == 0) allocate(v(n), stat = flag)
        if (flag == 0) allocate(y(n), stat = flag)
        if (flag == 0) allocate(bdx(n), stat = flag)
        if (flag == 0) allocate(gold(n), stat = flag)
        if (flag == 0) allocate(xnew(n), stat = flag)
        if (flag == 0) allocate(b(n,n), stat = flag)
        if (flag == 0) allocate(r(n,n), stat = flag)
        if (flag /= 0) then
            ! ERROR: Memory Error
            call errmgr%report_error("bfgs_solve", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Process
        fp = fcn%fcn(x, args)
        call fcn%gradient(x, g, fv = fp, args = args)
        neval = 1
        ngrad = 1

        ! Check for a "zero" gradient at the initial point
        gtest = norm2(g)
        if (gtest < gtol) then
            gcnvrg = .true.
        end if

        ! Main Loop
        flag = 0
        if (.not.gcnvrg) then
            do
                ! Update the iteration counter
                iter = iter + 1

                ! Define the initial direction, and a limit on the line search
                ! step
                if (iter == 1) then
                    dx = -g
                    stpmax = factor * max(norm2(x), real(n, real64))
                end if

                ! Perform the line search
                if (this%get_use_line_search()) then
                    call limit_search_vector(dx, stpmax)
                    call ls%search(fcn, x, g, dx, xnew, fp, fret, lib, &
                        args = args, err = errmgr)
                    neval = neval + lib%fcn_count
                    fp = fret
                else
                    xnew = x + dx
                    fp = fcn%fcn(xnew, args)
                    neval = neval + 1
                end if

                ! Update the gradient and line direction
                do i = 1, n
                    dx(i) = xnew(i) - x(i)
                    x(i) = xnew(i)
                    gold(i) = g(i)
                end do
                call fcn%gradient(x, g, fv = fp, args = args)
                ngrad = ngrad + 1

                ! Test for convergence on the change in X
                xtest = zero
                do i = 1, n
                    temp = abs(dx(i)) / max(abs(x(i)), one)
                    xtest = max(temp, xtest)
                end do
                if (xtest < xtol) then
                    xcnvrg = .true.
                    exit
                end if

                ! Test for convergence on the gradient
                gtest = norm2(g)
                if (gtest < gtol) then
                    gcnvrg = .true.
                    exit
                end if

                ! Perform the BFGS update
                y = g - gold
                ydx = dot_product(y, dx)

                ! Establish an initial approximation to the Hessian matrix
                if (iter == 1) then
                    temp = sqrt(dot_product(y, y) / ydx)
                    call dlaset('A', n, n, zero, temp, r, n)
                end if

                ! Compute: B = R**T * R
                call tri_mtx_mult(.true., one, r, zero, b)

                ! Compute bdx = B * dX (B is symmetric)
                call dsymv('u', n, one, b, n, dx, 1, zero, bdx, 1)

                ! Perform the actual update
                if (ydx > small) then
                    ! Compute the rank 1 update and downdate
                    u = y / sqrt(ydx)
                    v = bdx / sqrt(dot_product(dx, bdx))
                    call cholesky_rank1_update(r, u)
                    call cholesky_rank1_downdate(r, v)
                end if ! Else just skip the update

                ! Compute the solution to: B * dx = -g = (R**T * R) * dx
                dx = -g
                call solve_cholesky(.true., r, dx)

                ! Print iteration status
                if (this%get_print_status()) then
                    print *, ""
                    print 101, "Iteration: ", iter
                    print 101, "Function Evaluations: ", neval
                    print 102, "Function Value: ", fp
                    print 102, "Change in Variable: ", xtest
                    print 102, "Gradient: ", gtest
                end if

                ! Ensure we haven't made too many function evaluations
                if (neval >= maxeval) then
                    flag = 1
                    exit
                end if
            end do
        end if

        ! Report out iteration statistics
        if (present(ib)) then
            ib%iter_count = iter
            ib%fcn_count = neval
            ib%jacobian_count = 0
            ib%gradient_count = ngrad
            ib%converge_on_fcn = .false.
            ib%converge_on_chng = xcnvrg
            ib%converge_on_zero_diff = gcnvrg
        end if

        ! Get the function value at the computed minimum
        if (present(fout)) fout = fp

        ! Check for convergence issues
        if (flag /= 0) then
            write(errmsg, 103) &
                "The algorithm failed to converge." // new_line('c') // &
                "Function evaluations performed: ", neval, new_line('c') // &
                "Function Value: ", fp, new_line('c') // &
                "Change in Variable: ", xtest, new_line('c') // &
                "Gradient: ", gtest
            call errmgr%report_error("bfgs_solve", trim(errmsg), &
                NL_CONVERGENCE_ERROR)
        end if

        ! Formatting
100     format(A, I0, A, I0, A)
101     format(A, I0)
102     format(A, E10.3)
103     format(A, I0, A, E10.3, A, E10.3, A, E10.3)
    end subroutine

! ------------------------------------------------------------------------------
end module
