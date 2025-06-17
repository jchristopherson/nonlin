module nonlin_polynomials
    use iso_fortran_env
    use linalg, only : eigen, solve_least_squares
    use ferror, only : errors
    use nonlin_error_handling
    implicit none
    private
    public :: polynomial
    public :: assignment(=)
    public :: operator(+)
    public :: operator(-)
    public :: operator(*)

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface assignment(=)
        !! Defines polynomial assignment.
        module procedure :: poly_equals
        module procedure :: poly_dbl_equals
        module procedure :: poly_equals_array
    end interface

    interface operator(+)
        !! Defines polynomial addition.
        module procedure :: poly_poly_add
    end interface

    interface operator(-)
        !! Defines polynomial subtraction.
        module procedure :: poly_poly_subtract
    end interface

    interface operator(*)
        !! Defines polynomial multiplication
        module procedure :: poly_poly_mult
        module procedure :: poly_dbl_mult
        module procedure :: dbl_poly_mult
    end interface

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    type polynomial
        !! Defines a polynomial, and associated routines for performing
        !! polynomial operations.
        real(real64), private, allocatable, dimension(:) :: m_coeffs
            !! An array that contains the polynomial coefficients in ascending 
            !! order.
    contains
        generic, public :: initialize => init_poly, init_poly_coeffs
        procedure, public :: order => get_poly_order
        procedure, public :: fit => poly_fit
        procedure, public :: fit_thru_zero => poly_fit_thru_zero
        generic, public :: evaluate => evaluate_real, evaluate_complex
        procedure, public :: companion_mtx => poly_companion_mtx
        procedure, public :: roots => poly_roots
        procedure, public :: get => get_poly_coefficient
        procedure, public :: get_all => get_poly_coefficients
        procedure, public :: set => set_poly_coefficient

        procedure, private :: evaluate_real => poly_eval_double
        procedure, private :: evaluate_complex => poly_eval_complex
        procedure, private :: init_poly
        procedure, private :: init_poly_coeffs
    end type

contains
! ******************************************************************************
! POLYNOMIAL MEMBERS
! ------------------------------------------------------------------------------
    subroutine init_poly(this, order, err)
        !! Initializes the polynomial instance.
        class(polynomial), intent(inout) :: this
            !! The [[polynomial]] object.
        integer(int32), intent(in) :: order
            !! The order of the polynomial (must be >= 0).
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: n, istat
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = order + 1
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (order < 0) then
            ! ERROR: Negative order is not supported
            call errmgr%report_error("init_polynomial", &
                "A negative polynomial order is not supported.", &
                NL_INVALID_INPUT_ERROR)
            return
        end if

        ! Process
        if (allocated(this%m_coeffs)) deallocate(this%m_coeffs)
        allocate(this%m_coeffs(n), stat = istat)
        if (istat /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("init_polynomial", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if
        this%m_coeffs = zero
    end subroutine

! ------------------------------------------------------------------------------
    subroutine init_poly_coeffs(this, c, err)
        !! Initializes the polynomial instance.
        class(polynomial), intent(inout) :: this
            !! The [[polynomial]] object.
        real(real64), intent(in), dimension(:) :: c
            !! The array of polynomial coefficients. The coefficients are
            !! established as follows: c(1) + c(2) * x + c(3) * x**2 + ...
            !! c(n) * x**n-1.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Local Variables
        integer(int32) :: i, n
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(c)
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Initialize the polynomial
        call init_poly(this, n - 1, errmgr)
        if (errmgr%has_error_occurred()) return

        ! Populate the polynomial coefficients
        do i = 1, n
            call this%set(i, c(i))
        end do
    end subroutine

! ------------------------------------------------------------------------------
    pure function get_poly_order(this) result(n)
        !! Returns the order of the polynomial object.
        class(polynomial), intent(in) :: this
            !! The [[polynomial]] object.
        integer(int32) :: n
            !! The order of the polynomial.  Returns -1 in the event no
            !! polynomial coefficients have been defined.
        if (.not.allocated(this%m_coeffs)) then
            n = -1
        else
            n = size(this%m_coeffs) - 1
        end if
    end function

! ------------------------------------------------------------------------------
    subroutine poly_fit(this, x, y, order, err)
        !! Fits a polynomial of the specified order to the supplied data set.
        class(polynomial), intent(inout) :: this
            !! The [[polynomial]] object.
        real(real64), intent(in), dimension(:) :: x
            !! An N-element array containing the independent variable data
            !! points.  Notice, must be N > order.
        real(real64), intent(inout), dimension(:) :: y
            !! On input, an N-element array containing the dependent variable 
            !! data points.  On output, the contents are overwritten.
        integer(int32), intent(in) :: order
            !! The order of the polynomial (must be >= 1).
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Parameters
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32):: j, n, ncols, flag
        real(real64), pointer, dimension(:,:) :: a
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(x)
        ncols = order + 1
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(y) /= n) then
            ! ERROR: Array size mismatch
            call errmgr%report_error("polynomial_fit", "Array size mismatch.", &
                NL_ARRAY_SIZE_ERROR)
            return
        else if (order >= n .or. order < 1) then
            ! ERROR: Requested order does not make sense
            call errmgr%report_error("polynomial_fit", "The requested " // &
                "polynomial order is not valid for this data set.", &
                NL_INVALID_INPUT_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(a(n,ncols), stat = flag)
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("polynomial_fit", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Ensure the polynomial object is initialized and sized appropriately
        if (this%order() /= order) then
            call this%initialize(order, errmgr)
            if (errmgr%has_error_occurred()) return
        end if

        ! Populate A
        do j = 1, n
            a(j,1) = one
            a(j,2) = x(j)
        end do
        do j = 3, ncols
            a(:,j) = a(:,j-1) * x
        end do

        ! Solve: A * coeffs = y
        call solve_least_squares(a, y, err = errmgr)
        if (errmgr%has_error_occurred()) return

        ! Extract the coefficients from the first order+1 elements of Y
        this%m_coeffs = y(1:ncols)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine poly_fit_thru_zero(this, x, y, order, err)
        !! Fits a polynomial of the specified order that passes through zero
        !! to the supplied data set.
        class(polynomial), intent(inout) :: this
            !! The [[polynomial]] object.
        real(real64), intent(in), dimension(:) :: x
            !! An N-element array containing the independent variable data
            !! points.  Notice, must be N > order.
        real(real64), intent(inout), dimension(:) :: y
            !! On input, an N-element array containing the dependent
            !! variable data points.  On output, the contents are overwritten.
        integer(int32), intent(in) :: order
            !! The order of the polynomial (must be >= 1).
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32):: j, n, ncols, flag
        real(real64), pointer, dimension(:,:) :: a
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        n = size(x)
        ncols = order
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Check
        if (size(y) /= n) then
            ! ERROR: Array size mismatch
            call errmgr%report_error("polynomial_fit_thru_zero", &
                "Array size mismatch.", NL_ARRAY_SIZE_ERROR)
            return
        else if (order >= n .or. order < 1) then
            ! ERROR: Requested order does not make sense
            call errmgr%report_error("polynomial_fit_thru_zero", &
                "The requested polynomial order is not valid for this " // &
                "data set.", NL_INVALID_INPUT_ERROR)
            return
        end if

        ! Local Memory Allocation
        allocate(a(n,ncols), stat = flag)
        if (flag /= 0) then
            ! ERROR: Out of memory
            call errmgr%report_error("polynomial_fit_thru_zero", &
                "Insufficient memory available.", NL_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Ensure the polynomial object is initialized and sized appropriately
        if (this%order() /= order) then
            call this%initialize(order, errmgr)
            if (errmgr%has_error_occurred()) return
        end if

        ! Populate A
        a(:,1) = x
        do j = 2, ncols
            a(:,j) = a(:,j-1) * x
        end do

        ! Solve: A * coeffs = y
        call solve_least_squares(a, y, err = errmgr)
        if (errmgr%has_error_occurred()) return

        ! Extract the coefficients from the first order+1 elements of Y
        this%m_coeffs(1) = zero
        this%m_coeffs(2:ncols+1) = y(1:ncols)
    end subroutine

! ------------------------------------------------------------------------------
    elemental function poly_eval_double(this, x) result(y)
        !! Evaluates a polynomial at the specified points.
        class(polynomial), intent(in) :: this
            !! The [[polynomial]] object.
        real(real64), intent(in) :: x
            !! The value(s) at which to evaluate the polynomial.
        real(real64) :: y
            !! The value(s) of the polynomial at x.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: j, order, n

        ! Initialization
        order = this%order()
        n = order + 1
        if (order == -1) then
            y = zero
            return
        else if (order == 0) then
            y = this%m_coeffs(1)
            return
        end if

        ! Process
        y =this%m_coeffs(n) * x + this%m_coeffs(order)
        do j = n - 2, 1, -1
            y = y * x + this%m_coeffs(j)
        end do
    end function

! ------------------------------------------------------------------------------
    elemental function poly_eval_complex(this, x) result(y)
        !! Evaluates a polynomial at the specified points.
        class(polynomial), intent(in) :: this
            !! The [[polynomial]] object.
        complex(real64), intent(in) :: x
            !! The value(s) at which to evaluate the polynomial.
        complex(real64) :: y
            !! The value(s) of the polynomial at x.

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: j, order, n

        ! Initialization
        order = this%order()
        n = order + 1
        if (order == -1) then
            y = zero
            return
        else if (order == 0) then
            y = this%m_coeffs(1)
            return
        end if

        ! Process
        y =this%m_coeffs(n) * x + this%m_coeffs(order)
        do j = n - 2, 1, -1
            y = y * x + this%m_coeffs(j)
        end do
    end function

! ------------------------------------------------------------------------------
    pure function poly_companion_mtx(this) result(c)
        !! Returns the companion matrix for the polynomial.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Companion_matrix)
        !!
        !! - [Wolfram MathWorld](http://mathworld.wolfram.com/CompanionMatrix.html)
        class(polynomial), intent(in) :: this
            !! The [[polynomial]] object.
        real(real64), dimension(this%order(), this%order()) :: c
            !! The companion matrix.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Local Variables
        integer(int32) :: i, n

        ! Process
        n = this%order()
        if (n == -1) return
        c = zero
        do i = 1, n
            c(i,n) = -this%m_coeffs(i) / this%m_coeffs(n + 1)
            if (i < n) c(i+1,i) = one
        end do
    end function

! ------------------------------------------------------------------------------
    function poly_roots(this, err) result(z)
        !! Computes all the roots of a polynomial by computing the eigenvalues
        !! of the polynomial companion matrix.
        class(polynomial), intent(in) :: this
            !! The [[polynomial]] object.
        complex(real64), dimension(this%order()) :: z
            !! The roots of the polynomial.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Local Variables
        integer(int32) :: n
        real(real64), allocatable, dimension(:,:) :: c

        ! Initialization
        n = this%order()

        ! Quick Return
        if (n == 0) return

        ! Compute the companion matrix
        c = this%companion_mtx()

        ! Compute the eigenvalues of the companion matrix.  The eigenvalues are
        ! the roots of the polynomial.
        call eigen(c, z, err = err)
    end function

! ------------------------------------------------------------------------------
    function get_poly_coefficient(this, ind, err) result(c)
        !! Gets the requested polynomial coefficient by index.  The
        !! coefficient index is established as follows: c(1) + c(2) * x +
        !! c(3) * x**2 + ... c(n) * x**n-1.
        class(polynomial), intent(in) :: this
            !! The [[polynomial]] object.
        integer(int32), intent(in) :: ind
            !! The polynomial coefficient index.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.
        real(real64) :: c
            !! The requested coefficient.

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        c = 0.0d0
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Quick Return
        if (this%order() == -1) return

        ! Input Check
        if (ind <= 0 .or. ind > this%order() + 1) then
            ! ERROR: Index out of range
            call errmgr%report_error("get_polynomial_coefficient", &
                "The specified index is outside the bounds of the " // &
                "coefficient array.", NL_INVALID_INPUT_ERROR)
            return
        end if

        ! Get the coefficient
        c = this%m_coeffs(ind)
    end function

! ------------------------------------------------------------------------------
    pure function get_poly_coefficients(this) result(c)
        !! Gets an array containing all the coefficients of the polynomial.
        !! The coefficient index is established as follows: c(1) + c(2) * x +
        !! c(3) * x**2 + ... c(n) * x**n-1.
        class(polynomial), intent(in) :: this
            !! The [[polynomial]] object.
        real(real64), dimension(this%order() + 1) :: c
            !! The array of coefficients.

        ! Process
        if (this%order() == -1) return
        c = this%m_coeffs
    end function

! ------------------------------------------------------------------------------
    subroutine set_poly_coefficient(this, ind, c, err)
        !! Sets the requested polynomial coefficient by index.  The
        !! coefficient index is established as follows: c(1) + c(2) * x +
        !! c(3) * x**2 + ... c(n) * x**n-1.
        class(polynomial), intent(inout) :: this
            !! The [[polynomial]] object.
        integer(int32), intent(in) :: ind
            !! The polynomial coefficient index.
        real(real64), intent(in) :: c
            !! The polynomial coefficient.
        class(errors), intent(inout), optional, target :: err
            !! An error handling object.

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Quick Return
        if (this%order() == -1) return

        ! Input Check
        if (ind <= 0 .or. ind > this%order() + 1) then
            ! ERROR: Index out of range
            call errmgr%report_error("set_polynomial_coefficient", &
                "The specified index is outside the bounds of the " // &
                "coefficient array.", NL_INVALID_INPUT_ERROR)
            return
        end if

        ! Process
        this%m_coeffs(ind) = c
    end subroutine

! ******************************************************************************
! OPERATORS
! ------------------------------------------------------------------------------
    subroutine poly_equals(x, y)
        !! Assigns the contents of one polynomial to another.
        class(polynomial), intent(inout) :: x
            !! The assignee.
        class(polynomial), intent(in) :: y
            !! The item to copy.

        ! Local Variables
        integer(int32) :: i, ord

        ! Process
        ord = y%order()
        if (x%order() /= ord) call x%initialize(ord)
        do i = 1, ord + 1
            call x%set(i, y%get(i))
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine poly_dbl_equals(x, y)
        !! Assigns a number to each coefficient of the polynomial.
        class(polynomial), intent(inout) :: x
            !! The assignee.
        real(real64), intent(in) :: y
            !! The value to assign.

        ! Local Variables
        integer(int32) :: i, ord

        ! Process
        ord = x%order()
        do i = 1, ord + 1
            call x%set(i, y)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine poly_equals_array(x, y)
        !! Assigns the contents of an array as polynomial coefficients.
        class(polynomial), intent(inout) :: x
            !! The assignee.
        real(real64), intent(in), dimension(:) :: y
            !! The coefficient array.
        call x%initialize(y)
    end subroutine

! ------------------------------------------------------------------------------
    function poly_poly_add(x, y) result(z)
        !! Adds two polynomials.
        class(polynomial), intent(in) :: x
            !! The left-hand-side argument.
        class(polynomial), intent(in) :: y
            !! The right-hand-side argument.
        type(polynomial) :: z
            !! The resulting polynomial.

        ! Local Variables
        integer(int32) :: i, max_ord, x_ord, y_ord

        ! Initialization
        x_ord = x%order()
        y_ord = y%order()
        max_ord = max(x_ord, y_ord)
        call z%initialize(max_ord)

        ! Quick Return
        if (x_ord == -1 .and. y_ord == -1) return
        if (x_ord == -1 .and. y_ord /= -1) then
            do i = 1, max_ord + 1
                call z%set(i, y%get(i))
            end do
            return
        else if (x_ord /= -1 .and. y_ord == -1) then
            do i = 1, max_ord + 1
                call z%set(i, x%get(i))
            end do
            return
        end if

        ! Process
        if (x_ord > y_ord) then
            do i = 1, y_ord + 1
                call z%set(i, x%get(i) + y%get(i))
            end do
            do i = y_ord + 2, x_ord
                call z%set(i, x%get(i))
            end do
        else if (x_ord < y_ord) then
            do i = 1, x_ord + 1
                call z%set(i, x%get(i) + y%get(i))
            end do
            do i = x_ord + 2, y_ord + 1
                call z%set(i, y%get(i))
            end do
        else
            do i = 1, max_ord + 1
                call z%set(i, x%get(i) + y%get(i))
            end do
        end if
    end function

! ------------------------------------------------------------------------------
    function poly_poly_subtract(x, y) result(z)
        !! Subtracts two polynomials.
        class(polynomial), intent(in) :: x
            !! The left-hand-side argument.
        class(polynomial), intent(in) :: y
            !! The right-hand-side argument.
        type(polynomial) :: z
            !! The resulting polynomial.

        ! Local Variables
        integer(int32) :: i, max_ord, x_ord, y_ord

        ! Initialization
        x_ord = x%order()
        y_ord = y%order()
        max_ord = max(x_ord, y_ord)
        call z%initialize(max_ord)

        ! Quick Return
        if (x_ord == -1 .and. y_ord == -1) return
        if (x_ord == -1 .and. y_ord /= -1) then
            do i = 1, max_ord + 1
                call z%set(i, y%get(i))
            end do
            return
        else if (x_ord /= -1 .and. y_ord == -1) then
            do i = 1, max_ord + 1
                call z%set(i, x%get(i))
            end do
            return
        end if

        ! Process
        if (x_ord > y_ord) then
            do i = 1, y_ord + 1
                call z%set(i, x%get(i) - y%get(i))
            end do
            do i = y_ord + 2, x_ord
                call z%set(i, x%get(i))
            end do
        else if (x_ord < y_ord) then
            do i = 1, x_ord + 1
                call z%set(i, x%get(i) - y%get(i))
            end do
            do i = x_ord + 2, y_ord + 1
                call z%set(i, -y%get(i))
            end do
        else
            do i = 1, max_ord + 1
                call z%set(i, x%get(i) - y%get(i))
            end do
        end if
    end function

! ------------------------------------------------------------------------------
    function poly_poly_mult(x, y) result(z)
        !! Multiplies two polynomials.
        class(polynomial), intent(in) :: x
            !! The left-hand-side argument.
        class(polynomial), intent(in) :: y
            !! The right-hand-side argument.
        type(polynomial) :: z
            !! The resulting polynomial.

        ! Local Variables
        integer(int32) :: i, j, m, n
        real(real64) :: val

        ! Initialization
        n = x%order() + 1
        m = y%order() + 1
        call z%initialize(x%order() + y%order()) ! Sets z to all zeros

        ! Process
        do i = 1, n
            do j = 1, m
                val = z%get(i + j - 1) + x%get(i) * y%get(j)
                call z%set(i + j - 1, val)
            end do
        end do
    end function

! ------------------------------------------------------------------------------
    function poly_dbl_mult(x, y) result(z)
        !! Multiplies a polynomial by a scalar value.
        class(polynomial), intent(in) :: x
            !! The left-hand-side argument.
        real(real64), intent(in) :: y
            !! The right-hand-side argument.
        type(polynomial) :: z
            !! The resulting polynomial.

        ! Local Variables
        integer(int32) :: i, ord

        ! Process
        ord = x%order()
        call z%initialize(ord)
        do i = 1, ord + 1
            call z%set(i, x%get(i) * y)
        end do
    end function

! ------------------------------------------------------------------------------
    function dbl_poly_mult(x, y) result(z)
        !! Multiplies a polynomial by a scalar value.
        real(real64), intent(in) :: x
            !! The left-hand-side argument.
        class(polynomial), intent(in) :: y
            !! The right-hand-side argument.
        type(polynomial) :: z
            !! The resulting polynomial.

        ! Local Variables
        integer(int32) :: i, ord

        ! Process
        ord = y%order()
        call z%initialize(ord)
        do i = 1, ord + 1
            call z%set(i, y%get(i) * x)
        end do
    end function

! ------------------------------------------------------------------------------

! Example Polynomial Code (Coefficients go from lowest order to highest)
! src: https://github.com/JuliaMath/Polynomials.jl
!
! function *{T,S}(p1::Poly{T}, p2::Poly{S})
!     if p1.var != p2.var
!         error("Polynomials must have same variable")
!     end
!     R = promote_type(T,S)
!     n = length(p1)-1
!     m = length(p2)-1
!     a = zeros(R,m+n+1)

!     for i = 0:n
!         for j = 0:m
!             a[i+j+1] += p1[i] * p2[j]
!         end
!     end
!     Poly(a,p1.var)
! end

! ## older . operators, hack to avoid warning on v0.6
! dot_operators = quote
!     @compat Base.:.+{T<:Number}(c::T, p::Poly) = +(p, c)
!     @compat Base.:.+{T<:Number}(p::Poly, c::T) = +(p, c)
!     @compat Base.:.-{T<:Number}(p::Poly, c::T) = +(p, -c)
!     @compat Base.:.-{T<:Number}(c::T, p::Poly) = +(p, -c)
!     @compat Base.:.*{T<:Number,S}(c::T, p::Poly{S}) = Poly(c * p.a, p.var)
!     @compat Base.:.*{T<:Number,S}(p::Poly{S}, c::T) = Poly(p.a * c, p.var)
! end
! VERSION < v"0.6.0-dev" && eval(dot_operators)


! # are any values NaN
! hasnan(p::Poly) = reduce(|, (@compat isnan.(p.a)))

! function divrem{T, S}(num::Poly{T}, den::Poly{S})
!     if num.var != den.var
!         error("Polynomials must have same variable")
!     end
!     m = length(den)-1
!     if m == 0 && den[0] == 0
!         throw(DivideError())
!     end
!     R = typeof(one(T)/one(S))
!     n = length(num)-1
!     deg = n-m+1
!     if deg <= 0
!         return convert(Poly{R}, zero(num)), convert(Poly{R}, num)
!     end

!     aQ = zeros(R, deg)
!     # aR = deepcopy(num.a)
!     # @show num.a
!     aR = R[ num.a[i] for i = 1:n+1 ]
!     for i = n:-1:m
!         quot = aR[i+1] / den[m]
!         aQ[i-m+1] = quot
!         for j = 0:m
!             elem = den[j]*quot
!             aR[i-(m-j)+1] -= elem
!         end
!     end
!     pQ = Poly(aQ, num.var)
!     pR = Poly(aR, num.var)

!     return pQ, pR
! end

! div(num::Poly, den::Poly) = divrem(num, den)[1]
! rem(num::Poly, den::Poly) = divrem(num, den)[2]

! ==(p1::Poly, p2::Poly) = (p1.var == p2.var && p1.a == p2.a)
! ==(p1::Poly, n::Number) = (coeffs(p1) == [n])
! ==(n::Number, p1::Poly) = (p1 == n)
end module
