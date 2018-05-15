! nonlin_polynomials.f90


! TO DO:
! - Multiplication (convolution)
! - Division (deconvolution)
! - C interface

!> @brief \b polynomials
!!
!! @par Purpose
!! Provides a means of defining and operating on polynomials.
module nonlin_polynomials
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use linalg_core, only : eigen, solve_least_squares
    use ferror, only : errors
    use nonlin_constants, only : NL_INVALID_INPUT_ERROR, NL_ARRAY_SIZE_ERROR, &
        NL_OUT_OF_MEMORY_ERROR
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
!> @brief Defines polynomial assignment.
interface assignment(=)
    module procedure :: poly_equals
    module procedure :: poly_dbl_equals
    module procedure :: poly_equals_array
end interface

!> @brief Defines polynomial addition.
interface operator(+)
    module procedure :: poly_poly_add
end interface

!> @brief Defines polynomial subtraction.
interface operator(-)
    module procedure :: poly_poly_subtract
end interface

!> @brief Defines polynomial multiplication
interface operator(*)
    module procedure :: poly_poly_mult
    module procedure :: poly_dbl_mult
    module procedure :: dbl_poly_mult
end interface

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
!> @brief Defines a polynomial, and associated routines for performing
!! polynomial operations.
type polynomial
private
    !> An array that contains the polynomial coefficients in ascending order.
    real(real64), allocatable, dimension(:) :: m_coeffs
contains
    !> @brief Initializes the polynomial instance.
    generic, public :: initialize => init_poly, init_poly_coeffs
    !> @brief Returns the order of the polynomial object.
    procedure, public :: order => get_poly_order
    !> @brief Fits a polynomial of the specified order to a data set.
    procedure, public :: fit => poly_fit
    !> @brief Fits a polynomial of the specified order that passes through zero
    !! to a data set.
    procedure, public :: fit_thru_zero => poly_fit_thru_zero
    !> @brief Evaluates a polynomial at the specified points.
    generic, public :: evaluate => evaluate_real, evaluate_complex
    !> @brief Returns the companion matrix for the polynomial.
    procedure, public :: companion_mtx => poly_companion_mtx
    !> @brief Computes all the roots of a polynomial.
    procedure, public :: roots => poly_roots
    !> @brief Gets the requested polynomial coefficient.
    procedure, public :: get => get_poly_coefficient
    !> @brief Gets an array containing all the coefficients of the polynomial.
    procedure, public :: get_all => get_poly_coefficients
    !> @brief Sets the requested polynomial coefficient by index.
    procedure, public :: set => set_poly_coefficient

    procedure :: evaluate_real => poly_eval_double
    procedure :: evaluate_complex => poly_eval_complex
    procedure :: init_poly
    procedure :: init_poly_coeffs
end type

contains
! ******************************************************************************
! POLYNOMIAL MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Initializes the polynomial instance, and sets all coefficients
    !! to zero.
    !!
    !! @param[in,out] this The polynomial object.
    !! @param[in] order The order of the polynomial (must be >= 0).
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if a zero or negative polynomial order
    !!      was specified.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
    subroutine init_poly(this, order, err)
        ! Arguments
        class(polynomial), intent(inout) :: this
        integer(int32), intent(in) :: order
        class(errors), intent(inout), optional, target :: err

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
    !> @brief Initializes the polynomial instance.
    !!
    !! @param[in,out] this The polynomial object.
    !! @param[in] c The array of polynomial coefficients. The coefficients are
    !!  established as follows: c(1) + c(2) * x + c(3) * x**2 + ... c(n) * x**n-1.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if a zero or negative polynomial order
    !!      was specified.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
    subroutine init_poly_coeffs(this, c, err)
        ! Arguments
        class(polynomial), intent(inout) :: this
        real(real64), intent(in), dimension(:) :: c
        class(errors), intent(inout), optional, target :: err

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
    !> @brief Returns the order of the polynomial object.
    !!
    !! @param[in] this The polynomial object.
    !!
    !! @return The order of the polynomial.  Returns -1 in the event no
    !! polynomial coefficients have been defined.
    pure function get_poly_order(this) result(n)
        class(polynomial), intent(in) :: this
        integer(int32) :: n
        if (.not.allocated(this%m_coeffs)) then
            n = -1
        else
            n = size(this%m_coeffs) - 1
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Fits a polynomial of the specified order to a data set.
    !!
    !! @param[in,out] this The polynomial object.
    !! @param[in] x An N-element array containing the independent variable data
    !!  points.  Notice, must be N > @p order.
    !! @param[in,out] y On input, an N-element array containing the dependent
    !!  variable data points.  On output, the contents are overwritten.
    !! @param[in] order The order of the polynomial (must be >= 1).
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if a zero or negative polynomial order
    !!      was specified, or if order is too large for the data set.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
    !!  - NL_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are different sizes.
    !!
    !! @par Usage
    !! The following code provides an example of how to fit a polynomial to a
    !! set of data.
    !! @code{.f90}
    !! program example
    !!     use linalg_constants, only : dp, i32
    !!     use nonlin_polynomials
    !!
    !!     ! Local Variables
    !!     real(real64), dimension(21) :: xp, yp, yf, yc, err
    !!     real(real64) :: res
    !!     type(polynomial) :: p
    !!
    !!     ! Data to fit
    !!     xp = [0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, &
    !!         0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, &
    !!         1.8d0, 1.9d0, 2.0d0]
    !!     yp = [1.216737514d0, 1.250032542d0, 1.305579195d0, 1.040182335d0, &
    !!         1.751867738d0, 1.109716707d0, 2.018141531d0, 1.992418729d0, &
    !!         1.807916923d0, 2.078806005d0, 2.698801324d0, 2.644662712d0, &
    !!         3.412756702d0, 4.406137221d0, 4.567156645d0, 4.999550779d0, &
    !!         5.652854194d0, 6.784320119d0, 8.307936836d0, 8.395126494d0, &
    !!         10.30252404d0]
    !!
    !!     ! Create a copy of yp as it will be overwritten in the fit command
    !!     yc = yp
    !!
    !!     ! Fit the polynomial
    !!     call p%fit(xp, yp, 3)
    !!
    !!     ! Evaluate the polynomial at xp, and then determine the residual
    !!     yf = p%evaluate(xp)
    !!     err = abs(yf - yc)
    !!     res = maxval(err)
    !!
    !!     ! Print out the coefficients
    !!      print '(A)', "Polynomial Coefficients (c0 + c1*x + c2*x**2 + c3*x**3):"
    !!      do i = 1, 4
    !!          print '(AI0AF12.9)', "c", i - 1, " = ", p%get(i)
    !!      end do
    !!      print '(AE9.4)', "Residual: ", res
    !! end program
    !! @endcode
    !! The above program returns the following results.
    !! @code{.txt}
    !! Polynomial Coefficients (c0 + c1*x + c2*x**2 + c3*x**3):
    !! c0 =  1.186614186
    !! c1 =  0.446613631
    !! c2 = -0.122320499
    !! c3 =  1.064762822
    !! Residual: .5064E+00
    !! @endcode
    subroutine poly_fit(this, x, y, order, err)
        ! Arguments
        class(polynomial), intent(inout) :: this
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(inout), dimension(:) :: y
        integer(int32), intent(in) :: order
        class(errors), intent(inout), optional, target :: err

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
    !> @brief Fits a polynomial of the specified order that passes through zero
    !! to a data set.
    !!
    !! @param[in,out] this The polynomial object.
    !! @param[in] x An N-element array containing the independent variable data
    !!  points.  Notice, must be N > @p order.
    !! @param[in,out] y On input, an N-element array containing the dependent
    !!  variable data points.  On output, the contents are overwritten.
    !! @param[in] order The order of the polynomial (must be >= 1).
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - NL_INVALID_INPUT_ERROR: Occurs if a zero or negative polynomial order
    !!      was specified, or if order is too large for the data set.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
    !!  - NL_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are different sizes.
    subroutine poly_fit_thru_zero(this, x, y, order, err)
        ! Arguments
        class(polynomial), intent(inout) :: this
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(inout), dimension(:) :: y
        integer(int32), intent(in) :: order
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
    !> @brief Evaluates a polynomial at the specified points.
    !!
    !! @param[in] this The polynomial object.
    !! @param[in] x The value(s) at which to evaluate the polynomial.
    !!
    !! @return The value(s) of the polynomial at @p x.
    elemental function poly_eval_double(this, x) result(y)
        ! Arguments
        class(polynomial), intent(in) :: this
        real(real64), intent(in) :: x
        real(real64) :: y

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
    !> @brief Evaluates a polynomial at the specified points.
    !!
    !! @param[in] this The polynomial object.
    !! @param[in] x The value(s) at which to evaluate the polynomial.
    !!
    !! @return The value(s) of the polynomial at @p x.
    elemental function poly_eval_complex(this, x) result(y)
        ! Arguments
        class(polynomial), intent(in) :: this
        complex(real64), intent(in) :: x
        complex(real64) :: y

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
    !> @brief Returns the companion matrix for the polynomial.
    !!
    !! @param[in] this The polynomial object.
    !!
    !! @return The companion matrix.
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Companion_matrix)
    !! - [Wolfram MathWorld](http://mathworld.wolfram.com/CompanionMatrix.html)
    pure function poly_companion_mtx(this) result(c)
        ! Arguments
        class(polynomial), intent(in) :: this
        real(real64), dimension(this%order(), this%order()) :: c

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
    !> @brief Computes all the roots of a polynomial by computing the
    !! eigenvalues of the polynomial companion matrix.
    !!
    !! @param[in] this The polynomial object.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !!  - NL_OUT_OF_MEMORY_ERROR: Occurs if local memory must be allocated, and
    !!      there is insufficient memory available.
    !!  - NL_CONVERGENCE_ERROR: Occurs if the algorithm failed to converge.
    !!
    !! @par Usage
    !! The following code provides an example of how to compute the roots of a
    !! polynomial.  This examples uses a tenth order polynomial; however, this
    !! process is applicable to any order.
    !! @code{.f90}
    !! program example
    !!     use linalg_constants, only : dp, i32
    !!     use nonlin_polynomials
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: order = 10
    !!
    !!     ! Local Variables
    !!     integer(int32) :: i
    !!     type(polynomial) :: p
    !!     real(real64), dimension(order+1) :: coeff
    !!     complex(real64), allocatable, dimension(:) :: rts, sol
    !!
    !!     ! Define the polynomial
    !!     call random_number(coeff)
    !!     call p%initialize(order)
    !!     do i = 1, size(coeff)
    !!         call p%set(i, coeff(i))
    !!     end do
    !!
    !!     ! Compute the roots via the polynomial routine
    !!     rts = p%roots()
    !!
    !!     ! Compute the value of the polynomial at each root and ensure it
    !!     ! is sufficiently close to zero.
    !!     sol = p%evaluate(rts)
    !!     do i = 1, size(sol)
    !!         print '(AE9.3AE9.3A)', "(", real(sol(i)), ", ", aimag(sol(i)), ")"
    !!     end do
    !! end program
    !! @endcode
    !! The above program returns the following results.
    !! @code{.txt}
    !! (-.466E-14, -.161E-14)
    !! (-.466E-14, 0.161E-14)
    !! (-.999E-15, 0.211E-14)
    !! (-.999E-15, -.211E-14)
    !! (0.444E-15, 0.108E-14)
    !! (0.444E-15, -.108E-14)
    !! (-.144E-14, -.433E-14)
    !! (-.144E-14, 0.433E-14)
    !! (0.644E-14, -.100E-13)
    !! (0.644E-14, 0.100E-13)
    !! @endcode
    function poly_roots(this, err) result(z)
        ! Arguments
        class(polynomial), intent(in) :: this
        complex(real64), dimension(this%order()) :: z
        class(errors), intent(inout), optional, target :: err

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
    !> @brief Gets the requested polynomial coefficient by index.  The
    !! coefficient index is established as follows: c(1) + c(2) * x +
    !! c(3) * x**2 + ... c(n) * x**n-1.
    !!
    !! @param[in] this The polynomial.
    !! @param[in] ind The polynomial coefficient index (0 < ind <= order + 1).
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !! - NL_INVALID_INPUT_ERROR: Occurs if the requested index is less than or
    !!      equal to zero, or if the requested index exceeds the number of
    !!      polynomial coefficients.
    !!
    !! @return The requested coefficient.
    function get_poly_coefficient(this, ind, err) result(c)
        ! Arguments
        class(polynomial), intent(in) :: this
        integer(int32), intent(in) :: ind
        class(errors), intent(inout), optional, target :: err
        real(real64) :: c

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
    !> @brief Gets an array containing all the coefficients of the polynomial.
    !! The coefficient index is established as follows: c(1) + c(2) * x +
    !! c(3) * x**2 + ... c(n) * x**n-1.
    !!
    !! @param[in] this The polynomial object.
    !!
    !! @return The array of coefficients.
    pure function get_poly_coefficients(this) result(c)
        ! Arguments
        class(polynomial), intent(in) :: this
        real(real64), dimension(this%order() + 1) :: c

        ! Process
        if (this%order() == -1) return
        c = this%m_coeffs
    end function

! ------------------------------------------------------------------------------
    !> @brief Sets the requested polynomial coefficient by index.  The
    !! coefficient index is established as follows: c(1) + c(2) * x +
    !! c(3) * x**2 + ... c(n) * x**n-1.
    !!
    !! @param[in,out] this The polynomial.
    !! @param[in] ind The polynomial coefficient index (0 < ind <= order + 1).
    !! @param[in] c The polynomial coefficient.
    !! @param[out] err An optional errors-based object that if provided can be
    !!  used to retrieve information relating to any errors encountered during
    !!  execution.  If not provided, a default implementation of the errors
    !!  class is used internally to provide error handling.  Possible errors and
    !!  warning messages that may be encountered are as follows.
    !! - NL_INVALID_INPUT_ERROR: Occurs if the requested index is less than or
    !!      equal to zero, or if the requested index exceeds the number of
    !!      polynomial coefficients.
    subroutine set_poly_coefficient(this, ind, c, err)
        ! Arguments
        class(polynomial), intent(inout) :: this
        integer(int32), intent(in) :: ind
        real(real64), intent(in) :: c
        class(errors), intent(inout), optional, target :: err

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
    !> @brief Assigns the contents of one polynomial to another.
    !!
    !! @param[in,out] x The assignee.
    !! @param[in] y The polynomial to copy
    subroutine poly_equals(x, y)
        ! Arguments
        class(polynomial), intent(inout) :: x
        class(polynomial), intent(in) :: y

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
    !> @brief Assigns a number to each coefficient of the polynomial.
    !!
    !! @param[in,out] x The assignee.
    !! @param[in] y The value to assign.
    subroutine poly_dbl_equals(x, y)
        ! Arguments
        class(polynomial), intent(inout) :: x
        real(real64), intent(in) :: y

        ! Local Variables
        integer(int32) :: i, ord

        ! Process
        ord = x%order()
        do i = 1, ord + 1
            call x%set(i, y)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Assigns the contents of an array as polynomial coefficients.
    !!
    !! @param[in,out] x The assignee.
    !! @param[in] y The coefficient array.
    subroutine poly_equals_array(x, y)
        ! Arguments
        class(polynomial), intent(inout) :: x
        real(real64), intent(in), dimension(:) :: y
        call x%initialize(y)
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Adds two polynomials.
    !!
    !! @param[in] x The left-hand-side argument.
    !! @param[in] y The right-hand-side argument.
    !!
    !! @return The resulting polynomial.
    function poly_poly_add(x, y) result(z)
        ! Arguments
        class(polynomial), intent(in) :: x, y
        type(polynomial) :: z

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
    !> @brief Subtracts two polynomials.
    !!
    !! @param[in] x The left-hand-side argument.
    !! @param[in] y The right-hand-side argument.
    !!
    !! @return The resulting polynomial.
    function poly_poly_subtract(x, y) result(z)
        ! Arguments
        class(polynomial), intent(in) :: x, y
        type(polynomial) :: z

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
    !> @brief Multiplies two polynomials.
    !!
    !! @param[in] x The left-hand-side argument.
    !! @param[in] y The right-hand-side argument.
    !!
    !! @return The resulting polynomial.
    function poly_poly_mult(x, y) result(z)
        ! Arguments
        class(polynomial), intent(in) :: x, y
        type(polynomial) :: z

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
    !> @brief Multiplies a polynomial by a scalar value.
    !!
    !! @param[in] x The polynomial.
    !! @param[in] y The scalar value.
    !!
    !! @return The resulting polynomial.
    function poly_dbl_mult(x, y) result(z)
        ! Arguments
        class(polynomial), intent(in) :: x
        real(real64), intent(in) :: y
        type(polynomial) :: z

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
    !> @brief Multiplies a polynomial by a scalar value.
    !!
    !! @param[in] x The scalar value.
    !! @param[in] y The polynomial.
    !!
    !! @return The resulting polynomial.
    function dbl_poly_mult(x, y) result(z)
        ! Arguments
        real(real64), intent(in) :: x
        class(polynomial), intent(in) :: y
        type(polynomial) :: z

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
