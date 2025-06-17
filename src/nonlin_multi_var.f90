module nonlin_multi_var
    use iso_fortran_env
    use nonlin_types
    use ferror
    implicit none
    private
    public :: fcnnvar
    public :: gradientfcn
    public :: fcnnvar_helper
    public :: equation_optimizer
    public :: nonlin_optimize_fcn

    interface
        function fcnnvar(x, args) result(f)
            !! Describes a function of N variables.
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
                !! An N-element array containing the independent variables.
            class(*), intent(inout), optional :: args
                !! An optional argument to allow the user to communicate with
                !! the routine.
            real(real64) :: f
                !! The value of the function.
        end function

        subroutine gradientfcn(x, g, args)
            !! Describes a routine capable of computing the gradient vector
            !! of an equation of N variables.
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
                !! An N-element array containing the independent variables.
            real(real64), intent(out), dimension(:) :: g
                !! An N-element array where the gradient vector will be
                !! written as output.
            class(*), intent(inout), optional :: args
                !! An optional argument to allow the user to communicate with
                !! the routine.
        end subroutine
    end interface

    type fcnnvar_helper
        !! Defines a type capable of encapsulating an equation of N variables.
        private
        procedure(fcnnvar), pointer, nopass :: m_fcn => null()
            !! A pointer to the target fcnnvar routine.
        procedure(gradientfcn), pointer, nopass :: m_grad => null()
            !! A pointer to the gradient routine.
        integer(int32) :: m_nvar = 0
            !! The number of variables in m_fcn.
    contains
        procedure, public :: fcn => fnh_fcn
        procedure, public :: is_fcn_defined => fnh_is_fcn_defined
        procedure, public :: set_fcn => fnh_set_fcn
        procedure, public :: get_variable_count => fnh_get_nvar
        procedure, public :: set_gradient_fcn => fnh_set_grad
        procedure, public :: is_gradient_defined => fnh_is_grad_defined
        procedure, public :: gradient => fnh_grad_fcn
    end type

    type, abstract :: equation_optimizer
        !! A base class for optimization of an equation of multiple variables.
        integer(int32), private :: m_maxEval = 500
            !! The maximum number of function evaluations allowed.
        real(real64), private :: m_tol = 1.0d-12
            !! The error tolerance used to determine convergence.
        logical, private :: m_printStatus = .false.
            !! Set to true to print iteration status; else, false.
    contains
        procedure, public :: get_max_fcn_evals => oe_get_max_eval
        procedure, public :: set_max_fcn_evals => oe_set_max_eval
        procedure, public :: get_tolerance => oe_get_tol
        procedure, public :: set_tolerance => oe_set_tol
        procedure, public :: get_print_status => oe_get_print_status
        procedure, public :: set_print_status => oe_set_print_status
        procedure(nonlin_optimize_fcn), deferred, public, pass :: solve
    end type

    interface
        subroutine nonlin_optimize_fcn(this, fcn, x, fout, ib, args, err)
            !! Describes the interface of a routine for optimizing an
            !! equation of N variables.
            use, intrinsic :: iso_fortran_env, only : real64
            use nonlin_types, only : iteration_behavior
            use ferror, only : errors
            import equation_optimizer
            import fcnnvar_helper
            class(equation_optimizer), intent(inout) :: this
                !! The [[equation_optimizer]] object.
            class(fcnnvar_helper), intent(in) :: fcn
                !! The [[fcnnvar_helper]] object containing the equation to
                !! optimize.
            real(real64), intent(inout), dimension(:) :: x
                !! On input, the initial guess at the optimal point.  On 
                !! output, the updated optimal point estimate.
            real(real64), intent(out), optional :: fout
                !! An optional output, that if provided, returns the value of 
                !! the function at x.
            type(iteration_behavior), optional :: ib
                !! An optional output, that if provided, allows the caller to 
                !! obtain iteration performance statistics.
            class(*), intent(inout), optional :: args
                !! An optional argument to allow the user to communicate with
                !! the routine.
            class(errors), intent(inout), optional, target :: err
                !! An error handling object.
        end subroutine
    end interface
contains
! ******************************************************************************
! FCNNVAR_HELPER
! ------------------------------------------------------------------------------
    function fnh_fcn(this, x, args) result(f)
        !! Executes the routine containing the function to evaluate.
        class(fcnnvar_helper), intent(in) :: this
            !! The [[fcnnvar_helper]] object.
        real(real64), intent(in), dimension(:) :: x
            !! The value of the independent variables at which the function
            !! should be evaluated.
        class(*), intent(inout), optional :: args
            !! An optional argument to allow the user to communicate with
            !! the routine.
        real(real64) :: f
            !! The value of the function.
        if (associated(this%m_fcn)) then
            f = this%m_fcn(x, args)
        end if
    end function

! ------------------------------------------------------------------------------
    function fnh_is_fcn_defined(this) result(x)
        !! Tests if the pointer to the function has been assigned.
        class(fcnnvar_helper), intent(in) :: this
            !! The [[fcnnvar_helper]] object.
        logical :: x
            !! Returns true if the pointer has been assigned; else, false.
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    subroutine fnh_set_fcn(this, fcn, nvar)
        !! Establishes a pointer to the routine containing the function.
        class(fcnnvar_helper), intent(inout) :: this
            !! The [[fcnnvar_helper]] object.
        procedure(fcnnvar), intent(in), pointer :: fcn
            !! The function pointer.
        integer(int32), intent(in) :: nvar
            !! The number of variables in the function.
        this%m_fcn => fcn
        this%m_nvar = nvar
    end subroutine

! ------------------------------------------------------------------------------
    function fnh_get_nvar(this) result(n)
        !! Gets the number of variables in this system.
        class(fcnnvar_helper), intent(in) :: this
            !! The [[fcnnvar_helper]] object.
        integer(int32) :: n
            !! The number of variables.
        n = this%m_nvar
    end function

! ------------------------------------------------------------------------------
    subroutine fnh_set_grad(this, fcn)
        !! Establishes a pointer to the routine containing the gradient
        !! vector of the function.
        class(fcnnvar_helper), intent(inout) :: this
            !! The [[fcnnvar_helper]] object.
        procedure(gradientfcn), pointer, intent(in) :: fcn
            !! The pointer to the gradient routine.
        this%m_grad => fcn
    end subroutine

! ------------------------------------------------------------------------------
    function fnh_is_grad_defined(this) result(x)
        !! Tests if the pointer to the routine containing the gradient
        !! has been assigned.
        class(fcnnvar_helper), intent(in) :: this
            !! The [[fcnnvar_helper]] object.
        logical :: x
            !! Returns true if the pointer has been assigned; else, false.
        x = associated(this%m_grad)
    end function

! ------------------------------------------------------------------------------
    subroutine fnh_grad_fcn(this, x, g, fv, args, err)
        !! Computes the gradient of the function.
        class(fcnnvar_helper), intent(in) :: this
            !! The [[fcnnvar_helper]] object.
        real(real64), intent(inout), dimension(:) :: x
            !! An N-element array containing the independent variables defining 
            !! the point about which the derivatives will be calculated.  This
            !! array is restored upon output.
        real(real64), intent(out), dimension(:) :: g
            !! An N-element array where the gradient will be written upon
            !! output.
        real(real64), intent(in), optional :: fv
            !! An optional input providing the function value at x.
        class(*), intent(inout), optional :: args
            !! An optional argument to allow the user to communicate with
            !! the routine.
        integer(int32), intent(out), optional :: err
            !! An optional integer output that can be used to determine error 
            !! status.  If not used, and an error is encountered, the routine
            !! simply returns silently.  If used, the following error codes 
            !! identify error status:
            !!
            !!  - 0: No error has occurred.
            !!
            !!  - n: A positive integer denoting the index of an invalid input.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: j, n, flag
        real(real64) :: eps, epsmch, h, temp, f, f1

        ! Initialization
        if (present(err)) err = 0
        ! n = this%m_nvar
        n = this%get_variable_count()

        ! Input Checking
        flag = 0
        if (size(x) /= n) then
            flag = 2
        else if (size(g) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: Incorrectly sized input arrays
            if (present(err)) err = flag
            return
        end if

        ! Process
        if (.not.this%is_fcn_defined()) return
        if (this%is_gradient_defined()) then
            ! Call the user-defined gradient routine
            call this%m_grad(x, g, args)
        else
            ! Compute the gradient via finite differences
            if (present(fv)) then
                f = fv
            else
                f = this%fcn(x, args)
            end if

            ! Establish step size factors
            epsmch = epsilon(epsmch)
            eps = sqrt(epsmch)

            ! Compute the derivatives
            do j = 1, n
                temp = x(j)
                h = eps * abs(temp)
                if (h == zero) h = eps
                x(j) = temp + h
                f1 = this%fcn(x, args)
                x(j) = temp
                g(j) = (f1 - f) / h
            end do
        end if
    end subroutine

! ******************************************************************************
! EQUATION_OPTIMIZER
! ------------------------------------------------------------------------------
    pure function oe_get_max_eval(this) result(n)
        class(equation_optimizer), intent(in) :: this
            !! The [[equation_optimizer]] object.
        integer(int32) :: n
            !! The maximum number of function evaluations.
        n = this%m_maxEval
    end function

! --------------------
    subroutine oe_set_max_eval(this, n)
        class(equation_optimizer), intent(inout) :: this
            !! The [[equation_optimizer]] object.
        integer(int32), intent(in) :: n
            !! The maximum number of function evaluations.
        this%m_maxEval = n
    end subroutine

! ------------------------------------------------------------------------------
    pure function oe_get_tol(this) result(x)
        class(equation_optimizer), intent(in) :: this
            !! The [[equation_optimizer]] object.
        real(real64) :: x
            !! The tolerance.
        x = this%m_tol
    end function

! --------------------
    subroutine oe_set_tol(this, x)
        class(equation_optimizer), intent(inout) :: this
            !! The [[equation_optimizer]] object.
        real(real64), intent(in) :: x
            !! The tolerance.
        this%m_tol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure function oe_get_print_status(this) result(x)
        class(equation_optimizer), intent(in) :: this
            !! The [[equation_optimizer]] object.
        logical :: x
            !! True if the iteration status should be printed; else, false.
        x = this%m_printStatus
    end function

! --------------------
    subroutine oe_set_print_status(this, x)
        class(equation_optimizer), intent(inout) :: this
            !! The [[equation_optimizer]] object.
        logical, intent(in) :: x
            !! True if the iteration status should be printed; else, false.
        this%m_printStatus = x
    end subroutine

! ------------------------------------------------------------------------------
end module