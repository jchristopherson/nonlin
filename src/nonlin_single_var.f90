module nonlin_single_var
    use iso_fortran_env
    use nonlin_types
    use ferror
    implicit none
    private
    public :: fcn1var
    public :: fcn1var_helper
    public :: equation_solver_1var
    public :: nonlin_solver_1var

    interface
        function fcn1var(x, args) result(f)
            !! Describes a function of one variable.
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in) :: x
                !! The independent variable.
            class(*), intent(inout), optional :: args
                !! An optional argument to allow the user to communicate with
                !! the routine.
            real(real64) :: f
                !! The value of the function at x.
        end function
    end interface


    type fcn1var_helper
        !! Defines a type capable of encapsulating an equation of one
        !! variable of the form: f(x) = 0.
        procedure(fcn1var), private, pointer, nopass :: m_fcn => null()
            !! A pointer to the target fcn1var routine.
        procedure(fcn1var), private, pointer, nopass :: m_diff => null()
            !! A pointer to a function capable of computing the derivative of 
            !! m_fcn.
    contains
        procedure, public :: fcn => f1h_fcn
        procedure, public :: is_fcn_defined => f1h_is_fcn_defined
        procedure, public :: set_fcn => f1h_set_fcn
        procedure, public :: is_derivative_defined => f1h_is_diff_defined
        procedure, public :: diff => f1h_diff_fcn
        procedure, public :: set_diff => f1h_set_diff
    end type

    type, abstract :: equation_solver_1var
        !! A base class for various solvers of equations of one variable.
        integer(int32), private :: m_maxEval = 100
            !! The maximum number of function evaluations allowed per solve.
        real(real64), private :: m_fcnTol = 1.0d-8
            !! The convergence criteria on function value.
        real(real64), private :: m_xtol = 1.0d-12
            !! The convergence criteria on change in variable value.
        real(real64), private :: m_difftol = 1.0d-12
            !! The convergence criteria on the slope of the function
            !! (derivative).
        logical, private :: m_printStatus = .false.
            !! Set to true to print iteration status; else, false.
    contains
        procedure, public :: get_max_fcn_evals => es1_get_max_eval
        procedure, public :: set_max_fcn_evals => es1_set_max_eval
        procedure, public :: get_fcn_tolerance => es1_get_fcn_tol
        procedure, public :: set_fcn_tolerance => es1_set_fcn_tol
        procedure, public :: get_var_tolerance => es1_get_var_tol
        procedure, public :: set_var_tolerance => es1_set_var_tol
        procedure, public :: get_print_status => es1_get_print_status
        procedure, public :: set_print_status => es1_set_print_status
        procedure(nonlin_solver_1var), deferred, public, pass :: solve
        procedure, public :: get_diff_tolerance => es1_get_diff_tol
        procedure, public :: set_diff_tolerance => es1_set_diff_tol
    end type

    interface
        subroutine nonlin_solver_1var(this, fcn, x, lim, f, ib, args, err)
            !! Describes the interface of a solver for an equation of one
            !! variable.
            use, intrinsic :: iso_fortran_env, only : real64
            use nonlin_types, only : iteration_behavior, value_pair
            use ferror, only : errors
            import equation_solver_1var
            import fcn1var_helper
            class(equation_solver_1var), intent(inout) :: this
                !! The [[equation_solver_1var]] object.
            class(fcn1var_helper), intent(in) :: fcn
                !! The fcn1var_helper object containing the equation to solve.
            real(real64), intent(inout) :: x
                !! On input the initial guess at the solution.  On output the 
                !! solution.
            type(value_pair), intent(in) :: lim
                !! A value_pair object defining the search limits.
            real(real64), intent(out), optional :: f
                !! An optional parameter used to return the function residual 
                !! as computed at x.
            type(iteration_behavior), optional :: ib
                !! An optional output, that if provided, allows the
                !! caller to obtain iteration performance information.
            class(*), intent(inout), optional :: args
                !! An optional argument to allow the user to communicate with
                !! the routine.
            class(errors), intent(inout), optional, target :: err
                !! An error handling object.
        end subroutine
    end interface

contains
! ******************************************************************************
! FCN1VAR_HELPER
! ------------------------------------------------------------------------------
    function f1h_fcn(this, x, args) result(f)
        !! Executes the routine containing the function to evaluate.
        class(fcn1var_helper), intent(in) :: this
            !! The [[fcn1var_helper]] object.
        real(real64), intent(in) :: x
            !! The value of the independent variable at which the function
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
    function f1h_is_fcn_defined(this) result(x)
        !! Tests if the pointer to the function containing the equation
        !! to solve has been assigned.
        class(fcn1var_helper), intent(in) :: this
            !! The [[fcn1var_helper]] object.
        logical :: x
            !! Returns true if the pointer has been assigned; else, false.
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    subroutine f1h_set_fcn(this, fcn)
        !! Establishes a pointer to the routine containing the equations
        !! to solve.
        class(fcn1var_helper), intent(inout) :: this
            !! The [[fcn1var_helper]] object.
        procedure(fcn1var), intent(in), pointer :: fcn
            !! The function pointer.
        this%m_fcn => fcn
    end subroutine

! ------------------------------------------------------------------------------
    function f1h_is_diff_defined(this) result(x)
        !! Tests if the pointer to the function containing the derivative of 
        !! the function to solve is defined.
        class(fcn1var_helper), intent(in) :: this
            !! The [[fcn1var_helper]] object.
        logical :: x
            !! Returns true if the pointer has been assigned; else, false.
        x = associated(this%m_diff)
    end function

! ------------------------------------------------------------------------------
    function f1h_diff_fcn(this, x, f, args) result(df)
        !! Computes the derivative of the function.  If a routine for computing 
        !! the derivative is not defined, the derivative is estimated via 
        !! finite differences.
        class(fcn1var_helper), intent(in) :: this
            !! The [[fcn1var_helper]] object.
        real(real64), intent(in) :: x
            !! The value of the independent variable at which the derivative is 
            !! to be computed.
        real(real64), intent(in), optional :: f
            !! An optional input specifying the function value at x.  If 
            !! supplied, and the derivative is being estimated numerically, the 
            !! function will not be evaluated at x.
        class(*), intent(inout), optional :: args
            !! An optional argument to allow the user to communicate with
            !! the routine.
        real(real64) :: df
            !! The value of the derivative.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        real(real64) :: eps, epsmch, h, temp, f1, f0

        ! Initialization
        epsmch = epsilon(epsmch)
        eps = sqrt(epsmch)

        ! Process
        if (this%is_derivative_defined()) then
            ! Use the user-defined routine to compute the derivative
            df = this%m_diff(x, args)
        else
            ! Compute the derivative via a forward difference
            h = eps * abs(x)
            if (h < epsmch) h = eps
            temp = x + h
            f1 = this%fcn(temp, args)
            if (present(f)) then
                f0 = f
            else
                f0 = this%fcn(x, args)
            end if
            df = (f1 - f0) / h
        end if
    end function

! ------------------------------------------------------------------------------
    subroutine f1h_set_diff(this, diff)
        !! Establishes a pointer to the routine containing the derivative of the
        !! equations to solve.
        class(fcn1var_helper), intent(inout) :: this
            !! The [[fcn1var_helper]] object.
        procedure(fcn1var), pointer, intent(in) :: diff
            !! A pointer to the function for computing the first derivative.
        this%m_diff => diff
    end subroutine

! ******************************************************************************
! EQUATION_SOLVER_1VAR
! ------------------------------------------------------------------------------
    pure function es1_get_max_eval(this) result(n)
        !! Gets the maximum number of function evaluations allowed during
        !! a single solve.
        class(equation_solver_1var), intent(in) :: this
            !! The [[equation_solver_1var]] object.
        integer(int32) :: n
            !! The maximum number of function evaluations.
        n = this%m_maxEval
    end function

! --------------------
    subroutine es1_set_max_eval(this, n)
        !! Sets the maximum number of function evaluations allowed during
        !! a single solve.
        class(equation_solver_1var), intent(inout) :: this
            !! The [[equation_solver_1var]] object.
        integer(int32), intent(in) :: n
            !! The maximum number of function evaluations.
        this%m_maxEval = n
    end subroutine

! ------------------------------------------------------------------------------
    pure function es1_get_fcn_tol(this) result(x)
        !! Gets the convergence on function value tolerance.
        class(equation_solver_1var), intent(in) :: this
            !! The [[equation_solver_1var]] object.
        real(real64) :: x
            !! The tolerance value.
        x = this%m_fcnTol
    end function

! --------------------
    subroutine es1_set_fcn_tol(this, x)
        !! Sets the convergence on function value tolerance.
        class(equation_solver_1var), intent(inout) :: this
            !! The [[equation_solver_1var]] object.
        real(real64), intent(in) :: x
            !! The tolerance value.
        this%m_fcnTol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure function es1_get_var_tol(this) result(x)
        !! Gets the convergence on change in variable tolerance.
        class(equation_solver_1var), intent(in) :: this
            !! The [[equation_solver_1var]] object.
        real(real64) :: x
            !! The tolerance value.
        x = this%m_xtol
    end function

! --------------------
    subroutine es1_set_var_tol(this, x)
        !! Sets the convergence on change in variable tolerance.
        class(equation_solver_1var), intent(inout) :: this
            !! The [[equation_solver_1var]] object.
        real(real64), intent(in) :: x
            !! The tolerance value.
        this%m_xtol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure function es1_get_print_status(this) result(x)
        !! Gets a logical value determining if iteration status should be
        !! printed.
        class(equation_solver_1var), intent(in) :: this
            !! The [[equation_solver_1var]] object.
        logical :: x
            !! True if the iteration status should be printed; else, false.
        x = this%m_printStatus
    end function

! --------------------
    subroutine es1_set_print_status(this, x)
        !! Sets a logical value determining if iteration status should be
        !! printed.
        class(equation_solver_1var), intent(inout) :: this
            !! The [[equation_solver_1var]] object.
        logical, intent(in) :: x
            !! True if the iteration status should be printed; else, false.
        this%m_printStatus = x
    end subroutine

! ------------------------------------------------------------------------------
    pure function es1_get_diff_tol(this) result(x)
        !! Gets the convergence on slope of the function (derivative)
        !! tolerance.
        class(equation_solver_1var), intent(in) :: this
            !! The [[equation_solver_1var]] object.
        real(real64) :: x
            !! The tolerance value.
        x = this%m_difftol
    end function

! --------------------
    subroutine es1_set_diff_tol(this, x)
        !! Sets the convergence on slope of the function (derivative)
        !! tolerance.
        class(equation_solver_1var), intent(inout) :: this
            !! The [[equation_solver_1var]] object.
        real(real64), intent(in) :: x
            !! The tolerance value.
        this%m_difftol = x
    end subroutine

! ------------------------------------------------------------------------------
end module