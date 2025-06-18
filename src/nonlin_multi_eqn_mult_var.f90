module nonlin_multi_eqn_mult_var
    use iso_fortran_env
    use nonlin_types
    use ferror
    implicit none
    private
    public :: vecfcn
    public :: jacobianfcn
    public :: vecfcn_helper
    public :: equation_solver
    public :: nonlin_solver

    interface
        subroutine vecfcn(x, f)
            !! Describes an M-element vector-valued function of N-variables.
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
                !! An N-element array containing the independent variables.
            real(real64), intent(out), dimension(:) :: f
                !! An M-element array that, on output, contains the values
                !! of the M functions.
        end subroutine

        subroutine jacobianfcn(x, jac)
            !! Describes a routine capable of computing the Jacobian matrix
            !! of M functions of N unknowns.
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
                !! An N-element array containing the independent variables.
            real(real64), intent(out), dimension(:,:) :: jac
                !! An M-by-N matrix where the Jacobian will be written.
        end subroutine
    end interface

    type vecfcn_helper
        !! Defines a type capable of encapsulating a system of nonlinear
        !! equations of the form: F(X) = 0.  This type is used to establish the
        !! system of equations to solve, and provides a means for computing
        !! the Jacobian matrix for the system of equations, and any other
        !! ancillary operations that may be needed by the solver.
        procedure(vecfcn), private, pointer, nopass :: m_fcn => null()
            !! A pointer to the target vecfcn routine.
        procedure(jacobianfcn), private, pointer, nopass :: m_jac => null()
            !! A pointer to the jacobian routine - null if no routine is 
            !! supplied.
        integer(int32), private :: m_nfcn = 0
            !! The number of functions in m_fcn.
        integer(int32) :: m_nvar = 0
            !! The number of variables in m_fcn.
    contains
        procedure, public :: set_fcn => vfh_set_fcn
        procedure, public :: set_jacobian => vfh_set_jac
        procedure, public :: is_fcn_defined => vfh_is_fcn_defined
        procedure, public :: is_jacobian_defined => vfh_is_jac_defined
        procedure, public :: fcn => vfh_fcn
        procedure, public :: jacobian => vfh_jac_fcn
        procedure, public :: get_equation_count => vfh_get_nfcn
        procedure, public :: get_variable_count => vfh_get_nvar
    end type

    type, abstract :: equation_solver
        !! A base class for various solvers of nonlinear systems of equations.
        integer(int32), private :: m_maxEval = 100
            !! The maximum number of function evaluations allowed per solve.
        real(real64), private :: m_fcnTol = 1.0d-8
            !! The convergence criteria on function values.
        real(real64), private :: m_xtol = 1.0d-12
            !! The convergence criteria on change in variable values.
        real(real64), private :: m_gtol = 1.0d-12
            !! The convergence criteria for the slope of the gradient vector.
        logical, private :: m_printStatus = .false.
            !! Set to true to print iteration status; else, false.
    contains
        procedure, public :: get_max_fcn_evals => es_get_max_eval
        procedure, public :: set_max_fcn_evals => es_set_max_eval
        procedure, public :: get_fcn_tolerance => es_get_fcn_tol
        procedure, public :: set_fcn_tolerance => es_set_fcn_tol
        procedure, public :: get_var_tolerance => es_get_var_tol
        procedure, public :: set_var_tolerance => es_set_var_tol
        procedure, public :: get_gradient_tolerance => es_get_grad_tol
        procedure, public :: set_gradient_tolerance => es_set_grad_tol
        procedure, public :: get_print_status => es_get_print_status
        procedure, public :: set_print_status => es_set_print_status
        procedure(nonlin_solver), deferred, public, pass :: solve
    end type

    interface
        subroutine nonlin_solver(this, fcn, x, fvec, ib, err)
            !! Describes the interface of a nonlinear equation solver.
            use, intrinsic :: iso_fortran_env, only : real64
            use nonlin_types, only : iteration_behavior
            use ferror, only : errors
            import equation_solver
            import vecfcn_helper
            class(equation_solver), intent(inout) :: this
                !! The [[equation_solver]]-based object.
            class(vecfcn_helper), intent(in) :: fcn
                !! The [[vecfcn_helper]] object containing the equations to
                !! solve.
            real(real64), intent(inout), dimension(:) :: x
                !! On input, an N-element array containing an initial estimate 
                !! to the solution.  On output, the updated solution estimate.
                !! N is the number of variables.
            real(real64), intent(out), dimension(:) :: fvec
                !! An M-element array that, on output, will contain the values 
                !! of each equation as evaluated at the variable values given 
                !! in x.
            type(iteration_behavior), optional :: ib
                !! An optional output, that if provided, allows the caller to 
                !! obtain iteration performance statistics.
            class(errors), intent(inout), optional, target :: err
                !! An error handling object.
        end subroutine
    end interface

contains
! ******************************************************************************
! VECFCN_HELPER
! ------------------------------------------------------------------------------
    subroutine vfh_set_fcn(this, fcn, nfcn, nvar)
        !! Establishes a pointer to the routine containing the system of
        !! equations to solve.
        class(vecfcn_helper), intent(inout) :: this
            !! The [[vecfcn_helper]] object.
        procedure(vecfcn), intent(in), pointer :: fcn
            !! The function pointer.
        integer(int32), intent(in) :: nfcn
            !! The number of functions.
        integer(int32), intent(in) :: nvar
            !! The number of variables.
        this%m_fcn => fcn
        this%m_nfcn = nfcn
        this%m_nvar = nvar
    end subroutine

! ------------------------------------------------------------------------------
    subroutine vfh_set_jac(this, jac)
        !! Establishes a pointer to the routine for computing the
        !! Jacobian matrix of the system of equations.  If no routine is
        !! defined, the Jacobian matrix will be computed numerically (this is
        !! the default state).
        class(vecfcn_helper), intent(inout) :: this
            !! The [[vecfcn_helper]] object.
        procedure(jacobianfcn), intent(in), pointer :: jac
            !! The function pointer.
        this%m_jac => jac
    end subroutine

! ------------------------------------------------------------------------------
    function vfh_is_fcn_defined(this) result(x)
        !! Tests if the pointer to the subroutine containing the system
        !! of equations to solve has been assigned.
        class(vecfcn_helper), intent(in) :: this
            !! The [[vecfcn_helper]] object.
        logical :: x
            !! Returns true if the pointer has been assigned; else, false.
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    function vfh_is_jac_defined(this) result(x)
        !! Tests if the pointer to the Jacobian calculation routine has been
        !! defined.
        class(vecfcn_helper), intent(in) :: this
            !! The [[vecfcn_helper]] object.
        logical :: x
            !! Returns true if the pointer has been assigned; else, false.
        x = associated(this%m_jac)
    end function

! ------------------------------------------------------------------------------
    subroutine vfh_fcn(this, x, f)
        !! Executes the routine containing the system of equations to
        !! solve.  No action is taken if the pointer to the subroutine has not
        !! been defined.
        class(vecfcn_helper), intent(in) :: this
            !! The [[vecfcn_helper]] object.
        real(real64), intent(in), dimension(:) :: x
            !! An N-element array containing the independent variables.
        real(real64), intent(out), dimension(:) :: f
            !! An M-element array that, on output, contains the values
            !! of the M functions.
        if (this%is_fcn_defined()) then
            call this%m_fcn(x, f)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine vfh_jac_fcn(this, x, jac, fv, work, olwork, err)
        !! Executes the routine containing the Jacobian matrix if
        !! supplied.  If not supplied, the Jacobian is computed via finite
        !! differences.
        class(vecfcn_helper), intent(in) :: this
            !! The [[vecfcn_helper]] object.
        real(real64), intent(inout), dimension(:) :: x
            !! An N-element array containing the independent variables defining 
            !! the point about which the derivatives will be calculated.
        real(real64), intent(out), dimension(:,:) :: jac
            !! An M-by-N matrix where, on output, the Jacobian will
            !! be written.
        real(real64), intent(in), dimension(:), optional, target :: fv
            !! An optional M-element array containing the function values at x. 
            !! If not supplied, the function values are computed at x.
        real(real64), intent(out), dimension(:), optional, target :: work
            !! An optional input, that if provided, prevents any local memory 
            !! allocation.  If not provided, the memory required is allocated
            !! within.  If provided, the length of the array must be at least
            !! olwork.  Notice, a workspace array is only utilized if the user 
            !! does not provide a routine for computing the Jacobian.
        integer(int32), intent(out), optional :: olwork
            !! An optional output used to determine workspace size. If supplied, 
            !! the routine determines the optimal size for work, and returns 
            !! without performing any actual calculations.
        integer(int32), intent(out), optional :: err
            !! An optional integer output that can be used to determine
            !! error status.  If not used, and an error is encountered, the 
            !! routine simply returns silently.  If used, the following error 
            !! codes identify error status:
            !!
            !!  - 0: No error has occurred.
            !!
            !!  - n: A positive integer denoting the index of an invalid input.
            !!
            !!  - -1: Indicates internal memory allocation failed.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: j, m, n, lwork, flag
        real(real64) :: eps, epsmch, h, temp
        real(real64), pointer, dimension(:) :: fptr, f1ptr
        real(real64), allocatable, target, dimension(:) :: wrk

        ! Initialization
        if (present(err)) err = 0
        ! m = this%m_nfcn
        ! n = this%m_nvar
        m = this%get_equation_count()
        n = this%get_variable_count()

        ! Input Checking
        flag = 0
        if (size(x) /= n) then
            flag = 2
        else if (size(jac, 1) /= m .or. size(jac, 2) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: Incorrectly sized input arrays
            if (present(err)) err = flag
            return
        end if

        ! Process
        if (.not.this%is_fcn_defined()) return
        if (associated(this%m_jac)) then
            ! Workspace Query
            if (present(olwork)) then
                olwork = 0
                return
            end if

            ! Call the user-defined Jacobian routine
            call this%m_jac(x, jac)
        else
            ! Compute the Jacobian via finite differences
            if (present(fv)) then
                lwork = m
            else
                lwork = 2 * m
            end if

            if (present(olwork)) then
                ! The user is just making a workspace query.  Simply return the
                ! workspace length, and exit the routine.
                olwork = lwork
                return
            end if

            ! Local Memory Allocation
            if (present(work)) then
                if (size(work) < lwork) then
                    ! ERROR: Workspace is too small
                    if (present(err)) err = 5
                    return
                end if
                f1ptr => work(1:m)
                if (present(fv)) then
                    if (size(fv) < m) then
                        ! ERROR: Function vector too small
                        if (present(err)) err = 4
                        return
                    end if
                    fptr => fv(1:m)
                else
                    fptr => work(m+1:2*m)
                    call this%fcn(x, fptr)
                end if
            else
                allocate(wrk(lwork), stat = flag)
                if (flag /= 0) then
                    ! ERROR: Memory issues
                    if (present(err)) err = -1
                    return
                end if
                f1ptr => wrk(1:m)
                if (present(fv)) then
                    fptr => fv(1:m)
                else
                    fptr => wrk(m+1:2*m)
                    call this%fcn(x, fptr)
                end if
            end if

            ! Establish step size factors
            epsmch = epsilon(epsmch)
            eps = sqrt(epsmch)

            ! Compute the derivatives via finite differences
            do j = 1, n
                temp = x(j)
                h = eps * abs(temp)
                if (h == zero) h = eps
                x(j) = temp + h
                call this%fcn(x, f1ptr)
                x(j) = temp
                jac(:,j) = (f1ptr - fptr) / h
            end do
        end if
    end subroutine

! ------------------------------------------------------------------------------
    function vfh_get_nfcn(this) result(n)
        !! Gets the number of equations in this system.
        class(vecfcn_helper), intent(in) :: this
            !! The [[vecfcn_helper]] object.
        integer(int32) :: n
            !! The function count.
        n = this%m_nfcn
    end function

! ------------------------------------------------------------------------------
    function vfh_get_nvar(this) result(n)
        !! Gets the number of variables in this system.
        class(vecfcn_helper), intent(in) :: this
            !! The [[vecfcn_helper]] object.
        integer(int32) :: n
            !! The number of variables.
        n = this%m_nvar
    end function

! ******************************************************************************
! EQUATION_SOLVER
! ------------------------------------------------------------------------------
    pure function es_get_max_eval(this) result(n)
        !! Gets the maximum number of function evaluations allowed during
        !! a single solve.
        class(equation_solver), intent(in) :: this
            !! The [[equation_solver]] object.
        integer(int32) :: n
            !! The maximum number of function evaluations.
        n = this%m_maxEval
    end function

! --------------------
    subroutine es_set_max_eval(this, n)
        !! Sets the maximum number of function evaluations allowed during
        !! a single solve.
        class(equation_solver), intent(inout) :: this
            !! The [[equation_solver]] object.
        integer(int32), intent(in) :: n
            !! The maximum number of function evaluations.
        this%m_maxEval = n
    end subroutine

! ------------------------------------------------------------------------------
    pure function es_get_fcn_tol(this) result(x)
        !! Gets the convergence on function value tolerance.
        class(equation_solver), intent(in) :: this
            !! The [[equation_solver]] object.
        real(real64) :: x
            !! The tolerance value.
        x = this%m_fcnTol
    end function

! --------------------
    subroutine es_set_fcn_tol(this, x)
        !! Sets the convergence on function value tolerance.
        class(equation_solver), intent(inout) :: this
            !! The [[equation_solver]] object.
        real(real64), intent(in) :: x
            !! The tolerance value.
        this%m_fcnTol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure function es_get_var_tol(this) result(x)
        !! Gets the convergence on change in variable tolerance.
        class(equation_solver), intent(in) :: this
            !! The [[equation_solver]] object.
        real(real64) :: x
            !! The tolerance value.
        x = this%m_xtol
    end function

! --------------------
    subroutine es_set_var_tol(this, x)
        !! Sets the convergence on change in variable tolerance.
        class(equation_solver), intent(inout) :: this
            !! The [[equation_solver]] object.
        real(real64), intent(in) :: x
            !! The tolerance value.
        this%m_xtol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure function es_get_grad_tol(this) result(x)
        !! Gets the convergence on slope of the gradient vector
        !! tolerance.
        class(equation_solver), intent(in) :: this
            !! The [[equation_solver]] object.
        real(real64) :: x
            !! The tolerance value.
        x = this%m_gtol
    end function

! --------------------
    subroutine es_set_grad_tol(this, x)
        !! Sets the convergence on slope of the gradient vector tolerance.
        class(equation_solver), intent(inout) :: this
            !! The [[equation_solver]] object.
        real(real64), intent(in) :: x
            !! The tolerance value.
        this%m_gtol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure function es_get_print_status(this) result(x)
        !! Gets a logical value determining if iteration status should be
        !! printed.
        class(equation_solver), intent(in) :: this
            !! The [[equation_solver]] object.
        logical :: x
            !! True if the iteration status should be printed; else, false.
        x = this%m_printStatus
    end function

! --------------------
    subroutine es_set_print_status(this, x)
        !! Sets a logical value determining if iteration status should be
        !! printed.
        class(equation_solver), intent(inout) :: this
            !! The [[equation_solver]] object.
        logical, intent(in) :: x
            !! True if the iteration status should be printed; else, false.
        this%m_printStatus = x
    end subroutine

! ------------------------------------------------------------------------------
end module