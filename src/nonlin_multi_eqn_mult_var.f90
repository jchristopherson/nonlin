! MIT License
!
! Copyright (c) 2017-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
module nonlin_multi_eqn_mult_var
    use iso_fortran_env
    use nonlin_types
    use nonlin_error_handling
    use autodiff
    implicit none
    private
    public :: vecfcn
    public :: jacobianfcn
    public :: vecfcn_helper
    public :: equation_solver
    public :: nonlin_solver

    interface
        subroutine vecfcn(x, f, args)
            !! Describes an M-element vector-valued function of N-variables.
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
                !! An N-element array containing the independent variables.
            real(real64), intent(out), dimension(:) :: f
                !! An M-element array that, on output, contains the values
                !! of the M functions.
            class(*), intent(in), optional :: args
                !! An optional argument to allow the user to communicate with
                !! the routine.
        end subroutine

        subroutine jacobianfcn(x, jac, args)
            !! Describes a routine capable of computing the Jacobian matrix
            !! of M functions of N unknowns.
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
                !! An N-element array containing the independent variables.
            real(real64), intent(out), dimension(:,:) :: jac
                !! An M-by-N matrix where the Jacobian will be written.
            class(*), intent(in), optional :: args
                !! An optional argument to allow the user to communicate with
                !! the routine.
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
        procedure(vector_function), private, pointer, nopass :: m_autofcn => null()
            !! A pointer to an AD implementation of m_fcn that, if supplied, will be
            !! used in lieu of m_fcn, and will trigger an AD calculation of the gradient
            !! unless m_jac is supplied.  The pecking order is user-defined Jacobian
            !! as primary, AD Jacobian secondary, and the fallback is FD.
        procedure(jacobianfcn), private, pointer, nopass :: m_jac => null()
            !! A pointer to the jacobian routine - null if no routine is 
            !! supplied.
        integer, private :: m_nfcn = 0
            !! The number of functions in m_fcn.
        integer, private :: m_nvar = 0
            !! The number of variables in m_fcn.
        type(workspace), private :: m_workspace
            !! The AD workspace.
        integer, private :: m_threads = 1
            !! The number of threads to utilize when evaluating the Jacobian.
            !! If greater than 1, columns may be evaluated in parallel.
        real(real64), allocatable, dimension(:,:) :: m_jwork
            !! An FD workspace for the Jacobian calculations with one column
            !! available per thread.
        real(real64), allocatable, dimension(:,:) :: m_xjwork
            !! An FD workspace for the Jacobian calculations with one column
            !! available per thread.
        real(real64), allocatable, dimension(:) :: m_fwork
            !! A function evaluation workspace array.
    contains
        procedure, public :: set_fcn => vfh_set_fcn
        procedure, public :: set_auto_fcn => vfh_set_auto_fcn
        procedure, public :: set_jacobian => vfh_set_jac
        procedure, public :: is_fcn_defined => vfh_is_fcn_defined
        procedure, public :: is_jacobian_defined => vfh_is_jac_defined
        procedure, public :: fcn => vfh_fcn
        procedure, public :: jacobian => vfh_jac_fcn
        procedure, public :: get_equation_count => vfh_get_nfcn
        procedure, public :: get_variable_count => vfh_get_nvar
        procedure, public :: get_thread_count => vfh_get_nthreads
        procedure, public :: set_thread_count => vfh_set_nthreads
    end type

    type, abstract :: equation_solver
        !! A base class for various solvers of nonlinear systems of equations.
        integer, private :: m_maxEval = 100
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
        subroutine nonlin_solver(this, fcn, x, fvec, ib, args)
            !! Describes the interface of a nonlinear equation solver.
            use, intrinsic :: iso_fortran_env, only : real64
            use nonlin_types, only : iteration_behavior
            import equation_solver
            import vecfcn_helper
            class(equation_solver), intent(inout) :: this
                !! The [[equation_solver]]-based object.
            class(vecfcn_helper), intent(inout) :: fcn
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
            class(*), intent(in), optional :: args
                !! An optional argument to allow the user to communicate with
                !! the routine.
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
        integer, intent(in) :: nfcn
            !! The number of functions.
        integer, intent(in) :: nvar
            !! The number of variables.
        this%m_fcn => fcn
        this%m_nfcn = nfcn
        this%m_nvar = nvar
    end subroutine

! ------------------------------------------------------------------------------
    subroutine vfh_set_auto_fcn(this, fcn, nfcn, nvar)
        !! Establishes a pointer to the AD-compatible routine containing the
        !! function.  If utilized, this will allow for computation of the
        !! gradient, for applicable solvers, via AD.  However, a user-defined
        !! gradient does take precedence over an AD generated gradient.
        class(vecfcn_helper), intent(inout) :: this
            !! The [[vecfcn_helper]] object.
        procedure(vector_function), intent(in), pointer :: fcn
            !! The function pointer.
        integer, intent(in) :: nfcn
            !! The number of functions.
        integer, intent(in) :: nvar
            !! The number of variables.
        this%m_autofcn => fcn
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
        x = associated(this%m_fcn) .or. associated(this%m_autofcn)
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
    subroutine vfh_fcn(this, x, f, args)
        !! Executes the routine containing the system of equations to
        !! solve.  No action is taken if the pointer to the subroutine has not
        !! been defined.
        class(vecfcn_helper), intent(inout) :: this
            !! The [[vecfcn_helper]] object.
        real(real64), intent(in), dimension(:) :: x
            !! An N-element array containing the independent variables.
        real(real64), intent(out), dimension(:) :: f
            !! An M-element array that, on output, contains the values
            !! of the M functions.
        class(*), intent(in), optional :: args
            !! An optional argument to allow the user to communicate with
            !! the routine.
        if (associated(this%m_autofcn)) then
            if (.not.this%m_workspace%validate_size(size(x), size(f), &
                nthreads = this%get_thread_count())) &
            then
                call this%m_workspace%resize(size(x), size(f), &
                    nthreads = this%get_thread_count())
            end if
            this%m_workspace%xd(1:size(x))%x = x
            this%m_workspace%xd(1:size(x))%dx = 0.0d0
            call this%m_autofcn(this%m_workspace%xd(1:size(x)), &
                this%m_workspace%fd(1:size(f)), args = args)
            f = this%m_workspace%fd(1:size(f))%x
        else if (associated(this%m_fcn)) then
            call this%m_fcn(x, f, args)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine vfh_jac_fcn(this, x, jac, fv, args)
        !! Executes the routine containing the Jacobian matrix if
        !! supplied.  Otherwise, uses AD when an AD-compatible function is
        !! registered, falling back to finite differences.
!$ use omp_lib, only : omp_get_thread_num
        class(vecfcn_helper), intent(inout) :: this
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
        class(*), intent(in), optional :: args
            !! An optional argument to allow the user to communicate with
            !! the routine.

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer :: j, m, n, flag, nt, tid
        real(real64) :: eps, epsmch, h, temp

        ! Initialization
        m = this%get_equation_count()
        n = this%get_variable_count()
        nt = this%get_thread_count()

        ! Input Checking
        flag = 0
        if (size(x) /= n) then
            flag = 2
        else if (size(jac, 1) /= m .or. size(jac, 2) /= n) then
            flag = 3
        end if
        if (flag /= 0) error stop flag

        ! Process
        if (.not.this%is_fcn_defined()) error stop NL_UNDEFINED_FUNCTION_ERROR
        if (associated(this%m_jac)) then
            ! Call the user-defined Jacobian routine
            call this%m_jac(x, jac, args)
        else if (associated(this%m_autofcn)) then
            ! Compute the Jacobian via AD
            if (.not.this%m_workspace%validate_size(size(x), size(jac, 1), &
                nthreads = this%get_thread_count())) &
            then
                call this%m_workspace%resize(size(x), size(jac, 1), &
                    nthreads = this%get_thread_count())
            end if
            call jacobian(this%m_autofcn, x, jac, work = this%m_workspace, &
                args = args, nthreads = this%get_thread_count())
        else
            ! Compute the Jacobian via finite differences
            if (allocated(this%m_fwork)) then
                if (size(this%m_fwork) /= m) then
                    deallocate(this%m_fwork)
                    allocate(this%m_fwork(m))
                end if
            else
                allocate(this%m_fwork(m))
            end if
            
            if (allocated(this%m_jwork)) then
                if (size(this%m_jwork, 1) /= m .or. &
                    size(this%m_jwork, 2) /= nt) &
                then
                    deallocate(this%m_jwork)
                    allocate(this%m_jwork(m, nt))
                end if
            else
                allocate(this%m_jwork(m, nt))
            end if

            ! Allocate a workspace for x for each thread
            if (allocated(this%m_xjwork)) then
                if (size(this%m_xjwork, 1) /= n .or. &
                    size(this%m_xjwork, 2) /= nt) &
                then
                    deallocate(this%m_xjwork)
                    allocate(this%m_xjwork(n, nt))
                end if
            else
                allocate(this%m_xjwork(n, nt))
            end if
            
            ! Store the variables in the workspace
            do concurrent (j = 1:n)
                this%m_xjwork(j,:) = x(j)
            end do

            ! Define the function values at x
            if (present(fv)) then
                this%m_fwork = fv
            else
                call this%fcn(x, this%m_fwork, args = args)
            end if

            ! Establish step size factors
            epsmch = epsilon(epsmch)
            eps = sqrt(epsmch)

!$omp parallel num_threads(nt) default(shared) private(j, tid, h, temp)
            tid = 1
!$          tid = omp_get_thread_num() + 1

            ! Compute the Jacobian - each thread has it's own workspace column
!$omp do schedule(static)
            do j = 1, n
                temp = this%m_xjwork(j,tid)
                h = eps * abs(temp)
                if (h == zero) h = eps
                this%m_xjwork(j,tid) = temp + h
                call this%fcn(this%m_xjwork(:,tid), this%m_jwork(:,tid), args = args)
                this%m_xjwork(j,tid) = temp
                jac(:,j) = (this%m_jwork(:,tid) - this%m_fwork) / h
            end do
!$omp end do
!$omp end parallel
        end if
    end subroutine

! ------------------------------------------------------------------------------
    function vfh_get_nfcn(this) result(n)
        !! Gets the number of equations in this system.
        class(vecfcn_helper), intent(in) :: this
            !! The [[vecfcn_helper]] object.
        integer :: n
            !! The function count.
        n = this%m_nfcn
    end function

! ------------------------------------------------------------------------------
    function vfh_get_nvar(this) result(n)
        !! Gets the number of variables in this system.
        class(vecfcn_helper), intent(in) :: this
            !! The [[vecfcn_helper]] object.
        integer :: n
            !! The number of variables.
        n = this%m_nvar
    end function

! ------------------------------------------------------------------------------
    pure function vfh_get_nthreads(this) result(rst)
        !! Gets the maximum number of threads to use when evaluating the 
        !! Jacobian matrix.
        class(vecfcn_helper), intent(in) :: this
            !! The [[vecfcn_helper]] object.
        integer :: rst
            !! The number of threads allowed.
        rst = this%m_threads
    end function

! --------------------
    subroutine vfh_set_nthreads(this, n)
        !! Sets the maximum number of threads to use when evaluating the
        !! Jacobian matrix.
        !! With multiple threads, function callbacks must be thread-safe.
        !! Do not evaluate the same helper concurrently from separate callers.
        class(vecfcn_helper), intent(inout) :: this
            !! The [[vecfcn_helper]] object.
        integer, intent(in) :: n
            !! The maximum allowed number of threads.
        this%m_threads = max(1, n)
    end subroutine

! ******************************************************************************
! EQUATION_SOLVER
! ------------------------------------------------------------------------------
    pure function es_get_max_eval(this) result(n)
        !! Gets the maximum number of function evaluations allowed during
        !! a single solve.
        class(equation_solver), intent(in) :: this
            !! The [[equation_solver]] object.
        integer :: n
            !! The maximum number of function evaluations.
        n = this%m_maxEval
    end function

! --------------------
    subroutine es_set_max_eval(this, n)
        !! Sets the maximum number of function evaluations allowed during
        !! a single solve.
        class(equation_solver), intent(inout) :: this
            !! The [[equation_solver]] object.
        integer, intent(in) :: n
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