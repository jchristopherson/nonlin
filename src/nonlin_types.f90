! nonlin_types.f90

!> @mainpage
!!
!! @section intro_sec Introduction
!! NONLIN is a library that provides routines to compute the solutions to 
!! systems of nonlinear equations.
!!
!! @author Jason Christopherson
!! @version 1.1.5

!> @brief \b nonlin_types
!!
!! @par Purpose
!! To provide various types and constants useful in the solution of systems of
!! nonlinear equations.
module nonlin_types
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use, intrinsic :: iso_c_binding, only : c_bool
    use linalg_constants, only : LA_OUT_OF_MEMORY_ERROR, &
        LA_INVALID_OPERATION_ERROR, LA_CONVERGENCE_ERROR
    implicit none
    private
    public :: dp
    public :: i32
    public :: vecfcn
    public :: fcn1var
    public :: fcnnvar
    public :: jacobianfcn
    public :: gradientfcn
    public :: vecfcn_helper
    public :: fcn1var_helper
    public :: fcnnvar_helper
    public :: iteration_behavior
    public :: equation_solver
    public :: equation_solver_1var
    public :: equation_optimizer
    public :: value_pair
    public :: nonlin_solver
    public :: nonlin_solver_1var
    public :: nonlin_optimize
    public :: print_status
    public :: NL_INVALID_INPUT_ERROR
    public :: NL_ARRAY_SIZE_ERROR
    public :: NL_OUT_OF_MEMORY_ERROR
    public :: NL_INVALID_OPERATION_ERROR
    public :: NL_CONVERGENCE_ERROR
    public :: NL_DIVERGENT_BEHAVIOR_ERROR
    public :: NL_SPURIOUS_CONVERGENCE_ERROR
    public :: NL_TOLERANCE_TOO_SMALL_ERROR

! ******************************************************************************
! NUMERIC TYPE CONSTANTS
! ------------------------------------------------------------------------------
    !> @brief Defines a double-precision (64-bit) floating-point type.
    integer, parameter :: dp = real64
    !> @brief Defines a 32-bit signed integer type.
    integer, parameter :: i32 = int32

! ******************************************************************************
! ERROR FLAGS
! ------------------------------------------------------------------------------
    !> An error flag denoting an invalid input.
    integer, parameter :: NL_INVALID_INPUT_ERROR = 201
    !> An error flag denoting an improperly sized array.
    integer, parameter :: NL_ARRAY_SIZE_ERROR = 202
    !> An error denoting that there is insufficient memory available.
    integer, parameter :: NL_OUT_OF_MEMORY_ERROR = LA_OUT_OF_MEMORY_ERROR
    !> An error resulting from an invalid operation.
    integer, parameter :: NL_INVALID_OPERATION_ERROR = &
        LA_INVALID_OPERATION_ERROR
    !> An error resulting from a lack of convergence.
    integer, parameter :: NL_CONVERGENCE_ERROR = LA_CONVERGENCE_ERROR
    !> An error resulting from a divergent condition.
    integer, parameter :: NL_DIVERGENT_BEHAVIOR_ERROR = 206
    !> An error indicating a possible spurious convergence condition.
    integer, parameter :: NL_SPURIOUS_CONVERGENCE_ERROR = 207
    !> An error indicating the user-requested tolerance is too small to be
    !! practical for the problem at hand.
    integer, parameter :: NL_TOLERANCE_TOO_SMALL_ERROR = 208

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        !> @brief Describes a function of one variable.
        !!
        !! @param[in] x The independent variable.
        !! @return The value of the function at @p x.
        function fcn1var(x) result(f)
            use linalg_constants, only : dp
            real(dp), intent(in) :: x
            real(dp) :: f
        end function

        !> @brief Describes an M-element vector-valued function of N-variables.
        !!
        !! @param[in] x An N-element array containing the independent variables.
        !! @param[out] f An M-element array that, on output, contains the values
        !!  of the M functions.
        subroutine vecfcn(x, f)
            use linalg_constants, only : dp
            real(dp), intent(in), dimension(:) :: x
            real(dp), intent(out), dimension(:) :: f
        end subroutine

        !> @brief Describes a routine capable of computing the Jacobian matrix
        !! of M functions of N unknowns.
        !!
        !! @param[in] x An N-element array containing the independent variables.
        !! @param[out] jac An M-by-N matrix where the Jacobian will be written.
        subroutine jacobianfcn(x, jac)
            use linalg_constants, only : dp
            real(dp), intent(in), dimension(:) :: x
            real(dp), intent(out), dimension(:,:) :: jac
        end subroutine

        !> @brief Describes a function of N variables.
        !!
        !! @param[in] x An N-element array containing the independent variables.
        !! @return The value of the function at @p x.
        function fcnnvar(x) result(f)
            use linalg_constants, only : dp
            real(dp), intent(in), dimension(:) :: x
            real(dp) :: f
        end function

        !> @brief Describes a routine capable of computing the gradient vector
        !! of an equation of N variables.
        !!
        !! @param[in] x An N-element array containing the independent variables.
        !! @param[out] g An N-element array where the gradient vector will be
        !!  written as output.
        subroutine gradientfcn(x, g)
            use linalg_constants, only : dp
            real(dp), intent(in), dimension(:) :: x
            real(dp), intent(out), dimension(:) :: g
        end subroutine
    end interface

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief Defines a type capable of encapsulating a system of nonlinear 
    !! equations of the form: F(X) = 0.
    type vecfcn_helper
        private
        !> A pointer to the target vecfcn routine.
        procedure(vecfcn), pointer, nopass :: m_fcn => null()
        !> A pointer to the jacobian routine - null if no routine is supplied.
        procedure(jacobianfcn), pointer, nopass :: m_jac => null()
        !> The number of functions in m_fcn
        integer(i32) :: m_nfcn = 0
        !> The number of variables in m_fcn
        integer(i32) :: m_nvar = 0
    contains
        !> @brief Establishes a pointer to the routine containing the system of
        !!  equations to solve.
        procedure, public :: set_fcn => vfh_set_fcn
        !> @brief Establishes a pointer to the routine for computing the 
        !! Jacobian matrix of the system of equations.  If no routine is 
        !! defined, the Jacobian matrix will be computed numerically (this is 
        !! the default state).
        procedure, public :: set_jacobian => vfh_set_jac
        !> @brief Tests if the pointer to the subroutine containing the system
        !! of equations to solve has been assigned.
        procedure, public :: is_fcn_defined => vfh_is_fcn_defined
        !> @brief Tests if the pointer to the subroutine containing the system
        !! of equations to solve has been assigned.
        procedure, public :: is_jacobian_defined => vfh_is_jac_defined
        !> @brief Executes the routine containing the system of equations to
        !! solve.  No action is taken if the pointer to the subroutine has not 
        !! been defined.
        procedure, public :: fcn => vfh_fcn
        !> @brief Executes the routine containing the Jacobian matrix if 
        !! supplied.  If not supplied, the Jacobian is computed via finite 
        !! differences.
        procedure, public :: jacobian => vfh_jac_fcn
        !> @brief Gets the number of equations in this system.
        procedure, public :: get_equation_count => vfh_get_nfcn
        !> @brief Gets the number of variables in this system.
        procedure, public :: get_variable_count => vfh_get_nvar
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a type capable of encapsulating an equation of one 
    !! variable of the form: f(x) = 0.
    type fcn1var_helper
        private
        !> A pointer to the target fcn1var routine.
        procedure(fcn1var), pointer, nopass :: m_fcn => null()
    contains
        !> @brief Executes the routine containing the function to evaluate.
        procedure, public :: fcn => f1h_fcn
        !> @brief Tests if the pointer to the function containing the equation
        !! to solve has been assigned.
        procedure, public :: is_fcn_defined => f1h_is_fcn_defined
        !> @brief Establishes a pointer to the routine containing the equations
        !! to solve.
        procedure, public :: set_fcn => f1h_set_fcn
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a type capable of encapsulating an equation of N 
    !! variables.
    type fcnnvar_helper
        private
        !> A pointer to the target fcnnvar routine.
        procedure(fcnnvar), pointer, nopass :: m_fcn => null()
        !> A pointer to the gradient routine.
        procedure(gradientfcn), pointer, nopass :: m_grad => null()
        !> The number of variables in m_fcn
        integer(i32) :: m_nvar = 0
    contains
        !> @brief Executes the routine containing the function to evaluate.
        procedure, public :: fcn => fnh_fcn
        !> @brief Tests if the pointer to the function has been assigned.
        procedure, public :: is_fcn_defined => fnh_is_fcn_defined
        !> @brief Establishes a pointer to the routine containing the function.
        procedure, public :: set_fcn => fnh_set_fcn
        !> @brief Gets the number of variables in this system.
        procedure, public :: get_variable_count => fnh_get_nvar
        !> @brief Establishes a pointer to the routine containing the gradient
        !! vector of the function.
        procedure, public :: set_gradient_fcn => fnh_set_grad
        !> @brief Tests if the pointer to the routine containing the gradient 
        !! has been assigned.
        procedure, public :: is_gradient_defined => fnh_is_grad_defined
        !> @brief Computes the gradient of the function.
        procedure, public :: gradient => fnh_grad_fcn
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a set of parameters that describe the behavior of the
    !! iteration process.
    type, bind(c) :: iteration_behavior
        !> Specifies the number of iterations performed.
        integer(i32) :: iter_count
        !> Specifies the number of function evaluations performed.
        integer(i32) :: fcn_count
        !> Specifies the number of Jacobian evaluations performed.
        integer(i32) :: jacobian_count
        !> Specifies the number of gradient vector evaluations performed.
        integer(i32) :: gradient_count
        !> True if the solution converged as a result of a zero-valued
        !! function; else, false.
        logical(c_bool) :: converge_on_fcn
        !> True if the solution converged as a result of no appreciable
        !! change in solution points between iterations; else, false.
        logical(c_bool) :: converge_on_chng
        !> True if the solution appears to have settled on a stationary
        !! point such that the gradient of the function is zero-valued; else,
        !! false.
        logical(c_bool) :: converge_on_zero_diff
    end type

! ------------------------------------------------------------------------------
    !> @brief A base class for various solvers of nonlinear systems of 
    !! equations.
    type, abstract :: equation_solver
        private
        !> The maximum number of function evaluations allowed per solve.
        integer(i32) :: m_maxEval = 100
        !> The convergence criteria on function values.
        real(dp) :: m_fcnTol = 1.0d-8
        !> The convergence criteria on change in variable values.
        real(dp) :: m_xtol = 1.0d-12
        !> The convergence criteria for the slope of the gradient vector.
        real(dp) :: m_gtol = 1.0d-12
        !> Set to true to print iteration status; else, false.
        logical :: m_printStatus = .false.
    contains
        !> @brief Gets the maximum number of function evaluations allowed during
        !! a single solve.
        procedure, public :: get_max_fcn_evals => es_get_max_eval
        !> @brief Sets the maximum number of function evaluations allowed during
        !! a single solve.
        procedure, public :: set_max_fcn_evals => es_set_max_eval
        !> @brief Gets the convergence on function value tolerance.
        procedure, public :: get_fcn_tolerance => es_get_fcn_tol
        !> @brief Sets the convergence on function value tolerance.
        procedure, public :: set_fcn_tolerance => es_set_fcn_tol
        !> @brief Gets the convergence on change in variable tolerance.
        procedure, public :: get_var_tolerance => es_get_var_tol
        !> @brief Sets the convergence on change in variable tolerance.
        procedure, public :: set_var_tolerance => es_set_var_tol
        !> @brief Gets the convergence on slope of the gradient vector
        !! tolerance.
        procedure, public :: get_gradient_tolerance => es_get_grad_tol
        !> @brief Sets the convergence on slope of the gradient vector
        !! tolerance.
        procedure, public :: set_gradient_tolerance => es_set_grad_tol
        !> @brief Gets a logical value determining if iteration status should be
        !! printed.
        procedure, public :: get_print_status => es_get_print_status
        !> @brief Sets a logical value determining if iteration status should be
        !! printed.
        procedure, public :: set_print_status => es_set_print_status
        !> @brief Solves the system of equations.
        procedure(nonlin_solver), deferred, public, pass :: solve
    end type

! ------------------------------------------------------------------------------
    !> @brief A base class for various solvers of equations of one variable.
    type, abstract :: equation_solver_1var
        private
        !> The maximum number of function evaluations allowed per solve.
        integer(i32) :: m_maxEval = 100
        !> The convergence criteria on function value.
        real(dp) :: m_fcnTol = 1.0d-8
        !> The convergence criteria on change in variable value.
        real(dp) :: m_xtol = 1.0d-12
        !> Set to true to print iteration status; else, false.
        logical :: m_printStatus = .false.
    contains
        !> @brief Gets the maximum number of function evaluations allowed during
        !! a single solve.
        procedure, public :: get_max_fcn_evals => es1_get_max_eval
        !> @brief Sets the maximum number of function evaluations allowed during
        !! a single solve.
        procedure, public :: set_max_fcn_evals => es1_set_max_eval
        !> @brief Gets the convergence on function value tolerance.
        procedure, public :: get_fcn_tolerance => es1_get_fcn_tol
        !> @brief Sets the convergence on function value tolerance.
        procedure, public :: set_fcn_tolerance => es1_set_fcn_tol
        !> @brief Gets the convergence on change in variable tolerance.
        procedure, public :: get_var_tolerance => es1_get_var_tol
        !> @brief Sets the convergence on change in variable tolerance.
        procedure, public :: set_var_tolerance => es1_set_var_tol
        !> @brief Gets a logical value determining if iteration status should be
        !! printed.
        procedure, public :: get_print_status => es1_get_print_status
        !> @brief Sets a logical value determining if iteration status should be
        !! printed.
        procedure, public :: set_print_status => es1_set_print_status
        !> @brief Solves the equation.
        procedure(nonlin_solver_1var), deferred, public, pass :: solve
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a pair of numeric values.
    type, bind(c) :: value_pair
        !> Value 1.
        real(dp) :: x1
        !> Value 2.
        real(dp) :: x2
    end type

! ------------------------------------------------------------------------------
    !> @brief A base class for optimization of an equation of multiple 
    !! variables.
    type, abstract :: equation_optimizer
        private
        !> The maximum number of function evaluations allowed.
        integer(i32) :: m_maxEval = 500
        !> The error tolerance used to determine convergence.
        real(dp) :: m_tol = 1.0d-12
        !> Set to true to print iteration status; else, false.
        logical :: m_printStatus = .false.
    contains
        !> @brief Gets the maximum number of function evaluations allowed.
        procedure, public :: get_max_fcn_evals => oe_get_max_eval
        !> @brief Sets the maximum number of function evaluations allowed.
        procedure, public :: set_max_fcn_evals => oe_set_max_eval
        !> @brief Gets the tolerance on convergence.
        procedure, public :: get_tolerance => oe_get_tol
        !> @brief Sets the tolerance on convergence.
        procedure, public :: set_tolerance => oe_set_tol
        !> @brief Gets a logical value determining if iteration status should be
        !! printed.
        procedure, public :: get_print_status => oe_get_print_status
        !> @brief Sets a logical value determining if iteration status should be
        !! printed.
        procedure, public :: set_print_status => oe_set_print_status
        !> @brief Optimizes the equation.
        procedure(nonlin_optimize), deferred, public, pass :: solve
    end type

! ******************************************************************************
! ABSTRACT ROUTINE INTERFACES
! ------------------------------------------------------------------------------
    interface
        !> @brief Describes the interface of a nonlinear equation solver.
        !!
        !! @param[in,out] this The equation_solver-based object.
        !! @param[in] fcn The vecfcn_helper object containing the equations to
        !!  solve.
        !! @param[in,out] x On input, an N-element array containing an initial
        !!  estimate to the solution.  On output, the updated solution estimate.
        !!  N is the number of variables.
        !! @param[out] fvec An M-element array that, on output, will contain
        !!  the values of each equation as evaluated at the variable values
        !!  given in @p x.
        !! @param[out] ib An optional output, that if provided, allows the
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can
        !!  be used to retrieve information relating to any errors encountered
        !!  during execution.  If not provided, a default implementation of the
        !!  errors class is used internally to provide error handling.  The
        !!  possible error codes returned will likely vary from solver to
        !!  solver.
        subroutine nonlin_solver(this, fcn, x, fvec, ib, err)
            use linalg_constants, only : dp, i32
            use ferror, only : errors
            import equation_solver
            import vecfcn_helper
            import iteration_behavior
            class(equation_solver), intent(inout) :: this
            class(vecfcn_helper), intent(in) :: fcn
            real(dp), intent(inout), dimension(:) :: x
            real(dp), intent(out), dimension(:) :: fvec
            type(iteration_behavior), optional :: ib
            class(errors), intent(inout), optional, target :: err
        end subroutine

        !> @brief Describes the interface of a solver for an equation of one
        !! variable.
        !!
        !! @param[in,out] this The equation_solver_1var-based object.
        !! @param[in] fcn The fcn1var_helper object containing the equation
        !!  to solve.
        !! @param[in,out] x On input the initial guess at the solution.  On
        !!  output the updated solution estimate.
        !! @param[in] lim A value_pair object defining the search limits.
        !! @param[out] f An optional parameter used to return the function
        !!  residual as computed at @p x.
        !! @param[out] ib An optional output, that if provided, allows the
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can
        !!  be used to retrieve information relating to any errors encountered
        !!  during execution.  If not provided, a default implementation of the
        !!  errors class is used internally to provide error handling.  The
        !!  possible error codes returned will likely vary from solver to
        !!  solver.
        subroutine nonlin_solver_1var(this, fcn, x, lim, f, ib, err)
            use linalg_constants, only : dp, i32
            use ferror, only : errors
            import equation_solver_1var
            import fcn1var_helper
            import value_pair
            import iteration_behavior
            class(equation_solver_1var), intent(inout) :: this
            class(fcn1var_helper), intent(in) :: fcn
            real(dp), intent(inout) :: x
            type(value_pair), intent(in) :: lim
            real(dp), intent(out), optional :: f
            type(iteration_behavior), optional :: ib
            class(errors), intent(inout), optional, target :: err
        end subroutine

        !> @brief Describes the interface of a routine for optimizing an 
        !! equation of N variables.
        !!
        !! @param[in,out] this The equation_optimizer-based object.
        !! @param[in] fcn The fcnnvar_helper object containing the equation to
        !!  optimize.
        !! @param[in,out] x On input, the initial guess at the optimal point. 
        !!  On output, the updated optimal point estimate.
        !! @param[out] fout An optional output, that if provided, returns the
        !!  value of the function at @p x.
        !! @param[out] ib An optional output, that if provided, allows the
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can
        !!  be used to retrieve information relating to any errors encountered
        !!  during execution.  If not provided, a default implementation of the
        !!  errors class is used internally to provide error handling.  The
        !!  possible error codes returned will likely vary from solver to
        !!  solver.
        subroutine nonlin_optimize(this, fcn, x, fout, ib, err)
            use linalg_constants, only : dp, i32
            use ferror, only : errors
            import equation_optimizer
            import fcnnvar_helper
            import iteration_behavior
            class(equation_optimizer), intent(inout) :: this
            class(fcnnvar_helper), intent(in) :: fcn
            real(dp), intent(inout), dimension(:) :: x
            real(dp), intent(out), optional :: fout
            type(iteration_behavior), optional :: ib
            class(errors), intent(inout), optional, target :: err
        end subroutine

    end interface


contains
! ******************************************************************************
! VECFCN_HELPER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine containing the system of
    !!  equations to solve.
    !! 
    !! @param[in,out] this The vecfcn_helper object.
    !! @param[in] fcn The function pointer.
    !! @param[in] nfcn The number of functions.
    !! @param[in] nvar The number of variables.
    subroutine vfh_set_fcn(this, fcn, nfcn, nvar)
        class(vecfcn_helper), intent(inout) :: this
        procedure(vecfcn), intent(in), pointer :: fcn
        integer(i32), intent(in) :: nfcn, nvar
        this%m_fcn => fcn
        this%m_nfcn = nfcn
        this%m_nvar = nvar
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine for computing the Jacobian
    !! matrix of the system of equations.  If no routine is defined, the
    !! Jacobian matrix will be computed numerically (this is the default state).
    !!
    !! @param[in,out] this The vecfcn_helper object.
    !! @param[in] jac The function pointer.
    subroutine vfh_set_jac(this, jac)
        class(vecfcn_helper), intent(inout) :: this
        procedure(jacobianfcn), intent(in), pointer :: jac
        this%m_jac => jac
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the subroutine containing the system of
    !! equations to solve has been assigned.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function vfh_is_fcn_defined(this) result(x)
        class(vecfcn_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the subroutine containing the system of 
    !! equations to solve has been assigned.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function vfh_is_jac_defined(this) result(x)
        class(vecfcn_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_jac)
    end function

! ------------------------------------------------------------------------------
    !> @brief Executes the routine containing the system of equations to solve.
    !! No action is taken if the pointer to the subroutine has not been defined.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @param[in] x An N-element array containing the independent variables.
    !! @param[out] f An M-element array that, on output, contains the values
    !!  of the M functions.
    subroutine vfh_fcn(this, x, f)
        class(vecfcn_helper), intent(in) :: this
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: f
        if (this%is_fcn_defined()) then
            call this%m_fcn(x, f)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Executes the routine containing the Jacobian matrix if supplied.
    !! If not supplied, the Jacobian is computed via finite differences.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @param[in] x An N-element array containing the independent variabls 
    !!  defining the point about which the derivatives will be calculated.
    !! @param[out] jac An M-by-N matrix where, on output, the Jacobian will
    !!  be written.
    !! @param[in] fv An optional M-element array containing the function values
    !!  at @p x.  If not supplied, the function values are computed at @p x.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.  Notice, a workspace array is only utilized if the user does
    !!  not provide a routine for computing the Jacobian.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional integer output that can be used to determine
    !!  error status.  If not used, and an error is encountered, the routine
    !!  simply returns silently.  If used, the following error codes identify
    !!  error status:
    !!  - 0: No error has occurred.
    !!  - n: A positive integer denoting the index of an invalid input.
    !!  - -1: Indicates internal memory allocation failed.
    subroutine vfh_jac_fcn(this, x, jac, fv, work, olwork, err)
        ! Arguments
        class(vecfcn_helper), intent(in) :: this
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), dimension(:,:) :: jac
        real(dp), intent(in), dimension(:), optional, target :: fv
        real(dp), intent(out), dimension(:), optional, target :: work
        integer(i32), intent(out), optional :: olwork, err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: j, m, n, lwork, flag
        real(dp) :: eps, epsmch, h, temp
        real(dp), pointer, dimension(:) :: fptr, f1ptr
        real(dp), allocatable, target, dimension(:) :: wrk

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
    !> @brief Gets the number of equations in this system.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @return The function count.
    pure function vfh_get_nfcn(this) result(n)
        class(vecfcn_helper), intent(in) :: this
        integer(i32) :: n
        n = this%m_nfcn
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets the number of variables in this system.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @return The number of variables.
    pure function vfh_get_nvar(this) result(n)
        class(vecfcn_helper), intent(in) :: this
        integer(i32) :: n
        n = this%m_nvar
    end function

! ******************************************************************************
! FCN1VAR_HELPER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Executes the routine containing the function to evaluate.
    !!
    !! @param[in] this The fcn1var_helper object.
    !! @param[in] x The value of the independent variable at which the function
    !!  should be evaluated.
    !! @return The value of the function at @p x.
    function f1h_fcn(this, x) result(f)
        class(fcn1var_helper), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: f
        if (associated(this%m_fcn)) then
            f = this%m_fcn(x)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the function containing the equation to 
    !! solve has been assigned.
    !!
    !! @param[in] this The fcn1var_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function f1h_is_fcn_defined(this) result(x)
        class(fcn1var_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine containing the equations to 
    !! solve.
    !! 
    !! @param[in,out] this The fcn1var_helper object.
    !! @param[in] fcn The function pointer.
    subroutine f1h_set_fcn(this, fcn)
        class(fcn1var_helper), intent(inout) :: this
        procedure(fcn1var), intent(in), pointer :: fcn
        this%m_fcn => fcn
    end subroutine

! ******************************************************************************
! FCNNVAR_HELPER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Executes the routine containing the function to evaluate.
    !!
    !! @param[in] this The fcnnvar_helper object.
    !! @param[in] x The value of the independent variable at which the function
    !!  should be evaluated.
    !! @return The value of the function at @p x.
    function fnh_fcn(this, x) result(f)
        class(fcnnvar_helper), intent(in) :: this
        real(dp), intent(in), dimension(:) :: x
        real(dp) :: f
        if (associated(this%m_fcn)) then
            f = this%m_fcn(x)
        end if
    end function

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the function containing the equation to 
    !! solve has been assigned.
    !!
    !! @param[in] this The fcnnvar_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function fnh_is_fcn_defined(this) result(x)
        class(fcnnvar_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine containing the equations to 
    !! solve.
    !! 
    !! @param[in,out] this The fcnnvar_helper object.
    !! @param[in] fcn The function pointer.
    !! @param[in] nvar The number of variables in the function.
    subroutine fnh_set_fcn(this, fcn, nvar)
        class(fcnnvar_helper), intent(inout) :: this
        procedure(fcnnvar), intent(in), pointer :: fcn
        integer(i32), intent(in) :: nvar
        this%m_fcn => fcn
        this%m_nvar = nvar
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the number of variables in this system.
    !!
    !! @param[in] this The fcnnvar_helper object.
    !! @return The number of variables.
    pure function fnh_get_nvar(this) result(n)
        class(fcnnvar_helper), intent(in) :: this
        integer(i32) :: n
        n = this%m_nvar
    end function

! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine containing the gradient
    !! vector of the function.
    !!
    !! @param[in,out] this The fcnnvar_helper object.
    !! @param[in] fcn The pointer to the gradient routine.
    subroutine fnh_set_grad(this, fcn)
        class(fcnnvar_helper), intent(inout) :: this
        procedure(gradientfcn), pointer, intent(in) :: fcn
        this%m_grad => fcn
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the routine containing the gradient has
    !! been assigned.
    !!
    !! @param[in] this The fcnnvar_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function fnh_is_grad_defined(this) result(x)
        class(fcnnvar_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_grad)
    end function

! ------------------------------------------------------------------------------
    !> @brief Executes the routine containing the gradient, if supplied.  If not
    !! supplied, the gradient is computed via finite differences.
    !!
    !! @param[in] this The fcnnvar_helper object.
    !! @param[in,out] x An N-element array containing the independent variables
    !!  defining the point about which the derivatives will be calculated.  This
    !!  array is restored upon output.
    !! @param[out] g An N-element array where the gradient will be written upon
    !!  output.
    !! @param[in] fv An optional input providing the function value at @p x.
    !! @param[out] err An optional integer output that can be used to determine
    !!  error status.  If not used, and an error is encountered, the routine
    !!  simply returns silently.  If used, the following error codes identify
    !!  error status:
    !!  - 0: No error has occurred.
    !!  - n: A positive integer denoting the index of an invalid input.
    subroutine fnh_grad_fcn(this, x, g, fv, err)
        ! Arguments
        class(fcnnvar_helper), intent(in) :: this
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: g
        real(dp), intent(in), optional :: fv
        integer(i32), intent(out), optional :: err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: j, n, flag
        real(dp) :: eps, epsmch, h, temp, f, f1

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
            call this%m_grad(x, g)
        else
            ! Compute the gradient via finite differences
            if (present(fv)) then
                f = fv
            else
                f = this%fcn(x)
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
                f1 = this%fcn(x)
                x(j) = temp
                g(j) = (f1 - f) / h
            end do
        end if
    end subroutine

! ******************************************************************************
! EQUATION_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Gets the maximum number of function evaluations allowed during a
    !! single solve.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The maximum number of function evaluations.
    pure function es_get_max_eval(this) result(n)
        class(equation_solver), intent(in) :: this
        integer(i32) :: n
        n = this%m_maxEval
    end function

! --------------------
    !> @brief Sets the maximum number of function evaluations allowed during a
    !! single solve.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] n The maximum number of function evaluations.
    subroutine es_set_max_eval(this, n)
        class(equation_solver), intent(inout) :: this
        integer(i32), intent(in) :: n
        this%m_maxEval = n
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on function value tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    pure function es_get_fcn_tol(this) result(x)
        class(equation_solver), intent(in) :: this
        real(dp) :: x
        x = this%m_fcnTol
    end function

! --------------------
    !> @brief Sets the convergence on function value tolerance.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] x The tolerance value.
    subroutine es_set_fcn_tol(this, x)
        class(equation_solver), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_fcnTol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on change in variable tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    pure function es_get_var_tol(this) result(x)
        class(equation_solver), intent(in) :: this
        real(dp) :: x
        x = this%m_xtol
    end function

! --------------------
    !> @brief Sets the convergence on change in variable tolerance.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] x The tolerance value.
    subroutine es_set_var_tol(this, x)
        class(equation_solver), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_xtol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on slope of the gradient vector tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    pure function es_get_grad_tol(this) result(x)
        class(equation_solver), intent(in) :: this
        real(dp) :: x
        x = this%m_gtol
    end function

! --------------------
    !> @brief Sets the convergence on slope of the gradient vector tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    subroutine es_set_grad_tol(this, x)
        class(equation_solver), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_gtol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets a logical value determining if iteration status should be
    !! printed.
    !!
    !! @param[in] this The equation_solver object.
    !! @return True if the iteration status should be printed; else, false.
    pure function es_get_print_status(this) result(x)
        class(equation_solver), intent(in) :: this
        logical :: x
        x = this%m_printStatus
    end function

! --------------------
    !> @brief Sets a logical value determining if iteration status should be
    !! printed.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] x True if the iteration status should be printed; else, false.
    subroutine es_set_print_status(this, x)
        class(equation_solver), intent(inout) :: this
        logical, intent(in) :: x
        this%m_printStatus = x
    end subroutine

! ******************************************************************************
! EQUATION_SOLVER_1VAR MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Gets the maximum number of function evaluations allowed during a
    !! single solve.
    !!
    !! @param[in] this The equation_solver_1var object.
    !! @return The maximum number of function evaluations.
    pure function es1_get_max_eval(this) result(n)
        class(equation_solver_1var), intent(in) :: this
        integer(i32) :: n
        n = this%m_maxEval
    end function

! --------------------
    !> @brief Sets the maximum number of function evaluations allowed during a
    !! single solve.
    !!
    !! @param[in,out] this The equation_solver_1var object.
    !! @param[in] n The maximum number of function evaluations.
    subroutine es1_set_max_eval(this, n)
        class(equation_solver_1var), intent(inout) :: this
        integer(i32), intent(in) :: n
        this%m_maxEval = n
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on function value tolerance.
    !!
    !! @param[in] this The equation_solver_1var object.
    !! @return The tolerance value.
    pure function es1_get_fcn_tol(this) result(x)
        class(equation_solver_1var), intent(in) :: this
        real(dp) :: x
        x = this%m_fcnTol
    end function

! --------------------
    !> @brief Sets the convergence on function value tolerance.
    !!
    !! @param[in,out] this The equation_solver_1var object.
    !! @param[in] x The tolerance value.
    subroutine es1_set_fcn_tol(this, x)
        class(equation_solver_1var), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_fcnTol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on change in variable tolerance.
    !!
    !! @param[in] this The equation_solver_1var object.
    !! @return The tolerance value.
    pure function es1_get_var_tol(this) result(x)
        class(equation_solver_1var), intent(in) :: this
        real(dp) :: x
        x = this%m_xtol
    end function

! --------------------
    !> @brief Sets the convergence on change in variable tolerance.
    !!
    !! @param[in,out] this The equation_solver_1var object.
    !! @param[in] x The tolerance value.
    subroutine es1_set_var_tol(this, x)
        class(equation_solver_1var), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_xtol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets a logical value determining if iteration status should be
    !! printed.
    !!
    !! @param[in] this The equation_solver_1var object.
    !! @return True if the iteration status should be printed; else, false.
    pure function es1_get_print_status(this) result(x)
        class(equation_solver_1var), intent(in) :: this
        logical :: x
        x = this%m_printStatus
    end function

! --------------------
    !> @brief Sets a logical value determining if iteration status should be
    !! printed.
    !!
    !! @param[in,out] this The equation_solver_1var object.
    !! @param[in] x True if the iteration status should be printed; else, false.
    subroutine es1_set_print_status(this, x)
        class(equation_solver_1var), intent(inout) :: this
        logical, intent(in) :: x
        this%m_printStatus = x
    end subroutine

! ******************************************************************************
! EQUATION_OPTIMIZER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Gets the maximum number of function evaluations allowed.
    !!
    !! @param[in] this The equation_optimizer object.
    !! @return The maximum number of function evaluations.
    pure function oe_get_max_eval(this) result(n)
        class(equation_optimizer), intent(in) :: this
        integer(i32) :: n
        n = this%m_maxEval
    end function

! --------------------
    !> @brief Sets the maximum number of function evaluations allowed.
    !!
    !! @param[in,out] this The equation_optimizer object.
    !! @param[in] n The maximum number of function evaluations.
    subroutine oe_set_max_eval(this, n)
        class(equation_optimizer), intent(inout) :: this
        integer(i32), intent(in) :: n
        this%m_maxEval = n
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the tolerance on convergence.
    !!
    !! @param[in] this The equation_optimizer object.
    !! @return The convergence tolerance.
    pure function oe_get_tol(this) result(x)
        class(equation_optimizer), intent(in) :: this
        real(dp) :: x
        x = this%m_tol
    end function

! --------------------
    !> @brief Sets the tolerance on convergence.
    !!
    !! @param[in,out] this The equation_optimizer object.
    !! @param[in] x The convergence tolerance.
    subroutine oe_set_tol(this, x)
        class(equation_optimizer), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_tol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets a logical value determining if iteration status should be
    !! printed.
    !!
    !! @param[in] this The equation_optimizer object.
    !! @return True if the iteration status should be printed; else, false.
    pure function oe_get_print_status(this) result(x)
        class(equation_optimizer), intent(in) :: this
        logical :: x
        x = this%m_printStatus
    end function

! --------------------
    !> @brief Sets a logical value determining if iteration status should be
    !! printed.
    !!
    !! @param[in,out] this The equation_optimizer object.
    !! @param[in] x True if the iteration status should be printed; else, false.
    subroutine oe_set_print_status(this, x)
        class(equation_optimizer), intent(inout) :: this
        logical, intent(in) :: x
        this%m_printStatus = x
    end subroutine

! ******************************************************************************
! MISC. ROUTINES
! ------------------------------------------------------------------------------
    !> @brief Prints the iteration status.
    !!
    !! @param[in] iter The iteration number.
    !! @param[in] nfeval The number of function evaluations.
    !! @param[in] njaceval The number of Jacobian evaluations.
    !! @param[in] xnorm The change in variable value.
    !! @param[in] fnorm The residual.
    subroutine print_status(iter, nfeval, njaceval, xnorm, fnorm)
        ! Arguments
        integer(i32), intent(in) :: iter, nfeval, njaceval
        real(dp), intent(in) :: xnorm, fnorm

        ! Process
        print *, ""
        print '(AI0)', "Iteration: ", iter
        print '(AI0)', "Function Evaluations: ", nfeval
        if (njaceval > 0) print '(AI0)', "Jacobian Evaluations: ", njaceval
        print '(AE9.3)', "Change in Variable: ", xnorm
        print '(AE9.3)', "Residual: ", fnorm
    end subroutine

! ------------------------------------------------------------------------------
end module
