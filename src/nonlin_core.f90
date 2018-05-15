! nonlin_core.f90

!> @mainpage
!!
!! @section intro_sec Introduction
!! NONLIN is a library that provides routines for solving least squares
!! problems, systems of equations, and single equations. The library provides
!! function optimization routines.  There are also specific routines and data
!! types optimized for operation on polynomials.
!!
!! @par
!! The example below illustrates two techniques of determining the roots of a
!! polynomial using routines and data types provided by this library.
!! @image html Cubic_Roots_Example.png
!! @code{.f90}
!! program example
!!     use iso_fortran_env
!!     use nonlin_core
!!     use nonlin_solve
!!     implicit none
!!
!!     ! Local variables
!!     type(fcn1var_helper) :: obj
!!     procedure(fcn1var), pointer :: fcn
!!     type(newton_1var_solver) :: solver
!!     real(real64) :: x, f
!!     type(value_pair) :: limits
!!     type(iteration_behavior) :: tracking
!!
!!     ! Define the search limits
!!     limits%x1 = 2.0d0
!!     limits%x2 = -2.0d0
!!
!!     ! Establish the function
!!     fcn => fcn1
!!     call obj%set_fcn(fcn)
!!
!!     ! Allow the solver to print out updates at each iteration
!!     call solver%set_print_status(.true.)
!!
!!     ! Solve the equation
!!     call solver%solve(obj, x, limits, f, ib = tracking)
!!
!!     ! Print the output, residual, and details regarding the iteration
!!     print *, ""
!!     print '(AF7.5)', "The solution: ", x
!!     print '(AE10.3)', "The residual: ", f
!!     print '(AI0)', "Iterations: ", tracking%iter_count
!!     print '(AI0)', "Function Evaluations: ", tracking%fcn_count
!!     print '(AI0)', "Derivative Evaluations: ", tracking%jacobian_count
!!     print '(AL1)', "Converge on Function Value: ", tracking%converge_on_fcn
!!     print '(AL1)', "Converge on Change in Variable: ", tracking%converge_on_chng
!!     print '(AL1)', "Converge on Derivative: ", tracking%converge_on_zero_diff
!!
!! contains
!!     ! The function:
!!     ! f(x) = x**3 - 2 * x - 1
!!     function fcn1(x) result(f)
!!         real(real64), intent(in) :: x
!!         real(real64) :: f
!!         f = x**3 - 2.0d0 * x - 1.0d0
!!     end function
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! Iteration: 1
!! Function Evaluations: 4
!! Jacobian Evaluations: 2
!! Change in Variable: 0.500E+00
!! Residual: -.125E+00
!!
!! Iteration: 2
!! Function Evaluations: 5
!! Jacobian Evaluations: 3
!! Change in Variable: 0.125E+01
!! Residual: -.208E+01
!!
!! Iteration: 3
!! Function Evaluations: 6
!! Jacobian Evaluations: 4
!! Change in Variable: 0.625E+00
!! Residual: -.115E+01
!!
!! Iteration: 4
!! Function Evaluations: 7
!! Jacobian Evaluations: 5
!! Change in Variable: -.313E+00
!! Residual: 0.436E+00
!!
!! Iteration: 5
!! Function Evaluations: 8
!! Jacobian Evaluations: 6
!! Change in Variable: 0.665E-01
!! Residual: 0.221E-01
!!
!! Iteration: 6
!! Function Evaluations: 9
!! Jacobian Evaluations: 7
!! Change in Variable: 0.375E-02
!! Residual: 0.685E-04
!!
!! The solution: 1.61803
!! The residual:  0.665E-09
!! Iterations: 7
!! Function Evaluations: 11
!! Derivative Evaluations: 8
!! Converge on Function Value: T
!! Converge on Change in Variable: F
!! Converge on Derivative: F
!! @endcode
!! @par
!! An alternative for determining the roots of a polynomial is to use the
!! polynomial type from this library, and determine the roots using the
!! functionallity of that type.  The following example illustrates just such
!! an approach.
!! @code{.f90}
!! program example
!!     use iso_fortran_env
!!     use nonlin_polynomials
!!     implicit none
!!
!!     ! Local Variables
!!     type(polynomial) :: f
!!     real(real64) :: coeffs(4)
!!     complex(real64), allocatable :: rts(:)
!!     integer(int32) :: i
!!
!!     ! Define the polynomial (x**3 - 2 * x - 1)
!!     coeffs = [-1.0d0, -2.0d0, 0.0d0, 1.0d0]
!!     f = coeffs
!!
!!     ! Compute the polynomial roots
!!     rts = f%roots()
!!
!!     ! Display the results
!!     do i = 1, size(rts)
!!         print '(AI0AF9.6AF9.6A)', "Root ", i, " = (", real(rts(i), real64), &
!!             ", ", aimag(rts(i)), ")"
!!     end do
!! end program
!! @endcode
!! The above program produces the following output.
!! @code{.txt}
!! Root 1 = (-1.000000,  0.000000)
!! Root 2 = (-0.618034,  0.000000)
!! Root 3 = ( 1.618034,  0.000000)
!! @endcode

! ------------------------------------------------------------------------------

!> @brief \b nonlin_core
!!
!! @par Purpose
!! To provide various routines to solve equations of one or many variables.
!! Actual solvers are located in the following modules.
!! - Module: nonlin_solve
!!  - quasi_newton_solver
!!  - newton_solver
!!  - brent_solver
!!  - newton_1var_solver
!! - Module: nonlin_least_squares
!!  - least_squares_solver
!! - Module: nonlin_optimize
!!  - nelder_mead
!!  - bfgs
module nonlin_core
    use, intrinsic :: iso_fortran_env, only : real64, int32
    use nonlin_constants
    implicit none
    private
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
    public :: nonlin_optimize_fcn
    public :: print_status

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        !> @brief Describes a function of one variable.
        !!
        !! @param[in] x The independent variable.
        !! @return The value of the function at @p x.
        function fcn1var(x) result(f)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in) :: x
            real(real64) :: f
        end function

        !> @brief Describes an M-element vector-valued function of N-variables.
        !!
        !! @param[in] x An N-element array containing the independent variables.
        !! @param[out] f An M-element array that, on output, contains the values
        !!  of the M functions.
        subroutine vecfcn(x, f)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(out), dimension(:) :: f
        end subroutine

        !> @brief Describes a routine capable of computing the Jacobian matrix
        !! of M functions of N unknowns.
        !!
        !! @param[in] x An N-element array containing the independent variables.
        !! @param[out] jac An M-by-N matrix where the Jacobian will be written.
        subroutine jacobianfcn(x, jac)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(out), dimension(:,:) :: jac
        end subroutine

        !> @brief Describes a function of N variables.
        !!
        !! @param[in] x An N-element array containing the independent variables.
        !! @return The value of the function at @p x.
        function fcnnvar(x) result(f)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
            real(real64) :: f
        end function

        !> @brief Describes a routine capable of computing the gradient vector
        !! of an equation of N variables.
        !!
        !! @param[in] x An N-element array containing the independent variables.
        !! @param[out] g An N-element array where the gradient vector will be
        !!  written as output.
        subroutine gradientfcn(x, g)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(out), dimension(:) :: g
        end subroutine
    end interface

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief Defines a set of parameters that describe the behavior of the
    !! iteration process.
    type :: iteration_behavior
        !> Specifies the number of iterations performed.
        integer(int32) :: iter_count
        !> Specifies the number of function evaluations performed.
        integer(int32) :: fcn_count
        !> Specifies the number of Jacobian evaluations performed.
        integer(int32) :: jacobian_count
        !> Specifies the number of gradient vector evaluations performed.
        integer(int32) :: gradient_count
        !> True if the solution converged as a result of a zero-valued
        !! function; else, false.
        logical :: converge_on_fcn
        !> True if the solution converged as a result of no appreciable
        !! change in solution points between iterations; else, false.
        logical :: converge_on_chng
        !> True if the solution appears to have settled on a stationary
        !! point such that the gradient of the function is zero-valued; else,
        !! false.
        logical :: converge_on_zero_diff
    end type

! ------------------------------------------------------------------------------
    !> @brief Defines a pair of numeric values.
    type :: value_pair
        !> Value 1.
        real(real64) :: x1
        !> Value 2.
        real(real64) :: x2
    end type

! ******************************************************************************
! NONLIN_VECFCN_HELPER.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a type capable of encapsulating a system of nonlinear
    !! equations of the form: F(X) = 0.  This type is used to establish the
    !! system of equations to solve, and provides a means for computing
    !! the Jacobian matrix for the system of equations, and any other
    !! ancillary operations that may be needed by the solver.
    !!
    !! @par Example
    !! The following example illustrates the most basic use of this type to
    !! solve a system of 2 equations and 2 unknowns using Newton's method.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use nonlin_core
    !!     use nonlin_solve
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     type(vecfcn_helper) :: obj
    !!     procedure(vecfcn), pointer :: fcn
    !!     type(newton_solver) :: solver
    !!     real(real64) :: x(2), f(2)
    !!
    !!     ! Assign a pointer to the subroutine containing the equations to solve
    !!     fcn => fcns
    !!     call obj%set_fcn(fcn, 2, 2) ! There are 2 equations with 2 unknowns
    !!
    !!     ! Define an initial guess
    !!     x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]
    !!
    !!     ! Solve the system of equations
    !!     call solver%solve(obj, x, f)
    !!
    !!     ! Display the output
    !!     print '(AF7.5AF7.5A)', "Solution: (", x(1), ", ", x(2), ")"
    !!     print '(AE9.3AE9.3A)', "Residual: (", f(1), ", ", f(2), ")"
    !! contains
    !!     ! Define the routine containing the equations to solve.  The equations are:
    !!     ! x**2 + y**2 = 34
    !!     ! x**2 - 2 * y**2 = 7
    !!     subroutine fcns(x, f)
    !!         real(real64), intent(in), dimension(:) :: x
    !!         real(real64), intent(out), dimension(:) :: f
    !!         f(1) = x(1)**2 + x(2)**2 - 34.0d0
    !!         f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    !!     end subroutine
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! Solution: (5.00000, 3.00000)
    !! Residual: (0.000E+00, 0.000E+00)
    !! @endcode
    type vecfcn_helper
        private
        !> A pointer to the target vecfcn routine.
        procedure(vecfcn), pointer, nopass :: m_fcn => null()
        !> A pointer to the jacobian routine - null if no routine is supplied.
        procedure(jacobianfcn), pointer, nopass :: m_jac => null()
        !> The number of functions in m_fcn
        integer(int32) :: m_nfcn = 0
        !> The number of variables in m_fcn
        integer(int32) :: m_nvar = 0
    contains
        !> @brief Establishes a pointer to the routine containing the system of
        !!  equations to solve.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_fcn(class(vecfcn_helper) this, procedure(vecfcn) pointer fcn, integer(int32) nfcn, integer(int32) nvar)
        !! @endcode
        !!
        !! @param[in,out] this The vecfcn_helper object.
        !! @param[in] fcn The function pointer.
        !! @param[in] nfcn The number of functions.
        !! @param[in] nvar The number of variables.
        !!
        !! @par Example
        !! The following example illustrates how to define the function to
        !! solve.  Newton's method is being utilized via the newton_solver
        !! type.
        !! @code{.f90}
        !! program example
        !!     use iso_fortran_env
        !!     use nonlin_core
        !!     use nonlin_solve
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     type(vecfcn_helper) :: obj
        !!     procedure(vecfcn), pointer :: fcn
        !!     type(newton_solver) :: solver
        !!     real(real64) :: x(2), f(2)
        !!
        !!     ! Assign a pointer to the subroutine containing the equations to solve
        !!     fcn => fcns
        !!     call obj%set_fcn(fcn, 2, 2) ! There are 2 equations with 2 unknowns
        !!
        !!     ! Define an initial guess
        !!     x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]
        !!
        !!     ! Solve the system of equations
        !!     call solver%solve(obj, x, f)
        !!
        !!     ! Display the output
        !!     print '(AF7.5AF7.5A)', "Solution: (", x(1), ", ", x(2), ")"
        !!     print '(AE9.3AE9.3A)', "Residual: (", f(1), ", ", f(2), ")"
        !! contains
        !!     ! Define the routine containing the equations to solve.  The equations are:
        !!     ! x**2 + y**2 = 34
        !!     ! x**2 - 2 * y**2 = 7
        !!     subroutine fcns(x, f)
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64), intent(out), dimension(:) :: f
        !!         f(1) = x(1)**2 + x(2)**2 - 34.0d0
        !!         f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
        !!     end subroutine
        !! end program
        !! @endcode
        !! The above program produces the following output.
        !! @code{.txt}
        !! Solution: (5.00000, 3.00000)
        !! Residual: (0.000E+00, 0.000E+00)
        !! @endcode
        procedure, public :: set_fcn => vfh_set_fcn
        !> @brief Establishes a pointer to the routine for computing the
        !! Jacobian matrix of the system of equations.  If no routine is
        !! defined, the Jacobian matrix will be computed numerically (this is
        !! the default state).
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_jacobian(class(vecfcn_helper) this, procedure(jacobianfcn) pointer jac)
        !! @endcode
        !!
        !! @param[in,out] this The vecfcn_helper object.
        !! @param[in] jac The function pointer.
        !!
        !! @par Example
        !! The following example utilizes Newton's method to solve a system of
        !! 2 equations and 2 unknowns with a user-defined Jacobian.
        !! @code{.f90}
        !! program example
        !!     use iso_fortran_env
        !!     use nonlin_core
        !!     use nonlin_solve
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     type(vecfcn_helper) :: obj
        !!     procedure(vecfcn), pointer :: fcn
        !!     procedure(jacobianfcn), pointer :: jac
        !!     type(newton_solver) :: solver
        !!     real(real64) :: x(2), f(2)
        !!
        !!     ! Assign the function and Jacobian routines
        !!     fcn => fcns
        !!     jac => fcnjac
        !!     call obj%set_fcn(fcn, 2, 2)
        !!     call obj%set_jacobian(jac)
        !!
        !!     ! Define an initial guess
        !!     x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]
        !!
        !!     ! Solve the system of equations
        !!     call solver%solve(obj, x, f)
        !!
        !!     ! Display the output
        !!     print '(AF7.5AF7.5A)', "Solution: (", x(1), ", ", x(2), ")"
        !!     print '(AE9.3AE9.3A)', "Residual: (", f(1), ", ", f(2), ")"
        !! contains
        !!     ! The system of equations (source: https://www.mathworks.com/help/optim/ug/fsolve.html)
        !!     ! 2 * x1 - x2 = exp(-x1)
        !!     ! -x1 + 2 * x2 = exp(-x2)
        !!     subroutine fcns(x, f)
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64), intent(out), dimension(:) :: f
        !!         f(1) = 2.0d0 * x(1) - x(2) - exp(-x(1))
        !!         f(2) = -x(1) + 2.0d0 * x(2) - exp(-x(2))
        !!     end subroutine
        !!
        !!     ! The Jacobian matrix:
        !!     !     | exp(-x1) + 2          -1     |
        !!     ! J = |                              |
        !!     !     |     -1          exp(-x2) + 2 |
        !!     subroutine fcnjac(x, jac)
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64), intent(out), dimension(:,:) :: jac
        !!         jac(1,1) = exp(-x(1)) + 2.0d0
        !!         jac(2,1) = -1.0d0
        !!         jac(1,2) = -1.0d0
        !!         jac(2,2) = exp(-x(2)) + 2.0d0
        !!     end subroutine
        !! end program
        !! @endcode
        !! The above program produces the following output.
        !! @code{.txt}
        !! Solution: (0.56714, 0.56714)
        !! Residual: (-.693E-08, -.683E-08)
        !! @endcode
        procedure, public :: set_jacobian => vfh_set_jac
        !> @brief Tests if the pointer to the subroutine containing the system
        !! of equations to solve has been assigned.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function is_fcn_defined(class(vecfcn_helper) this)
        !! @endcode
        !!
        !! @param[in] this The vecfcn_helper object.
        !! @return Returns true if the pointer has been assigned; else, false.
        procedure, public :: is_fcn_defined => vfh_is_fcn_defined
        !> @brief Tests if the pointer to the subroutine containing the system
        !! of equations to solve has been assigned.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function is_jacobian_defined(class(vecfcn_helper) this)
        !! @endcode
        !!
        !! @param[in] this The vecfcn_helper object.
        !! @return Returns true if the pointer has been assigned; else, false.
        procedure, public :: is_jacobian_defined => vfh_is_jac_defined
        !> @brief Executes the routine containing the system of equations to
        !! solve.  No action is taken if the pointer to the subroutine has not
        !! been defined.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine fcn(class(vecfcn_helper) this, real(real64) x(:), real(real64) f(:))
        !! @endcode
        !!
        !! @param[in] this The vecfcn_helper object.
        !! @param[in] x An N-element array containing the independent variables.
        !! @param[out] f An M-element array that, on output, contains the values
        !!  of the M functions.
        procedure, public :: fcn => vfh_fcn
        !> @brief Executes the routine containing the Jacobian matrix if
        !! supplied.  If not supplied, the Jacobian is computed via finite
        !! differences.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine jacobian(class(vecfcn_helper) this real(real64) x(:), &
        !!  real(real64) jac(:), optional real(real64) fv(:), &
        !!  optional real(real64) work(:), optional integer(int32) lwork, &
        !!  optional integer(int32) err)
        !! @endcode
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
        procedure, public :: jacobian => vfh_jac_fcn
        !> @brief Gets the number of equations in this system.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) get_equation_count(class(vecfcn_helper) this)
        !! @endcode
        !!
        !! @param[in] this The vecfcn_helper object.
        !! @return The function count.
        procedure, public :: get_equation_count => vfh_get_nfcn
        !> @brief Gets the number of variables in this system.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) get_variable_count(class(vecfcn_helper) this)
        !! @endcode
        !!
        !! @param[in] this The vecfcn_helper object.
        !! @return The number of variables.
        procedure, public :: get_variable_count => vfh_get_nvar
    end type

! ------------------------------------------------------------------------------
    interface
        module subroutine vfh_set_fcn(this, fcn, nfcn, nvar)
            class(vecfcn_helper), intent(inout) :: this
            procedure(vecfcn), intent(in), pointer :: fcn
            integer(int32), intent(in) :: nfcn, nvar
        end subroutine

        module subroutine vfh_set_jac(this, jac)
            class(vecfcn_helper), intent(inout) :: this
            procedure(jacobianfcn), intent(in), pointer :: jac
        end subroutine

        pure module function vfh_is_fcn_defined(this) result(x)
            class(vecfcn_helper), intent(in) :: this
            logical :: x
        end function

        pure module function vfh_is_jac_defined(this) result(x)
            class(vecfcn_helper), intent(in) :: this
            logical :: x
        end function

        module subroutine vfh_fcn(this, x, f)
            class(vecfcn_helper), intent(in) :: this
            real(real64), intent(in), dimension(:) :: x
            real(real64), intent(out), dimension(:) :: f
        end subroutine

        module subroutine vfh_jac_fcn(this, x, jac, fv, work, olwork, err)
            class(vecfcn_helper), intent(in) :: this
            real(real64), intent(inout), dimension(:) :: x
            real(real64), intent(out), dimension(:,:) :: jac
            real(real64), intent(in), dimension(:), optional, target :: fv
            real(real64), intent(out), dimension(:), optional, target :: work
            integer(int32), intent(out), optional :: olwork, err
        end subroutine

        pure module function vfh_get_nfcn(this) result(n)
            class(vecfcn_helper), intent(in) :: this
            integer(int32) :: n
        end function

        pure module function vfh_get_nvar(this) result(n)
            class(vecfcn_helper), intent(in) :: this
            integer(int32) :: n
        end function
    end interface

! ******************************************************************************
! NONLIN_FCN1VAR_HELPER.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a type capable of encapsulating an equation of one
    !! variable of the form: f(x) = 0.
    !!
    !! @par Example
    !! The following example illustrates the use of this type to solve an
    !! equation using Brent's method.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use nonlin_core
    !!     use nonlin_solve
    !!     implicit none
    !!
    !!     ! Local variables
    !!     type(fcn1var_helper) :: obj
    !!     procedure(fcn1var), pointer :: fcn
    !!     type(brent_solver) :: solver
    !!     real(real64) :: x, f
    !!     type(value_pair) :: limits
    !!
    !!     ! Define the search limits
    !!     limits%x1 = 1.5d0
    !!     limits%x2 = 5.0d0
    !!
    !!     ! Establish the function
    !!     fcn => fcn1
    !!     call obj%set_fcn(fcn)
    !!
    !!     ! Solve the equation
    !!     call solver%solve(obj, x, limits, f)
    !!
    !!     ! Print the output and the residual
    !!     print '(AF7.5)', "The solution: ", x
    !!     print '(AE10.3)', "The residual: ", f
    !!
    !! contains
    !!     ! The function:
    !!     ! f(x) = sin(x) / x, solution: x = x * pi for n = 1, 2, 3, ...
    !!     function fcn1(x) result(f)
    !!         real(real64), intent(in) :: x
    !!         real(real64) :: f
    !!         f = sin(x) / x
    !!     end function
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! The solution: 3.14159
    !! The residual: -0.751E-11
    !! @endcode
    type fcn1var_helper
        private
        !> A pointer to the target fcn1var routine.
        procedure(fcn1var), pointer, nopass :: m_fcn => null()
        !> A pointer to a function capable of computing the derivative of m_fcn.
        procedure(fcn1var), pointer, nopass :: m_diff => null()
    contains
        !> @brief Executes the routine containing the function to evaluate.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function fcn(class(fcn1var_helper) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in] this The fcn1var_helper object.
        !! @param[in] x The value of the independent variable at which the function
        !!  should be evaluated.
        !! @return The value of the function at @p x.
        procedure, public :: fcn => f1h_fcn
        !> @brief Tests if the pointer to the function containing the equation
        !! to solve has been assigned.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function is_fcn_defined(class(fcn1var_helper) this)
        !! @endcode
        !!
        !! @param[in] this The fcn1var_helper object.
        !! @return Returns true if the pointer has been assigned; else, false.
        procedure, public :: is_fcn_defined => f1h_is_fcn_defined
        !> @brief Establishes a pointer to the routine containing the equations
        !! to solve.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_fcn(class(fcn1var_helper) this, procedure(fcn1var) pointer fcn)
        !! @endcode
        !!
        !! @param[in,out] this The fcn1var_helper object.
        !! @param[in] fcn The function pointer.
        !!
        !! @par Example
        !! The following example illustrates how to use this routine to
        !! inform the solver of the function to solve.  This particular
        !! example utilizes Brent's method to compute the solution.
        !! @code{.f90}
        !! program example
        !!     use iso_fortran_env
        !!     use nonlin_core
        !!     use nonlin_solve
        !!     implicit none
        !!
        !!     ! Local variables
        !!     type(fcn1var_helper) :: obj
        !!     procedure(fcn1var), pointer :: fcn
        !!     type(brent_solver) :: solver
        !!     real(real64) :: x, f
        !!     type(value_pair) :: limits
        !!
        !!     ! Define the search limits
        !!     limits%x1 = 1.5d0
        !!     limits%x2 = 5.0d0
        !!
        !!     ! Establish the function
        !!     fcn => fcn1
        !!     call obj%set_fcn(fcn)
        !!
        !!     ! Solve the equation
        !!     call solver%solve(obj, x, limits, f)
        !!
        !!     ! Print the output and the residual
        !!     print '(AF7.5)', "The solution: ", x
        !!     print '(AE10.3)', "The residual: ", f
        !!
        !! contains
        !!     ! The function:
        !!     ! f(x) = sin(x) / x, solution: x = x * pi for n = 1, 2, 3, ...
        !!     function fcn1(x) result(f)
        !!         real(real64), intent(in) :: x
        !!         real(real64) :: f
        !!         f = sin(x) / x
        !!     end function
        !! end program
        !! @endcode
        !! The above program produces the following output.
        !! @code{.txt}
        !! The solution: 3.14159
        !! The residual: -0.751E-11
        !! @endcode
        procedure, public :: set_fcn => f1h_set_fcn
        !> @brief Tests if the pointer to the function containing the
        !! derivative of the function to solve is defined.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function is_derivative_defined(class(fcn1var_helper) this)
        !! @endcode
        !!
        !! @param[in] this The fcn1var_helper object.
        !! @return Returns true if the pointer has been assigned; else, false.
        procedure, public :: is_derivative_defined => f1h_is_diff_defined
        !> @brief Computes the derivative of the function.  If a routine for
        !! computing the derivative is not defined, the derivative is estimated
        !! via finite differences.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function derivative(class(fcn1var_helper) this, real(real64) x, optional real(real64) f)
        !! @endcode
        !!
        !! @param[in] this The fcn1var_helper object.
        !! @param[in] x The value of the independent variable at which the derivative is to be computed.
        !! @param[in] f An optional input specifying the function value at @p x.  If supplied, and the
        !!  derivative is being estimated numerically, the function will not be evaluated at @p x.
        procedure, public :: diff => f1h_diff_fcn
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_diff(class(fcn1var_helper) this, procedure(fcn1var) pointer diff)
        !! @endcode
        !!
        !! @param[in,out] this The fcn1var_helper object.
        !! @param[in] diff A pointer to the function for computing the first
        !!  derivative.
        procedure, public :: set_diff => f1h_set_diff
    end type

! ------------------------------------------------------------------------------
    interface
        module function f1h_fcn(this, x) result(f)
            class(fcn1var_helper), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64) :: f
        end function

        pure module function f1h_is_fcn_defined(this) result(x)
            class(fcn1var_helper), intent(in) :: this
            logical :: x
        end function

        module subroutine f1h_set_fcn(this, fcn)
            class(fcn1var_helper), intent(inout) :: this
            procedure(fcn1var), intent(in), pointer :: fcn
        end subroutine

        pure module function f1h_is_diff_defined(this) result(x)
            class(fcn1var_helper), intent(in) :: this
            logical :: x
        end function

        module function f1h_diff_fcn(this, x, f) result(df)
            class(fcn1var_helper), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64), intent(in), optional :: f
            real(real64) :: df
        end function

        module subroutine f1h_set_diff(this, diff)
            class(fcn1var_helper), intent(inout) :: this
            procedure(fcn1var), pointer, intent(in) :: diff
        end subroutine
    end interface

! ******************************************************************************
! NONLIN_FCNNVAR_HELPER.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a type capable of encapsulating an equation of N
    !! variables.
    !!
    !! @par Example
    !! The following example illustrates how to use this type to find the
    !! minimum of a function of multiple variables.  In this instance the
    !! Nelder-Mead simplex method is utilized to minimize the Rosenbrock
    !! function in two variables.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use nonlin_optimize
    !!     use nonlin_core
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     type(nelder_mead) :: solver
    !!     type(fcnnvar_helper) :: obj
    !!     procedure(fcnnvar), pointer :: fcn
    !!     real(real64) :: x(2), f
    !!
    !!     ! Initialization
    !!     fcn => rosenbrock
    !!     call obj%set_fcn(fcn, 2)
    !!
    !!     ! Define an initial guess - the solution is (1, 1)
    !!     call random_number(x)
    !!
    !!     ! Call the solver
    !!     call solver%solve(obj, x, f)
    !!
    !!      ! Display the output
    !!      print '(AF7.5AF7.5A)', "Minimum: (", x(1), ", ", x(2), ")"
    !!      print '(AE9.3)', "Function Value: ", f
    !! contains
    !!     ! Rosenbrock's Function
    !!     function rosenbrock(x) result(f)
    !!         real(real64), intent(in), dimension(:) :: x
    !!         real(real64) :: f
    !!         f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    !!     end function
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! Minimum: (1.00000, 1.00000)
    !! Function Value: 0.213E-12
    !! @endcode
    type fcnnvar_helper
        private
        !> A pointer to the target fcnnvar routine.
        procedure(fcnnvar), pointer, nopass :: m_fcn => null()
        !> A pointer to the gradient routine.
        procedure(gradientfcn), pointer, nopass :: m_grad => null()
        !> The number of variables in m_fcn
        integer(int32) :: m_nvar = 0
    contains
        !> @brief Executes the routine containing the function to evaluate.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function fcn(class(fcnnvar_helper) this, real(real64) x(:))
        !! @endcode
        !!
        !! @param[in] this The fcnnvar_helper object.
        !! @param[in] x The value of the independent variable at which the function
        !!  should be evaluated.
        !! @return The value of the function at @p x.
        procedure, public :: fcn => fnh_fcn
        !> @brief Tests if the pointer to the function has been assigned.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function is_fcn_defined(class(fcnnvar_helper) this)
        !! @endcode
        !!
        !! @param[in] this The fcnnvar_helper object.
        !! @return Returns true if the pointer has been assigned; else, false.
        procedure, public :: is_fcn_defined => fnh_is_fcn_defined
        !> @brief Establishes a pointer to the routine containing the function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_fcn(class(fcnnvar_helper) this, procedure(fcnnvar) pointer fcn, integer(int32) nvar)
        !! @endcode
        !!
        !! @param[in,out] this The fcnnvar_helper object.
        !! @param[in] fcn The function pointer.
        !! @param[in] nvar The number of variables in the function.
        !!
        !! @par Example
        !! The following example illustrates how to use this routine to
        !! inform the solver of the function to optimize.  This particular
        !! example utilizes the Nelder-Mead simplex method to minimize the
        !! function.
        !! @code{.f90}
        !! program example
        !!     use iso_fortran_env
        !!     use nonlin_optimize
        !!     use nonlin_core
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     type(nelder_mead) :: solver
        !!     type(fcnnvar_helper) :: obj
        !!     procedure(fcnnvar), pointer :: fcn
        !!     real(real64) :: x(2), f
        !!
        !!     ! Initialization
        !!     fcn => rosenbrock
        !!     call obj%set_fcn(fcn, 2)
        !!
        !!     ! Define an initial guess - the solution is (1, 1)
        !!     call random_number(x)
        !!
        !!     ! Call the solver
        !!     call solver%solve(obj, x, f)
        !!
        !!      ! Display the output
        !!      print '(AF7.5AF7.5A)', "Minimum: (", x(1), ", ", x(2), ")"
        !!      print '(AE9.3)', "Function Value: ", f
        !! contains
        !!     ! Rosenbrock's Function
        !!     function rosenbrock(x) result(f)
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64) :: f
        !!         f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
        !!     end function
        !! end program
        !! @endcode
        !! The above program produces the following output.
        !! @code{.txt}
        !! Minimum: (1.00000, 1.00000)
        !! Function Value: 0.213E-12
        !! @endcode
        procedure, public :: set_fcn => fnh_set_fcn
        !> @brief Gets the number of variables in this system.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_variable_count(class(fcnnvar_helper) this)
        !! @endcode
        !!
        !! @param[in] this The fcnnvar_helper object.
        !! @return The number of variables.
        procedure, public :: get_variable_count => fnh_get_nvar
        !> @brief Establishes a pointer to the routine containing the gradient
        !! vector of the function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_gradient_fcn(class(fcnnvar_helper) this, procedure(gradientfcn) pointer fcn)
        !! @endcode
        !!
        !! @param[in,out] this The fcnnvar_helper object.
        !! @param[in] fcn The pointer to the gradient routine.
        !!
        !! @par Example
        !! The following example illustrates the use of a user-defined
        !! gradient function.  This particular example utilizes the BFGS
        !! method to minimize Beale's function.
        !! @code{.f90}
        !! program example
        !!     use iso_fortran_env
        !!     use nonlin_core
        !!     use nonlin_optimize
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     type(bfgs) :: solver
        !!     type(fcnnvar_helper) :: obj
        !!     procedure(fcnnvar), pointer :: fcn
        !!     procedure(gradientfcn), pointer :: grad
        !!     real(real64) :: x(2), f
        !!
        !!     ! Tell the solver where to find the function
        !!     fcn => beale
        !!     grad => bealegrad
        !!     call obj%set_fcn(fcn, 2)
        !!     call obj%set_gradient_fcn(grad)
        !!
        !!     ! Define an initial guess
        !!     x = 1.0d0
        !!
        !!     ! Compute the solution
        !!     call solver%solve(obj, x, f)
        !!
        !!     ! Display the output
        !!     print '(AF7.5AF7.5A)', "Minimum: (", x(1), ", ", x(2), ")"
        !!     print '(AE9.3)', "Function Value: ", f
        !! contains
        !!     ! The Beale function:
        !!     ! f(x) = (1.5 - x + xy)**2 + (2.25 - x + xy**2)**2 + (2.625 - x + xy**3)**2
        !!     ! The minimum is at x = 3, y = 0.5, and f(3, 0.5) = 0
        !!     function beale(x) result(f)
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64) :: f
        !!         f = (1.5d0 - x(1) + x(1) * x(2))**2 + &
        !!             (2.25d0 - x(1) + x(1) * x(2)**2)**2 + &
        !!             (2.625d0 - x(1) + x(1) * x(2)**3)**2
        !!     end function
        !!
        !!     ! The gradient
        !!     subroutine bealegrad(x, g)
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64), intent(out), dimension(:) :: g
        !!
        !!         g(1) = 2.0d0 * (x(2)**3 - 1.0d0) * (x(1) * x(2)**3 - x(1) + 2.625d0) + &
        !!             2.0d0 * (x(2)**2 - 1.0d0) * (x(1) * x(2)**2 - x(1) + 2.25d0) + &
        !!             2.0d0 * (x(2) - 1.0d0) * (x(1) * x(2) - x(1) + 1.5d0)
        !!
        !!         g(2) = 6.0d0 * x(1) * x(2)**2 * (x(1) * x(2)**3 - x(1) + 2.625d0) + &
        !!             4.0d0 * x(1) * x(2) * (x(1) * x(2)**2 - x(1) + 2.25d0) + &
        !!             2.0d0 * x(1) * (x(1) * x(2) - x(1) + 1.5d0)
        !!     end subroutine
        !! end program
        !! @endcode
        !! The above program produces the following output.
        !! @code{.txt}
        !! Minimum: (3.00000, 0.50000)
        !! Function Value: 0.999E-28
        !! @endcode
        procedure, public :: set_gradient_fcn => fnh_set_grad
        !> @brief Tests if the pointer to the routine containing the gradient
        !! has been assigned.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function is_gradient_defined(class(fcnnvar_helper) this)
        !! @endcode
        !!
        !! @param[in] this The fcnnvar_helper object.
        !! @return Returns true if the pointer has been assigned; else, false.
        procedure, public :: is_gradient_defined => fnh_is_grad_defined
        !> @brief Computes the gradient of the function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine gradient(class(fcnnvar_helper) this, real(real64) x(:), &
        !!  real(real64) g(:), optional real(real64) fv(:), &
        !!  optional integer(int32) err)
        !! @endcode
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
        procedure, public :: gradient => fnh_grad_fcn
    end type

! ------------------------------------------------------------------------------
    interface
        module function fnh_fcn(this, x) result(f)
            class(fcnnvar_helper), intent(in) :: this
            real(real64), intent(in), dimension(:) :: x
            real(real64) :: f
        end function

        pure module function fnh_is_fcn_defined(this) result(x)
            class(fcnnvar_helper), intent(in) :: this
            logical :: x
        end function

        module subroutine fnh_set_fcn(this, fcn, nvar)
            class(fcnnvar_helper), intent(inout) :: this
            procedure(fcnnvar), intent(in), pointer :: fcn
            integer(int32), intent(in) :: nvar
        end subroutine

        pure module function fnh_get_nvar(this) result(n)
            class(fcnnvar_helper), intent(in) :: this
            integer(int32) :: n
        end function

        module subroutine fnh_set_grad(this, fcn)
            class(fcnnvar_helper), intent(inout) :: this
            procedure(gradientfcn), pointer, intent(in) :: fcn
        end subroutine

        pure module function fnh_is_grad_defined(this) result(x)
            class(fcnnvar_helper), intent(in) :: this
            logical :: x
        end function

        module subroutine fnh_grad_fcn(this, x, g, fv, err)
            class(fcnnvar_helper), intent(in) :: this
            real(real64), intent(inout), dimension(:) :: x
            real(real64), intent(out), dimension(:) :: g
            real(real64), intent(in), optional :: fv
            integer(int32), intent(out), optional :: err
        end subroutine
    end interface

! ******************************************************************************
! NONLIN_EQUATION_SOLVER.F90
! ------------------------------------------------------------------------------
    !> @brief A base class for various solvers of nonlinear systems of
    !! equations.
    !!
    !! @par Example
    !! The following example illustrates the use of the newton_solver type,
    !! which includes this type in its inheritance chain, to solve a system
    !! of equations.  Several options are utilized in the example to illustrate
    !! much of the solver functionallity.  Notice, many of these options
    !! are not necessary to compute the solution; however, they do provide
    !! means of control over the solver.  For a more minimalistic example see
    !! the vecfcn_helper type.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use nonlin_core
    !!     use nonlin_solve
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     type(vecfcn_helper) :: obj
    !!     procedure(vecfcn), pointer :: fcn
    !!     type(newton_solver) :: solver
    !!     real(real64) :: x(2), f(2)
    !!     type(iteration_behavior) :: tracking
    !!
    !!     ! Assign a pointer to the subroutine containing the equations to solve
    !!     fcn => fcns
    !!     call obj%set_fcn(fcn, 2, 2) ! There are 2 equations with 2 unknowns
    !!
    !!     ! Set up solver parameters
    !!     call solver%set_max_fcn_evals(1000) ! Specify the maximum number of function evaluations before iteration termination
    !!     call solver%set_fcn_tolerance(1.0d-10)
    !!     call solver%set_var_tolerance(1.0d-10)
    !!     call solver%set_gradient_tolerance(1.0d-10)
    !!
    !!     ! Tell the solver to print out a status update at each iteration - the default behavior does not print updates
    !!     call solver%set_print_status(.true.)
    !!
    !!     ! Define an initial guess
    !!     x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]
    !!
    !!     ! Solve the system of equations, but include solver statistics tracking
    !!     call solver%solve(obj, x, f, tracking)
    !!
    !!     ! Display the output
    !!     print *, ""
    !!     print '(AF7.5AF7.5A)', "Solution: (", x(1), ", ", x(2), ")"
    !!     print '(AE9.3AE9.3A)', "Residual: (", f(1), ", ", f(2), ")"
    !!     print '(AI0)', "Iteration Count: ", tracking%iter_count
    !!     print '(AI0)', "Function Evaluations: ", tracking%fcn_count
    !!     print '(AI0)', "Jacobian Evaluations: ", tracking%jacobian_count
    !!     print '(AL1)', "Converge on Function Value: ", tracking%converge_on_fcn
    !!     print '(AL1)', "Converge on Change in Variable: ", tracking%converge_on_chng
    !!     print '(AL1)', "Converge on Zero Slope Gradient Vector: ", tracking%converge_on_zero_diff
    !! contains
    !!     ! Define the routine containing the equations to solve.  The equations are:
    !!     ! x**2 + y**2 = 34
    !!     ! x**2 - 2 * y**2 = 7
    !!     subroutine fcns(x, f)
    !!         real(real64), intent(in), dimension(:) :: x
    !!         real(real64), intent(out), dimension(:) :: f
    !!         f(1) = x(1)**2 + x(2)**2 - 34.0d0
    !!         f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    !!     end subroutine
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! Iteration: 1
    !! Function Evaluations: 3
    !! Jacobian Evaluations: 1
    !! Change in Variable: 0.545E+00
    !! Residual: 0.272E+02
    !!
    !! Iteration: 2
    !! Function Evaluations: 5
    !! Jacobian Evaluations: 2
    !! Change in Variable: 0.504E+00
    !! Residual: 0.743E+01
    !!
    !! Iteration: 3
    !! Function Evaluations: 6
    !! Jacobian Evaluations: 3
    !! Change in Variable: 0.132E+00
    !! Residual: 0.522E+00
    !!
    !! Iteration: 4
    !! Function Evaluations: 7
    !! Jacobian Evaluations: 4
    !! Change in Variable: 0.882E-02
    !! Residual: 0.199E-02
    !!
    !! Iteration: 5
    !! Function Evaluations: 8
    !! Jacobian Evaluations: 5
    !! Change in Variable: 0.389E-04
    !! Residual: 0.302E-07
    !!
    !! Solution: (5.00000, 3.00000)
    !! Residual: (0.000E+00, 0.000E+00)
    !! Iteration Count: 6
    !! Function Evaluations: 9
    !! Jacobian Evaluations: 6
    !! Converge on Function Value: T
    !! Converge on Change in Variable: F
    !! Converge on Zero Slope Gradient Vector: F
    !! @endcode
    type, abstract :: equation_solver
        private
        !> The maximum number of function evaluations allowed per solve.
        integer(int32) :: m_maxEval = 100
        !> The convergence criteria on function values.
        real(real64) :: m_fcnTol = 1.0d-8
        !> The convergence criteria on change in variable values.
        real(real64) :: m_xtol = 1.0d-12
        !> The convergence criteria for the slope of the gradient vector.
        real(real64) :: m_gtol = 1.0d-12
        !> Set to true to print iteration status; else, false.
        logical :: m_printStatus = .false.
    contains
        !> @brief Gets the maximum number of function evaluations allowed during
        !! a single solve.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_max_fcn_evals(class(equation_solver) this)
        !! @endcode
        !!
        !! @param[in] this The equation_solver object.
        !! @return The maximum number of function evaluations.
        procedure, public :: get_max_fcn_evals => es_get_max_eval
        !> @brief Sets the maximum number of function evaluations allowed during
        !! a single solve.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_max_fcn_evals(class(equation_solver) this, integer(int32) n)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver object.
        !! @param[in] n The maximum number of function evaluations.
        !!
        !! @par Example
        !! See the equation_solver type for an example on how to establish
        !! the limit on number of function evaluations.
        procedure, public :: set_max_fcn_evals => es_set_max_eval
        !> @brief Gets the convergence on function value tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) get_fcn_tolerance(class(equation_solver) this)
        !! @endcode
        !!
        !! @param[in] this The equation_solver object.
        !! @return The tolerance value.
        procedure, public :: get_fcn_tolerance => es_get_fcn_tol
        !> @brief Sets the convergence on function value tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_fcn_tolerance(class(equation_solver) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver object.
        !! @param[in] x The tolerance value.
        !!
        !! @par Example
        !! See the equation_solver type for an example on how to establish
        !! the function value tolerance.
        procedure, public :: set_fcn_tolerance => es_set_fcn_tol
        !> @brief Gets the convergence on change in variable tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_var_tolerance(class(equation_solver) this)
        !! @endcode
        !!
        !! @param[in] this The equation_solver object.
        !! @return The tolerance value.
        procedure, public :: get_var_tolerance => es_get_var_tol
        !> @brief Sets the convergence on change in variable tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_var_tolerance(class(equation_solver) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver object.
        !! @param[in] x The tolerance value.
        !!
        !! @par Example
        !! See the equation_solver type for an example on how to establish
        !! the change in variable tolerance.
        procedure, public :: set_var_tolerance => es_set_var_tol
        !> @brief Gets the convergence on slope of the gradient vector
        !! tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_gradient_tolerance(class(equation_solver))
        !! @endcode
        !!
        !! @param[in] this The equation_solver object.
        !! @return The tolerance value.
        procedure, public :: get_gradient_tolerance => es_get_grad_tol
        !> @brief Sets the convergence on slope of the gradient vector
        !! tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_gradient_tolerance(class(equation_solver) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver object.
        !! @param[in] x The tolerance value.
        !!
        !! @par Example
        !! See the equation_solver type for an example on how to establish
        !! the gradient vector slope tolerance.
        procedure, public :: set_gradient_tolerance => es_set_grad_tol
        !> @brief Gets a logical value determining if iteration status should be
        !! printed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical get_print_status(class(equation_solver) this)
        !! @endcode
        !!
        !! @param[in] this The equation_solver object.
        !! @return True if the iteration status should be printed; else, false.
        procedure, public :: get_print_status => es_get_print_status
        !> @brief Sets a logical value determining if iteration status should be
        !! printed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_print_status(class(equation_solver) this, logical x)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver object.
        !! @param[in] x True if the iteration status should be printed; else, false.
        !!
        !! @par Example
        !! See the equation_solver type for an example on how to enable
        !! iteration update printing.
        procedure, public :: set_print_status => es_set_print_status
        !> @brief Solves the system of equations.
        !!
        !! @par Example
        !! See the equation_solver type for an example on how to solve a
        !! system of equations.
        procedure(nonlin_solver), deferred, public, pass :: solve
    end type

! ------------------------------------------------------------------------------
    interface
        pure module function es_get_max_eval(this) result(n)
            class(equation_solver), intent(in) :: this
            integer(int32) :: n
        end function

        module subroutine es_set_max_eval(this, n)
            class(equation_solver), intent(inout) :: this
            integer(int32), intent(in) :: n
        end subroutine

        pure module function es_get_fcn_tol(this) result(x)
            class(equation_solver), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine es_set_fcn_tol(this, x)
            class(equation_solver), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function es_get_var_tol(this) result(x)
            class(equation_solver), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine es_set_var_tol(this, x)
            class(equation_solver), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function es_get_grad_tol(this) result(x)
            class(equation_solver), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine es_set_grad_tol(this, x)
            class(equation_solver), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function es_get_print_status(this) result(x)
            class(equation_solver), intent(in) :: this
            logical :: x
        end function

        module subroutine es_set_print_status(this, x)
            class(equation_solver), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine
    end interface

! ******************************************************************************
! NONLIN_EQUATION_SOLVER_1VAR.F90
! ------------------------------------------------------------------------------
    !> @brief A base class for various solvers of equations of one variable.
    !!
    !! @par Example
    !! The following example illustrates the solution of an equation of one
    !! variable by means of Brent's method.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use nonlin_core
    !!     use nonlin_solve
    !!     implicit none
    !!
    !!     ! Local variables
    !!     type(fcn1var_helper) :: obj
    !!     procedure(fcn1var), pointer :: fcn
    !!     type(brent_solver) :: solver
    !!     real(real64) :: x, f
    !!     type(value_pair) :: limits
    !!
    !!     ! Define the search limits
    !!     limits%x1 = 1.5d0
    !!     limits%x2 = 5.0d0
    !!
    !!     ! Establish the function
    !!     fcn => fcn1
    !!     call obj%set_fcn(fcn)
    !!
    !!     ! Solve the equation
    !!     call solver%solve(obj, x, limits, f)
    !!
    !!     ! Print the output and the residual
    !!     print '(AF7.5)', "The solution: ", x
    !!     print '(AE10.3)', "The residual: ", f
    !!
    !! contains
    !!     ! The function:
    !!     ! f(x) = sin(x) / x, solution: x = x * pi for n = 1, 2, 3, ...
    !!     function fcn1(x) result(f)
    !!         real(real64), intent(in) :: x
    !!         real(real64) :: f
    !!         f = sin(x) / x
    !!     end function
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! The solution: 3.15159
    !! The residual: -0.751E-11
    !! @endcode
    type, abstract :: equation_solver_1var
        private
        !> The maximum number of function evaluations allowed per solve.
        integer(int32) :: m_maxEval = 100
        !> The convergence criteria on function value.
        real(real64) :: m_fcnTol = 1.0d-8
        !> The convergence criteria on change in variable value.
        real(real64) :: m_xtol = 1.0d-12
        !> The convergence criteria on the slope of the function (derivative)
        real(real64) :: m_difftol = 1.0d-12
        !> Set to true to print iteration status; else, false.
        logical :: m_printStatus = .false.
    contains
        !> @brief Gets the maximum number of function evaluations allowed during
        !! a single solve.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_max_fcn_evals(class(equation_solver_1var) this)
        !! @endcode
        !!
        !! @param[in] this The equation_solver_1var object.
        !! @return The maximum number of function evaluations.
        procedure, public :: get_max_fcn_evals => es1_get_max_eval
        !> @brief Sets the maximum number of function evaluations allowed during
        !! a single solve.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_max_fcn_evals(class(equation_solver_1var) this, integer(int32) n)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver_1var object.
        !! @param[in] n The maximum number of function evaluations.
        procedure, public :: set_max_fcn_evals => es1_set_max_eval
        !> @brief Gets the convergence on function value tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_fcn_tolerance(class(equation_solver_1var) this)
        !! @endcode
        !!
        !! @param[in] this The equation_solver_1var object.
        !! @return The tolerance value.
        procedure, public :: get_fcn_tolerance => es1_get_fcn_tol
        !> @brief Sets the convergence on function value tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_fcn_tolerance(class(equation_solver_1var) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver_1var object.
        !! @param[in] x The tolerance value.
        procedure, public :: set_fcn_tolerance => es1_set_fcn_tol
        !> @brief Gets the convergence on change in variable tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_var_tolerance(class(equation_solver_1var) this)
        !! @endcode
        !!
        !! @param[in] this The equation_solver_1var object.
        !! @return The tolerance value.
        procedure, public :: get_var_tolerance => es1_get_var_tol
        !> @brief Sets the convergence on change in variable tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_var_tolerance(class(equation_solver_1var) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver_1var object.
        !! @param[in] x The tolerance value.
        procedure, public :: set_var_tolerance => es1_set_var_tol
        !> @brief Gets a logical value determining if iteration status should be
        !! printed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function get_print_status(class(equation_solver_1var) this)
        !! @endcode
        !!
        !! @param[in] this The equation_solver_1var object.
        !! @return True if the iteration status should be printed; else, false.
        procedure, public :: get_print_status => es1_get_print_status
        !> @brief Sets a logical value determining if iteration status should be
        !! printed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_print_status(class(equation_solver_1var) this, logical x)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver_1var object.
        !! @param[in] x True if the iteration status should be printed; else, false.
        procedure, public :: set_print_status => es1_set_print_status
        !> @brief Solves the equation.
        !!
        !! @par Example
        !! See the equation_solver_1var type for an example on how to solve a
        !! system of equations.
        procedure(nonlin_solver_1var), deferred, public, pass :: solve
        !> @brief Gets the convergence on slope of the function (derivative)
        !! tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_diff_tolerance(class(equation_solver_1var) this)
        !! @endcode
        !!
        !! @param[in] this The equation_solver object.
        !! @return The tolerance value.
        procedure, public :: get_diff_tolerance => es1_get_diff_tol
        !> @brief Sets the convergence on slope of the function (derivative)
        !! tolerance.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_diff_tolerance(class(equation_solver_1var) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The equation_solver object.
        !! @param[in] x The tolerance value.
        !!
        !! @par Example
        !! See the equation_solver type for an example on how to establish
        !! the gradient vector slope tolerance.
        procedure, public :: set_diff_tolerance => es1_set_diff_tol
    end type

! ------------------------------------------------------------------------------
    interface
        pure module function es1_get_max_eval(this) result(n)
            class(equation_solver_1var), intent(in) :: this
            integer(int32) :: n
        end function

        module subroutine es1_set_max_eval(this, n)
            class(equation_solver_1var), intent(inout) :: this
            integer(int32), intent(in) :: n
        end subroutine

        pure module function es1_get_fcn_tol(this) result(x)
            class(equation_solver_1var), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine es1_set_fcn_tol(this, x)
            class(equation_solver_1var), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function es1_get_var_tol(this) result(x)
            class(equation_solver_1var), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine es1_set_var_tol(this, x)
            class(equation_solver_1var), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function es1_get_print_status(this) result(x)
            class(equation_solver_1var), intent(in) :: this
            logical :: x
        end function

        module subroutine es1_set_print_status(this, x)
            class(equation_solver_1var), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine

        pure module function es1_get_diff_tol(this) result(x)
            class(equation_solver_1var), intent(in) :: this
            real(real64) :: x
        end function

    ! --------------------
        module subroutine es1_set_diff_tol(this, x)
            class(equation_solver_1var), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine
    end interface

! ******************************************************************************
! NONLIN_EQUATION_OPTIMIZER.F90
! ------------------------------------------------------------------------------
    !> @brief A base class for optimization of an equation of multiple
    !! variables.
    !!
    !! @par Example
    !! The following example illustrates how to find the minimum of a function
    !! of multipler variables by means of the BFGS method.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use nonlin_core
    !!     use nonlin_optimize
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     type(bfgs) :: solver
    !!     type(fcnnvar_helper) :: obj
    !!     procedure(fcnnvar), pointer :: fcn
    !!     procedure(gradientfcn), pointer :: grad
    !!     real(real64) :: x(2), f
    !!
    !!     ! Tell the solver where to find the function
    !!     fcn => beale
    !!     grad => bealegrad
    !!     call obj%set_fcn(fcn, 2)
    !!     call obj%set_gradient_fcn(grad)
    !!
    !!     ! Define an initial guess
    !!     x = 1.0d0
    !!
    !!     ! Compute the solution
    !!     call solver%solve(obj, x, f)
    !!
    !!     ! Display the output
    !!     print '(AF7.5AF7.5A)', "Minimum: (", x(1), ", ", x(2), ")"
    !!     print '(AE9.3)', "Function Value: ", f
    !! contains
    !!     ! The Beale function:
    !!     ! f(x) = (1.5 - x + xy)**2 + (2.25 - x + xy**2)**2 + (2.625 - x + xy**3)**2
    !!     ! The minimum is at x = 3, y = 0.5, and f(3, 0.5) = 0
    !!     function beale(x) result(f)
    !!         real(real64), intent(in), dimension(:) :: x
    !!         real(real64) :: f
    !!         f = (1.5d0 - x(1) + x(1) * x(2))**2 + &
    !!             (2.25d0 - x(1) + x(1) * x(2)**2)**2 + &
    !!             (2.625d0 - x(1) + x(1) * x(2)**3)**2
    !!     end function
    !!
    !!     ! The gradient
    !!     subroutine bealegrad(x, g)
    !!         real(real64), intent(in), dimension(:) :: x
    !!         real(real64), intent(out), dimension(:) :: g
    !!
    !!         g(1) = 2.0d0 * (x(2)**3 - 1.0d0) * (x(1) * x(2)**3 - x(1) + 2.625d0) + &
    !!             2.0d0 * (x(2)**2 - 1.0d0) * (x(1) * x(2)**2 - x(1) + 2.25d0) + &
    !!             2.0d0 * (x(2) - 1.0d0) * (x(1) * x(2) - x(1) + 1.5d0)
    !!
    !!         g(2) = 6.0d0 * x(1) * x(2)**2 * (x(1) * x(2)**3 - x(1) + 2.625d0) + &
    !!             4.0d0 * x(1) * x(2) * (x(1) * x(2)**2 - x(1) + 2.25d0) + &
    !!             2.0d0 * x(1) * (x(1) * x(2) - x(1) + 1.5d0)
    !!     end subroutine
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! Minimum: (3.00000, 0.50000)
    !! Function Value: 0.999E-28
    !! @endcode
    type, abstract :: equation_optimizer
        private
        !> The maximum number of function evaluations allowed.
        integer(int32) :: m_maxEval = 500
        !> The error tolerance used to determine convergence.
        real(real64) :: m_tol = 1.0d-12
        !> Set to true to print iteration status; else, false.
        logical :: m_printStatus = .false.
    contains
        !> @brief Gets the maximum number of function evaluations allowed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_max_fcn_evals(class(equation_optimizer) this)
        !! @endcode
        !!
        !! @param[in] this The equation_optimizer object.
        !! @return The maximum number of function evaluations.
        procedure, public :: get_max_fcn_evals => oe_get_max_eval
        !> @brief Sets the maximum number of function evaluations allowed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_max_fcn_evals(class(equation_optimizer) this, integer(int32) n)
        !! @endcode
        !!
        !! @param[in,out] this The equation_optimizer object.
        !! @param[in] n The maximum number of function evaluations.
        procedure, public :: set_max_fcn_evals => oe_set_max_eval
        !> @brief Gets the tolerance on convergence.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_tolerance(class(equation_optimizer) this)
        !! @endcode
        !!
        !! @param[in] this The equation_optimizer object.
        !! @return The convergence tolerance.
        procedure, public :: get_tolerance => oe_get_tol
        !> @brief Sets the tolerance on convergence.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_tolerance(class(equation_optimizer) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The equation_optimizer object.
        !! @param[in] x The convergence tolerance.
        procedure, public :: set_tolerance => oe_set_tol
        !> @brief Gets a logical value determining if iteration status should be
        !! printed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical get_print_status(class(equation_optimizer) this)
        !! @endcode
        !!
        !! @param[in] this The equation_optimizer object.
        !! @return True if the iteration status should be printed; else, false.
        procedure, public :: get_print_status => oe_get_print_status
        !> @brief Sets a logical value determining if iteration status should be
        !! printed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_print_status(class(equation_optimizer) this, logical x)
        !! @endcode
        !!
        !! @param[in,out] this The equation_optimizer object.
        !! @param[in] x True if the iteration status should be printed; else, false.
        procedure, public :: set_print_status => oe_set_print_status
        !> @brief Optimizes the equation.
        !!
        !! @par Example
        !! See the equation_optimizer type for an example on how to solve a
        !! system of equations.
        procedure(nonlin_optimize_fcn), deferred, public, pass :: solve
    end type

! ------------------------------------------------------------------------------
    interface
        pure module function oe_get_max_eval(this) result(n)
            class(equation_optimizer), intent(in) :: this
            integer(int32) :: n
        end function

        module subroutine oe_set_max_eval(this, n)
            class(equation_optimizer), intent(inout) :: this
            integer(int32), intent(in) :: n
        end subroutine

        pure module function oe_get_tol(this) result(x)
            class(equation_optimizer), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine oe_set_tol(this, x)
            class(equation_optimizer), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function oe_get_print_status(this) result(x)
            class(equation_optimizer), intent(in) :: this
            logical :: x
        end function

        module subroutine oe_set_print_status(this, x)
            class(equation_optimizer), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine
    end interface

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
            use, intrinsic :: iso_fortran_env, only : real64
            use ferror, only : errors
            import equation_solver
            import vecfcn_helper
            import iteration_behavior
            class(equation_solver), intent(inout) :: this
            class(vecfcn_helper), intent(in) :: fcn
            real(real64), intent(inout), dimension(:) :: x
            real(real64), intent(out), dimension(:) :: fvec
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
            use, intrinsic :: iso_fortran_env, only : real64
            use ferror, only : errors
            import equation_solver_1var
            import fcn1var_helper
            import value_pair
            import iteration_behavior
            class(equation_solver_1var), intent(inout) :: this
            class(fcn1var_helper), intent(in) :: fcn
            real(real64), intent(inout) :: x
            type(value_pair), intent(in) :: lim
            real(real64), intent(out), optional :: f
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
        subroutine nonlin_optimize_fcn(this, fcn, x, fout, ib, err)
            use, intrinsic :: iso_fortran_env, only : real64
            use ferror, only : errors
            import equation_optimizer
            import fcnnvar_helper
            import iteration_behavior
            class(equation_optimizer), intent(inout) :: this
            class(fcnnvar_helper), intent(in) :: fcn
            real(real64), intent(inout), dimension(:) :: x
            real(real64), intent(out), optional :: fout
            type(iteration_behavior), optional :: ib
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

contains
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
        integer(int32), intent(in) :: iter, nfeval, njaceval
        real(real64), intent(in) :: xnorm, fnorm

        ! Process
        print *, ""
        print '(AI0)', "Iteration: ", iter
        print '(AI0)', "Function Evaluations: ", nfeval
        if (njaceval > 0) print '(AI0)', "Jacobian Evaluations: ", njaceval
        print '(AE9.3)', "Change in Variable: ", xnorm
        print '(AE9.3)', "Residual: ", fnorm
    end subroutine
end module
