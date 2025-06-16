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
        print 100, "Iteration: ", iter
        print 100, "Function Evaluations: ", nfeval
        if (njaceval > 0) print 100, "Jacobian Evaluations: ", njaceval
        print 101, "Change in Variable: ", xnorm
        print 101, "Residual: ", fnorm

        ! Formatting
100     format(A, I0)
101     format(A, E10.3)
    end subroutine
end module
