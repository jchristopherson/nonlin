# nonlin
A library that provides routines to compute the solutions to systems of nonlinear equations.

## Example 1
This example solves a set of two equations of two unknowns using a Quasi-Newton type solver.  In this example, the solver is left to compute the derivatives numerically.

```fortran

    ! Define the routine containing the eqations to solve.  The equations are:
    ! x**2 + y**2 = 34
    ! x**2 - 2 * y**2 = 7
    subroutine fcns(x, f)
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: f
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    end subroutine

    ! ...

    ! Solver Code
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    type(quasi_newton_solver) :: solver
    real(dp) :: x(2), f(2)

    ! Define an initial guess
    x = [1.0d0, 1.0d0]

    ! Solve
    call solver%solve(obj, x, f)
    
```