module nonlin_simple
    use iso_fortran_env
    use nonlin_single_var, only : fcn1var_helper, fcn1var
    use nonlin_multi_eqn_mult_var
    use nonlin_multi_var
    use nonlin_types
    use nonlin_solve
    use autodiff
    implicit none
    private
    public :: solve_brent

contains
! ------------------------------------------------------------------------------
function solve_brent(fcn, bracket) result(x)
    !! Utilizes Brent's method to solve a scalar equation of the form 
    !! \(f(x) = 0\).  For more control over the solver it is recommended to
    !! utilize the [[brent_solver]] type.  This function is a convenience
    !! wrapper around this solver.
    procedure(fcn1var) :: fcn
        !! The function.
    type(value_pair), intent(in) :: bracket
        !! An initial estimate that brackets the solution.
    real(real64) :: x
        !! The value of \(x\) that satisfies \(f(x) = 0\) within the default
        !! solver tolerances.

    type(fcn1var_helper) :: eqn
    type(brent_solver) :: solver
    call eqn%set_fcn(fcn)
    call solver%solve(eqn, x, bracket)
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module