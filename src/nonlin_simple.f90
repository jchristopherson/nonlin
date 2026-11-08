module nonlin_simple
    use iso_fortran_env
    use nonlin_single_var, only : fcn1var_helper, fcn1var
    use nonlin_multi_eqn_mult_var, only : vecfcn_helper
    use nonlin_multi_var
    use nonlin_types
    use nonlin_solve
    use autodiff, only : scalar_function, vector_function
    implicit none
    private
    public :: solve_brent
    public :: solve_newton_1var
    public :: solve_newton
    public :: solve_quasi_newton

contains
! ------------------------------------------------------------------------------
function solve_brent(fcn, bracket, args) result(x)
    !! Utilizes Brent's method to solve a scalar equation of the form 
    !! \(f(x) = 0\).  For more control over the solver it is recommended to
    !! utilize the [[brent_solver]] type.  This function is a convenience
    !! wrapper around this solver.
    procedure(fcn1var), intent(in), pointer :: fcn
        !! The function.
    type(value_pair), intent(in) :: bracket
        !! An initial estimate that brackets the solution.
    class(*), intent(in), optional :: args
        !! An optional object that can be used to pass information to fcn.
    real(real64) :: x
        !! The value of \(x\) that satisfies \(f(x) = 0\) within the default
        !! solver tolerances.

    type(fcn1var_helper) :: eqn
    type(brent_solver) :: solver
    call eqn%set_fcn(fcn)
    call solver%solve(eqn, x, bracket, args = args)
end function

! ------------------------------------------------------------------------------
function solve_newton_1var(fcn, bracket, args) result(x)
    !! Utilizes Newton's method to solve a scalar equation of the form 
    !! \(f(x) = 0\).  For more control over the solver it is recommended to
    !! utilize the [[newton_1var_solver]] type.  This function is a convenience
    !! wrapper around this solver.
    procedure(fcn1var), intent(in), pointer :: fcn
        !! The function.
    type(value_pair), intent(in) :: bracket
        !! An initial estimate that brackets the solution.
    class(*), intent(in), optional :: args
        !! An optional object that can be used to pass information to fcn.
    real(real64) :: x
        !! The value of \(x\) that satisfies \(f(x) = 0\) within the default
        !! solver tolerances.

    type(fcn1var_helper) :: eqn
    type(newton_1var_solver) :: solver
    call eqn%set_fcn(fcn)
    call solver%solve(eqn, x, bracket, args = args)
end function

! ------------------------------------------------------------------------------
function solve_newton(fcn, xo, args) result(x)
    !! Utilizes Newton's method to solve a system of N equations and N unknowns.
    !! For more control over the solver it is recommended to utilize the 
    !! [[newton_solver]] type.  This function is a convenience wrapper around 
    !! this solver.
    procedure(vector_function), intent(in), pointer :: fcn
        !! The function.
    real(real64), intent(in), dimension(:) :: xo
        !! An N-element array containing an initial guess at the solution.
    class(*), intent(in), optional :: args
        !! An optional object that can be used to pass information to fcn.
    real(real64), allocatable, dimension(:) :: x
        !! An N-element array containing the values of \(x\) that satisfies 
        !! \(f(x) = 0\) within the default solver tolerances.

    integer :: nvar, neqn
    type(vecfcn_helper) :: eqns
    type(newton_solver) :: solver
    real(real64), allocatable, dimension(:) :: f

    nvar = size(x)
    neqn = nvar
    allocate(x(nvar), source = xo)
    allocate(f(neqn))
    call eqns%set_auto_fcn(fcn, neqn, nvar)
    call solver%solve(eqns, x, f, args = args)
end function

! ------------------------------------------------------------------------------
function solve_quasi_newton(fcn, xo, args) result(x)
    !! Utilizes a Quasi-Newton method to solve a system of N equations and N 
    !! unknowns.  For more control over the solver it is recommended to utilize
    !! the [[quasi_newton_solver]] type.  This function is a convenience wrapper 
    !! around this solver.
    procedure(vector_function), intent(in), pointer :: fcn
        !! The function.
    real(real64), intent(in), dimension(:) :: xo
        !! An N-element array containing an initial guess at the solution.
    class(*), intent(in), optional :: args
        !! An optional object that can be used to pass information to fcn.
    real(real64), allocatable, dimension(:) :: x
        !! An N-element array containing the values of \(x\) that satisfies 
        !! \(f(x) = 0\) within the default solver tolerances.

    integer :: nvar, neqn
    type(vecfcn_helper) :: eqns
    type(quasi_newton_solver) :: solver
    real(real64), allocatable, dimension(:) :: f

    nvar = size(x)
    neqn = nvar
    allocate(x(nvar), source = xo)
    allocate(f(neqn))
    call eqns%set_auto_fcn(fcn, neqn, nvar)
    call solver%solve(eqns, x, f, args = args)
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
end module