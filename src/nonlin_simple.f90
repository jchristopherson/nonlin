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
! nonlin_simple.f90

module nonlin_simple
    use iso_fortran_env
    use nonlin_single_var, only : fcn1var_helper, fcn1var
    use nonlin_multi_eqn_mult_var, only : vecfcn_helper
    use nonlin_multi_var
    use nonlin_types
    use nonlin_solve
    use nonlin_optimize
    use nonlin_continuation
    use autodiff, only : scalar_function, vector_function
    implicit none
    private
    public :: solve_brent
    public :: solve_newton_1var
    public :: solve_newton
    public :: solve_quasi_newton
    public :: solve_nelder_mead
    public :: solve_bfgs
    public :: solve_continuation
    public :: solve_continuation_pseudo_arclength

contains
! ------------------------------------------------------------------------------
function solve_brent(fcn, bracket, args) result(x)
    !! Utilizes Brent's method to solve a scalar equation of the form 
    !! \(f(x) = 0\).  For more control over the solver it is recommended to
    !! utilize the [[brent_solver]] type.  This function is a convenience
    !! wrapper around this solver.
    procedure(fcn1var) :: fcn
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
    procedure(fcn1var), pointer :: fptr

    fptr => fcn
    call eqn%set_fcn(fptr)
    call solver%solve(eqn, x, bracket, args = args)

end function

! ------------------------------------------------------------------------------
function solve_newton_1var(fcn, bracket, args) result(x)
    !! Utilizes Newton's method to solve a scalar equation of the form 
    !! \(f(x) = 0\).  For more control over the solver it is recommended to
    !! utilize the [[newton_1var_solver]] type.  This function is a convenience
    !! wrapper around this solver.
    procedure(fcn1var) :: fcn
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
    procedure(fcn1var), pointer :: fptr

    fptr => fcn
    call eqn%set_fcn(fptr)
    call solver%solve(eqn, x, bracket, args = args)

end function

! ------------------------------------------------------------------------------
function solve_newton(fcn, xo, args) result(x)
    !! Utilizes Newton's method to solve a system of N equations and N unknowns.
    !! For more control over the solver it is recommended to utilize the 
    !! [[newton_solver]] type.  This function is a convenience wrapper around 
    !! this solver.
    procedure(vector_function) :: fcn
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
    procedure(vector_function), pointer :: fptr
    real(real64), allocatable, dimension(:) :: f

    nvar = size(xo)
    neqn = nvar
    allocate(x(nvar), source = xo)
    allocate(f(neqn))
    fptr => fcn
    call eqns%set_auto_fcn(fptr, neqn, nvar)
    call solver%solve(eqns, x, f, args = args)

end function

! ------------------------------------------------------------------------------
function solve_quasi_newton(fcn, xo, args) result(x)
    !! Utilizes a Quasi-Newton method to solve a system of N equations and N 
    !! unknowns.  For more control over the solver it is recommended to utilize
    !! the [[quasi_newton_solver]] type.  This function is a convenience wrapper 
    !! around this solver.
    procedure(vector_function) :: fcn
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
    procedure(vector_function), pointer :: fptr
    real(real64), allocatable, dimension(:) :: f

    nvar = size(xo)
    neqn = nvar
    allocate(x(nvar), source = xo)
    allocate(f(neqn))
    fptr => fcn
    call eqns%set_auto_fcn(fptr, neqn, nvar)
    call solver%solve(eqns, x, f, args = args) 
end function

! ------------------------------------------------------------------------------
function solve_nelder_mead(fcn, xo, args) result(x)
    !! Utilizes the Nelder-Mead method for minimizing a function of N variables.
    !! For more control over the solver it is recommended to utilize the 
    !! [[nelder_mead]] type.  This function is a convenience wrapper around 
    !! this solver.
    procedure(scalar_function) :: fcn
        !! The function.
    real(real64), intent(in), dimension(:) :: xo
        !! An N-element array containing an initial guess at the solution.
    class(*), intent(in), optional :: args
        !! An optional object that can be used to pass information to fcn.
    real(real64), allocatable, dimension(:) :: x
        !! An N-element array containing the values of \(x\) that minimize 
        !! \(f(x)\) within the default solver tolerances.

    integer :: nvar
    procedure(scalar_function), pointer :: fptr
    type(fcnnvar_helper) :: eqn
    type(nelder_mead) :: solver

    nvar = size(xo)
    allocate(x(nvar), source = xo)
    fptr => fcn
    call eqn%set_auto_fcn(fptr, nvar)
    call solver%solve(eqn, x, args = args)

end function

! ------------------------------------------------------------------------------
function solve_bfgs(fcn, xo, args) result(x)
    !! Utilizes the BFGS method for minimizing a function of N variables.
    !! For more control over the solver it is recommended to utilize the 
    !! [[bfgs]] type.  This function is a convenience wrapper around 
    !! this solver.
    procedure(scalar_function) :: fcn
        !! The function.
    real(real64), intent(in), dimension(:) :: xo
        !! An N-element array containing an initial guess at the solution.
    class(*), intent(in), optional :: args
        !! An optional object that can be used to pass information to fcn.
    real(real64), allocatable, dimension(:) :: x
        !! An N-element array containing the values of \(x\) that minimize 
        !! \(f(x)\) within the default solver tolerances.

    integer :: nvar
    procedure(scalar_function), pointer :: fptr
    type(fcnnvar_helper) :: eqn
    type(bfgs) :: solver

    nvar = size(xo)
    allocate(x(nvar), source = xo)
    fptr => fcn
    call eqn%set_auto_fcn(fptr, nvar)
    call solver%solve(eqn, x, args = args)

end function

! ------------------------------------------------------------------------------
function solve_continuation(fcn, xo, args) result(x)
    !! Utilizes a continuation solver to solve a system of N equations and N 
    !! unknowns.  For more control over the solver it is recommended to utilize
    !! the [[continuation_solver]] type.  This function is a convenience wrapper 
    !! around this solver.
    procedure(vector_function) :: fcn
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
    type(continuation_solver) :: solver
    procedure(vector_function), pointer :: fptr
    real(real64), allocatable, dimension(:) :: f

    nvar = size(xo)
    neqn = nvar
    allocate(x(nvar), source = xo)
    allocate(f(neqn))
    fptr => fcn
    call eqns%set_auto_fcn(fptr, neqn, nvar)
    call solver%solve(eqns, x, f, args = args) 
end function

! ------------------------------------------------------------------------------
function solve_continuation_pseudo_arclength(fcn, xo, args) result(x)
    !! Utilizes a pseudo-arclength continuation solver to solve a system of N 
    !! equations and N unknowns.  For more control over the solver it is 
    !! recommended to utilize the [[continuation_solver]] type.  This function 
    !! is a convenience wrapper around this solver.
    procedure(vector_function) :: fcn
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
    type(continuation_solver) :: solver
    procedure(vector_function), pointer :: fptr
    real(real64), allocatable, dimension(:) :: f

    nvar = size(xo)
    neqn = nvar
    allocate(x(nvar), source = xo)
    allocate(f(neqn))
    fptr => fcn
    call eqns%set_auto_fcn(fptr, neqn, nvar)
    call solver%set_use_pseudo_arclength(.true.)
    call solver%solve(eqns, x, f, args = args) 
end function

! ------------------------------------------------------------------------------
end module