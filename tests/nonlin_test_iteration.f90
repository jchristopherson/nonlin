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
! nonlin_test_iteration.f90

! Tests that the iter_count, fcn_count, jacobian_count, and gradient_count
! fields of the iteration_behavior type agree with the number of times the
! solver actually evaluated the user-supplied function, Jacobian, derivative,
! or gradient.
module nonlin_test_iteration
    use iso_fortran_env
    use nonlin
    implicit none
    private
    public :: test_iteration_brent
    public :: test_iteration_newton_1var
    public :: test_iteration_newton
    public :: test_iteration_quasi_newton
    public :: test_iteration_least_squares
    public :: test_iteration_constrained_least_squares
    public :: test_iteration_nelder_mead
    public :: test_iteration_bfgs

    ! Number of times each user routine has been called since the last reset
    integer, save :: n_fcn = 0
    integer, save :: n_diff = 0
    integer, save :: n_jac = 0
    integer, save :: n_grad = 0
contains
! ******************************************************************************
! USER ROUTINES
! ------------------------------------------------------------------------------
    ! f(x) = x**3 - 2, SOLUTION: x = 2**(1/3)
    function f_cubic(x, args) result(f)
        real(real64), intent(in) :: x
        class(*), intent(in), optional :: args
        real(real64) :: f
        n_fcn = n_fcn + 1
        f = x**3 - 2.0d0
    end function

    ! f'(x) = 3x**2
    function df_cubic(x, args) result(df)
        real(real64), intent(in) :: x
        class(*), intent(in), optional :: args
        real(real64) :: df
        n_diff = n_diff + 1
        df = 3.0d0 * x**2
    end function

    ! Two-equation system, SOLUTION: (5, 3)
    subroutine f_system(x, f, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        n_fcn = n_fcn + 1
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - 2.0d0 * x(2)**2 - 7.0d0
    end subroutine

    ! Jacobian of f_system
    subroutine j_system(x, j, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: j
        class(*), intent(in), optional :: args
        n_jac = n_jac + 1
        j = 2.0d0 * reshape([x(1), x(1), x(2), -2.0d0 * x(2)], [2, 2])
    end subroutine

    ! Rosenbrock function, SOLUTION: (1, 1)
    function rosenbrock(x, args) result(f)
        real(real64), intent(in), dimension(:) :: x
        class(*), intent(in), optional :: args
        real(real64) :: f
        n_fcn = n_fcn + 1
        f = 1.0d2 * (x(2) - x(1)**2)**2 + (x(1) - 1.0d0)**2
    end function

    ! Gradient of the Rosenbrock function
    subroutine rosenbrock_grad(x, g, args)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: g
        class(*), intent(in), optional :: args
        n_grad = n_grad + 1
        g(1) = -4.0d2 * x(1) * (x(2) - x(1)**2) - 2.0d0 * (1.0d0 - x(1))
        g(2) = 2.0d2 * (x(2) - x(1)**2)
    end subroutine

! ******************************************************************************
! HELPER ROUTINES
! ------------------------------------------------------------------------------
    ! Resets all call counters.
    subroutine reset_counts()
        n_fcn = 0
        n_diff = 0
        n_jac = 0
        n_grad = 0
    end subroutine

    ! Checks the iteration_behavior fields against the number of calls made to
    ! the user routines.
    function counts_ok(ib, nfcn, njac, ngrad) result(ok)
        type(iteration_behavior), intent(in) :: ib
        integer, intent(in) :: nfcn, njac, ngrad
        logical :: ok
        ok = ib%fcn_count == nfcn .and. ib%jacobian_count == njac .and. &
            ib%gradient_count == ngrad .and. ib%iter_count >= 1
    end function

    subroutine report(name, ib)
        character(len = *), intent(in) :: name
        type(iteration_behavior), intent(in) :: ib
        print *, "Iteration count check failed: ", name
        print *, "  ib iter/fcn/jac/grad: ", ib%iter_count, ib%fcn_count, &
            ib%jacobian_count, ib%gradient_count
        print *, "  calls fcn/diff/jac/grad: ", n_fcn, n_diff, n_jac, n_grad
    end subroutine

! ******************************************************************************
! TESTS
! ------------------------------------------------------------------------------
    function test_iteration_brent() result(check)
        ! Local Variables
        type(brent_solver) :: solver
        type(fcn1var_helper) :: obj
        procedure(fcn1var), pointer :: fcn
        type(iteration_behavior) :: ib
        type(value_pair) :: limits
        real(real64) :: x, f
        logical :: check

        ! Parameters
        real(real64), parameter :: root = 1.2599210498948732d0
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        fcn => f_cubic
        call obj%set_fcn(fcn)
        limits%x1 = 0.0d0
        limits%x2 = 2.0d0
        call reset_counts()

        ! Solve
        call solver%solve(obj, x, limits, f, ib)

        ! Test
        check = abs(x - root) < tol .and. counts_ok(ib, n_fcn, n_jac, n_grad)
        if (.not.check) then
            call report("Brent", ib)
            print *, "  root: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_iteration_newton_1var() result(check)
        ! Local Variables
        type(newton_1var_solver) :: solver
        type(fcn1var_helper) :: obj
        procedure(fcn1var), pointer :: fcn, diff
        type(iteration_behavior) :: ib
        type(value_pair) :: limits
        real(real64) :: x, f
        logical :: check

        ! Parameters
        real(real64), parameter :: root = 1.2599210498948732d0
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        fcn => f_cubic
        diff => df_cubic
        call obj%set_fcn(fcn)
        call obj%set_diff(diff)
        limits%x1 = 0.0d0
        limits%x2 = 2.0d0
        call reset_counts()

        ! Solve
        call solver%solve(obj, x, limits, f, ib)

        ! Test - the derivative routine is the Jacobian for a 1-variable problem
        check = abs(x - root) < tol .and. counts_ok(ib, n_fcn, n_diff, n_grad)
        if (.not.check) then
            call report("Newton (1 variable)", ib)
            print *, "  root: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_iteration_newton() result(check)
        ! Local Variables
        type(newton_solver) :: solver
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2)
        logical :: check, good
        integer :: i

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        fcn => f_system
        jac => j_system
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)

        ! Solve with and without a line search
        do i = 1, 2
            call solver%set_use_line_search(i == 1)
            x = [4.0d0, 2.0d0]
            call reset_counts()
            call solver%solve(obj, x, f, ib)
            good = maxval(abs(x - [5.0d0, 3.0d0])) < tol
            good = good .and. counts_ok(ib, n_fcn, n_jac, n_grad)
            if (.not.good) then
                check = .false.
                call report("Newton, line search = " // merge("on ", "off", i == 1), ib)
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    function test_iteration_quasi_newton() result(check)
        ! Local Variables
        type(quasi_newton_solver) :: solver
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2)
        logical :: check, good
        integer :: i

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        fcn => f_system
        jac => j_system
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)

        ! Solve with and without a line search
        do i = 1, 2
            call solver%set_use_line_search(i == 1)
            x = [4.0d0, 2.0d0]
            call reset_counts()
            call solver%solve(obj, x, f, ib)
            good = maxval(abs(x - [5.0d0, 3.0d0])) < tol
            good = good .and. counts_ok(ib, n_fcn, n_jac, n_grad)
            if (.not.good) then
                check = .false.
                call report("Quasi-Newton, line search = " // merge("on ", "off", i == 1), ib)
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    function test_iteration_least_squares() result(check)
        ! Local Variables
        type(least_squares_solver) :: solver
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2)
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        fcn => f_system
        jac => j_system
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)
        x = [4.0d0, 2.0d0]
        call reset_counts()

        ! Solve
        call solver%solve(obj, x, f, ib)

        ! Test
        check = maxval(abs(x - [5.0d0, 3.0d0])) < tol
        check = check .and. counts_ok(ib, n_fcn, n_jac, n_grad)
        if (.not.check) then
            call report("Least squares", ib)
            print *, "  solution: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_iteration_constrained_least_squares() result(check)
        ! Local Variables
        type(constrained_least_squares_solver) :: solver
        type(vecfcn_helper) :: obj
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jac
        type(iteration_behavior) :: ib
        real(real64) :: x(2), f(2)
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        fcn => f_system
        jac => j_system
        call obj%set_fcn(fcn, 2, 2)
        call obj%set_jacobian(jac)
        call solver%set_lower_limits([0.0d0, 0.0d0])
        call solver%set_upper_limits([10.0d0, 10.0d0])
        x = [4.0d0, 2.0d0]
        call reset_counts()

        ! Solve
        call solver%solve(obj, x, f, ib)

        ! Test
        check = maxval(abs(x - [5.0d0, 3.0d0])) < tol
        check = check .and. counts_ok(ib, n_fcn, n_jac, n_grad)
        if (.not.check) then
            call report("Constrained least squares", ib)
            print *, "  solution: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_iteration_nelder_mead() result(check)
        ! Local Variables
        type(nelder_mead) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        type(iteration_behavior) :: ib
        real(real64) :: x(2), fout
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-3

        ! Initialization
        fcn => rosenbrock
        call obj%set_fcn(fcn, 2)
        x = [-1.2d0, 1.0d0]
        call reset_counts()

        ! Solve
        call solver%solve(obj, x, fout, ib)

        ! Test - Nelder-Mead uses no derivative information
        check = maxval(abs(x - [1.0d0, 1.0d0])) < tol
        check = check .and. counts_ok(ib, n_fcn, 0, 0)
        if (.not.check) then
            call report("Nelder-Mead", ib)
            print *, "  solution: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_iteration_bfgs() result(check)
        ! Local Variables
        type(bfgs) :: solver
        type(fcnnvar_helper) :: obj
        procedure(fcnnvar), pointer :: fcn
        procedure(gradientfcn), pointer :: grad
        type(iteration_behavior) :: ib
        real(real64) :: x(2), fout
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-5

        ! Initialization
        fcn => rosenbrock
        grad => rosenbrock_grad
        call obj%set_fcn(fcn, 2)
        call obj%set_gradient_fcn(grad)
        x = 0.0d0
        call reset_counts()

        ! Solve
        call solver%solve(obj, x, fout, ib)

        ! Test - BFGS uses gradient information but no Jacobian
        check = maxval(abs(x - [1.0d0, 1.0d0])) < tol
        check = check .and. counts_ok(ib, n_fcn, 0, n_grad)
        if (.not.check) then
            call report("BFGS", ib)
            print *, "  solution: ", x
        end if
    end function

! ------------------------------------------------------------------------------
end module
