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
! nonlin_test_simple.f90

! Tests for the convenience solver functions provided by the nonlin_simple module.
module nonlin_test_simple
    use iso_fortran_env
    use nonlin
    use nonlin_simple
    use autodiff, only : dual, vector_function, assignment(=), &
        scalar_function, operator(+), operator(-), operator(*), operator(**)
    implicit none
    private
    public :: test_solve_brent_1
    public :: test_solve_brent_2
    public :: test_solve_newton_1var_1
    public :: test_solve_newton_1var_2
    public :: test_solve_newton_1
    public :: test_solve_newton_2
    public :: test_solve_newton_3
    public :: test_solve_quasi_newton_1
    public :: test_solve_quasi_newton_2
    public :: test_solve_quasi_newton_3
    public :: test_solve_continuation
    public :: test_solve_continuation_pseudo_arclength
    public :: test_solve_nelder_mead
    public :: test_solve_bfgs
contains
! ******************************************************************************
! HELPER FUNCTIONS
! ------------------------------------------------------------------------------
    ! f(x) = sin(x) / x, SOLUTION: x = n * pi for n = 1, 2, 3, ...
    function f_sin_over_x(x, args) result(f)
        real(real64), intent(in) :: x
        class(*), intent(in), optional :: args
        real(real64) :: f
        f = sin(x) / x
    end function

    ! f(x) = x**3 - a, SOLUTION: x = a**(1/3), with a supplied via args
    function f_cubic_shifted(x, args) result(f)
        real(real64), intent(in) :: x
        class(*), intent(in), optional :: args
        real(real64) :: f
        real(real64) :: a
        a = 0.0d0
        if (present(args)) then
            select type (args)
            type is (real(real64))
                a = args
            end select
        end if
        f = x**3 - a
    end function

    ! Solution: x = +/-5, y = +/-3 when a = 2
    ! x**2 + y**2 = 34
    ! x**2 - a * y**2 = 7
    subroutine ad_system_1(x, f, args)
        type(dual), intent(in), dimension(:) :: x
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        real(real64) :: a
        a = 2.0d0
        if (present(args)) then
            select type (args)
            type is (real(real64))
                a = args
            end select
        end if
        f(1) = x(1)**2 + x(2)**2 - 34.0d0
        f(2) = x(1)**2 - a * x(2)**2 - 7.0d0
    end subroutine

    ! Solution: x(i) = sqrt(i) for i = 1, ..., size(x)
    subroutine ad_diagonal_system(x, f, args)
        type(dual), intent(in), dimension(:) :: x
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        integer :: i
        do i = 1, size(x)
            f(i) = x(i)**2 - real(i, real64)
        end do
    end subroutine

    function ad_shifted_bowl(x, args) result(f)
        type(dual), intent(in), dimension(:) :: x
        class(*), intent(in), optional :: args
        type(dual) :: f
        real(real64) :: center
        center = 1.0d0
        if (present(args)) then
            select type (args)
            type is (real(real64))
                center = args
            end select
        end if
        f = (x(1) - center)**2 + (x(2) + 2.0d0)**2
    end function

    ! Returns .true. if every element of x is within tol of the corresponding
    ! element of ref.
    pure function all_close(x, ref, tol) result(c)
        real(real64), intent(in), dimension(:) :: x, ref
        real(real64), intent(in) :: tol
        logical :: c
        c = all(abs(x - ref) <= tol)
    end function

! ******************************************************************************
! BRENT'S METHOD
! ------------------------------------------------------------------------------
    function test_solve_brent_1() result(check)
        ! Local Variables
        type(value_pair) :: bracket
        real(real64) :: x
        logical :: check

        ! Parameters
        real(real64), parameter :: pi = 3.141592653589793d0
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        bracket%x1 = 1.5d0
        bracket%x2 = 5.0d0

        ! Compute the solution
        x = solve_brent(f_sin_over_x, bracket)

        ! The solution on this interval should be: pi
        if (abs(x - pi) > tol) then
            check = .false.
            print 100, &
                "Test Failed: solve_brent Test 1.  Expected: ", pi, &
                ", Found: ", x
        end if

        ! Formatting
100     format(A, F12.8, A, F12.8)
    end function

! ------------------------------------------------------------------------------
    function test_solve_brent_2() result(check)
        ! Local Variables
        type(value_pair) :: bracket
        real(real64) :: x, a
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        a = 8.0d0
        bracket%x1 = 0.0d0
        bracket%x2 = 5.0d0

        ! Compute the solution.  The root is only reachable if args is
        ! forwarded to the function.
        x = solve_brent(f_cubic_shifted, bracket, args = a)

        ! The solution should be: 2
        if (abs(x - 2.0d0) > tol) then
            check = .false.
            print 100, &
                "Test Failed: solve_brent Test 2.  Expected: ", 2.0d0, &
                ", Found: ", x
        end if

        ! Formatting
100     format(A, F12.8, A, F12.8)
    end function

! ******************************************************************************
! NEWTON'S METHOD (1 VARIABLE)
! ------------------------------------------------------------------------------
    function test_solve_newton_1var_1() result(check)
        ! Local Variables
        type(value_pair) :: bracket
        real(real64) :: x
        logical :: check

        ! Parameters
        real(real64), parameter :: pi = 3.141592653589793d0
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        bracket%x1 = 1.5d0
        bracket%x2 = 5.0d0

        ! Compute the solution
        x = solve_newton_1var(f_sin_over_x, bracket)

        ! The solution on this interval should be: pi
        if (abs(x - pi) > tol) then
            check = .false.
            print 100, &
                "Test Failed: solve_newton_1var Test 1.  Expected: ", pi, &
                ", Found: ", x
        end if

        ! Formatting
100     format(A, F12.8, A, F12.8)
    end function

! ------------------------------------------------------------------------------
    function test_solve_newton_1var_2() result(check)
        ! Local Variables
        type(value_pair) :: bracket
        real(real64) :: x, a
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        a = 8.0d0
        bracket%x1 = 0.0d0
        bracket%x2 = 5.0d0

        ! Compute the solution.  The root is only reachable if args is
        ! forwarded to the function.
        x = solve_newton_1var(f_cubic_shifted, bracket, args = a)

        ! The solution should be: 2
        if (abs(x - 2.0d0) > tol) then
            check = .false.
            print 100, &
                "Test Failed: solve_newton_1var Test 2.  Expected: ", 2.0d0, &
                ", Found: ", x
        end if

        ! Formatting
100     format(A, F12.8, A, F12.8)
    end function

! ******************************************************************************
! NEWTON'S METHOD (N VARIABLES)
! ------------------------------------------------------------------------------
    function test_solve_newton_1() result(check)
        ! Local Variables
        real(real64), dimension(2) :: x, xo, exact
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        xo = [4.0d0, 2.0d0]
        exact = [5.0d0, 3.0d0]

        ! Compute the solution
        x = solve_newton(ad_system_1, xo)

        ! Test the solution
        if (size(x) /= 2) then
            check = .false.
            print *, "Test Failed: solve_newton Test 1.  Wrong size: ", size(x)
            return
        end if
        if (.not.all_close(x, exact, tol)) then
            check = .false.
            print *, "Test Failed: solve_newton Test 1.  Expected: ", exact, &
                ", Found: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_solve_newton_2() result(check)
        ! Local Variables
        real(real64), dimension(2) :: x, xo, exact
        real(real64) :: a
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        a = 3.0d0
        xo = [5.0d0, 2.0d0]

        ! With a = 3, the solution is x**2 = 27.25 and y**2 = 6.75.  The
        ! solution with a = 2 (x = 5, y = 3) would indicate that args was
        ! not forwarded to the function.
        exact = [sqrt(27.25d0), sqrt(6.75d0)]

        ! Compute the solution
        x = solve_newton(ad_system_1, xo, args = a)

        ! Test the solution
        if (.not.all_close(x, exact, tol)) then
            check = .false.
            print *, "Test Failed: solve_newton Test 2.  Expected: ", exact, &
                ", Found: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_solve_newton_3() result(check)
        ! Local Variables
        real(real64), dimension(3) :: x, xo, exact
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        xo = [2.0d0, 2.0d0, 2.0d0]
        exact = sqrt([1.0d0, 2.0d0, 3.0d0])

        ! Compute the solution
        x = solve_newton(ad_diagonal_system, xo)

        ! Test the solution
        if (.not.all_close(x, exact, tol)) then
            check = .false.
            print *, "Test Failed: solve_newton Test 3.  Expected: ", exact, &
                ", Found: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_solve_continuation() result(check)
        real(real64), allocatable :: x(:)
        real(real64), dimension(2) :: xo, exact
        real(real64) :: a
        logical :: check

        a = 3.0d0
        xo = [5.0d0, 2.0d0]
        exact = [sqrt(27.25d0), sqrt(6.75d0)]
        x = solve_continuation(ad_system_1, xo, args = a)

        check = size(x) == size(exact) .and. all_close(x, exact, 1.0d-6)
        if (.not.check) print *, &
            "Test Failed: solve_continuation. Expected:", exact, "Found:", x
    end function

! ------------------------------------------------------------------------------
    function test_solve_continuation_pseudo_arclength() result(check)
        real(real64), allocatable :: x(:)
        real(real64), dimension(3) :: xo, exact
        logical :: check

        xo = [2.0d0, 2.0d0, 2.0d0]
        exact = sqrt([1.0d0, 2.0d0, 3.0d0])
        x = solve_continuation_pseudo_arclength(ad_diagonal_system, xo)

        check = size(x) == size(exact) .and. all_close(x, exact, 1.0d-6)
        if (.not.check) print *, &
            "Test Failed: solve_continuation_pseudo_arclength. Expected:", &
            exact, "Found:", x
    end function

! ******************************************************************************
! QUASI-NEWTON METHOD
! ------------------------------------------------------------------------------
    function test_solve_quasi_newton_1() result(check)
        ! Local Variables
        real(real64), dimension(2) :: x, xo, exact
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        xo = [4.0d0, 2.0d0]
        exact = [5.0d0, 3.0d0]

        ! Compute the solution
        x = solve_quasi_newton(ad_system_1, xo)

        ! Test the solution
        if (size(x) /= 2) then
            check = .false.
            print *, "Test Failed: solve_quasi_newton Test 1.  Wrong size: ", &
                size(x)
            return
        end if
        if (.not.all_close(x, exact, tol)) then
            check = .false.
            print *, "Test Failed: solve_quasi_newton Test 1.  Expected: ", &
                exact, ", Found: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_solve_quasi_newton_2() result(check)
        ! Local Variables
        real(real64), dimension(2) :: x, xo, exact
        real(real64) :: a
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        a = 3.0d0
        xo = [5.0d0, 2.0d0]

        ! With a = 3, the solution is x**2 = 27.25 and y**2 = 6.75.  The
        ! solution with a = 2 (x = 5, y = 3) would indicate that args was
        ! not forwarded to the function.
        exact = [sqrt(27.25d0), sqrt(6.75d0)]

        ! Compute the solution
        x = solve_quasi_newton(ad_system_1, xo, args = a)

        ! Test the solution
        if (.not.all_close(x, exact, tol)) then
            check = .false.
            print *, "Test Failed: solve_quasi_newton Test 2.  Expected: ", &
                exact, ", Found: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_solve_quasi_newton_3() result(check)
        ! Local Variables
        real(real64), dimension(3) :: x, xo, exact
        logical :: check

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Initialization
        check = .true.
        xo = [2.0d0, 2.0d0, 2.0d0]
        exact = sqrt([1.0d0, 2.0d0, 3.0d0])

        ! Compute the solution
        x = solve_quasi_newton(ad_diagonal_system, xo)

        ! Test the solution
        if (.not.all_close(x, exact, tol)) then
            check = .false.
            print *, "Test Failed: solve_quasi_newton Test 3.  Expected: ", &
                exact, ", Found: ", x
        end if
    end function

! ------------------------------------------------------------------------------
    function test_solve_nelder_mead() result(check)
        real(real64), parameter :: tol = 1.0d-5
        real(real64), dimension(2) :: xo, x, expected
        logical :: check

        xo = [4.0d0, 3.0d0]
        expected = [2.5d0, -2.0d0]
        x = solve_nelder_mead(ad_shifted_bowl, xo, args = expected(1))
        check = all_close(x, expected, tol)
        if (.not.check) print *, &
            "Test Failed: solve_nelder_mead. Expected:", expected, "Found:", x
    end function

! ------------------------------------------------------------------------------
    function test_solve_bfgs() result(check)
        real(real64), parameter :: tol = 1.0d-6
        real(real64), dimension(2) :: xo, x, expected
        logical :: check

        xo = [4.0d0, 3.0d0]
        expected = [2.5d0, -2.0d0]
        x = solve_bfgs(ad_shifted_bowl, xo, args = expected(1))
        check = all_close(x, expected, tol)
        if (.not.check) print *, &
            "Test Failed: solve_bfgs. Expected:", expected, "Found:", x
    end function
end module
