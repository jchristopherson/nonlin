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
! nonlin_least_squares_example.f90

program example
    use iso_fortran_env
    use nonlin
    use example_problems
    implicit none

    ! Local Variables
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    type(least_squares_solver) :: solver
    real(real64) :: x(4), f(21) ! There are 4 coefficients and 21 data points

    ! Locate the routine containing the equations to solve
    fcn => lsq_poly_fit_fcn
    call obj%set_fcn(fcn, 21, 4)

    ! Define an initial guess
    x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0, 1.0d0, 1.0d0]

    ! Solve
    call solver%solve(obj, x, f)

    ! Display the output
    print 100, "c0: ", x(4)
    print 100, "c1: ", x(3)
    print 100, "c2: ", x(2)
    print 100, "c3: ", x(1)
    print 101, "Max Residual: ", maxval(abs(f))

    ! Formatting
100 format(A, F12.10)
101 format(A, F7.5)
end program
