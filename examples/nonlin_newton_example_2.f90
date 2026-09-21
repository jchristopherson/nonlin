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
! nonlin_newton_example_2.f90

program example
    use iso_fortran_env
    use nonlin
    use example_problems
    implicit none

    ! Variables
    type(vecfcn_helper) :: obj
    type(newton_solver) :: solver
    procedure(vecfcn), pointer :: fcn
    real(real64) :: x(2), f(2)

    ! Initialization
    fcn => powell_bad
    call obj%set_fcn(fcn, 2, 2)
    x = [0.0d0, 1.0d0]
    call solver%set_print_status(.true.)

    ! Solve the equations
    call solver%solve(obj, x, f)

    ! Display the output
    print *, ""
    print 100, "Solution: (", x(1), ", ", x(2), ")"
    print 100, "Residual: (", f(1), ", ", f(2), ")"

    ! Formatting
100 format(A, E12.6, A, E12.6, A)
end program
