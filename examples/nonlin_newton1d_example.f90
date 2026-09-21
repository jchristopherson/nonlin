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
! nonlin_newton1d_example.f90

program example
    use iso_fortran_env
    use nonlin
    use example_problems
    implicit none

    ! Local variables
    type(fcn1var_helper) :: obj
    procedure(fcn1var), pointer :: fcn
    type(newton_1var_solver) :: solver
    real(real64) :: x, f
    type(value_pair) :: limits
    type(iteration_behavior) :: tracking

    ! Define the search limits
    limits%x1 = 2.0d0
    limits%x2 = -2.0d0

    ! Establish the function
    fcn => fcn_1var
    call obj%set_fcn(fcn)

    ! Allow the solver to print out updates at each iteration
    call solver%set_print_status(.true.)

    ! Solve the equation
    call solver%solve(obj, x, limits, f, ib = tracking)

    ! Print the output and the residual
    print *, ""
    print 100, "The solution: ", x
    print 101, "The residual: ", f
    print 102, "Iterations: ", tracking%iter_count
    print 102, "Function Evaluations: ", tracking%fcn_count
    print 102, "Derivative Evaluations: ", tracking%jacobian_count
    print 103, "Converge on Function Value: ", tracking%converge_on_fcn
    print 103, "Converge on Change in Variable: ", tracking%converge_on_chng
    print 103, "Converge on Derivative: ", tracking%converge_on_zero_diff

    ! Formatting
100 format(A, F7.5)
101 format(A, E10.3)
102 format(A, I0)
103 format(A, L1)
end program
