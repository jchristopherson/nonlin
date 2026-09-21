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
! nonlin_newton_example.f90

program example
    use iso_fortran_env
    use nonlin
    use example_problems
    implicit none

    ! Local Variables
    type(vecfcn_helper) :: obj
    procedure(vecfcn), pointer :: fcn
    type(newton_solver) :: solver
    real(real64) :: x(2), f(2)
    type(iteration_behavior) :: tracking

    ! Assign a pointer to the subroutine containing the equations to solve
    fcn => misc_2fcn
    call obj%set_fcn(fcn, 2, 2) ! There are 2 equations with 2 unknowns

    ! Set up solver parameters
    call solver%set_max_fcn_evals(1000) ! Specify the maximum number of function evaluations before iteration termination
    call solver%set_fcn_tolerance(1.0d-10)
    call solver%set_var_tolerance(1.0d-10)
    call solver%set_gradient_tolerance(1.0d-10)

    ! Tell the solver to print out a status update at each iteration
    call solver%set_print_status(.true.)

    ! Define an initial guess
    x = 1.0d0 ! Equivalent to x = [1.0d0, 1.0d0]

    ! Solve the system of equations, but include solver statistics tracking
    call solver%solve(obj, x, f, tracking)

    ! Display the output
    print *, ""
    print 100, "Solution: (", x(1), ", ", x(2), ")"
    print 101, "Residual: (", f(1), ", ", f(2), ")"
    print 102, "Iteration Count: ", tracking%iter_count
    print 102, "Function Evaluations: ", tracking%fcn_count
    print 102, "Jacobian Evaluations: ", tracking%jacobian_count
    print 103, "Converge on Function Value: ", tracking%converge_on_fcn
    print 103, "Converge on Change in Variable: ", tracking%converge_on_chng
    print 103, "Converge on Zero Slope Gradient Vector: ", tracking%converge_on_zero_diff

    ! Formatting
100 format(A, F7.5, A, F7.5, A)
101 format(A, E9.3, A, E9.3, A)
102 format(A, I0)
103 format(A, L1)
end program
