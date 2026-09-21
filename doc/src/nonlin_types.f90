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
module nonlin_types
    use iso_fortran_env
    implicit none
    private
    public :: iteration_behavior
    public :: value_pair

    type iteration_behavior
        !! Defines a set of parameters that describe the behavior of the
        !! iteration process.
        integer(int32) :: iter_count
            !! Specifies the number of iterations performed.
        integer(int32) :: fcn_count
            !! Specifies the number of function evaluations performed.
        integer(int32) :: jacobian_count
            !! Specifies the number of Jacobian evaluations performed.
        integer(int32) :: gradient_count
            !! Specifies the number of gradient vector evaluations performed.
        logical :: converge_on_fcn
            !! True if the solution converged as a result of a zero-valued
            !! function; else, false.
        logical :: converge_on_chng
            !! True if the solution converged as a result of no appreciable
            !! change in solution points between iterations; else, false.
        logical :: converge_on_zero_diff
            !! True if the solution appears to have settled on a stationary
            !! point such that the gradient of the function is zero-valued; 
            !! else, false.
    end type

    type value_pair
        !! Defines a pair of numeric values.
        real(real64) :: x1
            !! Value 1.
        real(real64) :: x2
            !! Value 2.
    end type
end module