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
! nonlin_poly_roots_example.f90

program example
    use iso_fortran_env
    use nonlin
    implicit none

    ! Local Variables
    type(polynomial) :: f
    real(real64) :: coeffs(4)
    complex(real64), allocatable :: rts(:)
    integer(int32) :: i

    ! Define the polynomial (x**3 - 2 * x - 1)
    coeffs = [-1.0d0, -2.0d0, 0.0d0, 1.0d0]
    f = coeffs

    ! Compute the polynomial roots
    rts = f%roots()

    ! Display the results
    do i = 1, size(rts)
        print '(A,I0,A,F9.6,A,F9.6,A)', "Root ", i, " = (", real(rts(i), real64), &
            ", ", aimag(rts(i)), ")"
    end do
end program
