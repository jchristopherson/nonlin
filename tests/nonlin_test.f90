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
! nonlin_test.f90

! The testing application for the NONLIN library.
program main
    ! Imported Modules
    use nonlin_test_jacobian
    use nonlin_test_solve
    use nonlin_test_poly
    use nonlin_test_optimize

    ! Local Variables
    logical :: rst, overall

    ! Initialization
    overall = .true.

    ! Tests
    rst = test_jacobian_1()
    if (.not.rst) overall = .false.

    rst = test_jacobian_2()
    if (.not.rst) overall = .false.

    rst = test_quasinewton_1()
    if (.not.rst) overall = .false.

    rst = test_quasinewton_2()
    if (.not.rst) overall = .false.

    rst = test_newton_1()
    if (.not.rst) overall = .false.

    rst = test_newton_2()
    if (.not.rst) overall = .false.

    rst = test_least_squares_1()
    if (.not.rst) overall = .false.

    rst = test_least_squares_2()
    if (.not.rst) overall = .false.

    call test_least_squares_3()

    rst = test_brent_1()
    if (.not.rst) overall = .false.

    call test_poly_fit()
    rst = test_poly_roots()
    if (.not.rst) overall = .false.

    rst = test_poly_add()
    if (.not.rst) overall = .false.

    rst = test_poly_subtract()
    if (.not.rst) overall = .false.

    rst = test_poly_multiply()
    if (.not.rst) overall = .false.

    rst = test_poly_divide()
    if (.not.rst) overall = .false.

    rst = test_nelder_mead_1()
    if (.not.rst) overall = .false.

    rst = test_nelder_mead_2()
    if (.not.rst) overall = .false.

    rst = test_nelder_mead_3()
    if (.not.rst) overall = .false.

    rst = test_bfgs_1()
    if (.not.rst) overall = .false.

    rst = test_bfgs_2()
    if (.not.rst) overall = .false.

    rst = test_bfgs_3()
    if (.not.rst) overall = .false.

    rst = test_newton_3()
    if (.not.rst) overall = .false.

    rst = test_quasinewton_3()
    if (.not.rst) overall = .false.

    rst = test_quasinewton_4()
    if (.not.rst) overall = .false.

    rst = test_newton_4()
    if (.not.rst) overall = .false.

    rst = test_least_squares_4()
    if (.not.rst) overall = .false.

    rst = test_brent_2()
    if (.not.rst) overall = .false.

    rst = test_newton_1var_1()
    if (.not.rst) overall = .false.

    rst = test_newton_1var_2()
    if (.not.rst) overall = .false.

    rst = test_constrained_least_squares_1()
    if (.not.rst) overall = .false.

    rst = test_constrained_least_squares_2()
    if (.not.rst) overall = .false.

    rst = test_constrained_least_squares_3()
    if (.not.rst) overall = .false.

    rst = test_constrained_least_squares_4()
    if (.not.rst) overall = .false.

    rst = test_constrained_least_squares_bounds()
    if (.not.rst) overall = .false.

    ! End
    if (.not.overall) stop -1
end program
