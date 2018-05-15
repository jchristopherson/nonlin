! nonlin_newton_example_2.f90

program example
    use iso_fortran_env
    use nonlin_core
    use nonlin_solve

    ! Variables
    type(vecfcn_helper) :: obj
    type(newton_solver) :: solver
    procedure(vecfcn), pointer :: fcn
    real(real64) :: x(2), f(2)

    ! Initialization
    fcn => prblm
    call obj%set_fcn(fcn, 2, 2)
    x = [0.0d0, 1.0d0]
    call solver%set_print_status(.true.)

    ! Solve the equations
    call solver%solve(obj, x, f)

    ! Display the output
    print *, ""
    print '(AE12.6AE12.6A)', "Solution: (", x(1), ", ", x(2), ")"
    print '(AE12.6AE12.6A)', "Residual: (", f(1), ", ", f(2), ")"
contains
    ! REF: http://people.sc.fsu.edu/~jburkardt/f_src/test_nonlin/test_nonlin.f90
    ! Problem #3, Powell's badly scaled function
    !
    ! Solution:
    ! x(1) = 1.098159e-5
    ! x(2) = 9.106146
    subroutine prblm(x, f)
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        f(1) = 1.0d4 * x(1) * x(2) - 1.0d0
        f(2) = exp(-x(1)) + exp(-x(2)) - 1.0001d0
    end subroutine
end program
