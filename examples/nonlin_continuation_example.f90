program example
    use iso_fortran_env
    use nonlin
    implicit none

    integer :: i
    type(vecfcn_helper) :: obj
    procedure(vector_function), pointer :: fcn
    type(continuation_solver) :: solver
    real(real64) :: x(4), f(4)
    type(iteration_behavior) :: tracking
    type(convergence_history) :: history

    fcn => fcn_
    call obj%set_auto_fcn(fcn, 4, 4)

    x = 1.0d0   ! initial guess
    call solver%solve(obj, x, f, ib = tracking)

    ! Display the output
    print *, ""
    print "(A, F0.3, A, F0.3, A, F0.3, A, F0.3, A)", "Solution: (", x(1), &
        ", ", x(2), ", ", x(3), ", ", x(4), ")"
    print "(A, E0.3, A, E0.3, A, E0.3, A, E0.3, A)", "Residual: (", f(1), &
        ", ", f(2), ", ", f(3), ", ", f(4), ")"
    print "(A, I0)", "Iteration Count: ", tracking%iter_count

    history = tracking%get_convergence_history()
    print *, ""
    print "(A, A, A)", "Iteration", achar(9), "Residual"
    print "(A, I0, A, F0.3)", (achar(9), i, achar(9), history%log_residual(i), &
        i = 1, size(history%log_residual))

contains
    subroutine fcn_(x_, f_, args_)
        type(dual), intent(in), dimension(:) :: x_
        type(dual), intent(out), dimension(:) :: f_
        class(*), intent(in), optional :: args_

        f_(1) = x_(1) + 10.0d0 * x_(2)
        f_(2) = sqrt(5.0d0) * (x_(3) - x_(4))
        f_(3) = (x_(2) - 2.0d0 * x_(3))**2
        f_(4) = sqrt(10.0d0) * (x_(1) - x_(4))**2
    end subroutine
end program