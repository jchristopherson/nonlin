! nonlin_brent_example.f90

program example
    use iso_fortran_env
    use nonlin_core
    use nonlin_solve
    implicit none

    ! Local variables
    type(fcn1var_helper) :: obj
    procedure(fcn1var), pointer :: fcn
    type(brent_solver) :: solver
    real(real64) :: x, f
    type(value_pair) :: limits

    ! Define the search limits
    limits%x1 = 1.5d0
    limits%x2 = 5.0d0

    ! Establish the function
    fcn => fcn1
    call obj%set_fcn(fcn)

    ! Solve the equation
    call solver%solve(obj, x, limits, f)

    ! Print the output and the residual
    print '(AF7.5)', "The solution: ", x
    print '(AE10.3)', "The residual: ", f

contains
    ! The function:
    ! f(x) = sin(x) / x, solution: x = x * pi for n = 1, 2, 3, ...
    function fcn1(x) result(f)
        real(real64), intent(in) :: x
        real(real64) :: f
        f = sin(x) / x
    end function
end program
