! nonlin_nm_opt_example.f90

program example
    use iso_fortran_env
    use nonlin_optimize
    use nonlin_core
    use example_problems
    implicit none

    ! Local Variables
    type(nelder_mead) :: solver
    type(fcnnvar_helper) :: obj
    procedure(fcnnvar), pointer :: fcn
    real(real64) :: x(2), fout
    type(iteration_behavior) :: ib

    ! Initialization
    fcn => rosenbrock
    call obj%set_fcn(fcn, 2)

    ! Define an initial guess - the solution is (1, 1)
    call random_number(x)

    ! Call the solver
    call solver%solve(obj, x, fout, ib)

     ! Display the output
     print 100, "Minimum: (", x(1), ", ", x(2), ")"
     print 101, "Function Value: ", fout
     print 102, "Iterations: ", ib%iter_count
     print 102, "Function Evaluations: ", ib%fcn_count

    ! Formatting
100 format(A, F7.5, A, F7.5, A)
101 format(A, E9.3)
102 format(A, I0)
end program
