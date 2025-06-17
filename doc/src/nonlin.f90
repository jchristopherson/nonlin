module nonlin
    use nonlin_multi_eqn_mult_var
    use nonlin_single_var
    use nonlin_multi_var
    use nonlin_types
    use nonlin_least_squares
    use nonlin_linesearch
    use nonlin_optimize
    use nonlin_solve
    use nonlin_polynomials
    implicit none

    ! NONLIN_MULTI_EQN_MULTI_VAR
    public :: vecfcn
    public :: jacobianfcn
    public :: vecfcn_helper
    public :: equation_solver
    public :: nonlin_solver

    ! NONLIN_SINGLE_VAR
    public :: fcn1var
    public :: fcn1var_helper
    public :: equation_solver_1var
    public :: nonlin_solver_1var

    ! NONLIN_MULTI_VAR
    public :: fcnnvar
    public :: gradientfcn
    public :: fcnnvar_helper
    public :: equation_optimizer
    public :: nonlin_optimize_fcn

    ! NONLIN_TYPES
    public :: iteration_behavior
    public :: value_pair

    ! NONLIN_LEAST_SQUARES
    public :: least_squares_solver

    ! NONLIN_LINESEARCH
    public :: line_search
    public :: limit_search_vector

    ! NONLIN_OPTIMIZE
    public :: nelder_mead
    public :: line_search_optimizer
    public :: bfgs

    ! NONLIN_SOLVE
    public :: line_search_solver
    public :: quasi_newton_solver
    public :: newton_solver
    public :: brent_solver
    public :: newton_1var_solver

    ! NONLIN_POLYNOMIALS
    public :: polynomial
    public :: assignment(=)
    public :: operator(+)
    public :: operator(-)
    public :: operator(*)
end module