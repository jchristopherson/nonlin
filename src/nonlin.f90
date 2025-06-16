module nonlin
    use nonlin_multi_eqn_mult_var
    use nonlin_single_var
    use nonlin_multi_var
    use nonlin_types
    use nonlin_least_squares
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
end module