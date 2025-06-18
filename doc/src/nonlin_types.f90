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