! nonlin_solve.f90

!> @brief \b nonlin_solve
!!
!! @par Purpose
!! To provide various routines capapble of solving systems of nonlinear 
!! equations.
module nonlin_solve
    use linalg_constants, only : dp, i32
    use nonlin_types
    use nonlin_linesearch, only : line_search
    use ferror, only : errors
    implicit none
    private

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> A base class for various nonlinear equation solvers.
    type equation_solver
        private
        !> The maximum number of function evaluations allowed per solve.
        integer(i32) :: m_maxEval = 100
        !> The convergence criteria on function values.
        real(dp) :: m_fcnTol = 1.0d-8
        !> The convergence criteria on change in variable values.
        real(dp) :: m_xtol = 1.0d-8
        !> The convergence criteria for the slope of the gradient vector.
        real(dp) :: m_gtol = 1.0d-12
    end type

contains
! ******************************************************************************
! QUASI-NEWTON METHOD (BROYDEN'S METHOD)
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
