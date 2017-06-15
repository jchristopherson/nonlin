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
    public :: equation_solver
    public :: nonlin_solver

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
    contains
        !> @brief Gets the maximum number of function evaluations allowed during
        !! a single solve.
        procedure, public :: get_max_fcn_evals => es_get_max_eval
        !> @brief Sets the maximum number of function evaluations allowed during
        !! a single solve.
        procedure, public :: set_max_fcn_evals => es_set_max_eval
        !> @brief Gets the convergence on function value tolerance.
        procedure, public :: get_fcn_tolerance => es_get_fcn_tol
        !> @brief Sets the convergence on function value tolerance.
        procedure, public :: set_fcn_tolerance => es_set_fcn_tol
        !> @brief Gets the convergence on change in variable tolerance.
        procedure, public :: get_var_tolerance => es_get_var_tol
        !> @brief Sets the convergence on change in variable tolerance.
        procedure, public :: set_var_tolerance => es_set_var_tol
        !> @brief Gets the convergence on slope of the gradient vector 
        !! tolerance.
        procedure, public :: get_gradient_tolerance => es_get_grad_tol
        !> @brief Sets the convergence on slope of the gradient vector 
        !! tolerance.
        procedure, public :: set_gradient_tolerance => es_set_grad_tol
        !> @brief Solves the system of equations.
        procedure(nonlin_solver), deferred, public, pass :: solve
    end type

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        !> @brief Describes the interface of a nonlinear equation solver.
        !!
        !! @param[in,out] this The equation_solver-based object.
        !! @param[in] fcn The vecfcn_helper object containing the equations to
        !!  solve.
        !! @param[in,out] x On input, an N-element array containing an initial
        !!  estimate to the solution.  On output, the updated solution estimate.
        !!  N is the number of variables.
        !! @param[out] fvec An M-element array that, on output, will contain
        !!  the values of each equation as evaluated at the variable values
        !!  given in @p x.
        !! @param[out] ib An optional output, that if provided, allows the 
        !!  caller to obtain iteration performance statistics.
        !! @param[out] err An optional errors-based object that if provided can
        !!  be used to retrieve information relating to any errors encountered 
        !!  during execution.  If not provided, a default implementation of the
        !!  errors class is used internally to provide error handling.  The
        !!  possible error codes returned will likely vary from solver to 
        !!  solver.
        subroutine nonlin_solver(this, fcn, x, fvec, ib, err)
            use linalg_constants, only : dp, i32
            import equation_solver
            class(equation_solver), intent(inout) :: this
            class(vecfcn_helper), intent(in) :: fcn
            real(dp), intent(inout), dimension(:) :: x
            real(dp), intent(out), dimension(:) :: fvec
            type(iteration_behavior), optional :: ib
            class(errors), intent(in), optional, target :: err
        end subroutine
    end interface

contains
! ******************************************************************************
! EQUATION_SOLVER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Gets the maximum number of function evaluations allowed during a
    !! single solve.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The maximum number of function evaluations.
    pure function es_get_max_eval(this) result(n)
        class(equation_solver), intent(in) :: this
        integer(i32) :: n
        n = this%m_maxEval
    end function

! --------------------
    !> @brief Sets the maximum number of function evaluations allowed during a
    !! single solve.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] n The maximum number of function evaluations.
    subroutine es_set_max_eval(this, n)
        class(equation_solver), intent(inout) :: this
        integer(i32), intent(in) :: n
        this%m_maxEval = n
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on function value tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    pure function es_get_fcn_tol(this) result(x)
        class(equation_solver), intent(in) :: this
        real(dp) :: x
        x = this%m_fcnTol
    end function

! --------------------
    !> @brief Sets the convergence on function value tolerance.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] x The tolerance value.
    subroutine es_set_fcn_tol(this, x)
        class(equation_solver), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_fcnTol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on change in variable tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    pure function es_get_var_tol(this) result(x)
        class(equation_solver), intent(in) :: this
        real(dp) :: x
        x = this%m_xtol
    end function

! --------------------
    !> @brief Sets the convergence on change in variable tolerance.
    !!
    !! @param[in,out] this The equation_solver object.
    !! @param[in] x The tolerance value.
    subroutine es_set_var_tol(this, x)
        class(equation_solver), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_xtol = x
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Gets the convergence on slope of the gradient vector tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    pure function es_get_grad_tol(this) result(x)
        class(equation_solver), intent(in) :: this
        real(dp) :: x
        x = this%m_gtol
    end function

! --------------------
    !> @brief Sets the convergence on slope of the gradient vector tolerance.
    !!
    !! @param[in] this The equation_solver object.
    !! @return The tolerance value.
    subroutine es_set_grad_tol(this, x)
        class(equation_solver), intent(inout) :: this
        real(dp), intent(in) :: x
        this%m_gtol = x
    end subroutine

! ******************************************************************************
! QUASI-NEWTON METHOD (BROYDEN'S METHOD)
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
