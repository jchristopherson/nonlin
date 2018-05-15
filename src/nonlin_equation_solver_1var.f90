! nonlin_equation_solver_1var.f90

submodule (nonlin_core) nonlin_equation_solver_1var
contains
! ------------------------------------------------------------------------------
    pure module function es1_get_max_eval(this) result(n)
        class(equation_solver_1var), intent(in) :: this
        integer(int32) :: n
        n = this%m_maxEval
    end function

! --------------------
    module subroutine es1_set_max_eval(this, n)
        class(equation_solver_1var), intent(inout) :: this
        integer(int32), intent(in) :: n
        this%m_maxEval = n
    end subroutine

! ------------------------------------------------------------------------------
    pure module function es1_get_fcn_tol(this) result(x)
        class(equation_solver_1var), intent(in) :: this
        real(real64) :: x
        x = this%m_fcnTol
    end function

! --------------------
    module subroutine es1_set_fcn_tol(this, x)
        class(equation_solver_1var), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_fcnTol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function es1_get_var_tol(this) result(x)
        class(equation_solver_1var), intent(in) :: this
        real(real64) :: x
        x = this%m_xtol
    end function

! --------------------
    module subroutine es1_set_var_tol(this, x)
        class(equation_solver_1var), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_xtol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function es1_get_print_status(this) result(x)
        class(equation_solver_1var), intent(in) :: this
        logical :: x
        x = this%m_printStatus
    end function

! --------------------
    module subroutine es1_set_print_status(this, x)
        class(equation_solver_1var), intent(inout) :: this
        logical, intent(in) :: x
        this%m_printStatus = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function es1_get_diff_tol(this) result(x)
        class(equation_solver_1var), intent(in) :: this
        real(real64) :: x
        x = this%m_difftol
    end function

! --------------------
    module subroutine es1_set_diff_tol(this, x)
        class(equation_solver_1var), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_difftol = x
    end subroutine
end submodule
