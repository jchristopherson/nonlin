! nonlin_equation_optimizer.f90

submodule (nonlin_core) nonlin_equation_optimizer
contains
! ------------------------------------------------------------------------------
    pure module function oe_get_max_eval(this) result(n)
        class(equation_optimizer), intent(in) :: this
        integer(int32) :: n
        n = this%m_maxEval
    end function

! --------------------
    module subroutine oe_set_max_eval(this, n)
        class(equation_optimizer), intent(inout) :: this
        integer(int32), intent(in) :: n
        this%m_maxEval = n
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oe_get_tol(this) result(x)
        class(equation_optimizer), intent(in) :: this
        real(real64) :: x
        x = this%m_tol
    end function

! --------------------
    module subroutine oe_set_tol(this, x)
        class(equation_optimizer), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_tol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oe_get_print_status(this) result(x)
        class(equation_optimizer), intent(in) :: this
        logical :: x
        x = this%m_printStatus
    end function

! --------------------
    module subroutine oe_set_print_status(this, x)
        class(equation_optimizer), intent(inout) :: this
        logical, intent(in) :: x
        this%m_printStatus = x
    end subroutine
end submodule
