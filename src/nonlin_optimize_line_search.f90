! nonlin_optimize_line_search.f90

submodule (nonlin_optimize) nonlin_optimize_line_search
contains
! ------------------------------------------------------------------------------
    module subroutine lso_get_line_search(this, ls)
        class(line_search_optimizer), intent(in) :: this
        class(line_search), intent(out), allocatable :: ls
        if (allocated(this%m_lineSearch)) &
            allocate(ls, source = this%m_lineSearch)
    end subroutine

! ----------------------
    module subroutine lso_set_line_search(this, ls)
        class(line_search_optimizer), intent(inout) :: this
        class(line_search), intent(in) :: ls
        if (allocated(this%m_lineSearch)) deallocate(this%m_lineSearch)
        allocate(this%m_lineSearch, source = ls)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine lso_set_default(this)
        class(line_search_optimizer), intent(inout) :: this
        type(line_search) :: ls
        call this%set_line_search(ls)
    end subroutine

! ------------------------------------------------------------------------------
    pure module function lso_is_line_search_defined(this) result(x)
        class(line_search_optimizer), intent(in) :: this
        logical :: x
        x = allocated(this%m_lineSearch)
    end function

! ------------------------------------------------------------------------------
    pure module function lso_get_use_search(this) result(x)
        class(line_search_optimizer), intent(in) :: this
        logical :: x
        x = this%m_useLineSearch
    end function

! --------------------
    module subroutine lso_set_use_search(this, x)
        class(line_search_optimizer), intent(inout) :: this
        logical, intent(in) :: x
        this%m_useLineSearch = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function lso_get_var_tol(this) result(x)
        class(line_search_optimizer), intent(in) :: this
        real(real64) :: x
        x = this%m_xtol
    end function

! --------------------
    module subroutine lso_set_var_tol(this, x)
        class(line_search_optimizer), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_xtol = x
    end subroutine
end submodule
