! nonlin_solve_line_search.f90

submodule (nonlin_solve) nonlin_solve_line_search
contains
! ------------------------------------------------------------------------------
    module subroutine lss_get_line_search(this, ls)
        class(line_search_solver), intent(in) :: this
        class(line_search), intent(out), allocatable :: ls
        if (allocated(this%m_lineSearch)) &
            allocate(ls, source = this%m_lineSearch)
    end subroutine

! ----------------------
    module subroutine lss_set_line_search(this, ls)
        class(line_search_solver), intent(inout) :: this
        class(line_search), intent(in) :: ls
        if (allocated(this%m_lineSearch)) deallocate(this%m_lineSearch)
        allocate(this%m_lineSearch, source = ls)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine lss_set_default(this)
        class(line_search_solver), intent(inout) :: this
        type(line_search) :: ls
        call this%set_line_search(ls)
    end subroutine

! ------------------------------------------------------------------------------
    pure module function lss_is_line_search_defined(this) result(x)
        class(line_search_solver), intent(in) :: this
        logical :: x
        x = allocated(this%m_lineSearch)
    end function

! ------------------------------------------------------------------------------
    pure module function lss_get_use_search(this) result(x)
        class(line_search_solver), intent(in) :: this
        logical :: x
        x = this%m_useLineSearch
    end function

! --------------------
    module subroutine lss_set_use_search(this, x)
        class(line_search_solver), intent(inout) :: this
        logical, intent(in) :: x
        this%m_useLineSearch = x
    end subroutine
end submodule
