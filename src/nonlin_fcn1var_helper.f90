! nonlin_fcn1var_helper.f90

submodule (nonlin_core) nonlin_fcn1var_helper
contains
! ------------------------------------------------------------------------------
    module function f1h_fcn(this, x) result(f)
        class(fcn1var_helper), intent(in) :: this
        real(real64), intent(in) :: x
        real(real64) :: f
        if (associated(this%m_fcn)) then
            f = this%m_fcn(x)
        end if
    end function

! ------------------------------------------------------------------------------
    pure module function f1h_is_fcn_defined(this) result(x)
        class(fcn1var_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    module subroutine f1h_set_fcn(this, fcn)
        class(fcn1var_helper), intent(inout) :: this
        procedure(fcn1var), intent(in), pointer :: fcn
        this%m_fcn => fcn
    end subroutine
end submodule
