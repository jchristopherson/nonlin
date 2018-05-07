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

! ------------------------------------------------------------------------------
    pure module function f1h_is_diff_defined(this) result(x)
        class(fcn1var_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_diff)
    end function

! ------------------------------------------------------------------------------
    module function f1h_diff_fcn(this, x, f) result(df)
        ! Arguments
        class(fcn1var_helper), intent(in) :: this
        real(real64), intent(in) :: x
        real(real64), intent(in), optional :: f
        real(real64) :: df

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        real(real64) :: eps, epsmch, h, temp, f1, f0

        ! Initialization
        epsmch = epsilon(epsmch)
        eps = sqrt(epsmch)

        ! Process
        if (this%is_derivative_defined()) then
            ! Use the user-defined routine to compute the derivative
            df = this%m_diff(x)
        else
            ! Compute the derivative via a forward difference
            h = eps * abs(x)
            if (h < epsmch) h = eps
            temp = x + h
            f1 = this%fcn(temp)
            if (present(f)) then
                f0 = f
            else
                f0 = this%fcn(x)
            end if
            df = (f1 - f0) / h
        end if
    end function

! ------------------------------------------------------------------------------
    module subroutine f1h_set_diff(this, diff)
        class(fcn1var_helper), intent(inout) :: this
        procedure(fcn1var), pointer, intent(in) :: diff
        this%m_diff => diff
    end subroutine

! ------------------------------------------------------------------------------
end submodule
