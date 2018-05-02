! nonlin_fcnnvar_helper.f90

submodule (nonlin_core) nonlin_fcnnvar_helper
contains
! ------------------------------------------------------------------------------
    module function fnh_fcn(this, x) result(f)
        class(fcnnvar_helper), intent(in) :: this
        real(real64), intent(in), dimension(:) :: x
        real(real64) :: f
        if (associated(this%m_fcn)) then
            f = this%m_fcn(x)
        end if
    end function

! ------------------------------------------------------------------------------
    pure module function fnh_is_fcn_defined(this) result(x)
        class(fcnnvar_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    module subroutine fnh_set_fcn(this, fcn, nvar)
        class(fcnnvar_helper), intent(inout) :: this
        procedure(fcnnvar), intent(in), pointer :: fcn
        integer(int32), intent(in) :: nvar
        this%m_fcn => fcn
        this%m_nvar = nvar
    end subroutine

! ------------------------------------------------------------------------------
    pure module function fnh_get_nvar(this) result(n)
        class(fcnnvar_helper), intent(in) :: this
        integer(int32) :: n
        n = this%m_nvar
    end function

! ------------------------------------------------------------------------------
    module subroutine fnh_set_grad(this, fcn)
        class(fcnnvar_helper), intent(inout) :: this
        procedure(gradientfcn), pointer, intent(in) :: fcn
        this%m_grad => fcn
    end subroutine

! ------------------------------------------------------------------------------
    pure module function fnh_is_grad_defined(this) result(x)
        class(fcnnvar_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_grad)
    end function

! ------------------------------------------------------------------------------
    module subroutine fnh_grad_fcn(this, x, g, fv, err)
        ! Arguments
        class(fcnnvar_helper), intent(in) :: this
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: g
        real(real64), intent(in), optional :: fv
        integer(int32), intent(out), optional :: err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: j, n, flag
        real(real64) :: eps, epsmch, h, temp, f, f1

        ! Initialization
        if (present(err)) err = 0
        ! n = this%m_nvar
        n = this%get_variable_count()

        ! Input Checking
        flag = 0
        if (size(x) /= n) then
            flag = 2
        else if (size(g) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: Incorrectly sized input arrays
            if (present(err)) err = flag
            return
        end if

        ! Process
        if (.not.this%is_fcn_defined()) return
        if (this%is_gradient_defined()) then
            ! Call the user-defined gradient routine
            call this%m_grad(x, g)
        else
            ! Compute the gradient via finite differences
            if (present(fv)) then
                f = fv
            else
                f = this%fcn(x)
            end if

            ! Establish step size factors
            epsmch = epsilon(epsmch)
            eps = sqrt(epsmch)

            ! Compute the derivatives
            do j = 1, n
                temp = x(j)
                h = eps * abs(temp)
                if (h == zero) h = eps
                x(j) = temp + h
                f1 = this%fcn(x)
                x(j) = temp
                g(j) = (f1 - f) / h
            end do
        end if
    end subroutine


end submodule
