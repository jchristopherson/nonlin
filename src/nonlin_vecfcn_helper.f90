! nonlin_vecfcn_helper.f90

submodule (nonlin_core) nonlin_vecfcn_helper
contains
! ------------------------------------------------------------------------------
    module subroutine vfh_set_fcn(this, fcn, nfcn, nvar)
        class(vecfcn_helper), intent(inout) :: this
        procedure(vecfcn), intent(in), pointer :: fcn
        integer(int32), intent(in) :: nfcn, nvar
        this%m_fcn => fcn
        this%m_nfcn = nfcn
        this%m_nvar = nvar
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine vfh_set_jac(this, jac)
        class(vecfcn_helper), intent(inout) :: this
        procedure(jacobianfcn), intent(in), pointer :: jac
        this%m_jac => jac
    end subroutine

! ------------------------------------------------------------------------------
    pure module function vfh_is_fcn_defined(this) result(x)
        class(vecfcn_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    pure module function vfh_is_jac_defined(this) result(x)
        class(vecfcn_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_jac)
    end function

! ------------------------------------------------------------------------------
    module subroutine vfh_fcn(this, x, f)
        class(vecfcn_helper), intent(in) :: this
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        if (this%is_fcn_defined()) then
            call this%m_fcn(x, f)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine vfh_jac_fcn(this, x, jac, fv, work, olwork, err)
        ! Arguments
        class(vecfcn_helper), intent(in) :: this
        real(real64), intent(inout), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: jac
        real(real64), intent(in), dimension(:), optional, target :: fv
        real(real64), intent(out), dimension(:), optional, target :: work
        integer(int32), intent(out), optional :: olwork, err

        ! Parameters
        real(real64), parameter :: zero = 0.0d0

        ! Local Variables
        integer(int32) :: j, m, n, lwork, flag
        real(real64) :: eps, epsmch, h, temp
        real(real64), pointer, dimension(:) :: fptr, f1ptr
        real(real64), allocatable, target, dimension(:) :: wrk

        ! Initialization
        if (present(err)) err = 0
        ! m = this%m_nfcn
        ! n = this%m_nvar
        m = this%get_equation_count()
        n = this%get_variable_count()

        ! Input Checking
        flag = 0
        if (size(x) /= n) then
            flag = 2
        else if (size(jac, 1) /= m .or. size(jac, 2) /= n) then
            flag = 3
        end if
        if (flag /= 0) then
            ! ERROR: Incorrectly sized input arrays
            if (present(err)) err = flag
            return
        end if

        ! Process
        if (.not.this%is_fcn_defined()) return
        if (associated(this%m_jac)) then
            ! Workspace Query
            if (present(olwork)) then
                olwork = 0
                return
            end if

            ! Call the user-defined Jacobian routine
            call this%m_jac(x, jac)
        else
            ! Compute the Jacobian via finite differences
            if (present(fv)) then
                lwork = m
            else
                lwork = 2 * m
            end if

            if (present(olwork)) then
                ! The user is just making a workspace query.  Simply return the
                ! workspace length, and exit the routine.
                olwork = lwork
                return
            end if

            ! Local Memory Allocation
            if (present(work)) then
                if (size(work) < lwork) then
                    ! ERROR: Workspace is too small
                    if (present(err)) err = 5
                    return
                end if
                f1ptr => work(1:m)
                if (present(fv)) then
                    if (size(fv) < m) then
                        ! ERROR: Function vector too small
                        if (present(err)) err = 4
                        return
                    end if
                    fptr => fv(1:m)
                else
                    fptr => work(m+1:2*m)
                    call this%fcn(x, fptr)
                end if
            else
                allocate(wrk(lwork), stat = flag)
                if (flag /= 0) then
                    ! ERROR: Memory issues
                    if (present(err)) err = -1
                    return
                end if
                f1ptr => wrk(1:m)
                if (present(fv)) then
                    fptr => fv(1:m)
                else
                    fptr => wrk(m+1:2*m)
                    call this%fcn(x, fptr)
                end if
            end if

            ! Establish step size factors
            epsmch = epsilon(epsmch)
            eps = sqrt(epsmch)

            ! Compute the derivatives via finite differences
            do j = 1, n
                temp = x(j)
                h = eps * abs(temp)
                if (h == zero) h = eps
                x(j) = temp + h
                call this%fcn(x, f1ptr)
                x(j) = temp
                jac(:,j) = (f1ptr - fptr) / h
            end do
        end if
    end subroutine

! ------------------------------------------------------------------------------
    pure module function vfh_get_nfcn(this) result(n)
        class(vecfcn_helper), intent(in) :: this
        integer(int32) :: n
        n = this%m_nfcn
    end function

! ------------------------------------------------------------------------------
    pure module function vfh_get_nvar(this) result(n)
        class(vecfcn_helper), intent(in) :: this
        integer(int32) :: n
        n = this%m_nvar
    end function
end submodule
