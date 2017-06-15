! nonlin_types.f90

module nonlin_types
    use linalg_constants, only : dp, i32
    implicit none
    private
    public :: vecfcn
    public :: jacobianfcn
    public :: vecfcn_helper
    public :: NL_INVALID_INPUT_ERROR
    public :: NL_ARRAY_SIZE_ERROR
    public :: NL_OUT_OF_MEMORY_ERROR
    public :: NL_INVALID_OPERATION_ERROR
    public :: NL_CONVERGENCE_ERROR
    public :: NL_DIVERGENT_BEHAVIOR_ERROR

! ******************************************************************************
! ERROR FLAGS
! ------------------------------------------------------------------------------
    !> An error flag denoting an invalid input.
    integer, parameter :: NL_INVALID_INPUT_ERROR = 201
    !> An error flag denoting an improperly sized array.
    integer, parameter :: NL_ARRAY_SIZE_ERROR = 202
    !> An error denoting that there is insufficient memory available.
    integer, parameter :: NL_OUT_OF_MEMORY_ERROR = 203
    !> An error resulting from an invalid operation.
    integer, parameter :: NL_INVALID_OPERATION_ERROR = 204
    !> An error resulting from a lack of convergence.
    integer, parameter :: NL_CONVERGENCE_ERROR = 205
    !> An error resulting from a divergent condition.
    integer, parameter :: NL_DIVERGENT_BEHAVIOR_ERROR = 206

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        !> @brief Describes an M-element vector-valued function of N-variables.
        !!
        !! @param[in] x An N-element array containing the independent variables.
        !! @param[out] f An M-element array that, on output, contains the values
        !!  of the M functions.
        subroutine vecfcn(x, f)
            use linalg_constants, only : dp
            real(dp), intent(in), dimension(:) :: x
            real(dp), intent(out), dimension(:) :: f
        end subroutine

        !> @brief Describes a routine capable of computing the Jacobian matrix
        !! of M functions of N unknowns.
        !!
        !! @param[in] x An N-element array containing the independent variables.
        !! @param[out] jac An M-by-N matrix where the Jacobian will be written.
        subroutine jacobianfcn(x, jac)
            use linalg_constants, only : dp
            real(dp), intent(in), dimension(:) :: x
            real(dp), intent(out), dimension(:,:) :: jac
        end subroutine
    end interface

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    !> @brief Defines a type capable of encapsulating a system of nonlinear 
    !! equations of the form: F(X) = 0.
    type vecfcn_helper
        !> A pointer to the encapsulated vecfcn routine.
        procedure(vecfcn), pointer, nopass :: m_fcn => null()
        !> A pointer to the jacobian routine - null if no routine is supplied.
        procedure(jacobianfcn), pointer, nopass :: m_jac => null()
        !> The number of functions in m_fcn
        integer(i32) :: m_nfcn = 0
        !> The number of variables in m_fcn
        integer(i32) :: m_nvar = 0
    contains
        !> @brief Establishes a pointer to the routine containing the system of
        !!  equations to solve.
        procedure, public :: set_fcn => vfh_set_fcn
        !> @brief Establishes a pointer to the routine for computing the 
        !! Jacobian matrix of the system of equations.  If no routine is 
        !! defined, the Jacobian matrix will be computed numerically (this is 
        !! the default state).
        procedure, public :: set_jacobian => vfh_set_jac
        !> @brief Tests if the pointer to the subroutine containing the system
        !! of equations to solve has been assigned.
        procedure, public :: is_fcn_defined => vfh_is_fcn_defined
        !> @brief Tests if the pointer to the subroutine containing the system
        !! of equations to solve has been assigned.
        procedure, public :: is_jacobian_defined => vfh_is_jac_defined
        !> @brief Executes the routine containing the system of equations to
        !! solve.  No action is taken if the pointer to the subroutine has not 
        !! been defined.
        procedure, public :: fcn => vfh_fcn
        !> @brief Executes the routine containing the Jacobian matrix if 
        !! supplied.  If not supplied, the Jacobian is computed via finite 
        !! differences.
        procedure, public :: jacobian => vfh_jac_fcn
        !> @brief Gets the number of equations in this system.
        procedure, public :: get_equation_count => vfh_get_nfcn
        !> @brief Gets the number of variables in this system.
        procedure, public :: get_variable_count => vfh_get_nvar
    end type


contains
! ******************************************************************************
! VECFCN_HELPER MEMBERS
! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine containing the system of
    !!  equations to solve.
    !! 
    !! @param[in,out] this The vecfcn_helper object.
    !! @param[in] fcn The function pointer.
    !! @param[in] nfcn The number of functions.
    !! @param[in] nvar The number of variables.
    subroutine vfh_set_fcn(this, fcn, nfcn, nvar)
        class(vecfcn_helper), intent(inout) :: this
        procedure(vecfcn), intent(in), pointer :: fcn
        integer(i32), intent(in) :: nfcn, nvar
        this%m_fcn => fcn
        this%m_nfcn = nfcn
        this%m_nvar = nvar
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Establishes a pointer to the routine for computing the Jacobian
    !! matrix of the system of equations.  If no routine is defined, the
    !! Jacobian matrix will be computed numerically (this is the default state).
    !!
    !! @param[in,out] this The vecfcn_helper object.
    !! @param[in] jac The function pointer.
    subroutine vfh_set_jac(this, jac)
        class(vecfcn_helper), intent(inout) :: this
        procedure(jacobianfcn), intent(in), pointer :: jac
        this%m_jac => jac
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the subroutine containing the system of
    !! equations to solve has been assigned.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function vfh_is_fcn_defined(this) result(x)
        class(vecfcn_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    !> @brief Tests if the pointer to the subroutine containing the system of 
    !! equations to solve has been assigned.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @return Returns true if the pointer has been assigned; else, false.
    pure function vfh_is_jac_defined(this) result(x)
        class(vecfcn_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_jac)
    end function

! ------------------------------------------------------------------------------
    !> @brief Executes the routine containing the system of equations to solve.
    !! No action is taken if the pointer to the subroutine has not been defined.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @param[in] x An N-element array containing the independent variables.
    !! @param[out] f An M-element array that, on output, contains the values
    !!  of the M functions.
    subroutine vfh_fcn(this, x, f)
        class(vecfcn_helper), intent(in) :: this
        real(dp), intent(in), dimension(:) :: x
        real(dp), intent(out), dimension(:) :: f
        if (this%is_fcn_defined()) then
            call this%m_fcn(x, f)
        end if
    end subroutine

! ------------------------------------------------------------------------------
    !> @brief Executes the routine containing the Jacobian matrix if supplied.
    !! If not supplied, the Jacobian is computed via finite differences.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @param[in] x An N-element array containing the independent variabls 
    !!  defining the point about which the derivatives will be calculated.
    !! @param[out] jac An M-by-N matrix where, on output, the Jacobian will
    !!  be written.
    !! @param[in] fv An optional M-element array containing the function values
    !!  at @p x.  If not supplied, the function values are computed at @p x.
    !! @param[out] work An optional input, that if provided, prevents any local
    !!  memory allocation.  If not provided, the memory required is allocated
    !!  within.  If provided, the length of the array must be at least
    !!  @p olwork.  Notice, a workspace array is only utilized if the user does
    !!  not provide a routine for computing the Jacobian.
    !! @param[out] olwork An optional output used to determine workspace size.
    !!  If supplied, the routine determines the optimal size for @p work, and
    !!  returns without performing any actual calculations.
    !! @param[out] err An optional integer output that can be used to determine
    !!  error status.  If not used, and an error is encountered, the routine
    !!  simply returns silently.  If used, the following error codes identify
    !!  error status:
    !!  - 0: No error has occurred.
    !!  - n: A positive integer denoting the index of an invalid input.
    !!  - -1: Indicates internal memory allocation failed.
    subroutine vfh_jac_fcn(this, x, jac, fv, work, olwork, err)
        ! Arguments
        class(vecfcn_helper), intent(in) :: this
        real(dp), intent(inout), dimension(:) :: x
        real(dp), intent(out), dimension(:,:) :: jac
        real(dp), intent(in), dimension(:), optional, target :: fv
        real(dp), intent(out), dimension(:), optional, target :: work
        integer(i32), intent(out), optional :: olwork, err

        ! Parameters
        real(dp), parameter :: zero = 0.0d0

        ! Local Variables
        integer(i32) :: j, m, n, lwork, flag
        real(dp) :: eps, epsmch, h, temp
        real(dp), pointer, dimension(:) :: fptr, f1ptr
        real(dp), allocatable, target, dimension(:) :: wrk

        ! Initialization
        if (present(err)) err = 0
        m = this%m_nfcn
        n = this%m_nvar

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
    !> @brief Gets the number of equations in this system.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @return The function count.
    pure function vfh_get_nfcn(this) result(n)
        class(vecfcn_helper), intent(in) :: this
        integer(i32) :: n
        n = this%m_nfcn
    end function

! ------------------------------------------------------------------------------
    !> @brief Gets the number of variables in this system.
    !!
    !! @param[in] this The vecfcn_helper object.
    !! @return The number of variables.
    pure function vfh_get_nvar(this) result(n)
        class(vecfcn_helper), intent(in) :: this
        integer(i32) :: n
        n = this%m_nvar
    end function

! ------------------------------------------------------------------------------
end module
