! nonlin_constants.f90

!> @brief \b nonlin_constants
!!
!! @par Purpose
!! To provide various constants used by the NONLIN library.
module nonlin_constants
    use, intrinsic :: iso_fortran_env, only : int32
    use linalg_constants
    implicit none

! ******************************************************************************
! ERROR FLAGS
! ------------------------------------------------------------------------------
    !> A flag denoting no error.
    integer(int32), parameter :: NL_NO_ERROR = 0
    !> An error flag denoting an invalid input.
    integer(int32), parameter :: NL_INVALID_INPUT_ERROR = 201
    !> An error flag denoting an improperly sized array.
    integer(int32), parameter :: NL_ARRAY_SIZE_ERROR = 202
    !> An error denoting that there is insufficient memory available.
    integer(int32), parameter :: NL_OUT_OF_MEMORY_ERROR = LA_OUT_OF_MEMORY_ERROR
    !> An error resulting from an invalid operation.
    integer(int32), parameter :: NL_INVALID_OPERATION_ERROR = &
        LA_INVALID_OPERATION_ERROR
    !> An error resulting from a lack of convergence.
    integer(int32), parameter :: NL_CONVERGENCE_ERROR = LA_CONVERGENCE_ERROR
    !> An error resulting from a divergent condition.
    integer(int32), parameter :: NL_DIVERGENT_BEHAVIOR_ERROR = 206
    !> An error indicating a possible spurious convergence condition.
    integer(int32), parameter :: NL_SPURIOUS_CONVERGENCE_ERROR = 207
    !> An error indicating the user-requested tolerance is too small to be
    !! practical for the problem at hand.
    integer(int32), parameter :: NL_TOLERANCE_TOO_SMALL_ERROR = 208
end module
