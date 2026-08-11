module nonlin_error_handling
    use iso_fortran_env
    use linalg_errors, only : LA_OUT_OF_MEMORY_ERROR, LA_CONVERGENCE_ERROR, &
        LA_INVALID_OPERATION_ERROR
    implicit none

! ******************************************************************************
! ERROR FLAGS
! ------------------------------------------------------------------------------
    integer(int32), parameter :: NL_NO_ERROR = 0
        !! A flag denoting no error.
    integer(int32), parameter :: NL_INVALID_INPUT_ERROR = 201
        !! An error flag denoting an invalid input.
    integer(int32), parameter :: NL_ARRAY_SIZE_ERROR = 202
        !! An error flag denoting an improperly sized array.
    integer(int32), parameter :: NL_OUT_OF_MEMORY_ERROR = LA_OUT_OF_MEMORY_ERROR
        !! An error denoting that there is insufficient memory available.
    integer(int32), parameter :: NL_INVALID_OPERATION_ERROR = &
        LA_INVALID_OPERATION_ERROR
        !! An error resulting from an invalid operation.
    integer(int32), parameter :: NL_CONVERGENCE_ERROR = LA_CONVERGENCE_ERROR
        !! An error resulting from a lack of convergence.
    integer(int32), parameter :: NL_DIVERGENT_BEHAVIOR_ERROR = 206
        !! An error resulting from a divergent condition.
    integer(int32), parameter :: NL_SPURIOUS_CONVERGENCE_ERROR = 207
        !! An error indicating a possible spurious convergence condition.
    integer(int32), parameter :: NL_TOLERANCE_TOO_SMALL_ERROR = 208
        !! An error indicating the user-requested tolerance is too small to be
        !! practical for the problem at hand.
    integer(int32), parameter :: NL_INDEX_OUT_OF_RANGE_ERROR = 209
        !! An error flag denoting an array index is out of range.
    integer(int32), parameter :: NL_DIVIDE_BY_ZERO_ERROR = 210
        !! An error flag denoting a division by zero error.
    integer(int32), parameter :: NL_UNDEFINED_FUNCTION_ERROR = 211
        !! An error flag denoting that the routine containing the function to
        !! solve has not been properly defined.
    integer(int32), parameter :: NL_UNDERDEFINED_PROBLEM_ERROR = 212
        !! An error flag denoting an underdefined problem error.

end module