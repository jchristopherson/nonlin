! MIT License
!
! Copyright (c) 2017-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
module nonlin_types
    use iso_fortran_env
    use ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use nonlin_error_handling
    implicit none
    private
    public :: convergence_history
    public :: iteration_behavior
    public :: value_pair

    type convergence_history
        !! A collection of arrays providing solver-related data from each
        !! iteration.
        real(real64), allocatable, dimension(:) :: log_residual
            !! Base-10 logarithm of the residual norm at each iteration.
        real(real64), allocatable, dimension(:) :: log_step_change
            !! Base-10 logarithm of the norm of the direction vector for each
            !! iteration.
        real(real64), allocatable, dimension(:) :: log_gradient
            !! Base-10 logarithm of the norm of the gradient for each iteration.
    end type

    type iteration_behavior
        !! Defines a set of parameters that describe the behavior of the
        !! iteration process.
        integer :: iter_count
            !! Specifies the number of iterations performed.
        integer :: fcn_count
            !! Specifies the number of function evaluations performed.
        integer :: jacobian_count
            !! Specifies the number of Jacobian evaluations performed.
        integer :: gradient_count
            !! Specifies the number of gradient vector evaluations performed.
        logical :: converge_on_fcn
            !! True if the solution converged as a result of a zero-valued
            !! function; else, false.
        logical :: converge_on_chng
            !! True if the solution converged as a result of no appreciable
            !! change in solution points between iterations; else, false.
        logical :: converge_on_zero_diff
            !! True if the solution appears to have settled on a stationary
            !! point such that the gradient of the function is zero-valued; 
            !! else, false.
        type(convergence_history), private :: m_history
            !! Convergence history.
        integer, private :: m_history_count = 0
            !! Counter for the actual number of items in the m_history arrays.
    contains
        procedure, public :: initialize => ib_init
        procedure, private :: is_history_allocated => ib_is_history_allocated
        procedure, private :: allocate_history_arrays => ib_allocate_history
        procedure, private :: get_history_buffer_size => ib_get_history_buffer_size
        procedure, private :: reallocate_history_arrays => ib_reallocate_history_storage
        procedure, public :: log_iteration => ib_log_iteration
        procedure, public :: get_convergence_history => ib_get_history
    end type

    interface iteration_behavior
        module procedure :: ib_constructor
    end interface

    type value_pair
        !! Defines a pair of numeric values.
        real(real64) :: x1
            !! Value 1.
        real(real64) :: x2
            !! Value 2.
    end type

contains
! ******************************************************************************
! ITERATION_BEHAVIOR PROCEDURES
! ------------------------------------------------------------------------------
    subroutine ib_init(this, max_iters)
        !! Initializes the [[iteration_behavior]] object.  All fields are reset
        !! such that all counters are zero and all convergence flags are false.
        class(iteration_behavior), intent(inout) :: this
            !! The [[iteration_behavior]] object.
        integer, intent(in) :: max_iters
            !! The maximum number of iterations the solver is allowed to take.
            !! An NL_INVALID_INPUT_ERROR is thrown if this value is less than
            !! one.

        if (max_iters < 1) error stop NL_INVALID_INPUT_ERROR

        this%converge_on_chng = .false.
        this%converge_on_fcn = .false.
        this%converge_on_zero_diff = .false.
        this%fcn_count = 0
        this%gradient_count = 0
        this%iter_count = 0
        this%jacobian_count = 0

        this%m_history_count = 0
        call this%allocate_history_arrays(max_iters)
    end subroutine

! ------------------------------------------------------------------------------
    function ib_constructor(max_iters) result(rst)
        !! Initializes the [[iteration_behavior]] object.  All fields are set
        !! such that all counters are zero and all convergence flags are false.
        integer, intent(in) :: max_iters
            !! The maximum number of iterations the solver is allowed to take.
            !! An NL_INVALID_INPUT_ERROR is thrown if this value is less than
            !! one.
        type(iteration_behavior) :: rst
            !! The [[iteration_behavior]] object.

        call rst%initialize(max_iters)
    end function

! ------------------------------------------------------------------------------
    pure function ib_is_history_allocated(this) result(rst)
        !! Tests to see if the m_history arrays are allocated.
        class(iteration_behavior), intent(in) :: this
            !! The [[iteration_behavior]] object.
        logical :: rst
            !! True if all of the arrays are allocated; else, false.

        rst = allocated(this%m_history%log_gradient) .and. &
            allocated(this%m_history%log_residual) .and. &
            allocated(this%m_history%log_step_change)
    end function

! ------------------------------------------------------------------------------
    subroutine ib_allocate_history(this, n)
        !! Allocates the history array to the requested size.
        class(iteration_behavior), intent(inout) :: this
            !! The [[iteration_behavior]] object.
        integer, intent(in) :: n
            !! The size.

        real(real64) :: nan
        nan = ieee_value(nan, ieee_quiet_nan)

        if (allocated(this%m_history%log_gradient)) deallocate(this%m_history%log_gradient)
        if (allocated(this%m_history%log_residual)) deallocate(this%m_history%log_residual)
        if (allocated(this%m_history%log_step_change)) deallocate(this%m_history%log_step_change)
        allocate( &
            this%m_history%log_gradient(n), &
            this%m_history%log_residual(n), &
            this%m_history%log_step_change(n), &
            source = nan &
        )
    end subroutine

! ------------------------------------------------------------------------------
    pure function ib_get_history_buffer_size(this) result(rst)
        !! Gets the size of the history buffers.
        class(iteration_behavior), intent(in) :: this
            !! The [[iteration_behavior]] object.
        integer :: rst
            !! The size.

        if (.not.this%is_history_allocated()) then
            rst = 0
        else
            rst = min( &
                size(this%m_history%log_gradient), &
                size(this%m_history%log_residual), &
                size(this%m_history%log_step_change) &
            )
        end if
    end function

! ------------------------------------------------------------------------------
    subroutine ib_reallocate_history_storage(this, new_size)
        !! Reallocates storage for the history buffers.
        class(iteration_behavior), intent(inout) :: this
            !! The [[iteration_behavior]] object.
        integer :: new_size
            !! The new size of the storage arrays.

        integer :: n
        real(real64), allocatable, dimension(:) :: g, r, x

        if (.not.this%is_history_allocated()) then
            ! Simply allocate as nothing is allocated yet
            call this%allocate_history_arrays(new_size)
        else
            ! Resize each array
            call move_alloc(this%m_history%log_gradient, g)
            call move_alloc(this%m_history%log_residual, r)
            call move_alloc(this%m_history%log_step_change, x)
            call this%allocate_history_arrays(new_size)
            n = size(g)
            this%m_history%log_gradient(1:n) = g
            this%m_history%log_residual(1:n) = r
            this%m_history%log_step_change(1:n) = x
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine ib_log_iteration(this, resid, dx, grad)
        !! Registers an iteration given the appropriate metrics.
        class(iteration_behavior), intent(inout) :: this
            !! The [[iteration_behavior]] object.
        real(real64), intent(in) :: resid
            !! The residual value.
        real(real64), intent(in) :: dx
            !! The norm of the solution direction vector.
        real(real64), intent(in) :: grad
            !! The norm of the gradient vector.

        this%m_history_count = this%m_history_count + 1
        if (this%m_history_count > this%get_history_buffer_size()) then
            call this%reallocate_history_arrays(this%m_history_count + 100)
        end if

        this%m_history%log_residual(this%m_history_count) = log10(resid)
        this%m_history%log_step_change(this%m_history_count) = log10(dx)
        this%m_history%log_gradient(this%m_history_count) = log10(grad)
    end subroutine

! ------------------------------------------------------------------------------
    pure function ib_get_history(this) result(rst)
        !! Gets the convergence history.
        class(iteration_behavior), intent(in) :: this
            !! The [[iteration_behavior]] object.
        type(convergence_history) :: rst
            !! The [[convergence_history]] object containing convergence 
            !! tracking information.

        integer :: n

        if (this%m_history_count == 0) then
            allocate(rst%log_gradient(0), rst%log_residual(0), rst%log_step_change(0))
        else
            n = this%m_history_count
            allocate(rst%log_gradient(n), source = this%m_history%log_gradient(1:n))
            allocate(rst%log_residual(n), source = this%m_history%log_residual(1:n))
            allocate(rst%log_step_change(n), source = this%m_history%log_step_change(1:n))
        end if
    end function

! ------------------------------------------------------------------------------
end module