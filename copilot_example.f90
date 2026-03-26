module dls_solver
   !! Damped nonlinear least-squares with simple bound constraints
   !!
   !! Minimize  0.5 * || r(x) ||_2^2
   !! subject to  lb <= x <= ub
   !!
   !! Levenberg–Marquardt / damped Gauss–Newton with box projection.

   implicit none
   private

   !=================================================================
   ! Public types
   !=================================================================
   type :: dls_options
      real(real64) :: tol_grad   = 1.0e-8_real64
      real(real64) :: tol_step   = 1.0e-10_real64
      real(real64) :: tol_f      = 1.0e-12_real64
      integer       :: max_iter  = 100
      real(real64) :: lambda0    = 1.0e-3_real64
      real(real64) :: lambda_inc = 10.0_real64
      real(real64) :: lambda_dec = 0.3_real64
      logical       :: verbose   = .false.
   end type dls_options

   type :: dls_status
      logical       :: converged = .false.
      integer       :: iter      = 0
      real(real64) :: fnorm      = huge(1.0_real64)
      integer       :: info      = 0   !! 0=OK, >0=stopped, <0=error
   end type dls_status

   public :: dls_options, dls_status, dls_solve

   !=================================================================
   ! Procedure interfaces for user callbacks
   !=================================================================
   interface
      subroutine dls_residual(x, r)
         import :: real64
         real(real64), intent(in)  :: x(:)
         real(real64), intent(out) :: r(:)
      end subroutine dls_residual

      subroutine dls_jacobian(x, J)
         import :: real64
         real(real64), intent(in)  :: x(:)
         real(real64), intent(out) :: J(:,:)
      end subroutine dls_jacobian
   end interface

contains

   !=================================================================
   ! Main solver
   !=================================================================
   subroutine dls_solve(x, residual, jacobian, lb, ub, opts, stat)
      use iso_fortran_env, only: real64
      implicit none

      !------------------------------
      ! Arguments
      !------------------------------
      real(real64), intent(inout)          :: x(:)
      procedure(dls_residual)              :: residual
      procedure(dls_jacobian), optional    :: jacobian
      real(real64), intent(in),  optional  :: lb(:), ub(:)
      type(dls_options), intent(in), optional :: opts
      type(dls_status),  intent(out), optional :: stat

      !------------------------------
      ! Local variables
      !------------------------------
      type(dls_options) :: o
      type(dls_status)  :: s

      integer           :: n, m, k
      real(real64), allocatable :: r(:), r_new(:)
      real(real64), allocatable :: J(:,:), JTJ(:,:), g(:), p(:), x_new(:)
      real(real64) :: fnorm, fnorm_new, lambda, rho, pred_red, act_red
      real(real64) :: grad_inf, step_inf
      logical      :: have_bounds

      ! For finite-difference Jacobian
      real(real64) :: eps_fd
      integer      :: i

      !------------------------------
      ! Initialization
      !------------------------------
      if (present(opts)) then
         o = opts
      else
         o = dls_options()
      end if

      n = size(x)
      call residual(x, r = [0.0_real64])  ! dummy to get size? Not allowed.
      ! Better: require user to ensure residual size via a workspace query
      ! but for simplicity we do a one-time call with temporary alloc.

      ! Determine residual dimension
      allocate(r(1))
      call residual(x, r)
      m = size(r)
      deallocate(r)

      allocate(r(m), r_new(m))
      allocate(J(m,n), JTJ(n,n), g(n), p(n), x_new(n))

      eps_fd = sqrt(epsilon(1.0_real64))
      lambda = o.lambda0

      ! Initial evaluation
      call residual(x, r)
      fnorm = norm2(r)

      have_bounds = present(lb) .or. present(ub)

      if (o.verbose) then
         write(*,'(a,1x,i4,1x,es12.4)') 'Iter  F-norm:', 0, fnorm
      end if

      !------------------------------
      ! Main iteration loop
      !------------------------------
      do k = 1, o.max_iter

         ! Jacobian: analytic if provided, else finite differences
         if (present(jacobian)) then
            call jacobian(x, J)
         else
            call fd_jacobian(residual, x, r, J, eps_fd)
         end if

         ! Compute JTJ and gradient g = J^T r
         call compute_jtj_g(J, r, JTJ, g)

         grad_inf = maxval(abs(g))

         ! Check gradient-based stopping
         if (grad_inf < o.tol_grad) then
            s%converged = .true.
            s%iter      = k-1
            s%fnorm     = fnorm
            s%info      = 1
            exit
         end if

         ! Solve (JTJ + lambda I) p = -g
         call damped_step(JTJ, g, lambda, p)

         ! Trial step and projection onto bounds
         x_new = x + p
         if (have_bounds) call project_bounds(x_new, lb, ub)

         step_inf = maxval(abs(x_new - x))
         if (step_inf < o.tol_step) then
            s%converged = .true.
            s%iter      = k-1
            s%fnorm     = fnorm
            s%info      = 2
            exit
         end if

         ! Evaluate new residual and objective
         call residual(x_new, r_new)
         fnorm_new = norm2(r_new)

         act_red = 0.5_real64*(fnorm**2 - fnorm_new**2)

         ! Predicted reduction: 0.5 * p^T (lambda p - g)
         pred_red = 0.5_real64 * dot_product(p, lambda*p - g)

         if (pred_red <= 0.0_real64) then
            ! Numerical issue; increase damping and retry
            lambda = lambda * o.lambda_inc
            cycle
         end if

         rho = act_red / pred_red

         if (rho > 0.0_real64) then
            ! Accept step
            x     = x_new
            r     = r_new
            fnorm = fnorm_new

            ! Decrease damping
            lambda = max(lambda * o.lambda_dec, 1.0e-16_real64)
         else
            ! Reject step, increase damping
            lambda = lambda * o.lambda_inc
         end if

         if (o.verbose) then
            write(*,'(a,1x,i4,1x,es12.4,1x,es10.3,1x,es10.3)') &
                 'Iter F-norm rho lambda:', k, fnorm, rho, lambda
         end if

         ! Function-based stopping
         if (abs(act_red) < o.tol_f) then
            s%converged = .true.
            s%iter      = k
            s%fnorm     = fnorm
            s%info      = 3
            exit
         end if

      end do

      if (.not. s%converged) then
         s%iter  = o.max_iter
         s%fnorm = fnorm
         s%info  = 10   ! max iter
      end if

      if (present(stat)) stat = s

      ! Clean up
      deallocate(r, r_new, J, JTJ, g, p, x_new)

   end subroutine dls_solve

   !=================================================================
   ! Finite-difference Jacobian
   !=================================================================
   subroutine fd_jacobian(residual, x, r_base, J, eps_fd)
      use iso_fortran_env, only: real64
      implicit none
      procedure(dls_residual) :: residual
      real(real64), intent(in)  :: x(:)
      real(real64), intent(out) :: r_base(:)
      real(real64), intent(out) :: J(:,:)
      real(real64), intent(in)  :: eps_fd

      integer :: n, m, j
      real(real64), allocatable :: x_pert(:), r_pert(:)
      real(real64) :: h

      n = size(x)
      call residual(x, r_base)
      m = size(r_base)

      if (size(J,1) /= m .or. size(J,2) /= n) then
         stop 'fd_jacobian: dimension mismatch in J'
      end if

      allocate(x_pert(n), r_pert(m))

      x_pert = x
      do j = 1, n
         h = eps_fd * max(1.0_real64, abs(x(j)))
         x_pert(j) = x(j) + h
         call residual(x_pert, r_pert)
         J(:,j) = (r_pert - r_base) / h
         x_pert(j) = x(j)
      end do

      deallocate(x_pert, r_pert)
   end subroutine fd_jacobian

   !=================================================================
   ! Compute JTJ and gradient g = J^T r
   !=================================================================
   subroutine compute_jtj_g(J, r, JTJ, g)
      use iso_fortran_env, only: real64
      implicit none
      real(real64), intent(in)  :: J(:,:)
      real(real64), intent(in)  :: r(:)
      real(real64), intent(out) :: JTJ(:,:)
      real(real64), intent(out) :: g(:)

      integer :: m, n, i, j, k

      m = size(J,1)
      n = size(J,2)

      if (size(JTJ,1) /= n .or. size(JTJ,2) /= n) stop 'compute_jtj_g: JTJ size'
      if (size(g) /= n) stop 'compute_jtj_g: g size'

      JTJ = 0.0_real64
      g   = 0.0_real64

      do i = 1, m
         do j = 1, n
            g(j) = g(j) + J(i,j)*r(i)
         end do
      end do

      do j = 1, n
         do k = j, n
            do i = 1, m
               JTJ(j,k) = JTJ(j,k) + J(i,j)*J(i,k)
            end do
            JTJ(k,j) = JTJ(j,k)
         end do
      end do
   end subroutine compute_jtj_g

   !=================================================================
   ! Solve (JTJ + lambda I) p = -g via Cholesky
   !=================================================================
   subroutine damped_step(JTJ, g, lambda, p)
      use iso_fortran_env, only: real64
      implicit none
      real(real64), intent(in)  :: JTJ(:,:)
      real(real64), intent(in)  :: g(:)
      real(real64), intent(in)  :: lambda
      real(real64), intent(out) :: p(:)

      integer :: n, i, j, k, info
      real(real64), allocatable :: A(:,:), L(:,:), y(:)

      n = size(g)
      if (size(JTJ,1) /= n .or. size(JTJ,2) /= n) stop 'damped_step: JTJ size'
      if (size(p) /= n) stop 'damped_step: p size'

      allocate(A(n,n), L(n,n), y(n))

      A = JTJ
      do i = 1, n
         A(i,i) = A(i,i) + lambda
      end do

      ! Cholesky factorization A = L L^T
      L = 0.0_real64
      info = 0
      do j = 1, n
         do i = j, n
            L(i,j) = A(i,j)
            do k = 1, j-1
               L(i,j) = L(i,j) - L(i,k)*L(j,k)
            end do
            if (i == j) then
               if (L(i,j) <= 0.0_real64) then
                  info = i
                  exit
               end if
               L(i,j) = sqrt(L(i,j))
            else
               L(i,j) = L(i,j) / L(j,j)
            end if
         end do
         if (info /= 0) exit
      end do

      if (info /= 0) then
         ! Fallback: steepest descent step
         p = -g
         deallocate(A, L, y)
         return
      end if

      ! Solve L y = -g
      y = -g
      do i = 1, n
         do k = 1, i-1
            y(i) = y(i) - L(i,k)*y(k)
         end do
         y(i) = y(i) / L(i,i)
      end do

      ! Solve L^T p = y
      p = y
      do i = n, 1, -1
         do k = i+1, n
            p(i) = p(i) - L(k,i)*p(k)
         end do
         p(i) = p(i) / L(i,i)
      end do

      deallocate(A, L, y)
   end subroutine damped_step

   !=================================================================
   ! Projection onto bound constraints
   !=================================================================
   subroutine project_bounds(x, lb, ub)
      use iso_fortran_env, only: real64
      implicit none
      real(real64), intent(inout)         :: x(:)
      real(real64), intent(in), optional  :: lb(:), ub(:)

      integer :: n, i

      n = size(x)

      if (present(lb)) then
         if (size(lb) /= n) stop 'project_bounds: lb size mismatch'
         do i = 1, n
            if (x(i) < lb(i)) x(i) = lb(i)
         end do
      end if

      if (present(ub)) then
         if (size(ub) /= n) stop 'project_bounds: ub size mismatch'
         do i = 1, n
            if (x(i) > ub(i)) x(i) = ub(i)
         end do
      end if
   end subroutine project_bounds

   !=================================================================
   ! 2-norm helper (Fortran 2008 has norm2, but we keep a local one)
   !=================================================================
   pure function norm2(v) result(nrm)
      use iso_fortran_env, only: real64
      implicit none
      real(real64), intent(in) :: v(:)
      real(real64)             :: nrm
      integer :: i

      nrm = 0.0_real64
      do i = 1, size(v)
         nrm = nrm + v(i)*v(i)
      end do
      nrm = sqrt(nrm)
   end function norm2

end module dls_solver