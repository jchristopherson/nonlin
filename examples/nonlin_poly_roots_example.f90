! nonlin_poly_roots_example.f90

program example
    use iso_fortran_env
    use nonlin
    implicit none

    ! Local Variables
    type(polynomial) :: f
    real(real64) :: coeffs(4)
    complex(real64), allocatable :: rts(:)
    integer(int32) :: i

    ! Define the polynomial (x**3 - 2 * x - 1)
    coeffs = [-1.0d0, -2.0d0, 0.0d0, 1.0d0]
    f = coeffs

    ! Compute the polynomial roots
    rts = f%roots()

    ! Display the results
    do i = 1, size(rts)
        print '(AI0AF9.6AF9.6A)', "Root ", i, " = (", real(rts(i), real64), &
            ", ", aimag(rts(i)), ")"
    end do
end program
