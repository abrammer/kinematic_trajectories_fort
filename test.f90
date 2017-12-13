    REAL FUNCTION cubic_interpolate(p, w)
    IMPLICIT NONE

    ! Input Variables
    REAL, INTENT(IN) :: w
    REAL, INTENT(IN), DIMENSION(4) :: p

    ! Local Variables
!    REAL :: cubic_interpolate

	print*, p
	print*,w
    cubic_interpolate = p(2) + 0.5*w*(p(3) - p(1) + w*(2.0*p(1) - 5.0*p(2) + 4.0*p(3) - p(4) + w*(3.0*(p(2) - p(3)) + p(4) - p(1))))

  END FUNCTION cubic_interpolate
