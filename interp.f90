  REAL FUNCTION cubic_interp(p, w)
    REAL, INTENT(IN) :: w
    REAL, INTENT(IN), DIMENSION(4) :: p
   ! cubic_interp = p(2) + 0.5*w*(p(3) - p(1) + w*(2.0*p(1) - 5.0*p(2) + 4.0*p(3) - p(4) + w*(3.0*(p(2) - p(3)) + p(4) - p(1))))
 cubic_interp=((-0.5*p(1)+1.5*p(2)-1.5*p(3)+0.5*p(4))*(w**3))+((p(1)-2.5*p(2)+2*p(3)-0.5*p(4))*(w**2))+((-0.5*p(1)+0.5*p(3))*w)+p(2)

  END FUNCTION cubic_interp

  REAL FUNCTION bicubic_interp(p, y, x)
    ! Input Variables
    REAL, INTENT(IN) :: y, x
    REAL, INTENT(IN), DIMENSION(4,4) :: p
    ! Local Variables
    REAL, DIMENSION(4) :: p_temp
	
    p_temp(1) = cubic_interp(p(1,:), y)
    p_temp(2) = cubic_interp(p(2,:), y)
    p_temp(3) = cubic_interp(p(3,:), y)
    p_temp(4) = cubic_interp(p(4,:), y)
    bicubic_interp = cubic_interp(p_temp, x)
  END FUNCTION bicubic_interp
	
  REAL FUNCTION linear_interp(x1, x2, xi)
        IMPLICIT NONE
        REAL, INTENT(IN) :: x1, x2, xi
        linear_interp = x1*(1.0-xi) + (x2*xi)
   return
   END FUNCTION linear_interp
 
    REAL FUNCTION neville_interp(psub, prof, inp, n)  !! neville interpolation New Version
    ! Input Variable
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n  			!! number of values in array.
    REAL, INTENT(IN) :: inp
    REAL, INTENT(IN), DIMENSION(n) :: psub, prof
    ! Local Variables
    INTEGER :: i, j
    REAL, DIMENSION(n) :: f
	
    f = prof
     do i = 1,n-1
      do j = 1,(n-i)
         f(j) = ( (inp-psub(j))*f(j+1) - (inp-psub(j+i))*f(j) )  / (psub(j+i)-psub(j) );
      end do
    end do
    neville_interp = f(1)
  END FUNCTION neville_interp


