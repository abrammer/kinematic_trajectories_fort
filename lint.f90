  FUNCTION abcd(x1, x2, xi,xx)
        IMPLICIT NONE
  	REAL, INTENT(IN) :: x1, x2, xi
  	real :: abcd, xx
        abcd = x1*(1.0-xi) + (x2*xi)
        xx = abcd
         print*,xx
   return
   END FUNCTION abcd

  FUNCTION lbcd(x1, x2, xi,xx)
        IMPLICIT NONE
        REAL, INTENT(IN) :: x1, x2, xi
        real :: lbcd, xx
        lbcd = x1*(1.0-xi) + (x2*xi)
        xx = lbcd
         print*,xx
   return
   END FUNCTION lbcd


