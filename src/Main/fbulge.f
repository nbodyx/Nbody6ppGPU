      SUBROUTINE FBULGE(XI,XIDOT,FM,FD)
*
*
*       Gamma/eta potential for the bulge.
*       ----------------------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      INCLUDE 'galaxy.h'
      REAL*8  XI(3),XIDOT(3),FM(3),FD(3)
*
*
*     Obtain force and first derivative for gamma/eta model (Dehnen 1993).
*     
*     phi(r) = -GM/[a(2-g)] * [ 1 - (1+a/r)**(g-2) ]  (g # 2)
*     phi(r) = -GM/a*ln(1 + a/r)                      (g = 2)
*
      RR = SQRT(XI(1)**2 + XI(2)**2 + XI(3)**2)
      RRD = (XI(1)*XIDOT(1) + XI(2)*XIDOT(2) + XI(3)*XIDOT(3))/(RR*RR)
      RRD = RRD*(AR*GAM + 3.0*RR)/(RR + AR)
      H3 = (GMB/RR**3)*(1.0 + AR/RR)**(GAM-3.0)
*
      DO 10 K = 1,3
         FM(K) = -H3*XI(K)
         FD(K) = -H3*(XIDOT(K) - RRD*XI(K))
 10   CONTINUE
*
      RETURN
*
      END
