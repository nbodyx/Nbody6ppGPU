      REAL*8 FUNCTION TSTEP2(XDOT,F,FDOT,ETA)
*
*
*       General time-step criterion by xdot, f and fdot.
*       ----------------------------
*
      IMPLICIT NONE
      REAL*8  XDOT(3),F(3),FDOT(3)
      REAL*8 F2,FDOT2,V2,ETA
*
*       Obtain new integration step using composite expression.
*       STEP = (ETA*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2))**0.5.
*
      V2 = XDOT(1)**2 + XDOT(2)**2 + XDOT(3)**2
      F2 = F(1)**2 + F(2)**2 + F(3)**2
      FDOT2 = FDOT(1)**2 + FDOT(2)**2 + FDOT(3)**2
*
      TSTEP2 = 0.25*ETA*SQRT(MIN(V2/F2,F2/FDOT2))
*
      RETURN
*
      END
