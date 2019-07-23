*   NFW potential of i (Galpy; Bovy 2015)
*                   amp
*     rho =    ------------
*              r/A (1 + r/A)^2
      SUBROUTINE F_NFW(XI,XIDOT,FM,FD)
*
*
*       ----------------------------------
*
      implicit none
      INCLUDE 'MWpotential.h'
      REAL*8  XI(3),XIDOT(3),FM(3),FD(3)
      REAL*8  r2,DP,r,drdv,ra,r3,ddp
*
*
*       Obtain force and first derivative for logarithmic potential.
      r2 = XI(1)**2 + XI(2)**2 + XI(3)**2
      r = sqrt(r2)
      r3 = r*r2
      ra = H_A + r
      drdv = XI(1)*XIDOT(1) + XI(2)*XIDOT(2) + XI(3)*XIDOT(3)
      DP  = H_AMP*(1.0D0/(r2*ra) - log(1.0D0+r/H_A)/r3)
      DDP = H_AMP*(-2.0D0/(r*ra) - 1.0D0/(ra*ra) + 3*log(1+r/H_A)/r2)
     &     *drdv/r3
*
*     Note softening in the logarithmic potential (cf. Binney & Tremaine).
      FM = DP*XI
      FD = DP*XIDOT - DDP*XI
*
      RETURN
*
      END

*********************
      REAL*8 FUNCTION Pot_NFW(r,amp,a)
**    potential

      implicit none
      REAL*8 r,amp,a

      Pot_NFW = -amp*log(1.0D0 + r/a)/r

      return 

      end

