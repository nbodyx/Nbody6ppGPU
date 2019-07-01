*  Bulge: spherical power-law potential w/ cutoff (Galpy; Bovy 2015)
*
*                                     amp
*                          rho(r)= ---------   e^{-(r/rc)^2}
*                                   r^\alpha

      SUBROUTINE F_PowCut(XI,XIDOT,FM,FD)
*
*
*       ----------------------------------
*
      implicit none
      INCLUDE 'MWpotential.h'
      REAL*8  XI(3),XIDOT(3),FM(3),FD(3),r,mr,rvr2,r2,dmr
      REAL*8 M_PowCut,Den_PowCut
*
*
*     Radial distance
      r = SQRT(XI(1)**2 + XI(2)**2 + XI(3)**2)
      r2 = r*r
      rvr2 = (XI(1)*XIDOT(1) + XI(2)*XIDOT(2) + XI(3)*XIDOT(3))/r2
      mr = -M_PowCut(r,B_AMP,B_ALPHA,B_RC,PI)/(r2*r)
      dmr = -4.0*pi*Den_PowCut(r,B_AMP,B_ALPHA,B_RC)

      FM(1:3) = mr*XI(1:3)
      FD = mr*(XIDOT -3.0*rvr2*XI) + dmr*rvr2/r2*XI

*
      RETURN
*
      END
      
*********************
      REAL*8 FUNCTION Den_PowCut(r,amp,alpha,rc)
**    density profile

      REAL*8 r,alpha,rc,amp

      DEN_PowCut = amp/r**alpha*exp(-(r/rc)**2.0)

      return 
      end

*********************
      REAL*8 FUNCTION M_PowCut(r,amp,alpha,rc,pi)
**    cumulative mass

      REAL*8 R,amp,alpha,rc,pi,a,xup,ginc,gammds,gamma
      INTEGER ifault

      a = 1.5-alpha/2.0
      xup = (r/rc)**2.0
      ginc =gammds(xup,a,ifault) 
      M_PowCut = amp*2.0*pi*rc*(3.0-alpha)*ginc*gamma(a)

      return 
      end
      
********************
      REAL*8 FUNCTION Pot_PowCut(r,amp,alpha,rc,pi)
**    potential

      REAL*8 r,alpha,rc,amp,pi,a,a2,xup,ginc,ginc2,gammds
      REAL*8 gamma,g,g2
      INTEGER ifault
      
      xup = (r/rc)**2.0
      a = 1.0D0 - alpha/2.0D0
      a2= 1.5D0 - alpha/2.0D0
      ginc = gammds(xup, a,ifault) 
      g = gamma(a)
      if(ifault>0) write(6,*) 'Fault',a,xup

      ginc2 =gammds(xup,a2,ifault) 
      g2 = gamma(a2)
      if(ifault>0) write(6,*) 'Fault',a2,xup

      Pot_PowCut = amp*2.0*pi*rc**(3.0-alpha)*(ginc*g/rc-ginc2*g2/r)
      return 
      end

