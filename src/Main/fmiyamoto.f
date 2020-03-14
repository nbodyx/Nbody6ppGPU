*   MiyamotoNagaiPotential.py: class that implements the Miyamoto-Nagai (Galpy; Bovy 2015)
*                              potential
*                                                           GM
*                              phi(R,z) = -  ---------------------------------
*                                             \sqrt(R^2+(a+\sqrt(z^2+b^2))^2)
*
      SUBROUTINE F_Miyamoto(XI,XIDOT,FM,FD)
*
*
*       ----------------------------------
*
      implicit none
      INCLUDE 'MWpotential.h'
      REAL*8  XI(3),XIDOT(3),FM(3),FD(3)
      REAL*8  BZ,R2,AZ,AZ3,AZDOT,Y1,Y2
*
*
*       Obtain force & derivative for global variables XI & XIDOT.
      R2 = XI(1)**2 + XI(2)**2
      BZ = SQRT(D_B**2 + XI(3)**2)
      AZ = SQRT(R2 + (D_A + BZ)**2)
      AZ3 = D_AMP/AZ**3
*       Note missing square root sign in (b^2 + z^2) of Book eq. (8.52).
      AZDOT = XI(1)*XIDOT(1) + XI(2)*XIDOT(2) +
     &     (D_A + BZ)*XI(3)*XIDOT(3)/BZ
      FM(1) = -AZ3*XI(1)
      FM(2) = -AZ3*XI(2)
      FM(3) = -AZ3*XI(3)*(D_A + BZ)/BZ
*     RDOT = (XI(1)*XIDOT(1) + XI(2)*XIDOT(2))/SQRT(R2)
      FD(1) = -AZ3*(XIDOT(1) - 3.0*AZDOT*XI(1)/AZ**2)
      FD(2) = -AZ3*(XIDOT(2) - 3.0*AZDOT*XI(2)/AZ**2)
*       Note wrong sign in first term of eq. (8.52) (see printer errata).
      Y1 = 3.0*(D_A + BZ)*XI(3)*AZDOT/(AZ**2*BZ)
*       Note printer mistake of dZ/dt in denominator (see errata).
      Y2 = (D_A*D_B**2 + BZ**3)*XIDOT(3)/BZ**3
      FD(3) = AZ3*(Y1 - Y2)
*       Refer to Miyamoto & Nagai (PASJ 27, 533) and Book eq. (8.52).
*
      RETURN
*
      END

*********************
      REAL*8 FUNCTION Pot_Miyamoto(r,z,amp,a,b)
**    potential

      implicit none
      REAL*8 r,z,amp,a,b

      Pot_Miyamoto= -amp/sqrt(r*r + (a + sqrt(z*z + b*b))**2.0D0)

      return 

      end

