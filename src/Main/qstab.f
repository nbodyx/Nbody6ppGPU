      REAL*8 FUNCTION QSTAB(e,e_out,zi,m1,m2,m3)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 m1,m2,m3
*
*     Three-body stability function (Valtonen, AAS 2015).
*     ---------------------------------------------------
*
*     Inner and outer eccentricity: e & e_out.
*     Inclination in radians: zi.
*     Masses: inner, m1 & m2, outer, m3; m1 > m2.
*
*     Adopt 10,000 outer orbits for random walk time-scale.
      ZN = 10000.0
      zz = 1.0/6.0D0
*
      F = 1.0 - 2.0*e/3.0*(1.0 - 0.5*e**2) - 0.3*cos(zi)*
     &   (1.0 - 0.5*e + 2.0*cos(zi)*(1.0 - 2.5*e**1.5 - cos(zi)))
*
      G = SQRT(MAX(m1,m2) /(m1 + m2))*(1.0 + m3/(m1 + m2))
*
      QSTAB = 1.52*(SQRT(ZN)/(1.0 - e_out))**zz*(F*G)**0.3333
*
      RETURN
      END
