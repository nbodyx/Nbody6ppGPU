      PROGRAM MWPotential

      implicit none
      include 'MWpotential.h'
      REAL*8 M_PowCut,Pot_PowCut,ran2
      REAL*8 rxy,xi(3),vi(3),theta
      REAL*8 M,F(3),FD(3)
      REAL*8 rmin,rmax,dr,ri,fr,pot
      INTEGER nbin,I,idum

      B_AMP = 1.0D0
      B_ALPHA = 1.8D0
      B_RC = 1.9D0/8.0D0
      B_NORM = 0.05D0
      PI = 4.0D0*ATAN(1.0D0)

      write(6,*) 'Input: Rmin, Rmax, Nbin, RandSeed'
      READ(5,*) Rmin,Rmax,nbin,idum

      xi(1)=1.0
      xi(2)=0.0
      xi(3)=0.0
      vi(1:3) = 0
      F(1:3) = 0
      call F_PowCut(xi,vi,F,FD)
      B_AMP = B_AMP*B_NORM/sqrt(F(1)*F(1)+F(2)*F(2))

      write(6,*) 'Normalize ',B_AMP

      pi = 4.0D0*ATAN(1.0D0)
      
      dr = (Rmax/Rmin)**(1.0/REAL(nbin))

      ri = Rmin
      DO I = 1, nbin
         xi(3) = (1.0 - 2.0*ran2(idum))*ri
         theta = 2.0*pi*ran2(idum)
         rxy = sqrt(ri*ri-xi(3)*xi(3))
         xi(1) = rxy*cos(theta)
         xi(2) = rxy*sin(theta)
         vi(1:3) = 0

         M = M_PowCut(RI,B_AMP,B_ALPHA,B_RC,PI)
         pot = Pot_PowCut(RI,B_AMP,B_ALPHA,B_RC,PI)
         F(1:3) = 0.0
         FD(1:3) = 0.0
         call F_PowCut(xi,vi,F,FD)
         fr = (F(1)*xi(1)+F(2)*xi(2))/rxy
         WRITE(10,*) xi,vi,rxy,M,pot,F,FD,fr
         ri = ri*dr
      END DO 

      write(6,*) 'Write to fort.10'

      RETURN

      END

      
