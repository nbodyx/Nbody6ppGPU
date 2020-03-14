      PROGRAM MWPotential

      implicit none
      include 'MWpotential.h'
      REAL*8 M_PowCut,Pot_PowCut,Pot_Miyamoto,Pot_NFW
      REAL*8 rxy,xi(3),vi(3),theta,ran2,pi
      REAL*8 B_M,B_F(3),B_FD(3),B_pot
      REAL*8 D_F(3),D_FD(3),D_pot
      REAL*8 H_F(3),H_FD(3),H_pot
      REAL*8 F(3),FD(3),F_SCALE
      REAL*8 rmin,rmax,dr,ri,fr,pot
      INTEGER nbin,I,idum

      write(6,*) 'Input: Rmin, Rmax, Nbin, RandSeed'
      READ(5,*) Rmin,Rmax,nbin,idum

      pi = 3.1415926535897932D0

*     Bulge 
      B_AMP = 1.0 
      B_ALPHA = 1.8D0
      B_RC = 1.9D0          ! kpc
      B_NORM = 0.05D0


*     Disk
      D_AMP = 1.0
      D_A = 3.0D0         ! kpc
      D_B = 0.28D0        ! kpc
      D_NORM = 0.6D0

*     Halo
      H_AMP = 1.0
      H_A = 16.0D0        ! kpc
      H_NORM = 0.35D0

*     scale
      F_Scale = 6.32793804994D0  ! pc/Myr^2
*     Potential unit: kpc*pc/Myr^2

*     Normalization
      xi(1)=8.0D0  ! kpc
      xi(2)=0.0
      xi(3)=0.0
      vi(1:3) = 0
      F(1:3) = 0
      call F_PowCut(xi,vi,F,FD)
      FR = sqrt(F(1)*F(1)+F(2)*F(2))
      B_AMP = B_NORM/FR*F_SCALE
      call F_Miyamoto(xi,vi,F,FD)
      FR = sqrt(F(1)*F(1)+F(2)*F(2))
      D_AMP = D_NORM/FR*F_SCALE
      call F_NFW(xi,vi,F,FD)
      FR = sqrt(F(1)*F(1)+F(2)*F(2))
      H_AMP = H_NORM/FR*F_SCALE
      
      write(6,*) 'Normalize ',B_AMP, D_AMP, H_AMP

*     loop for calculation
      dr = (Rmax/Rmin)**(1.0/REAL(nbin))

      ri = Rmin
      DO I = 1, nbin
         xi(3) = (1.0 - 2.0*ran2(idum))*ri
         theta = 2.0*pi*ran2(idum)
         rxy = sqrt(ri*ri-xi(3)*xi(3))
         xi(1) = rxy*cos(theta)
         xi(2) = rxy*sin(theta)
         vi(1:3) = 0

         B_M = M_PowCut(RI,B_AMP,B_ALPHA,B_RC)
         B_pot = Pot_PowCut(RI,B_AMP,B_ALPHA,B_RC)
         B_F(1:3) = 0.0
         B_FD(1:3) = 0.0
         call F_PowCut(xi,vi,B_F,B_FD)

         D_pot = Pot_Miyamoto(rxy,xi(3),D_AMP,D_A,D_B)
         D_F(1:3) = 0.0
         D_FD(1:3) = 0.0
         call F_Miyamoto(xi,vi,D_F,D_FD)

         H_pot = Pot_NFW(ri, H_AMP, H_A)
         H_F(1:3) = 0.0
         H_FD(1:3) = 0.0
         call F_NFW(xi,vi,H_F,H_FD)

*         fr = (F(1)*xi(1)+F(2)*xi(2))/rxy
         WRITE(10,*) xi,vi,rxy,B_M,B_pot,B_F,B_FD,D_pot,D_F,D_FD,
     &        H_pot,H_F,H_FD
         ri = ri*dr
      END DO 

      write(6,*) 'Write to fort.10'

      RETURN

      END

      
