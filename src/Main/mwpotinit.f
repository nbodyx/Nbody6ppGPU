      SUBROUTINE mwpotinit(RSCALE, VSCALE, TSCALE)
*     
*     
*     initial MWpotential 
*     ------------------------------
*     
      implicit none
      INCLUDE 'MWpotential.h'
      REAL*8 F_MW_UNIT, XI(3), VI(3), FI(3),FDI(3), FRI
      REAL*8 RSCALE, VSCALE, TSCALE

*     Bulge 
      B_AMP = 1.0 
      B_ALPHA = 1.8D0
      B_RC = 1.9D0              ! kpc
      B_NORM = 0.05D0

*     Disk
      D_AMP = 1.0
      D_A = 3.0D0               ! kpc
      D_B = 0.28D0              ! kpc
      D_NORM = 0.6D0

*     Halo
      H_AMP = 1.0
      H_A = 16.0D0              ! kpc
      H_NORM = 0.35D0

*     Normalization
      F_MW_UNIT = 6.32793804994D0 ! pc/Myr^2
      XI(1)=8.0D0               ! kpc
      XI(2)=0.0
      XI(3)=0.0
      VI(1:3)= 0.0
      FDI(1:3)= 0.0
      FI(1:3) = 0.0
      call F_PowCut(xi,vi,FI,FDI)
      FRI = sqrt(FI(1)*FI(1)+FI(2)*FI(2))
      B_AMP = B_NORM/FRI*F_MW_UNIT
      call F_Miyamoto(xi,vi,FI,FDI)
      FRI = sqrt(FI(1)*FI(1)+FI(2)*FI(2))
      D_AMP = D_NORM/FRI*F_MW_UNIT
      call F_NFW(xi,vi,FI,FDI)
      FRI = sqrt(FI(1)*FI(1)+FI(2)*FI(2))
      H_AMP = H_NORM/FRI*F_MW_UNIT

*     scale from NB to Astronomical units
      R_KPC = RSCALE/1000.0D0   ! kpc
      V_PcMyr = VSCALE*1.02269032D0 ! km/s->pc/Myr
      F_PcMyr2 = RSCALE/(TSCALE*TSCALE) ! pc/Myr^2
      FD_PcMyr3= F_PcMyr2/TSCALE*1000.0D0 ! pc/Myr^3
      P_kpcPcMyr2= F_PcMyr2*R_KPC !pc^2/Myr^2


      RETURN

      END

