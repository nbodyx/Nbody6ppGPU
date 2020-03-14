      SUBROUTINE fmwpot(RG, VG, FG, FGD)
*     
*     
*     Calculate external force of MW pontetial for one star
*     --------------------------
      implicit none
      INCLUDE 'MWpotential.h'
      REAL*8 RG(3), VG(3), FG(3), FGD(3)
      REAL*8 RGKPC(3),VGPCMYR(3),FS(3),FSD(3)

*     convert to astronomical unit
      RGKPC = RG*R_KPC
      VGPCMYR=VG*V_PcMyr

      FG = 0
      FGD =0
*     bulge
      call F_PowCut(RGKPC,VGPCMYR,FS,FSD)
      FG = FG  + FS
      FGD= FGD + FSD
*     disk
      call F_Miyamoto(RGKPC,VGPCMYR,FS,FSD)
      FG = FG  + FS
      FGD= FGD + FSD
*     halo
      call F_NFW(RGKPC,VGPCMYR,FS,FSD)
      FG = FG  + FS
      FGD= FGD + FSD

*     debug
*      write(111,*) RGKPC, VGPCMYR, FG, FGD

*     convert to NB unit
      FG = FG/F_PcMyr2
*      FGD =0.0
      FGD= FGD/FD_PcMyr3


      RETURN

      END
      
