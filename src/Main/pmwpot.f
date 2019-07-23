      SUBROUTINE pmwpot(RG, POTG)
*     
*     
*     Calculate external force of MW pontetial for one star
*     --------------------------
      implicit none
      INCLUDE 'MWpotential.h'
      REAL*8 RG(3), VG(3), POTG
      REAL*8 Pot_PowCut, Pot_Miyamoto, Pot_NFW
      REAL*8 RGKPC(3),rxy2, rxy, r

*     convert to astronomical unit
      RGKPC = RG*R_KPC

      rxy2  = RGKPC(1)*RGKPC(1) + RGKPC(2)*RGKPC(2)
      rxy = SQRT(rxy2)
      r   = SQRT(rxy2 + RGKPC(3)*RGKPC(3))

      POTG = 0.0D0
*     bulge
      POTG = POTG + Pot_PowCut(r,B_AMP,B_ALPHA,B_RC)
*     disk
      POTG = POTG + Pot_Miyamoto(rxy,RGKPC(3),D_AMP, D_A, D_B)
*     halo
      POTG = POTG + Pot_NFW(r, H_AMP, H_A)

*     convert to NB unit
      POTG = POTG/P_kpcPcMyr2

      RETURN

      END
      
