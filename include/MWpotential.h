*     Galpy
*     Bulge spherical power-law potential w/ cutoff
* 
*                                     amp
*                          rho(r)= ---------   e^{-(r/rc)^2}
*                                    r^\alpha
* 
*
*     Tidal force for Milky-Way potential
      COMMON/MWPOT/ B_AMP, B_ALPHA, B_RC, B_NORM,
     &     D_AMP, D_A, D_B, D_NORM,
     &     H_AMP, H_A, H_NORM, 
     &     F_PcMyr2, FD_PcMyr3, V_PcMyr, R_KPC, P_kpcPcMyr2

      REAL*8 B_AMP, B_ALPHA, B_RC, B_NORM
      REAL*8 D_AMP, D_A, D_B, D_NORM  
      REAL*8 H_AMP, H_A, H_NORM
      REAL*8 F_PcMyr2, FD_PcMyr3, V_PcMyr, R_KPC, P_kpcPcMyr2
      
      
      
