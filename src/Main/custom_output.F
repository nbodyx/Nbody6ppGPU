      subroutine custom_output(NXTLEN,NXTLST,NXTLIMIT,NGHOSTS)
*
*
*     Custom output with HDF5/binary format  (Long Wang)
*-----------------------------------------------------------------------------------
*     IF(KZ(19).EQ.0.AND.KZ(12).EQ.-1) The output is not scaled
*     single: NAME, M[M*], X(1:3)[pc], V(1:3)[km/s], POT[NB],
*             RS[R*], L[L*], Teff[K], MCORE[M*], RSCORE[R*], K*

*     binary: NAME1, NAME2, NAME(ICM), M1[M*], M2[M*], XCM(1:3)[pc], 
*             VCM(1:3)[km/s], XREL(1:3)[AU], VREL(1:3)[km/s], POT[NB],
*             semi[AU], ecc, P[days], Gamma,
*             RS1[R*], RS2[R*], L1[L*], L2[L*], Teff1[K], Teff2[K],
*             MCORE1[M*], MCORE2[M*], RSCORE1[R*], RSCORE2[R*],
*             K*1, K*2, K*(ICM)

*     merger: NAME1,NAME2,NAME3,NAME(INCM), M1[M*], M2[M*], M3[*],
*             XCM(1:3)[pc], VCM(1:3)[km/s], XREL0(1:3)[AU], VREL0(1:3)[km/s],
*             XREL1(1:3)[AU], VREL1(1:3)[km/s],POT[NB],
*             semi0[AU], ecc0, P0[days], semi1[AU], ecc1, P1[days],
*             RS1[R*], RS2[R*], RS3[R*], L1[L*], L2[L*], L3[L*], Teff1[K], Teff2[K], Teff3[K],
*             MCORE1[M*], MCORE2[M*], MCORE3[M*], RSCORE1[R*], RSCORE2[R*], RSCORE3[R*]
*             K*1, K*2, K*3, K*(INCM), 
*     (binary parameter: 0 means inner binary, 1 means outer component )
*     (single parameter: 1/2, means inner binary components, 3 means outer component )
*     (merger parameter: 0 inner binary, 1 outer binary)
*     (POT: potential, sum of - G m_j/r_j)
*-----------------------------------------------------------------------------------
*

      include 'common6.h'
      include 'timing.h'
      include 'galaxy.h'
      include 'output.h'
      include 'scale_out.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2,GALIMIT=0.1D0)

      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CMC(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
      REAL*8  M,MASS,MC,MIJ,MKK
      REAL*8 TMASK(NMAX),VIREL(3),XIREL(3)
      INTEGER NLOOP

      
C*     single
C      COMMON/OUTSINGLE/ S_M(NMAX),S_X1(NMAX), S_X2(NMAX), S_X3(NMAX),
C     &                  S_V1(NMAX), S_V2(NMAX), S_V3(NMAX),
CC     &                  S_F1(NMAX), S_F2(NMAX), S_F3(NMAX),
C     &                  S_RS(NMAX), S_L(NMAX), S_TE(NMAX),
C     &                  S_RC(NMAX), S_MC(NMAX),
C     &                  NS_KW(NMAX),NS_NAM(NMAX)
C      REAL*4 S_M, S_X1, S_X2, S_X3, S_V1, S_V2, S_V3
CC      REAL*4 S_F1, S_F2, S_F3
C      REAL*4 S_RS, S_L, S_TE, S_RC, S_MC
C      INTEGER NS_KW,NS_NAM
C*     Binary
C      COMMON/OUTBINARY/ B_M1(KMAX), B_M2(KMAX),
C     &                  B_XC1(KMAX), B_XC2(KMAX), B_XC3(KMAX),
C     &                  B_VC1(KMAX), B_VC2(KMAX), B_VC3(KMAX), 
C     &                  B_XR1(KMAX), B_XR2(KMAX), B_XR3(KMAX), 
C     &                  B_VR1(KMAX), B_VR2(KMAX), B_VR3(KMAX), 
CC     &                  B_FC1(KMAX), B_FC2(KMAX), B_FC3(KMAX), 
C     &                  B_RS1(KMAX), B_L1(KMAX), B_TE1(KMAX),
C     &                  B_RS2(KMAX), B_L2(KMAX), B_TE2(KMAX),
C     &                  B_RC1(KMAX), B_MC1(KMAX), B_RC2(KMAX), 
C     &                  B_MC2(KMAX), B_A(KMAX), B_ECC(KMAX), B_P(KMAX),
C     &                  NB_KW1(KMAX), NB_NAM1(KMAX), NB_KW2(KMAX), 
C     &                  NB_NAM2(KMAX), NB_KWC(KMAX), NB_NAMC(KMAX)
C      REAL*4 B_M1, B_M2, B_XC1, B_XC2, B_XC3, B_VC1, B_VC2, B_VC3
C      REAL*4 B_XR1, B_XR2, B_XR3, B_VR1, B_VR2, B_VR3
CC      REAL*4 B_FC1, B_FC2, B_FC3
C      REAL*4 B_RS1, B_L1, B_TE1, B_RS2, B_L2, B_TE2
C      REAL*4 B_RC1, B_MC1, B_RC2, B_MC2,B_A, B_ECC, B_P
C      INTEGER NB_KW1, NB_NAM1, NB_KW2, NB_NAM2, NB_KWC, NB_NAMC
C*     Merger
C      COMMON/OUTMERGER/ M_M1(MMAX), M_M2(MMAX), M_M3(MMAX),
C     &                  M_XC1(MMAX), M_XC2(MMAX), M_XC3(MMAX),
C     &                  M_VC1(MMAX), M_VC2(MMAX), M_VC3(MMAX),
C     &                  M_XR01(MMAX), M_XR02(MMAX), M_XR03(MMAX), 
C     &                  M_VR01(MMAX), M_VR02(MMAX), M_VR03(MMAX), 
C     &                  M_XR11(MMAX), M_XR12(MMAX), M_XR13(MMAX), 
C     &                  M_VR11(MMAX), M_VR12(MMAX), M_VR13(MMAX), 
CC     &                  M_FC1(MMAX), M_FC2(MMAX), M_FC3(MMAX),
C     &                  M_RS1(MMAX), M_L1(MMAX), M_TE1(MMAX), 
C     &                  M_RS2(MMAX), M_L2(MMAX), M_TE2(MMAX), 
C     &                  M_RS3(MMAX), M_L3(MMAX), M_TE3(MMAX), 
C     &                  M_RC1(MMAX), M_MC1(MMAX), M_RC2(MMAX), 
C     &                  M_MC2(MMAX), M_RC3(MMAX), M_MC3(MMAX),
C     &                  M_A0(MMAX), M_ECC0(MMAX), M_P0(MMAX),
C     &                  M_A1(MMAX), M_ECC1(MMAX), M_P1(MMAX),
C     &                  NM_KW1(MMAX), NM_NAM1(MMAX), NM_KW2(MMAX),
C     &                  NM_NAM2(MMAX), NM_KW3(MMAX), NM_NAM3(MMAX),
C     &                  NM_KWC(MMAX), NM_NAMC(MMAX)
C      REAL*4 M_M1, M_M2, M_M3, M_XC1, M_XC2, M_XC3, M_VC1, M_VC2, M_VC3
C      REAL*4 M_XR01, M_XR02, M_XR03, M_VR01, M_VR02, M_VR03
C      REAL*4 M_XR11, M_XR12, M_XR13, M_VR11, M_VR12, M_VR13
CC      REAL*4 M_FC1, M_FC2, M_FC3
C      REAL*4 M_RS1, M_L1, M_TE1, M_RS2, M_L2, M_TE2, M_RS3, M_L3, M_TE3
C      REAL*4 M_RC1, M_MC1, M_RC2, M_MC2, M_RC3, M_MC3
C      REAL*4 M_A0, M_ECC0, M_P0, M_A1, M_ECC1, M_P1
C      INTEGER NM_KW1, NM_NAM1, NM_KW2, NM_NAM2
C      INTEGER NM_KW3, NM_NAM3, NM_KWC, NM_NAMC

      INTEGER NXTLST(NMAX),IWBINC
      INTEGER N_SINGLE,N_BINARY,N_MERGER, IMLST(4,30),NIMLST
      REAL*8 MMLST(30),MTMP,M3
      SAVE TMASK
      DATA TMASK/NMAX*-999.0D0/
      
      TTOT = TIME+TOFF
      IWBINC = NTOT
      N_SINGLE = 0
      N_BINARY = 0
      N_MERGER = 0

*     Try no scaling output if KZ(46)
      IF (KZ(19).EQ.0.AND.KZ(12).EQ.-1) THEN
         RSCALE_OUT=1.0
         MSCALE_OUT=1.0
         DSCALE_OUT=1.0
         VSCALE_OUT=1.0
         RAU_OUT=1.0
         TSCALE_OUT=1.0
         DAYS_OUT=1.0
         FSCALE_OUT=1.0
         POTSCALE_OUT=1.0
         G_OUT=1.0
      else
         RSCALE_OUT=RBAR
         MSCALE_OUT=ZMBAR
         VSCALE_OUT=VSTAR
         DSCALE_OUT=ZMBAR/RBAR**3
         RAU_OUT=RAU
         TSCALE_OUT=TSTAR
         DAYS_OUT=DAYS
         FSCALE_OUT=1.0
         G_OUT=4.302E-3
*         POTSCALE_OUT=G_OUT*ZMBAR/RSCALE_OUT
         POTSCALE_OUT=1.0
      END IF

*     choose full or active particles output
      NLOOP = NXTLEN
      IF(KZ(46).EQ.3.OR.KZ(46).EQ.4.OR.TTOT.EQ.0.0D0) THEN
         NLOOP = NXTLIMIT
      END IF
      
      DO 1 J = 1, NLOOP
         I = NXTLST(J)
*
         S = TIME - T0(I)
         S1 = 1.5*S
         S2 = 2.0*S
*     Note X may become corrected value X0 for zero interval.
         X(1,I) = ((FDOT(1,I)*S + F(1,I))*S + X0DOT(1,I))*S + X0(1,I)
         X(2,I) = ((FDOT(2,I)*S + F(2,I))*S + X0DOT(2,I))*S + X0(2,I)
         X(3,I) = ((FDOT(3,I)*S + F(3,I))*S + X0DOT(3,I))*S + X0(3,I)
         XDOT(1,I) = (FDOT(1,I)*S1 + F(1,I))*S2 + X0DOT(1,I)
         XDOT(2,I) = (FDOT(2,I)*S1 + F(2,I))*S2 + X0DOT(2,I)
         XDOT(3,I) = (FDOT(3,I)*S1 + F(3,I))*S2 + X0DOT(3,I)
*     
*     --10/18/14 18:04-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         IF (NXTLEN.GT.60.and.TTOT.NE.0.D0) THEN
c$$$            RRMIN = 1.E6
c$$$            I2 = IMINR(I)
c$$$            DO L = 2, LIST(1,I)+1
c$$$               K = LIST(L,I)
c$$$               RRC = 0.D0
c$$$               call jpred(K,time,time)
c$$$               DO KK =1,3
c$$$                  RRC = RRC + (X(KK,K)-X(KK,I))**2
c$$$               END DO
c$$$               IF (RRC.LT.RRMIN) THEN
c$$$                  IIC = K
c$$$                  RRMIN = RRC
c$$$               END IF
c$$$               IF (K.EQ.I2) RRI2 = RRC
c$$$            END DO
c$$$            IF(IIC.NE.I2.and.LIST(1,I).NE.0) THEN
c$$$            print*,'J',J,'I1',I,'I2',I2,'IIC',IIC,'RM',RRMIN,'RI2',RRI2,
c$$$     &              'NNB',LIST(1,I),'TIME',TIME,'RS',RS(I),'X',X(1:3,I)
c$$$            call flush(6)
c$$$            END IF
c$$$            IF (J.GE.NXTLEN) call abort()
c$$$         END IF
*     --10/18/14 18:04-lwang-end----------------------------------------*
*
*     --K.S. binary case--
         IF(I.GT.N.AND.NAME(I).GT.0) THEN
            N_BINARY = N_BINARY + 1
            IK = N_BINARY
            IPAIR = I - N
*     Resolve unperturbed binary
            CALL KSRES_OP(IPAIR,J1,J2,XIREL,VIREL,0)
*     Name
            if (GAMMA(IPAIR).LT.GALIMIT.AND.H(IPAIR).LT.0) THEN
               NB_NAMC(IK) = NAME(I)
               NB_NAM1(IK) = NAME(J1)
               NB_NAM2(IK) = NAME(J2)
*     Dynamic information
               B_M1(IK) = REAL(BODY(J1)*MSCALE_OUT)
               B_M2(IK) = REAL(BODY(J2)*MSCALE_OUT)
               B_XC1(IK) = REAL(X(1,I)*RSCALE_OUT) 
               B_XC2(IK) = REAL(X(2,I)*RSCALE_OUT) 
               B_XC3(IK) = REAL(X(3,I)*RSCALE_OUT) 
               B_VC1(IK) = REAL(XDOT(1,I)*VSCALE_OUT) 
               B_VC2(IK) = REAL(XDOT(2,I)*VSCALE_OUT) 
               B_VC3(IK) = REAL(XDOT(3,I)*VSCALE_OUT) 
               B_XR1(IK) = REAL(XIREL(1)*RAU_OUT) 
               B_XR2(IK) = REAL(XIREL(2)*RAU_OUT) 
               B_XR3(IK) = REAL(XIREL(3)*RAU_OUT) 
               B_VR1(IK) = REAL(VIREL(1)*VSCALE_OUT) 
               B_VR2(IK) = REAL(VIREL(2)*VSCALE_OUT) 
               B_VR3(IK) = REAL(VIREL(3)*VSCALE_OUT) 
               B_FC1(IK) = REAL(F(1,I)*FSCALE_OUT)
               B_FC2(IK) = REAL(F(2,I)*FSCALE_OUT)
               B_FC3(IK) = REAL(F(3,I)*FSCALE_OUT)
               B_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
               B_G(IK)   = REAL(GAMMA(IPAIR))
C            B_FDC1(IK) = REAL(FDOT(1,I)*FDSCALE_OUT)
C            B_FDC2(IK) = REAL(FDOT(2,I)*FDSCALE_OUT)
C            B_FDC3(IK) = REAL(FDOT(3,I)*FDSCALE_OUT)
*     Stellar evolution of two components
               IF (KZ(12).GT.0) THEN
                  call sev_one_star(J1,NB_KW1(IK),B_RS1(IK),B_L1(IK),
     &                 B_TE1(IK),B_MC1(IK),B_RC1(IK),B_M1(IK))
                  call sev_one_star(J2,NB_KW2(IK),B_RS2(IK),B_L2(IK),
     &                 B_TE2(IK),B_MC2(IK),B_RC2(IK),B_M2(IK))
               END IF
               NB_KWC(IK) = KSTAR(I)
*     Binary parameters
               SEMI = -0.5*BODY(I)/H(IPAIR)
               B_ECC(IK) = REAL(SQRT((1.0D0 - R(IPAIR)/SEMI)**2+
     &              TDOT2(IPAIR)**2/(BODY(I)*SEMI)))
               B_P(IK) = REAL(DAYS_OUT*SEMI*SQRT(ABS(SEMI)/BODY(I)))
               B_A(IK) = REAL(SEMI*RAU_OUT)
            ELSE
*     If perturbation is too strong (Gamma>GALIMIT) or the binding energy is positive, the K.S. may not be a binary.
               TMASK(J1) = TTOT
               TMASK(J2) = TTOT
               TMASK(I)  = TTOT
               N_SINGLE = N_SINGLE + 1
               IK = N_SINGLE
               call custom_store_single(J1,IK)
*     correction for the force and potential
               DR3I = BODY(J2)/R(IPAIR)**3
               S_F1(IK) = S_F1(IK) - REAL(XIREL(1)*DR3I*FSCALE_OUT)
               S_F2(IK) = S_F2(IK) - REAL(XIREL(2)*DR3I*FSCALE_OUT)
               S_F3(IK) = S_F3(IK) - REAL(XIREL(3)*DR3I*FSCALE_OUT)
               S_POT(IK) = S_POT(IK) -
     &              REAL(BODY(J2)/R(IPAIR)*POTSCALE_OUT)

               N_SINGLE = N_SINGLE + 1
               IK = N_SINGLE
               call custom_store_single(J2,IK)
*     correction for the force and potential
               DR3I = -BODY(J1)/R(IPAIR)**3
               S_F1(IK) = S_F1(IK) - REAL(XIREL(1)*DR3I*FSCALE_OUT) 
               S_F2(IK) = S_F2(IK) - REAL(XIREL(2)*DR3I*FSCALE_OUT) 
               S_F3(IK) = S_F3(IK) - REAL(XIREL(3)*DR3I*FSCALE_OUT) 
               S_POT(IK) = S_POT(IK) -
     &              REAL(BODY(J1)/R(IPAIR)*POTSCALE_OUT)
               
            END IF
*     CONTINUE
            GO TO 1
         END IF
*     
*     --Merger case--
         IF(I.GT.N.AND.NAME(I).LT.0) THEN
            NIMLST = 0
            NAMEI = NAME(I)
            ICM = I
            IPAIR = ICM - N
            I2 = 2*IPAIR
            I1 = I2 - 1
*     Check possible quartet [[B,S],S] and quintet [[B,S],B], save two components seperately
 4          IF (NAMEI.LT.-2*NZERO) THEN
*     Merger index
               IM = 0
*     Ghost index
               IG = 0
*     Notice, here assume merger number is small
               DO K = 1, NMERGE
                  IF (NAMEM(K).EQ.NAMEI) THEN
                     IM = K
                     GO TO 10
                  END IF
               END DO
*     error messsage
               print*,'Error!: Cannot find merger index for name =',
     &              NAMEI, ', NMERGER =',NMERGE
               call flush(6)
               call abort()
*
 10            DO K = NXTLIMIT+1,NXTLIMIT+NGHOSTS
                  IF (NAME(NXTLST(K)).EQ.NAMEG(IM)) THEN
                     IG = NXTLST(K)
                     GO TO 20
                  END IF
               END DO
*     error messsage
               print*,'Error!: Cannot find ghost of merger IM=',
     &              IM, ', NAMEG(IM) =',NAMEG(IM)
               call flush(6)
               call abort()
*
 20            CALL KSRES_OP(IPAIR,J1,J2,XIREL,VIREL,1)
*     Outer component
               IF (NAME(IG).LT.0) THEN
                  NIMLST = NIMLST + 1
                  IMLST(1,NIMLST) = J2
                  IMLST(2,NIMLST) = IG
                  IMLST(3,NIMLST) = NAME(IG)
                  MMLST(NIMLST) = CM(4,IM)
               ELSE IF(NAME(IG).GT.NZERO) THEN
                  N_BINARY = N_BINARY + 1
                  JPAIR = IG - N
                  CALL KSRES_OP(JPAIR,JP1,JP2,XIREL,VIREL,0)
                  call custom_store_binary(N_BINARY,J2,NAME(IG),
     &                 JP1,JP2,CM(3,IM),CM(4,IM),H(JPAIR),R(JPAIR),
     &                 TDOT2(JPAIR),XIREL,VIREL,KSTAR(IG),GAMMA(JPAIR))
               ELSE
                  N_SINGLE = N_SINGLE + 1
                  IK = N_SINGLE
*     Name
                  NS_NAM(IK) = NAME(IG)
*     Dynamics
                  S_M(IK)  = REAL(BODY(J2)*MSCALE_OUT)
                  S_X1(IK) = REAL(X(1,J2)*RSCALE_OUT) 
                  S_X2(IK) = REAL(X(2,J2)*RSCALE_OUT) 
                  S_X3(IK) = REAL(X(3,J2)*RSCALE_OUT) 
                  S_V1(IK) = REAL(XDOT(1,J2)*VSCALE_OUT) 
                  S_V2(IK) = REAL(XDOT(2,J2)*VSCALE_OUT) 
                  S_V3(IK) = REAL(XDOT(3,J2)*VSCALE_OUT) 
                  S_F1(IK) = REAL(F(1,I)*FSCALE_OUT) 
                  S_F2(IK) = REAL(F(2,I)*FSCALE_OUT) 
                  S_F3(IK) = REAL(F(3,I)*FSCALE_OUT) 
                  S_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
C                  S_FD1(IK) = REAL(FDOT(1,I)*FDSCALE_OUT) 
C                  S_FD2(IK) = REAL(FDOT(2,I)*FDSCALE_OUT) 
C                  S_FD3(IK) = REAL(FDOT(3,I)*FDSCALE_OUT) 
*     Stellar evolution of single star
                  IF (KZ(12).GT.0) THEN
                     MTMP = BODY(IG)
                     BODY(IG) = BODY(J2)
                     call sev_one_star(IG,NS_KW(IK),S_RS(IK),S_L(IK),
     &                    S_TE(IK),S_MC(IK),S_RC(IK),S_M(IK))
                     BODY(IG) = MTMP
                  END IF
               END IF
*     Inner merger
               NIMLST = NIMLST + 1
               IMLST(1,NIMLST) = J1
               IMLST(2,NIMLST) = -IM
               IMLST(3,NIMLST) = NAMEI + 2*NZERO
               MMLST(NIMLST) = CM(2,IM)
               GO TO 8
            END IF
*     For normal merger
            N_MERGER = N_MERGER + 1
            IK = N_MERGER
            DO K = 1, NMERGE
               IF (NAMEM(K).EQ.NAMEI) THEN
                  IM = K
                  GO TO 6
               END IF
            END DO
*     error messsage
            print*,'Error!: Cannot find merger index for name =',
     &           NAMEI, ', NMERGER =',NMERGE
            call flush(6)
            call abort()
*     Check B-B system, save as two binaries
 6          IF (NAME(I2).GT.NZERO) THEN
               DO K = NXTLIMIT+1,NXTLIMIT+NGHOSTS
                  IF (NAME(NXTLST(K)).EQ.NAMEG(IM)) THEN
                     IG = NXTLST(K)
                     GO TO 17
                  END IF
               END DO
*     error messsage
               print*,'Error!: Cannot find ghost of merger IM=',
     &              IM, ', NAMEG(IM) =',NAMEG(IM)
               call flush(6)
               call abort()
 17            CALL KSRES_OP(IPAIR,J1,J2,XIREL,VIREL,1)
*     Inner binary
               N_BINARY = N_BINARY + 1
               RIP = SQRT(XREL(1,IM)**2+XREL(2,IM)**2+XREL(3,IM)**2)
               RDOT = XREL(1,IM)*VREL(1,IM) + XREL(2,IM)*VREL(2,IM) + 
     &              XREL(3,IM)*VREL(3,IM)
               call custom_store_binary(N_BINARY,J1,NAMEM(IM),J1,IG,
     &              CM(1,IM),CM(2,IM),HM(IM),RIP,RDOT,
     &              XREL(1,IM),VREL(1,IM),KSTARM(IM),GALIMIT)
*     Outer binary
               N_BINARY = N_BINARY + 1
               JPAIR = IG - N
               call ksres_op(JPAIR,JP1,JP2,XIREL,VIREL,0)
               call custom_store_binary(N_BINARY,J2,NAME(I2),JP1,JP2,
     &              CM(3,IM),CM(4,IM),H(JPAIR),R(JPAIR),TDOT2(JPAIR),
     &              XIREL,VIREL,KSTAR(I2),GALIMIT)
               GO TO 8
            END IF
*     Resolve merger
            M_XC1(IK) = REAL(X(1,ICM)*RSCALE_OUT) 
            M_XC2(IK) = REAL(X(2,ICM)*RSCALE_OUT) 
            M_XC3(IK) = REAL(X(3,ICM)*RSCALE_OUT) 
            M_VC1(IK) = REAL(XDOT(1,ICM)*VSCALE_OUT) 
            M_VC2(IK) = REAL(XDOT(2,ICM)*VSCALE_OUT) 
            M_VC3(IK) = REAL(XDOT(3,ICM)*VSCALE_OUT)
            M_FC1(IK) = REAL(F(1,ICM)*FSCALE_OUT)
            M_FC2(IK) = REAL(F(2,ICM)*FSCALE_OUT)
            M_FC3(IK) = REAL(F(3,ICM)*FSCALE_OUT)
            M_POT(IK) = REAL(-PHIDBL(ICM)*POTSCALE_OUT)

*     Check whether ipair is avaiable
            IF(IPAIR.GT.0) THEN
               CALL KSRES_OP(IPAIR,J1,J2,XIREL,VIREL,0)
               M_XR11(IK) = REAL(XIREL(1)*RAU_OUT)
               M_XR12(IK) = REAL(XIREL(2)*RAU_OUT)
               M_XR13(IK) = REAL(XIREL(3)*RAU_OUT)
               M_VR11(IK) = REAL(VIREL(1)*VSCALE_OUT)
               M_VR12(IK) = REAL(VIREL(2)*VSCALE_OUT)
               M_VR13(IK) = REAL(VIREL(3)*VSCALE_OUT)
            ELSE
               IM2 = -IPAIR
               M_XR11(IK) = REAL(XREL(1,IM2)*RAU_OUT) 
               M_XR12(IK) = REAL(XREL(2,IM2)*RAU_OUT) 
               M_XR13(IK) = REAL(XREL(3,IM2)*RAU_OUT) 
               M_VR11(IK) = REAL(VREL(1,IM2)*VSCALE_OUT) 
               M_VR12(IK) = REAL(VREL(2,IM2)*VSCALE_OUT) 
               M_VR13(IK) = REAL(VREL(3,IM2)*VSCALE_OUT) 
            END IF
            M_XR01(IK) = REAL(XREL(1,IM)*RAU_OUT) 
            M_XR02(IK) = REAL(XREL(2,IM)*RAU_OUT) 
            M_XR03(IK) = REAL(XREL(3,IM)*RAU_OUT) 
            M_VR01(IK) = REAL(VREL(1,IM)*VSCALE_OUT) 
            M_VR02(IK) = REAL(VREL(2,IM)*VSCALE_OUT) 
            M_VR03(IK) = REAL(VREL(3,IM)*VSCALE_OUT) 
            M_M1(IK) = REAL(CM(1,IM)*MSCALE_OUT) 
            M_M2(IK) = REAL(CM(2,IM)*MSCALE_OUT) 
            M_M3(IK) = REAL(BODY(J2)*MSCALE_OUT)
            
            NM_NAMC(IK) = NAMEI
            NM_NAM1(IK) = NAME(J1)
            NM_NAM2(IK) = NAMEG(IM)
            NM_NAM3(IK) = NAME(J2)

            SEMI = 0.D0
            IF (IPAIR.GT.0) THEN
               SEMI = -0.5*BODY(ICM)/H(IPAIR)
               M_ECC1(IK) = REAL(SQRT((1.0D0 - R(IPAIR)/SEMI)**2 + 
     &              TDOT2(IPAIR)**2/(BODY(ICM)*SEMI)))
            ELSE
               IM2 = -IPAIR
               SEMI = -0.5*BODY(ICM)/HM(IM2)
               RJ = SQRT(XREL(1,IM2)**2+XREL(2,IM2)**2+XREL(3,IM2)**2)
C               RDOT = XREL(1,IM2)*VREL(1,IM2) + XREL(2,IM2)*VREL(2,IM2)
C     &              + XREL(3,IM2)*VREL(3,IM2)
               RDOT = 2.0*(UM(1,IM2)*UMDOT(1,IM2)+UM(2,IM2)*UMDOT(2,IM2)
     &              +UM(3,IM2)*UMDOT(3,IM2)+UM(4,IM2)*UMDOT(4,IM2))
               M_ECC1(IK) = REAL(SQRT((1.0D0 - RJ/SEMI)**2 + 
     &              RDOT**2/(BODY(ICM)*SEMI)))
            END IF
            M_P1(IK) = REAL(DAYS_OUT*SEMI*SQRT(ABS(SEMI)/BODY(ICM)))
            M_A1(IK) = REAL(SEMI*RAU_OUT)
            SEMI0 = -0.5*BODY(J1)/HM(IM)
            RJ = SQRT(XREL(1,IM)**2 + XREL(2,IM)**2 + XREL(3,IM)**2)
C            RDOT = XREL(1,IM)*VREL(1,IM) + XREL(2,IM)*VREL(2,IM)
C     &           + XREL(3,IM)*VREL(3,IM)
            RDOT = 2.0*(UM(1,IM)*UMDOT(1,IM) + UM(2,IM)*UMDOT(2,IM) +
     &           UM(3,IM)*UMDOT(3,IM) + UM(4,IM)*UMDOT(4,IM))
            M_ECC0(IK) = REAL(SQRT((1.0D0 - RJ/SEMI0)**2 +
     &           RDOT**2/(BODY(J1)*SEMI0)))
            M_P0(IK) = REAL(DAYS_OUT*SEMI0*SQRT(ABS(SEMI0)/BODY(J1)))
            M_A0(IK) = REAL(SEMI0*RAU_OUT)
            
*     stellar evolution
            IF (KZ(12).GT.0) THEN
               DO K = NXTLIMIT+1,NXTLIMIT+NGHOSTS
                  IF (NAME(NXTLST(K)).EQ.NAMEG(IM)) THEN
                     IG = NXTLST(K)
                     GO TO 7
                  END IF
               END DO
*     error messsage
               print*,'Error!: Cannot find ghost of merger IM=',
     &              IM, ', NAMEG(IM) =',NAMEG(IM)
               call flush(6)
               call abort()
*     first star
 7             MTMP = BODY(J1)
               BODY(J1) = CM(1,IM)
               call sev_one_star(J1,NM_KW1(IK),M_RS1(IK),M_L1(IK),
     &              M_TE1(IK),M_MC1(IK),M_RC1(IK),M_M1(IK))
               BODY(J1) = MTMP
*     second star
               MTMP = BODY(IG)
               BODY(IG) = CM(2,IM)
               call sev_one_star(IG,NM_KW2(IK),M_RS2(IK),M_L2(IK),
     &              M_TE2(IK),M_MC2(IK),M_RC2(IK),M_M2(IK))
               BODY(IG) = MTMP
*     third star
               call sev_one_star(J2,NM_KW3(IK),M_RS3(IK),M_L3(IK),
     &              M_TE3(IK),M_MC3(IK),M_RC3(IK),M_M3(IK))
               NM_KWC(IK) = KSTAR(ICM)
            END IF
 8          IF(NIMLST.GT.0) THEN
*     ICM (M,X,XD), IPAIR(XIREL,VIREL), NAMEI, M3
               ICM = IMLST(1,NIMLST)
               IPAIR = IMLST(2,NIMLST)
               NAMEI = IMLST(3,NIMLST)
               M3 = MMLST(NIMLST)
               NIMLST = NIMLST - 1
               GO TO 4
            END IF
            GO TO 1
         END IF
*
*     --Chain case--
         IF(NAME(I).EQ.0) THEN
            DO L = 1,NCH
               N_SINGLE = N_SINGLE + 1
               IK = N_SINGLE
               NS_NAM(IK) = NAMEC(L)
               S_M(IK)  = REAL(BODYC(L)*MSCALE_OUT)
               S_X1(IK) = REAL((XCH(3*L-2)+X(1,ICH))*RSCALE_OUT)
               S_X2(IK) = REAL((XCH(3*L-1)+X(2,ICH))*RSCALE_OUT)
               S_X3(IK) = REAL((XCH(3*L)+X(3,ICH))*RSCALE_OUT)
               S_V1(IK) = REAL((VCH(3*L-2)+XDOT(1,ICH))*VSCALE_OUT)
               S_V2(IK) = REAL((VCH(3*L-1)+XDOT(2,ICH))*VSCALE_OUT)
               S_V3(IK) = REAL((VCH(3*L)+XDOT(3,ICH))*VSCALE_OUT)
               S_F1(IK) = REAL(F(1,I)*FSCALE_OUT) 
               S_F2(IK) = REAL(F(2,I)*FSCALE_OUT) 
               S_F3(IK) = REAL(F(3,I)*FSCALE_OUT) 
               S_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
*     stellar evolution         
               IF(KZ(12).GT.0) THEN
                  IF (L.EQ.1) THEN
                     MTMP = BODY(I)
                     BODY(I) = BODYC(L)
                     call sev_one_star(I,NS_KW(IK),S_RS(IK),S_L(IK),
     &                    S_TE(IK),S_MC(IK),S_RC(IK),S_M(IK))
                     BODY(I) = MTMP
                  ELSE
                     IG = 0
                     DO K = NXTLIMIT+1,NXTLIMIT+NGHOSTS
                        IF(NAME(NXTLST(K)).EQ.NAMEC(L)) THEN
                           IG = NXTLST(K)
                           GO TO 11
                        END IF
                     END DO
*     error message
                     print*,'Error!: Chain ghost name =',NAMEC(L),
     &                    'not found.'
                     call flush(6)
                     call abort()
 11                  MTMP = BODY(IG)
                     BODY(IG) = BODYC(L)
                     call sev_one_star(IG,NS_KW(IK),S_RS(IK),S_L(IK),
     &                    S_TE(IK),S_MC(IK),S_RC(IK),S_M(IK))
                     BODY(IG) = MTMP
                  END IF
               END IF
            END DO
*     continue
            GO TO 1
         END IF
*
*     Nearest neighbor index
         I2 = IMINR(I)
         IF(I2.LE.0) I2 = I + 1
 12      IF(BODY(I2).EQ.0.D0.AND.I2.LE.N) THEN
            I2 = I2 + 1
            GO TO 12
         END IF
            
*     --Single star case--
*     escape the registered wide binary index or second component is multiple system
         IF(TMASK(I).EQ.TTOT.OR.TMASK(I2).EQ.TTOT
     &        .OR.I2.GT.N.OR.NAME(I2).LE.0) THEN
            IF(TMASK(I).NE.TTOT) THEN
               TMASK(I) = TTOT
               N_SINGLE = N_SINGLE + 1
               call custom_store_single(I,N_SINGLE)
            END IF
*     CONTINUE
            GO TO 1
*
         ELSE
*     --Wide binary case--
            call JPRED_int(I2,TIME)
            RIJ2 = 0.D0
            VIREL2 = 0.D0
            RDOT = 0.D0
            DO K = 1,3
               XIREL(K) = X(K,I) - X(K,I2)
               VIREL(K) = XDOT(K,I) - XDOT(K,I2)
               RIJ2 = RIJ2 + XIREL(K)**2
               VIREL2 = VIREL2 + VIREL(K)**2
               RDOT = RDOT + XIREL(K)*VIREL(K)
            END DO
            RIJ = SQRT(RIJ2)
            ZMB = BODY(I) + BODY(I2)
            SEMI = 2.0/RIJ - VIREL2/ZMB
            SEMI = 1.0/SEMI
*     Check whether satisfy wide binary criterion
            IF (SEMI.GT.0.0.AND.SEMI.LT.0.1*RSCALE) THEN
               IWBINC = IWBINC + 1
               N_BINARY = N_BINARY + 1
               IK = N_BINARY
*     set time flag to avoid duplicate case
               TMASK(I) = TTOT
               TMASK(I2) = TTOT
*     Name
               NB_NAMC(IK) = IWBINC
               NB_NAM1(IK) = NAME(I)
               NB_NAM2(IK) = NAME(I2)
*     Dynamic information
               B_M1(IK) = REAL(BODY(I)*MSCALE_OUT)
               B_M2(IK) = REAL(BODY(I2)*MSCALE_OUT)

               B_XC1(IK) = REAL((BODY(I)*X(1,I)+BODY(I2)*X(1,I2))
     &              /ZMB*RSCALE_OUT)
               B_XC2(IK) = REAL((BODY(I)*X(2,I)+BODY(I2)*X(2,I2))
     &              /ZMB*RSCALE_OUT)
               B_XC3(IK) = REAL((BODY(I)*X(3,I)+BODY(I2)*X(3,I2))
     &              /ZMB*RSCALE_OUT)
               B_VC1(IK) = REAL((BODY(I)*XDOT(1,I)+BODY(I2)*XDOT(1,I2))
     &              /ZMB*VSCALE_OUT)
               B_VC2(IK) = REAL((BODY(I)*XDOT(2,I)+BODY(I2)*XDOT(3,I2))
     &              /ZMB*VSCALE_OUT)
               B_VC3(IK) = REAL((BODY(I)*XDOT(3,I)+BODY(I2)*XDOT(3,I2))
     &              /ZMB*VSCALE_OUT)
               B_XR1(IK) = REAL(XIREL(1)*RAU_OUT) 
               B_XR2(IK) = REAL(XIREL(2)*RAU_OUT) 
               B_XR3(IK) = REAL(XIREL(3)*RAU_OUT) 
               B_VR1(IK) = REAL(VIREL(1)*VSCALE_OUT) 
               B_VR2(IK) = REAL(VIREL(2)*VSCALE_OUT) 
               B_VR3(IK) = REAL(VIREL(3)*VSCALE_OUT) 
               DR3I = BODY(I2)/(RIJ2*RIJ)
*     Only estimate the c.m. force at the first component position
               B_FC1(IK) = REAL((F(1,I) + XIREL(1)*DR3I)*FSCALE_OUT) 
               B_FC2(IK) = REAL((F(2,I) + XIREL(2)*DR3I)*FSCALE_OUT) 
               B_FC3(IK) = REAL((F(3,I) + XIREL(3)*DR3I)*FSCALE_OUT) 
               B_POT(IK) = REAL((-PHIDBL(I) +BODY(I2)/RIJ)*POTSCALE_OUT)
               B_G(IK) = REAL(2*SQRT((F(1,I)-F(1,I2))**2
     &              +(F(2,I)-F(2,I2))**2
     &              +(F(3,I)-F(3,I2))**2)*RIJ2/(BODY(I)+BODY(I2)))
*     Stellar evolution of two components
               IF (KZ(12).GT.0) THEN
                  call sev_one_star(I,NB_KW1(IK),B_RS1(IK),B_L1(IK),
     &                 B_TE1(IK),B_MC1(IK),B_RC1(IK),B_M1(IK))
                  call sev_one_star(I2,NB_KW2(IK),B_RS2(IK),B_L2(IK),
     &                 B_TE2(IK),B_MC2(IK),B_RC2(IK),B_M2(IK))
               END IF
               NB_KWC(IK) = -1
*     Binary parameters
               B_ECC(IK) = REAL(SQRT((1.0D0 - RIJ/SEMI)**2 + 
     &              RDOT**2/(ZMB*SEMI)))
               BODYI = (B_M1(IK) + B_M2(IK))/MSCALE_OUT
               B_P(IK) = REAL(DAYS_OUT*SEMI*SQRT(ABS(SEMI)/BODYI))
               B_A(IK) = REAL(SEMI*RAU_OUT)
*     Continue
               GO TO 1
            ELSE
*     Treat as single star
               TMASK(I) = TTOT
               N_SINGLE = N_SINGLE + 1
               call custom_store_single(I,N_SINGLE)
*     CONTINUE
               GO TO 1
            END IF
         END IF
*
 1    CONTINUE
      
*     output
*     --10/25/14 22:16-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      if(rank.eq.0) print*,'CUSTOM OUTPUT: TTOT',TTOT,'N_SINGLE',
     &     N_SINGLE,'N_BINARY',N_BINARY,'N_MERGER',N_MERGER
c$$$      print*,'NS_NAM',NS_NAM(1:5),'S_M',S_M(1:5),'S_X1',S_X1(1:5),
c$$$     &     'S_V1',S_V1(1:5),'NS_KW',NS_KW(1:5),'S_RS',S_RS(1:5),
c$$$     &     'S_L',S_L(1:5),'S_TE',S_TE(1:5),'S_MC',S_MC(1:5),'S_RC',
c$$$     &     S_RC(1:5)
*     --10/25/14 22:16-lwang-end----------------------------------------*

      if(rank.eq.0) call output_single(TTOT,N_SINGLE,KZ(12),KZ(46))
c$$$     &     NS_NAM,
c$$$     &     S_M,S_X1,S_X2,S_X2,S_V1,S_V2,S_V3,
c$$$     &     S_F1,S_F2,S_F3,
c$$$     &     NS_KW,S_RS,S_L,S_TE,S_MC,S_RC)

      if(rank.eq.0) call output_binary(TTOT,N_BINARY,KZ(12),KZ(46))
c$$$     &     NS_NAMC,NS_NAM1,NS_NAM2,
c$$$     &     B_M1,B_M2,B_XC1,B_XC2,B_XC3,B_VC1,B_VC2,B_VC3,
c$$$     &     B_XR1,B_XR2,B_XR3,B_VR1,B_VR2,B_VR3,
c$$$C     &     B_FC1,B_FC2,B_FC3,
c$$$     &     B_A,B_ECC,B_P,
c$$$     &     NB_KWC,NB_KW1,NB_KW2,B_RS1,B_RS2,B_L1,B_L2,
c$$$     &     B_TE1,B_TE2,B_MC1,B_MC2,B_RC1,B_RC2)

      if(rank.eq.0) call output_merger(TTOT,N_MERGER,KZ(12),KZ(46))
c$$$     &     NM_NAMC,NM_NAM1,NM_NAM2,NM_NAM3,
c$$$     &     M_M1,M_M2,M_M3,M_XC1,M_XC2,M_XC3,M_VC1,M_VC2.M_VC3,
c$$$     &     M_XR01,M_XR02,M_XR03,M_VR01,M_VR02,M_VR03,
c$$$     &     M_XR11,M_XR12,M_XR13,M_VR11,M_VR12,M_VR13,
c$$$C     &     M_FC1,M_FC2,M_FC3,
c$$$     &     M_A0,M_A2,M_ECC0,M_ECC1,M_P0,M_P1,
c$$$     &     NM_KWC,NM_KW1,NM_KW2,NM_KW3,
c$$$     &     M_RS1,M_RS2,M_RS3,M_L1,M_L2,M_L3,M_TE1,M_TE2,M_TE3,
c$$$     &     M_MC1,M_MC2,M_MC3,M_RC1,M_RC2,M_RC3)

      if (kz(14) .eq. 3) then
         SX = RBAR/1000.0
         GZ = RG(1)*VG(2) - RG(2)*VG(1)
         write (60,*) TTOT, (RG(K)*SX,K=1,3), (VG(K)*VSTAR,K=1,3),
     &        GZ, ETIDE
         call flush(60)
      end if

      if (NXTLEN.eq.NXTLIMIT.or.TTOT.eq.0.D0) then
         ZSM0 = ZSMASS0*MSCALE_OUT
         ZBM0 = ZBMASS0*MSCALE_OUT
         call global_output(ZSM0,ZBM0,RC,
     &        N_SINGLE,N_BINARY,N_MERGER,KZ(7),KZ(19))
      END IF
      
      RETURN

      END
      
*----------------------------------------------------------------------------
      subroutine sev_one_star(I,KW,RM,LUM,TE,MC,RCC,M1)
*
*
*     Get stellar evolution parameters for one star
*     ---------------------------------------------
*
      include 'common6.h'
      include 'scale_out.h'
      INTEGER I,KW
      REAL*8  LUMS(10),TSCLS(20),GB(10)
      REAL*8  RM8,LUM8,MC8,RCC8,M18,M0
      REAL*4  RM,LUM,TE,MC,RCC,M1

*
      KW = KSTAR(I)
      AGE = MAX(TIME,TEV0(I))*TSCALE_OUT - EPOCH(I)
      M0 = BODY0(I)*MSCALE_OUT
      M18 = BODY(I)*MSCALE_OUT
      CALL STAR(KW,M0,M18,TM,TN,TSCLS,LUMS,GB,ZPARS)
      CALL HRDIAG(M0,AGE,M18,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &     RM8,LUM8,KW,MC8,RCC8,ME,RE,K2)
*     Temperature (Solar temperature got from Williams, D. R. (1 July 2013). "Sun Fact Sheet". NASA. Retrieved 12 August 2013.)
      TE = REAL(5778*(LUM8/(RM8*RM8))**0.25)

      RM = REAL(RM8)
      LUM = REAL(LUM8)
      MC = REAL(MC8)
      RCC = REAL(RCC8)
      M1 = REAL(M18)

      RETURN

      END

*----------------------------------------------------------------------------
      subroutine custom_store_single(I,IK)
*
*
*     Save one single star in vectors
*     ---------------------------------
*     
      include 'common6.h'
      include 'output_single.h'
      include 'scale_out.h'
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
C*     single
C      COMMON/OUTSINGLE/ S_M(NMAX),S_X1(NMAX), S_X2(NMAX), S_X3(NMAX),
C     &                  S_V1(NMAX), S_V2(NMAX), S_V3(NMAX),
CC     &                  S_F1(NMAX), S_F2(NMAX), S_F3(NMAX),
C     &                  S_RS(NMAX), S_L(NMAX), S_TE(NMAX),
C     &                  S_RC(NMAX), S_MC(NMAX),
C     &                  NS_KW(NMAX),NS_NAM(NMAX)
C      REAL*4 S_M, S_X1, S_X2, S_X3, S_V1, S_V2, S_V3
CC      REAL*4 S_F1, S_F2, S_F3
C      REAL*4 S_RS, S_L, S_TE, S_RC, S_MC
C      INTEGER NS_KW,NS_NAM

*     Name
      NS_NAM(IK) = NAME(I)
*     Dynamics
      S_M(IK)  = REAL(BODY(I)*MSCALE_OUT)
      S_X1(IK) = REAL(X(1,I)*RSCALE_OUT) 
      S_X2(IK) = REAL(X(2,I)*RSCALE_OUT) 
      S_X3(IK) = REAL(X(3,I)*RSCALE_OUT) 
      S_V1(IK) = REAL(XDOT(1,I)*VSCALE_OUT) 
      S_V2(IK) = REAL(XDOT(2,I)*VSCALE_OUT) 
      S_V3(IK) = REAL(XDOT(3,I)*VSCALE_OUT) 
      S_F1(IK) = REAL(F(1,I)*FSCALE_OUT) 
      S_F2(IK) = REAL(F(2,I)*FSCALE_OUT) 
      S_F3(IK) = REAL(F(3,I)*FSCALE_OUT) 
      S_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
*     Stellar evolution of single star
      IF (KZ(12).GT.0) THEN
         call sev_one_star(I,NS_KW(IK),S_RS(IK),S_L(IK),
     &        S_TE(IK),S_MC(IK),S_RC(IK),S_M(IK))
      END IF

      RETURN

      END

*----------------------------------------------------------------------------
      subroutine custom_store_binary(IK,I,NCM,J1,J2,M1,M2,HIP,RIP,RDOT,
     &     XIREL,VIREL,KC,GI)
*
*
*     Save special binary system
*     ---------------------------------
*     
      include 'common6.h'
      include 'output_binary.h'
      include 'scale_out.h'
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)

C*     Binary
C      COMMON/OUTBINARY/ B_M1(KMAX), B_M2(KMAX),
C     &                  B_XC1(KMAX), B_XC2(KMAX), B_XC3(KMAX),
C     &                  B_VC1(KMAX), B_VC2(KMAX), B_VC3(KMAX), 
C     &                  B_XR1(KMAX), B_XR2(KMAX), B_XR3(KMAX), 
C     &                  B_VR1(KMAX), B_VR2(KMAX), B_VR3(KMAX), 
CC     &                  B_FC1(KMAX), B_FC2(KMAX), B_FC3(KMAX), 
C     &                  B_RS1(KMAX), B_L1(KMAX), B_TE1(KMAX),
C     &                  B_RS2(KMAX), B_L2(KMAX), B_TE2(KMAX),
C     &                  B_RC1(KMAX), B_MC1(KMAX), B_RC2(KMAX), 
C     &                  B_MC2(KMAX), B_A(KMAX), B_ECC(KMAX), B_P(KMAX),
C     &                  NB_KW1(KMAX), NB_NAM1(KMAX), NB_KW2(KMAX), 
C     &                  NB_NAM2(KMAX), NB_KWC(KMAX), NB_NAMC(KMAX)
C      REAL*4 B_M1, B_M2, B_XC1, B_XC2, B_XC3, B_VC1, B_VC2, B_VC3
C      REAL*4 B_XR1, B_XR2, B_XR3, B_VR1, B_VR2, B_VR3
CC      REAL*4 B_FC1, B_FC2, B_FC3
C      REAL*4 B_RS1, B_L1, B_TE1, B_RS2, B_L2, B_TE2
C      REAL*4 B_RC1, B_MC1, B_RC2, B_MC2,B_A, B_ECC, B_P
C      INTEGER NB_KW1, NB_NAM1, NB_KW2, NB_NAM2, NB_KWC, NB_NAMC

      REAL*8 XIREL(3), VIREL(3), M1, M2, MTMP,MCM,RDOT,GI

*     Name
      NB_NAMC(IK) = NCM
      NB_NAM1(IK) = NAME(J1)
      NB_NAM2(IK) = NAME(J2)
*     Dynamic information
      B_M1(IK) = REAL(M1*MSCALE_OUT) 
      B_M2(IK) = REAL(M2*MSCALE_OUT) 
      B_XC1(IK) = REAL(X(1,I)*RSCALE_OUT) 
      B_XC2(IK) = REAL(X(2,I)*RSCALE_OUT) 
      B_XC3(IK) = REAL(X(3,I)*RSCALE_OUT) 
      B_VC1(IK) = REAL(XDOT(1,I)*VSCALE_OUT) 
      B_VC2(IK) = REAL(XDOT(2,I)*VSCALE_OUT) 
      B_VC3(IK) = REAL(XDOT(3,I)*VSCALE_OUT) 
      B_XR1(IK) = REAL(XIREL(1)*RAU_OUT) 
      B_XR2(IK) = REAL(XIREL(2)*RAU_OUT) 
      B_XR3(IK) = REAL(XIREL(3)*RAU_OUT) 
      B_VR1(IK) = REAL(VIREL(1)*VSCALE_OUT) 
      B_VR2(IK) = REAL(VIREL(2)*VSCALE_OUT) 
      B_VR3(IK) = REAL(VIREL(3)*VSCALE_OUT) 
      B_FC1(IK) = REAL(F(1,I)*FSCALE_OUT)
      B_FC2(IK) = REAL(F(2,I)*FSCALE_OUT)
      B_FC3(IK) = REAL(F(3,I)*FSCALE_OUT)
      B_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
      B_G(IK)   = REAL(GI)
*     Stellar evolution of two components
      IF (KZ(12).GT.0) THEN
         MTMP = BODY(J1)
         BODY(J1) = M1
         call sev_one_star(J1,NB_KW1(IK),B_RS1(IK),B_L1(IK),
     &        B_TE1(IK),B_MC1(IK),B_RC1(IK),B_M1(IK))
         BODY(J1) = MTMP
         MTMP = BODY(J2)
         BODY(J2) = M2
         call sev_one_star(J2,NB_KW2(IK),B_RS2(IK),B_L2(IK),
     &        B_TE2(IK),B_MC2(IK),B_RC2(IK),B_M2(IK))
         BODY(J2) = MTMP
      END IF
      NB_KWC(IK) = KC
*     Binary parameters
      MCM = M1 + M2
      SEMI = -0.5*MCM/HIP
      B_ECC(IK) = REAL(SQRT((1.0D0 - RIP/SEMI)**2 + 
     &     RDOT**2/(MCM*SEMI)))
      B_P(IK) = REAL(DAYS_OUT*SEMI*SQRT(ABS(SEMI)/MCM))
      B_A(IK) = REAL(SEMI*RAU_OUT)

      RETURN

      END
