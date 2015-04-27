      SUBROUTINE KSREL(IPAIR,I1,M1,M2,XI,VI,DT,FRELP,FRELD,VWHOLE,APPLY,
     &     KSAK,KSBK)
*     
*     Relativistic perturbation on KS pair.
*     -------------------------------------
*     
*     Coded by Pau Amaro-Seoane 2006.
*     -------------------------------
*     
*     3PN,3.5PN and Spin terms added by Patrick Brem 2011.
*     ----------------------------------------------------
*     
C      INCLUDE 'common6.h'
      include 'params.h'
      include 'postnewton.h'
C     COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
C     &                LISTC(LMAX)
      
      REAL*8 XI(6),VI(6),FRELP(3),FRELD(3),DT
      REAL*8 RIJ, RIJ2,A(3),AD(3),KSAK,KSBK,ADK,BDK,PI2
      REAL*8 NX,NY,NZ,RCSCH
      REAL*8 KSV(3),KSX(3)
      REAL*8 CL, CL3, CL2, CL4, CL5, CL6, CL7
      REAL*8 A1,A2,A2_5,A3,A3_5,B1,B2,B2_5,B3,B3_5,ETA,MASAS,RP,MOR,RPP
      REAL*8 A1D,A2D,A2_5D,A3D,A3_5D,B1D,B2D,B2_5D,B3D,B3_5D,VA
C     REAL*8 KSECC,KSH,V1_V22,V12,V22,VWHOLE
      REAL*8 V1_V22,VWHOLE
C     REAL*8 XCOMPIMP,YCOMPIMP,ZCOMPIMP,FRP,GREL
*     (c) patrick
      REAL*8 KSS(3),KSSIG(3),NVEC(3),AT(3)
      REAL*8 NCV(3),NCS(3),NCSIG(3),VCS(3),VCSIG(3)
      REAL*8 SDNCV,SIGDNCV,NDV,DM,S1(3),S2(3)
      REAL*8 C1_5(3),C2(3),C2_5(3),C1_5D(3),C2_5D(3),C2D(3)
      REAL*8 VDS,VDSIG,NDS,NDSIG,S1DLU,S2DLU
      REAL*8 SU1(3),SU2(3),SU3(3),SV(3),SV1(3),SV2(3)
      REAL*8 L(3),LU(3) 
      REAL*8 LABS,SS1(3),SS2(3)
      REAL*8 NVDOT, NDOT(3), SNVDOT, SIGNVDOT
      REAL*8 NDOTCV(3),NCA(3),NDOTCS(3),NCSU(3)
      REAL*8 NDOTCSIG(3),NCSV(3),ACS(3),VCSU(3)
      REAL*8 ACSIG(3),VCSV(3),NSDOT,NSIGDOT,VSDOT,VSIGDOT
      REAL*8 XS(3),XA(3),NXA,NXS,XS2,XA2,XAD(3),XSD(3)
      REAL*8 NXADOT, NXSDOT
C     REAL*8 NXADOT, NXSDOT, ESPIN, ANEWTON, SUM,E_CORRECT
*     For energy correction
C      REAL*8 SQM(3), CHI12,CHI22
*     variables such as NVEC and RIJ are recycled in the perturber sum at 101
      REAL*8 M1,M2

      INTEGER I1,I2,APPLY,IPAIR
      INTEGER SPNI1,SPNI2

      PI2 = 9.869604401089359
      I2 = I1 + 1
*     refer to the individual spins by the name of the individual particles rather than I1 and I2
      SPNI1 = I1
      SPNI2 = I2
*      IF (SPNI1.GT.N) THEN
*        SPN(1,SPNI1) = 0.0
*        SPN(2,SPNI1) = 0.0
*        SPN(3,SPNI1) = 0.0
*      END IF !Because treating the spin of a KS pair makes no sense
*      IF (SPNI2.GT.N) THEN
*        SPN(1,SPNI2) = 0.0
*        SPN(2,SPNI2) = 0.0
*        SPN(3,SPNI2) = 0.0
*      END IF !Because treating the spin of a KS pair makes no sense

*     calculate physical spin from 0<SPIN<1 in the loop
*     Define relative coordinates and velocity components.
      DO 2 K = 1,3
         S1(K) = SPN(K,SPNI1)*M1*M1/CLIGHTPN
         S2(K) = SPN(K,SPNI2)*M2*M2/CLIGHTPN
         KSX(K) = XI(K)-XI(K+3)
         KSV(K) = VI(K)-VI(K+3)
         KSS(K) = S1(K)+S2(K)
         KSSIG(K) = (M1+M2)*(S2(K)/M2-
     &         S1(K)/M1)
         XS(K)  = 0.5*(SPN(K,SPNI1)+SPN(K,SPNI2))
         XA(K)  = 0.5*(SPN(K,SPNI1)-SPN(K,SPNI2))
    2 CONTINUE
*      IF(rank.eq.0)THEN  
*         WRITE(*,*)'S1(K) ,S2(K), S3(K)BUG',S1(1) ,S2(1), 
*     &   S1(2) ,S2(2),S1(3) ,S2(3)
*      END IF
*     diagnostics
*      IF (I1.GT.N.OR.I2.GT.N) WRITE (*,*) 'ALARM!'
      
*     Set speed of light CL and its powers.
      CL = CLIGHTPN
      CL2 = CL*CL
      CL3 = CL2*CL
      CL4 = CL2*CL2
      CL5 = CL4*CL
      CL6 = CL5*CL
      CL7 = CL6*CL
      
      MASAS = M1+M2
      RCSCH = 2.0*MASAS/CL2
      ETA = M1*M2/(MASAS*MASAS)
      DM = M1-M2

      RIJ2 = KSX(1)**2+KSX(2)**2+KSX(3)**2
      RIJ = SQRT(RIJ2)
      MOR = MASAS/RIJ
*     lighter body over heavier body
C      RATIO = M2/M1
C      IF (M2.GT.M1) THEN RATIO = 1/RATIO END IF
      
      V1_V22 = KSV(1)*KSV(1)+KSV(2)*KSV(2)+KSV(3)*KSV(3)
      VWHOLE = SQRT(V1_V22)
*      IF (rank.eq.0) then
*         WRITE(*,*) 'looking for V1_V22:', V1_V22
*      END IF 
      NX = KSX(1)/RIJ
      NY = KSX(2)/RIJ
      NZ = KSX(3)/RIJ

      NVEC(1) = NX
      NVEC(2) = NY
      NVEC(3) = NZ

C      CHI12 = SPN(1,SPNI1)**2+SPN(2,SPNI1)**2+SPN(3,SPNI1)**2
C      CHI22 = SPN(1,SPNI2)**2+SPN(2,SPNI2)**2+SPN(3,SPNI2)**2

*     Dot and Crossproducts needed for Spin-Orbit contribution
      NDV = NX*KSV(1) + NY*KSV(2) + NZ*KSV(3)
      VDS = KSV(1)*KSS(1)+KSV(2)*KSS(2)+KSV(3)*KSS(3)
      VDSIG = KSV(1)*KSSIG(1)+KSV(2)*KSSIG(2)+
     &     KSV(3)*KSSIG(3)
      NDS = NX*KSS(1)+NY*KSS(2)+NZ*KSS(3)
      NDSIG = NX*KSSIG(1)+NY*KSSIG(2)+NZ*KSSIG(3)
      CALL CROSSP(NVEC,KSV,NCV)
      CALL CROSSP(NVEC,KSS,NCS)
      CALL CROSSP(NVEC,KSSIG,NCSIG)
      CALL CROSSP(KSV,KSS,VCS)
      CALL CROSSP(KSV,KSSIG,VCSIG)
      SDNCV = KSS(1)*NCV(1)+KSS(2)*NCV(2)+KSS(3)*NCV(3)
      SIGDNCV = KSSIG(1)*NCV(1)+KSSIG(2)*NCV(2)+KSSIG(3)*NCV(3)
      CALL CROSSP(KSX,KSV,L)
      L(1) = L(1)*MASAS*ETA
      L(2) = L(2)*MASAS*ETA
      L(3) = L(3)*MASAS*ETA
      LABS = SQRT(L(1)*L(1)+L(2)*L(2)+L(3)*L(3))
      LU(1) = L(1)/LABS
      LU(2) = L(2)/LABS
      LU(3) = L(3)/LABS 
      S1DLU = S1(1)*LU(1)+S1(2)*LU(2)+S1(3)*LU(3)
      S2DLU = S2(1)*LU(1)+S2(2)*LU(2)+S2(3)*LU(3)
      NXA = NVEC(1)*XA(1)+NVEC(2)*XA(2)+NVEC(3)*XA(3)
      NXS = NVEC(1)*XS(1)+NVEC(2)*XS(2)+NVEC(3)*XS(3)
      XA2 = XA(1)*XA(1)+XA(2)*XA(2)+XA(3)*XA(3)
      XS2 = XS(1)*XS(1)+XS(2)*XS(2)+XS(3)*XS(3)

      
      
*     RP = RDOT
      RP = (KSX(1)*KSV(1)+KSX(2)*KSV(2)+KSX(3)*KSV(3))/RIJ
      DIFFS(6) = RIJ
      DIFFS(7) = V1_V22
      DIFFS(8) = RP
      IF (APPLY.EQ.0.OR.APPLY.EQ.1) THEN
************1PN*********
         A1 = 2.0*(2.0+ETA)*MOR-(1.0+3.0*ETA)*V1_V22 +1.5*ETA*RP*RP
         B1 = 2.0*(2.0-ETA)*RP
************1PN*********
************2PN*********
         A2 = -0.75*(12.0+29.0*ETA)*MOR*MOR-
     &        ETA*(3.0-4.0*ETA)*V1_V22*V1_V22-
     &        1.875*ETA*(1.0-3.0*ETA)*RP*RP*RP*RP+
     &        0.5*ETA*(13.0-4.0*ETA)*MOR*V1_V22+
     &        (2.0+25.0*ETA+2.0*ETA*ETA)*MOR*RP*RP+
     &        1.5*ETA*(3.0-4.0*ETA)*V1_V22*RP*RP
         B2 = -0.5*RP*((4.0+41.0*ETA+8.0*ETA*ETA)*MOR-
     &        ETA*(15.0+4.0*ETA)*V1_V22+3.0*ETA*(3.0+2.0*ETA)*RP*RP)
************2PN*********
************2.5PN********
         A2_5 = 1.6*ETA*MOR*RP*(17.0*MOR/3.0+3.0*V1_V22)
         B2_5 = -1.6*ETA*MOR*(3.0*MOR+V1_V22)
************2.5PN********
************3PN********
         A3 = MOR*MOR*MOR*(16.0+(1399.0/12.0-41.0*PI2/16.0)*ETA+
     &        71.0*ETA*ETA/2.0)+ETA*(20827.0/840.0+123.0*PI2/64.0-
     &        ETA*ETA)*MOR*MOR*V1_V22-(1.0+(22717.0/168.0+615.0*
     &        PI2/64.0)*ETA+11.0*ETA*ETA/8.0-7.0*ETA*ETA*ETA)*MOR*
     &        MOR*RP*RP-0.25*ETA*(11.0-49.0*ETA+52.0*ETA*ETA)*V1_V22*
     &        V1_V22*V1_V22+35.0*ETA*(1.0-5.0*ETA+5.0*ETA*ETA)*RP*RP*
     &        RP*RP*RP*RP/16.0-0.25*ETA*(75.0+32.0*ETA-40.0*ETA*ETA)*
     &        MOR*V1_V22*V1_V22-0.5*ETA*(158.0-69.0*ETA-60.0*ETA*ETA)*
     &        MOR*RP*RP*RP*RP+ETA*(121.0-16.0*ETA-20.0*ETA*ETA)*MOR*
     &        V1_V22*RP*RP+3.0*ETA*(20.0-79.0*ETA+60.0*ETA*ETA)*V1_V22*
     &        V1_V22*RP*RP/8.0-15.0*ETA*(4.0-18.0*ETA+17.0*ETA*ETA)*
     &        V1_V22*RP*RP*RP*RP/8.0

         B3 = RP*((4.0+(5849.0/840.0+123.0*PI2/32.0)*ETA-25.0*ETA*ETA-
     &        8.0*ETA*ETA*ETA)*MOR*MOR+ETA*(65.0-152.0*ETA-48.0*ETA*
     &        ETA)*V1_V22*V1_V22/8.0+15.0*ETA*(3.0-8.0*ETA-2.0*ETA*
     &        ETA)*RP*RP*RP*RP/8.0+
     &        ETA*(15.0+27.0*ETA+10.0*ETA*ETA)*MOR*V1_V22-ETA*(329.0+
     &        177.0*ETA+108.0*ETA*ETA)*MOR*RP*RP/6.0-
     &        3.0*ETA*(16.0-37.0*ETA-16.0*ETA*ETA)*V1_V22*RP*RP/4.0)
************3PN********
************3.5PN********
         A3_5 = MOR*ETA*(V1_V22*V1_V22*(-366./35.-12.*ETA)+
     &        V1_V22*RP*RP*(114.+12.*ETA)-112.*RP*RP*RP*RP+
     &        MOR*(V1_V22*(-692./35.+724.*ETA/15.)+RP*RP*
     &        (-294./5.-376.*ETA/5.)+MOR*(-3956./35.-184.*ETA/5.)))
         
         B3_5 = 8.0*ETA*MOR*((1325.0+546.0*ETA)*MOR*MOR/42.0+
     &        (313.0+42.0*ETA)*V1_V22*V1_V22/28.0+75.0*RP*RP*RP*RP
     &        -(205.0+777.0*ETA)*MOR*V1_V22/42.0+(205.0+
     &        424.0*ETA)*MOR*RP*RP/12.0-3.0*(113.0+2.0*ETA)*
     &        V1_V22*RP*RP/4.0)/5.
************3.5PN********
*******1.5PN SPIN-ORBIT****
         DO 1 K = 1,3
            C1_5(K) = (NVEC(K)*(12.*SDNCV+6.*DM*SIGDNCV/MASAS)+
     &           9.*NDV*NCS(K)+3.*DM*NDV*NCSIG(K)/MASAS -
     &           7.*VCS(K)-3.*DM*VCSIG(K)/MASAS)/(RIJ2*RIJ)
*******1.5PN SPIN-ORBIT****
*******2 PN SPIN-SPIN******
            C2(K) = -MOR*MOR*MOR/RIJ*3.*ETA*(NVEC(K)*(XS2-XA2-
     &           5.*NXS*NXS+5.*NXA*NXA)+2.*(XS(K)*NXS-XA(K)*NXA))
*******2 PN SPIN-SPIN******
*******2.5PN SPIN-ORBIT****
            C2_5(K) = (NVEC(K)*(SDNCV*(-30.*ETA*NDV*NDV+
     &           24.*ETA*V1_V22-MOR*(38.+25.*ETA))+DM/MASAS*SIGDNCV*
     &           (-15.*ETA*NDV*NDV+12.*ETA*V1_V22-MOR*(18.+14.5*
     &           ETA)))+NDV*KSV(K)*(SDNCV*(-9.+9.*ETA)+DM/MASAS*
     &           SIGDNCV*(-3.+6.*ETA))+NCV(K)*(NDV*VDS*(-3.+3.*ETA)-
     &           8.*MOR*ETA*NDS-DM/MASAS*(4.*MOR*ETA*NDSIG+3.*NDV*
     &           VDSIG))+NDV*NCS(K)*(-22.5*ETA*NDV*NDV+21.*ETA*
     &           V1_V22-MOR*(25.+15.*ETA))+DM/MASAS*NDV*NCSIG(K)*(-15.*
     &           ETA*NDV*NDV+12.*ETA*V1_V22-MOR*(9.+8.5*ETA))+
     &           VCS(K)*(16.5*ETA*NDV*NDV+MOR*(21.+9.*ETA)-14.*ETA*
     &           V1_V22)+DM/MASAS*VCSIG(K)*(9.*ETA*NDV*NDV-7.*ETA*V1_V22
     &           +MOR*(9.+4.5*ETA)))/(RIJ2*RIJ)
*******2.5PN SPIN-ORBIT****
*******2 PN SPIN QUADRUPOLE
C     IF (CHI12.GT.0.AND.CHI22.GT.0.) THEN
C     SQM(K) = -3.*MOR*MOR*MOR/RIJ*ETA/2.*(CHI12/RATIO*((1.-
C     & 5.*(NVEC(1)*SPN(1,I1)+NVEC(2)*SPN(2,I1)+NVEC(3)*
C     & SPN(3,I1))**2/CHI12)*NVEC(K)+2.*(NVEC(1)*SPN(1,I1)+
C     & NVEC(2)*SPN(2,I1)+NVEC(3)*
C     & SPN(3,I1))/CHI12*SPN(K,I1))+
C     & RATIO*CHI22*((1.-
C     & 5.*(NVEC(1)*SPN(1,I2)+NVEC(2)*SPN(2,I2)+NVEC(3)*
C     & SPN(3,I2))**2/CHI22)*NVEC(K)+2.*(NVEC(1)*SPN(1,I2)+
C     & NVEC(2)*SPN(2,I2)+NVEC(3)*
C     & SPN(3,I2))/CHI12*SPN(K,I2)))
C     ELSE
C     SQM(K) = 0.0
C     END IF
*******2 PN SPIN QUADRUPOLE
 1       CONTINUE
         A3_5 = A3_5/CL7
         B3_5 = B3_5/CL7
         A3 = A3/CL6
         B3 = B3/CL6
         A1 = A1/CL2
         B1 = B1/CL2
         A2_5 = A2_5/CL5
         B2_5 = B2_5/CL5
         A2 = A2/CL4
         B2 = B2/CL4        
         KSAK = A2_5+A1+A2+A3+A3_5
         KSBK = B2_5+B1+B2+B3+B3_5
*     GAMREL for switch-on criterion
         GAMREL(1,IPAIR) = SQRT((A1*NX + B1*KSV(1))**2+(A1*NY+B1
     &        *KSV(2))**2+(A1*NZ+B1*KSV(3))**2)
         GAMREL(2,IPAIR) = SQRT((A2*NX + B2*KSV(1))**2+(A2*NY
     &        +B2*KSV(2))**2+(A2*NZ+B2*KSV(3))**2)
         GAMREL(3,IPAIR) = SQRT((A2_5*NX + B2_5*KSV(1))**2+(A2_5*NY
     &        +B2_5*KSV(2))**2+(A2_5*NZ+B2_5*KSV(3))**2)
         GAMREL(4,IPAIR) = SQRT((A3*NX + B3*KSV(1))**2+(A3*NY+B3*KSV(2)
     &        )**2+(A3*NZ+B3*KSV(3))**2)
         GAMREL(5,IPAIR) = SQRT((A3_5*NX + B3_5*KSV(1))**2+(A3_5*NY
     &        +B3_5*KSV(2))**2+(A3_5*NZ+B3_5*KSV(3))**2)
         GAMREL(6,IPAIR) = SQRT(C1_5(1)*C1_5(1)+C1_5(2)*C1_5(2)+
     &        C1_5(3)*C1_5(3))/CL2/MOR*RIJ
         GAMREL(7,IPAIR) = SQRT(C2(1)*C2(1)+C2(2)*C2(2)+C2(3)*C2(3)
     &        )/CL4/MOR*RIJ
         GAMREL(8,IPAIR) = SQRT(C2_5(1)*C2_5(1)+C2_5(2)*C2_5(2)+
     &        C2_5(3)*C2_5(3))/CL4/MOR*RIJ

         A(1) = MOR*(KSAK*NX + KSBK*KSV(1))/RIJ
C     A(1) = A(1)+ C1_5(1)/CL2 + (C2(1)+C2_5(1))/CL4
!     
         A(2) = MOR*(KSAK*NY + KSBK*KSV(2))/RIJ
C     A(2) = A(2)+ C1_5(2)/CL2 + (C2(2)+C2_5(2))/CL4

         A(3) = MOR*(KSAK*NZ + KSBK*KSV(3))/RIJ
C     A(3) = A(3)+ C1_5(3)/CL2 + (C2(3)+C2_5(3))/CL4

         DO 345 K = 1,3
            FRELP(K) = A(K)
 345     CONTINUE
c$$$         WRITE(1123,*)'TIME', TIME, 'A(1:3)', A(1:3), 'MOR', MOR,
c$$$     &   'KSAK', KSAK ,'KSBK', KSBK, 'NX', NX, 'KSV', KSV(1),
c$$$     &   'RIJ', RIJ
c$$$         CALL FLUSH(1123)
      END IF
      IF (APPLY.EQ.1.OR.APPLY.EQ.2) THEN

*     Total acceleration (+ newtonian)
         AT(1) = FRELP(1) - MOR*NX/RIJ
         AT(2) = FRELP(2) - MOR*NY/RIJ
         AT(3) = FRELP(3) - MOR*NZ/RIJ
         
c$$$       IF(RANK.EQ.0)THEN
c$$$         WRITE(*,*) 'AT(1)',AT(1), 'AT(2)',AT(2), 'AT(3)',AT(3)
c$$$       END IF
         RPP = V1_V22/RIJ + AT(1)*NX+AT(2)*NY + AT(3)*NZ - RP*RP/RIJ
***   VA = KSV*VDOT needed when differentiating KSV**2
         VA = AT(1)*KSV(1) + AT(2)*KSV(2) + AT(3)*KSV(3)

*     *For Spin
         DO 4 K = 1,3
            NDOT(K) = (KSV(K)-NVEC(K)*RP)/RIJ
 4       CONTINUE
         NVDOT = NDOT(1)*KSV(1)+NDOT(2)*KSV(2)+NDOT(3)*KSV(3)+
     &        NVEC(1)*AT(1)+NVEC(2)*AT(2)+NVEC(3)*AT(3)

         CALL CROSSP(NDOT,KSV,NDOTCV)
         CALL CROSSP(NVEC,AT,NCA)


*     * FDOT'S

************1PN*********
         A1D = -2.0*(2.0+ETA)*MOR*RP/RIJ - 2.0*(1.0+3.0*ETA)*VA +
     &        3.0*ETA*RP*RPP
         B1D = 2.0*(2.0-ETA)*RPP 
************1PN*********
************2PN*********
         A2D = 1.5*(12.0+29.0*ETA)*MOR*MOR*RP/RIJ - 
     &        ETA*(3.0-4.0*ETA)*4.0*V1_V22*VA - 7.5*ETA*(1.0-3.0*ETA)*
     &        RPP -
     &        0.5*ETA*(13.0-4.0*ETA)*MOR*RP*V1_V22/RIJ+
     &        ETA*(13.0-4.0*ETA)*MOR*VA -
     &        (2.0+25.0*ETA+2.0*ETA*ETA)*MOR*RP*RP*RP/RIJ+
     &        2.0*(2.0+25.0*ETA+2.0*ETA*ETA)*MOR*RP*RPP + 
     &        3.0*ETA*(3.0-4.0*ETA)*VA*RP*RP +
     &        3.0*ETA*(3.0-4.0*ETA)*V1_V22*RP*RPP

         B2D = -0.5*RPP*((4.0+41.0*ETA+8.0*ETA*ETA)*MOR -
     &        ETA*(15.0+4.0*ETA)*V1_V22+3.0*ETA*(3.0+2.0*ETA)*RP*RP) -
     &        0.5*RP*(-(4.0+41.0*ETA+8.0*ETA*ETA)*MOR*RP/RIJ - 
     &        2.0*ETA*(15.0+4.0*ETA)*VA + 6.0*ETA*(3.0+2.0*ETA)*RP*RPP)
************2PN*********
************2.5PN********
         A2_5D = -1.6*ETA*MOR*RP*RP*(17.0/3.0*MOR+3.0*V1_V22)/RIJ +
     &        1.6*ETA*MOR*RPP*(17.0/3.0*MOR+3.0*V1_V22)+
     &        1.6*ETA*MOR*RP*(-17.0*MOR*RP/3.0/RIJ+6.0*VA)
         
         B2_5D = 1.6*ETA*MOR*RP*(3.0*MOR+V1_V22)/RIJ - 
     &        1.6*ETA*MOR*(-3.0*MOR*RP/RIJ+2.0*VA) 
************2.5PN********
************3PN********
         A3D = 6.*ETA*RP*RP*RP*RP*RP*RPP*(35.-175.*ETA+
     &        175.*ETA*ETA)/16. + ETA*(4.*RP*RP*RP*RPP*V1_V22 + 
     &        2.*RP*RP*RP*RP*VA)*(-15.+135.*ETA/2.-255.*ETA*ETA/4.)/2. +
     &        ETA*(2.*RP*RPP*V1_V22*V1_V22+4.*RP*RP*V1_V22*VA)/2.*
     &        (15.-237.*ETA/2.+45.*ETA*ETA) + 6.*V1_V22*V1_V22*VA*ETA*
     &        (-11./4.-49.*ETA/4.-13.*ETA*ETA) + MOR*(4.*RP*RP*RP*RPP*
     &        ETA*(-79.+69./2.*ETA+30.*ETA*ETA) + ETA*(2.*RP*RPP*V1_V22+
     &        2.*RP*RP*VA)*(121.-16.*ETA-20.*ETA*ETA)+4.*V1_V22*VA*ETA*
     &        (-75./4.-8.*ETA+10.*ETA*ETA)) - MOR*RP*((-79.+69.*ETA/
     &        2.+30.*ETA*ETA)*RP*RP*RP*RP*ETA+ETA*RP*RP*V1_V22*(121.-
     &        16.*ETA-20.*ETA*ETA)+ETA*V1_V22*V1_V22*(-75./4.-8.*ETA+
     &        10.*ETA*ETA))/RIJ - 2.*MOR*MOR*RP*(RP*RP*((-1.-615.*PI2*
     &        ETA/64.)-22717.*ETA/168.-11.*ETA*ETA/8.+7.*ETA*ETA*ETA)+
     &        ETA*V1_V22*((20827./840.+123.*PI2/64.)-ETA*ETA))/RIJ + 
     &        MOR*MOR*(2.*RP*RPP*((-1.-615*PI2*ETA/64.)-22717.*ETA/
     &        168.-11.*ETA*ETA/8.+7*ETA*ETA*ETA)+2.*ETA*VA*
     &        ((20827./840. +
     &        123.*PI2/64.)-ETA*ETA)) - 3.*MOR*MOR*MOR*RP*(16.+
     &        (1399./12.-41.*PI2/16.)*ETA+71.*ETA*ETA/2.)/RIJ

         B3D = 75.*RP*RP*RP*RP*RPP*ETA*(3./8.-ETA-.25*ETA*ETA)+
     &        ETA*(3.*RP*RP*RPP*V1_V22+2.*RP*RP*RP*VA)*(-12.+111.*
     &        ETA/4.+12.*ETA*ETA)+ETA*(RPP*V1_V22*V1_V22+4.*RP*V1_V22*
     &        VA)*(65./8.-19.*ETA-6.*ETA*ETA)-MOR*RP*(RP*RP*RP*
     &        ETA*(-329./6.-59.*ETA/2.-18.*ETA*ETA)+RP*V1_V22*ETA*
     &        (15.+27.*ETA+10.*ETA*ETA))/RIJ+MOR*(3.*RP*RP*RPP*ETA*
     &        (-329./6.-59.*ETA/2.-18.*ETA*ETA)+ETA*(RPP*V1_V22+
     &        2.*RP*VA)*(15.+27.*ETA+10.*ETA*ETA))-2.*MOR*MOR*RP*
     &        (RP*((4.+123.*PI2*ETA/32.)+5849.*ETA/840.-25.*ETA*
     &        ETA-8.*ETA*ETA*ETA))/RIJ+MOR*MOR*(RPP*((4.+123.*PI2*
     &        ETA/32.)+5849./840.*ETA-25.*ETA*ETA-8.*ETA*ETA*ETA))
************3PN********
************3.5PN********
         A3_5D = MOR*ETA*(-RP*(V1_V22*V1_V22*(-366./35.-12.*ETA)+
     &        V1_V22*RP*RP*(114.+12.*ETA)+RP*RP*RP*RP*(-112.))/RIJ+
     &        4*V1_V22*VA*(-366./35.-12.*ETA)+2.*(VA*RP*RP+RP*RPP*
     &        V1_V22)*(114.+12.*ETA)+4.*RP*RP*RP*RPP*(-112.)+MOR*
     &        (2*VA*(-692./35.+724.*ETA/15.)+2.*RP*RPP*(-294./5.-
     &        376.*ETA/5.)-2.*RP*(V1_V22*(-692./35.+724.*ETA/15.)+
     &        RP*RP*(-294./5.-376.*ETA/5.))/RIJ-3.*MOR*RP*(-3956./
     &        35.-184.*ETA/5.)/RIJ))

         B3_5D = MOR*ETA*(4.*V1_V22*VA*(626./35.+12.*ETA/5.)+
     &        2.*(VA*RP*RP+V1_V22*RP*RPP)*(-678./5.-12.*ETA/5.)+
     &        4.*RP*RP*RP*RPP*120.-RP*(V1_V22*V1_V22*(626./35.+
     &        12.*ETA/5.)+V1_V22*RP*RP*(-678./5.-12.*ETA/5.)+
     &        120.*RP*RP*RP*RP)/RIJ+MOR*(2.*VA*(-164./21.-148.*ETA/
     &        5.)+2*RP*RPP*(82./3.+848.*ETA/15.)-2.*RP*(V1_V22*
     &        (-164./21-148.*ETA/5.)+RP*RP*(82./3.+848.*ETA/15.))/
     &        RIJ-3.*MOR*RP*(1060./21.+104.*ETA/5.)/RIJ))

************3.5PN********

***********************

*     Simple Spin Euler evolution
         DO 3 K = 1,3
**********1PN**********
            SU1(K) = MOR*ETA*(NVEC(K)*(-4.*VDS-2.*DM/MASAS*VDSIG)+
     &           KSV(K)*(3.*NDS+DM/MASAS*NDSIG)+NDV*(2.*KSS(K)+
     &           DM/MASAS*KSSIG(K)))/RIJ
            
            SV1(K) = MOR*(NVEC(K)*(VDSIG*(-2.+4.*ETA)-2.*DM/MASAS*VDS)+
     &           KSV(K)*(NDSIG*(1.-ETA)+DM/MASAS*NDS)+NDV*(KSSIG(K)*(1.-
     &           2.*ETA)+
     &           DM/MASAS*KSS(K)))/RIJ
**********1PN**********
****  1.5PN SPIN SPIN****
            SS1(K) = 0.5*(L(K)*(4.+3.*(M2/M1))+
     &           (S2(K)-3.*S2DLU*LU(K)))/(RIJ2*RIJ)

            SS2(K) = 0.5*(L(K)*(4.+3.*(M1/M2))+
     &           (S1(K)-3.*S1DLU*LU(K)))/(RIJ2*RIJ)

****  1.5PN SPIN SPIN****
**********2PN**********

            SU2(K) = MOR*ETA/RIJ*(NVEC(K)*(VDS*(-2.*V1_V22+3.*NDV*NDV-
     &           6.*ETA*NDV*NDV+7.*MOR-8.*ETA*MOR)-14.*MOR*NDS*NDV+
     &           DM/MASAS*VDSIG*ETA*(-3.*NDV*NDV-4.*MOR)+DM/MASAS*MOR
     &           *NDSIG*NDV*
     &           (2.-ETA/2.))+KSV(K)*(NDS*(2.*V1_V22-4.*ETA*V1_V22-3.*
     &           NDV*NDV+7.5*ETA*NDV*NDV+4.*MOR-6.*ETA*MOR)+VDS*NDV*(2.-
     &           6.*ETA)+
     &           DM/MASAS*NDSIG*(-1.5*ETA*V1_V22+3.*ETA*NDV*NDV-MOR-
     &           3.5*ETA*
     &           MOR)-3.*DM/MASAS*VDSIG*NDV*ETA)+KSS(K)*NDV
     &           *(V1_V22-2.*ETA*
     &           V1_V22-1.5*NDV*NDV+3.*ETA*NDV*NDV-MOR+2.*ETA*MOR)+
     &           DM/MASAS*KSSIG(K)*NDV*(-ETA*V1_V22+1.5*ETA*NDV*NDV+
     &           (ETA-1.)*MOR))

            SV2(K) = MOR/RIJ*(NVEC(K)*(VDSIG*ETA*(-2.*V1_V22+6.*ETA*NDV*
     &           NDV+(3.+8.*ETA)*MOR)+MOR*NDSIG*NDV*(2.-22.5*ETA+2.*
     &           ETA*ETA)+
     &           DM/MASAS*VDS*ETA*(-3.*NDV*NDV-4.*MOR)+DM/MASAS*
     &           MOR*NDS*NDV*(2.-
     &           0.5*ETA))+KSV(K)*(NDSIG*(0.5*ETA*V1_V22+2.*ETA*ETA*
     &           V1_V22-4.5*ETA*ETA*NDV*NDV+(4.5*ETA-1.+8.*ETA*ETA)*
     &           MOR)+VDSIG*NDV*ETA*(6.*ETA-1.)-3.*DM/MASAS*VDS*NDV*ETA+
     &           DM/MASAS*NDS*(-1.5*ETA*V1_V22+
     &           3.*ETA*NDV*NDV-(1.+3.5*ETA)*MOR))+KSSIG(K)*NDV*(2.*
     &           ETA*ETA*V1_V22-3.*ETA*ETA*NDV*NDV+(-1.+4.*ETA-2.*ETA*
     &           ETA)*MOR)+DM/MASAS*KSS(K)*NDV*(-ETA*V1_V22+1.5*ETA*NDV*
     &           NDV+(-1.+ETA)*MOR))

**********2PN**********
 3       CONTINUE
*     *belongs to 1.5 SPIN SPIN**
         CALL CROSSP(SS1,S1,SS1)
         CALL CROSSP(SS2,S2,SS2)
**************
         DO 7 K = 1,3
            SU3(K) = SU1(K)/CL2 + SU2(K)/CL4 + (SS1(K) + SS2(K))/CL2
            SV(K) = SV1(K)/CL2 + SV2(K)/CL4+MASAS*
     &      (SS2(K)/M2-SS1(K)/
     &           M1)/CL2
            KSS(K) = KSS(K) + SU3(K)*DT
            KSSIG(K) = KSSIG(K) + SV(K)*DT
C     * calculate actual spins from the total S and Sigma
            SPN(K,SPNI1) = M1*(MASAS*KSS(K)-M2*KSSIG(K))
     &           /MASAS/MASAS/
     &           M1/M1*CLIGHTPN

            SPN(K,SPNI2) = M2*(MASAS*KSS(K)+M1*KSSIG(K))
     &           /MASAS/MASAS/
     &           M2/M2*CLIGHTPN
C     derivatives of tagoshi spin variables
            XAD(K) = 0.5/(MASAS*MASAS*M1*M2)*(-SU3(K)
     &           *MASAS*DM-SV(K)*
     &           (M1*M1+M2*M2))
            XSD(K) = 0.5/(MASAS*MASAS*M1*M2)*(SU3(K)
     &            *MASAS*MASAS+SV(K)*
     &           (M1*M1-M2*M2))
 7       CONTINUE

*     *DERIVATIVES OF SPIN TERMS**
*     *Auxiliary terms**
         CALL CROSSP(NDOT,KSS,NDOTCS)
         CALL CROSSP(NVEC,SU3,NCSU)
         CALL CROSSP(NDOT,KSSIG,NDOTCSIG)
         CALL CROSSP(NVEC,SV,NCSV)
         CALL CROSSP(AT,KSS,ACS)
         CALL CROSSP(KSV,SU3,VCSU)
         CALL CROSSP(AT,KSSIG,ACSIG)
         CALL CROSSP(KSV,SV,VCSV)

         SNVDOT = SU3(1)*NCV(1)+SU3(2)*NCV(2)+SU3(3)*NCV(3)+
     &        KSS(1)*NDOTCV(1)+KSS(2)*NDOTCV(2)+KSS(3)*NDOTCV(3)+
     &        KSS(1)*NCA(1)+KSS(2)*NCA(2)+KSS(3)*NCA(3)

         SIGNVDOT = SV(1)*NCV(1)+SV(2)*NCV(2)+SV(3)*NCV(3)+
     &        KSSIG(1)*NDOTCV(1)+KSSIG(2)*NDOTCV(2)+KSSIG(3)*NDOTCV(3)+
     &        KSSIG(1)*NCA(1)+KSSIG(2)*NCA(2)+KSSIG(3)*NCA(3)

         NSDOT = NDOT(1)*KSS(1)+NDOT(2)*KSS(2)+NDOT(3)*KSS(3)+
     &        NVEC(1)*SU3(1)+NVEC(2)*SU3(2)+NVEC(3)*SU3(3)
         NSIGDOT = NDOT(1)*KSSIG(1)+NDOT(2)*KSSIG(2)+NDOT(3)*KSSIG(3)+
     &        NVEC(1)*SV(1)+NVEC(2)*SV(2)+NVEC(3)*SV(3)
         VSDOT = AT(1)*KSS(1)+AT(2)*KSS(2)+AT(3)*KSS(3)+
     &        KSV(1)*SU3(1)+KSV(2)*SU3(2)+KSV(3)*SU3(3)
         VSIGDOT = AT(1)*KSSIG(1)+AT(2)*KSSIG(2)+AT(3)*KSSIG(3)+
     &        KSV(1)*SV(1)+KSV(2)*SV(2)+KSV(3)*SV(3)

         NXSDOT = NDOT(1)*XS(1)+NDOT(2)*XS(2)+NDOT(3)*XS(3)+
     &        NVEC(1)*XSD(1)+NVEC(2)*XSD(2)+NVEC(3)*XSD(3)
         NXADOT = NDOT(1)*XA(1)+NDOT(2)*XA(2)+NDOT(3)*XA(3)+
     &        NVEC(1)*XAD(1)+NVEC(2)*XAD(2)+NVEC(3)*XAD(3)
******end of auxiliary terms
********1.5 PN********
         DO 6 K = 1,3
            C1_5D(K) = -3.*RP/RIJ*C1_5(K)+(NDOT(K)*
     &           (12.*SDNCV+6.*DM/MASAS*
     &           SIGDNCV)+NVEC(K)*(12.*SNVDOT+6.*DM/MASAS*SIGNVDOT)+9.*
     &           NVDOT*NCS(K)+9.*NDV*(NDOTCS(K)+NCSU(K))+3.*DM/MASAS*
     &           (NVDOT*NCSIG(K)+NDV*(NDOTCSIG(K)+NCSV(K)))-7.*(ACS(K)
     &           +VCSU(K))-3.*DM/MASAS*(ACSIG(K)+VCSV(K)))/(RIJ2*RIJ)
*******1.5 PN*********
******2PN SS**********
            C2D(K) = -4.*RP/RIJ*C2(K)-MOR*MOR*MOR*3.*ETA/RIJ*(NDOT(K)*
     &           (XS2-XA2-5.*NXS*NXS+5.*NXA*NXA)+NVEC(K)*(2.*(XS(1)*
     &           XSD(1)+XS(2)*XSD(2)+XS(3)*XSD(3)-XA(1)*XAD(1)-XA(2)*
     &           XAD(2)-XA(3)*XAD(3))-10.*NXS*NXSDOT+10.*NXA*NXADOT)+
     &           2.*(XSD(K)*NXS+XS(K)*NXSDOT-XAD(K)*NXA-XA(K)*NXADOT))
******2PN SS**********

*******2.5 PN*********
            C2_5D(K) = -3.*RP/RIJ*C2_5(K)+(NDOT(K)*(SDNCV*(-30.*ETA*
     &           NDV*NDV+24.*ETA*V1_V22-MOR*(38.+25.*ETA))+DM/MASAS*
     &            SIGDNCV*(-15.*ETA*NDV*NDV+12.*ETA*V1_V22-MOR
     &           *(18.+14.5*ETA)))+
     &           NVEC(K)*(SNVDOT*(-30.*ETA*NDV*NDV+24.*ETA*V1_V22-MOR*
     &           (38.+25.*ETA))+SDNCV*(-60.*ETA*NDV*NVDOT+48.*ETA*VA+
     &           MOR*RP/RIJ*(38.+25.*ETA))+DM/MASAS*SIGNVDOT
     &           *(-15.*ETA*NDV*
     &           NDV+12.*ETA*V1_V22-MOR*(18.+14.5*ETA))
     &           +DM/MASAS*SIGDNCV*
     &           (-30.*ETA*NDV*NVDOT+24.*ETA*VA+MOR*RP/RIJ*(18.+
     &           14.5*ETA)))+(NVDOT*KSV(K)+NDV*AT(K))*
     &           (SDNCV*(-9.+9.*ETA)+DM/MASAS*SIGDNCV*
     &           (-3.+6.*ETA))+NDV*KSV(K)*(SNVDOT*(-9.+9.*ETA)+DM/MASAS*
     &           SIGNVDOT*(-3.+6.*ETA))+(NDOTCV(K)+NCA(K))*(NDV*VDS*
     &           (-3.+3.*ETA)-8.*MOR*ETA*NDS-DM/MASAS*(4.*MOR*ETA*
     &           NDSIG+3.*NDV*VDSIG))+NCV(K)*((NVDOT*VDS+NDV*VSDOT)*
     &           (-3.+3.*ETA)-8.*ETA*MOR*(NSDOT-RP/RIJ*NDS)-
     &           DM/MASAS*(4.*ETA*MOR*(NSIGDOT-RP/RIJ*NDSIG)+
     &           3.*(NVDOT*VDSIG+NDV*VSIGDOT)))+(NVDOT*NCS(K)+NDV*
     &           (NDOTCS(K)+NCSU(K)))*(-22.5*ETA*NDV*NDV+21.*ETA*V1_V22-
     &           MOR*(25.+15.*ETA))+NDV*NCS(K)*(-45.*ETA*NDV*NVDOT+
     &           42.*ETA*VA+MOR*RP/RIJ*(25.+15.*ETA))+DM/MASAS*(NVDOT*
     &           NCSIG(K)+NDV*(NDOTCSIG(K)+NCSV(K)))*(-15.*ETA*NDV*
     &           NDV+12.*ETA*V1_V22-MOR*(9.+8.5*ETA))+DM/MASAS*NDV*
     &           NCSIG(K)*(-30.*ETA*NDV*NVDOT+
     &           24.*ETA*VA+MOR*RP/RIJ*(9.+8.5*ETA))+(ACS(K)+VCSU(K))*
     &           (16.5*ETA*NDV*NDV+MOR*(21.+9.*ETA)-14.*ETA*V1_V22)+
     &           VCS(K)*(33.*ETA*NDV*NVDOT-MOR*RP/RIJ*(21.+9.*ETA)-
     &           28.*ETA*VA)+DM/MASAS*(ACSIG(K)+VCSV(K))
     &           *(9.*ETA*NDV*NDV-
     &           7.*ETA*V1_V22+MOR*(9.+4.5*ETA))+DM/MASAS*VCSIG(K)*(18.*
     &           ETA*NDV*NVDOT-14.*ETA*VA-MOR*RP/RIJ*(9.+4.5*ETA)))/
     &           (RIJ2*RIJ)
*******2.5 PN*********
 6       CONTINUE
*         DEBUGC2D = SQRT(C2D(1)**2+C2D(2)**2+C2D(3)**2)
         A3_5D = A3_5D/CL7
         B3_5D = B3_5D/CL7
         A3D = A3D/CL6
         B3D = B3D/CL6
         A1D = A1D/CL2
         B1D = B1D/CL2
         A2_5D = A2_5D/CL5
         B2_5D = B2_5D/CL5
         A2D = A2D/CL4
         B2D = B2D/CL4

*     Sum up the Adots
         ADK = A2_5D+A1D+A2D+A3D+A3_5D
         BDK = B2_5D+B1D+B2D+B3D+B3_5D

         AD(1) = -2.0*MOR*RP*(KSAK*NX+KSBK*KSV(1))/RIJ2 + 
     &        MOR*(ADK*NX+BDK*KSV(1))/RIJ + MOR*(KSAK*(KSV(1)-NX*RP)
     &        /RIJ+KSBK*AT(1))/RIJ
C     &   + C1_5D(1)/CL2 + C2D(1)/CL4 +
C     &       C2_5D(1)/CL4

         AD(2) = -2.0*MOR*RP*(KSAK*NY+KSBK*KSV(2))/RIJ2 + 
     &        MOR*(ADK*NY+BDK*KSV(2))/RIJ + MOR*(KSAK*(KSV(2)-NY*RP)
     &        /RIJ+KSBK*AT(2))/RIJ
C     & + C1_5D(2)/CL2 + C2D(2)/CL4 +
C     &       C2_5D(2)/CL4

         AD(3) = -2.0*MOR*RP*(KSAK*NZ+KSBK*KSV(3))/RIJ2 + 
     &        MOR*(ADK*NZ+BDK*KSV(3))/RIJ + MOR*(KSAK*(KSV(3)-NZ*RP)
     &        /RIJ+KSBK*AT(3))/RIJ
C     & + C1_5D(3)/CL2 + C2D(3)/CL4 +
C     &       C2_5D(3)/CL4
***********************
         DO 346 K = 1,3
            FRELD(K) = AD(K)
 346     CONTINUE
      END IF

      RETURN

      END
