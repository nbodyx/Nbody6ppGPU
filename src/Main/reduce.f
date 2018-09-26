      SUBROUTINE REDUCE(IESC,JESC,ISUB)
*
*
*       Reduction of chain.
*       -------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,XCM(3),VCM(3),DXC(3),DVC(3),CG(6)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      REAL*8  FIRR(3),FD(3)
      INTEGER JNEW(NMX)
*
*
*       Define discrete time for new polynomials (DT2 < 0 is OK).
      TIME0 = TIME
      DT2 = T0S(ISUB) + TIMEC - TPREV
      DT8 = (TBLOCK - TPREV)/8.0D0
      TIME = TPREV + NINT(DT2/DT8)*DT8
      TIME = MIN(TBLOCK,TIME)
C      TIME = TBLOCK
*       Re-define initial epoch for consistency (ignore phase error).
      T0S(ISUB) = TIME - TIMEC
*
      IF (rank.eq.0.and.KZ(30).GT.2) THEN
          WRITE (6,1)  TIME0+TOFF, TIME+TOFF, DT2, DT8
    1     FORMAT (' REDUCE:   TIME0 TIME DT2 DT8 ',2F12.6,1P,2E10.2)
      END IF
*
C      CALL CONST(XCH,VCH,M,NN,ENERGY,ANG,ALAG)
C      write(6,*) 'ECH ', ENERGY
      
*       Improve coordinates & velocities of c.m. body to order F3DOT.
      CALL XVPRED(ICH,-1)
*
*       Predict X & XDOT of pertubers to order FDOT.
      NB1 = LISTC(1) + 1
      DO 15 L = 2,NB1
          J = LISTC(L)
          CALL XVPRED(J,-2)
   15 CONTINUE

*
C      if(time.gt.0) call XVPRED(IFIRST,NTOT)
C      call adjust
C 
C      DO I =1,NTOT
C         write(333,*) I,NAME(I),BODY(I),X(1:3,I),XDOT(1:3,I)
C      END DO
C      DO L = 1,NCH
C         write(6,*) L, BODYC(L),
C     &        X4(1:3,L)+X(1:3,ICH), XDOT4(1:3,L)+XDOT(1:3,ICH)
C      END DO
C      RI = 0
C      RI2 = 0
C      BCM = BODYC(IESC) + BODYC(JESC)
C      DO K=1,3
C         XCM(K) = (BODYC(IESC)*X4(K,IESC) + BODYC(JESC)*X4(K,JESC))/BCM
C         RI2 = RI2 + XCM(K)**2
C      END DO
C      RI = SQRT(RI2)
C      write(6,*) 'IESC',IESC,'JESC',JESC,'RI',RI,'XC',XCM+X(1:3,ICH)
C      XCM = X4(1:3,IESC)-X4(1:3,JESC)
C      RI2 = XCM(1)**2 + XCM(2)**2 + XCM(3)**2
C      write(6,*) 'POIb',BODYC(IESC)*BODYC(JESC)/SQRT(RI2)
C      EBINBK=0
C      DO K=1,3
C         XCM(K) =
C     &    (BODYC(IESC)*XDOT4(K,IESC) + BODYC(JESC)*XDOT4(K,JESC))/BCM
C         EBINBK = EBINBK + 0.5*((XDOT4(K,IESC)-XCM(K))**2*BODYC(IESC)+
C     &    (XDOT4(K,JESC)-XCM(K))**2*BODYC(JESC))
C      END DO
C      write(6,*) 'Ekinb',EBINBK,'EBB',
C     &     EBINBK-BODYC(IESC)*BODYC(JESC)/SQRT(RI2)
C      
C      XCM(1:3) = 0
C      BCM = 0
C      DO KK =1,3
C         XCM(1:3) = XCM(1:3) + BODYC(KK)*X4(1:3,KK)
C         BCM  = BCM + BODYC(KK)
C      END DO
C      XCM(1:3) = XCM(1:3)/BCM
C      RI2 = XCM(1)**2 + XCM(2)**2 + XCM(3)**2
C      write(6,*) 'Other RI',SQRT(RI2),'XC',XCM+X(1:3,ICH)
C      
CC      write(6,*) 'LISTC ',NB1, NAME(LISTC(2)), NAME(LISTC(3))
C      CALL CONST(XCH,VCH,M,NN,ENERGY,ANG,ALAG)
C      EKINCC = 0.5*BODY(ICH)*(XDOT(1,ICH)**2+XDOT(2,ICH)**2
C     &     +XDOT(3,ICH)**2)
C      write(6,*) 'ECH ', ENERGY,'EKIN',EKINCC,'ET',ENERGY+EKINCC
C 
C      call adjust
*       Make new chain from remaining members.
    2 LK = 0
C      DPOT = 0.0
      DO 10 L = 1,NCH
          IF (L.EQ.IESC) GO TO 10
          RIJ2 = 0.0
          DO 5 K = 1,3
              LK = LK + 1
              XCH(LK) = X4(K,L)
              VCH(LK) = XDOT4(K,L)
              RIJ2 = RIJ2 + (XCH(LK) - X4(K,IESC))**2
 5         CONTINUE
C          DPOT = DPOT + BODYC(L)*BODYC(IESC)/SQRT(RIJ2)
   10 CONTINUE
*
*       Reduce chain membership and mass (global & local COMMON).
      NCH = NCH - 1
      NN = NCH
      MASS = MASS - M(IESC)
*
*
*       Set new c.m. for reduced system and save old c.m. variables.
      DO 20 K = 1,3
          DXC(K) = -BODYC(IESC)*X4(K,IESC)/(BODY(ICH) - BODYC(IESC))
          DVC(K) = -BODYC(IESC)*XDOT4(K,IESC)/(BODY(ICH) - BODYC(IESC))
          XCM(K) = X(K,ICH) + DXC(K)
          VCM(K) = XDOT(K,ICH) + DVC(K)
          CM(K) = X(K,ICH)
          CM(K+3) = XDOT(K,ICH)
C          write(6,*) 'K',K,' XCM',XCM(K),' VCM',VCM(K),
C     &         ' XCM[old]',CM(K),' VCM[old]',CM(K+3)
   20 CONTINUE
*
*       Re-define new chain variables w.r. to modified c.m. (NB! retain X4).
      LK = 0
      DO 30 L = 1,NCH
          DO 25 K = 1,3
              LK = LK + 1
              XCH(LK) = XCH(LK) - DXC(K)
              VCH(LK) = VCH(LK) - DVC(K)
   25     CONTINUE
   30 CONTINUE
*
*       Save original mass & neighbour radius of c.m. body.
      BODYCH = BODY(ICH)
C      write(6,*) 'BODY CM', BODY(ICH)
      RS0 = RS(ICH)
*
*       Search for global index of escaper.
      DO 40 J = IFIRST,NTOT
          IF (NAME(J).EQ.NAMEC(IESC)) THEN
              I = J
              IF (rank.eq.0.and.BODY(J).GT.0.0D0) 
     &        WRITE (6,35)  I, IESC, NAMEC(IESC)
   35         FORMAT (' WARNING!   NON-ZERO GHOST    I IESC NAMEC ',3I5)
              GO TO 55
          END IF
   40 CONTINUE
*
*       Switch to another reference body since #ICH is escaping (NAME = 0).
      I = ICH
      IF (IESC.GT.1) THEN
          NEW = 1
      ELSE
          NEW = 2
      END IF
*
*       Identify global index of new reference body (including binary c.m.).
      DO 45 J = IFIRST,NTOT
          IF (NAME(J).EQ.NAMEC(NEW)) THEN
              ICH = J
C*     Replace old c.m. index to new index in NXTLST
C              call exchange_tlist(I,ICH,STEP,DTK)
              GO TO 50
          END IF
   45 CONTINUE
*
*       Include warning if no reference body (this should not occur).
      if(rank.eq.0) WRITE (6,48)  IESC, NAMEC(NEW)
   48 FORMAT (' REDUCE:   DANGER!   NO REFERENCE BODY    IESC NAME',2I5)
      NCH = NCH + 1
      NN = NCH
      MASS = MASS + M(IESC)
      GO TO 100
*
*       Copy neighbour list of old c.m. body for routine NBREST.
   50 NNB = LIST(1,I)
      DO 52 L = 2,NNB+1
          JPERT(L-1) = LIST(L,I)
   52 CONTINUE
*
*       Restore ghost to lists of neighbours (body #I will be skipped).
C      JLIST(1) = ICH
      JLIST(1) = I
      JNEW(1) = I
      JNEW(2) = ICH
C      CALL NBREST(I,1,NNB)
      CALL NBCHANGE(JLIST,1,JNEW,2,LIST,1,NTOT)
      JLIST(1) = ICH
*
      IF (rank.eq.0.and.KZ(30).GT.1) THEN
          WRITE (6,53)  NAME0, NAME(ICH), ICH
   53     FORMAT (' REDUCE:    SWITCH C.M.    NAME0 NAMECH ICH ',3I5)
      END IF
*
*       Exchange name of reference body and initialize new c.m. name.
      NAME(I) = NAME0
      NAME0 = NAME(ICH)
      NAME(ICH) = 0
*
*       Update total mass and initialize new c.m. body variables.
   55 BODY(ICH) = BODYCH - BODYC(IESC)
      CM(7) = BODY(ICH)
      T0(ICH) = TIME
      T0R(ICH) = TIME
      DO 60 K = 1,3
          X(K,ICH) = XCM(K)
          X0(K,ICH) = XCM(K)
          XDOT(K,ICH) = VCM(K)
          X0DOT(K,ICH) = VCM(K)
   60 CONTINUE
*
      
*       Restore the mass and transform to global coordinates & velocities.
      BODY(I) = BODYC(IESC)
      T0(I) = TIME
      T0R(I) = TIME
C      DKP = 0.0
C      DECM = 0.0
C      DKE = 0.0
      DO 65 K = 1,3
          X(K,I) = X4(K,IESC) + CM(K)
          XDOT(K,I) = XDOT4(K,IESC) + CM(K+3)
          X0(K,I) = X(K,I)
          X0DOT(K,I) = XDOT(K,I)
*       Form kinetic energy terms (definition of DVC needs negative sign).
C          DKP = DKP - BODYC(IESC)*DVC(K)*XDOT4(K,IESC)
C          DECM = DECM + 0.5*BODY(ICH)*DVC(K)**2
C          DKE = DKE + 0.5*BODYC(IESC)*XDOT4(K,IESC)**2
   65 CONTINUE
*
*       Remove chain (and clump) mass & reference name of escaper.
      DO 70 L = IESC,NCH
          M(L) = M(L+1)
          BODYC(L) = BODYC(L+1)
          NAMEC(L) = NAMEC(L+1)
          SIZE(L) = SIZE(L+1)
          ISTAR(L) = ISTAR(L+1)
          BODYS(L,ISUB) = BODYS(L+1,ISUB)
          NAMES(L,ISUB) = NAMES(L+1,ISUB)
   70 CONTINUE
*
*       Initialize neighbour list using current radius of (old) #ICH.
      CALL NBLIST(ICH,RS0)
*
*       Set dominant F & FDOT on body #I for #ICH in FPOLY2.
      JLIST(1) = ICH
      JLIST(2) = JCLOSE
      NP = 2
      IF (JCLOSE.EQ.0) NP = 1
      CALL FCLOSE(I,NP)
*
*       Perform re-initialization of c.m. polynomials & perturber list.
*     Avoid call twice when binary escape
      IF(JESC.LE.0) CALL REINIT(ISUB)
*
*       Copy new chain coordinates & velocities to standard variables.
      LK = 0
      DO 80 L = 1,NCH
          DO 75 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
   75     CONTINUE
   80 CONTINUE
*
*       Update total energy using regular expression (checked OK).
C      ENERGY = ENERGY + DECM - DKP + DPOT - DKE
C      WRITE (6,81)  DECM, DKP, DPOT, DKE, EN0-ENERGY, ENERGY
C   81 FORMAT (' CHECK!   DECM DKP DPOT DKE DE ENER  ',6F10.6)
      
*       Copy neighbour list of cm. body for routine NBREST.
      NNB = LIST(1,ICH)
      DO 82 L = 2,NNB+1
          JPERT(L-1) = LIST(L,ICH)
   82 CONTINUE
*
*       Restore ghost to lists of neighbours (body #ICH will be skipped).
C      JLIST(1) = I
C      CALL NBREST(ICH,1,NNB)
      JLIST(1) = ICH
      JNEW(1) = ICH
      JNEW(2) = I
      CALL NBCHANGE(JLIST,1,JNEW,2,LIST,1,NTOT)
      JLIST(1) = I
*
*       Make new neighbour list.
      CALL NBLIST(I,RS0)
*
*       Distinguish between single particle and binary (JESC = 0 & > 0).
      IF (JESC.EQ.0) THEN
*       Initialize force polynomials & time-steps (add differential F & FD).
          IPHASE = 8
*     Remove from ghost list
          call delay_remove_tlist(I,STEP,DTK)
          CALL FPOLY1(I,I,0)
          DO 85 K = 1,3
              FIRR(K) = 0.0
              FD(K) = 0.0
   85     CONTINUE
          CALL FCHAIN(I,0,X(1,I),XDOT(1,I),FIRR,FD)
          RIJ2 = 0.0
          VIJ2 = 0.0
          RDOT = 0.0
          DO 86 K = 1,3
              F(K,I) = F(K,I) + FIRR(K)
              FDOT(K,I) = FDOT(K,I) + FD(K)
              RIJ2 = RIJ2 + (X(K,I) - X(K,ICH))**2
              VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,ICH))**2
              RDOT = RDOT + (X(K,I) - X(K,ICH))*(XDOT(K,I)-XDOT(K,ICH))
   86     CONTINUE
          CALL FPOLY2(I,I,0)
*     Add escaper into NLSTDELAY
          call delay_store_tlist(I,STEP,DTK)
*       Check high-velocity ejection for outwards hyperbolic motion.
          HI = 0.5*VIJ2 - (BODY(I) + BODY(ICH))/SQRT(RIJ2)
          IF (HI.GT.0.0.AND.RDOT.GT.0.0.AND.KZ(37).GT.0) THEN
               CALL HIVEL(I)
          END IF
          IPHASE = -1
*       Include case of escaping binary (set JESC < 0 for new KS).
      ELSE IF (JESC.GT.0) THEN
          ICLOSE = I
          IF (JESC.GT.IESC) JESC = JESC - 1
          IESC = JESC
          JESC = -1
          GO TO 2
      ELSE
*       Initialize KS regularization after second reduction.
          ICOMP = MIN(ICLOSE,I)
          JCOMP = MAX(ICLOSE,I)
C          write(6,*) 'T', TIME
C          DO L = 1,NCH
C             write(6,*) L, BODYC(L),X4(1:3,L)+X(1:3,ICH),
C     &            XDOT4(1:3,L)+XDOT(1:3,ICH)
C          END DO
C          write(6,*) 4, BODY(ICOMP), X(1:3,ICOMP), XDOT(1:3,ICOMP)
C     &         , T0(ICOMP), T0R(ICOMP)
C          write(6,*) 5, BODY(JCOMP), X(1:3,JCOMP), XDOT(1:3,JCOMP)
C     &         , T0(JCOMP), T0R(JCOMP)
C          write(6,*) 'LISTC ',NB1, NAME(LISTC(2)), NAME(LISTC(3))
C          CALL CONST(XCH,VCH,M,NN,ENERGY,ANG,ALAG)
C          write(6,*) 'ECH ', ENERGY
         CALL KSREG
C          CALL CONST(XCH,VCH,M,NN,ENERGY2,ANG,ALAG)
C          EBINNEW= BODY(7755)*BODY(7756)*H(3878)/BODY(21190)
C          EKINCC2 = 0.5*BODY(ICH)*(XDOT(1,ICH)**2+XDOT(2,ICH)**2
C     &         +XDOT(3,ICH)**2)
C          ICC = 21190
C          EKINBB = 0.5*BODY(ICC)*(XDOT(1,ICC)**2+XDOT(2,ICC)**2
C     &         +XDOT(3,ICC)**2)
C          write(6,*) 'ECH ', ENERGY2,'EB2',EBINNEW,'EKINC',EKINCC2,
C     &         'EKINB',EKINBB,
C     &         'ENEW',ENERGY2+EBINNEW+EKINCC2+EKINBB,
C     &         'DE',ENERGY2+EBINNEW+EKINCC2+EKINBB-ENERGY-EKINCC
C            DO kk =1,NTOT
C               write(334,*) kk,NAME(kk),BODY(kk),X(1:3,kk),XDOT(1:3,kk)
C            END DO
C          write(6,*) 'LISTC ',NB1, NAME(LISTC(2)), NAME(LISTC(3))
C          call adjust
      END IF
*

*       Re-activate any dormant binary.
      IF (I.GT.N.AND.JESC.EQ.0) THEN
          CALL RENEW(I)
      END IF
*
*       Prepare c.m. check.
      DO 88 K = 1,6
          CG(K) = 0.0
   88 CONTINUE
*
      LK = 0
      DO 95 L = 1,NCH
          DO 90 K = 1,3
              LK = LK + 1
              CG(K) = CG(K) + BODYC(L)*XCH(LK)
              CG(K+3) = CG(K+3) + BODYC(L)*VCH(LK)
   90     CONTINUE
   95 CONTINUE
*
      DO 96 K = 1,6
          CG(K) = CG(K)/BODY(ICH)
   96 CONTINUE
*
      IF (rank.eq.0.and.KZ(30).GT.2) THEN
          WRITE (6,97)  NCH,TIME+TOFF, (CG(K),K=1,6), XDOT(1:3,ICH)
   97     FORMAT (' REDUCE: NCH  T CG Vcm',I4,F12.5,1P,9E9.1)
      END IF
*
*       Ensure current coordinates & velocities for chain components.
      CALL XCPRED(1)
*
  100 RETURN
*
      END

