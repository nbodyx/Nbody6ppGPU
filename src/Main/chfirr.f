      SUBROUTINE CHFIRR(I,IR,XI,XIDOT,FIRR,FD)
*
*
*       Irregular force & derivative on chain c.m.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  M,MASS,MC,MIJ,MKK
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XJ(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  XI(3),XIDOT(3),DX(3),DV(3),FIRR(3),FD(3),FP(3*NCMAX),
     &        FPSUM(3),FPD(3*NCMAX),FPDSUM(3)
*
*
*     --07/31/14 23:23-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      print*,'CHFIRR I FIRR',I,FIRR(1:3)
*     --07/31/14 23:23-lwang-end----------------------------------------*
*       Check updating perturber list at frequent intervals.
*      IF (IR.GT.0) THEN
*          CALL CHLIST(ICH)
*       Note XC & UC for chain members predicted at new block-time.
*      END IF
*
*       First subtract all perturber forces with respect to single body #ICH.
      NPC = LISTC(1) + 1
      JDUM = 0

      RPERT2 = CMSEP2*(0.5*RSUM)**2
      DO 30 L = 2,NPC
          J = LISTC(L)
          call jpred_int(J,TIME)
          DX(1) = X(1,J) - XI(1)
          DX(2) = X(2,J) - XI(2)
          DX(3) = X(3,J) - XI(3)
          RIJ2 = DX(1)**2 + DX(2)**2 + DX(3)**2
*       Skip correcting all interactions outside RPERT2.
          IF (RIJ2.GT.RPERT2) GO TO 30
*     --03/13/14 21:07-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          write(130+rank,*),'J',J,'N',NAME(J),'M',BODY(J),'X',X(1:3,J),
c$$$     &         'XD',XDOT(1:3,J),'I',I,'M',BODY(I),'X',XI(1:3),'XD',
c$$$     &         XIDOT(1:3),'T',TIME
c$$$          if(rank.eq.0) then
c$$$          end if
c$$$          if(time.eq.29.047119140625000) then
c$$$             write(130+rank,*),'j',j,'n',name(j),'xd',xdot(1,j),
c$$$     &            'npc',npc,'t',time,'x0d',x0dot(1,j),'f',f(1,j),
c$$$     &            'fd',fdot(1,j)
c$$$          call flush(130+rank)
c$$$          end if
*     --03/13/14 21:07-lwang-end----------------------------------------*
          IF (J.LE.N) GO TO 5
          JP = J - N
          J1 = 2*JP - 1
*     --07/25/14 23:34-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          print*,rank,'CM J1',J1,'X1',X(1:3,J1),'XD',XDOT(1:3,J1),
c$$$     &         'J2',X(1:3,J1+1),'XD',XDOT(1:3,J1+1),'XI',XI(1:3)
*     --07/25/14 23:34-lwang-end----------------------------------------*
*       See whether c.m. approximation applies (skip unperturbed case).
          IF (RIJ2.GT.CMSEP2*R(JP)**2.OR.LIST(1,J1).EQ.0) GO TO 5
*       Resolve components of pair #JP (suppressed 2/9/98).
*         CALL KSRES2(JP,J1,J2,RIJ2)
*       Note danger of low-order prediction (high-order KSRES2 in INTGRT).
          J = J1
          JDUM = J1
*       Sum over individual KS components.
    5     DR2 = 0.0
*     --03/13/14 21:07-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          print*,'CHFIRR I FIRR',I,FIRR(1:3),'J',J,'N',NAME(J),
c$$$     &         'NP',LIST(1,2*(J-N)-1)
c$$$          if(time.eq.29.047119140625000) then
c$$$             write(130+rank,*),'j',j,'n',name(j),'xd',xdot(1,j),
c$$$     &            'npc',npc,'t',time,'x0d',x0dot(1,j),'f',f(1,j),
c$$$     &            'fd',fdot(1,j)
c$$$             call flush(130+rank)
c$$$          end if
*     --03/13/14 21:07-lwang-end----------------------------------------*
          DRDV = 0.0
          DO 10 K = 1,3
              DX(K) = X(K,J) - XI(K)
              DV(K) = XDOT(K,J) - XIDOT(K)
              DR2 = DR2 + DX(K)**2
              DRDV = DRDV + DX(K)*DV(K)
   10     CONTINUE
          DR2I = 1.0/DR2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*DRDV*DR2I
*
*       Subtract force & first derivative.
          DO 20 K = 1,3
              FIRR(K) = FIRR(K) - DX(K)*DR3I
              FD(K) = FD(K) - (DV(K) - DX(K)*DRDV)*DR3I
   20     CONTINUE
*     --07/25/14 23:44-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          if (rank.eq.0) then
c$$$          print*,'J',J,'DX',DX(1:3),'DR3I',DR3I,'FIRR',FIRR(1:3),
c$$$     &         'M',BODY(J),'MI',BODY(I)
c$$$          call flush(6)
c$$$          end if
c$$$          if(rank.eq.0.and.time.gt.23.17335927486419678) stop
*     --07/25/14 23:44-lwang-end----------------------------------------*
          IF (J.EQ.JDUM) THEN
              J = J + 1
              GO TO 5
          END IF
   30 CONTINUE
*
*       Initialize the perturbing force components and vectorial sum.
      DO 35 K = 1,3*NCH
          FP(K) = 0.0D0
          FPD(K) = 0.0D0
   35 CONTINUE
      DO 40 K = 1,3
          FPSUM(K) = 0.0D0
          FPDSUM(K) = 0.0D0
   40 CONTINUE
*
*       Set perturber distance and calculate force for chain c.m.
      IM1 = 0
      DO 70 LL = 2,NPC
*       Reset dummy index after use (otherwise bug with two KS pairs).
          KDUM = 0
          K = LISTC(LL)
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
*       Skip perturbers consistent with c.m. approximation test above.
          IF (RIJ2.GT.RPERT2) GO TO 70
*     --07/31/14 23:24-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          print*,'70 CHFIRR I FIRR',I,FIRR(1:3),'K',K,'N',NAME(K)
*     --07/31/14 23:24-lwang-end----------------------------------------*
*       Decide appropriate summation (c.m. approximation or components).
          IF (K.GT.N) then
             if (RIJ2.LT.CMSEP2*R(K-N)**2) THEN
                KDUM = 2*(K - N) - 1
                K = KDUM
                GO TO 55
             end if
          END IF
*
*       Obtain perturbation on each component.
          DO 50 IM = 1,NCH
              IM1 = 3*(IM - 1)
   46         DR2 = 0.0
              DRDV = 0.0
              DO 48 L = 1,3
                  DX(L) = X(L,K) - XC(L,IM)
                  DV(L) = XDOT(L,K) - UC(L,IM)
                  DR2 = DR2 + DX(L)**2
                  DRDV = DRDV + DX(L)*DV(L)
   48         CONTINUE
              DR2I = 1.0/DR2
              DR3I = BODY(K)*DR2I*SQRT(DR2I)
              DRDV = 3.0*DRDV*DR2I
*     --08/01/14 19:40-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$              print*,'48 DR2I',DR2I,'IM',IM,'K',K,'KDUM',KDUM
*     --08/01/14 19:40-lwang-end----------------------------------------*
*
*       Accumulate perturbing force & derivative.
              DO 49 L = 1,3
                  FP(IM1+L) = FP(IM1+L) + DX(L)*DR3I
                  FPD(IM1+L) = FPD(IM1+L) + (DV(L) - DX(L)*DRDV)*DR3I
   49         CONTINUE
*
              IF (K.EQ.KDUM) THEN
                  K = K + 1
                  GO TO 46
              ELSE
                  IF (KDUM.GT.0) K = KDUM
              END IF
   50     CONTINUE
          GO TO 70
*
*       Sum over individual components of pair #J using c.m. approximation.
   55     DR2 = 0.0
          DRDV = 0.0
          DO 60 L = 1,3
              DX(L) = X(L,K) - XI(L)
              DV(L) = XDOT(L,K) - XIDOT(L)
              DR2 = DR2 + DX(L)**2
              DRDV = DRDV + DX(L)*DV(L)
   60     CONTINUE
          DR2I = 1.0/DR2
          DR3I = BODY(K)*DR2I*SQRT(DR2I)
          DRDV = 3.0*DRDV*DR2I
*
*       Add force & first derivative.
          DO 65 L = 1,3
              FIRR(L) = FIRR(L) + DX(L)*DR3I
              FD(L) = FD(L) + (DV(L) - DX(L)*DRDV)*DR3I
   65     CONTINUE
          IF (K.EQ.KDUM) THEN
              K = K + 1
              GO TO 55
          END IF
   70 CONTINUE
*
*     --07/31/14 23:39-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          print*,'IM1 FIRR',FIRR(1:3)
*     --07/31/14 23:39-lwang-end----------------------------------------*
*       Add perturbations on the components to the c.m. force & derivative.
      IF (IM1.GT.0) THEN
*       Sum individual perturbations.
          DO 80 L = 1,NCH
              L1 = 3*(L - 1)
*     --07/31/14 23:48-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$              print*,'80 FPSUM',FPSUM(1:3),'BODYC',BODYC(L),'FP',
c$$$     &             (FP(L1+K),K=1,3)
*     --07/31/14 23:48-lwang-end----------------------------------------*
              DO 75 K = 1,3
                  FPSUM(K) = FPSUM(K) + BODYC(L)*FP(L1+K)
                  FPDSUM(K) = FPDSUM(K) + BODYC(L)*FPD(L1+K)
   75         CONTINUE
   80     CONTINUE
*
*       Include mass-weighted contributions.
          BODYIN = 1.0/BODY(I)
          DO 90 K = 1,3
              FIRR(K) = FPSUM(K)*BODYIN + FIRR(K)
              FD(K) = FPDSUM(K)*BODYIN + FD(K)
   90     CONTINUE
      END IF
*
*     --07/25/14 23:22-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(rank.eq.0) then
c$$$      print*,'F - CHFIRR I FIRR',I,FIRR(1:3),'J',J,'N',NAME(J)
c$$$      end if
*     --07/25/14 23:22-lwang-end----------------------------------------*

      RETURN
*
      END
