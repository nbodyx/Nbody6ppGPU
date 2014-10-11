      SUBROUTINE KSREG
*
*
*       New KS regularization.
*       ----------------------
*
      INCLUDE 'common6.h'
      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
      REAL*8 TPRED,TRES,DSAVE(17)
      LOGICAL ipredall
*      EXTERNAL RENAME
*
*
*       Replace #JCOMP by arbitrary body in case it is the only neighbour.
      NNB = LIST(1,ICOMP)
*     --03/05/14 16:21-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      print*,rank,'KSREG ICOMP',ICOMP,'JCOMP',JCOMP,'NB',NNB,'t',time
c$$$      call flush(6)
*     --03/05/14 16:21-lwang-end----------------------------------------*
      IF (NNB.EQ.1.AND.LIST(2,ICOMP).EQ.JCOMP) THEN
          LIST(2,ICOMP) = JCOMP + 1
          IF (JCOMP + 1.GT.NTOT) THEN
              LIST(2,ICOMP) = MAX(ICOMP-1,IFIRST+2)
          END IF
      END IF
*
*       Copy neighbour list of #ICOMP without JCOMP.
      NNB1 = 1
      DO 1 L = 1,NNB
          IF (LIST(L+1,ICOMP).EQ.JCOMP) GO TO 1
          NNB1 = NNB1 + 1
          ILIST(NNB1) = LIST(L+1,ICOMP)
    1 CONTINUE
      ILIST(1) = NNB1 - 1
*
*       Save basic variables for components unless in correct location.
      DO 10 KCOMP = 1,2
*       Treat the first & second component in turn.
          IF (KCOMP.EQ.1) THEN
              I = ICOMP
          ELSE
              I = JCOMP
          END IF
          J = 2*NPAIRS + KCOMP
          IF (I.EQ.J) GO TO 10
*     
*     Exchange NXTLST index also
          IF(IPHASE.EQ.8) THEN
             call delay_remove_tlist(J,STEP,DTK)
             call delay_store_tlist(I)
          ELSE
             IF(TIME+TOFF.NE.0.D0) call exchange_tlist(I,J,STEP,DTK)
          END IF
*
          call jpred(I,TIME,TIME)
          DO 2 K = 1,3
              DSAVE(K) = X(K,I)
              DSAVE(K+3) = X0DOT(K,I)
    2     CONTINUE
          DSAVE(7) = BODY(I)
          DSAVE(8) = RS(I)
          DSAVE(9) = RADIUS(I)
          DSAVE(10) = TEV(I)
          DSAVE(11) = BODY0(I)
          DSAVE(12) = TEV0(I)
          DSAVE(13) = EPOCH(I)
          DSAVE(14) = SPIN(I)
          DSAVE(15) = ZLMSTY(I)
          DSAVE(16) = TPRED(I)
          DSAVE(17) = STEP(I)
          NAMEI = NAME(I)
          KSI = KSTAR(I)
*
*       Exchange first & second single particle with ICOMP & JCOMP.
          DO 4 K = 1,3
              X(K,I) = X(K,J)
              X0(K,I) = X0(K,J)
              X0DOT(K,I) = X0DOT(K,J)
              XDOT(K,I) = XDOT(K,J)
              F(K,I) = F(K,J)
              FDOT(K,I) = FDOT(K,J)
              FI(K,I) = FI(K,J)
              FIDOT(K,I) = FIDOT(K,J)
              D0(K,I) = D0(K,J)
              D1(K,I) = D1(K,J)
              D2(K,I) = D2(K,J)
              D3(K,I) = D3(K,J)
              FR(K,I) = FR(K,J)
              FRDOT(K,I) = FRDOT(K,J)
              D0R(K,I) = D0R(K,J)
              D1R(K,I) = D1R(K,J)
              D2R(K,I) = D2R(K,J)
              D3R(K,I) = D3R(K,J)
              X(K,J) = DSAVE(K)
              X0DOT(K,J) = DSAVE(K+3)
    4     CONTINUE
*
          TPRED(I) = TPRED(J)
          BODY(I) = BODY(J)
          RS(I) = RS(J)
          RADIUS(I) = RADIUS(J)
          TEV(I) = TEV(J)
          TEV0(I) = TEV0(J)
          BODY0(I) = BODY0(J)
          EPOCH(I) = EPOCH(J)
          SPIN(I) = SPIN(J)
          ZLMSTY(I) = ZLMSTY(J)
          NAME(I) = NAME(J)
          KSTAR(I) = KSTAR(J)
          STEP(I) = STEP(J)
          STEPR(I) = STEPR(J)
          T0(I) = T0(J)
          TIMENW(I) = T0(I) + STEP(I)
          T0R(I) = T0R(J)
          K = LIST(1,J) + 1
          DO 5 L = 1,K
              LIST(L,I) = LIST(L,J)
    5     CONTINUE
          BODY(J) = DSAVE(7)
          RS(J) = DSAVE(8)
          RADIUS(J) = DSAVE(9)
          TEV(J) = DSAVE(10)
          BODY0(J) = DSAVE(11)
          TEV0(J) = DSAVE(12)
          EPOCH(J) = DSAVE(13)
          SPIN(J) = DSAVE(14)
          ZLMSTY(J) = DSAVE(15)
          TPRED(J) = DSAVE(16)
          NAME(J) = NAMEI
          KSTAR(J) = KSI
          STEP(J) = DSAVE(17)

#ifdef SIMD
*     Update Particle data in AVX/SSE library
          IF(TIME+TOFF.NE.0.0D0) THEN
             CALL IRR_SIMD_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),
     &            FDOT(1,I),BODY(I),T0(I))
             CALL IRR_SIMD_SET_LIST(I,LIST(1,I))
          END IF
#endif

 10    CONTINUE

*     Reset TRES for ksres
      TRES(NPAIRS+1) = 0.D0
      
*       Save neighbour list of first component for RENAME & FPOLY.
      DO 15 L = 1,NNB1
          LIST(L,IFIRST) = ILIST(L)
   15 CONTINUE
*
*       Increase pair index, total number & single particle index.
      NPAIRS = NPAIRS + 1
      NTOT = N + NPAIRS
      IFIRST = 2*NPAIRS + 1
*
*       Update all relevant COMMON list arrays.
      CALL RENAME_KS

*      Update the tlist array in intgrt
      IF(IPHASE.NE.8.AND.time+TOFF.ne.0.0D0) then
         call delay_remove_tlist(IFIRST-2,STEP,DTK)
         call delay_remove_tlist(IFIRST-1,STEP,DTK)
         call shrink_tlist
      end if
*
*       Copy neighbour list for second component & c.m. (NNB1 = LIST(1,I)+1).
      I = 2*NPAIRS - 1
      DO 40 L = 1,NNB1
          LIST(L,I+1) = LIST(L,I)
          LIST(L,NTOT) = LIST(L,I)
   40 CONTINUE
*
*       Initialize the regularized solution.
      CALL KSINIT
*
*       Check optional binary analysis after merger or multiple collision.
C      IF (KZ(4).GT.0.AND.IPHASE.GT.3) THEN
C          CALL EVOLVE(NPAIRS,-1)
C      END IF
*
*       Check updating of global index for chain c.m.
      IF (NCH.GT.0) THEN
          CALL CHFIND
      END IF
*
      RETURN
*
      END
