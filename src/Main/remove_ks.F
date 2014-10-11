      SUBROUTINE REMOVE_KS(IPAIR)
*
*
*     Remove KS pair, components or c.m. 
*       -----------------
*
      INCLUDE 'common6.h'
      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
      REAL*8 TPRED,DSAVE
      LOGICAL ipredall
*
      ICM = IPAIR + N
*     Move last cm data to I (escaper or old c.m. & KS comps).
*     Notice NTOT is already reduced by 1, the last cm should be NTOT+1
*     If I is last cm, escape movement

      IF (ICM.GT.NTOT) THEN
*     Remove ICM from NXTLST
         IF(TIME+TOFF.NE.0.D0) then
            call delay_remove_tlist(ICM,STEP,DTK)
            call shrink_tlist
         END IF
         RETURN
      END IF
*
      J = NTOT + 1
*     Exchange ICM and Last c.m. in NXTLST
      IF(TIME+TOFF.NE.0.D0) call exchange_tlist(ICM,J,STEP,DTK)

*     save step for remove icm from NXTLST
      DSAVE = STEP(ICM)

      X(1:3,ICM) = X(1:3,J)
      X0(1:3,ICM) = X0(1:3,J)
      X0DOT(1:3,ICM) = X0DOT(1:3,J)
      XDOT(1:3,ICM) = XDOT(1:3,J)
      F(1:3,ICM) = F(1:3,J)
      FDOT(1:3,ICM) = FDOT(1:3,J)
      FI(1:3,ICM) = FI(1:3,J)
      FIDOT(1:3,ICM) = FIDOT(1:3,J)
      D0(1:3,ICM) = D0(1:3,J)
      D1(1:3,ICM) = D1(1:3,J)
      D2(1:3,ICM) = D2(1:3,J)
      D3(1:3,ICM) = D3(1:3,J)
      FR(1:3,ICM) = FR(1:3,J)
      FRDOT(1:3,ICM) = FRDOT(1:3,J)
      D0R(1:3,ICM) = D0R(1:3,J)
      D1R(1:3,ICM) = D1R(1:3,J)
      D2R(1:3,ICM) = D2R(1:3,J)
      D3R(1:3,ICM) = D3R(1:3,J)
*
      TPRED(ICM) = TPRED(J)
      BODY(ICM) = BODY(J)
      RS(ICM) = RS(J)
      RADIUS(ICM) = RADIUS(J)
      TEV(ICM) = TEV(J)
      TEV0(ICM) = TEV0(J)
      BODY0(ICM) = BODY0(J)
      EPOCH(ICM) = EPOCH(J)
      SPIN(ICM) = SPIN(J)
      ZLMSTY(ICM) = ZLMSTY(J)
      KSTAR(ICM) = KSTAR(J)
      NAME(ICM) = NAME(J)
      STEP(ICM) = STEP(J)
      STEPR(ICM) = STEPR(J)
      T0(ICM) = T0(J)
      T0R(ICM) = T0R(J)
      TIMENW(ICM) = T0(J) + STEP(J)
*
*     Transfer unmodified neighbour list.
      NNB = LIST(1,J) + 1
*     Include flag of 2nd component (note new IFIRST if escaping pair).
      IF (ICM.LE.IFIRST.AND.NNB.EQ.1) NNB = 2
      LIST(1:NNB,ICM) = LIST(1:NNB,J)

*     Remove terminated c.m. in NXTLST
      STEP(J) = DSAVE
      IF(TIME+TOFF.NE.0.D0) then
         call delay_remove_tlist(J,STEP,DTK)
         call shrink_tlist
      end if

*
#ifdef SIMD
*     Update Particle data in AVX/SSE library
      IF(TIME.NE.0.0D0) THEN
         CALL IRR_SIMD_SET_JP(ICM,X0(1,ICM),X0DOT(1,ICM),F(1,ICM),
     &        FDOT(1,ICM),BODY(ICM),T0(ICM))
         CALL IRR_SIMD_SET_LIST(ICM,LIST(1,ICM))
      END IF
#endif
*      
*     Move KS data of last pair to I (Notice here npairs is already reduced by 1, thus the last pair index is npairs+1)
      LPAIR = NPAIRS + 1
      U(1:4,IPAIR) = U(1:4,LPAIR)
      U0(1:4,IPAIR) = U0(1:4,LPAIR)
      UDOT(1:4,IPAIR) = UDOT(1:4,LPAIR)
      FU(1:4,IPAIR) = FU(1:4,LPAIR)
      FUDOT(1:4,IPAIR) = FUDOT(1:4,LPAIR)
      FUDOT2(1:4,IPAIR) = FUDOT2(1:4,LPAIR)
      FUDOT3(1:4,IPAIR) = FUDOT3(1:4,LPAIR)
      SF(1:4,IPAIR) = SF(1:4,LPAIR)
      FP0(1:4,IPAIR) = FP0(1:4,LPAIR)
      FD0(1:4,IPAIR) = FD0(1:4,LPAIR)
*
      R(IPAIR) = R(LPAIR)
      R0(IPAIR) = R0(LPAIR)
      DTAU(IPAIR) = DTAU(LPAIR)
      TDOT2(IPAIR) = TDOT2(LPAIR)
      TDOT3(IPAIR) = TDOT3(LPAIR)
      GAMMA(IPAIR) = GAMMA(LPAIR)
      H(IPAIR) = H(LPAIR)
      HDOT(IPAIR) = HDOT(LPAIR)
      HDOT2(IPAIR) = HDOT2(LPAIR)
      HDOT3(IPAIR) = HDOT3(LPAIR)
      HDOT4(IPAIR) = HDOT4(LPAIR)
      KSLOW(IPAIR) = KSLOW(LPAIR)
      SF(5,IPAIR) = SF(5,LPAIR)
      SF(6,IPAIR) = SF(6,LPAIR)
      SF(7,IPAIR) = SF(7,LPAIR)
      H0(IPAIR) = H0(LPAIR)
*
*     Move TRES of last pair to ipair
      TRES(IPAIR) = TRES(LPAIR)

      RETURN
*
      END
