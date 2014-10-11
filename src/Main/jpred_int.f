      SUBROUTINE JPRED_INT(I,TTIME)
*
*
*       Prediction of single particle, used in intgrt.F
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
      REAL*8 TPRED,TTIME
      LOGICAL ipredall
*
      IF(I.LT.ifirst.or.ipredall) return
*      IF (TPRED(I).ne.TTIME) THEN
*     Adopt low order prediction for standard single particle.
         S = TTIME - T0(I)
         S1 = 1.5*S
         S2 = 2.0*S
*     Note X may become corrected value X0 for zero interval.
         X(1,I) = ((FDOT(1,I)*S + F(1,I))*S + X0DOT(1,I))*S + X0(1,I)
         X(2,I) = ((FDOT(2,I)*S + F(2,I))*S + X0DOT(2,I))*S + X0(2,I)
         X(3,I) = ((FDOT(3,I)*S + F(3,I))*S + X0DOT(3,I))*S + X0(3,I)
         XDOT(1,I) = (FDOT(1,I)*S1 + F(1,I))*S2 + X0DOT(1,I)
         XDOT(2,I) = (FDOT(2,I)*S1 + F(2,I))*S2 + X0DOT(2,I)
         XDOT(3,I) = (FDOT(3,I)*S1 + F(3,I))*S2 + X0DOT(3,I)
*         NBPRED = NBPRED + 1
*         TPRED(I) = TTIME
*      END IF
      IF(I.GT.N) THEN
         JPAIR = I - N
         IF(TRES(JPAIR).NE.TTIME.AND.LIST(1,2*JPAIR-1).GT.0) THEN
            ZZ = 1.0
*     Distinguish between low and high-order prediction of U & UDOT.
            IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
            CALL KSRES2(JPAIR,J1,J2,ZZ,TTIME)
            TRES(JPAIR) = TTIME
         END IF
      END IF
*
      RETURN
*
      END
