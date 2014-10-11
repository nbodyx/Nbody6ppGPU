      SUBROUTINE KCPERT_COR(II,I1,FIRR,FD)
*
*
*       Differential force correction on active KS from chain.
*       ------------------------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  FIRR(3),FD(3),DX(3),DV(3),FP(3),FPD(3),FCM(3),FCMD(3)
*
*
*       See whether perturber list contains chain c.m. body #ICH. 
      J = 0
      NNB1 = LIST(1,I1) + 1
      DO 1 L = 2,NNB1
          JJ = LIST(L,I1)
          IF (JJ.EQ.ICH) J = JJ
    1 CONTINUE
*       Exit if chain c.m. not identified.
      IF (J.EQ.0) GO TO 30
*
*       Obtain chain c.m. force on each KS component for diff correction.
      BODYIN = 1.0/BODY(II)
      I = I1
    5 A1 = X(1,J) - X(1,I)
      A2 = X(2,J) - X(2,I)
      A3 = X(3,J) - X(3,I)
      RIJ2 = A1*A1 + A2*A2 + A3*A3
*
      DV(1) = XDOT(1,J) - XDOT(1,I)
      DV(2) = XDOT(2,J) - XDOT(2,I)
      DV(3) = XDOT(3,J) - XDOT(3,I)
      DR2I = 1.0/RIJ2
      DR3I = BODY(J)*DR2I*SQRT(DR2I)
      DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
*       Save c.m. contribution (to be subtracted).
      FCM(1) = A1*DR3I
      FCM(2) = A2*DR3I
      FCM(3) = A3*DR3I
      FCMD(1) = (DV(1) - A1*DRDV)*DR3I
      FCMD(2) = (DV(2) - A2*DRDV)*DR3I
      FCMD(3) = (DV(3) - A3*DRDV)*DR3I
*
      DO 10 K = 1,3
          FP(K) = 0.0
          FPD(K) = 0.0
   10 CONTINUE
*
*       Form contributions from each chain member for each component.
      DO 20 JJ = 1,NCH
          DR2 = 0.0
          DRDV = 0.0
          DO 15 L = 1,3
              DX(L) = XC(L,JJ) - X(L,I)
              DV(L) = UC(L,JJ) - XDOT(L,I)
              DR2 = DR2 + DX(L)**2
              DRDV = DRDV + DX(L)*DV(L)
   15     CONTINUE
          DR2I = 1.0/DR2
          DR3I = BODYC(JJ)*DR2I*SQRT(DR2I)
          DRDV = 3.0*DRDV*DR2I
*
*       Sum force & first derivative.
          DO 18 L = 1,3
              FP(L) = FP(L) + DX(L)*DR3I
              FPD(L) = FPD(L) + (DV(L) - DX(L)*DRDV)*DR3I
   18     CONTINUE
   20 CONTINUE
*
*       Add differential corrections (chain contribution minus c.m. term).
      DO 25 K = 1,3
          FIRR(K) = FIRR(K) + BODY(I)*(FP(K) - FCM(K))*BODYIN
          FD(K) = FD(K) + BODY(I)*(FPD(K) - FCMD(K))*BODYIN
   25 CONTINUE
*
*       Treat the second component in a similar way.
      IF (I.EQ.I1) THEN
          I = I1 + 1
          GO TO 5
      END IF
*
   30 RETURN
*
      END
