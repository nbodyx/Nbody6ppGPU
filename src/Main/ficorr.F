      SUBROUTINE FICORR(I,DM)
*
*
*       Local force corrections due to mass loss.
*       -----------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDe 'omp_lib.h'
      REAL*8  A(6)
*      REAL*8  XBACK(3,LMAX)
*
#ifdef SIMD
*     Reset AVX/SSE particle data
      call IRR_SIMD_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &     BODY(I),T0(I))
#endif
*       Correct force & first derivative for neighbours of body #I.
      call xbpredall
*      call jpred(I,TIME,TIME)
      NNB1 = LIST(1,I) + 1
      DO 20 L = 2,NNB1
          J = LIST(L,I)
          NNB2 = LIST(1,J) + 1
*
*       Ensure that body #I is a neighbour of body #J.
          DO 1 K = 2,NNB2
              IF (LIST(K,J).EQ.I) GO TO 5
    1     CONTINUE
          GO TO 20
*
    5     RIJ2 = 0.0D0
          A7 = 0.0D0
*
*          call jpred(j,time,time)
          DO 10 K = 1,3
              A(K) = X(K,I) - X(K,J)
              A(K+3) = XDOT(K,I) - XDOT(K,J)
              RIJ2 = RIJ2 + A(K)**2
              A7 = A7 + A(K)*A(K+3)
   10     CONTINUE
*
          A7 = 3.0*A7/RIJ2
          A8 = DM/(RIJ2*SQRT(RIJ2))
*
          DO 15 K = 1,3
              A(K+3) = (A(K+3) - A7*A(K))*A8
              F(K,J) = F(K,J) - 0.5*A(K)*A8
              FI(K,J) = FI(K,J) - A(K)*A8
              FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
              D1(K,J) = D1(K,J) - A(K+3)
              FIDOT(K,J) = FIDOT(K,J) - A(K+3)
   15     CONTINUE
#ifdef SIMD
*     Reset AVX/SSE particle data
          call IRR_SIMD_SET_JP(J,X0(1,J),X0DOT(1,J),F(1,J),FDOT(1,J),
     &         BODY(J),T0(J))
#endif
   20 CONTINUE
*
*     Obtain the potential energy due to all particles.
      POTJ = 0.0D0

!$omp parallel do private(J,K,RIJ2)
!$omp& reduction(+:POTJ)
      DO 30 J = IFIRST,NTOT
          IF (J.EQ.I) GO TO 30
          RIJ2 = 0.0D0
          DO 25 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   25     CONTINUE
          POTJ = POTJ + BODY(J)/SQRT(RIJ2)
   30 CONTINUE
!$omp end parallel do
*
*     Update the potential energy loss and include kinetic energy correction.
      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
      EMDOT = EMDOT + DM*(0.5*VI2 - POTJ)
*
*       See whether linearized tidal terms should be included.
      IF (KZ(14).GT.0.AND.KZ(14).LT.3) THEN
          EMDOT = EMDOT - 0.5*DM*(TIDAL(1)*X(1,I)**2 +
     &                            TIDAL(3)*X(3,I)**2)
      END IF
*
*       Check optional Plummer potential.
      IF (KZ(14).EQ.4.OR.KZ(14).EQ.3) THEN
          RI2 = AP2
          DO 50 K = 1,3
              RI2 = RI2 + X(K,I)**2
  50      CONTINUE
          EMDOT = EMDOT - DM*MP/SQRT(RI2)
      END IF
*
      RETURN
*
      END
