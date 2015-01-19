      SUBROUTINE SHRINK(TMIN)
*
*
*       Shrinking of large time-steps.
*       ------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (DTM = 0.03125D0)
      REAL*8  RIJ(3),VIJ(3)
*
*
*       Determine impact to each high-velocity particle for large STEPR.
      NNB = LISTV(1)
      DO 20 L = 1,NNB
          I = LISTV(L+1)
C!$omp parallel do private(J,NNB1,K,RV,VV,RIJ,DT,IT,R2,FJ2,FI2)
          DO 10 J = IFIRST,NTOT
*       Form relative quantities for large regular time-steps.
              IF (STEPR(J).LT.DTM.OR.J.EQ.I) GO TO 10
*       Check whether body #I is already a neighbour (accept LIST > I).
              NNB1 = LIST(1,J) + 1
              DO 1 K = 2,NNB1
                  IF (LIST(K,J).GT.I) GO TO 2
                  IF (LIST(K,J).EQ.I) GO TO 10
    1         CONTINUE
    2         RV = 0.0
              VV = 0.0
              DO 3 K = 1,3
                  RIJ(K) = X(K,I) - X(K,J)
                  VIJ(K) = XDOT(K,I) - XDOT(K,J)
                  RV = RV + RIJ(K)*VIJ(K)
                  VV = VV + VIJ(K)**2
    3         CONTINUE
*
*       Skip treatment for receding particles.
              IF (RV.GE.0.0D0) GO TO 10
*
*       Evaluate time of minimum approach and truncate to next time.
              DT = -RV/VV
              IT = 0
    5         DT = MIN(T0R(J) + STEPR(J) - TIME,DT)
*
*       Obtain minimum impact parameter for straight-line orbit.
              R2 = 0.0
              FJ2 = 0.0
              DO 6 K = 1,3
                  R2 = R2 + (RIJ(K) + VIJ(K)*DT)**2
                  FJ2 = FJ2 + F(K,J)**2
    6         CONTINUE
*
*       Compare force at minimum distance with total force on body #J (> 0).
              FI2 = (BODY(I)/R2)**2
              IF (FI2.LT.0.04*FJ2.OR.BODY(J).EQ.0.0D0) GO TO 10
*
*       See whether the regular time-step can be shortened.
              IF (T0R(J) + 0.5*STEPR(J).GT.TMIN.AND.IT.LT.5) THEN
                  if(rank.eq.0)
     &            WRITE (14,8)  J, SQRT(R2), SQRT(FI2/FJ2), DT, STEPR(J)
    8             FORMAT (' SHRINK    I RN FI/FJ DT STEPR ',
     &                                I6,F7.3,F6.2,1P,2E10.2)
                  CALL FLUSH(14)
                  STEPR(J) = 0.5*STEPR(J)
                  IF (STEP(J).GT.STEPR(J)) THEN
*     Update NXTLST also
C                     call shift_tlist(J,STEP(J),1,DTK)
C!$omp critical           
                      call delay_remove_tlist(J,STEP,DTK)
C!$omp end critical           
                      STEP(J) = 0.5*STEP(J)
C!$omp critical           
                      call delay_store_tlist(J)
C!$omp end critical           
                      TIMENW(J) = T0(J) + STEP(J)
*                      TMIN = MIN(TIMENW(J),TMIN)
                  END IF
                  IT = IT + 1
                  GO TO 5
              END IF
   10     CONTINUE
C!$omp end parallel do
   20 CONTINUE

*
      RETURN
*
      END
