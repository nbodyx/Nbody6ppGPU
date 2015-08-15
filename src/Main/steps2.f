      SUBROUTINE STEPS2(I1,I2)
*
*
*       Initialization of time-steps by Xdot, F and Fdot & prediction variables.
*       ----------------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  FDUM(3)
*
*
*       Set new steps and initialize prediction variables.
!$omp parallel do private(I,K,DT,DTR,DTN,DTRN,FDUM)
      DO 40 I = I1,I2
*
*       Include regular force in the irregular step (cf. Makino & SJA 1992).
          DO K = 1,3
             FDUM(K) = FI(K,I) + FR(K,I)
          END DO
*
*       Determine irregular and regular steps by the general criterion.
          DT = TSTEP2(XDOT(1,I),FDUM,D1(1,I),ETAI)
          DTR = TSTEP2(XDOT(1,I),FR(1,I),D1R(1,I),ETAR)
*       Avoid the case when velocity is zero, time step is zero
          DT = MAX(DT,DTK(64))
          DTR = MAX(DT,DTR)
*
*       Initialize the times and obtain discrete steps (block-step version).
          T0(I) = TIME
          T0R(I) = TIME
*
*     Convert predicted step to nearest block time-step (truncated down).
          CALL STEPK(DT,DTN)
          CALL STEPK(DTR,DTRN)
*     --11/22/13 20:26-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         IF(I.EQ.12534.and.rank.eq.0) then
c$$$             print*,'I',I,'T',time,'TB',tblock,'DT',dt,'DTN',dtn
c$$$          end if
*     --11/22/13 20:26-lwang-end----------------------------------------*
          STEP(I) = DTN
          STEPR(I) = DTRN

*       Reduce regular step if STEPR > SMAX
          STEPR(I) = MIN(STEPR(I),SMAX)
*          
*       Reduce irregular step if STEPR < STEP.
 25       IF (STEPR(I).LT.STEP(I)) THEN
             STEP(I) = 0.5D0*STEP(I)
             GO TO 25
          END IF
*
*       Initialize or update array for new block times.
          TIMENW(I) = T0(I) + STEP(I)
*
*       Set prediction variables (X0DOT set by START, KSREG or KSTERM).
          DO 30 K = 1,3
              X0(K,I) = X(K,I)
              F(K,I) = 0.5D0*F(K,I)
              FDOT(K,I) = ONE6*FDOT(K,I)
   30     CONTINUE
   40 CONTINUE
!$omp end parallel do
      
*
      RETURN
*
      END
