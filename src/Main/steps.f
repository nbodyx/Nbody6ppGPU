      SUBROUTINE STEPS(I1,I2)
*
*
*       Initialization of time-steps & prediction variables.
*       ----------------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  FDUM(3)
*
*
*       Set new steps and initialize prediction variables.
      DO 40 I = I1,I2
*
*       Include regular force in the irregular step (cf. Makino & SJA 1992).
          DO 5 K = 1,3
              FDUM(K) = FI(K,I) + FR(K,I)
    5     CONTINUE
*
*       Determine irregular and regular steps by the general criterion.
          DT = TSTEP(FDUM,D1(1,I),D2(1,I),D3(1,I),ETAI)
          DTR = TSTEP(FR(1,I),D1R(1,I),D2R(1,I),D3R(1,I),ETAR)
*
*       Reduce irregular step for case of triple, quad, chain or merger.
          IF (IPHASE.GE.4) DT = 0.5*DT
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
          IF (TIME.LE.0.0D0) THEN
              STEP(I) = DTN
              STEPR(I) = DTRN
          ELSE 
*       Reduce steps by factor 2 until commensurate with current time.
              STEP(I) = DTN
              STEPR(I) = DTRN
              ITER = 0
   10         IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
                  STEP(I) = 0.5D0*STEP(I)
                  ITER = ITER + 1
                  IF (ITER.LT.16.OR.STEP(I).GT.DTK(40)) GO TO 10
                  STEP(I) = DTK(40)
                  if(rank.eq.0)
     &            WRITE (6,15) I, ITER, TIME/STEP(I), DT, STEP(I)
   15             FORMAT (' WARNING!   I ITER T/STEP DT STEP ',
     &                                 I5,I4,E16.4,1P,2E9.1)
              END IF
              ITER = 0
   18         IF (DMOD(TIME,STEPR(I)).NE.0.0D0) THEN
                  STEPR(I) = 0.5D0*STEPR(I)
                  ITER = ITER + 1
                  IF (ITER.LT.16.OR.STEPR(I).GT.DTK(40)) GO TO 18
                  STEPR(I) = DTK(40)
                  if(rank.eq.0)
     &            WRITE (6,20) I, ITER, TIME/STEPR(I), DTR, STEPR(I)
   20             FORMAT (' WARNING!   I ITER T/STEPR DTR STEPR ',
     &                                 I5,I4,E16.4,1P,2E9.1)
              END IF
          END IF
*
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
*
      RETURN
*
      END
