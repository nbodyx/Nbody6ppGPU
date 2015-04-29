      SUBROUTINE INSERT(I,LI)
*
*
*       Insert particle index in KS time-step list.
*       -------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  TI,TI1,TI2,TJ
*
*
*       Avoid increasing KBLIST if body #I can be exchanged with next member.
      LI2 = MIN(LI + 2,NNTB)
      I2 = KBLIST(LI2)
      TI = T0(I) + STEP(I)
      TI2 = T0(I2) + STEP(I2)
*
*       Check swapping condition (#I due before #I2 but after #I1 at LI + 1).
      IF (TI.LE.TI2) THEN
          LI1 = MIN(LI+1,NNTB)
          I1 = KBLIST(LI1)
          TI1 = T0(I1) + STEP(I1)
          IF (TI.GT.TI1) THEN
              KBLIST(LI) = KBLIST(LI1)
              KBLIST(LI1) = I
          END IF
*       Reduce pointer so that current location will be selected again.
          LI = LI - 1
          GO TO 50
*       Also swap if body #I is among the two last members and TI > TI2.
      ELSE IF (LI.GT.NNTB - 2) THEN
          KBLIST(LI) = KBLIST(LI2)
          KBLIST(LI2) = I
          LI = LI - 1
          GO TO 50
      END IF
*
*       Estimate the insert index from the remaining time interval.
      FAC = MIN(10.0,STEP(I)/(TBLIST - TIME))
      LSTAR = LI + FLOAT(NNTB - LI)*FAC
*     --12/21/13 22:36-lwang-print--------------------------------------*
***** Note:------------------------------------------------------------**
#ifdef PRINT      
      PRINT*,' I,LSTAR,LI,FAC,FLOAT(NNTB - LI),STEP=',
     *       I,LSTAR,LI,FAC,FLOAT(NNTB - LI),STEP(I)
      call flush(6)
#endif      
*     --12/21/13 22:36-lwang-end----------------------------------------*
*     
      IF(LSTAR.LE.NNTB) THEN
*       Improve insert index by another iteration (check LI < LSTAR <= NNTB).
         J = KBLIST(LSTAR)
         TJ = T0(J) + STEP(J)
*       Avoid division by zero (unperturbed steps may also be quantized).
         IF (TBLIST.NE.TJ) THEN
            FAC = (TI - TJ)/(TBLIST - TJ)
            LSTAR = LSTAR + FLOAT(NNTB - LSTAR)*FAC
            LSTAR = MAX(LSTAR,LI+1)
            LSTAR = MIN(LSTAR,NNTB)
            J = KBLIST(LSTAR)
         END IF
      ELSE
         LSTAR = MAX(NNTB-1,LI+1)
         J = KBLIST(LSTAR)
*
      END IF
*
*      --12/21/13 22:37-lwang-print--------------------------------------*
***** Note:------------------------------------------------------------**
#ifdef PRINT
      J1 = KBLIST(LI)
      JNNTB = KBLIST(NNTB)
*     PRINT*,' KBL ',J1,JNNTB,T0(J1)+STEP(J1),T0(JNNTB)+STEP(JNNTB)
#endif      
*     --12/21/13 22:37-lwang-end----------------------------------------*
*       Determine correct index by comparing neighbouring KBLIST members.
      IF (TI.GE.T0(J) + STEP(J)) THEN
          L1 = LSTAR + 1
          LSTAR = L1
          DO 10 L = L1,NNTB
              J = KBLIST(L)
              TJ = T0(J) + STEP(J)
*       Advance index until TI < TJ or last member.
              IF (TI.LT.TJ) GO TO 30 
              LSTAR = L + 1
   10     CONTINUE
      ELSE
   20     LSTAR = LSTAR - 1
*     --12/30/13 20:23-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          IF(LSTAR.EQ.0) THEN
c$$$             do k=1,NNTB
c$$$             print*,'DEBUG',k,KBLIST(k),T0(KBLIST(K))+STEP(KBLIST(K))
c$$$             end do
c$$$             print*,'DEBUG, I',I,'TI',TI,'LI',LI,T0(LI)+STEP(LI)
c$$$             call flush(6)
c$$$          end if
*     --12/30/13 20:23-lwang-end----------------------------------------*
          J = KBLIST(LSTAR)
          IF (J.EQ.I) THEN
              LI = LI - 1
              GO TO 50
          END IF
          IF (TI.LT.T0(J) + STEP(J).AND.LSTAR.GT.1) GO TO 20
          LSTAR = LSTAR + 1
*       Ensure that current index LI is not chosen due to round-off.
          LSTAR = MAX(LSTAR,LI+1)
      END IF 
*
*       Create free location at LSTAR by moving all subsequent members down.
   30 DO 40 L = NNTB,LSTAR,-1
          KBLIST(L+1) = KBLIST(L)
   40 CONTINUE 
*
*       Insert body #I and update memberships.
      KBLIST(LSTAR) = I
      NNTB = NNTB + 1
*
   50 RETURN
*
      END
