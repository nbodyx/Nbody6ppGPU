      subroutine delay_add_tlist(T0,STEP,DTK)
*
*
*     Add particle in NXTLST due to T0+STEP
*
      include 'params.h'
      include 'tlist.h'
      REAL*8 STEP(NMAX),DTK(64),T0(NMAX)
      INTEGER L,J
      
      IF(NLSTDELAY(1).EQ.0) RETURN
      L = 2
 1    J = NLSTDELAY(L)
      IF(T0(J)+STEP(J).LE.T0(NXTLST(1))+DTK(NDTMAX)) THEN
*     --07/15/14 23:23-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         print*,'DELAY ADD',J,'STEP',STEP(J),'T0',T0(J),'TIME',TIME,
c$$$     &        'TMIN',T0(NXTLST(1))+ DTK(NDTMAX),
c$$$     &        'CURRENT LIST',NLSTDELAY(1:NLSTDELAY(1)+1)
*     --07/15/14 23:23-lwang-end----------------------------------------*
         call add_tlist(J,STEP,DTK)
         IF(L.LT.NLSTDELAY(1)+1) THEN
            NLSTDELAY(L) = NLSTDELAY(NLSTDELAY(1)+1)
            NLSTDELAY(1) = NLSTDELAY(1) - 1
            GO TO 1
         ELSE
            NLSTDELAY(1) = NLSTDELAY(1) - 1
         END IF
      ELSE
         L = L + 1
         IF(L.LE.NLSTDELAY(1)+1) GO TO 1
      END IF

      RETURN

      END
