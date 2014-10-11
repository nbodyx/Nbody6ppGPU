      subroutine delay_remove_tlist(J,STEP,DTK)
*
*
*     remove particle in NLSTDELAY, if not exist, remove from NXTLST
*
      include 'params.h'
      include 'tlist.h'
      REAL*8 STEP(NMAX),DTK(64)
      INTEGER L

      IF(NLSTDELAY(1).GT.0) THEN
         L = 2
 1       IF (NLSTDELAY(L).EQ.J) THEN
*     --07/16/14 11:56-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            print*,'DELAY REMOVE J',J
*     --07/16/14 11:56-lwang-end----------------------------------------*
            NLSTDELAY(L) = NLSTDELAY(NLSTDELAY(1)+1)
            NLSTDELAY(1) = NLSTDELAY(1) - 1
            RETURN
         ELSE      
            L = L + 1
            IF(L.LE.NLSTDELAY(1)+1) GO TO 1
         END IF
      END IF
*     If not find in NLSTDELAY, try to remove from NXTLST
      call remove_tlist(J,STEP,DTK)

      RETURN

      END

