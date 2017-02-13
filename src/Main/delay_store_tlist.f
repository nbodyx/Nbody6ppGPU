      subroutine delay_store_tlist(J)
*
*
*     store particle in NLSTDELAY
*
      include 'params.h'
      include 'tlist.h'
      INTEGER L
      
      IF(NLSTDELAY(1).LT.NDELAYMAX-1) THEN
         DO L = 2, NLSTDELAY(1)+1
            IF (NLSTDELAY(L).EQ.J) RETURN
         END DO
         DO L = NXTLIMIT+1,NXTLIMIT+NGHOSTS
            IF (NXTLST(L).EQ.J) THEN
               print*,'Error!: J particle',J,' already in ghost list!'
               call abort()
            END IF
         END DO
*     --07/16/14 11:53-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*         print*,'DELAY STORE J',J,'NLSTDELAY NUM',NLSTDELAY(1)
*     &        NLSTDELAY(1:NLSTDELAY(1)+1)
*     --07/16/14 11:53-lwang-end----------------------------------------*
         NLSTDELAY(1) = NLSTDELAY(1) + 1
         NLSTDELAY(NLSTDELAY(1)+1) = J
      ELSE
         write(6,*) 'Error: NLSTDELAY overflow!'
         call abort()
      END IF

      RETURN

      END

