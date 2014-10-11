      subroutine replace_tlist(I,J,STEP,DTK)
*
*
*     Replace index I by J in NXTLST
*     --------------------------------------------------------
*
      include 'params.h'
      include 'tlist.h'
      REAL*8 STEP(NMAX),DTK(64)
      INTEGER I,J,LI,LL,L

*     --07/15/14 17:16-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      print*,'REPLACE',I,J
c$$$      call flush(6)
*     --07/15/14 17:16-lwang-end----------------------------------------*
*     Get level
      LI = k_step(STEP(I),DTK)
*
*     Check ghost
      IF(LI.EQ.-1) THEN
         IF(NGHOSTS.GT.0) THEN
            L = NXTLIMIT + 1
 2          IF(NXTLST(L).EQ.I) THEN
               NXTLST(L) = J
               RETURN
            ELSE
               L = L + 1
               IF(L-NXTLIMIT.LE.NGHOSTS) GO TO 2
            END IF
         END IF
         write(6,*) 'Error: No ghost particle J',J,'STEP',STEP(J)
         call flush(6)
         call abort()
      ELSE
         LL = NDTK(LI+1)+1
 1       IF(NXTLST(LL).EQ.I) THEN
            NXTLST(LL) = J
         ELSE
            IF(LL.LT.NDTK(LI)) THEN
               LL = LL + 1
               GO TO 1
            ELSE
               IF(NLSTDELAY(1).GT.0) THEN
                  DO L = 2, NLSTDELAY(1)+1
                     IF(NLSTDELAY(L).EQ.I) THEN
                        NLSTDELAY(L) = J
                        RETURN
                     END IF
                  END DO
               END IF
               write(6,*) 'Error: Index ',I,' not found in step level ',
     &              LI,'!'
               call flush(6)
               call abort()
            END IF
         END IF
      END IF

      RETURN

      END
