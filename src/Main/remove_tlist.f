      subroutine remove_tlist(J,STEP,DTK)
*
*
*     Remove particle in NXTLST
*
      include 'params.h'
      include 'tlist.h'
      REAL*8 STEP(NMAX),DTK(64)
      INTEGER L,J,LL,K

*     Get step level
      K = k_step(STEP(J),DTK)
*     Check ghost
      IF(K.EQ.-1) THEN
         IF(NGHOSTS.GT.0) THEN
            L = NXTLIMIT + 1
 2          IF(NXTLST(L).EQ.J) THEN
               IF(L-NXTLIMIT.LT.NGHOSTS) 
     &              NXTLST(L) =  NXTLST(NXTLIMIT+NGHOSTS)
               NGHOSTS = NGHOSTS - 1
*     --07/16/14 20:23-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$               print*,'REMOVE GHOST J',J,'STEP',STEP(J),'NGHOSTS',
c$$$     &              NGHOSTS,'LIST',NXTLST(NXTLIMIT+1:NXTLIMIT+NGHOSTS)
*     --07/16/14 20:23-lwang-end----------------------------------------*
               RETURN
            ELSE
               L = L + 1
               IF(L-NXTLIMIT.LE.NGHOSTS) GO TO 2
            END IF
         END IF
         write(6,*) 'Error: No ghost particle J',J,'STEP',STEP(J),
     &        'NGHOSTS',NGHOSTS,'LIST',
     &        NXTLST(NXTLIMIT+1:NXTLIMIT+NGHOSTS)
         call flush(6)
         call abort()
      END IF

*     decrease nxtlst ending point
      IF(NXTLIMIT.LE.1) THEN
         write(6,*) 'Error: No particle in NXTLST!'
         call flush(6)
         call abort()
      END IF
      NXTLIMIT = NXTLIMIT - 1
*     Reduce all NDTK outside NDTMIN by one
      DO L = 1,NDTMIN-1
         NDTK(L) = NDTK(L) - 1
      END DO
*     Replace removing particle index by the end index of step level K
      LL = NDTK(K+1)+1
 1    IF(NXTLST(LL).EQ.J) THEN
         NXTLST(LL) = NXTLST(NDTK(K))
         NDTK(K) = NDTK(K) - 1
      ELSE
         IF(LL.LT.NDTK(K)) THEN
            LL = LL + 1
            GO TO 1
         ELSE
            write(6,*) 'Error: Index ',J,' not found in step level ',
     &           K,'!'
            call flush(6)
            call abort()
         END IF
      END IF
      DO L = K-1,NDTMIN,-1
*     Shift last index position to beginning of level L
         NXTLST(NDTK(L+1)+1) = NXTLST(NDTK(L))
*     Reduce step level L position by one
         NDTK(L) = NDTK(L) - 1
      END DO
*     Shift ghost particles
      IF(NGHOSTS.GT.0) NXTLST(NXTLIMIT+1) = NXTLST(NXTLIMIT+NGHOSTS+1)

      CALL shrink_tlist

*     --07/08/14 16:49-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      print*,'REMOVE J',J,'NXTLIMIT',NXTLIMIT,'K',K,'L',LL,
c$$$     &     'NDTK',NDTK(K),'NDTK(K+1)',NDTK(K+1)
*     --07/08/14 16:49-lwang-end----------------------------------------*

      RETURN

      END
