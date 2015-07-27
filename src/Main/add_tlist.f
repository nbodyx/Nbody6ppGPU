      subroutine add_tlist(J,STEP,DTK)
*
*
*     Add particle in NXTLST
*
      include 'params.h'
      include 'tlist.h'
      REAL*8 STEP(NMAX),DTK(64)
      INTEGER L,J,K

*     Get step level
      K = k_step(STEP(J),DTK)
*     Check ghost
      IF(K.EQ.-1) THEN
         DO L = 2, NLSTDELAY(1)+1
            IF (NLSTDELAY(L).EQ.J) then
               print*,'Error!: ghost partile ',J,' already exist in ',
     &              'delay list'
               call abort()
            END IF
         END DO
         NGHOSTS = NGHOSTS + 1
         NXTLST(NXTLIMIT+NGHOSTS) = J
*     --07/16/14 20:21-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         print*,'ADD GHOST J',J,'STEP',STEP(J),' NGHOSTS', NGHOSTS
*     --07/16/14 20:21-lwang-end----------------------------------------*
         RETURN
      END IF
*     Safe check
      IF(K.LT.1) then
         write(6,*) 'Error: New particle step level too large K=',K
         call flush(6)
         call abort()
      END IF
      IF(K.GT.64) then
         write(6,*) 'Error: New particle step level too small K=',K
         call flush(6)
         call abort()
      END IF
*     Update smallest step level indicator
      IF(NDTMAX.LT.K) NDTMAX = K
      IF(NDTMIN.GT.K) NDTMIN = K
*     Safe check
      IF(NXTLIMIT.GE.NMAX) THEN
         write(6,*) 'Error: NXTLST maximum reached!'
         call flush(6)
         call abort()
      END IF
*     Increase nxtlst ending point
      NXTLIMIT = NXTLIMIT + 1
      IF(NGHOSTS.GT.0) NXTLST(NXTLIMIT+NGHOSTS) = NXTLST(NXTLIMIT)
*     Increase all NDTK outside NDTMIN by one
      DO L = 1,NDTMIN-1
         NDTK(L) = NDTK(L) + 1
      END DO
      DO L = NDTMIN,K-1
*     Increase step level L position by one
         NDTK(L) = NDTK(L) + 1
*     Shift first index position to end of level L
         NXTLST(NDTK(L)) = NXTLST(NDTK(L+1)+1)
      END DO
*     Increase step level K by one
      NDTK(K) = NDTK(K) + 1
*     Add new particle at the end of step level K
      NXTLST(NDTK(K)) = J

*     --07/08/14 16:48-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      print*,'ADD J',J,'STEP',STEP(J),'NXTLIMIT',NXTLIMIT,'K',K,
c$$$     &     'NDTMAX',NDTMAX
*     --07/08/14 16:48-lwang-end----------------------------------------*
      RETURN

      END

