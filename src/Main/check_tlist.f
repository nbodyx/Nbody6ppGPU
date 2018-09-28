      subroutine check_tlist()
*
*
*     Remove particle in NXTLST for escape
*
*     B_FLAG: binary remove flag

      include 'common6.h'
      include 'tlist.h'

      DO L = NDTMIN,NDTMAX-1
         I1 = NDTK(L)
         I2 = NDTK(L+1)+1
         DTCHECK = DTK(L)
         IF(I1.GT.I2) THEN
         DO K = I1,I2,-1
            IF(STEP(NXTLST(K)).NE.DTCHECK) THEN
               write(6,*) 'ERROR: STEP(',NXTLST(K),
     &              ') != DTK(',L,')',STEP(NXTLST(K)),DTCHECK,
     &              'NXTLST:K ',K
               call flush(6)
               call abort
            END IF 
         END DO
         END IF
      END DO

C      TMIN = T0(NXTLST(1)) + DTK(NDTMAX)
C 
C      
C      DO L = 1, NXTLEN
C         J=NXTLST(L)
C         IF(T0(J)+STEP(J).NE.TMIN) then
C            print*,'ERROR: T0+STEP .NE. TMIN for particle: ',
C     &           'I',J,'N',NAME(J),'TMIN',TMIN,'T0+STEP',T0(J)+STEP(J),
C     &           'TIME',TIME,'TPREV',TPREV,'T0',T0(J),'STEP',STEP(J),
C     &           'T0(1)',T0(NXTLST(1)),'NDTMAX',NDTMAX,DTK(NDTMAX),
C     &           'NXTLEVEL',NXTLEVEL,'NXTLEN',NXTLEN,'I_NXTLST',L,
C     &           'NXTLIMIT',NXTLIMIT,'T0(2)',T0(NXTLST(2))
C            print*,'NDTK',(K,NDTK(K),DTK(K),K=1,NDTMAX)
C            call flush(6)
C            call abort()
C         END IF
C         IF(BODY(J).EQ.0.D0) THEN
C            print*,'ERROR: Ghost particle in nxtlst: L',L,
C     &           'I',J,'N',NAME(J),'T0+STEP',T0(J)+STEP(J),
C     &           'TIME',TIME,'T0',T0(J),'STEP',STEP(J),
C     &           'T0(1)',T0(NXTLST(1)),'NDTMAX',NDTMAX,DTK(NDTMAX),
C     &           'NXTLIMIT',NXTLIMIT,'M0',BODY0(J),'X',X(1,J),
C     &           'XD',XDOT(1,J),'K',KSTAR(J),'TEV',TEV(J)
C            print*,'NGHOST',NXTLST(NXTLIMIT+1:NXTLIMIT+NGHOSTS)
C            call flush(6)
C            call abort()
C         END IF
C      END DO
C      DO L = NXTLEN+1,NXTLIMIT
C         J = NXTLST(L)
C         IF(T0(J)+STEP(J).LE.TMIN) THEN
C            print*,'ERROR: Lost integrating particle: ',
C     &           'I',J,'N',NAME(J),'TMIN',TMIN,'T0+STEP',T0(J)+STEP(J),
C     &           'T0',T0(J),'STEP',STEP(J),'NDTMAX',NDTMAX,DTK(NDTMAX),
C     &           'NXTLEVEL',NXTLEVEL,'NXTLEN',NXTLEN,'I_NXTLST',L
C            print*,'NDTK',(K,NDTK(K),DTK(K),K=1,NDTMAX)
C            call flush(6)
C            call abort()
C         END IF
C      END DO

      RETURN

      END
      
