      subroutine star_trace(NXTLEN,NXTLST)
*
*
*     Trace one star orbit and output every irregular step
*-------------------------------------------------------------------------------
*     
*     NAMEST: start name
*     NAMEED: end name
      include 'common6.h'
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
      INTEGER NXTLST(NMAX),NXTLEN
      LOGICAL istart
      DATA istart /.true./
      SAVE istart
      CHARACTER(LEN=20) :: FORM
      CHARACTER(LEN=64) :: FNAME
*      REAL*8 XREL(3),VREL(3)
      
      IF (istart) then
         DO K=NAMEST,NAMEED
            KUNIT = K - NAMEST + 300
            NNN = NAMEST
            ip = 0
 3          IF(NNN.gt.1) then
               NNN = NNN/10
               ip = ip + 1
               go to 3
            END IF
            IF(IP.gt.10) then
               write(FORM,'(A2,I2,A1)') '(I',ip,')'
            ELSE
               write(FORM,'(A2,I1,A1)') '(I',ip,')'
            END IF
               
            write(FNAME,FORM) K
            FNAME='star.'//trim(FNAME)
            write(6,*) 'Open trace file: ',FNAME
            call flush(6)
            OPEN (UNIT=KUNIT,STATUS='UNKNOWN',FORM='FORMATTED',
     &           FILE=FNAME)
         END DO
         istart = .FALSE.
      END IF 
      
      DO J = 1, NXTLEN
         I = NXTLST(J)
         IF(I.GT.N) THEN
            IPAIR = I - N
            I1 = 2*IPAIR - 1
            I2 = I1 + 1
            IRES2 = 0
            IF (NAME(I1).GE.NAMEST.AND.NAME(I1).LE.NAMEED) THEN
*               call ksres_op(IPAIR,I1,I2,Q,RDOT,1)
               call ksres2(IPAIR,I1,I2,0.0,TIME)
               IRES2 = 1
               KUNIT = NAME(I1) - NAMEST + 300
               write(KUNIT,*) TIME*TSTAR,BODY(I1)*ZMBAR, X(1:3,I1)*RBAR, 
     &              XDOT(1:3,I1)*VSTAR, I1, LIST(1,I1),GAMMA(IPAIR),
     &              H(IPAIR)
               call flush(KUNIT)
            END IF
            IF (NAME(I2).GE.NAMEST.AND.NAME(I2).LE.NAMEED) THEN
               IF(IRES2.EQ.0) call ksres2(IPAIR,I1,I2,0.0,TIME)
               KUNIT = NAME(I2) - NAMEST + 300
               write(KUNIT,*) TIME*TSTAR,BODY(I2)*ZMBAR, X(1:3,I2)*RBAR, 
     &              XDOT(1:3,I2)*VSTAR, I2, LIST(1,I1),GAMMA(IPAIR),
     &              H(IPAIR)
               call flush(KUNIT)
            END IF 
         ELSE
            IF (NAME(I).GE.NAMEST.AND.NAME(I).LE.NAMEED) THEN
               KUNIT = NAME(I) - NAMEST + 300
               FVALUE = (F(1,I)*F(1,I) + F(2,I)*F(2,I) +
     &              F(3,I)*F(3,I))**0.5
               write(KUNIT,*) TIME*TSTAR,BODY(I)*ZMBAR, X(1:3,I)*RBAR, 
     &              XDOT(1:3,I)*VSTAR, I, LIST(1,I), FVALUE, -PHIDBL(I)
               call flush(KUNIT)
            END IF
         END IF
      END DO
      
      RETURN

      END

