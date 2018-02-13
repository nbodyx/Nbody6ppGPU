      subroutine star_trace(NXTLEN,NXTLST)
*
*
*     Trace one star orbit and output every irregular step
*-------------------------------------------------------------------------------
*     
*     NAMEST: start name
*     NAMEED: end name
      include 'common6.h'
      LOGICAL istart
      DATA istart /.true./
      SAVE istart
      CHARACTER(LEN=20) :: TCHAR
      CHARACTER(LEN=64) :: FNAME
      
      IF (istart) then
         DO K=NAMEST,NAMEED
            KUINT = K - NAMEST + 300
            call string_left(TCHAR,REAL(K),1.0)
            FNAME='star.'//trim(TCHAR)
            FNAME=trim(FNAME)
            OPEN (UNIT=KUNIT,STATUS='UNKNOWN',FORM='FORMATTED',
     &           FILE=FNAME)
            write(6,*) 'Open trace file: ',FNAME
         END DO
         istart = .FALSE.
      END IF 
      
      DO J = 1, NXTLEN
         I = NXTLST(J)
         IF(I.GT.N) THEN
            IPAIR = I - N
            I1 = 2*IPAIR + 1
            I2 = I1 + 1
            IRES2 = 0
            IF (NAME(I1).GE.NAMEST.AND.NAME(I1).LE.NAMEED) THEN
               call ksres2(IPAIR,I1,I2,0.0,TIME)
               IRES2 = 1
               KUNIT = NAME(I1) - NAMEST + 300
               write(KUNIT,*) BODY(I1)*ZMBAR, X(1:3,I1)*RBAR, 
     &              XDOT(1:3,I1)*VSTAR, I1
            END IF
            IF (NAME(I2).GE.NAMEST.AND.NAME(I2).LE.NAMEED) THEN
               IF(IRES2.EQ.0) call ksres2(IPAIR,I1,I2,0.0,TIME)
               KUNIT = NAME(I1) - NAMEST + 300
               write(KUNIT,*) BODY(I2)*ZMBAR, X(1:3,I2)*RBAR, 
     &              XDOT(1:3,I2)*VSTAR, I2
            END IF 
         ELSE
            IF (NAME(I).GE.NAMEST.AND.NAME(I).LE.NAMEED) THEN
               KUNIT = NAME(I) - NAMEST + 300
               write(KUNIT,*) BODY(I)*ZMBAR, X(1:3,I)*RBAR, 
     &              XDOT(1:3,I)*VSTAR, I
            END IF
         END IF 
      END DO
      
      RETURN

      END

