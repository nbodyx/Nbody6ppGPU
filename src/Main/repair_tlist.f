      SUBROUTINE REPAIR_TLIST(IADD,ADDL,IRM,RML,STEP,DTK,IGHOST)
*
*
*     Modify TLIST due to the change of data list
*     --------------------------------------------------------
*
*     IADD: Number of adding particles
*     ADDL: Adding particle index list
*     IRM:  Number of removing particles
*     RML:  Removing particle index list
*     IGHOST: Whether shift all removed particles into ghost list
      include 'params.h'
      include 'tlist.h'
      REAL*8 STEP(NMAX),DTK(64)
      INTEGER ADDL(10),RML(10)
      INTEGER I,J
      LOGICAL IGHOST

*     Adding particle
*     Removing particle
*     Safe check
      IF(IRM.GT.10) THEN
         write(6,*) 'Error: Removing particles exceed maximum 10!'
         call flush(6)
         call abort()
      END IF
*     Loop all removing particles
      DO I = 1, IRM
         J = RML(I)
         IF(IGHOST) THEN
            IF(STEP(J).LE.DTK(1)) THEN
               call delay_remove_tlist(J,STEP,DTK)
               STEP(J) = 2*DTK(1)
               call add_tlist(J,STEP,DTK)
            END IF
         ELSE
            call delay_remove_tlist(J,STEP,DTK)
         END IF
      END DO

*     Safe check
      IF(IADD.GT.10) THEN
         write(6,*) 'Error: Adding particles exceed maximum 10!'
         call flush(6)
         call abort()
      END IF
*     Loop all adding particles
      DO I = 1, IADD
         J = ADDL(I)
         call add_tlist(J,STEP,DTK)
      END DO
      
      CALL shrink_tlist

      RETURN

      END

