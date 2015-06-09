      SUBROUTINE TTCAL
*
*
*       Computation of tidal elements by quadratic interpolation (mode A)
*       --------------------------------------------------------
*** FlorentR - new subroutine

      INCLUDE 'common6.h'
      INCLUDE 'tt.h'
      REAL*8 DT21,DT32,DT31,TTA,TTB,T1,T2,T3, TNOW
      integer nrotjac
      LOGICAL FIRST

      SAVE FIRST
      DATA FIRST /.TRUE./

* INDTT saves the id of the tensor whose timestamp is immediatly after
*   the current time, set in TNOW
      TNOW = TIME+TOFF
      
      IF(FIRST.AND.TNOW.EQ.0.0D0) THEN
        INDTT=2
        TNOWSAVE = TNOW
        FIRST=.FALSE.
      ELSE
* If the time has not been updated, the tensor is up-to-date: do nothing and exit
        IF(TNOW.EQ.TNOWSAVE) RETURN
      ENDIF

      TNOWSAVE = TNOW

      IF(TNOW.GT.TTTIME(NBTT)) THEN
*  When TNOW is after the timestamp of the last tensor:
*  Terminate after optional COMMON save.
        WRITE(6,*) 'ERROR: the tidal tensor is not defined anymore'
        IF (KZ(1).GT.0) CALL MYDUMP(1,1)
        STOP
      ELSE
* update indtt
        DO WHILE(TNOW.GT.TTTIME(INDTT))
          INDTT = INDTT + 1
        END DO
      ENDIF

* Compute the tensor and its derivatives by quadratic interpolation
      DT21 = TTTIME(INDTT) - TTTIME(INDTT-1)
      DT32 = TTTIME(INDTT+1) - TTTIME(INDTT)
      DT31 = TTTIME(INDTT+1) - TTTIME(INDTT-1)
      DO I=1,3
        DO J=1,3
          T1= TTENS(I,J,INDTT-1)
          T2= TTENS(I,J,INDTT)
          T3= TTENS(I,J,INDTT+1)
          TTA = (T2-T1)/DT21
          TTB = ((T3-T1)*DT21-(T2-T1)*DT31) / (DT31*DT32*DT21)

          TTEFF(I,J) = T1 + TTA * (TNOW-TTTIME(INDTT-1)) 
     &        + TTB * (TNOW-TTTIME(INDTT-1))*(TNOW-TTTIME(INDTT))
          DTTEFF(I,J) = TTA + TTB *(2.0*TNOW-TTTIME(INDTT-1)
     &        -TTTIME(INDTT))
        END DO
      END DO
*
*     Get the eigenvalues and eigenvectors for current tensor
      call jacobi_transform(tteff,3,3,tteigen,tteigenv,nrotjac)
* 
      RTIDE = 0.0
      NK = 0
      DO K=1,3
         TTRTIDE(K) = (ZMASS/abs(TTEIGEN(k)))**ONE3
         IF (TTEIGEN(K).LT.0) THEN
            TTRTIDE(K) = -TTRTIDE(K)
         ELSE
            RTIDE = RTIDE + TTRTIDE(K)
            NK = NK + 1
         END IF
      END DO

      IF(NK.GT.0) THEN
         RTIDE = RTIDE/float(NK)
      ELSE
*     If tidal force not expand cluster, just use huge value of RTIDE
         RTIDE = 1.0E6*RSCALE
      END IF

      RETURN

      END
