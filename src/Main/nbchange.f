      SUBROUTINE NBCHANGE(LSOLD,NOLD,LSNEW,NNEW,LIST,I1,I2)
*
*
*     Exchange neighbors from old to new (LISTNEW member do not add listnew member)
*     --------------------------
*
      include 'params.h'
      INTEGER NOLD,NNEW
      INTEGER LISTOLD(NOLD),LISTNEW(NNEW)
      INTEGER LSOLD(NOLD),ISO(NOLD),LSNEW(NNEW),ISN(NNEW)
      INTEGER LIST(LMAX,NMAX),I1,I2,I,J,IJ,K
      LOGICAL CFLAG

*     Copy the array
      LISTOLD(1:NOLD) = LSOLD(1:NOLD)
      LISTNEW(1:NNEW) = LSNEW(1:NNEW)
      DO I=1,NOLD
         ISO(I) = I
      END Do
      DO I=1,NNEW
         ISN(I) = I
      END DO
      
*     sort list
      call SORT1I(NOLD,LISTOLD,ISO)
      call SORT1I(NNEW,LISTNEW,ISN)

C*     for DEBUGING
C      DO I = I1,I2
C         NNB = LIST(1,I)
C         write(6,*) 'I ',I,' NNB ',NNB,' LIST',LIST(2:NNB+1,I)
C      END DO

C*     For debuging
C      write(6,*) 'OLD: N ',NOLD,' LIST ',LISTOLD(1:NOLD)
C      write(6,*) 'NEW: N ',NNEW,' LIST ',LISTNEW(1:NNEW)
      
*      K = -1
*      DO I = 1,NOLD
*         IF(LISTOLD(I).LE.K) then
*            write(6,*) 'Error! LISTOLD is not in increasing order!',
*     &           'LIST:',LISTOLD(1:NOLD)
*         END IF
*      END DO
* 
*      K = -1
*      DO I = 1,NNEW
*         IF(LISTNEW(I).LE.K) then
*            write(6,*) 'Error! LISTOLD is not in increasing order!',
*     &           'LIST:',LISTOLD(1:NOLD)
*            stop
*         END IF
*      END DO
      
!$omp parallel do private(I,CFLAG,NNB,K,IJ)
      DO I = I1, I2
         NNB = LIST(1,I)+1
         IF (NNB.GT.1) THEN
*     Escape the case when the first or last in neighbor list is already out of the range
            IF(.NOT.(LIST(2,I).GT.LISTOLD(NOLD).OR.
     &           LIST(NNB,I).LT.LISTOLD(1))) THEN
*     Check each neighbor
*     Remove the old if NNOLD > 0
               CFLAG=.false.
               IF(NOLD.GT.0) THEN
                  K=1
                  DO J = 2,NNB
 1                   IJ = LIST(J,I)
                     IF(IJ.EQ.LISTOLD(K)) THEN
*     find one, then shift the array back by one
                        LIST(J:NNB-1,I) = LIST(J+1:NNB,I)
                        NNB = NNB-1
                        CFLAG = .true.
                        K = K + 1
                        IF (K.GT.NOLD) GO TO 2
                        GO TO 1
                     ELSE IF(IJ.GT.LISTOLD(K)) THEN
                        K = K + 1
                        IF (K.GT.NOLD) GO TO 2
                        GO TO 1
                     END IF
                  END DO
               END IF
               
*     Do not add new member if I is in listnew
 2             DO J=1,NNEW
                  IF(I.EQ.LISTNEW(J)) CFLAG=.false.
               END DO

*     Add new if CFLAG is true and NNEW > 0
               IF(CFLAG.AND.NNEW.GT.0) THEN
                  K=1
                  DO J = 2,NNB
                     IJ = LIST(J,I)
                     IF(IJ.GT.LISTNEW(K)) THEN
                        IF(NNB+1.GT.LMAX) THEN
                           write(6,*) 'Error! LIST overflow!, I',I,
     &                          'LMAX',LMAX
                           STOP
                        END IF 
                        LIST(J+1:NNB+1,I) = LIST(J:NNB,I)
                        LIST(J,I) = LISTNEW(K)
                        NNB = NNB + 1
                        K = K + 1
                        IF (K.GT.NNEW) GO TO 3
                     ELSE IF(IJ.EQ.LISTNEW(K)) THEN
                        K = K + 1
                        IF (K.GT.NNEW) GO TO 3
                     END IF 
                  END DO
*     Add remaining
                  IF(K.LT.NNEW) THEN
                     IF(NNB+NNEW-K+1.GT.LMAX) THEN
                        write(6,*) 'Error! LIST overflow!, I',I,
     &                       'LMAX',LMAX
                        STOP
                     END IF
                     LIST(NNB+1:NNB+NNEW-K+1,I) = LISTNEW(K:NNEW)
                     NNB = NNB + NNEW-K+1
                  END IF 
               END IF
*     Renew neighbor number
 3             LIST(1,I) = NNB - 1
            END IF 
         END IF
      END DO
!$omp end parallel do


C*     for DEBUGING
C      DO I = I1,I2
C         NNB = LIST(1,I)
C         write(6,*) 'I ',I,' NNB ',NNB,' LIST',LIST(2:NNB+1,I)
C      END DO
      
      RETURN

      END
      
