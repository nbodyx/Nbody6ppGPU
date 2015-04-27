      subroutine ksrel_finspin(TIME,SPNI1,M1,M2,XI,VI)
*
*     * final spin determination after rezzolla
*---------------------------------------------
*
      include 'params.h'
      include 'mpi_base.h'
      include 'postnewton.h'
      REAL*8 TIME,M1,M2,ETA,MASAS,XI(6),VI(6),KSX(3),KSV(3)
      REAL*8 RATIO, ABSL, ABSA1, ABSA2, COSA, COSB, COSG
      REAL*8 L(3),LABS,S1(3),S2(3),KSS(3)
      REAL*8 AFIN(3),AFINABS,ROOTHE,JVEC(3),JABS
      INTEGER SPNI1,SPNI2

      SPNI2 = SPNI1 + 1

      MASAS = M1+M2
      ETA = M1*M2/(MASAS*MASAS)
      
      DO K = 1,3
         KSX(K) = XI(K)-XI(K+3)
         KSV(K) = VI(K)-VI(K+3)
      END DO

      CALL CROSSP(KSX,KSV,L)
      L(1) = L(1)*MASAS*ETA
      L(2) = L(2)*MASAS*ETA
      L(3) = L(3)*MASAS*ETA
      LABS = SQRT(L(1)*L(1)+L(2)*L(2)+L(3)*L(3))

*     lighter body over heavier body
      RATIO = M2/M1
      IF (M2.GT.M1) THEN RATIO = 1/RATIO END IF

      DO K =1,3
         S1(K) = SPN(K,SPNI1)*M1*M1/CLIGHTPN
         S2(K) = SPN(K,SPNI2)*M2*M2/CLIGHTPN
         KSS(K) = S1(K)+S2(K)
      END DO

*     WRITE(*,*) 'SPN(1,SPNI1) SPN(2,SPNI1) SPN(3,SPNI1)%%%',
*     & SPN(1,SPNI1),
*     & SPN(2,SPNI1), SPN(3,SPNI1), SPN(1,SPNI2),
*     & SPN(2,SPNI2), SPN(3,SPNI2)
*     WRITE(*,*) 'SPNI1 SPNI2%%', SPNI1, SPNI2
      ABSA1 = SQRT(SPN(1,SPNI1)**2+SPN(2,SPNI1)**2+SPN(3,SPNI1)**2)
      ABSA2 = SQRT(SPN(1,SPNI2)**2+SPN(2,SPNI2)**2+SPN(3,SPNI2)**2)
      WRITE (*,*) 'INDIVIDUAL SPINS', ABSA1,ABSA2
      IF (ABSA1.LE.1.0E-10.OR.ABSA2.LE.1.0E-10) THEN
         COSA = 0.
         COSB = 0.
         WRITE (*,*) TIME, 'PROBLEM COSA/B = 0'
      ELSE
         COSA = (SPN(1,SPNI1)*SPN(1,SPNI2)+SPN(2,SPNI1)*SPN(2,SPNI2)+
     &        SPN(3,SPNI1)*SPN(3,SPNI2))/ABSA1/ABSA2
         COSB = (SPN(1,SPNI1)*L(1)+SPN(2,SPNI1)*L(2)
     &        +SPN(3,SPNI1)*L(3))/
     &        ABSA1/LABS
      END IF
      IF (ABSA2.LE.1.0E-10.OR.LABS.LE.1.0E-10) THEN
         COSG = 0.
         WRITE (*,*) TIME, 'PROBLEM COSG = 0'
      ELSE
         COSG = (SPN(1,SPNI2)*L(1)+SPN(2,SPNI2)*L(2)
     &        +SPN(3,SPNI2)*L(3))/
     &        ABSA2/LABS
      END IF
      ABSL = -0.129/(1+RATIO*RATIO)**2*(ABSA1**2+ABSA2**2*RATIO**4+
     &     2.*ABSA1*ABSA2*RATIO**2*COSA)+(-3.84*ETA-2.686+2.)/
     &     (1+RATIO**2)*(ABSA1*COSB+ABSA2*RATIO**2*COSG)+
     &     3.4641-3.454*ETA+2.353*ETA*ETA
      ROOTHE = ABSA1**2+ABSA2**2*RATIO**4+2.*
     &     ABSA2*ABSA1*RATIO**2*COSA+2.*(ABSA1*COSB+ABSA2*RATIO**2*
     &     COSG)*ABSL*RATIO+ABSL*ABSL*RATIO*RATIO
      IF (ROOTHE.GT.0.) THEN
         AFINABS = 1/((1+RATIO)**2)*SQRT(ROOTHE)
      ELSE
         WRITE (*,*) TIME, 'PROBLEM negative root'
         IF (ABSA1.GT.ABSA2) THEN
            AFINABS = ABSA1
         ELSE
            AFINABS = ABSA2
         END IF
      END IF
      IF (AFINABS.GT.1.0) THEN
         AFINABS = 1.0
         WRITE (*,*) TIME, 'PROBLEM BIGGER THAN 1'
      END IF
      JVEC(1) = KSS(1) + L(1)
      JVEC(2) = KSS(2) + L(2)
      JVEC(3) = KSS(3) + L(3)
      JABS = SQRT(JVEC(1)**2+JVEC(2)**2+JVEC(3)**2)
      AFIN(1) = AFINABS*JVEC(1)/JABS
      AFIN(2) = AFINABS*JVEC(2)/JABS
      AFIN(3) = AFINABS*JVEC(3)/JABS
!     if we are here, we are going to merge it so give the new spin now
!     since we dont know for sure which one is the new body, just give the value to both
      SPN(1,SPNI1) = AFIN(1)
      SPN(2,SPNI1) = AFIN(2)
      SPN(3,SPNI1) = AFIN(3)
      SPN(1,SPNI2) = AFIN(1)
      SPN(2,SPNI2) = AFIN(2)
      SPN(3,SPNI2) = AFIN(3)
      IF(rank.eq.0) WRITE (51,558) TIME,(M1+M2),AFINABS,SPNI1,SPNI2
      IF(rank.eq.0) WRITE (*,*) 'L AND J', L(1), L(2), L(3), JVEC(1),
     &     JVEC(2), JVEC(3)
 558  FORMAT (3E15.6,2I5)
      CALL FLUSH(51)

      RETURN

      END
      
