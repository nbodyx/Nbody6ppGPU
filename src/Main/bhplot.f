      SUBROUTINE BHPLOT
*
*
*       Black hole plotting data.
*       ----------------------------------
*
      INCLUDE 'common6.h'
*      COMMON/IMBH/ BIMBH,DTBH,NIMBH
      INTEGER IBH
      REAL*8 TBH,RREL(3),VREL(3),RAVE(3),VAVE(3),VSIGMA(3)
      REAL*8 MT,DEN,RIJMAX
      LOGICAL LOOP_FLAG
      CHARACTER*8 WHICH
      SAVE TBH,IBH
      DATA TBH,IBH/0.0,-1/
*
*     Reset TBH after restart.
      IF (TBH.EQ.0.D0.AND.TIME+TOFF.NE.0.D0.AND.KSTART.GT.1) 
     &     TBH = INT(TIME/DTBH)*DTBH + TOFF
*     --07/02/14 17:25-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
C      TBH = 8.25
*     --07/02/14 17:25-lwang-end----------------------------------------*

*     IF not reach next output time, return
      IF (TBH.GT.TIME+TOFF) RETURN
*     
*     Output labels
      IF(IBH.EQ.-1.and.KSTART.EQ.1.and.rank.eq.0) then
         write(45,*)'     STAT         Time[Myr]                 IBH',
     &        '         X(1)[PC]                   X(2)[PC]        ',
     &        '         X(3)[PC]                  V(1)[km/s]       ',
     &        '        V(2)[km/s]                 V(3)[km/s]       ',
     &        '          NB        XAVE(1)[AU]       ',
     &        '        XAVE(2)[AU]               XAVE(3)[AU]       ',
     &        '       VAVE(1)[km/s]             VAVE(2)[km/s]      ',
     &        '       VAVE(3)[km/s]             DEN[M*/PC^3]       ',
     &        '       RIJMAX[PC]            VSIGMA(1)[km^2/s^2]    ',
     &        '   VSIGMA(2)[km^2/s^2]       VSIGMA(3)[km^2/s^2]    '
      END IF
*
*     Search for massive black hole
      IF(IBH.EQ.-1) IBH = NIMBH
      IF(IBH.GE.1.AND.NAME(IBH).EQ.NIMBH) GO TO 1 
*     
*     Search index of IMBH
      IF(ABS(N-IBH).GT.ABS(IFIRST-IBH)) THEN
         IBH = 0
         DO K = 1,N
            IF (NAME(K).EQ.NIMBH) THEN
               IBH = K
               GO TO 1
            END IF
         END Do
      ELSE
         IBH = 0
         DO K = N,1,-1
            IF (NAME(K).EQ.NIMBH) THEN
               IBH = K
               GO TO 1
            END IF
         END DO
      END IF
*
c$$$      if(rank.eq.0) then
c$$$         print*,'BH not found! N',NIMBH
c$$$         call flush(6)
c$$$         stop
c$$$      end if
      
*     Output data
 1    IF(IBH.GT.0) THEN
         IBL = IBH
         RAVE(1:3) = 0.D0
         VAVE(1:3) = 0.D0
         RIJMAX = 0.D0
         MT = 0.D0
         IF (IBH.LT.IFIRST) THEN
*     In case massive black hole is in binary system
            LOOP_FLAG = .true.
            IBL = KVEC(IBH)+N
*     Escape the chain case
            IF(NAME(IBL).LE.0) RETURN
            CALL JPRED_INT(IBL,TIME)
*     Do KS transformation even perturber list is zero
            JPAIR = 2*KVEC(IBH)-1 
            IF(LIST(1,JPAIR).EQ.0) THEN
               ZZ = 1.0
*     Distinguish between low and high-order prediction of U & UDOT.
               IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
               CALL KSRES2(JPAIR,J1,J2,ZZ,TTIME)
            END IF
*     --07/02/14 17:14-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*            IF(TIME.GE.8.2578) THEN
c$$$            write(120,*),'TIME',TIME,'IBL',IBL,'X',X(1,IBL),
c$$$     &           'XD',XDOT(1,IBL)*VSTAR,'M',BODY(IBL)*ZMBAR,
c$$$     &              'F',F(1,IBL),'FD',FDOT(1,IBL),
c$$$     &              'STEP',STEP(IBL),'STEPR',STEPR(IBL),
c$$$     &              'FI',FI(1,IBL),'FID',FIDOT(1,IBL),
c$$$     &              'T0',T0(IBL),'T0R',T0R(IBL),'NB',LIST(1,IBL),
c$$$     &              'FR',FR(1,IBL),'FRD',FRDOT(1,IBL)
C     &           'XBH',X(1,IBH),'XDBH',XDOT(1,IBH)*VSTAR,
C     &           'X2',X(1,IBH+1),'XD2',XDOT(1,IBH+1)*VSTAR,
C     &           'MBH',BODY(IBH)*ZMBAR,'M2',BODY(IBH+1)*ZMBAR
c$$$            call flush(120)
*            END IF
c$$$            IF(XDOT(1,IBL)*VSTAR.GT.100) call abort()
*     --07/02/14 17:14-lwang-end----------------------------------------*
         ELSE
            call JPRED_INT(IBL,TIME)
         END IF

*     Output neighbor information into file mbhnb.46
         NNB = LIST(1,IBL)
         NNBT = NNB

*     Print binary companion as first neighbor
         IF(IBH.LT.IFIRST) THEN
            NNBT = NNBT + 1
            IF(rank.eq.0) 
     &           WRITE (46,*) 'Time[Myr]',(TIME+TOFF)*TSTAR,'NB',NNBT
            WHICH = 'BINARY'
            IBC = IBH + 1
            IF(MOD(IBH,2).EQ.0) IBC = IBH - 1
            DO K = 1, 3
               RREL(K) = X(K,IBC) - X(K,IBH)
               VREL(K) = XDOT(K,IBC) - XDOT(K,IBH)
               RAVE(K) = BODY(IBC)*RREL(K)
               VAVE(K) = BODY(IBC)*VREL(K)
               RIJ = RREL(K)*RREL(K)
            END DO
            RIJ = SQRT(RIJ)
            RIJMAX = RIJ
            MT = BODY(IBC)
            IF(rank.eq.0)
     &           write (46,*) NAME(IBC),BODY(IBC)*ZMBAR,
     &           RREL(1:3)*RAU,RIJ*RAU,VREL(1:3)*VSTAR,KSTAR(IBC)
         ELSE
            IF(rank.eq.0)
     &           WRITE (46,*) 'Time[Myr]',(TIME+TOFF)*TSTAR,'NB',NNBT
            WHICH = 'SINGLE'
         END IF

         NNB = NNB + 1
         DO K = 2, NNB
            I = LIST(K,IBL)
            call JPRED(I,TIME,TIME)
            RIJ = 0.D0
            DO KJ = 1, 3
               RREL(KJ) = X(KJ,I) - X(KJ,IBH)
               VREL(KJ) = XDOT(KJ,I) - XDOT(KJ,IBH)
               RAVE(KJ) = RAVE(KJ) + BODY(I)*RREL(KJ)
               VAVE(KJ) = VAVE(KJ) + BODY(I)*VREL(KJ)
               RIJ = RREL(KJ)*RREL(KJ)
            END DO
            RIJ = SQRT(RIJ)
            RIJMAX = MAX(RIJMAX,RIJ)
            MT = MT + BODY(I)
            IF(rank.eq.0) 
     &           write (46,*) NAME(I),BODY(I)*ZMBAR,
     &           RREL(1:3)*RAU,RIJ*RAU,VREL(1:3)*VSTAR,KSTAR(I)
         END DO
*     Write two blank lines 
         write (46,*)
         write (46,*)
*     
*     Density center, average velocity and density of black hole neighbors
         RAVE(1:3) = RAVE(1:3)/MT
         VAVE(1:3) = VAVE(1:3)/MT
         DEN = MT*ZMBAR/(4.0/3.0*3.1415926*(RIJMAX*RBAR)**3)
*     Calculate velocity dispersion

         VSIGMA(1:3) = 0.D0
         IF(IBH.LT.IFIRST) THEN
            IBC = IBH + 1
            IF(MOD(IBH,2).EQ.0) IBC = IBH - 1
            DO K = 1, 3
               VREL(K) = XDOT(K,IBC) - XDOT(K,IBH) - VAVE(K)
               VSIGMA(K) = BODY(I)*VREL(K)*VREL(K)
            END DO
         END IF
         
         DO K = 2, NNB
            I = LIST(K,IBL)
            DO KJ = 1, 3
               VREL(KJ) = XDOT(KJ,I) - XDOT(KJ,IBH) - VAVE(KJ)
               VSIGMA(KJ) = VSIGMA(KJ) + BODY(I)*VREL(KJ)*VREL(KJ)
            END DO
         END DO

         VSIGMA(1:3) = VSIGMA(1:3)/MT

         if(rank.eq.0) THEN
            WRITE (45,*) WHICH,(TIME+TOFF)*TSTAR,IBH,X(1:3,IBH)*RBAR,
     &           XDOT(1:3,IBH)*VSTAR,NNBT,RAVE(1:3)*RAU,VAVE(1:3)*VSTAR,
     &           DEN,RIJMAX*RBAR,VSIGMA(1:3)*VSTAR*VSTAR
         end if

         TBH = TBH + DTBH
      END IF
         
      RETURN

      END
