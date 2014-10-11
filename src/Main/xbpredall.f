      subroutine xbpredall
*
*
*     Predict x and xdot. (L.WANG)

      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
      REAL*8 TPRED
      LOGICAL iPREDALL

      IF (IPREDALL) RETURN

      NNPRED = NNPRED + 1
!$omp parallel do private(J,S,S1,S2,JPAIR,J1,J2,ZZ)
      DO 40 J = IFIRST,NTOT
*     IF(TPRED(J).NE.TIME) THEN
         S = TIME - T0(J)
         S1 = 1.5*S
         S2 = 2.0*S
         X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
         X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
         X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
         XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
         XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
         XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
         TPRED(J) = TIME
         IF (J.GT.N) THEN
            JPAIR = J - N
            IF (LIST(1,2*JPAIR - 1).GT.0) THEN
               ZZ = 1.0
               IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
               CALL KSRES2(JPAIR,J1,J2,ZZ,TIME)
            END IF
*     --03/07/14 21:22-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            IF(J.EQ.12195) THEN
c$$$               print*,rank,'PA I',J,'N',NAME(J),'X',X(1,J),'X1',X(1,J1),
c$$$     &              'X2',X(1,J2),'T',TIME,'TPR',TPRED(J),ipredall,
c$$$     &              'T0',T0(J),'x0',x0(1,j),'F',F(1,j),'FD',FDOT(1,j)
c$$$               call flush(6)
c$$$            end if
*     --03/07/14 21:22-lwang-end----------------------------------------*
         END IF
 40   CONTINUE
!$omp end parallel do
*     --05/12/14 17:26-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      do KK=ifirst,ntot
c$$$         write(102+rank,*) 'I',KK,'N',NAME(KK),'X',X(1:3,KK),
c$$$     &        'XD',XDOT(1:3,KK),'X0',X0(1:3,KK),'X0D',X0DOT(1:3,kk),
c$$$     &        'T0',T0(KK)
c$$$         if(KK.GT.N) then
c$$$            jpair=kk-n
c$$$            k1=2*(kk-n)-1
c$$$            k2=k1+1
c$$$            if(list(1,k1).gt.0) then
c$$$            write(102+rank,*) 'P',K1,K2,'N',NAME(K1),NAME(K2),'X1',
c$$$     &           X(1:3,K1),'X2',X(1:3,K2),'XD1',XDOT(1:3,K1),
c$$$     &              'XD2',XDOT(1:3,K2),'U0',U0(1:4,KK-n),'UDOT',
c$$$     &              UDOT(1:4,jpair),'FU',FU(1:4,jpair),'FUDOT',
c$$$     &              FUDOT(1:4,jpair),'FUDOT2',FUDOT2(1:4,jpair),
c$$$     &              'T',TIME
c$$$            end if
c$$$         end if
c$$$      end do
*     --05/12/14 17:26-lwang-end----------------------------------------*
      iPREDALL = .true.

      return

      end
      
