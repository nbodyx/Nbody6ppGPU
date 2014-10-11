      SUBROUTINE JPRED(I,TTIME,TTRES)
*
*
*       Neighbour prediction of single particle.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
      REAL*8 TPRED,TTIME,TTRES
      LOGICAL ipredall
*
*
*       Skip prediction if done at current time.
*      IF (predall.or.I.LT.ifirst) return
*     --03/13/14 21:22-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(ttime.le.
c$$$     &     29.047119140625000.and.name(i).eq.12503) then
c$$$        write(120+rank,*),'t',ttime,'tpr',tpred(i),'xd',xdot(1,2*(i-n)),
c$$$     &        xdot(1,2*(i-n)-1),'x0d',x0dot(1,i),'x0',x0(1,i),
c$$$     &        'f',f(1,i),'fd',fdot(1,i),'nb',list(1,2*(i-n)-1),
c$$$     &        'U',U(1,i-n),'UD',UDOT(1,i-n),'i1',i-n,'t0',t0(i-n)
c$$$         call flush(120+rank)
c$$$         if(ttime.ge.29.047092579805405.and.
c$$$     &        tpred(i).ge.29.047090694399618) call abort()
c$$$      end if
*     --03/13/14 21:22-lwang-end----------------------------------------*
      IF(I.LT.ifirst) return
*     --03/15/14 20:05-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(x0(1,I).ge.1000.and.body(i).gt.0) call abort()
*     --03/15/14 20:05-lwang-end----------------------------------------*
      IF (TPRED(I).EQ.TIME) THEN
         IF (TPRED(I).EQ.T0(I).and.X(1,I).NE.X0(1,I)) THEN
            X(1:3,I) = X0(1:3,I)
            XDOT(1:3,I) = X0DOT(1:3,I)
*     Resolve the components of any perturbed pair.
            IF (I.GT.N) THEN
               JPAIR = I - N
               IF (LIST(1,2*JPAIR-1).GT.0) THEN
*!$omp critical
                  ZZ = 1.0
*     Distinguish between low and high-order prediction of U & UDOT.
                  IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
                  CALL KSRES2(JPAIR,J1,J2,ZZ,TTIME)
*!$omp end critical
               END IF
            END IF
         END IF
         
         RETURN
      ELSE
*     Adopt low order prediction for standard single particle.
         S = TTIME - T0(I)
         S1 = 1.5*S
         S2 = 2.0*S
*     Note X may become corrected value X0 for zero interval.
         X(1,I) = ((FDOT(1,I)*S + F(1,I))*S + X0DOT(1,I))*S + X0(1,I)
         X(2,I) = ((FDOT(2,I)*S + F(2,I))*S + X0DOT(2,I))*S + X0(2,I)
         X(3,I) = ((FDOT(3,I)*S + F(3,I))*S + X0DOT(3,I))*S + X0(3,I)
         XDOT(1,I) = (FDOT(1,I)*S1 + F(1,I))*S2 + X0DOT(1,I)
         XDOT(2,I) = (FDOT(2,I)*S1 + F(2,I))*S2 + X0DOT(2,I)
         XDOT(3,I) = (FDOT(3,I)*S1 + F(3,I))*S2 + X0DOT(3,I)
*     
*     Resolve the components of any perturbed pair.
         IF (I.GT.N) THEN
            JPAIR = I - N
            IF (LIST(1,2*JPAIR-1).GT.0) THEN
*!$omp critical
               ZZ = 1.0
*     Distinguish between low and high-order prediction of U & UDOT.
               IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
               CALL KSRES2(JPAIR,J1,J2,ZZ,TTIME)
*!$omp end critical
            END IF
         END IF
*
*     --03/07/14 21:11-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         IF(I.EQ.12195) THEN
c$$$            print*,rank,'JP I',I,'N',NAME(I),'X',X(1,I),'X1',X(1,J1),
c$$$     &           'X2',X(1,J2),'T',TIME,'TPR',TPRED(I),'T0',T0(I),
c$$$     &           'x0',x0(1,i)
c$$$            call flush(6)
c$$$         end if
*     --03/07/14 21:11-lwang-end----------------------------------------*
      TPRED(I) = TTIME
      END IF
*
      RETURN
*
      END
