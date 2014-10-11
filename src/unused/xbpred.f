      subroutine xbpred(IR,N0,NL,nxtlst)
*
*
*     Predict x and xdot. (L.WANG)
*     ------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      INCLUDE 'omp_lib.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
      COMMON/SLOW0/  RANGE,ISLOW(10)
*      REAL*8 UI(4),RDOT(3),V(4),A1(3,4)
      REAL*8 TPRED,TRES
      LOGICAL ipfirst,ipredall
      INTEGER N0,NL,NXTLST(NMAX)
      SAVE ipfirst
      DATA ipfirst/.true./

*     Initialize logical array only once
      IF (ipfirst) then
         TPRED(1:NMAX) = 0.D0
         TRES(1:KMAX) = 0.D0
         ipfirst = .false.
*         NLIMIT = MAX(10,NTOT/(5*LMAX))
      end if

*     DEBUGING-----------------------------------------------
*      DO I = 1, NTOT
*         if (.NOT.NOPRED(I)) print*,'I',I,' MISS'
*      END DO

*      if (time.ge.0.35) then
*      print*,rank,'1409,x,xdot,f,fdot',x(1,1409),xdot(1,1409),f(1,1409)
*     &        ,fdot(1,1409),time,timenw(1409),t0(1409),step(1409)
*      print*,rank,'705+N,x,xdot,f,fdot',x(1,705+N),xdot(1,705+N),
*     &   f(1,705+N),fdot(1,705+N),time,timenw(705+N),t0(705+N),
*     &     step(705+N)
*      end if
*     DEBUGING-----------------------------------------------      

*     only do nnb prediction if iphase .ne 2 and no regular step
*      IF (IR.eq.0.and.IQ.ne.2) THEN
      IF (IR.eq.0) THEN
         call cputim(tttprea)
         ipredall=.false.
*         IF (NPAIRS.GT.0) THEN
**!$omp parallel do if(NPAIRS.GE.ITHREAD) default(shared) private(I)
*            DO I = 1, NPAIRS
*               NOPRED(N+I) = .true.
*            END DO
**!$omp end parallel do
*         END IF

!$omp parallel do private(LL,I,S,K,NNB1,J,S1,S2,JPAIR,ZZ,J1,J2)
!$omp& reduction(+:NBPRED)
         DO LL = N0, NL
*     Predict current state vector of body #I to order FDOT.
            I = NXTLST(LL)
            IF (TPRED(I).ne.TIME) THEN
               S = TIME - T0(I)
               DO K = 1,3
                  X(K,I) = ((FDOT(K,I)*S + F(K,I))*S + X0DOT(K,I))*S
     *                 + X0(K,I)
                  XDOT(K,I) = (3.0*FDOT(K,I)*S + 2.0*F(K,I))*S
     *                 + X0DOT(K,I)
               END DO
               NBPRED = NBPRED + 1
               TPRED(I) = TIME
            END IF
            IF (I.GT.N) THEN
               JPAIR = I - N
               IF (tres(jpair).ne.time.and.
     &              LIST(1,2*JPAIR - 1).GT.0) THEN
                  ZZ = 1.0
                  IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
                  CALL KSRES2(JPAIR,J1,J2,ZZ,TIME)
                  tres(jpair) = time
               END IF
            END IF
*            
            NNB1 = LIST(1,I) + 1
            DO L = 2,NNB1
               J = LIST(L,I)
               IF (TPRED(J).ne.TIME) THEN
                  S = TIME - T0(J)
                  S1 = 1.5*S
                  S2 = 2.0*S
                  X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S
     *                 +X0(1,J)
                  X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S
     *                 +X0(2,J)
                  X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S
     *                 +X0(3,J)
                  XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
                  XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
                  XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
                  NBPRED = NBPRED + 1
                  TPRED(J) = TIME
               END IF
               IF (J.GT.N) THEN
                  JPAIR = J - N
                  IF (tres(jpair).ne.time.and.
     &                 LIST(1,2*JPAIR - 1).GT.0) THEN
                     ZZ = 1.0
                     IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
                     CALL KSRES2(JPAIR,J1,J2,ZZ,TIME)
                     tres(jpair) = time
                  END IF
               END IF
            END DO
         END DO
!$omp end parallel do
         call cputim(tttpreb)
         ttprei = ttprei + (tttpreb-tttprea)*60
      ELSE
         call cputim(tttprea)
         NNPRED = NNPRED + 1
         iPREDALL=.true.
!$omp parallel do private(J,S,S1,S2,JPAIR,ZZ,J1,J2)
         DO J = IFIRST,NTOT
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
                  TRES(JPAIR) = TIME
               END IF
            END IF
         END DO
!$omp end parallel do
         call cputim(tttpreb)
         ttpreall = ttpreall + (tttpreb-tttprea)*60
      END IF

*     Resolve any KS coordinates & velocities using most recent c.m.
c$$$      IF (NPAIRS.GT.0) THEN
c$$$         call cputim(tttresa)
c$$$*       Resolve perturbed KS pairs with c.m. prediction after NBSORT.
c$$$*      JJ = -1
c$$$**!$omp parallel do IF(NPAIRS.GE.ITHREAD) default(shared)
c$$$**!$omp& private(J,JPAIR,S,J1,J2,A1,A2,A3,A4,IMOD,DTU,DTU1,DTU2,
c$$$**!$omp&        Q1,Q2,Q3,UI,RDOT,V,RI,RINV)
c$$$**!$omp& reduction(+:NBPRED)
c$$$         DO 45 JPAIR = 1,NPAIRS
c$$$*      IF (.NOT.NOPRED(JJ).OR..NOT.NOPRED(JJ+1)) THEN
c$$$            IF (LIST(1,2*JPAIR - 1).GT.0) THEN
c$$$*       Ignore c.m. prediction after full N loop (all active KS needed).
c$$$               J = N + JPAIR
c$$$               IF(TPRED(J).eq.TIME) THEN
c$$$*     Predict ALL binaries even unperturbed ones for parallel code (R.Sp.)
c$$$                  ZZ = 1.0
c$$$*       Distinguish between low and high-order prediction of U & UDOT.
c$$$                  IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
c$$$                  CALL KSRES2(JPAIR,J1,J2,ZZ,TIME)
c$$$               END IF
c$$$            END IF
c$$$ 45      CONTINUE
c$$$**!$omp end parallel do
c$$$*
c$$$         call cputim(tttresb)
c$$$         ttres = ttres + (tttresb-tttresa)*60
c$$$      END IF
c$$$*       Resolve Chain if it is in block or in neighbour lists (RS Nov. 03).
c$$$      IF (NCH.GT.0) THEN
c$$$          ICHPR = 0
c$$$*       First check whether chain c.m. is in current block (ICHPR = 2).
c$$$          IF (TIMENW(ICH).EQ.TIME) ICHPR = 2
c$$$*       Second check whether neighbour lists contain chain c.m. (ICHPR = 1).
c$$$          IF (ICHPR.EQ.0) THEN
c$$$              DO 47 L = 1,NXTLEN
c$$$                  I = NXTLST(L)
c$$$                  NNB1 = LIST(1,I) + 1
c$$$*       Second check whether neighbour lists contain chain c.m. (ICHPR = 1).
c$$$                  DO 48 K = 2,NNB1
c$$$                      J = LIST(K,I)
c$$$                      IF (J.GT.ICH) GO TO 47
c$$$                      IF (J.EQ.ICH) ICHPR = 1
c$$$ 48               CONTINUE
c$$$ 47           CONTINUE
c$$$          END IF
c$$$          IF (ICHPR.GT.1) CALL CHLIST(ICH)
c$$$          IF (ICHPR.GT.0) CALL XCPRED(0)
c$$$      END IF
c$$$*
      RETURN

      END
      
