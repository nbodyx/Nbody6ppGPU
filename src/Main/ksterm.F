      SUBROUTINE KSTERM
*
*
*       Termination of KS regularization.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*8  DSAVE(16)
      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
      LOGICAL ipredall
      REAL*8 TPRED
*
*
*       Copy pair index from COMMON save and define KS components & c.m.
      IPAIR = KSPAIR
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      ICM = N + IPAIR
c$$$      JMIN = 0
      call JPRED(ICM,TIME,TIME)
*
*     --02/27/14 21:23-lwang-modify-------------------------------------*
***** Note: Suppress this to avoid new chain happened in impact, which is
*****       called by ksint here
c$$$*       Prepare termination at block time (KS, triple, quad, chain or merge).
c$$$      IF (TIME.LE.TBLOCK.AND.IPHASE.LT.9) THEN
c$$$          TIME0 = TBLOCK
c$$$*       Skip KS integration for unperturbed orbit or T0(I1) at end of block.
c$$$          IF (LIST(1,I1).EQ.0.OR.T0(I1).EQ.TIME0) GO TO 3
c$$$*
c$$$*       See whether the time interval should be modified by KSLOW procedure.
c$$$          IF (KSLOW(IPAIR).GT.1) THEN
c$$$              IMOD = KSLOW(IPAIR)
c$$$              ZMOD = FLOAT(ISLOW(IMOD))
c$$$          ELSE
c$$$              ZMOD = 1.0
c$$$          END IF
c$$$*
c$$$    1     DT0 = TIME0 - T0(I1)
c$$$*       Integrate up to current block time in case interval is too large.
c$$$          IF (DT0.GT.STEP(I1)) THEN
c$$$              TIME = T0(I1) + STEP(I1)
c$$$              H0(IPAIR) = H(IPAIR)
c$$$              Z = -0.5D0*H(IPAIR)*DTAU(IPAIR)**2
c$$$              CALL STUMPF(IPAIR,Z)
c$$$              CALL KSINT(I1)
c$$$              DTU = DTAU(IPAIR)
c$$$              STEP(I1) = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU
c$$$     &                                                   + R(IPAIR))*DTU
c$$$              STEP(I1) = ZMOD*STEP(I1)
c$$$*       Restrict increase of R for superfast particles in one block-step.
c$$$              IF (H(IPAIR).GT.100.0.AND.R(IPAIR).GT.RMIN) GO TO 3
c$$$              GO TO 1
c$$$          END IF
c$$$*
c$$$*       Determine the last regularized step by Newton-Raphson iteration.
c$$$          DTU = DT0/(R(IPAIR)*ZMOD)
c$$$          DTU = MIN(DTU,DTAU(IPAIR))
c$$$*       Include rare case of zero interval due to subtracting large values.
c$$$          DTU = MAX(DTU,1.0D-10)
c$$$          ITER = 0
c$$$    2     Y0 = DT0 - ZMOD*((ONE6*TDOT3(IPAIR)*DTU +
c$$$     &                             0.5*TDOT2(IPAIR))*DTU + R(IPAIR))*DTU
c$$$          YPR = -((0.5*TDOT3(IPAIR)*DTU + TDOT2(IPAIR))*DTU + R(IPAIR))
c$$$          YPR = ZMOD*YPR
c$$$          DTU = DTU - Y0/YPR
c$$$          DT1 = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU +
c$$$     &                                                     R(IPAIR))*DTU
c$$$          DT1 = ZMOD*DT1
c$$$          ITER = ITER + 1
c$$$          IF (ABS(DT0 - DT1).GT.1.0E-10*STEP(I1).AND.ITER.LT.10) GO TO 2
c$$$*
c$$$*       Advance the KS solution to next block time and terminate at TIME0.
c$$$          DTAU(IPAIR) = DTU
c$$$          STEP(I1) = DT1
c$$$          TIME = T0(I1) + DT1
c$$$          H0(IPAIR) = H(IPAIR)
c$$$          Z = -0.5D0*H(IPAIR)*DTU**2
c$$$          CALL STUMPF(IPAIR,Z)
c$$$          CALL KSINT(I1)
c$$$    3     TIME = TIME0
c$$$*
c$$$*       Predict X & XDOT for body #JCOMP (note TIME = TBLOCK if second call).
c$$$          IF (JCOMP.GE.IFIRST) THEN
c$$$              CALL XVPRED(JCOMP,-1)
c$$$              IF (GAMMA(IPAIR).GT.0.2.AND.JCOMP.LE.N) THEN
c$$$                  JMIN = JCOMP
c$$$*       Initialize T0, X0 & X0DOT for XVPRED & FPOLY on large perturbation.
c$$$                  T0(JCOMP) = TIME
c$$$                  DO 4 K = 1,3
c$$$                      X0(K,JCOMP) = X(K,JCOMP)
c$$$                      X0DOT(K,JCOMP) = XDOT(K,JCOMP)
c$$$    4             CONTINUE
c$$$              END IF
c$$$          END IF
c$$$      END IF
*     --02/27/14 21:24-lwang-end-modify---------------------------------*
*
*       Predict coordinates and evaluate potential energy w.r.t. perturbers.
*      CALL KSRES(IPAIR,J1,J2,0.0D0)
      NP = LIST(1,I1)
      DO 5 L = 1,NP
          JPERT(L) = LIST(L+1,I1)
    5 CONTINUE
      JLIST(1) = I1
      JLIST(2) = I2
      CALL NBPOT(2,NP,POT1)
*
*       Rectify the orbit to yield U & UDOT consistent with binding energy.
      CALL KSRECT(IPAIR)
*
*       Retain final KS variables for explicit restart at merge termination.
      IF (TIME.LE.TBLOCK.AND.IPHASE.EQ.6) THEN
          HM(NMERGE) = H(IPAIR)
          DO 6 K = 1,4
              UM(K,NMERGE) = U(K,IPAIR)
              UMDOT(K,NMERGE) = UDOT(K,IPAIR)
    6     CONTINUE
      END IF
*
*       Check optional diagnostic output for disrupted new hard binary.
      IF (KZ(8).EQ.0) GO TO 10
      IF (LIST(2,I1+1).NE.0.OR.H(IPAIR).GT.0.0) GO TO 10
      IF (GAMMA(IPAIR).GT.0.5.AND.JCOMP.GT.0.OR.IPHASE.EQ.7) THEN
          IF (JCOMP.EQ.0.OR.IPHASE.EQ.7) JCOMP = I1
          K = 0
          IF (JCOMP.GT.N) THEN
              J2 = 2*(JCOMP - N)
              K = LIST(2,J2)
          END IF
          SEMI = -0.5*BODY(ICM)/H(IPAIR)
          EB = -0.5*BODY(I1)*BODY(I2)/SEMI
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
          if(rank.eq.0)
     &    WRITE (8,8)  TIME+TOFF, NAME(I1), NAME(I2), K, BODY(I1),
     &                 BODY(I2), EB, SEMI, R(IPAIR), 
     &                 GAMMA(IPAIR), RI
 8        FORMAT ('END BINARY   Time[NB] ',1P,E25.16,0P,'  NAME(I1) '
     &         ,I11,'  NAME(I2) ',I11,'  FLAG-RB ',I3,'  M(I1)[NB] ',
     &         1P,E13.4,'  M(I2)[NB] ',E13.4,'  EB[NB] ',E20.11,
     &         '  SEMI[NB] ',E15.6,'  R12[NB] ',E15.6,'  GAMMA[NB] ',
     &         E13.4,'  RI[NB] ',E15.6,0P)
      END IF
*
   10 IF (KZ(10).GT.1) THEN
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
          if(rank.eq.0)
     &    WRITE (6,15)  TIME+TOFF, NAME(I1),NAME(I1+1),
     &         BODY(I1)*ZMBAR, BODY(I1+1)*ZMBAR,
     &          DTAU(IPAIR),R(IPAIR), RI, H(IPAIR), IPAIR, NAME(ICM),
     &         GAMMA(IPAIR), STEP(I1), LIST(1,I1), LIST(1,ICM)
   15     FORMAT (/,' END KSREG   TIME[NB]',1P,E17.10,0P,' NAME(I1)',
     &         I10,' NAME(I2)',I10,' M(I1)[M*]',
     &         F9.3,' M(I2)[M*]',F9.3,' DTAU',F8.3,' R12[NB]',1P,E9.1
     &         ,0P,' RI[NB]',F7.2,' H',F9.2,' IPAIR',I9,' NAME(ICM)'
     &         ,I10,' GAMMA',F8.3,' STEP(I1)',1P,E11.1,0P,' NPERT',I5
     &         ,' NB(ICM)',I5)
*     --04/15/14 21:40-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          if(rank.eq.0.and.name(icm).eq.2567) call abort()
c$$$          if(rank.eq.0.and.name(icm).eq.11569)
c$$$     &         print*,'IFIRST',IFIRST,'N',NAME(IFIRST),'IPAIR',IPAIR,'N',
c$$$     &         ,NAME(2*IPAIR-1),'LIST',LIST(1:6,I1),
c$$$     &     'N',(NAME(LIST(2+K,I1)),K=1,5)
*     --04/15/14 21:40-lwang-end----------------------------------------*
          call flush(6)
      END IF
*
*       Obtain global coordinates & velocities.
      CALL RESOLV(IPAIR,2)
*
*       Correct for differential potential energy due to rectification.
      CALL NBPOT(2,NP,POT2)
*       Add correction term with opposite sign for conservation.
*     ECOLL = ECOLL + (POT2 - POT1)
*     if(rank.eq.0)then
*     IF (ABS(POT1-POT2).GT.0.0001) WRITE (6,16)  POT1,BE(3),POT1-POT2
*  16 FORMAT (' CORRECT:    POT1 BE3 POT1-POT2  ',2F10.6,F10.6)
*     end if
*
*       Modify c.m. neighbour radius by density contrast and set new values.
      NNB = LIST(1,ICM) + 1
*       Check predicted neighbour number and form volume ratio.
*      NBP = MIN(ALPHA*SQRT(FLOAT(NNB)*RS(ICM))/(RS(ICM)**2),ZNBMAX)
      NBP = INT(MIN(NNB + 0.1*(FLOAT(NNBOPT) - FLOAT(NNB)),ZNBMAX))
      NBP = MAX(NBP,INT(ZNBMIN))
      A0 = FLOAT(NBP)/FLOAT(NNB)
*       Re-determine neighbour list on zero membership for distant binary.
      IF (LIST(1,ICM).EQ.0) THEN
          RS0 = 0.1*(ABS(X(1,ICM)) + ABS(X(2,ICM)))
          call xbpredall
          CALL NBLIST(ICM,RS0)
          NNB = LIST(1,ICM) + 1
      END IF
*       Copy all neighbours in the case of merger.
      IF (IPHASE.EQ.6) THEN
          A0 = 1.0
          NBP = NNB - 1
      END IF
      IF (RS(ICM).GT.-100.0*BODY(ICM)/H(IPAIR)) A0 = 1.0
*       Accept old c.m. values for small length scale ratio or H > 0.
      RS(I1) = RS(ICM)*A0**0.3333
      RS(I1+1) = RS(I1)
      RS2 = RS(I1)**2
*
*       Select neighbours for components inside the modified c.m. sphere.
   20 NNB1 = 1
      DO 25 L = 2,NNB
          J = LIST(L,ICM)
          call jpred(J,time,time)
          RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
     &                                    (X(3,ICM) - X(3,J))**2
*       Ensure that at least the predicted neighbour number is reached.
          IF (RIJ2.LT.RS2.OR.L + NBP.GT.NNB1 + NNB) THEN
              NNB1 = NNB1 + 1
              ILIST(NNB1) = J
          END IF
   25 CONTINUE
*
*       Check that there is space for adding dominant component later.
      IF (NNB1.GE.NNBMAX.AND.IPHASE.NE.6) THEN
          RS2 = 0.9*RS2
          GO TO 20
      END IF
*
*       Reduce pair index, total number & single particle index.
      NPAIRS = NPAIRS - 1
      NTOT = N + NPAIRS
      IFIRST = 2*NPAIRS + 1
*
*       Save name of components & flag for modifying LISTD in UPDATE.
      JLIST(1) = NAME(I1)
      JLIST(2) = NAME(I1+1)
      JLIST(3) = LIST(2,I1+1)
*
*       Skip adjustment of tables if last or only pair being treated.
      IF (IPAIR.EQ.NPAIRS+1) GO TO 60
*
*     Exchange the last pair and terminated pair
      DO 50 KCOMP = 0,1
*     First save ipair data
          I = 2*IPAIR - KCOMP
*
          DSAVE(1:3) = X(1:3,I)
          DSAVE(4:6) = X0DOT(1:3,I)
*       Current velocity has been set in routine RESOLV.
          DSAVE(7) = BODY(I)
          DSAVE(8) = RS(I)
          DSAVE(9) = RADIUS(I)
          DSAVE(10) = TEV(I)
          DSAVE(11) = BODY0(I)
          DSAVE(12) = TEV0(I)
          DSAVE(13) = EPOCH(I)
          DSAVE(14) = SPIN(I)
          DSAVE(15) = ZLMSTY(I)
          DSAVE(16) = TPRED(I)
          NAMEI = NAME(I)
          KSI = KSTAR(I)
*
*       Copy last pair to ipair
          J = 2*NPAIRS + 2 - KCOMP
          X(1:3,I) = X(1:3,J)
*       Copy latest X & X0DOT (= 0) of single components for predictor
          X0(1:3,I) = X(1:3,I)
          X0DOT(1:3,I) = X0DOT(1:3,J)
          XDOT(1:3,I) = XDOT(1:3,J)
          BODY(I) = BODY(J)
          RS(I) = RS(J)
          RADIUS(I) = RADIUS(J)
          TEV(I) = TEV(J)
          BODY0(I) = BODY0(J)
          TEV0(I) = TEV0(J)
          EPOCH(I) = EPOCH(J)
          SPIN(I) = SPIN(J)
          ZLMSTY(I) = ZLMSTY(J)
          TPRED(I) = TPRED(J)
          NAME(I) = NAME(J)
          KSTAR(I) = KSTAR(J)
          STEP(I) = STEP(J)
          T0(I) = T0(J)
          K = LIST(1,J) + 1
          IF (K.EQ.1) K = 2
*       Transfer unmodified neighbour lists (include flag in 2nd comp).
          LIST(1:K,I) = LIST(1:K,J)
*     
*       Set new component index and copy basic variables.
          X(1:3,J) = DSAVE(1:3)
          X0DOT(1:3,J) = DSAVE(4:6)
          XDOT(1:3,J) = DSAVE(4:6)
          BODY(J) = DSAVE(7)
          RS(J) = DSAVE(8)
          RADIUS(J) = DSAVE(9)
          TEV(J) = DSAVE(10)
          BODY0(J) = DSAVE(11)
          TEV0(J) = DSAVE(12)
          EPOCH(J) = DSAVE(13)
          SPIN(J) = DSAVE(14)
          ZLMSTY(J) = DSAVE(15)
          TPRED(J) = DSAVE(16)
          NAME(J) = NAMEI
          KSTAR(J) = KSI

   50 CONTINUE
*
*       Include removal of the circularization name NAMEC from chaos table.
      IF (KSTAR(ICM).GE.10.AND.NCHAOS.GT.0.AND.IPHASE.NE.6) THEN
*       Note that NAMEC may remain dormant for hierarchical systems.
          II = -ICM
          CALL SPIRAL(II)
      END IF
*
*       Remove regularized variables and c.m. of ipair and move last pair data to ipair position.
 60   CALL REMOVE_KS(IPAIR)
*
*       Set new global index of first & second component.
      ICOMP = 2*NPAIRS + 1
      JCOMP = ICOMP + 1
*
*       Save c.m. neighbour list for routine FPOLY1/2 (may be renamed below).
      ILIST(1) = NNB1 - 1
      LIST(1:NNB1,ICOMP) = ILIST(1:NNB1)
*
*       Modify all relevant COMMON list arrays.
      CALL UPDATE(IPAIR)
*
*      Update the TLIST in intgrt
C      CALL repair_tlist(ICOMP,JCOMP,ICM,2)
*     --04/15/14 21:40-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          if(rank.eq.0)
c$$$     &     print*,'NEW IFIRST',IFIRST,'N',NAME(IFIRST),'IPAIR',IPAIR,
c$$$     &     'N',NAME(2*IPAIR-1)
*     --04/15/14 21:40-lwang-end----------------------------------------*
*       Check replacing of single KS component by corresponding c.m.
   70 IF (LIST(2,ICOMP).LT.ICOMP) THEN
          J2 = LIST(2,ICOMP)
          J = KVEC(J2) + N
          DO 80 L = 2,NNB1
              IF (L.LT.NNB1.AND.LIST(L+1,ICOMP).LT.J) THEN
                  LIST(L,ICOMP) = LIST(L+1,ICOMP)
              ELSE
                  LIST(L,ICOMP) = J
              END IF
   80     CONTINUE
*       Check again until first neighbour > ICOMP.
          GO TO 70
      END IF
*
*       Make space for dominant component and copy members to JCOMP list.
      DO 90 L = NNB1,2,-1
          LIST(L+1,ICOMP) = LIST(L,ICOMP)
          LIST(L+1,JCOMP) = LIST(L,ICOMP)
   90 CONTINUE
*
*       Set dominant component in first location and specify membership.
      LIST(2,ICOMP) = JCOMP
      LIST(2,JCOMP) = ICOMP
      LIST(1,ICOMP) = NNB1
      LIST(1,JCOMP) = NNB1
*
*       Initialize T0, T0R, X0 & X0DOT for both components.
      T0(ICOMP) = TIME
      T0(JCOMP) = TIME
      T0R(ICOMP) = TIME
      T0R(JCOMP) = TIME
      DO 95 K = 1,3
          X0(K,ICOMP) = X(K,ICOMP)
          X0(K,JCOMP) = X(K,JCOMP)
          X0DOT(K,ICOMP) = XDOT(K,ICOMP)
          X0DOT(K,JCOMP) = XDOT(K,JCOMP)
   95 CONTINUE
*
*       Form new force polynomials (skip triple, quad, merge & chain).
      IF (IPHASE.LT.4) THEN
*       Predict current coordinates & velocities for the neighbours.
          CALL XVPRED(ICOMP,NNB1)
*     --09/25/13 19:34-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(time.ge.12.920379638671875) then
c$$$         do i=ifirst,ntot
c$$$            if(name(i).eq.560) then
c$$$           print*,rank,'fp1',i,name(i),x0(1,i),x0dot(1,i),
c$$$     *      step(i),t0(i),f(1,i),fdot(1,i),time
c$$$           end if
c$$$        end do
c$$$      end if
*     --09/25/13 19:35-lwang-end----------------------------------------*
*
*       Obtain new polynomials & steps.
          CALL FPOLY1_KS(ICOMP,JCOMP,2)
*     --09/09/13 22:05-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*          if (name(icomp).eq.1281.and.time.ge.0.6) then
*            print *,rank,'ks1 1281',icomp,fdot(1,icomp),frdot(1,icomp),
*     &            time,step(icomp)
*            call xbpredall
*            CALL FPOLY2(ICOMP,JCOMP,2)
*            print *,rank,'ks 1281',icomp,fdot(1,icomp),frdot(1,icomp),
*     &           time,step(icomp)
*            stop
*          end if
*     --09/09/13 22:05-lwang-end----------------------------------------*
          CALL FPOLY2_KS(ICOMP,JCOMP,2)
*     --09/09/13 22:05-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(time.ge.12.920379638671875) then
c$$$         do i=ifirst,ntot
c$$$            if(name(i).eq.560) then
c$$$           print*,rank,'fp2',i,name(i),x0(1,i),x0dot(1,i),
c$$$     *      step(i),t0(i),f(1,i),fdot(1,i),time
c$$$           end if
c$$$        end do
c$$$      end if
*     --09/09/13 22:05-lwang-end----------------------------------------*
*
c$$$*       Improve force polynomials of strong perturber after rectification.
c$$$          IF (JMIN.GE.IFIRST) THEN
c$$$              NNB1 = LIST(1,JMIN)
c$$$              call xvpred(JMIN, NNB1)
c$$$              CALL FPOLY1_KS(JMIN,JMIN,0)
c$$$              CALL FPOLY2_KS(JMIN,JMIN,0)
c$$$          END IF
       ELSE
*     Set artifical step for icomp and jcomp thus can be removed later in NXTLST
          STEP(ICOMP) = 2*SMAX
          STEP(JCOMP) = 2*SMAX
      END IF

*     Add two particles into NXTLST
      IF(TIME+TOFF.NE.0.D0) THEN
         call add_tlist(ICOMP,STEP,DTK)
         call add_tlist(JCOMP,STEP,DTK)
         call shrink_tlist
      END IF
*
#ifdef SIMD
*     Update Particle data in AVX/SSE library
      IF(TIME+TOFF.NE.0.0D0) THEN
         CALL IRR_SIMD_SET_JP(ICOMP,X0(1,ICOMP),X0DOT(1,ICOMP),
     &        F(1,ICOMP),FDOT(1,ICOMP),BODY(ICOMP),T0(ICOMP))
         CALL IRR_SIMD_SET_LIST(ICOMP,LIST(1,ICOMP))
         CALL IRR_SIMD_SET_JP(JCOMP,X0(1,JCOMP),X0DOT(1,JCOMP),
     &        F(1,JCOMP),FDOT(1,JCOMP),BODY(JCOMP),T0(JCOMP))
         CALL IRR_SIMD_SET_LIST(JCOMP,LIST(1,JCOMP))
      END IF
#endif
*     --09/25/13 18:02-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if (rank.eq.0) then
c$$$         print*,'N',name(ICOMP),name(JCOMP),'S',stepr(icomp),step(icomp)
c$$$     &        ,stepr(jcomp),step(jcomp),time,'t0',t0(Icomp),t0(jcomp),
c$$$     &        'iphase',iphase
c$$$         call flush(6)
c$$$         if(stepr(jcomp).ge.100) then
c$$$            call abort()
c$$$         end if
c$$$      end if
c$$$      j = 15967
c$$$      print*,rank,'ksterm ',j,'n',name(j),'x',x(1,j),'xdot',xdot(1,j),
c$$$     &     'fdot',fdot(1,j),'time',time,'iphase',iphase
c$$$      call flush(6)
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)
c$$$      if (time.ge.0.73070144653320313) then
c$$$         itest = 0
c$$$         print*,body(itest)
c$$$      end if
*     --09/25/13 18:02-lwang-end----------------------------------------*
*       Check updating of global index for chain c.m.
      IF (NCH.GT.0) THEN
          CALL CHFIND
      END IF
*
      RETURN
*
      END
