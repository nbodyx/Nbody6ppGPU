      SUBROUTINE XCPRED(KCASE)
*
*
*       Prediction of global chain coordinates.
*       ---------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  M,MASS,MC,MIJ,MKK
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
*     Safe for parallel
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
*     Unsafe for parallel! update variables: xc,uc,listc
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
*     --03/15/14 21:27-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      INTEGER ICALLP
c$$$      SAVE ICALLP
c$$$      DATA ICALLP/0/
c$$$
c$$$      ICALLP = ICALLP + 1
*     --03/15/14 21:27-lwang-end----------------------------------------*
*
*
*       Check indicator for prediction of perturbers & c.m.
      IF (KCASE.EQ.0) GO TO 4
*
*       Check adding chain c.m. #ICH to perturber list for prediction.
      IF (KCASE.EQ.1) THEN
          NNB2 = LISTC(1) + 2
          LISTC(NNB2) = ICH
      ELSE
*      Note KCASE = 2 for chain c.m. prediction at new block-time.
          NNB2 = LISTC(1) + 1
      END IF
*
*       Predict coordinates of perturbers & possibly c.m. to order FDOT.
      DO 1 L = 2,NNB2
          J = LISTC(L)
          call jpred(j,time,time)
c$$$          S = TIME - T0(J)
c$$$*       Do not allow prediction outside range (NB! No bad effects in DIFSY1).
c$$$*         S = MIN(S,STEP(J))
c$$$          X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + X0DOT(1,J))*S + X0(1,J)
c$$$          X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + X0DOT(2,J))*S + X0(2,J)
c$$$          X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + X0DOT(3,J))*S + X0(3,J)
*       Note most recent velocities are used for perturber integration.
    1 CONTINUE
*
*       Obtain global coordinates & velocities from current chain & c.m.
   4  LK = 0
      call jpred(ich,time,time)
      DO 10 I = 1,NN
          DO 5 K = 1,3
              LK = LK + 1
              XC(K,I) = XCH(LK) + X(K,ICH)
              UC(K,I) = VCH(LK) + XDOT(K,ICH)
    5     CONTINUE
*     --03/14/14 15:45-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*          if(I.eq.1605) then
c$$$             write(120+rank,*)'xcpred i',i,'xc', xc(1,i),'xch',xch(lk),
c$$$     &            'uc',uc(1,i),'x',x(1,ich),'t',time,'icall',icallp
c$$$             call flush(120+rank)
c$$$             call mpi_barrier(MPI_COMM_WORLD,ierr)
*             if(icallp.eq.13548) call abort()
*          end if
*     --03/14/14 15:45-lwang-end----------------------------------------*
          
   10 CONTINUE
*
      RETURN
*
      END
