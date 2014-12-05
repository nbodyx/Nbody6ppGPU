      SUBROUTINE KSPREG
*     
*     
*     New KS regularization.
*     ----------------------
*     
      INCLUDE 'common6.h'
      REAL*8 DSAVE(17)
      REAL*8  Q(3),RDOT(3),UI(4),A1(3,4)
*     EXTERNAL RENAME
*     
*     
*     Save basic variables for components unless in correct location.
      DO 10 KCOMP = 1,2
*     Treat the first & second component in turn.
         IF (KCOMP.EQ.1) THEN
            I = ICOMP
         ELSE
            I = JCOMP
         END IF
         J = 2*NPAIRS + KCOMP
         IF (I.EQ.J) GO TO 10
*     
         DO 2 K = 1,3
            DSAVE(K) = X(K,I)
            DSAVE(K+3) = XDOT(K,I)
 2       CONTINUE
         DSAVE(7) = BODY(I)
         DSAVE(8) = RS(I)
         DSAVE(9) = RADIUS(I)
         DSAVE(10) = TEV(I)
         DSAVE(11) = BODY0(I)
         DSAVE(12) = TEV0(I)
         DSAVE(13) = EPOCH(I)
         DSAVE(14) = SPIN(I)
         DSAVE(15) = ZLMSTY(I)
         DSAVE(17) = STEP(I)
         NAMEI = NAME(I)
         KSI = KSTAR(I)
*     
*     Exchange first & second single particle with ICOMP & JCOMP.
         DO 4 K = 1,3
            X(K,I) = X(K,J)
            X0(K,I) = X0(K,J)
            X0DOT(K,I) = X0DOT(K,J)
            XDOT(K,I) = XDOT(K,J)
            F(K,I) = F(K,J)
            FDOT(K,I) = FDOT(K,J)
            FI(K,I) = FI(K,J)
            FIDOT(K,I) = FIDOT(K,J)
            D0(K,I) = D0(K,J)
            D1(K,I) = D1(K,J)
            D2(K,I) = D2(K,J)
            D3(K,I) = D3(K,J)
            FR(K,I) = FR(K,J)
            FRDOT(K,I) = FRDOT(K,J)
            D0R(K,I) = D0R(K,J)
            D1R(K,I) = D1R(K,J)
            D2R(K,I) = D2R(K,J)
            D3R(K,I) = D3R(K,J)
            X(K,J) = DSAVE(K)
            X0(K,J) = X(K,J)
            XDOT(K,J) = DSAVE(K+3)
            X0DOT(K,J) = XDOT(K,J)
 4       CONTINUE
*     
         BODY(I) = BODY(J)
         RS(I) = RS(J)
         RADIUS(I) = RADIUS(J)
         TEV(I) = TEV(J)
         TEV0(I) = TEV0(J)
         BODY0(I) = BODY0(J)
         EPOCH(I) = EPOCH(J)
         SPIN(I) = SPIN(J)
         ZLMSTY(I) = ZLMSTY(J)
         NAME(I) = NAME(J)
         KSTAR(I) = KSTAR(J)
         STEP(I) = STEP(J)
         STEPR(I) = STEPR(J)
         T0(I) = T0(J)
         TIMENW(I) = T0(I) + STEP(I)
         T0R(I) = T0R(J)
         K = LIST(1,J) + 1
         DO 5 L = 1,K
            LIST(L,I) = LIST(L,J)
 5       CONTINUE
         BODY(J) = DSAVE(7)
         RS(J) = DSAVE(8)
         RADIUS(J) = DSAVE(9)
         TEV(J) = DSAVE(10)
         BODY0(J) = DSAVE(11)
         TEV0(J) = DSAVE(12)
         EPOCH(J) = DSAVE(13)
         SPIN(J) = DSAVE(14)
         ZLMSTY(J) = DSAVE(15)
         NAME(J) = NAMEI
         KSTAR(J) = KSI
         STEP(J) = DSAVE(17)

 10   CONTINUE

*     Increase pair index, total number & single particle index.
      NPAIRS = NPAIRS + 1
      NTOT = N + NPAIRS
      IFIRST = 2*NPAIRS + 1
*     
*     Set new global indices of the components and current pair index.
      ICOMP = 2*NPAIRS - 1
      JCOMP = ICOMP + 1
      IPAIR = NPAIRS
*     
*     Specify mass, neighbour radius, name, radius & type for new c.m.
      BODY(NTOT) = BODY(ICOMP) + BODY(JCOMP)
      RS(NTOT) = RS(ICOMP)
      NAME(NTOT) = NZERO + NAME(ICOMP)
      RADIUS(NTOT) = 0.0
      TEV(NTOT) = 1.0E+10
      TEV0(NTOT) = 1.0E+10
      BODY0(NTOT) = BODY(NTOT)
      EPOCH(NTOT) = TIME*TSTAR
      KSTAR(NTOT) = 0
*     
*     Define c.m. coordinates & velocities and set XDOT for components.
      DO K = 1,3
         X(K,NTOT) = (BODY(ICOMP)*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))/
     &        BODY(NTOT)
         X0(K,NTOT) = X(K,NTOT)
         X0DOT(K,NTOT) = (BODY(ICOMP)*X0DOT(K,ICOMP) + BODY(JCOMP)*
     &        X0DOT(K,JCOMP))/BODY(NTOT)
         XDOT(K,NTOT) = X0DOT(K,NTOT)
      END DO
*     
*     Define relative coordinates and velocities in physical units.
      DO 20 K = 1,3
         Q(K) = X(K,ICOMP) - X(K,JCOMP)
         RDOT(K) = XDOT(K,ICOMP) - XDOT(K,JCOMP)
 20   CONTINUE
*     
*     Introduce regularized variables using definition of 1985 paper.
      R(IPAIR) = SQRT(Q(1)**2 + Q(2)**2 + Q(3)**2)
*     
*     Initialize the regularized coordinates according to sign of Q(1).
      IF (Q(1).LE.0.0D0) THEN
         UI(3) = 0.0D0
         UI(2) = SQRT(0.5D0*(R(IPAIR) - Q(1)))
         UI(1) = 0.5D0*Q(2)/UI(2)
         UI(4) = 0.5D0*Q(3)/UI(2)
      ELSE
         UI(4) = 0.0D0
         UI(1) = SQRT(0.5D0*(R(IPAIR) + Q(1)))
         UI(2) = 0.5D0*Q(2)/UI(1)
         UI(3) = 0.5D0*Q(3)/UI(1)
      END IF
*     
*     Set current transformation matrix.
      CALL MATRIX(UI,A1)
*     
*     Form regularized velocity and set initial KS coordinates.
      TDOT2(IPAIR) = 0.0
      DO 30 K = 1,4
         UDOT(K,IPAIR) = 0.50D0*(A1(1,K)*RDOT(1) + A1(2,K)*RDOT(2) +
     &        A1(3,K)*RDOT(3))
*     Note that A1(J,K) is the transpose of A1(K,J).
         U(K,IPAIR) = UI(K)
         U0(K,IPAIR) = U(K,IPAIR)
         TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*UI(K)*UDOT(K,IPAIR)
 30   CONTINUE
*     
*     Evaluate initial binding energy per unit mass.
      H(IPAIR) = (2.0D0*(UDOT(1,IPAIR)**2 + UDOT(2,IPAIR)**2 +
     &     UDOT(3,IPAIR)**2 + UDOT(4,IPAIR)**2) -
     &     BODY(NTOT))/R(IPAIR)
      EB = H(IPAIR)*BODY(ICOMP)*BODY(JCOMP)/BODY(NTOT)
*     
*     Specify zero membership and large step for second component).
      LIST(1,JCOMP) = 0
*     Set large step for second component to avoid detection.
      STEP(JCOMP) = 1.0E+06
      STEPR(JCOMP) = 1.0E+06
*     
*     Obtain apocentre distance.
      SEMI = -0.5*BODY(NTOT)/H(IPAIR)
      ECC2 = (1.0-R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(NTOT)*SEMI)
      RAP = SEMI*(1.0 + SQRT(ECC2))
*     
*     Include suggestion for monitoring hyperbolic encounters (suppressed).
*     IF (SEMI.LT.0.0) THEN
*     PMIN = SEMI*(1.0 - SQRT(ECC2))
*     if(rank.eq.0)
*     &    WRITE (6,56)  TIME+TOFF, NAME(ICOMP), NAME(JCOMP), PMIN
*     56     FORMAT (' HYPERBOLIC    T NAM PMIN ',F8.2,2I6,1P,E10.2)
*     END IF
*     
*     Set 2*SEMI as termination scale for hard binary if 2*SEMI < RS/20.
      IF (EB.LT.EBH.AND.RAP.LT.0.05*RS(NTOT)) THEN 
         R0(IPAIR) = MAX(RMIN,-BODY(NTOT)/H(IPAIR))
      ELSE
         R0(IPAIR) = R(IPAIR)
      END IF
*     
*     Increase regularization counters (#9 & NKSHYP for hyperbolic orbits).
      NKSREG = NKSREG + 1
      IF (H(IPAIR).GT.0.0) NKSHYP = NKSHYP + 1
*     
*     Check optional output for new KS.
      IF (KZ(10).GT.0) THEN
         RI = SQRT((X(1,NTOT) - RDENS(1))**2 +
     &        (X(2,NTOT) - RDENS(2))**2 +
     &        (X(3,NTOT) - RDENS(3))**2)
         if(rank.eq.0)
     &        WRITE (6,60)  NAME(ICOMP), NAME(JCOMP), 
     &        BODY(ICOMP)*ZMBAR,BODY(JCOMP)*ZMBAR,DTAU(IPAIR),
     &        R(IPAIR), R(IPAIR)*RAU, RI, H(IPAIR), IPAIR, NAME(N+IPAIR)
 60      FORMAT (' NEW KSREG   NAME(I1)',
     &        I10,' NAME(I2)',I10,' M(I1)[M*]',
     &        F9.3,' M(I2)[M*]',F9.3,' DTAU',F8.3,' R12[NB,AU]',1P,2E9.1
     &        ,0P,' RI[NB]',F7.2,' H',F9.2,' IPAIR',I9,' NAME(ICM)',I10)
         call flush(6)
      END IF
*     
*     Modify the termination criterion according to value of NPAIRS.
      IF (NPAIRS.GT.KMAX - 3) GMAX = 0.8*GMAX
      IF (NPAIRS.LT.KMAX - 5.AND.GMAX.LT.0.001) GMAX = 1.2*GMAX
      IF (NPAIRS.EQ.KMAX.and.rank.eq.0) WRITE (6,70)  NPAIRS, TIME+TOFF
 70   FORMAT (5X,'WARNING!   MAXIMUM KS PAIRS   NPAIRS TIME',I5,F9.2)
*     
*     Initialize prediction variables to avoid skipping KS components.
      DO 75 KCOMP = 1,2
         JDUM = 2*NPAIRS - 2 + KCOMP
         DO 72 K = 1,3
            X0(K,JDUM) = X(K,JDUM)
            X0DOT(K,JDUM) = 0.0D0
            F(K,JDUM) = 0.0D0
            FDOT(K,JDUM) = 0.0D0
            D2(K,JDUM) = 0.0D0
            D3(K,JDUM) = 0.0D0
            D2R(K,JDUM) = 0.0D0
            D3R(K,JDUM) = 0.0D0
 72      CONTINUE
 75   CONTINUE
*     
*     Set flag to distinguish between existing and new binaries (K = -1/0).
      KSLOW(IPAIR) = 1
*     
      RETURN
*     
      END
