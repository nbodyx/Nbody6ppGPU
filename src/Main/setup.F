      SUBROUTINE SETUP
*
*
*       Generation of initial coordinates & velocities.
*       -----------------------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER (A0=1.0D0)
      REAL*8  RAN2,A(8)
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif
*
*       Choose between uniform density and Plummer model.
      KDUM = IDUM1
      IF (KZ(5).GT.0) GO TO 20
*
*       Set up a uniform spherical system.
      DO 10 I = 1,N
    1     A(1) = 0.0D0
          DO 2 K = 1,3
              A(K+1) = 2.0*RAN2(KDUM) - 1.0
              A(1) = A(1) + A(K+1)**2
    2     CONTINUE
          IF (A(1).GT.1.0) GO TO 1
*
    4     A(5) = 0.0D0
          DO 5 K = 1,3
              A(K+5) = 2.0*RAN2(KDUM) - 1.0
              A(5) = A(5) + A(K+5)**2
    5     CONTINUE
          IF (A(5).GT.1.0) GO TO 4
*
          DO 8 K = 1,3
*             X(K,I) = A(1)*A(K+1)
*       Density proportional to 1/R**2.
              X(K,I) = A(K+1)
*       Constant density.
              XDOT(K,I) = A(K+5)
*       Isotropic velocities (magnitude randomized; no radial dependence).
    8     CONTINUE
   10 CONTINUE
*
      GO TO 200 
*
*       Initialize centre of mass terms.
   20 DO 25 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   25 CONTINUE
*
*       Treat Plummer and lower options before Jaffe and Zhao models.
      IF (KZ(5).GE.5) GO TO 52
*
*       Generate initial conditions from Plummer model (A & A 37, 183).
      DO 40 I = 1,N
   30     A(1) = RAN2(KDUM)
          IF (A(1).LT.1.0D-10) GO TO 30
          RI = (A(1)**(-0.6666667) - 1.0)**(-0.5)
*       Reject distant particles.
          IF (RI.GT.10.0) GO TO 30
*
          A(2) = RAN2(KDUM)
          A(3) = RAN2(KDUM)
          X(3,I) = (1.0 - 2.0*A(2))*RI
          X(1,I) = SQRT(RI**2 - X(3,I)**2)*COS(TWOPI*A(3))
          X(2,I) = SQRT(RI**2 - X(3,I)**2)*SIN(TWOPI*A(3))
   32     A(4) = RAN2(KDUM)
          A(5) = RAN2(KDUM)
          A(6) = A(4)**2*(1.0 - A(4)**2)**3.5
          IF (0.1*A(5).GT.A(6)) GO TO 32
*
          A(8) = A(4)*SQRT(2.0)/(1.0 + RI**2)**0.25
          A(6) = RAN2(KDUM)
          A(7) = RAN2(KDUM)
          XDOT(3,I) = (1.0 - 2.0*A(6))*A(8)
          XDOT(1,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*COS(TWOPI*A(7))
          XDOT(2,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*SIN(TWOPI*A(7))
*
*       Accumulate c.m. terms.
          DO 35 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   35     CONTINUE
   40 CONTINUE
*
*       Scale coordinates & velocities to analytical expectation values.
      SX = 1.5*TWOPI/16.0
      SV = SQRT(ZMASS/SX)
   42 DO 50 I = 1,N
          DO 45 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
              X(K,I) = SX*X(K,I)
              XDOT(K,I) = SV*XDOT(K,I)
   45     CONTINUE
   50 CONTINUE
*
*       Terminate standard cases (#5 <=1) or Jaffee/BH (#5 >= 5 label 52).
      IF (KZ(5).LE.1.OR.KZ(5).GE.5) GO TO 200
*
*       Save membership of first system for colour plot (N2 = NZERO - N1).
   52 N1 = N
*       Check optional generation of two orbiting Plummer models.
      IF (KZ(5).EQ.2) THEN
          READ (5,*)  APO, ECC, N2, SCALE
#if MPIINIT
          CALL MPI_BCAST(APO,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ECC,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(N2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(SCALE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
          N2 = MIN(N,N2)
          SEMI = APO/(1.0 + ECC)
          SEMI = MIN(SEMI,50.0D0)
          SEMI = MAX(SEMI,2.0D0)
          ECC = MIN(ECC,0.999D0)
          ECC = MAX(ECC,0.0D0)
          ZM2 = 0.0
          KSKIP = N1/N2
          DO 54 I = 1,N2
              J = KSKIP*I
              ZM2 = ZM2 + BODY(J)
   54     CONTINUE
          FAC1 = ZM2/(ZMASS + ZM2)
          FAC2 = ZMASS/(ZMASS + ZM2)
*       Restrict volume ratio to 125 (i.e. unreasonable density contrast).
          IF (SCALE.LE.0.2D0) SCALE = 0.2
*       Increase total mass.
          ZMASS = ZMASS + ZM2
*       Set apocentre velocity for new combined mass.
          VAP = SQRT(ZMASS/SEMI)*SQRT((1.0 - ECC)/(1.0 + ECC))
          DO 55 I = 1,N
              IF (I.LE.N2) THEN
*       Copy members from first system by uniform skipping (N2 <= N1).
                  J = KSKIP*I
                  BODY(I+N) = BODY(J)
                  X(1,I+N) = SCALE*X(1,J) + FAC2*APO
                  X(2,I+N) = SCALE*X(2,J)
                  X(3,I+N) = SCALE*X(3,J)
                  XDOT(1,I+N) = XDOT(1,J)/SQRT(SCALE)
                  XDOT(2,I+N) = XDOT(2,J)/SQRT(SCALE) + FAC2*VAP
                  XDOT(3,I+N) = XDOT(3,J)/SQRT(SCALE)
              END IF
              X(1,I) = X(1,I) - FAC1*APO
              XDOT(2,I) = XDOT(2,I) - FAC1*VAP
   55     CONTINUE
      ELSE IF (KZ(5).EQ.3) THEN
*       Prepare case of accretion disk with massive perturber.
          READ (5,*)  APO, ECC, SCALE
#if MPIINIT
          CALL MPI_BCAST(APO,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ECC,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(SCALE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
          RIN = 0.5
          ROUT = 1.0
          ZMASS = 1.0
          BODY(1) = ZMASS
          DO 58 K = 1,3
              X(K,1) = 0.0
              XDOT(K,1) = 0.0
   58     CONTINUE
*       Generate a thin disk population in circular orbits.
          DO 60 I = 2,N
              BODY(I) = 1.0D-03/FLOAT(N)
              SEMI = RIN + (ROUT - RIN)*FLOAT(I)/FLOAT(N)
              VCIRC = SQRT((BODY(1) + BODY(I))/SEMI)
              PHASE = TWOPI*RAN2(KDUM)
              X(1,I) = SEMI*COS(PHASE)
              X(2,I) = SEMI*SIN(PHASE)
              X(3,I) = 0.01*(2.0*RAN2(KDUM) - 1.0)
              XDOT(1,I) = -VCIRC*SIN(PHASE)
              XDOT(2,I) = VCIRC*COS(PHASE)
              XDOT(3,I) = 0.01*VCIRC*(2.0*RAN2(KDUM) - 1.0)
   60     CONTINUE
*       Define membership of perturber and ensure no external tide.
          N2 = 1
          KZ(14) = 0
*       Redefine solar mass unit and astronomical length scale in AU.
          ZMBAR = 1.0
          RBAR = 1.0/2.05D+05
          BODY(N+1) = SCALE*BODY(1)
*       Set appropriate mass ratios for transforming to new c.m. frame.
          FAC1 = BODY(N+1)/(ZMASS + BODY(N+1))
          FAC2 = ZMASS/(ZMASS + BODY(N+1))
          ZMASS = ZMASS + BODY(N+1)
*       Form orbital elements for massive perturber (avoid ECC = 1).
          IF (ABS(ECC - 1.0).GT.1.0D-05) THEN
              SEMI = APO/(1.0 + ECC)
              PERI = SEMI*(1.0 - ECC)
          ELSE
              SEMI = -1.0D+05
              PERI = 2.0*APO
          END IF
          VM2 = ZMASS*(2.0/PERI - 1.0/SEMI)
          VAP2 = ZMASS*(2.0/APO - 1.0/SEMI)
*       Determine initial y-velocity from angular momentum conservation.
          VY = SQRT(VM2)*PERI/APO
          IF(rank.eq.0.and.VAP2.LT.VY**2) THEN
             write(6,*) 'Error! Bad APO ',APO, 'VAP2',VAP2,'VY**2',VY**2
             call flush(6)
             call abort()
          end if
          VX = SQRT(VAP2 - VY**2)
*       Place perturber on the Y-axis with appropriate velocities.
          X(1,N+1) = APO*FAC2
          X(2,N+1) = 0.0
          X(3,N+1) = 0.0
          XDOT(1,N+1) = -VX*FAC2
          XDOT(2,N+1) = VY*FAC2
          XDOT(3,N+1) = 0.0
*       Displace the disk members and include negative y-velocity.
          DO 70 I = 1,N
              X(1,I) = X(1,I) - FAC1*APO
              XDOT(1,I) = XDOT(1,I) + FAC1*VX
              XDOT(2,I) = XDOT(2,I) - FAC1*VY
   70     CONTINUE
      ELSE IF (KZ(5).EQ.4) THEN
*       Include two massive bodies (ECC > 1: NAME = 1 & 2 free floating).
          N2 = 0
          READ (5,*)  SEMI, ECC, ZM1, ZM2
#if MPIINIT
          CALL MPI_BCAST(APO,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ECC,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ZM1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ZM2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
          if(rank.eq.0) WRITE (6,72)  SEMI, ECC, ZM1, ZM2
   72     FORMAT (/,12X,'MASSIVE BODIES    A =',1P,E9.1,
     &            '  E =',0P,F6.2,'  M1/<M> =',F6.2,'  M2/<M> =',F6.2)
          BODY(1) = ZM1
          BODY(2) = ZM2
          IF (ECC.LT.1.0) THEN
*       Set apocentre velocity for new combined mass (using NAME = 1 & 2).
              VAP = SQRT((ZM1 + ZM2)/SEMI)*SQRT((1.0 - ECC)/(1.0 + ECC))
              FAC1 = ZM2/(ZM1 + ZM2)
              FAC2 = ZM1/(ZM1 + ZM2)
              DO 75 K = 1,3
                  X(K,1) = 0.0
                  X(K,2) = 0.0
                  XDOT(K,2) = 0.0
   75         CONTINUE
*       Initialize binary with c.m. at rest (elements change in SCALE).
              X(1,1) = -FAC1*SEMI*(1.0 + ECC)
              X(1,2) = FAC2*SEMI*(1.0 + ECC)
              XDOT(2,1) = -FAC1*VAP
              XDOT(2,2) = FAC2*VAP
          END IF
      ELSE IF (KZ(5).EQ.5) THEN
*       Generate Jaffe model with length scale A0.
          ZMH = 1.0/SQRT(FLOAT(2*N))
          DO 100 I = 1,N
   80         ZM = RAN2(KDUM)
              RI = ZM*A0/(1.0 - ZM)
*       Reject distant particles.
              IF (RI.GT.10.0.OR.RI.LT.5.0D-03) GO TO 80
              A(2) = RAN2(KDUM)
              A(3) = RAN2(KDUM)
              X(3,I) = (1.0 - 2.0*A(2))*RI
              X(1,I) = SQRT(RI**2 - X(3,I)**2)*COS(TWOPI*A(3))
              X(2,I) = SQRT(RI**2 - X(3,I)**2)*SIN(TWOPI*A(3))
              VI2 = 1.0/(RI + A0) + 2.0*ZMH/RI
              DO 90 K = 1,3
                  VK2 = VI2*RAN2(KDUM)/3.0
                  ZZ = 2.0*RAN2(KDUM) - 1.0
                  SS = 1.0
                  IF (ZZ.LT.0.0) SS = -1.0
                  XDOT(K,I) = SS*SQRT(VK2)
   90         CONTINUE
*       Accumulate c.m. terms.
              DO 95 K = 1,3
                  CMR(K) = CMR(K) + BODY(I)*X(K,I)
                  CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   95         CONTINUE
  100     CONTINUE 
          SX = 1.0
          SV = 1.0
*       Use existing procedure for c.m. correction.
          GO TO 42
      ELSE IF (KZ(5).GE.6) THEN
*       Generate dwarf galaxy model with cusp & BH (Zhao, MN 278, 488).
          Y0 = 5.0*1.25/TWOPI/SQRT(FLOAT(N))
*       Calculate RH and ZMH by setting mass inside RH to 5*N^{1/2}*<M>.
          RH = (EXP(Y0) - 1.0D0)**0.4
          ZMH = 4.0*TWOPI/5.0*LOG(1.0D0 + RH**2.5)
*       Include possibility of using actual mass (read by SCALE).
          RCUT = 0.0
          IF (KZ(24).LT.0) THEN
              READ (5,*) ZMH, RCUT
#if MPIINIT
              CALL MPI_BCAST(ZMH,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(RCUT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
              IF (ZMH.GT.0.0) THEN
                  Y0 = 1.25/TWOPI*(ZMH/FLOAT(N))
                  RH = (EXP(Y0) - 1.0D0)**0.4
              END IF
          END IF
*
*       Add experimental choice of favourite value.
          IF (KZ(28).EQ.3) THEN
              ZMH = 1.0/SQRT(FLOAT(2*N))
              Y0 = 1.25/TWOPI*(ZMH/FLOAT(N))
              RH = (EXP(Y0) - 1.0D0)**0.4
          END IF
*
*       Define BH mass inside radius of influence (unscaled units).
          REJECT = 100.0
          CUTM = 4.0*TWOPI/5.0*LOG(1.0D0 + REJECT**2.5)
*       Specify mass fraction using cut-off value and total in M_sun.
          ZMH = ZMH/(FLOAT(N)*ZMBAR)*CUTM
          Y0 = 1.25/TWOPI*ZMH
          RH = (EXP(Y0) - 1.0D0)**0.4
          R2 = 100.0
          JDUM = -1
          RANJ = GASDEV(JDUM)
          IF (KZ(5).EQ.7) THEN
              RH = 0.006
              ZMH = LOG(1.0D0 + RH**0.75)
              REJECT = 200.0
              CUTM = LOG(1.0D0 + REJECT**0.75)
              R2 = 200.0
          END IF
*
          ICUT = 0
          DO 120 I = 1,N
  110         ZM = RAN2(KDUM)
              IF (KZ(5).EQ.6) THEN
                  ZX = 1.25/TWOPI*ZM
                  ZX = ZX*CUTM
                  ZY = EXP(ZX)
                  RI = (ZY - 1.0D0)**0.4
              ELSE
                  ZX = ZM*CUTM
                  ZY = EXP(ZX)
                  RI = (ZY - 1.0D0)**(4.0/3.0)
              END IF
*       Reject large and small distances before scaling.
              IF (RI.GT.REJECT.OR.RI.LT.1.0D-04) GO TO 110
              A(2) = RAN2(KDUM)
              A(3) = RAN2(KDUM)
              X(3,I) = (1.0 - 2.0*A(2))*RI
              X(1,I) = SQRT(RI**2 - X(3,I)**2)*COS(TWOPI*A(3))
              X(2,I) = SQRT(RI**2 - X(3,I)**2)*SIN(TWOPI*A(3))
  115         IF (KZ(5).EQ.6) THEN
                  CALL VDISP(RI,R2,ZMH,VI2)
              ELSE
                  CALL VDISP2(RI,R2,ZMH,VI2)
              END IF
              SIGMA = SQRT(VI2)
              VX = SIGMA*GASDEV(JDUM)
              VY = SIGMA*GASDEV(JDUM)
              VZ = SIGMA*GASDEV(JDUM)
              IF (MAX(ABS(VX),ABS(VY),ABS(VZ)).GT.3.0*SIGMA) GO TO 115
              XDOT(1,I) = VX
              XDOT(2,I) = VY
              XDOT(3,I) = VZ
              IF (RI.GT.RCUT) GO TO 120
              SEMI = 2.0/RI - (VX**2 + VY**2 + VZ**2)/ZMH
              SEMI = 1.0/SEMI
              RRDOT = X(1,I)*VX + X(2,I)*VY + X(3,I)*VZ
              ECC2 = (1.0 - RI/SEMI)**2 + RRDOT**2/(SEMI*ZMH)
              ECC = SQRT(ECC2)
              IF (SEMI*(1.0 - ECC).LT.0.5*RCUT) THEN
                  ICUT = ICUT + 1
                  GO TO 110
              END IF
  120     CONTINUE
*
          if(rank.eq.0) WRITE (6,125)  RH, ZMH, CUTM, RCUT, ICUT
  125     FORMAT (/,12X,'BLACK HOLE    RH =',F6.3,'  MBH =',F8.4,
     &                '  CUTM =',F7.3,'  RCUT =',F9.5,'  ICUT =',I5)
*
          SX = 1.0
          SV = 1.0
*       Substitute c.m. value (hardly matters with zero velocity).
          IF (KZ(24).LT.0) BODY(1) = ZMH
          DO 130 K = 1,3
              X(K,1) = 0.0
              XDOT(K,1) = 0.0
              CMR(K) = 0.0
              CMRDOT(K) = 0.0
  130     CONTINUE
          DO 140 I = 1,N
              DO 135 K = 1,3
                  CMR(K) = CMR(K) + BODY(I)*X(K,I)
                  CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
  135         CONTINUE
  140     CONTINUE
          GO TO 42
      END IF
*
*       Specify new membership.
      N = N + N2
      NZERO = N
      NTOT = N
      IF (N.GE.NMAX-10) THEN
          if(rank.eq.0) WRITE (6,150)  N, NMAX
  150     FORMAT (' DANGER!    LIMIT EXCEEDED   N =',I6,'  NMAX =',I6)
          STOP
      END IF
*
      IF (KZ(5).EQ.2) THEN
          if(rank.eq.0) WRITE (6,155)  SEMI, ECC, N1, N2, SCALE
  155     FORMAT (/,12X,'PLUMMER BINARY    A =',F6.2,'  E =',F6.2,
     &                  '  N1 =',I6,'  N2 =',I6,'  SCALE =',F6.2)
      ELSE IF (KZ(5).EQ.3) THEN
          if(rank.eq.0) WRITE (6,160)  APO, ECC, PERI, SCALE
  160     FORMAT (/,12X,'MASSIVE PERTURBER    APO =',F6.2,'  E =',F6.2,
     &                  '  PERI =',F6.2,'  MP/M1 =',F6.2)
      END IF
*
*       Re-initialize centre of mass terms.
      DO 162 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
  162 CONTINUE
      ZMASS = 0.0
      DO 170 I = 1,N
          ZMASS = ZMASS + BODY(I)
          DO 165 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
  165     CONTINUE
  170 CONTINUE
      DO 180 I = 1,N
          DO 175 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
  175     CONTINUE
  180 CONTINUE
*
*       Save random number sequence in COMMON for future use.
  200 IDUM1 = KDUM
*
      RETURN
*
      END
      
*----------------------------------------------------------------
      REAL*8 FUNCTION GASDEV(IDUM)
C
C          GAUSSIAN DISTRIBUTION (PRESS P. 203).
C          -------------------------------------
C
      COMMON/RANDG/  GSET,ISET
CC      DATA ISET/0/
      REAL*8 RAN2,V1,V2,R,FAC,GSET
      IF (ISET.EQ.0) THEN
1       V1 = 2.*RAN2(IDUM)-1.
        V2 = 2.*RAN2(IDUM)-1.
        R = V1**2+V2**2
        IF(R.GE.1.)GO TO 1
        FAC = SQRT(-2.*LOG(R)/R)
        GSET = V1*FAC
        GASDEV = V2*FAC
        ISET = 1
      ELSE
        GASDEV = GSET
        ISET = 0
      ENDIF

      RETURN

      END
      
*----------------------------------------------------------------
      SUBROUTINE VDISP(RI,R2,ZMH,VI2)
*
*       Velocity dispersion of dwarf galaxy model.
*       ------------------------------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
*
      ONE3 = 1.0/3.0D0
      SUM = 0.0
      RHO = 1.0/(SQRT(RI)*(1.0 + RI**2*SQRT(RI)))
      NCALL = 1000 + INT(100.0/RI**2)
      DR = 2.0*(R2 - RI)/FLOAT(NCALL)
      U = RI
*       Solve by Simpson's rule (NB! DR is twice the interval).
      DO 30 J = 1,NCALL
      U1 = U
      CALL VFUNC(ZMH,U1,F1)
      U2 = U + 0.5D0*DR
      CALL VFUNC(ZMH,U2,F2)
      U3 = U + DR
      CALL VFUNC(ZMH,U3,F3)
      SUM = SUM + ONE3*(F1 + 4.0D0*F2 + F3)*DR/2.0D0
      U = U + DR
      IF (U.LT.1.0) THEN
          DR = 1.4*DR
      END IF
      IF (U.GT.100.0) GO TO 40
   30 CONTINUE
*
   40 VI2 = SUM/RHO
*
      RETURN

      END

*----------------------------------------------------------------
      SUBROUTINE VFUNC(ZMH,U,F)
      IMPLICIT REAL*8 (A-H,M,O-Z)
*
      U2 = U**2
      ZMU = 0.8*6.2831853*LOG(1.0D0 + U2*SQRT(U))
      RHO = 1.0/(SQRT(U)*(1.0 + U2*SQRT(U)))
      F = RHO*(ZMU + ZMH)/U2
*
      RETURN
      END

      SUBROUTINE VDISP2(RI,R2,ZMH,VI2)
*
*       Velocity dispersion of r^{-9/4} galaxy model.
*       ---------------------------------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
*
      ONE3 = 1.0/3.0D0
      SUM = 0.0
      RHO = 1.0/(RI**2.25*(1.0 + RI**0.75))
      NCALL = 1000 + INT(100.0/RI**2)
      DR = 2.0*(R2 - RI)/FLOAT(NCALL)
      U = RI
*       Solve by Simpson's rule (NB! DR is twice the interval).
      DO 30 J = 1,NCALL
      U1 = U
      CALL VFUNC2(ZMH,U1,F1)
      U2 = U + 0.5D0*DR
      CALL VFUNC2(ZMH,U2,F2)
      U3 = U + DR
      CALL VFUNC2(ZMH,U3,F3)
      SUM = SUM + ONE3*(F1 + 4.0D0*F2 + F3)*DR/2.0D0
      U = U + DR
      IF (U.LT.1.0) THEN
          DR = 1.4*DR
      END IF
      IF (U.GT.200.0) GO TO 40
   30 CONTINUE
*
   40 VI2 = SUM/RHO
*
      RETURN
      END

*----------------------------------------------------------------
      SUBROUTINE VFUNC2(ZMH,U,F)
      IMPLICIT REAL*8 (A-H,M,O-Z)
*
      U2 = U**2
      ZMU = LOG(1.0D0 + U**0.75)
      RHO = 1.0/(U**2.25*(1.0 + U**0.75))
      F = RHO*(ZMU + ZMH)/U2
*
      RETURN
      END
