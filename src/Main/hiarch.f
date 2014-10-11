      SUBROUTINE HIARCH(IPAIR)
*
*
*       Hierarchical system diagnostics.
*       --------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A1(3),A2(3),XREL(3),VREL(3),EI(3),HI(3),HO(3),
     &        TK0(100),TLAST(100),TTERM(100)
      INTEGER  LAST(100),IPRINT(100)
      LOGICAL  FIRST
      SAVE  FIRST,TLAST,TTERM,RAP,TK0,LAST,IPRINT,NL
      DATA  FIRST /.TRUE./
      DATA  LAST,IPRINT  /200*0/
*
*
*       Open unit #10 the first time.
      IF (FIRST) THEN
*          OPEN (UNIT=10,STATUS='UNKNOWN',FORM='FORMATTED',FILE='HIARCH')
          FIRST = .FALSE.
*
*       Print cluster scaling parameters at start of the run.
          IF (rank.eq.0.AND.KSTART.EQ.1) THEN
              WRITE (12,1)  RBAR, BODYM*ZMBAR, BODY1*ZMBAR, TSCALE,
     &                      NBIN0, NZERO
 1            FORMAT (/,6X,1P, 'MODEL:    RBAR =',E26.17,'  <M>[M*] =',
     &            E26.17,'  M1[M*] =',E26.17,'  TSCALE =',E26.17,0P,
     &            '  NB0 =',I12,'  N0 =',I12,//)
          END IF
*
*       Define the saved variables in case of COMMON dump during merger.
          TK0(1) = 1.0
          RAP = 1.0
          NL = 0
          TLAST(1) = 0.0
          TTERM(1) = 0.0
C          if(rank.eq.0.AND.KSTART.EQ.1)then
C             WRITE (12,2)
C 2           FORMAT ('      TIME   A         A1        E1    PMIN',
C     &            '    P1/P0  Q     PCR   MB    Q1     I      NAME',
C     &            '        KC',/)
C             WRITE (12,3)
C 3           FORMAT ('             r/Rc                           ',
C     &            '           Rc   GAM   IT        N0   N1   KS  NAME',
C     &            /)
C          end if
      END IF
*
*       Set c.m. and component indices.
      I = N + IPAIR
      call jpred(I,TIME,TIME)
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      TTOT = TIME + TOFF
*
*       Skip printing termination at line #2 until after first new merger.
      IF (IPHASE.EQ.7.AND.NL.EQ.0) THEN
          GO TO 30
      END IF
*
*       Determine index of current merger (may have been set previously).
      IL = 0
      DO 4 L = 1,NL
          IF (NAME(I1).EQ.LAST(L)) THEN
              IL = L
          END IF
    4 CONTINUE
*
*       Increase membership and set identifier = 0 temporarily for new case.
      IF (IL.EQ.0) THEN
          NL = NL + 1
          IL = NL
          LAST(IL) = 0
      END IF
*
*       Decide whether new merger or termination should be printed.
      IF (IPHASE.EQ.6) THEN
          IF (NAME(I1).EQ.LAST(IL)) THEN
              TLAST(IL) = TTOT
*       Avoid repeated output of same system if terminated recently (< TCR).
              IF (TTOT - TTERM(IL).LT.TCR) THEN
                  IPRINT(IL) = 0
                  GO TO 30
              ELSE
                  IPRINT(IL) = 1
              END IF
          ELSE
*       Save identity and merge time and set print indicator for line #2.
              LAST(IL) = NAME(I1)
              TLAST(IL) = TTOT
              TTERM(IL) = 0.0D0
              IPRINT(IL) = 1
          END IF
      ELSE
          TTERM(IL) = TTOT
*       Skip output of line #2 if line #1 was not printed.
          IF (IPRINT(IL).EQ.0) THEN
              GO TO 30
          END IF
          IPRINT(IL) = 0
      END IF
*
*       Form semi-major axis and eccentricity (skip hyperbolic orbits).
      SEMI = -0.5*BODY(I)/H(IPAIR)
      IF (SEMI.LT.0.0) GO TO 30
      ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
*       Skip highly eccentric fly-by's which are not formally stable.
      IF (ECC.GT.0.99) GO TO 30
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
*
*       Decide between new and terminated merger.
      IF (IPHASE.EQ.6) THEN
*       Ensure that unperturbed binary is resolved.
         call jpred(Jcomp,time,time)
          IF (LIST(1,I1).EQ.0.OR.X(1,I1).EQ.X(1,I2)) THEN
              CALL RESOLV(IPAIR,1)
          END IF
          Q = BODY(I)/BODY(JCOMP)
          Q1 = MAX(BODY(I2)/BODY(I1),BODY(I1)/BODY(I2))
          RAP = SEMI*(1.0D0 + ECC)
          RIJ2 = 0.0
          VIJ2 = 0.0
          RDOT = 0.0
          A12 = 0.0
          A22 = 0.0
          A1A2 = 0.0
          RI2 = 0.0
          VI2 = 0.0
          RVI = 0.0
          DO 5 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,JCOMP))**2
              VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,JCOMP))**2
              RDOT = RDOT + (X(K,I) - X(K,JCOMP))*
     &                                       (XDOT(K,I) - XDOT(K,JCOMP))
              K1 = K + 1
              IF (K1.GT.3) K1 = 1
              K2 = K1 + 1
              IF (K2.GT.3) K2 = 1
              A1(K) = (X(K1,I1) - X(K1,I2))*(XDOT(K2,I1) - XDOT(K2,I2))
     &              - (X(K2,I1) - X(K2,I2))*(XDOT(K1,I1) - XDOT(K1,I2))
              A2(K) = (X(K1,JCOMP) - X(K1,I))*
     &                                     (XDOT(K2,JCOMP) - XDOT(K2,I))
     &              - (X(K2,JCOMP) - X(K2,I))*
     &                                     (XDOT(K1,JCOMP) - XDOT(K1,I))
              A12 = A12 + A1(K)**2
              A22 = A22 + A2(K)**2
              A1A2 = A1A2 + A1(K)*A2(K)
*       Form relative vectors and scalars for inner binary.
              XREL(K) = X(K,I1) - X(K,I2)
              VREL(K) = XDOT(K,I1) - XDOT(K,I2)
              RI2 = RI2 + XREL(K)**2
              VI2 = VI2 + VREL(K)**2
              RVI = RVI + XREL(K)*VREL(K)
    5     CONTINUE
*
*      Evaluate orbital parameters for outer orbit (skip hyperbolic case).
          RIJ = SQRT(RIJ2)
          ZMB = BODY(I) + BODY(JCOMP)
          SEMI1 = 2.0/RIJ - VIJ2/ZMB
          SEMI1 = 1.0/SEMI1
          IF (SEMI1.LT.0.0) GO TO 30
          ECC1 = SQRT((1.0 - RIJ/SEMI1)**2 + RDOT**2/(SEMI1*ZMB))
          PMIN = SEMI1*(1.0D0 - ECC1)/RAP
          TK1 = TWOPI*SEMI1*SQRT(SEMI1/ZMB)/TK
C          PMIN = MIN(PMIN/RAP,999.0D0)
C          TK1 = MIN(TK1/TK,99999.0D0)
*       Determine inclination (8 bins of 22.5 degrees).
          FAC = A1A2/SQRT(A12*A22)
          FAC = ACOS(FAC)
*         II = 1 + FAC*360.0/(TWOPI*22.5)
          DEG = 360.0*FAC/TWOPI
*
*       Construct Runge-Lenz vector (Heggie & Rasio, 1995, IAU174, Eq.(5)).
          EI2 = 0.0
          DO 6 K = 1,3
              EI(K) = (VI2*XREL(K) - RVI*VREL(K))/BODY(I) -
     &                                                 XREL(K)/SQRT(RI2)
              EI2 = EI2 + EI(K)**2
    6     CONTINUE
*
*       Define unit vectors for inner eccentricity and angular momenta.
          COSJ = 0.0
          SJSG = 0.0
          DO 8 K = 1,3
              EI(K) = EI(K)/SQRT(EI2)
              HI(K) = A1(K)/SQRT(A12)
              HO(K) = A2(K)/SQRT(A22)
              COSJ = COSJ + HI(K)*HO(K)
              SJSG = SJSG + EI(K)*HO(K)
    8     CONTINUE
*
*       Form the expressions A & Z.
          A = COSJ*SQRT(1.0 - EI2)
          Z = (1.0 - EI2)*(2.0 - COSJ**2) + 5.0*EI2*SJSG**2
*
*       Obtain maximum inner eccentricity (Douglas Heggie, Sept. 1995).
          Z2 = Z**2 + 25.0 + 16.0*A**4 - 10.0*Z - 20.0*A**2 - 8.0*A**2*Z
          EMAX = ONE6*(Z + 1.0 - 4.0*A**2 + SQRT(Z2))
          EMAX = SQRT(EMAX)
          KCM = KSTAR(I)
          IF (NAME(I).LT.0) KCM = -10
          NEWHI = NEWHI + 1
*
*
          RBIG = MAX(RADIUS(I1),RADIUS(I2))
          PMIN = SEMI*(1.0 - ECC)
          PMIN2 = SEMI*(1.0 - EMAX)
          IF (KZ(19).GE.3) THEN
              R1 = PMIN/RBIG
              R2 = PMIN2/RBIG
          ELSE
              R1 = 0.0
              R2 = 0.0
          END IF
          ZI = FAC
          EM = SQRT(SIN(ZI)**2 + EI2*COS(ZI)**2)
          if(rank.eq.0) then
          WRITE (12,10)  TTOT, SEMI, SEMI1, ECC1, PMIN, PMIN2, TK1, Q,
     &            PCRIT/SEMI, BODY(I)/BODYM, Q1, DEG,
     &            NAME(I1), NAME(I2), NAME(JCOMP), KCM,
     &            SQRT(EI2), EM, EMAX, KSTAR(I1), KSTAR(I2), RBIG
   10     FORMAT (/,'NEW HIAR  Time SEMI0 SMEI1 ECC1 ',
     &         'PERI0 PERI0M P1/P0 M(INCM)/M(I3) PCR/SEMI0 ',
     &         'M(INCM)/<M> MR INA NAME(I1) NAME(I2) NAME(I3) ', 
     &         'K*(INCM) ECC0 ECCMIN ECCMAX K*(I1) K*(I2) RSM ',
     &         1P,E20.11,2E14.5,0P,F6.2,1P,8E13.4,0P,4I12,3F6.2,2I12,1P,
     &         E14.5,0P)
          end if
          CALL FLUSH(12)
*       Save parameters for termination diagnostics.
          TK0(IL) = TK
*
      ELSE IF (IPHASE.EQ.7) THEN
          PMIN = SEMI*(1.0D0 - ECC)/RAP
C          PMIN = MIN(PMIN/RAP,999.0D0)
          IF (TK0(IL).LE.0.0D0) TK0(IL) = TK
          TK1 = TK/TK0(IL)
C          TK1 = MIN(TK/TK0(IL),99999.0D0)
          NK = INT((TTOT - TLAST(IL))/TK0(IL))
          NK1 = INT((TTOT - TLAST(IL))/TK)
C         NK = MIN(NK,99999999)
C          NK1 = MIN(NK1,9999)
          NK = MAX(NK,0)
          RI = 0.0
          DO 15 K = 1,3
              RI = RI + (X(K,I) - RDENS(K))**2
   15     CONTINUE
          RR = RC/RSCALE
*
*       Define type index (0, 1 or KSTAR; note: #I1 is inner c.m.).
C          IT = 0
C          IF (TIME.GT.MIN(TEV(I1),TEV(I2),TEV(I))) IT = 1
C          IF (KSTAR(I1).NE.0.AND.IT.EQ.0) IT = KSTAR(I1)
*
          if(rank.eq.0)then
          WRITE (12,20)  TTOT, SQRT(RI)/RC, SEMI, ECC, PMIN, TK1,
     &                   RR, GAMMA(IPAIR), NK, NK1, NPAIRS, NAME(I2)
   20     FORMAT ('END HIAR  Time RI/RC SEMI0 ECC0 PERI0 P0F/P0I ',
     &         'RC/RSCALE GAMMA(INCM) NKI NKF NPAIRS NAME(I2)',
     &         1P,E20.11,2E14.5,0P,F6.2,1P,4E13.4,0P,2I4,2I12)
          CALL FLUSH(12)
          end if
      END IF
*
*       Restrict memory of previous cases by removing the oldest one.
   30 IF (NL.GE.100) THEN
          LT = 1
*       Skip any current mergers (defined by TTERM = 0.0).
   32     IF (TTERM(LT).LE.0.0D0) THEN
              LT = LT + 1
              GO TO 32
          END IF
          DO 35 L = LT,NL-1
              TLAST(L) = TLAST(L+1)
              TTERM(L) = TTERM(L+1)
              TK0(L) = TK0(L+1)
              LAST(L) = LAST(L+1)
              IPRINT(L) = IPRINT(L+1)
   35     CONTINUE
          NL = NL - 1
      END IF
*
      RETURN
*
      END
