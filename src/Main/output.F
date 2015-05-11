
      SUBROUTINE OUTPUT
*
*
*       Output and data save.
*       ---------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
#ifdef TT
      INCLUDE 'tt.h'
#endif
      INCLUDE 'galaxy.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      COMMON/ECHAIN/  ECH
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
*      COMMON/NNBCOT/  nnbave
      REAL*8  X1(3,4),V1(3,4),UI(4),VI(4),XREL2(3),VREL2(3)
      REAL*4  XS(3,NMAX),VS(3,NMAX),BODYS(NMAX),RHOS(NMAX),AS(20)
      REAL*4  XJ(3,6),VJ(3,6),BODYJ(6)
      REAL*4  XNS(NMAX),PHI(NMAX)
      CHARACTER*27 OUTFILE
      CHARACTER*27 OUTTAIL
      CHARACTER*20 TCHAR
*
*     Call Computation of Moments of Inertia (with Chr. Theis)
      IF(KZ(49).GT.0) CALL ELLAN

*     Obtain energy error in case routine ADJUST not called recently.
      IF (TIME.GE.TADJ.OR.TIME.LE.0.0D0) GO TO 10
*
*       Predict X & XDOT for all particles (except unperturbed pairs).
      CALL XVPRED(IFIRST,NTOT)
*
*       Obtain the total energy at current time (resolve all KS pairs).
#ifdef PARALLEL
      CALL ENERGY_MPI(.false.)
#else
      CALL ENERGY(.false.)
#endif
*
*       Include KS pairs, triple, quad, mergers, collisions & chain.
      ETOT = ZKIN - POT + ETIDE + EBIN + ESUB + EMERGE + ECOLL + EMDOT
     &                                                         + ECDOT
      IF (NCH.GT.0) THEN
          ETOT = ETOT + ECH
      END IF
*
*       Update energies and form the relative error (divide by ZKIN or ETOT).
      BE(2) = BE(3)
      BE(3) = ETOT
      DE = BE(3) - BE(2)
      DETOT = DETOT + DE
      DE = DE/MAX(ZKIN,ABS(ETOT))
*       Save sum of relative energy error for main output and accumulate DE.
      ERROR = ERROR + DE
*
*       Find density centre & core radius (Casertano & Hut, Ap.J. 298, 80).
      IF (N.GE.20.AND.KZ(29).EQ.0) THEN
          CALL CORE
      END IF
*
*       Check optional sorting of Lagrangian radii & half-mass radius.
      IF (KZ(7).GT.0) THEN
          CALL LAGR(RDENS)
      END IF
*
*       Initialize diagnostic variables.
   10 NP = 0
      IUNP = 0
      AMIN = 100.0
      MULT = 0
*
*       Find smallest semi-major axis and count unperturbed KS pairs.
      DO 20 IPAIR = 1,NPAIRS
          NP = NP + LIST(1,2*IPAIR-1)
          SEMI = -0.5*BODY(N+IPAIR)/H(IPAIR)
          IF (SEMI.GT.0.0) AMIN = MIN(AMIN,SEMI)
          IF (LIST(1,2*IPAIR-1).EQ.0) IUNP = IUNP + 1
          IF (NAME(N+IPAIR).LT.-2*NZERO) MULT = MULT + 1
   20 CONTINUE
*
*       Include search of any hierarchical triples.
      DO 25 IM = 1,NMERGE
          ZMB = CM(1,IM) + CM(2,IM)
          SEMI = -0.5*ZMB/HM(IM)
          AMIN = MIN(AMIN,SEMI)
   25 CONTINUE
*
*       Perform time-step & neighbour statistics (NS is # single stars).
      DTI = 0.0
      DTRI = 0.0
      CNNB = 0.0
      CMAX = 0.0
      NNB = 0
      NS = 0
      SUM = 0.0
      DO 30 I = IFIRST,NTOT
          DTI = DTI + 1.0/STEP(I)
          DTRI = DTRI + 1.0/STEPR(I)
          CNNB = CNNB + LIST(1,I)/STEP(I)
          RHON = LIST(1,I)/RS(I)**3
          CMAX = MAX(CMAX,RHON)
          NNB = NNB + LIST(1,I)
          IF (I.LE.N.AND.BODY(I).GT.0.0D0) NS = NS + 1
          SUM = SUM + BODY(I)**2
   30 CONTINUE
      NS = NS - NSUB
*
*       Estimate relative cost & effective neighbour number of AC scheme.
      COST = CNNB/(FLOAT(N - NPAIRS)*DTRI)
      CNNB = CNNB/DTI
*       Scale maximum particle density contrast by the mean value.
      CMAX = 2.0*CMAX*RSCALE**3/FLOAT(N)
*
*       Set average neighbour number & density centre displacement.
      NNB = INT(FLOAT(NNB)/FLOAT(N - NPAIRS))
*      nnbave = nnb
      RD = SQRT(RDENS(1)**2 + RDENS(2)**2 + RDENS(3)**2)
*
*       Check print frequency indicator & optional model counter.
      NPRINT = NPRINT + 1
      IF (NPRINT.GT.NFIX.OR.TIME.LE.0.0D0) THEN
          NPRINT = 1
          IF (KZ(3).GT.0) MODEL = MODEL + 1
      END IF
*
*       Form binary & merger energy ratios.
      EB = EBIN/(ZKIN - POT)
      EM = EMERGE/(ZKIN - POT)
      IF (KZ(21).GT.1) THEN
          CALL JACOBI(NESC)
      ELSE
          NESC = 0
      END IF
*
*       Print main output diagnostics.
      I6 = INT(TSCALE*TTOT)
*
      IF(rank.EQ.0)THEN
      WRITE (6,40) TTOT, N, NNB, NPAIRS, NMERGE, MULT, NS, NSTEPI,
     &         NSTEPB,NSTEPR, NSTEPU, ERROR, BE(3)
   40 FORMAT (//,' TIME[NB] ',1P,E15.6,0P,'  N ',I8,'  <NB> ',I5,
     &     '  NPAIRS ',I7,'  NMERGE ',I5,' MULT ',I4,' NS ',I8,
     &     '  NSTEP(I,B,R,U) ',4I13,'  DE ',E15.6,'  E ',F17.6)
      END IF
*
      IF (KZ(21).GT.0) THEN
          CALL CPUTIM(ttotb)
          TCOMP = ttotb-ttota
          IF (VC.EQ.0.0D0) VC = RSCALE/TCR
          TRC = 1.02*FLOAT(NC)**2*BODYM/(VC**3*LOG(FLOAT(NC)))
          DMIN1 = MIN(DMIN1, DMIN2, DMIN3, DMIN4, DMINC)
          NEFF = INT(ZMASS**2/SUM)
*
          IF(rank.EQ.0)THEN
             WRITE (6,45)  NRUN, MODEL, TCOMP, TRC, DMIN1, DMIN2, DMIN3,
     &            DMIN4, AMIN, RMAX, RSMIN, NEFF
 45          FORMAT (/,' NRUN =',I3,'  M# =',I3,'  CPU =',1P,E12.5,0P,
     &            '  TRC =',F5.1, '  DMIN =',1P,4E8.1,'  AMIN =',E8.1,
     &            '  RMAX =',E8.1,'  RSMIN =',0P,F5.2,'  NEFF =',I6)
          END IF
      END IF
      VRMS = SQRT(0.5*ZMASS/RSCALE)*VSTAR
*
      IF(rank.EQ.0)THEN
      WRITE (6,50)
   50 FORMAT (/,'    <R>  RTIDE  RDENS   RC      NC   MC   RHOD   RHOM',
     &          '  CMAX   <Cn>  Ir/R      UN      NP    RCM    VCM',
     &          '         AZ     EB/E   EM/E   TCR     T6  NESC',
     &                   '  VRMS')
*
      WRITE (6,55)  RSCALE, RTIDE, RD, RC, NC, ZMC, RHOD, RHOM, CMAX,
     &              CNNB, COST, IUNP, NP, CMR(4), CMRDOT(4), AZ, EB, EM,
     &              TCR, I6, NESC, VRMS
   55 FORMAT (' #1',F5.2,F6.1,F7.2,F6.2,I7,F7.3,F6.0,F7.0,F6.0,F6.1,
     &                   F6.2,2I8,F9.3,F8.4,F11.6,2F7.3,F6.2,2I6,F6.1)
*
      WRITE (6,60)
   60 FORMAT (/,7X,'   NIRRF      NBCORR  NNPRED  NBFULL  NBVOID ',
     &     '  NICONV  NLSMIN NBSMIN   NBDIS  NBDIS2  NCMDER  ',
     &     'NFAST NBFAST      NBLOCK      NBLCKR    NBFLUX')
      WRITE (6,65)  NIRRF, NBCORR, NNPRED, NBFULL, NBVOID, NICONV,
     &              NLSMIN, NBSMIN, NBDIS, NBDIS2, NCMDER,  NFAST,
     &              NBFAST, NBLOCK, NBLCKR,  NBFLUX
*
   65 FORMAT (' #2',2I12,9I8,2I7,2I12,I10)
*
      WRITE (6,70)
   70 FORMAT (/,5X,'   NKSTRY     NKSREG     NKSHYP     NKSPER  ',
     &     'NKSMOD   NTTRY   NTRIP    NQUAD  NCHAIN   ',
     &     '  NMERG   NEWHI  NSTEPT  NSTEPQ  NSTEPC')
      WRITE (6,75)  NKSTRY, NKSREG,  NKSHYP, NKSPER,
     &              NKSMOD, NTTRY, NTRIP, NQUAD, NCHAIN, NMERG, NEWHI,
     &              NSTEPT, NSTEPQ, NSTEPC
   75 FORMAT (' #3',3I11,I11,3I8,3I9,I8,I8,3I8,2I10)
      END IF

*     Output global features and counters into global.30
      IF(rank.eq.0) then
         IF(kstart.eq.1.and.ttot.eq.0.0) then
         write (30,76) 
 76      format('TIME[NB} TIME[Myr] TCR[Myr] DE BE(3) ',
     &           'RSCALE[PC] RTIDE[PC] RDENS[PC] RC[PC]  RHOD[PC] ',
     &           'RHOM[PC] MC[M*] CMAX <Cn> Ir/R RCM VCM AZ ',
     &           'EB/E EM/E VRMS ',
     &           'N NS NPAIRS NUPKS NPKS NMERGE MULT <NB> NC NESC ',
     &           'NSTEPI NSTEPB NSTEPR NSTEPU NSTEPT NSTEPQ NSTEPC ',
     &           'NBLOCK NBLCKR NNPRED ',
     &           'NIRRF NBCORR NBFLUX NBFULL NBVOID NICONV NLSMIN ',
     &           'NBSMIN NBDIS NBDIS2 NCMDER ',
     &           'NFAST NBFAST ',
     &           'NKSTRY NKSREG NKSHYP NKSPER NKSMOD',
     &           'NTTRY NTRIP NQUAD NCHAIN NMERG NEWHI')
         END IF
         write (30,77) TTOT,TTOT*TSTAR,TCR*TSTAR,ERROR,BE(3),
     &        RSCALE*RBAR,RTIDE*RBAR,RD*RBAR,RC*RBAR,RHOD*RBAR,
     &        RHOM*RBAR,ZMC*ZMBAR,CMAX,CNNB,COST,CMR(4),CMRDOT(4),AZ,
     &        EB,EM,VRMS,
     &        N,NS,NPAIRS,IUNP,NP,NMERGE,MULT,NNB,NC,NESC,
     &        NSTEPI,NSTEPB,NSTEPR,NSTEPU,NSTEPT,NSTEPQ,NSTEPC,
     &        NBLOCK,NBLCKR,NNPRED,
     &        NIRRF,NBCORR,NBFLUX,NBFULL,NBVOID,NICONV,NLSMIN,NBSMIN,
     &        NBDIS,NBDIS2,NCMDER,
     &        NFAST,NBFAST,
     &        NKSTRY,NKSREG,NKSHYP,NKSPER,NKSMOD,
     &        NTTRY,NTRIP,NQUAD,NCHAIN,NMERG,NEWHI
 77      format(21E26.17,44I12)
         call flush(30)
      END IF
*
*       Check output for mass loss or tidal capture.
      IF (KZ(19).GT.0.OR.KZ(27).GT.0) THEN
          CALL EVENTS
      END IF

*     Reset counter to zero
      call counter_reset

*     Obtain half-mass radii for two groups (NAME <= NZERO/5 & > NZERO/5). 
      IF (KZ(7).GE.6) THEN
          CALL LAGR2(RDENS)
      END IF
*
*       Include diagnostics about cluster orbit in general external field.
#ifdef TT
*** FlorentR - add KZ(14)=9 in case of mode B
      IF (KZ(14).EQ.3.OR.(KZ(14).EQ.9.AND.TTMODE.EQ.0)) THEN
*** FRenaud
#else
      IF (KZ(14).EQ.3) THEN
#endif
          GZ = RG(1)*VG(2) - RG(2)*VG(1)
          SX = RBAR/1000.0
          WRITE (6,78)  (RG(K)*SX,K=1,3), (VG(K)*VSTAR,K=1,3),
     &                  GZ, ETIDE
   78     FORMAT (/,5X,'CLUSTER ORBIT    RG VG JZ ET ',
     &                                 3F7.2,2X,3F7.1,1P,E16.8,E10.2)
      END IF
      IF (KZ(14).EQ.4) THEN
          WRITE (6,79)  TTOT, N, RSCALE, ZMASS, MP, DETOT
   79     FORMAT (/,5X,'GAS EXPULSION    T N <R> M MP DETOT ',
     &                                   F7.1,I7,3F7.3,1P,E10.2)
      END IF
*
*       Reset minimum encounter distances & maximum apocentre separation.
      DMIN2 = 100.0
      DMIN3 = 100.0
      DMIN4 = 100.0
      DMINC = 100.0
      RSMIN = 100.0
      RMAX = 0.0
*
*       Check integer overflows (2^{32} or 2.1 billion).
C      IF (NSTEPI.GT.2000000000.OR.NSTEPI.LT.0) THEN
C          NSTEPI = 0
C          NIRECT = NIRECT + 1
C      END IF
C      IF (NSTEPR.GT.2000000000.OR.NSTEPR.LT.0) THEN
C          NSTEPR = 0
C          NRRECT = NRRECT + 1
C      END IF
C      IF (NSTEPU.GT.2000000000.OR.NSTEPU.LT.0) THEN
C          NSTEPU = 0
C          NURECT = NURECT + 1
C      END IF
C      IF (NBPRED.GT.2000000000.OR.NBPRED.LT.0) THEN
C          NBPRED = 0
C      END IF
C      IF (NBFLUX.GT.2000000000.OR.NBFLUX.LT.0) THEN
C          NBFLUX = 0
C          NBRECT = NBRECT + 1
C      END IF
C      IF (NBCORR.GT.2000000000.OR.NBCORR.LT.0) THEN
C          NBCORR = 0
C      END IF
*     IF (NBLOCK.GT.2000000000.OR.NBLOCK.LT.0) THEN
*         NBLOCK = 0
*     END IF
*
*       Exit if error exceeds restart tolerance (TIME < TADJ means no CHECK).
      IF (ABS(ERROR).GT.5.0*QE.AND.TIME.LT.TADJ) GO TO 100
*
*       Check optional analysis & output of KS binaries.
      IF (KZ(8).GT.0.AND.NPAIRS.GT.0) THEN
          CALL BINOUT
      END IF
*
*       Include optional diagnostics of block-steps.
      IF (KZ(33).GT.0) THEN
          CALL LEVELS
      END IF
*
*       Check optional output of single bodies & binaries.
      IF (KZ(6).GT.0) THEN
          CALL BODIES
      END IF
*
*       See whether to write data bank of binary diagnostics on unit 9.
      IF (KZ(8).GE.2) THEN
          CALL BINDAT
          IF (KZ(18).GE.4) THEN
              CALL HIDAT
          END IF
      END IF
*
*       Check optional diagnostics of evolving stars.
      IF (KZ(12).GT.0.AND.TIME.GE.TPLOT) THEN
*       Update stellar evolution
c$$$              IF (KZ(19).GE.3) THEN
c$$$                  CALL MDOT
c$$$              ELSE
c$$$                  CALL MLOSS
c$$$              END IF
          CALL HRPLOT
*       Update next output time.
          TPLOT = TPLOT + DTPLOT
      END IF
*
*       Check optional writing of data on unit 3 (frequency NFIX).
      IF (KZ(3).EQ.0.OR.NPRINT.NE.1) GO TO 100
*
      DO 80 IPAIR = 1,NPAIRS
          RHO(2*IPAIR-1) = RHO(N+IPAIR)
          RHO(2*IPAIR) = RHO(N+IPAIR)
*
*       Include prediction of unperturbed binaries (except ghosts).
          IF (LIST(1,2*IPAIR-1).EQ.0.AND.BODY(N+IPAIR).GT.0.0D0) THEN
             CALL RESOLV(IPAIR,1)
          END IF
 80    CONTINUE

*       Convert masses, coordinates & velocities to single precision.
      DO 90 I = 1,NTOT
          BODYS(I) = REAL(BODY(I))
          RHOS(I) = REAL(RHO(I))
          XNS(I) = REAL(XNDBL(I))
          PHI(I) = REAL(PHIDBL(I))
          DO 85 K = 1,3
              XS(K,I) = REAL(X(K,I))
              VS(K,I) = REAL(XDOT(K,I))
   85     CONTINUE
   90 CONTINUE
*
*       Replace any ghosts by actual M, R & V (including 2 binaries).
      DO 95 JPAIR = 1,NPAIRS
         J2 = 2*JPAIR
         J1 = J2 - 1
         ICM = N + JPAIR
*     Determine merger & ghost index for negative c.m. name.
         IF (NAME(ICM).LT.0.AND.BODY(ICM).GT.0.0) THEN
*     Include possible quartet [[B,S],S] and quintet [[B,S],B] first.
            IF (NAME(ICM).LT.-2*NZERO) THEN
*     Find ghost and merger index at previous hierarchical level.
               IM = 1
               DO 92 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAME(ICM) + 2*NZERO) IM = K
 92            CONTINUE
               JG = N
               DO 94 K = IFIRST,NTOT
                  IF (NAMEG(IM).EQ.NAME(K)) JG = K
 94            CONTINUE
*     Determine the current ghost and merger index.
               CALL FINDJ(J1,JG2,IM2)
*     Distinguish netween quartet and quintet.
               IF (NAME(J2).LE.NZERO) THEN
                  IF (JG.LE.N) THEN
                     BODYS(JG) = REAL(CM(2,IM))
                     BODYS(J1) = REAL(CM(1,IM2))
                     BODYS(JG2) = REAL(CM(2,IM2))
                  ELSE
                     JP = JG - N
                     BODYS(2*JP-1) = REAL(CM(3,IM))
                     BODYS(2*JP) = REAL(CM(4,IM))
                     BODYS(JG2) = REAL(CM(2,IM2))
                  END IF
               ELSE
                  IF (JG.LE.N) THEN
                     BODYS(JG) = REAL(CM(2,IM))
                     BODYS(JG2) = REAL(CM(2,IM2))
                  ELSE
                     JP = JG - N
                     BODYS(2*JP-1) = REAL(CM(3,IM2))
                     BODYS(2*JP) = REAL(CM(4,IM2))
                  END IF
               END IF
               GO TO 95
            END IF
*     
            CALL FINDJ(J1,J,IM)
*     Note: J is ghost index and IM is merger index.
            IF (J.LE.0) GO TO 95
            BODYS(J1) = REAL(CM(1,IM))
            BODYS(J) = REAL(CM(2,IM))
            ZMB = REAL(CM(1,IM) + CM(2,IM))
*     Form global coordinates and velocities from c.m. with XREL & VREL.
            DO K = 1,3
               X1(K,1) = X(K,J1) + CM(2,IM)*XREL(K,IM)/ZMB
               X1(K,2) = X(K,J1) - CM(1,IM)*XREL(K,IM)/ZMB
               V1(K,1) = XDOT(K,J1) + CM(2,IM)*VREL(K,IM)/ZMB
               V1(K,2) = XDOT(K,J1) - CM(1,IM)*VREL(K,IM)/ZMB
*     
               XS(K,J1) = REAL(X1(K,1))
               XS(K,J)  = REAL(X1(K,2))
               VS(K,J1) = REAL(V1(K,1))
               VS(K,J)  = REAL(V1(K,2))
            END DO
*     Look for ghosts of possible second (i.e. outer) merged binary.
c$$$            if(rank.eq.0)print*,' after findj t,j1,j2,jpair,j,im=',j1,
c$$$     &           j2,jpair,j,im
c$$$            if(rank.eq.0)print*,' after findj nmerge names=',
c$$$     &           nmerge,name(j1),name(j2),name(icm)
*     See whether ghost is second merged binary (QUAD) instead of triple.
            IF (J.GT.N) THEN
               ICM2 = J
*     Save outer binary masses as well as second KS component & ghost c.m.
               IPAIR = ICM2 - N
               I1 = 2*IPAIR - 1
               I2 = I1 + 1
               BODYS(I1) = REAL(CM(3,IM))
               BODYS(I2) = REAL(CM(4,IM))
               BODYS(J2) = REAL(CM(2,IM))
               BODYS(J) = REAL(CM(3,IM) + CM(4,IM))
*     Copy KS variables to local scalars.
               DO K = 1,4
                  UI(K) = U(K,IPAIR)
                  VI(K) = UDOT(K,IPAIR)
               END DO
*     Transform to physical variables and multiply by 4 (momentum formula).
               CALL KSPHYS(UI,VI,XREL2,VREL2)
               ZM = CM(3,IM) + CM(4,IM)
               DO K = 1,3
                  VREL2(K) = 4.0*VREL2(K)
                  X1(K,3) = X(K,J2) + CM(4,IM)*XREL2(K)/ZM
                  X1(K,4) = X(K,J2) - CM(3,IM)*XREL2(K)/ZM
                  V1(K,3) = XDOT(K,J2) + CM(4,IM)*VREL2(K)/ZM
                  V1(K,4) = XDOT(K,J2) - CM(3,IM)*VREL2(K)/ZM
*     
                  XS(K,I1) = REAL(X1(K,3))
                  XS(K,I2)  = REAL(X1(K,4))
                  VS(K,I1) = REAL(V1(K,3))
                  VS(K,I2)  = REAL(V1(K,4))
                  XS(K,ICM2) = REAL(X(K,J2))
                  VS(K,ICM2) = REAL(XDOT(K,J2))
               END DO
            END IF
         END IF
 95   CONTINUE
*
*       Check modification for chain regularization (case NAME(ICM) = 0).
      IF (NCH.GT.0) THEN
          CALL CHDATA(XJ,VJ,BODYJ)
          DO 98 L = 1,NCH
*       Copy global address from common JLIST (set in CHDATA).
              J = JLIST(L)
              BODYS(J) = REAL(BODYJ(L))
              DO 97 K = 1,3
                  XS(K,J) = REAL(XJ(K,L))
                  VS(K,J) = REAL(VJ(K,L))
   97         CONTINUE
   98     CONTINUE
      END IF
*
          AS(1) = REAL(TTOT)
          AS(2) = REAL(FLOAT(NPAIRS))
          AS(3) = REAL(RBAR)
          AS(4) = REAL(ZMBAR)
          AS(5) = REAL(RTIDE)
          AS(6) = REAL(TIDAL(4))
          AS(7) = REAL(RDENS(1))
          AS(8) = REAL(RDENS(2))
          AS(9) = REAL(RDENS(3))
          AS(10) = REAL(TTOT/TCR0)
          AS(11) = REAL(TSCALE)
          AS(12) = REAL(VSTAR)
          AS(13) = REAL(RC)
          AS(14) = REAL(NC)
          AS(15) = REAL(VC)
          AS(16) = REAL(RHOM)
          AS(17) = REAL(CMAX)
          AS(18) = REAL(RSCALE)
          AS(19) = REAL(RSMIN)
          AS(20) = REAL(DMIN1)
          NK = 20
*
*        Write all data in binary format on unit 3 = conf.3.
*
      if(rank.eq.0)then
*     Split the conf.3 files by time
         call string_left(TCHAR,TTOT,DELTAT)
         write(OUTFILE,118) TCHAR
 118     format('conf.3_',A20)

         IF(KZ(46).GT.0) THEN
            if(rank.eq.0.AND.TTOT.GT.0.0D0) THEN
               call custom_update_file(TTOT,DELTAT)
            END IF
         END IF

         OPEN (UNIT=3,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=OUTFILE)
         
         WRITE (3)  NTOT, MODEL, NRUN, NK
         WRITE (3)  (AS(K),K=1,NK),
C         DO J = 1, NTOT
C            WRITE (3,*) BODYS(J)*ZMBAR,XS(1:3,J)*RBAR,VS(1:3,J)*VSTAR,
C     &           NAME(J)
C         END DO
     &         (BODYS(J),J=1,NTOT),(RHOS(J),J=1,NTOT),(XNS(J),J=1,NTOT),
     &         ((XS(K,J),K=1,3),J=1,NTOT), ((VS(K,J),K=1,3),J=1,NTOT),
     &         (PHI(J),J=1,NTOT),(NAME(J),J=1,NTOT)
         CLOSE(3)

*     Tidal tail data
         IF(KZ(23).GE.3) THEN
            write(OUTTAIL,119) TCHAR
 119        format('tail.5_',A20)
            OPEN (UNIT=5,STATUS='UNKNOWN',FILE=OUTTAIL)
            
            NK = 13
            WRITE (5) NTAIL, MODEL, NK
*     Include cluster centre just in case.
            DO K = 1,3
               AS(K) = REAL(RG(K))
               AS(K+3) = REAL(VG(K))
               AS(K+6) = REAL(RDENS(K))
            END DO
            AS(10) = REAL(TTOT)
            AS(11) = REAL(RBAR)
            AS(12) = REAL(TSTAR)
            AS(13) = REAL(VSTAR)
            NK = 13
            WRITE (5)  (AS(K),K=1,NK)
            
            IF (NTAIL.GT.0) THEN
*     Individual stars
               DO I = ITAIL0,NTTOT
                  BODYS(I) = REAL(BODY(I))
                  DO K = 1,3
                     XS(K,I) = REAL(X(K,I) - RG(K))
                     VS(K,I) = REAL(XDOT(K,I) - VG(K))
                  END DO
                  WRITE (5) NAME(I),BODYS(I),XS(1:3,I),VS(1:3,I)
               END DO
            END IF
            CLOSE(5)
         END IF

      end if
*
*       Update next output interval and initialize the corresponding error.
  100 TNEXT = TNEXT + DELTAT
      ERROR = 0.0D0
*
      RETURN
*
      END
