      subroutine global_output(tsmass0,tbmass0,RC8,
     &     N_SINGLE,N_BINARY,N_MERGER,KZ7,KZ19)
*
*
*     Global parameters output
*     -----------
      implicit none
      INTEGER NMAX,KMAX,LMAX,MMAX,MLD,MLR,MLV,MCL,NCMAX,NTMAX
      INTEGER maxpe,ithread
      INTEGER NN,NN1,NN14,NN42,NN98,NTP
      include 'params.h'
      parameter(NN=9,NN1=NN-1,NN14=14*NN,NN42=42*NN,NN98=98*NN,NTP=11)
      include 'output.h'
      include 'scale_out.h'
*      include 'omp_lib.h'

*     Input arguments
      INTEGER N_SINGLE, N_BINARY,N_MERGER, KZ7, KZ19
      REAL*8 TSMASS0,TBMASS0,RC8

***   Global parameters
*     Energy
*     ETOT, EKIN, EPOT, EBIN, ETIDE, EMDOT, ECOLL, ECDOT
*     MMAX: maximum stellar mass
      REAL EPARS(12), EKIN, EPOT, EBIN
      REAL FPARS(16), SFPARS(58) , MM
      INTEGER NPARS(5), SNPARS(72)

*     Flag
      LOGICAL fshell,fse

***   Lagragian radii
*     for VLAGR and SLAGR, 7 members store different components:
*     1. Total
*     2. X
*     3. Y
*     4. Z
*     5. R
*     6. T
*     7. ROT
*     Number counts for Total Lagrangian resolved all binaries and mergers to get correct average mass
*     Number counts for binaris use center-of-mass.
*     Velocity dispersion is mass weighted.
*     Velocity dispersion of binaries (also for different stellar types) use the central mass velocity.
*     But for stellar types, the mass use only the component mass.

      REAL FLAGR(NN)
      INTEGER NLAGR(NN),NSLAGR(NN),NBLAGR(NN),NPLAGR(NN,NTP)
      REAL RLAGR(NN),RSLAGR(NN),RBLAGR(NN),RPLAGR(NN,NTP)
      REAL MLAGR(NN),MSLAGR(NN),MBLAGR(NN),MPLAGR(NN,NTP)
      REAL VLAGR(NN,7),VSLAGR(NN,7),VBLAGR(NN,7),VPLAGR(NN,7,NTP)
      REAL SLAGR(NN,7),SSLAGR(NN,7),SBLAGR(NN,7),SPLAGR(NN,7,NTP)
      REAL ELAGR(NN),ESLAGR(NN),EBLAGR(NN),EPLAGR(NN,NTP)
      REAL VTAVE(3,NN),VTBAVE(3,NN),VTSAVE(3,NN),VTPAVE(3,NN,NTP)

***   Binary fraction
*     NSBLAGR: Total number of binaries (resolved) inside different global R_lagr
*     NSPBLAGR: Total number of primordial binaries (resolved) inside different global R_lagr
*     EBINL:   Binary binding energy inside global R_lagr
*     EBINBL:  Binary binding energy inside binary R_lagr
*     EPBINL:   Primordial binary binding energy inside global R_lagr (last is core radius)
*     EPBINBL:  Primordial binary binding energy inside binary R_lagr
      REAL MSBLAGR(NN1),MSPBLAGR(NN)
      INTEGER NSBLAGR(NN1),NSPBLAGR(NN)
      REAL EBINL(NN),EBINBL(NN1), EPBINL(NN), EPBINBL(NN1)
      REAL A(NN,3)


*** COUNTER:
*     NC: Total number of stars (binaries are resolved)
*     NCB: Total number of binaries (binaries and mergers are counted once)
*     NCS: Total number of single stars

      REAL TMASS,TSMASS,TBMASS,TPMASS(NTP),HFMASS
      REAL TMASS_F,TSMASS_F,TBMASS_F
      REAL R2(NMAX),RH,RC
      REAL MCPREV, MCSPREV, MCBPREV, MCPPREV(NTP)
      REAL C1,C2,C3, CMASS, BMASS, SMASS, PMASS(NTP)
      INTEGER NPB, NP, NPS, NPSB, NPR, NPK(NTP)
      INTEGER IALIST(NMAX), NC, NCS, NCB, NCPREV, NCSPREV, NCBPREV
      INTEGER NCP(NTP), NCPPREV(NTP)
      INTEGER I,J,K,L,kk,kks,kkb,kkp(NTP),kpre

***   Storage array
      REAL  VR(NMAX), VT(3,NMAX), VROT(NMAX), EI(NMAX), EBI(NMAX)

***   Template
      REAL r(3),v(3),mi,pot,ri,rvxy,rxy2,vrot1,vrot2,n1,n2,v2
      REAL eki,ai(3)
*     SE temp
      REAL mki(3)
      INTEGER KN,KI,kwi(3)

***   Scale
      REAL*8 ebscale,ekscale,ascale

***   stellar evoltuion
*     NKW: number of stars (single, binary with only one component, binary with two components)
*     MKW: mass of stars (single, binary with only one component (only count this component mass, binary with two components)
      INTEGER KTab(-1:15)
      INTEGER NKW(3,-1:15)
      REAL    MKW(3,-1:15)


***** Initialization ----------------------      
***   Lagrangian radii fraction of total mass
*     Last is reserved for core
      DATA FLAGR/0.001,0.01,0.1,0.3,0.5,0.7,0.9,1.0,0.0/

***   stelalr evolution table
*     1. Low mass main sequence (M < 0.7) (0)
*     2. High mass main sequence  (1)
*     3. Hertzsprung gap (HG). (2)
*     4. Red giant. (3)
*     5. Core Helium burning. (HB) (4)
*     6. AGB (5-6)
*     7. Helium types (7-9)
*     8. White dwarf (10-12)
*     9. Neutron star (13)
*     10.Black hole (14)
*     11.Pre main sequence (-1)

*     KW        -1  0  1  2  3  4  5  6  7  8  9 10 11 12 13  14  15
      DATA KTab/11, 1, 2, 3, 4, 5, 6, 6, 7, 7, 7, 8, 8, 8, 9, 10, -1/

      NKW = 0
      MKW = 0

*     escale
      ebscale = RAU_OUT/MSCALE_OUT**2
      ekscale = 1.0/(MSCALE_OUT*VSCALE_OUT**2)
      ascale  = 1.0/(MSCALE_OUT*RSCALE_OUT*VSCALE_OUT)

***   Lagragian radii
      RLAGR  = 0
      RSLAGR = 0
      RBLAGR = 0
      RPLAGR = 0
*     Number
      NLAGR  = 0
      NSLAGR = 0
      NBLAGR = 0
      NPLAGR = 0
*     Mass
      MLAGR  = 0
      MSLAGR = 0 
      MBLAGR = 0
      MPLAGR = 0
*     Velocity
      VLAGR  = 0
      VSLAGR = 0
      VBLAGR = 0
      VPLAGR = 0
*     Velocity dispersion
      SLAGR  = 0
      SSLAGR = 0
      SBLAGR = 0
      SPLAGR = 0
*     Energy
      ELAGR  = 0
      ESLAGR = 0
      EBLAGR = 0
      EPLAGR = 0
*     Average velocity
      VTAVE  = 0
      VTBAVE = 0
      VTSAVE = 0
      VTPAVE = 0
*     binary fraction
      MSBLAGR  = 0
      MSPBLAGR = 0
      NSBLAGR  = 0
      NSPBLAGR = 0
*     energy fraction
      EBINL   = 0
      EBINBL  = 0
      EPBINL  = 0
      EPBINBL = 0
      
C      DATA RLAGR, RSLAGR, RBLAGR, RPLAGR /NN14 * 0.0/ 
C      DATA NLAGR, NSLAGR,  /NN14 * 0/ 
C      DATA MLAGR, MSLAGR, MBLAGR, MPLAGR /NN14 * 0.0/ 
C      DATA VLAGR, VSLAGR, VBLAGR, VPLAGR /NN98 * 0.0/ 
C      DATA SLAGR, SSLAGR, SBLAGR, SPLAGR /NN98 * 0.0/ 
C      DATA ELAGR, ESLAGR, EBLAGR, EPLAGR /NN14 * 0.0/
C      DATA VTAVE, VTBAVE, VTSAVE, VTPAVE  /NN42 * 0.0/
C      DATA MSBLAGR, MSPBLAGR /NN1*0.0, NN1*0.0/
C      DATA NSBLAGR, NSPBLAGR /NN1*0.0, NN1*0.0/
C      DATA EBINL,   EBINBL   /NN*0.0, NN1*0.0/
C      DATA EPBINL,  EPBINBL  /NN*0.0, NN1*0.0/

***   Counters
      C1 = 0
      C2 = 0
      C3 = 0
*     total mass cumulator
      CMASS = 0
      BMASS = 0
      SMASS = 0
      PMASS = 0
*     previous mass
      MCPREV  = 0
      MCSPREV = 0
      MCBPREV = 0
      MCPPREV = 0
*     previous number
      NCPREV  = 0
      NCSPREV = 0
      NCBPREV = 0
      NCPPREV = 0
*     number counter
      NC      = 0
      NCS     = 0
      NCB     = 0
      NCP     = 0
      NPK     = 0
*     bin counter
      KK  = 1
      KKS = 1
      KKB = 1
      KKP = 1
*     total mass
      TMASS  = 0
      TSMASS = 0
      TBMASS = 0
      TPMASS = 0
      
C      DATA C1,C2,C3,CMASS,BMASS,SMASS,PMASS /6*0.0, NTP*0.0/
C      DATA MCPREV, MCSPREV, MCBPREV, MCPPREV /3*0.0, NTP*0.0/
C      DATA NC, NCS, NCB, NCPREV, NCSPREV, NCBPREV /6*0/
C      DATA NCP, NCPPREV, NPK /NTP*0, NTP*0, NTP*0/
C      DATA I,J,K,kk,kks,kkb,kkp /3*0,3*1,NTP*1/
C      DATA TMASS, TSMASS, TBMASS /3*0.0/

***   Global parameters
*     energy
      EPARS = 0
      EKIN  = 0
      EPOT  = 0
      EBIN  = 0
      MM    = 0
      RH    = 0
C      DATA EPARS,EKIN,EPOT,EBIN,MM,RH /17*0.0/

***   Flag
      IF (KZ7.EQ.2.OR.KZ7.EQ.3) THEN
         fshell = .true.
      ELSE
         fshell = .false.
      END IF

***   stellar evolution flag
      IF (KZ19.GE.3) then
         fse=.true.
      else
         fse=.false.
      end if

***   Particle number counts
      NPB   = N_BINARY + N_MERGER
      NP    = N_SINGLE + NPB
      NPS   = N_SINGLE
      NPSB  = N_SINGLE + N_BINARY
      NPR   = N_SINGLE + 2*N_BINARY + 3*N_MERGER

***   Core radius
      RC = REAL(RC8)
      RLAGR(NN) = RC

***** Calculation------------------------
***     Determine the center-of-mass
*     Set square radii of resolved binaries
      C1 = 0
      C2 = 0
      C3 = 0
C!$omp parallel do private(I) 
C!$omp reduction(+:C1) reduction(+:C2) reduction(+:C3)
      DO I = 1, NPS
         C1 = C1 + S_M(I)*S_X1(I)
         C2 = C2 + S_M(I)*S_X2(I)
         C3 = C3 + S_M(I)*S_X3(I)
         TMASS = TMASS + S_M(I)
         TSMASS = TSMASS + S_M(I)
*     stellar evolution count
         IF (fse) then
            kwi(1) = NS_KW(I)
            NKW(1,kwi(1)) = NKW(1,kwi(1)) + 1
            MKW(1,kwi(1)) = MKW(1,kwi(1)) + S_M(I)
            KI = KTab(KWI(1))
            NPK(KI) = NPK(KI) + 1
            TPMASS(KI) = TPMASS(KI) + S_M(I)
         END IF
*     Maximum mass
         MM = MAX(MM,S_M(I))
      END DO
C!$omp end parallel do

C!$omp parallel do private(I,mi) 
C!$omp reduction(+:C1) reduction(+:C2) reduction(+:C3)
      DO I = 1, N_BINARY
         mi = B_M1(I) + B_M2(I)
         C1 = C1 + mi*B_XC1(I)
         C2 = C2 + mi*B_XC2(I)
         C3 = C3 + mi*B_XC3(I)
         TMASS = TMASS + mi
         TBMASS = TBMASS + mi
*     stellar evolution count
         IF (fse) then
            kwi(1) = NB_KW1(I)
            kwi(2) = NB_KW2(I)
            mki(1) = B_M1(I)
            mki(2) = B_M2(I)
            IF (kwi(1).eq.kwi(2)) then
               NKW(3,kwi(1)) = NKW(3,kwi(1)) + 1
               DO L = 1,2
                  MKW(3,kwi(L)) = MKW(3,kwi(L)) + mki(L)
               END DO
            ELSE
               DO L = 1,2
                  NKW(2,kwi(L)) = NKW(2,kwi(L)) + 1
                  MKW(2,kwi(L)) = MKW(2,kwi(L)) + mki(L)
               END DO
            END IF
            DO L = 1,2
               KI = KTab(KWI(L))
               NPK(KI) = NPK(KI) + 1
               TPMASS(KI) = TPMASS(KI) + mki(L)
            END DO
         END IF
*     Maximum mass
         MM = MAX(MM,B_M1(I))
         MM = MAX(MM,B_M2(I))
      END DO
C!$omp end parallel do

      DO I = 1, N_MERGER
         mi = M_M1(I) + M_M2(I) + M_M3(I)
         C1 = C1 + mi*M_XC1(I)
         C2 = C2 + mi*M_XC2(I)
         C3 = C3 + mi*M_XC3(I)
         TMASS = TMASS + mi
         TBMASS = TBMASS + mi
*     stellar evolution count
         IF (fse) then
            kwi(1) = NM_KW1(I)
            kwi(2) = NM_KW2(I)
            kwi(3) = NM_KW3(I)
            mki(1) = M_M1(I)
            mki(2) = M_M2(I)
            mki(3) = M_M3(I)
            IF (kwi(1).eq.kwi(2)) then
               NKW(3,kwi(1)) = NKW(3,kwi(1)) + 1
               DO L = 1,2
                  MKW(3,kwi(L)) = MKW(3,kwi(L)) + mki(L)
               END DO
            ELSE
               DO L = 1,2
                  NKW(2,kwi(L)) = NKW(2,kwi(L)) + 1
                  MKW(2,kwi(L)) = MKW(2,kwi(L)) + mki(L)
               END DO
            END IF
            NKW(1,kwi(3)) = NKW(1,kwi(3)) + 1
            MKW(1,kwi(3)) = MKW(1,kwi(3)) + mki(3)
            DO L = 1,3
               KI = KTab(KWI(L))
               NPK(KI) = NPK(KI) + 1
               TPMASS(KI) = TPMASS(KI) + mki(L)
            END DO
         END IF
*     Maximum mass
         MM = MAX(MM,M_M1(I))
         MM = MAX(MM,M_M2(I))
         MM = MAX(MM,M_M3(I))
      END DO

      C1 = C1/TMASS
      C2 = C2/TMASS
      C3 = C3/TMASS

***   Obtain the index and radius
*     Single stars
C!$omp parallel do private(I)
      DO I = 1, NPS
         R2(I) = (S_X1(I) - C1)**2 + (S_X2(I) - C2)**2 +
     &        (S_X3(I) - C3)**2
         IALIST(I) = I
      END DO
C!$omp end parallel do
*
*     binaries
      IF (N_BINARY.GT.0) THEN
C!$omp parallel do private(I)
         DO I = 1, N_BINARY
            R2(I+NPS) = (B_XC1(I) - C1)**2 + (B_XC2(I) - C2)**2 +
     &           (B_XC3(I) - C3)**2
            IALIST(I+NPS) = I + NPS
         END DO
C!$omp end parallel do
      END IF
      
*     Mergers
      IF (N_MERGER.GT.0) THEN
         DO I = 1, N_MERGER
            J = I + NPSB
            R2(J) = (M_XC1(I) - C1)**2 + (M_XC2(I) - C2)**2 
     &           + (M_XC3(I) - C3)**2
            IALIST(J) = J
         END DO
      END IF
         
***   Sort square distances of all particles with respect to the centre C.
      CALL SORT1F(NP,R2,IALIST)

***   Total mass
      IF (KZ7.EQ.2.OR.KZ7.EQ.4) THEN
*     Get initial total mass
         TMASS_F = REAL(TSMASS0+TBMASS0)
         TSMASS_F = REAL(TSMASS0)
         TBMASS_F = REAL(TBMASS0)
      ELSE
         TMASS_F = TMASS
         TSMASS_F = TSMASS
         TBMASS_F = TBMASS
      END IF
*     half mass
      HFMASS = 0.5*TMASS

***   Mass array initialization
      MLAGR  = TMASS_F  *FLAGR
      MSLAGR = TSMASS_F *FLAGR
      MBLAGR = TBMASS_F *FLAGR
      DO I = 1, NTP
         MPLAGR(:,I) = TPMASS(I) * FLAGR
      END DO

***** Lagragian calculation loop -----------------
***   Loop all stars
      DO 1 I = 1, NP
*     Get index of sorted particles
         J = IALIST(I)
         RI = SQRT(R2(I))
*     Increase total counter
         NC = NC + 1

*     Check the type of stars (single/binary/merger)
         IF (J.GT.NPS) THEN
            NCB = NCB + 1
            NC = NC + 1

            IF (J.GT.NPSB) THEN
               K = J - NPSB
*     Mergers:
               NC = NC + 1
*     cumulative binary bins for binary fraction
               IF (kk.lt.NN)  NSBLAGR(kk) = NSBLAGR(kk) + 1
               mi = M_M1(K) + M_M2(K) + M_M3(K)
               r(1) = M_XC1(K)
               r(2) = M_XC2(K)
               r(3) = M_XC3(K)
               v(1) = M_VC1(K)
               v(2) = M_VC2(K)
               v(3) = M_VC3(K)
               pot  = M_POT(K)
*     Binding energy
               EBI(I) = REAL(0.5*(M_M1(K)*M_M2(K)/M_A0(K) 
     &              + (M_M1(K)+M_M2(K))*M_M3(K)/M_A1(K))*ebscale)
            ELSE
*     Binaries:
*     Get origin index of dataset
               K = J - NPS
*     Get data
               mi = B_M1(K) + B_M2(K)
               n1 = NB_NAM1(K) 
               n2 = NB_NAM2(K) 
               r(1) = B_XC1(K)
               r(2) = B_XC2(K)
               r(3) = B_XC3(K)
               v(1) = B_VC1(K)
               v(2) = B_VC2(K)
               v(3) = B_VC3(K)
               pot  = B_POT(K)
*     Binding energy
               EBI(I) = REAL(0.5*(B_M1(K)*B_M2(K)/B_A(K))*ebscale)
*     primordial binaries:
               if (abs(n1-n2).eq.1) then
                  if (kk.lt.NN) then
                     mspblagr(kk) = mspblagr(kk) + mi
                     nspblagr(kk) = nspblagr(kk) + 2
                     EPBINL(KK)   = EPBINL(KK)   - EBI(I)
                     EPBINBL(KKB) = EPBINL(KKB)  - EBI(I)
                  end if
                  if (ri.le.rc) then
                     EPBINL(NN)  = EPBINL(NN)  - EBI(I)
                     mspblagr(NN) = mspblagr(NN) + mi
                     nspblagr(NN) = nspblagr(NN) + 2
                  end if
               end if
            END IF
*     Get binding energy
            EBIN = EBIN - EBI(I)
*     cumulative binary bins for binary fraction
            IF (kk.lt.NN) THEN
               MSBLAGR(kk) = MSBLAGR(kk) + mi
               NSBLAGR(kk) = NSBLAGR(kk) + 2
               EBINL(KK)   = EBINL(KK)   - EBI(I)
               EBINBL(KKB) = EBINBL(KKB)  - EBI(I)
            END IF
            if (ri.le.rc) EBINL(NN)  = EBINL(NN)   - EBI(I)
         ELSE
*     Singles:
            NCS  = NCS + 1
            mi = S_M(J)
            r(1) = S_X1(J)
            r(2) = S_X2(J)
            r(3) = S_X3(J)
            v(1) = S_V1(J)
            v(2) = S_V2(J)
            v(3) = S_V3(J)
            pot  = S_POT(J)
         END IF

*     Get velocity information
*     x-y plane radial velocity * rxy
         rvxy = r(1)*v(1) + r(2)*v(2)
*     radial velocity value
         vr(I) = (rvxy + r(3)*v(3))/ri;
*     tangential velocity vector
         vt(1:3,I) = v(1:3) - vr(I)*r(1:3)/ri;
*     x-y plane radial position square
         rxy2 = r(1)*r(1) + r(2)*r(2);
*     rotational velocity
         vrot1 = v(1) - rvxy*r(1)/rxy2;
         vrot2 = v(2) - rvxy*r(2)/rxy2;
         vrot(I) = sqrt(vrot1*vrot1+vrot2*vrot2);
*     rotational direction sign
         IF ((vrot1*r(2) - vrot2*r(1)).LT.0.0) vrot(I) = - vrot(I)
*     Get energy
         v2 = (v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
         eki = REAL(0.5*mi*v2*ekscale)
         EI(I) = eki + REAL(0.5*pot*mi/mscale_out)
         EKIN = EKIN + eki
         EPOT = EPOT + REAL(0.5*pot*mi/mscale_out)
*     Angular momentum
         ai(1) = REAL(mi*(r(2)*v(3)-r(3)*v(2))*ascale)
         ai(2) = REAL(mi*(r(3)*v(1)-r(1)*v(3))*ascale)
         ai(3) = REAL(mi*(r(1)*v(2)-r(2)*v(1))*ascale)

*     Fill angular momentum bins
         A(kk,1:3) = A(kk,1:3) + ai(1:3)
         
*     Fill the bins
*     Total:
         kpre = kk
         call lagr_fill(kk,CMASS,NPR,MCPREV,NCPREV,RC,NC,
     &        ri,mi,v,vr(I),vt(1,I),vrot(I),EI(I),
     &        NLAGR,RLAGR,MLAGR,VLAGR,VTAVE,ELAGR,fshell)
*     binary fraction
         IF (kk.lt.NN.and.kpre.lt.kk.and..not.fshell) THEN
            IF (NSBLAGR(KK).EQ.0) THEN
               MSBLAGR(kk) = MSBLAGR(kk-1)
               NSBLAGR(kk) = NSBLAGR(kk-1)
            END IF
            IF (NSPBLAGR(KK).EQ.0) THEN
               MSPBLAGR(kk) = MSPBLAGR(kk-1)
               NSPBLAGR(kk) = NSPBLAGR(kk-1)
            END IF
            IF (A(1,KK).eq.0)  A(:,KK) = A(:,KK-1)
         END IF

*     half mass radius
         IF (CMASS.GE.HFMASS.and.RH.eq.0.0) RH = ri

         IF (J.GT.NPS) THEN
*     Binary:
            call lagr_fill(kkb,BMASS,NPB,MCBPREV,NCBPREV,RC,NCB,
     &           ri,mi,v,vr(I),vt(1,I),vrot(I),EI(I),
     &           NBLAGR,RBLAGR,MBLAGR,VBLAGR,VTBAVE,EBLAGR,fshell)
         ELSE
*     Single:
            call lagr_fill(kks,SMASS,NPS,MCSPREV,NCSPREV,RC,NCS,
     &           ri,mi,v,vr(I),vt(1,I),vrot(I),EI(I),
     &           NSLAGR,RSLAGR,MSLAGR,VSLAGR,VTSAVE,ESLAGR,fshell)
         END IF

*     Stellar types
         IF (fse) then
            IF (J.GT.NPS) THEN
               IF (J.GT.NPSB) THEN
                  KN = 3
                  K = J - NPSB
                  kwi(1) = NM_KW1(K)
                  mki(1) = M_M1(K)
                  kwi(2) = NM_KW2(K)
                  mki(2) = M_M2(K)
                  kwi(3) = NM_KW3(K)
                  mki(3) = M_M3(K)
               ELSE
                  KN = 2
                  K = J - NPS
                  kwi(1) = NB_KW1(K)
                  mki(1) = B_M1(K)
                  kwi(2) = NB_KW2(K)
                  mki(2) = B_M2(K)
               END IF
            ELSE
               KN = 1
               kwi(1) = NS_KW(J)
               mki(1) = S_M(J)
            END IF
*     fill coresponding stellar types
            DO 3 L = 1, KN
               KI = KTab(KWI(L))
               IF (KI.LT.0) GO TO 3
               NCP(KI) = NCP(KI) + 1
               call lagr_fill(kkp(KI),PMASS(KI),NPK(KI),MCPPREV(KI),
     &              NCPPREV(KI),RC,NCP(KI),
     &              ri,mki(L),v,vr(I),vt(1,I),vrot(I),EI(I),
     &              NPLAGR(1,KI),RPLAGR(1,KI),MPLAGR(1,KI),
     &              VPLAGR(1,1,KI),VTPAVE(1,1,KI),EPLAGR(1,KI),fshell)
 3          CONTINUE
         END IF

 1    CONTINUE

*     remember to get the mplagr in first loop

      IF (.not.fshell) THEN
*     Fill emtpy bins
*     total:
         call lagr_empty_fill(kk,NLAGR,RLAGR,MLAGR,VLAGR,VTAVE,ELAGR)
*     single:
         call lagr_empty_fill(kks,NSLAGR,RSLAGR,MSLAGR,VSLAGR,VTSAVE,
     &        ESLAGR)
*     binary:
         call lagr_empty_fill(kkb,NBLAGR,RBLAGR,MBLAGR,VBLAGR,VTBAVE,
     &        EBLAGR)
*     stellar types:
         if (fse) THEN
            DO K=1, NTP
               call lagr_empty_fill(kkp(K),NPLAGR(1,K),RPLAGR(1,K),
     &              MPLAGR(1,K),VPLAGR(1,1,K),VTPAVE(1,1,K),EPLAGR(1,K))
            END DO
         END IF
*     binary fraction:
         IF (KK.LT.NN) THEN
            MSBLAGR(KK:NN-1) = MSBLAGR(KK-1)
            NSBLAGR(KK:NN-1) = NSBLAGR(KK-1)
            MSPBLAGR(KK:NN-1) = MSPBLAGR(KK-1)
            NSPBLAGR(KK:NN-1) = NSPBLAGR(KK-1)
         END IF
      END IF
      
         
*     Average bins
*     total:
      call lagr_average(MLAGR,VLAGR,VTAVE,ELAGR)
*     single:
      call lagr_average(MSLAGR,VSLAGR,VTSAVE,ESLAGR)
*     binary:
      call lagr_average(MBLAGR,VBLAGR,VTBAVE,EBLAGR)
*     stellar types:
      IF (fse) THEN
         DO K=1, NTP
            call lagr_average(MPLAGR(1,K),VPLAGR(1,1,K),
     &           VTPAVE(1,1,K),EPLAGR(1,K))
         END DO
      END IF

*     Velocity dispersion
*     initialization
      NC  = 0
      NCS = 0
      NCB = 0
      NCP = 0

      kk  = 1
      kks = 1
      kkb = 1
      kkp = 1

      MCPREV  = 0.0
      MCSPREV = 0.0
      MCBPREV = 0.0
      MCPPREV = 0.0

      NCPREV  = 0
      NCSPREV = 0
      NCBPREV = 0
      NCPPREV = 0
      
*     Loop all stars
      DO 2 I = 1, NP
*     Get index of sorted particles
         J = IALIST(I)
         RI = SQRT(R2(I))
*     Increase total counter
         NC = NC + 1

*     Check the type of stars (single/binary/merger)
         IF (J.GT.NPS) THEN
            NCB = NCB + 1
            NC = NC + 1

            IF (J.GT.NPSB) THEN
*     Mergers:
               K = J - NPSB
               NC = NC + 1
               mi = M_M1(K) + M_M2(K) + M_M3(K)
               v(1)   = M_VC1(K)
               v(2)   = M_VC2(K)
               v(3)   = M_VC3(K)
            ELSE
*     Binaries:
*     Get origin index of dataset
               K = J - NPS
*     Get data
               mi = B_M1(K) + B_M2(K)
               v(1)  = B_VC1(K)
               v(2)  = B_VC2(K)
               v(3)  = B_VC3(K)
            END IF

*     Calculate velocity dispersion
            call lagr_sig_fill(kkb,NCB,NCBPREV,ri,rc,
     &           mi,v,vr(I),vt(1,I),vrot(I),
     &           NBLAGR,VBLAGR,VTBAVE,SBLAGR,fshell)

         ELSE
*     Singles:
            NCS  = NCS + 1
            mi = S_M(J)
            v(1)  = S_V1(J)
            v(2)  = S_V2(J)
            v(3)  = S_V3(J)

*     Calculate velocity dispersion
            call lagr_sig_fill(kks,NCS,NCSPREV,ri,rc,
     &           mi,v,vr(I),vt(1,I),vrot(I),
     &           NSLAGR,VSLAGR,VTSAVE,SSLAGR,fshell)
         END IF
         
*     Calculate velocity dispersion
         call lagr_sig_fill(kk,NC,NCPREV,ri,rc,
     &        mi,v,vr(I),vt(1,I),vrot(I),
     &        NLAGR,VLAGR,VTAVE,SLAGR,fshell)

*     Stellar types
         IF (fse) then
            IF (J.GT.NPS) THEN
               IF (J.GT.NPSB) THEN
                  KN = 3
                  K = J - NPSB
                  kwi(1) = NM_KW1(K)
                  mki(1) = M_M1(K)
                  kwi(2) = NM_KW2(K)
                  mki(2) = M_M2(K)
                  kwi(3) = NM_KW3(K)
                  mki(3) = M_M3(K)
               ELSE
                  KN = 2
                  K = J - NPS
                  kwi(1) = NB_KW1(K)
                  mki(1) = B_M1(K)
                  kwi(2) = NB_KW2(K)
                  mki(2) = B_M2(K)
               END IF
            ELSE
               KN = 1
               kwi(1) = NS_KW(J)
               mki(1) = S_M(J)
            END IF
*     fill coresponding stellar types
            DO 4 L = 1, KN
               KI = KTab(KWI(L))
               IF (KI.LT.0) GO TO 4
               NCP(KI) = NCP(KI) + 1
*     Calculate velocity dispersion
               call lagr_sig_fill(kkp(KI),NCP(KI),NCPPREV(KI),ri,rc,
     &              mki(L),v,vr(I),vt(1,I),vrot(I),
     &              NPLAGR(1,KI),VPLAGR(1,1,KI),VTPAVE(1,1,KI),
     &              SPLAGR(1,1,KI),fshell)
 4          CONTINUE
         END IF
         
 2    CONTINUE

*     Fill empty bins for non-shell mode
      if (.not.fshell) then
*     total:
         call lagr_sig_empty_fill(kk,SLAGR)
*     single:
         call lagr_sig_empty_fill(kks,SSLAGR)
*     binary:
         call lagr_sig_empty_fill(kkb,SBLAGR)
*     stellar types:
         IF (fse) THEN
            DO I = 1,NTP
               call lagr_sig_empty_fill(kkp(I),SPLAGR(1,1,I))
            END DO
         END IF
      end if

*     Get mass average
*     total:
      call lagr_sig_average(NLAGR,MLAGR,SLAGR)
*     single:
      call lagr_sig_average(NSLAGR,MSLAGR,SSLAGR)
*     binary:
      call lagr_sig_average(NBLAGR,MBLAGR,SBLAGR)
*     stellar types:
      IF (fse) THEN
         DO I = 1,NTP
            call lagr_sig_average(NPLAGR(1,I),MPLAGR(1,I),SPLAGR(1,1,I))
         END DO
      END IF

***** Get global parameters for output ----------------
***   Energy
      EPARS(2) = EKIN
      EPARS(3) = EPOT
      EPARS(4) = EBIN
***   Mass
*     5.   TMass: total mass (Msun)
      FPARS(5) = TMASS
*     6.   TSMass: single mass (Msun)
      FPARS(6) = TSMASS
*     7.   TBMass: binary mass (Msun)
      FPARS(7) = TBMASS
***   Radius
*     9.   Rh: half-mass radius
      FPARS(9) = Rh
*     11-13: Rden: center-of-mass center position
      FPARS(11) = C1
      FPARS(12) = C2
      FPARS(13) = C3
*     Maximum stellar mass
      FPARS(16) = MM

***   Number
*     1. Total number (resolved)
      NPARS(1) = NPR
*     2. Total single number
      NPARS(2) = NCS
*     3. Total binary number + merger number (unresolved)
      NPARS(3) = NPB
*     4. Total Merger number (unresolved)
      NPARS(4) = N_MERGER
*     5. Total number (unresolved)
      NPARS(5) = NP


      if (fse) then
***   Stellar evolution
         DO I = 0,16
*     Mass
            SFPARS(8+3*I:10+3*I) = MKW(:,I-1)
*     Number
            SNPARS(22+3*I:24+3*I) = NKW(:,I-1)
         END DO
      end if

      call global_params_gether(EPARS,FPARS,NPARS,SFPARS,SNPARS)

***** print-----------------
*     FPARS, EPARS, NPARS, LAGR(all,Single.bianry), Binary fraction (all, primordial), binary binding energy (all, primordial), Angular momentum (unresolved binareis), SFPARS, SNPARS, LAGR(stellar types)
*     Each Lagr group (9) R,  N,  M, V, Vx, Vy, Vz, Vr, Vt, Vrot, S, Sx, Sy, Sz, Sr, St, Srot, e
*     Fraction: 0.001,0.01,0.1,0.3,0.5,0.7,0.9,1.0, RC
            
      write(36,*) 
     &     FPARS, EPARS, NPARS,
     &     RLAGR , NLAGR , MLAGR , VLAGR , SLAGR , ELAGR,
     &     RSLAGR, NSLAGR, MSLAGR, VSLAGR, SSLAGR, ESLAGR,
     &     RBLAGR, NBLAGR, MBLAGR, VBLAGR, SBLAGR, EBLAGR,
     &     MSBLAGR, NSBLAGR, MSPBLAGR, NSPBLAGR,
     &     EBINL,   EBINBL, EPBINL, EPBINBL, A,
     &     SFPARS, SNPARS,
     &     (RPLAGR(:,K), NPLAGR(:,K), MPLAGR(:,K), VPLAGR(:,:,K), 
     &     SPLAGR(:,:,K), EPLAGR(:,K), K=1,NTP)

      RETURN
      
      END
      
**************************************
      subroutine lagr_fill(kb,mass,ntot,mpre,npre,rc,
     &     n,ri,mi,v,vr,vt,vrot,e,
     &     nl,rl,ml,vl,vta,el,fshell)
*
*      
*     Lagr fill function
*     -----------------------
*
      implicit none
      INTEGER NN
      parameter (NN=9)
      INTEGER kb,n,ntot,npre,nl(NN)
      REAL mass,mi,ri,v(3),vr,vt(3),vrot,e,mpre,rc
      REAL rl(NN),ml(NN),vl(NN,7),vta(3,NN),el(NN)
      Logical fshell
      
      mass = mass + mi

*     average velocity
      vl(kb,2:4) = vl(kb,2:4) + mi*v(1:3) 
      vl(kb,5)   = vl(kb,5)   + mi*vr 
      vl(kb,7)   = vl(kb,7)   + mi*vrot
      vta(1:3,kb)= vta(1:3,kb)+ mi*vt(1:3)
      el(kb)     = el(kb)     + e

*     Core:
      IF (ri.LE.rc) THEN
         ml(NN) = ml(NN) + mi
         nl(NN) = nl(NN) + 1
         vl(NN,2:4) = vl(NN,2:4) + mi*v(1:3) 
         vl(NN,5)   = vl(NN,5)   + mi*vr 
         vl(NN,7)   = vl(NN,7)   + mi*vrot
         vta(1:3,NN)= vta(1:3,NN)+ mi*vt(1:3)
         el(NN)     = el(NN)     + e
      END IF

      
*     Go to next bin if mass reach the R_lagr limit
      if (kb.lt.NN) THEN
         if ((mass.ge.ml(kb)).or.(kb.eq.NN-1.and.n.eq.ntot)) then
*     update mass
            ml(kb) = mass
*     Get R_lagr for single
            rl(kb) = ri
*     Get number for single
            nl(kb) = n
*     For shell cases:
            if (fshell) then
               ml(kb) = ml(kb) - mpre
               nl(kb) = nl(kb) - npre
            end if
*     Update previous mass and number
            mpre = mass
            npre = n
*     increase bin index
            kb = kb + 1
*     initial next bins
            if (kb.lt.NN) then
               if(.not.fshell.and.nl(kb).eq.0) then
                  vl(kb,2:5) = vl(kb-1,2:5)
                  vl(kb,7)   = vl(kb-1,7)
                  vta(1:3,kb)= vta(1:3,kb-1)
                  el(kb)     = el(kb-1)
               end if
            end if
         end if
      end if
      
      return

      end

*************************************
      subroutine lagr_empty_fill(kb,nl,rl,ml,vl,vta,el)
*
*
*     Fill empty bins of lagragian parameters
*     ----------------------------
*
      INTEGER NN
      parameter(NN=9)
      INTEGER kb,nl(NN),I,kn
      REAL rl(NN),ml(NN),vl(NN,7),vta(3,NN),el(NN)

      IF (KB.Lt.NN) then
         kn = kb -1
         rl(kb:NN-1) = rl(kn) 
         nl(kb:NN-1) = nl(kn) 
         ml(kb:NN-1) = ml(kn)
         DO I = 2,5
            vl(kb:NN-1,I)= vl(kn,I)
         END DO
         vl(kb:NN-1,7) = vl(kn,7)
         el(kb:NN-1) = el(kn)
         DO I = kb, NN-1
            vta(1:3,I) = vta(1:3,kn)
         END DO
      end if
      
      return

      end

*************************************
      subroutine lagr_average(ml,vl,vta,el)
*
*
*     Average bins for lagrangian parameters
*     ----------------------------
*
      INTEGER NN
      parameter(NN=9)
      INTEGER I
      REAL ml(NN),vl(NN,7),vta(3,NN),el(NN)
      
      DO I = 1, NN
         if (ml(I).gt.0.0) then
            vl(I,2:5)  = vl(I,2:5) / ml(I)
            vta(1:3,I) = vta(1:3,I) / ml(I)
            vl(I,7)    = vl(I,7) / ml(I)
            el(I)      = el(I) / ml(I)
            vl(I,1) = SQRT(vl(I,2)**2+vl(I,3)**2+vl(I,4)**2)
            vl(I,6) = SQRT(vta(1,I)**2+vta(2,I)**2+vta(3,I)**2)
         end if
      END DO
      
      return

      end

*************************************
      subroutine lagr_sig_fill(kb,n,npre,ri,rc,mi,v,vr,vt,vrot,
     &     nl,vl,vta,sl,fshell)
*
*
*     Average bins for lagrangian parameters
*     ----------------------------
*
      parameter(NN=9)
      INTEGER nl(NN),n,npre
      REAL mi,v(3),vr,vt(3),vrot
      REAL vl(NN,7),vta(3,NN),sl(NN,7)
      logical fshell
      
      REAL dr(6),dt(3),dr2(6)
      if(kb.gt.NN) print*, n,npre
*     x,y,z
      dr(1:3) = v(1:3) - vl(kb,2:4)
      dr(4)   = vr     - vl(kb,5)
      dt(1:3) = vt(1:3)- vta(1:3,kb)
      dr(6)   = vrot   - vl(kb,7)

*     mass weighted square
      dr2(1:4) = mi*dr(1:4)**2
      dr2(5)   = mi*(dt(1)*dt(1) + dt(2)*dt(2) + dt(3)*dt(3))
      dr2(6)   = mi*dr(6)**2

*     dispersion
      sl(kb,2:7) = sl(kb,2:7) + dr2(1:6)

      IF (ri.le.rc) sl(NN,2:7) = sl(NN,2:7) + dr2(1:6)

*     check whether need to increase kk
      if (n-npre.eq.nl(kb)) then
         if (fshell) npre = npre + nl(kb)
         kb = kb + 1
         if (.not.fshell.and.kb.lt.NN) sl(kb,2:7) = sl(kb-1,2:7)
      end if

      return

      end

*************************************
      subroutine lagr_sig_empty_fill(kb,sl)
*
*
*     Fill empty bins of sig. lagragian parameters
*     ----------------------------
*
      INTEGER NN
      parameter(NN=9)
      INTEGER kb,I,kn
      REAL sl(NN,7)

      IF (KB.Lt.NN) THEN
         kn = kb - 1
         DO I = 2,7
            sl(kb:NN-1,I) = sl(kn,I)
         END DO
      END IF

      return

      end

*************************************
      subroutine lagr_sig_average(nl,ml,sl)
*
*
*     Average bins for sig. lagrangian parameters
*     ----------------------------
*
      parameter(NN=9)
      INTEGER nl(NN),I
      REAL ml(NN),sl(NN,7)
      
      DO I = 1, NN
         if (nl(I).gt.0) then
            sl(I,2:7) = sl(I,2:7) / ml(I)
            sl(I,1)   = sl(I,2) + sl(I,3) + sl(I,4)
            ml(I)  = ml(I) / nl(I)
         end if
      END DO
      
      return

      end
