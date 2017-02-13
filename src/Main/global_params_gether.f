      subroutine global_params_gether(epars,fpars,npars,sfpars,snpars)
*
*
*     Gether global parameters for output
*     -----------------------------------
*     
      include 'common6.h'
      include 'scale_out.h'
*     TMASS: total mass (NB)
      REAL epars(12), FPARS(16), SFPARS(58)
      INTEGER NPARS(5),SNPARS(72)

      REAL TMASS, TSMASS, TBMASS, EKM

*     Notice [A] means defautled is Astronomical Unit, otherwise it is NB unit if KZ(12)=-1 & KZ(19)=3

***** Energy ----------------------------
*     2. EKIN: kinetic energy
      ZKIN = EPARS(2)
*     3. EPOT: Potential energy
      POT = -EPARS(3)
*     4. EBIN: Binary binding energy
*     5. ETIDE: Tidal energy
*       Obtain the tidal potential energy for linearized external field. 
      IF (KZ(14).EQ.0) THEN
*       Note: EPL holds accumulated tidal energy if KZ(14) = 3.
         VIR = 0.0
         EPL = 0.0D0
      ELSE
*       Employ general expression sum {m*r*F} for virial energy.
         CALL XTRNLV(1,N)
*       Form tidal energy with Plummer potential (note EPL use for #14=3).
         IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
            EPL = 0.0
            DO 50 I = 1,N
               RI2 = AP2
               DO 45 K = 1,3
                  RI2 = RI2 + X(K,I)**2
 45            CONTINUE
               EPL = EPL - BODY(I)*MP/SQRT(RI2)
 50         CONTINUE
         END IF
      END IF
      EPARS(5) = REAL(EPL)
      IF(KZ(14).NE.3.and.KZ(14).NE.9) EPARS(5) = EPARS(5) + REAL(ETIDE)

*     6. EMDOT: mass loss energy
      EPARS(6)  = REAL(EMDOT)
*     7. ECOLL: collision energy
      EPARS(7)  = REAL(ECOLL)
*     8. ECMDOT: common envelope energy
      EPARS(8)  = REAL(ECDOT)
*     1. ETOT: total energy without escaping energy
      EPARS(1)  = SUM(EPARS(2:8))
*     9. EKICK: kick energy
      EPARS(9)  = REAL(EKICK)
*     10. ESESC: escaper kinetic and potential energy
      EPARS(10) = REAL(ESESC)
*     11. EBESC: binary escaper binding energy
      EPARS(11) = REAL(EBESC)
*     12. EMESC: triple escaper binding energy
      EPARS(12) = REAL(EMESC)

*     Mechanical energy
      EKM = SUM(EPARS(2:4))
      
***** Global parameters      
*     1.   Time[N]: time of simulation
      FPARS(1) = REAL(TTOT)

*     2.   Time[A]: time of simulation (Myr)
      FPARS(2) = REAL(TTOT*TSCALE_OUT)

*     5.   TMass: total mass (Msun)
      TMASS = REAL(FPARS(5)/MSCALE_OUT)

*     6.   TSMass: single mass (Msun)
      TSMASS = REAL(FPARS(6)/MSCALE_OUT)

*     7.   TBMass: binary mass (Msun)
      TBMASS = REAL(FPARS(7)/MSCALE_OUT)

*     8.  Virial ratio
*     Include non-zero virial energy for Plummer potential and/or 3D case.
      VIR =  POT - VIR
*     Note angular momentum term is contained in virial energy (#14=1/2).
      Q =  ZKIN/VIR
      FPARS(8) = REAL(Q)
      
*     3. Tcr[A]:  half-mass radius crossing time (Myr)
*     Define crossing time using single particle energy (cf. option 14).
      TCR = TMASS**2.5/(2.0*ABS(EKM))**1.5
*     Define crossing time and save single particle energy.
      IF (Q.GT.1.0.AND.KZ(14).LT.3) THEN
         TCR = TCR*SQRT(2.0*Q)
      END IF
      FPARS(3) = REAL(TCR*TSCALE_OUT)

*     9.   Rh[A]:   half mass radius
      Rh = FPARS(9)/RSCALE_OUT
      
*     4.   Trh[A]:  half-mass radius relaxation time 
*     Unresolved total number
      ALLN = FLOAT(NPARS(5))
*     Form equilibrium rms velocity (temporarily defined as VC).
      VC = SQRT(2.0D0*ABS(EKM)/TMASS)
      TRH = 4.0*TWOPI/3.0*(VC*Rh)**3/(15.4*TMASS**2*LOG(ALLN)/ALLN)
      FPARS(4) = REAL(TRH*TSCALE_OUT)

*     10.   Rtid[A]: Tidal radius
      FPARS(10) = REAL(RTIDE*RSCALE_OUT)

*     11-13. R_den[A]: Density center position
      
*     14.  Rho_d[A]: Density weighted average density ΣRHO2/ΣRHO 
*                    (RHO:Mass density of individual star calculated by nearest 5 neighbors.
*                     only avaiable for particles inside core radius)
      FPARS(14) = REAL(RHOD*DSCALE_OUT)

*     15.  Rho_m[A]: Maximum mass density / half mass mean value
      FPARS(15) = REAL(RHOM*DSCALE_OUT)

***** Number ---------------------------
C*     6. NESC: Escaper number
C      NPARS(6) = NESC

      if (KZ(19).ge.3) then
***** Stellar evolution-----------------
*     Number:
*     1. NDISS: Tidal dissipation at pericenter (#27 >0)
         SNPARS(1) = NDISS
*     2. NTIDE: Tidal captures (#27 >0)
         SNPARS(2) = NTIDE
*     3. NSYNC: Synchronous binaries (#27 >0)
         SNPARS(3) = NSYNC
*     4. NCOLL: stellar collision
         SNPARS(4) = NCOLL
*     5. NCOAL: stellar coalescence
         SNPARS(5) = NCOAL
*     6. NCIRC: circularized binaries (#27 >0)
         SNPARS(6) = NCIRC
*     7. NROCHE: Roche stage triggered times
         SNPARS(7) = NROCHE
*     8. NRO: Roche binary events
         SNPARS(8) = NRO
*     9. NCE: Common envelope binaries
         SNPARS(9) = NCE
*     10. NHYP: Hyperbolic collision
         SNPARS(10) = NHYP
*     11. NHYPC: Hyperbolic common envelope binaries
         SNPARS(11) = NHYPC
*     12. NKICK: WD/NS/BH kick
         SNPARS(12) = NKICK
*     13. NMDOT: Stellar mass loss event
         SNPARS(13) = NMDOT
*     14. NRG: New red giants
         SNPARS(14) = NRG
*     15. NHE: New helium stars
         SNPARS(15) = NHE
*     16. NRS: New red supergiants
         SNPARS(16) = NRS
*     17. NNH: New naked helium stars
         SNPARS(17) = NNH
*     18. NWD: New white dwarfs
         SNPARS(18) = NWD
*     19. NSN: New neutron stars
         SNPARS(19) = NSN
*     20. NBH: New black holes
         SNPARS(20) = NBH
*     21. NBS: New blue stragglers
         SNPARS(21) = NBS
*     22-72: NKW: single, binary (one), binary (two) groups (3 number together) from KW = -1 to 15
         
*     mass
*     1. NMDOT: Stellar mass loss 
         SFPARS(1) = REAL(ZMDOT)
*     2. NRG: New red giant mass
         SFPARS(2) = REAL(ZMRG)
*     3. NHE: New helium star mass
         SFPARS(3) = REAL(ZMHE)
*     4. NRS: New red supergiant mass
         SFPARS(4) = REAL(ZMRS)
*     5. NNH: New naked helium star mass
         SFPARS(5) = REAL(ZMNH)
*     6. NWD: New white dwarf mass
         SFPARS(6) = REAL(ZMWD)
*     7. NSN: New neutron star mass
         SFPARS(7) = REAL(ZMSN)
*     8-58: MKW: single, binary (one), binary (two) groups (3 mass together) fro KW = -1 to 15
      end if

      return

      end

