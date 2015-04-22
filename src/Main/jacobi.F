      SUBROUTINE JACOBI(NESC)
*
*
*       Jacobi escape criterion.
*       ------------------------
*
      INCLUDE 'common6.h'
*
#ifdef TT
*** FlorentR - The evaluation of the number of stars and the mass is
***    done iteratively, according to an energy criterion. 
***    Idea from Mark Gieles.

      INTEGER NUNBOUND(NMAX)
      NUNBOUND = 0
      ZMBOUND = ZMASS
      NOLDESC2 = 0
      NITER = 0

 80   CONTINUE
      NITER = NITER + 1
      NESC = 0 
      NESC2 = 0

      ZMESC = 0
      ZMESC2 = 0

*** FRenaud
#endif
*
*       Specify escape energy (tidal field or isolated system).
      IF (KZ(14).GT.0.AND.KZ(14).LT.2) THEN
#ifdef TT
*** FlorentR - use ZMBOUND instead of ZMASS
          ECRIT = -1.5*(TIDAL(1)*ZMBOUND**2)**0.333
#else
          ECRIT = -1.5*(TIDAL(1)*ZMASS**2)**0.333
#endif
      ELSE
          ECRIT = 0.0
      END IF
*       Define current mean square velocity (= 0.5 initially).
#ifdef TT
*** FlorentR - use ZMBOUND instead of ZMASS
      V2M = 0.5*ZMBOUND/RSCALE
#else
      V2M = 0.5*ZMASS/RSCALE
#endif
*
*       Count all escapers.
#ifndef TT
      NESC = 0
#endif
      DO 60 I = IFIRST,NTOT
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
          VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
          IF (RI2.GT.RSCALE**2) THEN
              VC2 = ZMASS/SQRT(RI2) + ECRIT
              IF (VI2.LT.VC2) GO TO 60
          ELSE
*       Skip velocities below sqrt(2) times equilibrium value.
              IF (VI2.LT.2.0*V2M + ECRIT) GO TO 60
          END IF
          POTI = 0.0
          DO 55 J = IFIRST,NTOT
              IF (J.EQ.I) GO TO 55
              RIJ2 = 0.0
              DO 52 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   52         CONTINUE
#ifdef TT
*** FlorentR - skip contribution of stars that have been identified as unbound in first pass
            IF (NUNBOUND(J).EQ.0) POTI = POTI + BODY(J)/SQRT(RIJ2)
*** FRenaud
#else
              POTI = POTI + BODY(J)/SQRT(RIJ2)
#endif
   55     CONTINUE
          EI = 0.5*VI2 - POTI
#ifdef TT
*** FlorentR - compute the escaping mass
*          IF (EI.GT.ECRIT) NESC = NESC + 1
         IF (EI.GT.ECRIT) THEN
            NESC = NESC + 1
            ZMESC = ZMESC + BODY(I)
            IF (I.GT.N) THEN
               NESC = NESC + 1
            ENDIF
         ENDIF
 
         IF (EI.GT.0) THEN
            NESC2 = NESC2 + 1
            ZMESC2 = ZMESC2 + BODY(I)
            NUNBOUND(I) = 1
            IF (I.GT.N) THEN
               NESC2 = NESC2 + 1
            ENDIF
         ENDIF
*** FRenaud
#else
          IF (EI.GT.ECRIT) NESC = NESC + 1
#endif
   60 CONTINUE
*
#ifdef TT
*** FlorentR - iterative scheme and output
      ZMBOUND = ZMASS - ZMESC
      
      IF(NITER.LT.5) THEN
        IF(ABS(NESC2-NOLDESC2).GT.NESC2/10.0) THEN
* re-compute but with updated ZMBOUND
          NOLDESC2 = NESC2
          GO TO 80
        ENDIF
      ELSE
        WRITE(6,*) "Warning: Max nb of iterations reached in Jacobi.f"
      ENDIF

      WRITE(6,90) TTOT,ZMASS,ZMASS-ZMESC,ZMASS-ZMESC2,N,N-NESC,N-NESC2
 90   FORMAT(" EBOUND   TTOT ZMASS MBOUND MBOUND2   N NBOUND NBOUND2 ",
     &     f10.3,3(f10.4),3(i7))

C 70   CONTINUE
*** FRenaud
#endif
      RETURN
*
      END
