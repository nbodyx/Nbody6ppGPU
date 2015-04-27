      SUBROUTINE UNITS
*
*
*       Initialization of units & scaling factors.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
#ifdef POSTNEWTON
      include 'postnewton.h'
#endif
*
*
*       Define GM & PC in cgs units and #AU in pc (2008 IAU values).
      GM = 6.6743D-08*1.9884D+33
      PC = 3.0856776D+18
      AU = PC/1.4959787D+13  ! bug fix from D+12 4/7/11 (was OK in 2009).
*
*       Form scaling factors for binary periods A*SQRT(A/M) to yrs and days.
      YRS = (RBAR*AU)**1.5/SQRT(ZMBAR)
      DAYS = 365.24*YRS
*
*       Specify conversion factors for lengths to solar radii & AU.
      SU = PC/(AU*6.955D+10)*RBAR*AU
      RAU = RBAR*AU
*
*       Copy solar mass scaling to new variable (M = BODY*<M>).
      SMU = ZMBAR
*
*       Form time scale in seconds and velocity scale in km/sec.
      TSTAR = SQRT(PC/GM)*PC
      VSTAR = 1.0D-05*SQRT(GM/PC)
*
*       Convert time scale from units of seconds to million years.
      TSTAR = TSTAR/(3.15D+07*1.0D+06)
*
*       Ensure ZMBAR & RBAR > 0 (=0: assume <M>/Sun = 1, RBAR = 1 pc).
      IF (ZMBAR.LE.0.0D0) ZMBAR = FLOAT(N)/ZMASS
      IF (RBAR.LE.0.0D0) RBAR = 1.0
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      TSTAR = TSTAR*SQRT(RBAR**3/(ZMASS*ZMBAR))
      VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
*
*       Copy TSTAR to secondary time-scale factor.
      TSCALE = TSTAR
*
#ifdef POSTNEWTON
*       Scale clight into NB unit
      CLIGHTPN = CLIN/VSTAR
#endif
*       Physical scaling: X, M, V, T from RBAR*X, ZMBAR*M, VSTAR*V, TSTAR*T.
      if(rank.eq.0)
     &WRITE (6,10)  RBAR, ZMBAR, VSTAR, TSTAR, BODYM*ZMBAR, SU
 10   FORMAT (/,12X,1P,'PHYSICAL SCALING:    R* =',E21.12,'  M* =',
     &     E21.12,'  V* =',E21.12,'  T* =',E21.12,'  <M> =',E21.12,
     &     '  SU =',E21.12)
#ifdef POSTNEWTON
      if(rank.eq.0) THEN
*     Real light speed in km/s
         creal=299792.458
*     Input light speed / real light speed
         cratio=clin/creal
         WRITE (6,15) CLIGHTPN,RBAR*cratio**2,TSCALE*cratio**3,
     &        VSTAR/cratio,ZMBAR
 15      FORMAT (/,12X,1P,'PHYSICAL SCALING FOR POSTNEWTON:',
     &        '  Clight[NB] =',E21.12, '  R* =',E21.12,
     &        '  T* =',E21.12,'  V* =',E21.12,'  M* =',E21.12)
      end if
#endif
*
*       Define relevant parameter for the GR case (RZ = 6*<m>/c^2).
      IF (KZ(27).EQ.3.OR.KZ(28).GT.0) THEN
          CLIGHT = 3.0D+05/VSTAR
          RZ = 6.0*ZMASS/(FLOAT(N)*CLIGHT**2)
          if(rank.eq.0)
     &    WRITE (6,20)  VSTAR, CLIGHT, RZ
   20     FORMAT (/,12X,'GR SCALING:    V* =',1P,E10.2,'  C =',E10.2,
     &                                  '  RZ =',E10.2)
      END IF
*
      RETURN
*
      END
