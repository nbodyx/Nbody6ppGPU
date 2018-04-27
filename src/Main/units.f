      SUBROUTINE UNITS
*
*
*       Initialization of units & scaling factors.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Define GM & PC in cgs units and #AU in pc (2008 IAU values).
      GM = 1.32712442099D+26
      AU = 1.49597870700D+13
      PC = 1296000.0D0/TWOPI*AU !3.0856776e+18
*
*       Form scaling factors for binary periods A*SQRT(A/M) to yrs and days.
      YRS = (RBAR*1296000.0D0/TWOPI)**1.5/SQRT(ZMBAR)
      DAYS = 365.25*YRS
*
*       Specify conversion factors for lengths to solar radii & AU.
      RSUN = 6.957D+10
      SU = PC/RSUN*RBAR
      RAU = PC/AU*RBAR
*
*       Copy solar mass scaling to new variable (M = BODY*<M>).
      SMU = ZMBAR
*
*     Form time scale in seconds and velocity scale in km/sec.
      TSTAR = SQRT(PC/GM)*PC

*     Convert time scale from units of seconds to million years.
      TSTAR = TSTAR/(3.15576D+07*1.0D+06)

      IF (KZ(22).NE.10) THEN
         VSTAR = 1.0D-05*SQRT(GM/PC)
*
*       Ensure ZMBAR & RBAR > 0 (=0: assume <M>/Sun = 1, RBAR = 1 pc).
         IF (ZMBAR.LE.0.0D0) ZMBAR = FLOAT(N)/ZMASS
         IF (RBAR.LE.0.0D0) RBAR = 1.0
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
         TSTAR = TSTAR*SQRT(RBAR**3/(ZMASS*ZMBAR))
         VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
      ELSE
*     USE scale factor from original input data instead
         TSTAR = TSTAR*SQRT(RBAR**3/ZMBAR)
      END IF
*
*       Copy TSTAR to secondary time-scale factor.
      TSCALE = TSTAR
*
*       Physical scaling: X, M, V, T from RBAR*X, ZMBAR*M, VSTAR*V, TSTAR*T.
      if(rank.eq.0)
     &WRITE (6,10)  RBAR, ZMBAR, VSTAR, TSTAR, BODYM*ZMBAR, SU
 10   FORMAT (/,12X,1P,'PHYSICAL SCALING:    R* =',E21.12,'  M* =',
     &     E21.12,'  V* =',E21.12,'  T* =',E21.12,'  <M> =',E21.12,
     &     '  SU =',E21.12)
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
