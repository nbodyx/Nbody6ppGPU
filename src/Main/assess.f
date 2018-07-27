      SUBROUTINE ASSESS(IPAIR,IM,ECC,SEMI,ITERM)
*
*
*       Assessment of hierarchical stability.
*       -------------------------------------
*
      INCLUDE 'common6.h'
*     Safe for parallel
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8  XX(3,3),VV(3,3)
*     for diagnostics, safe for parallel
      SAVE ITIME
      DATA ITIME /0/
*
*
*       Define indices for components and c.m.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      ICM = N + IPAIR
*
*       Determine inclination between inner relative motion and outer orbit.
      RV = 0.0
      DO 10 K = 1,3
          XX(K,1) = XREL(K,IM)
          XX(K,2) = 0.0
          XX(K,3) = X(K,I2)
          VV(K,1) = VREL(K,IM)
          VV(K,2) = 0.0
          VV(K,3) = XDOT(K,I2)
          RV = RV + XREL(K,IM)*VREL(K,IM)
   10 CONTINUE
*
      CALL INCLIN(XX,VV,X(1,ICM),XDOT(1,ICM),ANGLE)
*
*       Form inner eccentricity (neglect radial velocity for minimum).
      RIN = SQRT(XREL(1,IM)**2 + XREL(2,IM)**2 + XREL(3,IM)**2)
      SEMI0 = -0.5*BODY(I1)/HM(IM)
      ECC2 = (1.0 - RIN/SEMI0)**2 + RV**2/(BODY(I1)*SEMI0)
      ECC0 = SQRT(ECC2)
      PMIN = SEMI*(1.0 - ECC)
*
*       Evaluate the general stability function.
      IF (ECC.LT.1.0) THEN
*       Include de/dt corrections inside generous stability boundary.
          EOUT = ECC
          IF (PMIN.LT.6.0*SEMI0.OR.EOUT.GT.0.9) THEN
              NNB1 = LIST(1,I1) + 1
              DE = 0.0
              DO 20 L = 2,NNB1      ! Sum positive perturber contributions.
                  JP = LIST(L,I1)
                  CALL EDOT(I1,I2,JP,SEMI,ECC,ECCDOT)
                  IF (ECCDOT.GT.0.0) THEN
                      TK = TWOPI*SEMI*(SQRT(SEMI)/BODY(ICM))
                      DE = DE - MIN(ECCDOT*TK,0.001D0)
                  END IF
   20         CONTINUE
              EOUT = EOUT - DE
              EOUT = MAX(EOUT,0.0D0)
          END IF
*       Obtain new stability boundary from Valtonen's criterion (MN 2017).
          QST = QSTAB(ECC0,EOUT,ANGLE,CM(1,IM),CM(2,IM),BODY(I2))
*
          ITIME = ITIME + 1
          IF (ITIME.GT.2000000000) ITIME = 0
          IF (MOD(ITIME,1000).EQ.0) THEN
              ALPH = 360.0*ANGLE/TWOPI
              PCRIT0 = QST*SEMI0
              WRITE (6,30) ECC0, ECC, ALPH, SEMI, PCRIT0, PMIN
   30         FORMAT (' ASSESS    E0 E1 INC A1 PCR PMIN ',
     &                            2F8.4,F7.1,1P,3E9.1)
          END IF
*       Define termination by QSTAB.
          IF (QST*SEMI0.LT.PMIN) THEN
              ITERM = 0
          ELSE
              ITERM = 1
          END IF
      END IF
      
C          EOUT = ECC
C*       Modify outer eccentricity for consistency with acceptance.
C          IF (EOUT.GT.0.80) THEN
C              DE = 0.5*(1.0 - EOUT)
C              DE = MIN(DE,0.01D0)
C              EOUT = ECC - DE
C              PMIN = SEMI*(1.0 - EOUT)
C          END IF
C          NST = NSTAB(SEMI0,SEMI,ECC0,EOUT,ANGLE,CM(1,IM),
C     &                                    CM(2,IM),BODY(I2))
C          IF (NST.EQ.0) THEN
C              PCRIT = 0.98*PMIN
C              PCR = stability(CM(1,IM),CM(2,IM),BODY(I2),
C     &                                          ECC0,ECC,ANGLE)*SEMI0
C*        Reduce termination distance if old criterion < PMIN/2.
C              IF (PCR.LT.0.5*PMIN) THEN
C                  PCRIT = 0.75*PMIN
C              END IF
C              ITIME = ITIME + 1
C              IF (ITIME.GT.2000000000) ITIME = 0
C              IF (MOD(ITIME,1000).EQ.0) THEN
C                  ALPH = 360.0*ANGLE/TWOPI
C                  WRITE (6,20)  ECC0, ECC, ALPH, SEMI, PCRIT, PCR,
C     &                          R0(IPAIR)
C   20             FORMAT (' ASSESS    E0 E1 INC A1 PCR PC0 R0 ',
C     &                                2F7.3,F7.1,1P,4E9.1)
C              END IF
C          ELSE
C              PCRIT = 1.01*PMIN
C          END IF
C      ELSE
C          PCRIT = stability(CM(1,IM),CM(2,IM),BODY(I2),
C     &                                        ECC0,ECC,ANGLE)*SEMI0
C      END IF
C*
C*       Set new stability distance or define termination.
C      IF (PCRIT.LT.PMIN) THEN
C          ITERM = 0
C          R0(IPAIR) = PCRIT
C      ELSE
C          ITERM = 1
C      END IF
*
      RETURN
*
      END
