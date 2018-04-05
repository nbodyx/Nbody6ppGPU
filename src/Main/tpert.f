      SUBROUTINE TPERT(IPAIR,GA,DT)
*
*
*       Perturbation time scale.
*       ------------------------
*
      INCLUDE 'common6.h'
*
*
*       Set c.m. index and initialize scalars.
      I = N + IPAIR
      FMAX = 0.0
      JCL = -1
      DTIN = 1.0E+20
      NNB1 = LIST(1,I) + 1
*
*       Check rare case of no neighbours.
C      IF (LIST(1,I).EQ.0) THEN
      IF (NNB1.LE.1) THEN
          DT = 2.0*STEPR(I)
          GO TO 20
      END IF
*
      IF (KZ(50).EQ.0) then
*       Find the most likely perturbers (first approach & maximum force).
         DO 10 L = 2,NNB1
            J = LIST(L,I)
            call jpred(J,TIME,TIME)
            RIJ2 = 0.0
            RDOT = 0.0
*     
            DO 6 K = 1,3
               XREL = X(K,J) - X(K,I)
               VREL = XDOT(K,J) - XDOT(K,I)
               RIJ2 = RIJ2 + XREL**2
               RDOT = RDOT + XREL*VREL
 6          CONTINUE
*     
            VR = RDOT/RIJ2
            IF (VR.LT.DTIN) THEN
               DTIN = VR
*     Note DTIN is inverse travel time to include case of no RDOT < 0.
               RCRIT2 = RIJ2
               JCRIT = J
            END IF
            FIJ = (BODY(I) + BODY(J))/RIJ2
            IF (FIJ.GT.FMAX) THEN
               FMAX = FIJ
               RJMIN2 = RIJ2
               JCL = J
            END IF
 10      CONTINUE
*
      ELSE
         JCRIT = IMINR(I)
         JCL = JCRIT
      END IF
      
      IF(JCL.LE.0) THEN
         DT = STEP(I)
         GO TO 20
      END IF
      CALL JPRED(JCRIT,TIME,TIME)
      DO K = 1,3
         XREL = X(K,JCRIT) - X(K,I)
         VREL = XDOT(K,JCRIT) - XDOT(K,I)
         RIJ2 = RIJ2 + XREL**2
         RDOT = RDOT + XREL*VREL
      END DO
      DTIN = RDOT/RIJ2
      RCRIT2 = RIJ2
*
*       Form radial velocity of body with shortest approach time (if any).
      RCRIT = SQRT(RCRIT2)
      RDOT = RCRIT*ABS(DTIN)
      A1 = 2.0/(BODY(I)*GA)
      SEMI = -0.5*BODY(I)/H(IPAIR)
*
*       Use the actual apocentre for unperturbed travel time.
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      RI = SEMI*(1.0 + SQRT(ECC2))
*
*       Estimate time interval to reach tidal perturbation of GA.
      DT = (RCRIT - RI*(BODY(JCRIT)*A1)**0.3333)/RDOT
*
*       Compare the travel time based on acceleration only.
      DTMAX = SQRT(2.0D0*ABS(DT)*RDOT*RCRIT2/(BODY(I) + BODY(JCRIT)))
      DT = MIN(DT,DTMAX)
*
*       Skip dominant force test if there is only one critical body.
C      IF (JCRIT.NE.JCL) THEN
C*       Form the return time of the dominant body and choose the minimum.
C          DR = SQRT(RJMIN2) - RI*(BODY(JCL)*A1)**0.3333
C          DTMAX = SQRT(2.0D0*ABS(DR)/FMAX)
C          DT = MIN(DT,DTMAX)
C      END IF
*
*       Apply safety test in case background force dominates c.m. motion.
      DT = MIN(DT,4.0D0*STEP(I))
      DT = MAX(DT,0.0D0)
*
   20 RETURN
*
      END
