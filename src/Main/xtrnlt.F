      SUBROUTINE XTRNLT(XI,XIDOT,FREG,FDR)
*
*
*       Galactic force & first derivative.
*       ----------------------------------
*
      INCLUDE 'common6.h'
#ifdef TT
      INCLUDE 'tt.h'
#endif
      INCLUDE 'galaxy.h'
      REAL*8  XI(3),XIDOT(3),FREG(3),FDR(3),FM(3),FMD(3)
*
*
*       Consider point-mass, disk and/or logarithmic halo model.
      IF (KZ(14).EQ.3) THEN
*       Employ global instead of linearized forms for better accuracy.
          IF (GMG.GT.0.0D0) THEN
              CALL FNUC(XI,XIDOT,FM,FMD)
              DO 10 K = 1,3
                  FREG(K) = FREG(K) + FM(K)
                  FDR(K) = FDR(K) + FMD(K)
   10         CONTINUE
          END IF
*
*       Check bulge force.
      IF (GMB.GT.0.0D0) THEN
          CALL FBULGE(XI,XIDOT,FM,FMD)
          DO 15 K = 1,3
              FREG(K) = FREG(K) + FM(K)
              FDR(K) = FDR(K) + FMD(K)
   15     CONTINUE
      END IF
*
*       Include Miyamoto disk for positive disk mass.
          IF (DISK.GT.0.0D0) THEN
              CALL FDISK(XI,XIDOT,FM,FMD)
              DO 20 K = 1,3
                  FREG(K) = FREG(K) + FM(K)
                  FDR(K) = FDR(K) + FMD(K)
   20         CONTINUE
          END IF
*
*       Check addition of logarithmic halo potential to regular force.
          IF (V02.GT.0.0D0) THEN
              CALL FHALO(XI,XIDOT,FM,FMD)
              DO 30 K = 1,3
                  FREG(K) = FREG(K) + FM(K)
                  FDR(K) = FDR(K) + FMD(K)
   30         CONTINUE
          END IF
      END IF
#ifdef TT
*** FlorentR - Compute the force from the user-definition of the pot.
      IF( (KZ(14).EQ.9).AND.TTMODE.EQ.0) THEN
        CALL TTFORCE(XI,XIDOT,FM,FMD,0)
          DO 40 K = 1,3
            FREG(K) = FREG(K) + FM(K)
            FDR(K) = FDR(K) + FMD(K)
 40      CONTINUE
      END IF
*** FRenaud
#endif
      RETURN
*
      END
