      SUBROUTINE GCINT
*
*
*       Integration of 3D cluster orbit.
*       --------------------------------
*
      INCLUDE 'common6.h'
#ifdef TT
      INCLUDE 'tt.h'
#endif
      INCLUDE 'galaxy.h'
      REAL*8  FM(3),FD(3),FS(3),FSD(3),FD2(3),FD3(3)
*
*     Check whether it is necessary to integrate
      DT = TIME + TOFF - TGNEXT
      IF (DT.LT.0.0) RETURN

*       Predict coordinates and velocities to order FDOT.
      DT = TIME + TOFF - TG
*       Note: integration step may exceed STEPX (depends on block-step).
      DT2 = 0.5*DT
      DT3 = ONE3*DT
      RG2 = 0.0
      RGVG = 0.0
      DO 10 K = 1,3
          RG(K) = ((FGD(K)*DT3 + FG(K))*DT2 + VG(K))*DT + RG(K)
          VG(K) = (FGD(K)*DT2 + FG(K))*DT + VG(K)
          RG2 = RG2 + RG(K)**2
          RGVG = RGVG + RG(K)*VG(K)
   10 CONTINUE
*

      FM = 0.0
      FD = 0.0
*       Obtain force and first derivative of point-mass galaxy.
      IF (KZ(14).EQ.3) THEN
         IF (GMG.GT.0.0D0) THEN
            CALL FNUC(RG,VG,FM,FD)
         END IF
*
*       Check bulge force.
         IF (GMB.GT.0.0D0) THEN
            CALL FBULGE(RG,VG,FS,FSD)
            DO 20 K = 1,3
               FM(K) = FM(K) + FS(K)
               FD(K) = FD(K) + FSD(K)
 20         CONTINUE
         END IF
*
*     Include optional Miyamoto disk component (Book eq. 8.52).
         IF (DISK.GT.0.0D0) THEN
            CALL FDISK(RG,VG,FS,FSD)
            DO 25 K = 1,3
               FM(K) = FM(K) + FS(K)
               FD(K) = FD(K) + FSD(K)
 25         CONTINUE
         END IF
*
*     Check addition of logarithmic galaxy potential (not linearized).
         IF (V02.GT.0.0D0) THEN
            CALL FHALO(RG,VG,FS,FSD)
            DO 30 K = 1,3
               FM(K) = FM(K) + FS(K)
               FD(K) = FD(K) + FSD(K)
 30         CONTINUE
         END IF
      END IF
*
*     MWpotential14
      IF (KZ(14).EQ.5) THEN
         call fmwpot(RG, VG, FM ,FD)
      END IF

#ifdef TT
*** FlorentR - Compute the force from the user-definition of the pot.
*      equivalent to the FNUC call below  
      IF( (KZ(14).EQ.9).AND.TTMODE.EQ.0) THEN
         CALL TTFORCE(RG,VG,FM,FD,DT)
      END IF
*** FRenaud
#endif
*       Set time factors for corrector.
      DT13 = ONE3*DT
      DTSQ = DT**2
      DTSQ12 = ONE12*DTSQ
      DT2 = 2.0/DTSQ
      DT6 = 6.0/(DT*DTSQ)
      TG = TG + DT
*
*       Include the Hermite corrector and update F & FDOT.
      DO 40 K = 1,3
          DF = FG(K) - FM(K)
          SUM = FGD(K) + FD(K)
          AT3 = 2.0*DF + SUM*DT
          BT2 = -3.0*DF - (SUM + FGD(K))*DT
          RG(K) = RG(K) + (0.6*AT3 + BT2)*DTSQ12
          VG(K) = VG(K) + (0.75*AT3 + BT2)*DT13
          FG(K) = FM(K)
          FGD(K) = FD(K)
          FD2(K) = (3.0*AT3+BT2)*DT2
          FD3(K) = AT3*DT6
   40 CONTINUE

      TTMP = TSTEP(FG,FGD,FD2,FD3,ETAI)
      TGNEXT = TG + MIN(TTMP,STEPX)

C      write(6,*) 'GC int: RG ' ,RG,' VG ',VG, 
C     &     ' DT ',TTMP,' TGNEXT ',TGNEXT, ' TG ',TG
      
*
*       Update angular velocity in case of non-circular orbit.
*     OM1 = (RG(2)*VG(3) - RG(3)*VG(2))/RG2
*     OM2 = (RG(3)*VG(1) - RG(1)*VG(3))/RG2
      OM3 = (RG(1)*VG(2) - RG(2)*VG(1))/RG2
*     OMEGA2 = OM1**2 + OM2**2 + OM3**2
      OMEGA2 = OM3**2
      OMEGA = SQRT(OMEGA2)
      TIDAL(4) = 2.0*OMEGA
*
      RETURN
*
      END
