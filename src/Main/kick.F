      SUBROUTINE KICK(I,ICASE,KW,DM)
*
*
*       Velocity kick for WD, neutron stars or black holes.
*       ---------------------------------------------------
*
      Include 'kspars.h'
      INCLUDE 'common6.h'
      REAL*8  RAN2,VK(4)
*     Unsafe for parallel
      SAVE  IPAIR, KC, VDIS, RI
      DATA  IPAIR,KC /0,0/
*       Reinstall the traditional kicks (March 2012)
      PARAMETER (VFAC=-1.0D0)
*       Get fallback from hrdiag.f
*     Safe for parallel
      common /fall/fallback
*
*     --03/07/14 11:30-lwang-ks-parallel--------------------------------*
***** Note: Here when icase eq. -147 means for ks communication
*****       I will be the processer id to share data
      IF (ICASE.EQ.-147) THEN
#ifdef PARALLEL         
         CALL MPI_BCAST(IPAIR,1,MPI_INTEGER,I,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(KC,1,MPI_INTEGER,I,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(VDIS,1,MPI_REAL8,I,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(RI,1,MPI_REAL8,I,MPI_COMM_WORLD,ierr)
#endif         
         RETURN
      END IF
*     --03/07/14 11:30-lwang-end-ks-parallel----------------------------*
      
*       Save orbital parameters in case of KS binary (called from KSAPO).
      IF (ICASE.EQ.0) THEN
          IPAIR = I
          KC = KSTAR(N+IPAIR)
*       Identify the correct component (KSTAR reversed in MDOT).
          I1 = 2*IPAIR - 1
          IF (KSTAR(I1).LT.0) THEN
              IN = I1
          ELSE IF (KSTAR(I1+1).LT.0) THEN
              IN = I1 + 1
          END IF
          KSTAR(IN) = -KSTAR(IN)
*     ks MPI communication KSTAR
          call ksparmpi(K_store,K_int,K_KSTAR,IN,0,KSTAR(IN))
*
*     When call expel in cmbody and WD/NS binary form after Binary CE, 
*     There is big energy error when DM is large. It seems here the DM 
*     will be set to zero and then cause the issue. Thus suppress now.
C      Determine mass loss and actual disruption velocity.
C          DM = BODY(IN) - 1.4/ZMBAR
C          IF (KW.LT.13) DM = 0.0
C          VD2 = 2.0*(BODY(N+IPAIR) - DM)/R(IPAIR)
          VD2 = 2.0*(BODY(N+IPAIR))/R(IPAIR)
          VDIS = SQRT(VD2)*VSTAR
*       Set cluster escape velocity (add twice central potential).
C          VP2 = 2.0*BODY(N+IPAIR)/R(IPAIR)
          VESC = SQRT(VD2 + 4.0)*VSTAR
          SEMI = -0.5*BODY(N+IPAIR)/H(IPAIR)
          ZM1 = BODY(I1)*SMU
          ZM2 = BODY(I1+1)*SMU
          EB = BODY(I1)*BODY(I1+1)/BODY(N+IPAIR)*H(IPAIR)
          RI = R(IPAIR)
*       Skip on #25 = 0/1 for consistent net WD modification of EKICK.
          IF ((KW.LT.13.AND.KZ(25).EQ.0).OR.
     &        (KW.EQ.12.AND.KZ(25).NE.2)) GO TO 30
*       Sum whole binding energy (used by BINOUT for net change).
          EKICK = EKICK + EB
          EGRAV = EGRAV + EB
*     ks MPI communicaton EKICK EGRAV
          call ksparmpi(K_store,K_real8,K_EGRAV,0,0,EB)
          call ksparmpi(K_store,K_real8,K_EKICK,0,0,EB)
          I2 = I1 + 1
          if(rank.eq.0)
     &    WRITE (6,1)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2), ZM1,
     &                 ZM2, VESC, VDIS, R(IPAIR)/SEMI, EB, R(IPAIR)
    1     FORMAT (' BINARY KICK:    NAME(I1) NAME(I2) K*(I1) K*(I2) ',
     &         'M(I1)[M*] M(I2)[M*] VESC[km/s] VDIS[km/s] R12/SEMI ',
     &         'EB R12[NB] ',2I10,2I4,4F7.1,F6.2,F9.4,1P,E10.2)
          NBKICK = NBKICK + 1
*       Remove any circularized KS binary from the CHAOS/SYNCH table.
          IF (KSTAR(N+IPAIR).GE.10.AND.NCHAOS.GT.0) THEN
              II = -(N + IPAIR)
              CALL SPIRAL(II)
              KSTAR(N+IPAIR) = 0
*     ks MPI communication KSTAR
              call ksparmpi(K_store,K_int,K_KSTAR,N+IPAIR,0,
     &             KSTAR(N+IPAIR))
          END IF
          GO TO 30
      END IF
*
*       Generate velocity kick for neutron star (Gordon Drukier Tokyo paper).
*     IT = 0
*     V0 = 330.0
*     VCUT = 1000.0
*   2 VT = VCUT/V0*RAN2(IDUM1)
*     VP = VT*(2.0*RAN2(IDUM1) - 1.0)
*     VN = SQRT(VT**2 - VP**2)
*     FAC = 1.0/0.847*VN**0.3/(1.0 + VN**3.3)
*     IF (FAC.LT.RAN2(IDUM1).AND.IT.LT.10) GO TO 2
*     VKICK = V0*VT
*
C*       Adopt the Maxwellian of Hansen & Phinney (MNRAS 291, 569, 1997).
C      DISP = 190.0
*
*       Adopt the Maxwellian of Hobbs (MNRAS 360, 974, 2005)
      DISP = 265.0
*
*       Include velocity dispersion in terms of VSTAR (Parameter statement!).
      IF (VFAC.GT.0.0D0) THEN
          DISP = VFAC*VSTAR
      END IF
*
*       Allow for optional type-dependent WD kick.
      IF (KW.LT.13) THEN
          IF (KZ(25).EQ.1.AND.(KW.EQ.10.OR.KW.EQ.11)) THEN
              DISP = 5.0
          ELSE IF (KZ(25).GT.1.AND.KW.EQ.12) THEN
              DISP = 5.0
          ELSE
              DISP = 0.0
          END IF
      END IF
*
*       Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
      DO 2 K = 1,2
          X1 = RAN2(IDUM1)
          X2 = RAN2(IDUM1)
*       Generate two velocities from polar coordinates S & THETA.
          S = DISP*SQRT(-2.0*LOG(1.0 - X1))
          THETA = TWOPI*X2
          VK(2*K-1) = S*COS(THETA)
          VK(2*K) = S*SIN(THETA)
    2 CONTINUE
      VKICK = SQRT(VK(1)**2 + VK(2)**2 + VK(3)**2)
      VK(4) = VKICK
      IF(KW.EQ.14) VKICK = VKICK*(1.D0-fallback)
*
*       Limit kick velocity to VDIS+10*VSTAR/10*VST for binary/single stars.
C      IF (IPAIR.GT.0) THEN
C          VBF = SQRT(VDIS**2 + 100.0*VSTAR**2)
C          VKICK = MIN(VKICK,VBF)
C      ELSE
C          VKICK = MIN(VKICK,10.0D0*VSTAR)
*       Ensure escape of massless star.
      IF (IPAIR.LE.0.AND.BODY(I).EQ.0.0D0) VKICK = 10.0*VSTAR
C      END IF
      VKICK = VKICK/VSTAR
*
      IF (VKICK.EQ.0.0) GO TO 30
*
*       Randomize the velocity components.
*     A(4) = 0.0
*     DO 5 K = 1,3
*         A(K) = 2.0*RAN2(IDUM1) - 1.0
*         A(4) = A(4) + A(K)**2
*   5 CONTINUE
*
*       Add truncated/full kick velocity and initialize X0DOT.
      VI2 = 0.0
      VF2 = 0.0
      CALL JPRED(I,TIME,TIME)
      DO 10 K = 1,3
          VI2 = VI2 + XDOT(K,I)**2
*         XDOT(K,I) = XDOT(K,I) + VKICK*A(K)/SQRT(A(4))
          XDOT(K,I) = XDOT(K,I) + VKICK*VK(K)/VK(4)
          X0DOT(K,I) = XDOT(K,I)
*     ks MPI communication
          call ksparmpi(K_store,K_real8,K_X0DOT,k,i,x0dot(k,i))
          VF2 = VF2 + XDOT(K,I)**2
   10 CONTINUE
*
*       Modify energy loss due to increased velocity of single particle.
      DETMP = - 0.5*BODY(I)*(VF2 - VI2)
      ECDOT = ECDOT + DETMP
*     ks MPI communication
      call ksparmpi(K_store,K_real8,K_ECDOT,0,0,DETMP)
      NKICK = NKICK + 1
*
*       Replace final velocity by relative velocity for binary kick.
      IF (IPAIR.GT.0) THEN
          JP = KVEC(I)
          J = I + 1
          IF (I.EQ.2*JP) J = I - 1
          VF2 = 0.0
          CALL JPRED(J,TIME,TIME)
          DO 15 K = 1,3
              VF2 = VF2 + (XDOT(K,I) - XDOT(K,J))**2
   15     CONTINUE
          HNEW = 0.5*VF2 - (BODY(I) + BODY(J))/RI
          EB1 = BODY(I)*BODY(J)/(BODY(I) + BODY(J))*HNEW
          IF (EB1.LT.0.0) THEN
              EKICK = EKICK - EB1
              EGRAV = EGRAV - EB1
*     ks MPI communicaton EKICK EGRAV
              call ksparmpi(K_store,K_real8,K_EGRAV,0,0,-EB1)
              call ksparmpi(K_store,K_real8,K_EKICK,0,0,-EB1)
          END IF
          IPAIR = 0
      END IF
*
*     IF (NKICK.LT.50.OR.NAME(I).LE.2*NBIN0.OR.
*    &    (KW.GE.13.AND.TTOT*TSTAR.GT.100.0)) THEN
          ZM = BODY(I)*ZMBAR
          if(rank.eq.0)
     &    WRITE (6,20)  I, NAME(I), (TIME+TOFF)*TSTAR, KSTAR(I), KW, KC,
     &         BODY0(I)*ZMBAR,ZM,SQRT(VI2)*VSTAR, VKICK*VSTAR, 
     &         SQRT(VF2)*VSTAR, VK(4),fallback
   20     FORMAT (' VELOCITY KICK: I',I10'  NAME',I10,
     &         '  Time[Myr]',E12.5,'  K*0',I4,
     &         '  K*',I4,'  K*(ICM)',I4'  M0[M*]',F9.4,'  MN[M*]',F9.4,
     &         '  VI[km/s]',F9.4,'  VK[km/s]',F9.4'  VF[km/s]',F9.4,
     &         '  VK0[km/s]',F9.4,'  FB',F9.4)
          KC = 0
*     END IF
*
*       Highlight BH/NS velocities below 4 times rms velocity.
      IF (VKICK.LT.4.0*SQRT(0.5).AND.KW.GE.13) THEN
          if(rank.eq.0)
     &    WRITE (6,25)  I, NAME(I), KW, VKICK*VSTAR, SQRT(VF2)*VSTAR
   25     FORMAT (' LOW KICK:    I',I10,'  NAME',I10,'  K*',I4,
     &         '  VK[km/s]',F7.2,'  VF[km/s]',F7.2)
      END IF
*
*       Include optional list of high-velocity particles (double counting!).
*     IF (KZ(37).GT.0) THEN
*         CALL HIVEL(I)
*     END IF
*
   30 RETURN
*
      END
