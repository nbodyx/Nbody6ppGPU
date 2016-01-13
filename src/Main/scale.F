      SUBROUTINE SCALE
*
*
*       Scaling to new units.
*       ---------------------
*
      INCLUDE 'common6.h'
      LOGICAL LSCALE
*      Scaling if initial model is constructed or King model is read
      LSCALE = KZ(22).LT.2.OR.KZ(22).GE.6
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#ifdef ENSEMBLE
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
#endif
*
*       Read virial ratio, rotation scaling factors, boundary radius & tidal radius.
      if(rank.eq.0)then
      READ (5,*)  Q, VXROT, VZROT, RTIDE

      end if
*      SMAX = 0.5
#if MPIINIT
*     print*,' MPIINIT selected to be 1'
      CALL MPI_BCAST(Q,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(VXROT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(VZROT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RTIDE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
*
*       Note RTIDE should be non-zero for isolated systems (cf. CALL LAGR).
      RSPH2 = RTIDE
      QVIR = Q
      
      ZMASS = 0.0D0
      DO 10 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   10 CONTINUE
*
*       Form total mass and centre of mass displacements.
      DO 30 I = 1,N
          ZMASS = ZMASS + BODY(I)
          DO 25 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   25     CONTINUE
   30 CONTINUE
*
*       Adjust coordinates and velocities to c.m. rest frame.
*       (if start model data were read from other source do NOT scale)(R.Sp.)
      IF (LSCALE) THEN
         DO 40 I = 1,N
            DO 35 K = 1,3
               X(K,I) = X(K,I) - CMR(K)/ZMASS
               XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
 35         CONTINUE
 40      CONTINUE
*     
*     Scale masses to standard units of <M> = 1/N.
*     Skip scaling of masses for unscaled upload or planetesimal disk.
          IF(KZ(22).GE.2.OR.KZ(5).NE.3) THEN
            DO 50 I = 1,N
               BODY(I) = BODY(I)/ZMASS
 50         CONTINUE
         END IF
*
*     Determine astronomical scale from original data
*     Assume Mass use the unit of M_sun 
         IF(KZ(22).EQ.10) ZMBAR = ZMASS/FLOAT(N)
         
         ZMASS = 1.D0
*
      END IF
*
      GRAV = 1.0D0
*     Astronomical unit in non-scaled input data
      IF (LSCALE) GRAV = 4.30145D-3*ZMBAR*FLOAT(N)
      
*     --10/22/13 17:12-lwang-binary-scale-------------------------------*
***** Note:------------------------------------------------------------**
      IF (KZ(8).EQ.2) THEN
*       Sum the total energy of pairs.
         NPAIRS = NBIN0
         IFIRST = 2*NPAIRS + 1
         NTOT = N + NPAIRS
*
         if(rank.eq.0)
     *        PRINT*,' Start Energy Scale: N,NTOT,NPAIRS,IFIRST,GRAV ',
     *        N,NTOT,NPAIRS,IFIRST,GRAV
*
         XMBTOT = 0.D0
         IF (NBIN0.GT.0) THEN
            EBIN = 0.0D0
            DO IPAIR = 1,NBIN0
               ICOMP = 2*IPAIR - 1
               JCOMP = 2*IPAIR
               RIJ2 = 0.0D0
               RDOT = 0.0D0
               VIJ2 = 0.0D0
               J = N + IPAIR
               XMB = BODY(ICOMP) + BODY(JCOMP)
               BODY(J) = XMB
               XMBTOT = XMBTOT + XMB
               DO K = 1,3
                  RIJ2 = RIJ2 + (X(K,ICOMP) - X(K,JCOMP))**2
                  RDOT = RDOT + (X(K,ICOMP) - X(K,JCOMP))*
     &                 (XDOT(K,ICOMP) - XDOT(K,JCOMP))
                  VIJ2 = VIJ2 + (XDOT(K,ICOMP) - XDOT(K,JCOMP))**2
                  X(K,J) = (BODY(ICOMP)*X(K,ICOMP) + 
     &                 BODY(JCOMP)*X(K,JCOMP))/XMB 
                  XDOT(K,J) = (BODY(ICOMP)*XDOT(K,ICOMP) + 
     &                 BODY(JCOMP)*XDOT(K,JCOMP))/XMB
               END DO
               RIJ = DSQRT(RIJ2)
*     
               SEMI = 1.D0/(2.0D0/RIJ - VIJ2/(GRAV*XMB))
               ECC2 = (1.0D0 - RIJ/SEMI)**2 + RDOT**2/(SEMI*GRAV*XMB)
               ECC = DSQRT(ECC2)
               EBIN = EBIN + 0.5D0*GRAV*BODY(ICOMP)*BODY(JCOMP)/SEMI
*      if(rank.eq.0)WRITE(67,111)IPAIR,SEMI,ECC,XMB,(X(K,J),K=1,3),
*    &    (XDOT(K,J),K=1,3)
*       IF(RIJ.LT.RMIN.and.rank.eq.0)
*    &     WRITE(68,112)IPAIR,RIJ,RDOT,DSQRT(VIJ2),
*    &     GRAV*XMB/SEMI,(1.0D0 - RIJ/SEMI)**2,
*    &     RDOT**2/(SEMI*GRAV*XMB)
            END DO
*111  FORMAT(1X,I5,1P,12D15.5)
*112  FORMAT(1X,I5,1P,12D15.5)
         END IF
      END IF

*     --10/22/13 17:12-lwang-end----------------------------------------*
*       Obtain the total kinetic & potential energy.
#ifdef PARALLEL
      CALL ENERGY_MPI(.true.)
#else
      CALL ENERGY(.true.)
#endif
*
*     --10/22/13 18:40-lwang-binary-scale-------------------------------*
***** Note: Reset value to zero ---------------------------------------**
      if(KZ(8).EQ.2) THEN
         if(rank.eq.0)then
            PRINT*,' Energies Absolute Scale ',ZKIN,POT,EBIN
            PRINT*,' Energies per mass ',ZKIN/ZMASS,POT/ZMASS,EBIN/ZMASS
            PRINT*,' Total Mass in Binaries ',XMBTOT
            call flush(6)
         end if
         EBIN = 0
         NPAIRS = 0
         IFIRST = 1
         NTOT = N
      end if
*     --10/22/13 18:41-lwang-end----------------------------------------*
*       Use generalized virial theorem for external tidal field.
      IF (KZ(14).GT.0) THEN
          AZ = 0.0D0
          DO 55 I = 1,N
              AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
   55     CONTINUE
          IF (KZ(14).EQ.1) THEN
*       Use Chandrasekhar eq. (5.535) for virial ratio (rotating frame only).
              VIR = POT - 2.0*(ETIDE + 0.5*TIDAL(4)*AZ)
          ELSE
              VIR = POT - 2.0*ETIDE
          END IF
      ELSE
          VIR = POT
      END IF
*
*     Initial scaling factor
      SX = 1.0
      
*       Allow two optional ways of skipping standard velocity scaling.
      IF (KZ(22).LE.1.AND.(KZ(5).EQ.2.OR.KZ(5).EQ.3)) THEN
         QV = SQRT(Q*VIR/ZKIN)
         E0 = ZKIN*QV**2 - POT + ETIDE
*     Rescale velocities to new masses for two Plummer spheres.
         IF (KZ(5).EQ.2) THEN
            ZKIN = 0.0
            DO 57 I = 1,N
               DO 56 K = 1,3
                  XDOT(K,I) = XDOT(K,I)*QV
                  ZKIN = ZKIN + 0.5*BODY(I)*XDOT(K,I)**2
 56            CONTINUE
 57         CONTINUE
            E0 = ZKIN - POT + ETIDE
            ETOT = E0
            Q = ZKIN/POT
            if(rank.eq.0) WRITE (6,59)  E0, ZKIN/POT
 59         FORMAT (/,12X,'UNSCALED ENERGY    E =',F10.6,
     &           '  Q =',F6.2)
         ELSE
            IF (KZ(5).EQ.3) E0 = ZKIN - POT
            if(rank.eq.0) WRITE (6,54)  E0, ZKIN, POT
 54         FORMAT (/,12X,'UNSCALED ENERGY    E =',F10.6,
     &           ' ZKIN =',F10.6,' POT =',F10.6 )
         END IF

*       Scale non-zero velocities by virial theorem ratio.
*       (if start model data were read from other source do NOT scale)(R.Sp.)
      ELSE IF (LSCALE) THEN
*     Determine Q from original data and initialize VSTAR
         IF (KZ(22).EQ.10) THEN
            Q = ZKIN/(GRAV*VIR)
            QVIR = Q
            VSTAR = 1.0D0
         END IF
         IF (ZKIN.GT.0.0D0) THEN
            QV = SQRT(Q*VIR/ZKIN)
            DO 60 I = 1,N
               DO 58 K = 1,3
                  XDOT(K,I) = XDOT(K,I)*QV
 58            CONTINUE
 60         CONTINUE
            IF (KZ(22).EQ.10) VSTAR = 1.0D0/QV
         END IF
*
*     Scale total energy to standard units (E = -0.25 for Q < 1).
         E0 = -0.25
         ETOT = (Q - 1.0)*POT
*     Include case of hot system inside reflecting boundary.
         IF (KZ(29).GT.0.AND.Q.GT.1.0) THEN
            E0 = ETOT
         END IF

*     Define scaling factor (set E0 = ETOT if energy scaling not desired).
         SX = E0/ETOT

*     Safe check
         IF (SX.LE.0) then
            print*, 'Error!: Scaling factor <= 0. Virial ratio Q >= 1.',
     &           ' Scaled energy E0 = ',E0,
     &           ' Original energy ETOT = ',ETOT
            call abort()
         END IF
*
*     Scale coordinates & velocities to the new units.
         DO 70 I = 1,N
            DO 68 K = 1,3
               X(K,I) = X(K,I)/SX
               XDOT(K,I) = XDOT(K,I)*SQRT(SX)
 68         CONTINUE
 70      CONTINUE
*     Determine astronomical scale from original data
*     Assume original data use M_sun, PC and km/s as units
         IF (KZ(22).EQ.10) THEN
            VSTAR = VSTAR/SQRT(SX)
            RBAR = SX
         END IF
         
      ELSE
         ETOT = ZKIN - POT
         E0 = ETOT
      END IF

*     Print scaling information
      if(rank.eq.0)
     *     WRITE (6,65)  SX, ETOT, BODY(1), BODY(N), ZMASS/FLOAT(N),
     &     Q
 65   FORMAT (//,12X,'SCALING:   SX =',1P,D13.5,'  E =',E10.2,
     &     '  M(1) =',E9.2,'  M(N) =',E9.2,'  <M> =',E9.2,
     &     '  Q =',0P,F10.5)
*     
      call flush(6)

*     In case of no tidal field choose very large RTIDE (R.Sp.)
C      IF (KZ(14).EQ.0.AND.KZ(23).LE.2) RTIDE = 1.D8
*       In case of King model scale initial tidal radius
C      IF(KZ(23).GE.3)THEN
C         RTIDE = RTIDE/SX
c$$$      if(rank.eq.0)PRINT*,' RTIDE =',RTIDE,' ETID=',ZMASS/RTIDE
C      END IF
*
*       Check whether to include rotation (VXROT = 0 in standard case).
      IF (VXROT.GT.0.0D0) THEN
*
*       Set angular velocity for retrograde motion (i.e. star clusters).
          OMEGA = -SX*SQRT(ZMASS*SX)
       if(rank.eq.0)
     *    WRITE (6,75)  VXROT, VZROT, OMEGA
   75     FORMAT (/,12X,'VXROT =',F6.2,'  VZROT =',F6.2,
     &                                                 '  OMEGA =',F7.2)
*
*       Add solid-body rotation about Z-axis (reduce random velocities).
          DO 80 I = 1,N
              XDOT(1,I) = XDOT(1,I)*VXROT - X(2,I)*OMEGA
              XDOT(2,I) = XDOT(2,I)*VXROT + X(1,I)*OMEGA
              XDOT(3,I) = XDOT(3,I)*VZROT
   80     CONTINUE
      END IF
*
*       Set initial crossing time in scaled units.
      TCR = ZMASS**2.5/(2.0D0*ABS(E0))**1.5
      TCR0 = TCR
*
*       Obtain approximate half-mass radius after scaling.
      RSCALE = 0.5*ZMASS**2/(SX*POT)
*       Set square radius of reflecting sphere.
      RSPH2 = (RSPH2*RSCALE)**2
*       Form equilibrium rms velocity (temporarily defined as VC).
      VC = SQRT(2.0D0*ABS(E0)/ZMASS)
*
C*       Check for general binary search of initial condition.
C      IF (KZ(4).GT.0) THEN
C          CALL EVOLVE(0,0)
C      END IF
*
*       Print half-mass relaxation time & equilibrium crossing time.
      A1 = FLOAT(N)
      TRH = 4.0*TWOPI/3.0*(VC*RSCALE)**3/(15.4*ZMASS**2*LOG(A1)/A1)
      if(rank.eq.0)WRITE (6,95)  TRH, TCR, 2.0*RSCALE/VC
   95 FORMAT (/,12X,'TIME SCALES:   TRH =',1PE8.1,'  TCR =',E8.1,
     &                                            '  2<R>/<V> =',E8.1,/)
*
      RETURN
*
      END
