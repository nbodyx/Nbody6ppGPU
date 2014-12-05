      SUBROUTINE INPUT
*
*
*       Parameter input.
*       ----------------
*
      INCLUDE 'common6.h'
      EXTERNAL VERIFY
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
*       Make a formal call to define input parameters & counters.
C      CALL DEFINE
*
      IF(rank.eq.0)THEN
*       Read & print the main input parameters.
         READ (5,*)  N, NFIX, NCRIT, NRAND, NNBOPT, NRUN
C Termination time in physical units, TCRITP, read in nbody6.F
         READ (5,*)  ETAI, ETAR, RS0, DTADJ, DELTAT, TCRIT,
     &               QE, RBAR, ZMBAR
         READ (5,*)  (KZ(J),J=1,50)
         READ (5,*)  DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX, SMAX
      END IF

*     Check SMAX to make sure it have correct value
      if(rank.eq.0) THEN
         IF(SMAX.GT.1) THEN
            print*, 'Warning! SMAX > 1.0, reduce to 1.0.'
            SMAX = 1.0
         else
            DTNSMAX = 1.0
 1          IF(SMAX/DTNSMAX.LE.0.75) THEN
               DTNSMAX = 0.5D0*DTNSMAX
               IF(DTNSMAX.GT.1E-19) GO TO 1
            END IF
            SMAX = DTNSMAX
         END IF
      END IF

*
#if MPIINIT
      CALL MPI_BCAST(N,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NFIX,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NCRIT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NRAND,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NNBOPT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NRUN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
*
      CALL MPI_BCAST(KZ(1),50,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
*
      CALL MPI_BCAST(ETAI,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ETAR,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RS0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DTADJ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DELTAT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TCRIT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(QE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RBAR,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ZMBAR,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DTMIN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RMIN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      RMIN2 = RMIN**2 
      CALL MPI_BCAST(ETAU,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ECLOSE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(GMIN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(GMAX,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(SMAX,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)      
#endif
*

#ifdef ENSEMBLE
      NRAND = NRAND + 17*rank
      WRITE (6,11)rank,NRAND
   11 FORMAT (//,' Run configured for ensemble average PE ',I5,
     *   ' NRAND=',I9)
#endif

      if(rank.eq.0)then
         WRITE (6,10)
   10    FORMAT (/////,15X,'N  NFIX  NCRIT  NRAND  NNBOPT  NRUN')
         WRITE (6,12)  N, NFIX, NCRIT, NRAND, NNBOPT, NRUN
   12    FORMAT (/,I16,I6,2I7,I8,I6)
*
C New: (Aug.1998, P.Kroupa)
         WRITE(6,15)
   15    FORMAT (//,12X,' ETAI      ETAR      RS0       DTADJ',
     &                  '     DELTAT',
     &                  '     TCRITP    TCRIT     QE', 
     &                  '        RBAR       ZMBAR')
         WRITE (6,20)  ETAI, ETAR, RS0, DTADJ, DELTAT, TCRITP, TCRIT, 
     &              QE, RBAR,
     &              ZMBAR
   20    FORMAT (/,10X,1P10E10.1)
*
         WRITE (6,22)
   22    FORMAT (//,12X,'OPTIONS')
         WRITE (6,24)  (J,J=1,50)
   24    FORMAT (/,9X,50I3)
         WRITE (6,26)  (KZ(J),J=1,50)
   26    FORMAT (/,9X,50I3)
         WRITE (6,28)
   28    FORMAT (//,12X,'DTMIN     RMIN      ETAU      ECLOSE    GMIN',
     &        '      GMAX     SMAX')
         WRITE (6,30)  DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX, SMAX
   30    FORMAT (/,9X,1P7E10.1)
      end if
      call flush(6)
*
*       Define total particle number & neighbour membership range.
      NTOT = N
      NZERO = N
      NNBMAX = MIN(N/2,LMAX - 3)
      ZNBMIN = MAX(0.01*FLOAT(NNBMAX),1.0)
      ZNBMAX = 0.9*FLOAT(NNBMAX)
*       Save initial ETAI.
      ETA0 = ETAI
      RSMIN = RS0
      RC = RS0
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
      GPRINT(1) = 0.0
      DELTAS = 0.0
*     Suppress this KZ(4) since the output need special analysis tool
C      IF (KZ(4).GT.0) THEN
C*       Read parameters for binary evolution analysis.
C          K = KZ(4)
C          if(rank.eq.0)then
C          READ (5,*)  DELTAS, ORBITS(1), (GPRINT(J),J=1,K)
C          end if
C*
C#if MPIINIT
C      CALL MPI_BCAST(DELTAS,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C      CALL MPI_BCAST(ORBITS(1),9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C      CALL MPI_BCAST(GPRINT(1),9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C#endif
C*
C      if(rank.eq.0)WRITE (6,40)  DELTAS, ORBITS(1), (GPRINT(J),J=1,K)
C   40     FORMAT (//,12X,'DELTAS =',F6.2,'  ORBITS(1) =',F6.2,
C     &                                            '  GPRINT(J) =',9F7.3)
C*       Modify binary output factor by perturbation at different levels.
C          DO 50 L = 2,K
C              ORBITS(L) = ORBITS(1)*(GPRINT(1)/GPRINT(L))**0.3333
C   50     CONTINUE
C      END IF
*
C Old version:
*       Set random number skip for routine DATA.
c      IDUM1 = NRAND
C NEW version (14.08.98, P.Kroupa):
C*       Set random number SEED for routine DATA.
      IDUM1 = -1*NRAND
c+++ Notify others of this change on log file:
C      if(rank.eq.0)then
C      write(6,*)
C      write(6,*)' ****** NOTE: new random number seed initialisation!'
C      write(6,*)' ****** AND new ran2 from new ed. of Press et al.'
C      write(6,*)
C      end if
*
*
*       Save square of c.m. approximation parameter (coupled to GMIN).
      CMSEP2 = GMIN**(-0.666666667)
*
      RETURN
*
      END
