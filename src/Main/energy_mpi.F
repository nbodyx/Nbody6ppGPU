      SUBROUTINE ENERGY_MPI(iscale)
*     
*
*       Total energy.
*       -------------
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      INCLUDE 'omp_lib.h'
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
*     --11/10/13 23:32-lwang-check--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      logical ipstart
c$$$      INTEGER NPLIST(NMAX)
c$$$      data ipstart/.true./
c$$$      save ipstart,NPLIST
*     --11/10/13 23:32-lwang-end----------------------------------------*
#ifndef GPU
      REAL*8 SHPOT(NMAX)
#endif
      integer inum(maxpe),ista(maxpe)
      logical iscale
*
*     --11/11/13 12:56-lwang-check--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if (ipstart) then
c$$$         DO I = 1, NMAX
c$$$            NPLIST(I) = I
c$$$         END DO
c$$$         ipstart=.false.
c$$$      end if
*     --11/11/13 12:56-lwang-end----------------------------------------*
*       Sum the total energy of regularized pairs.
      if(.not.iscale) then
         EBIN = 0.0D0
         DO 1 IPAIR = 1,NPAIRS
*     Skip pairs with zero mass of c.m. particle (merged binary ghost).
            IF (BODY(N+IPAIR).GT.0.0D0) THEN
*     Predict coordinates, velocities & binding energy.
               CALL RESOLV(IPAIR,1)
               EBIN = EBIN + BODY(2*IPAIR-1)*BODY(2*IPAIR)*HT/
     &              BODY(N+IPAIR)
            END IF
 1       CONTINUE
      end if
*
*       Calculate the potential energy.
      ZKIN = 0.D00
      VIR = 0.0
*
      call cputim(tt998)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call cputim(tt999)
      ttbar = ttbar + (tt999-tt998)*60
      ibarcount=ibarcount+1

#ifdef GPU
*
      nl = ntot-ifirst+1
      if (nl.gt.5E4) then
#else
      nl = ntot
#endif         
*
      inl = nl/isize
      idiff = nl - isize*inl
      irun = 0
*
      do 1003 ix = 1,isize
      inum(ix)=inl
      if(ix.le.idiff)inum(ix) = inum(ix) + 1
      ista(ix) = irun+1
 1003 irun = irun + inum(ix)
*
      istart = ista(rank+1)
      iend = ista(rank+1) + inum(rank+1) - 1
*
#ifdef GPU      
*     --05/17/13 9:13-lwang-gpupot-------------------------------------*
***** Note:Use gpupot here---------------------------------------------**
*       Evaluate individual potentials on GPU (including all c.m.).
C      VIR = 0.0
      POT = 0.0
      NN = iend - istart + 1
      NNT = NTOT - IFIRST + 1
*     --11/10/13 23:30-lwang-check--------------------------------------*
***** Note:Loop for each 1024 particle---------------------------------**
c$$$      DO IJ = istart, iend, 1024
c$$$         if (iend-IJ+1.lt.1024) then
c$$$            NN = iend-ij+1
c$$$         else
c$$$            NN = 1024
c$$$         end if
c$$$      CALL GPUPOT(rank,ij,NN,NNT,BODY(IFIRST),X(1,IFIRST)
c$$$     &        ,phidbl(Ij+ifirst-1))
c$$$      end do
***** Note:Finish once-------------------------------------------------**      
c$$$      CALL GPUPOT(rank,NN,NNT,NPLIST(Istart),BODY(IFIRST),X(1,IFIRST),
c$$$     *     phidbl(istart+ifirst-1))
      CALL GPUPOT(rank,istart,NN,NNT,BODY(IFIRST),X(1,IFIRST)
     &     ,phidbl(Istart+ifirst-1))
*     --11/10/13 23:30-lwang-end----------------------------------------*
#else
!$omp parallel do private(I,JMIN,IPAIR,POTI,POTJ,A1,A2,A3,A4)
      DO 220 I = istart,iend
         JMIN = I + 1
         IF (I.LE.2*NPAIRS) THEN
*     Binding energy of regularized pairs is included explicitly above.
            IPAIR = KVEC(I)
            JMIN = 2*IPAIR + 1
         END IF
C     
         IPAIR = 0
         IF (I.GT.N)  THEN
*     Binding energy at center of mass position without binary members
            IPAIR = I - N
         END IF
*     
         POTJ = 0.D00
         POTI = 0.D00
*     POTI contains potential at particles position to be stored later (R.Sp.)
*     
         DO 30 J = 1,N
            IF (J.EQ.I .OR. J.EQ.2*IPAIR-1 .OR. J.EQ.2*IPAIR .OR.
     *           BODY(J).EQ.0.0D0 .OR. BODY(I).EQ.0.0D0)  GO TO 30
            A1 = X(1,I) - X(1,J)
            A2 = X(2,I) - X(2,J)
            A3 = X(3,I) - X(3,J)
            A4 = BODY(J)/DSQRT (A1*A1 + A2*A2 + A3*A3)
            POTI = POTI - A4
*     also J.LT.N?
            IF(J.GE.JMIN)POTJ = POTJ + A4
 30      CONTINUE
*     Store potential in shared vector first (R.Sp.)
         PHIDBL(I) = -POTI
         SHPOT(I) = BODY(I)*POTJ
 220  CONTINUE
!$omp end parallel do
#endif      
*
*        Distribute variables into private vectors again T3D (R.Sp.)
      
      isend = rank + 1
      if(isend.eq.isize)isend = 0
      irecv = rank - 1
      if(irecv.eq.-1)irecv = isize - 1
*
      do 1001 ir = 0,isize-2
*
      irank = rank - ir
      if(irank.lt.0)irank=irank+isize
*
#ifdef GPU      
      istart=ista(irank+1)+ifirst-1
      icnt = inum(irank+1)
*
      if(irank.eq.0)irank=isize
      istrec = ista(irank)+ifirst-1
      icnt2 = inum(irank)
#else
      istart=ista(irank+1)
      icnt = inum(irank+1)
*
      if(irank.eq.0)irank=isize
      istrec = ista(irank)
      icnt2 = inum(irank)
#endif      
*
*      print*,' ENERGY: bef rank,irank=',rank,irank,npairs
*      print*,' ENERGY: bef rank ',rank,' phidbl(',istrec,')=',
*     *          phidbl(istrec),NN,NNT,ifirst,iend-istart,icnt2
*      print*,' Send: phidbl(',istart,')=', phidbl(istart),icnt
*      print*,' ENERGY: bef rank ',rank,' shpot(',istrec,')=',
*     *          shpot(istrec)
*
*     call mpi_barrier(MPI_COMM_WORLD,ierr)
      call cputim(ttae)
      CALL MPI_SENDRECV(PHIDBL(istart),icnt,MPI_REAL8,isend,rank,
     *                  PHIDBL(istrec),icnt2,MPI_REAL8,irecv,irecv,
     *                  MPI_COMM_WORLD,status,ierr)
#ifndef GPU
      CALL MPI_SENDRECV(SHPOT(istart),icnt,MPI_REAL8,isend,rank,
     *                  SHPOT(istrec),icnt2,MPI_REAL8,irecv,irecv,
     *                  MPI_COMM_WORLD,status,ierr)
#endif      
      call cputim(ttab)
      ttsube = ttsube + (ttab-ttae)*60
*
      call cputim(tt998)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call cputim(tt999)
      ttbar = ttbar + (tt999-tt998)*60
      ibarcount=ibarcount+1

 1001  continue
#ifdef GPU       
      else
C         VIR = 0.0
         POT = 0.0
         NNT = NTOT - IFIRST + 1
*     --11/10/13 23:30-lwang-check--------------------------------------*
***** Note:Loop for each 1024 particle---------------------------------**
c$$$         DO IJ = 1, NNT, 1024
c$$$            if (NNT-IJ+1.lt.1024) then
c$$$               NN = NNT-ij+1
c$$$            else
c$$$               NN = 1024
c$$$            end if
c$$$     CALL GPUPOT(rank,ij,NN,NNT,BODY(IFIRST),X(1,IFIRST)
c$$$     &           ,phidbl(Ij+ifirst-1))
c$$$         end do
***** Note:Finish once-------------------------------------------------**
c$$$         CALL GPUPOT(rank,NNT,NNT,NPLIST(ifirst),BODY(IFIRST),
c$$$     *        X(1,IFIRST),phidbl(ifirst))
         CALL GPUPOT(rank,1,NNT,NNT,BODY(IFIRST),X(1,IFIRST)
     &        ,phidbl(ifirst))
c$$$         print*, 'phi: ',phidbl(ifirst:ifirst+5)
*     --11/11/13 12:53-lwang-end----------------------------------------*
      endif
*
*
*       Move the table entries down to give room for any KS components.
      I2 = 2*NPAIRS
      IF (NPAIRS.GT.0) THEN
*       Copy c.m. potential to the components.
          DO 10 IPAIR = 1,NPAIRS
              I1 = 2*IPAIR - 1
              PHIDBL(I1) = PHIDBL(N+IPAIR)
              PHIDBL(I1+1) = PHIDBL(N+IPAIR)
   10     CONTINUE
      END IF
*
*       Sum individual contributions after differential correction.
      DO 20 I = IFIRST,NTOT
          CALL PHICOR(I,DPHI1,DPHI2)
          IF (I.LE.N) THEN
              PHIDBL(I) = PHIDBL(I) + DPHI1
              POT = POT + BODY(I)*PHIDBL(I)
          ELSE
              I1 = 2*(I - N) - 1
              PHIDBL(I1) = PHIDBL(I1) + DPHI1
              PHIDBL(I1+1) = PHIDBL(I1+1) + DPHI2
              POT = POT + BODY(I1)*PHIDBL(I1) + BODY(I1+1)*PHIDBL(I1+1)
          END IF
   20 CONTINUE
*
*       Take half the value because of double counting.
      POT = 0.5*POT
#else
      POT = 0.D0
      DO 444 I = 1,NTOT
         POT = POT + SHPOT(I)
 444  CONTINUE
#endif      
*       Sum the kinetic energy (include c.m. bodies but not components).
      DO 40 I = IFIRST,NTOT
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   40 CONTINUE
*
      ZKIN = 0.5D0*ZKIN
*
*       Obtain the tidal potential energy for linearized external field. 
      IF (KZ(14).EQ.0) THEN
*       Note: EPL holds accumulated tidal energy if KZ(14) = 3.
         EPL = 0.0D0
      ELSE
         CALL XTRNLV(1,N)
*       Form tidal energy with Plummer potential (note EPL use for #14=3).
         IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
            EPL = 0.0
            DO 50 I = 1,N
               RI2 = AP2
               DO 45 K = 1,3
                  RI2 = RI2 + X(K,I)**2
 45            CONTINUE
               EPL = EPL - BODY(I)*MP/SQRT(RI2)
 50         CONTINUE
         END IF
      END IF
*
*       Check differential potential energy due to chain subsystem.
      IF (NCH.GT.0) THEN
          CALL CHPOT(DP)
          POT = POT + DP
      END IF
*
*       Total energy = ZKIN - POT + EPL + EBIN + ESUB + EMERGE + ECOLL.
*
      RETURN
*
      END
