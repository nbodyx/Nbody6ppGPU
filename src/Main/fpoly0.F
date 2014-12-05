      SUBROUTINE FPOLY0(RS0)
*
*
*       GPU initialization of neighbour lists and forces.
*       -------------------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      INCLUDE 'omp_lib.h'
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
      PARAMETER (NIMAX=1024)
      REAL*8   GPUACC(3,NIMAX),GPUJRK(3,NIMAX)
#ifdef GPU
      REAL*8   H2I(NIMAX),GPUPHI(NIMAX),DTR(NIMAX)
#endif
      REAL*8   DV(3)
      INTEGER  LISTGP(LMAX,NIMAX)
      SAVE  LISTGP
      INTEGER  M_FLAG
#ifdef PARALLEL
      integer inum(NIMAX),ista(NIMAX)
#endif
#ifndef GPU
      integer IREG(NMAX)
#endif
*
*
#ifdef GPU
*       Open the GPU libraries on new run (note nnbmax = NN is printed).
      NN = NTOT + 10
      CALL GPUNB_OPEN(NN,rank)
#endif

      M_FLAG = 0
      NB_FLAG = 0
*     Use mass ratio as neighbor member criterion
      IF(KZ(39).EQ.2.OR.KZ(39).EQ.3) M_FLAG = 1
      IF(KZ(39).EQ.0.OR.KZ(39).EQ.2) NB_FLAG = 1

*     Set larger value for GPUIRR (note further possible increase of NTOT).
*      NNN = N + 10
*      CALL GPUIRR_OPEN(NNN,LMAX)
*
*       Set provisional neighbour radius and regular time-step for GPUNB.
!$omp parallel do private(I,K,RI2)
      DO 10 I = IFIRST,NTOT
#ifndef GPU
          IREG(I-IFIRST+1) = I
#endif
          T0(I) = 0.0
*       Modify neighbour radius according to NBLIST procedure.
          RI2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
          RS(I) = RS0*SQRT(1.0 + RI2)
*       Set an estimated r-dependent regular time-step.
          STEPR(I) = SMAX/8.0D0*SQRT(1.0 + RI2)
          STEPR(I) = MIN(STEPR(I),SMAX)
          DO 5 K = 1,3
              F(K,I) = 0.0
              FDOT(K,I) = 0.0
    5     CONTINUE
   10 CONTINUE
!$omp end parallel do
*
*       Send all particles (X0 & X0DOT) to the GPU for prediction.
**!$omp parallel do private(I)
*      DO 15 I = 1,N
*          CALL GPUIRR_SET_JP(I,X(1,I),XDOT(1,I),F(1,I),FDOT(1,I),
*     &                                          BODY(I),T0(I))
*   15 CONTINUE
**!$omp end parallel do
*
      NN = NTOT - IFIRST + 1
*
#ifdef GPU
*       Send all single particles to the GPU.
      CALL GPUNB_SEND(NN,BODY(IFIRST),X(1,IFIRST),XDOT(1,IFIRST))
*
*       Define maximum GPU neighbour number and initialize counters.
      NBMAX = MIN(NNBMAX + 150,LMAX-5)
      NOFL2 = 0
#endif
*
*       Loop over all particles split into NIMAX blocks.
#ifdef PARALLEL
      nl = nn
      inl = nl/isize
      jsize = isize*inl
      idiff = nl - jsize
      irun = 0
*
      do 1103 ix = 1,isize
      inum(ix)=inl
      if(ix.le.idiff)inum(ix) = inum(ix) + 1
      ista(ix) = irun+1
      if(ista(ix).gt.nl)inum(ix) = 0
 1103 irun = irun + inum(ix)
*
      istart = ista(rank+1)
      iend = ista(rank+1) + inum(rank+1) - 1
      JNEXT = istart - 1      
#else
      JNEXT = 0
      istart = 1
      iend = NN
#endif 
      DO 100 II = istart,iend,NIMAX
         NI = MIN(iend-JNEXT,NIMAX)
*       Copy neighbour radius, STEPR and state vector for each block.
#ifdef GPU
*!$omp parallel do private(LL, I)
          DO 30 LL = 1,NI
              I = JNEXT + LL + IFIRST - 1
              IF(M_FLAG.EQ.1) THEN
*     Use mass ratio as neighbor member criterion
                 H2I(LL) = RS(I)**2/BODYM
              ELSE 
                 H2I(LL) = RS(I)**2
              END IF
              DTR(LL) = STEPR(I)
   30     CONTINUE
*!$omp end parallel do
*
*       Evaluate forces, first derivatives and neighbour lists for new block.
 31       I = JNEXT + IFIRST
          CALL GPUNB_REGF(NI,H2I,DTR,X(1,I),XDOT(1,I),GPUACC,GPUJRK,
     &         GPUPHI,LMAX,NNBMAX,LISTGP,M_FLAG)
*       Check neighbour lists for overflow or zero membership (NNB = 1).
*!$omp parallel do private(LL,NNB,I,RI2) reduction(+:NOFL2)
          DO 50 LL = 1,NI
              NNB = LISTGP(1,LL)
*       Repeat last block with reduced RS(I) on NNB < 2 (at end of loop).
              IF (NNB.LT.2) THEN
                  I = JNEXT + LL + IFIRST - 1
                  RI2 = (X(1,I)-RDENS(1))**2 + (X(2,I)-RDENS(2))**2 +
     &                                         (X(3,I)-RDENS(3))**2
                  IF (NNB.LT.0) THEN
                     WRITE (41,40)  NSTEPR, NAME(I), -NNB,
     &                    RS(I), SQRT(RI2)
 40                  FORMAT (' NB OVERFLOW! NSTEPR NAME NB RNB[NB] ',
     &                    'RI[NB] ',I11,I8,I7,2F8.3)
                     CALL FLUSH(41)
*     For smooth density system, use NNBOPT/NNB to estimate the neighbor radius
                     IF(NB_FLAG.EQ.1) THEN
                        IF(-NNB.GT.NNBOPT) THEN
                          RS(I)=(-FLOAT(NNBOPT)/FLOAT(NNB))**0.333*RS(I)
                        ELSE
                          RS(I)=(FLOAT(NNBOPT)/FLOAT(NNBMAX))**0.4*RS(I)
                        END IF
                     ELSE
                        RS(I) = 0.9*RS(I)
                     END IF
                  ELSE
                     WRITE (41,41)  NSTEPR, NAME(I), NNB,
     &                    RS(I), SQRT(RI2)
 41                  FORMAT (' ZERO NB!   #R NAME NB RS ri ',
     &                    I11,I7,I5,2F8.3)
                     CALL FLUSH(41)
                     RS(I) = 1.5*RS(I)
                  END IF
                  IF(M_FLAG.EQ.1) THEN
*     Use mass ratio as neighbor member criterion
                     H2I(LL) = RS(I)**2/BODYM
                  ELSE
                     H2I(LL) = RS(I)**2
                  END IF
                  NOFL2 = NOFL2 + 1
              END IF
   50     CONTINUE
*!$omp end parallel do
*
*       Repeat the last block for rare exceptions (NNB < 0 and = 1).
          IF (NOFL2.GT.0) THEN
             NOFL2 = 0
             GO TO 31
          END IF
*
*      Copy regular force and neighbour list from GPU.
!$omp parallel do private(LL, I, ITEMP, NNB, L1, L)
          DO 70 LL = 1,NI
              I = JNEXT + LL + IFIRST - 1
              DO 55 K = 1,3
                  FR(K,I) = GPUACC(K,LL)
                  D1R(K,I) = GPUJRK(K,LL)
   55         CONTINUE
              NNB = LISTGP(1,LL)
              L1 = 1
              DO 60 L = 2,NNB+1
*       Note GPU address starts from 0 (hence add IFIRST to neighbour list).
                  ITEMP = LISTGP(L,LL) + IFIRST
                  IF (ITEMP.NE.I) THEN
                      L1 = L1 + 1
                      LISTGP(L1,LL) = ITEMP
                  END IF
   60         CONTINUE
              LISTGP(1,LL) = L1 - 1
*              CALL GPUIRR_SET_LIST(I,LISTGP(1,LL))
              DO 65 L = 1,L1
                  LIST(L,I) = LISTGP(L,LL)
   65         CONTINUE
   70     CONTINUE
!$omp end parallel do
*
#else
*     OpemMP version of regular force, fdot, potential and neighbor list
          IOFF = JNEXT + 1
*     --11/30/14 16:46-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          print*,rank,'IOFF',IOFF,'NI',NI,'IS',ISTART,'IE',IEND,
c$$$     &         'IREG',IREG(IOFF:IOFF+2),'JNEXT',JNEXT
*     --11/30/14 16:46-lwang-end----------------------------------------*
          CALL regint(IOFF,NI,IREG,IFIRST,NTOT,X,XDOT,BODY,STEPR,BODYM,
     &         M_FLAG,NB_FLAG,NNBOPT,RS,GPUACC,GPUJRK,LISTGP,PHIDBL)
*     Check for zero neighbor
          NINEW = NI
 42       NICHECK = NINEW
          NINEW = 0
          DO LL = 1, NICHECK
             I = IREG(LL+IOFF-1)
             NNB = LISTGP(1,LL) + 1
             IF(NNB.EQ.1) THEN
                RI2 = (X(1,I)-RDENS(1))**2 + (X(2,I)-RDENS(2))**2 +
     &               (X(3,I)-RDENS(3))**2
                WRITE (41,41)  LL, I, NAME(I), LISTGP(1,LL),
     &               RS(I), SQRT(RI2)
 41             FORMAT (' ZERO NB!  L I NAME NB RS ri ',
     &               3I11,I5,2F8.3)
                CALL FLUSH(41)
                RS(I) = 1.5*RS(I)
                IREG(NINEW+IOFF) = I
                NINEW = NINEW + 1
             ELSE
                FR(1:3,I) = GPUACC(1:3,LL)
                D1R(1:3,I) = GPUJRK(1:3,LL)
                LIST(1:NNB,I) = LISTGP(1:NNB,LL)
             END IF
          END DO
          IF (NINEW.GT.0) THEN
*     --11/30/14 16:46-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$             print*,rank,'IOFF',IOFF,'NINEW',NINEW,
c$$$     &            'IREG',IREG(IOFF:2+IOFF)
*     --11/30/14 16:46-lwang-end----------------------------------------*
             CALL regint(IOFF,NINEW,IREG,IFIRST,NTOT,X,XDOT,BODY,STEPR,
     &            BODYM,M_FLAG,NB_FLAG,NNBOPT,RS,GPUACC,GPUJRK,LISTGP,
     &            PHIDBL)
             GO TO 42
          END IF
#endif

*       Evaluate current irregular forces
!$omp parallel do private(LL,I,NNB,L,K,A1,A2,A3,DV,RIJ2,DR2I,
!$omp& DR3I,DRDV)
          DO 334 LL = 1,NI
             I = JNEXT + LL + IFIRST - 1
             FI(1:3,I) = 0
             D1(1:3,I) = 0
             NNB = LIST(1,I) + 1
             DO L = 2,NNB
                K = LIST(L,I)
                A1 = X(1,K) - X(1,I)
                A2 = X(2,K) - X(2,I)
                A3 = X(3,K) - X(3,I)
                DV(1) = XDOT(1,K) - XDOT(1,I)
                DV(2) = XDOT(2,K) - XDOT(2,I)
                DV(3) = XDOT(3,K) - XDOT(3,I)
                RIJ2 = A1*A1 + A2*A2 + A3*A3
*     
                DR2I = 1.0/RIJ2
                DR3I = BODY(K)*DR2I*SQRT(DR2I)
                DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
                FI(1,I) = FI(1,I) + A1*DR3I
                FI(2,I) = FI(2,I) + A2*DR3I
                FI(3,I) = FI(3,I) + A3*DR3I

                D1(1,I) = D1(1,I) + (DV(1) - A1*DRDV)*DR3I
                D1(2,I) = D1(2,I) + (DV(2) - A2*DRDV)*DR3I
                D1(3,I) = D1(3,I) + (DV(3) - A3*DRDV)*DR3I
             END DO
 334      CONTINUE
!$omp end parallel do
*
          JNEXT = JNEXT + NI
  100 CONTINUE

*       Check option for external force.
      IF (KZ(14).GT.0) THEN
          CALL XTRNLD(istart+ifirst-1,iend+ifirst-1,1)
      END IF
*
*
#ifdef PARALLEL
      isend = rank + 1
      if(isend.eq.isize)isend = 0
      irecv = rank - 1
      if(irecv.eq.-1)irecv = isize - 1
*
      do 1002 ir = 0,isize-2
*
      irank = rank - ir
      if(irank.lt.0)irank=irank+isize
*
      istsen=ista(irank+1)+ifirst-1
      icnt = inum(irank+1)

      if(irank.eq.0)irank=isize
      istrec = ista(irank)+ifirst-1
      icnt2 = inum(irank)

c$$$      print*, 'rank ',rank,' istsen ',istsen,' istrec-1 ',istrec-1,
c$$$     &     ' before fi(1) ',fi(1,istsen), ' fi(e) ',fi(1,istrec-1),
c$$$     &     ' size',icnt,' isend',isend,' irecv',irecv
*

#ifdef PUREMPI
      call cputim(tta)
      CALL MPI_SENDRECV(FI(1,istsen),3*icnt,MPI_REAL8,isend,
     *         rank,FI(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(FR(1,istsen),3*icnt,MPI_REAL8,isend,
     *         rank,FR(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(D1(1,istsen),3*icnt,MPI_REAL8,isend,
     *         rank,D1(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(D1R(1,istsen),3*icnt,MPI_REAL8,isend,
     *         rank,D1R(1,istrec),3*icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(RS(istsen),icnt,MPI_REAL8,isend,
     *         rank,RS(istrec),icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      CALL MPI_SENDRECV(LIST(1,istsen),LMAX*icnt,MPI_INTEGER,isend,
     *         rank,LIST(1,istrec),LMAX*icnt2,MPI_REAL8,irecv,irecv,
     *         MPI_COMM_WORLD,status,ierr)
      call cputim(ttb)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call cputim(tt999)
      ibarcount=ibarcount+1
      ttbar = ttbar + (tt999-ttb)*60
      xtsub2 = xtsub2 + dble((13*8+lmax*8)*(icnt+icnt2))
      ttsub2 = ttsub2 + (ttb-tta)*60.
c$$$      print*, 'rank ',rank,' istsen ',istsen,' istrec ',istrec,
c$$$     &     ' after fi(1) ',fi(1,istsen), ' fi(e) ',fi(1,istrec)
#endif
*
#ifdef SHMEM
       call barrier()
       call shmem_get(FI(1,istren),FI(1,istrec),3*icnt2,irecv)
       call shmem_get(FR(1,istren),FR(1,istrec),3*icnt2,irecv)
       call shmem_get(D1(1,istren),D1(1,istrec),3*icnt2,irecv)
       call shmem_get(D1R(1,istren),D1R(1,istrec),3*icnt2,irecv)
       call shmem_get(RS(1,istren),RS(1,istrec),icnt2,irecv)
       call shmem_get(LIST(1,istren),LIST(1,istrec),lmax*icnt2,irecv)
#endif
*
 1002 continue
#endif
      
*     Form total force & force derivative and extra variables for XVPRED.
!$omp parallel do private(I,K)
      DO 110 I = IFIRST,NTOT
          DO 105 K = 1,3
              F(K,I) = FI(K,I) + FR(K,I)
              FDOT(K,I) = D1(K,I) + D1R(K,I)
              D0(K,I) = FI(K,I)
              D0R(K,I) = FR(K,I)
              FIDOT(K,I) = D1(K,I)
              FRDOT(K,I) = D1R(K,I)
  105     CONTINUE
  110 CONTINUE
!$omp end parallel do
*

#ifdef GPU
*       Close the GPU libraries (limits change in INTGRT).
      CALL GPUNB_CLOSE
#endif
*      CALL GPUIRR_CLOSE
*
      RETURN
*
      END
