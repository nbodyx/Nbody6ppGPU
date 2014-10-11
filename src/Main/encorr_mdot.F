      subroutine encorr_mdot(NPNUM,NDMLIST,DMLIST,GPUDM,GPUSW)
*
*     
*     Correct the energy due to mass loss
*     -----------------------------------

      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      INCLUDE 'omp_lib.h'
      INTEGER NDMLIST(NMAX),NPNUM
*    &  ,NDCLIST(NMAX/10),NCNUM
      REAL*8 DMLIST(NMAX),DMPHI(NMAX)
#ifdef PARALLEL
      REAL*8 DMPHI_L(NMAX)
#endif
      LOGICAL GPUDM,GPUSW
*     &  ,DCLIST(NMAX/10)
      
*     --11/11/13 18:27-lwang-check--------------------------------------*
***** Note: Energy correction for mass loss stars----------------------**
      if (NPNUM.GT.0) THEN
#ifdef NGPU
         call cputim(tttpota)
         NNT = NTOT - IFIRST + 1
*     Check whether need to send all particles to gpu
         if(GPUSW) then
            call gpupot_send(rank,NNT,BODY(IFIRST),X(1,IFIRST))
            GPUSW = .false.
         end if
*     Calculation potential
         IF (NPNUM.LE.2048) THEN
            IF(GPUDM) THEN 
               CALL GPUPOT_DM(rank,NPNUM,NNT,1-IFIRST,NDMLIST(1),
     *              DMLIST(1),DMPHI(1))
            ELSE
               CALL GPUPOT_FLOAT(rank,NPNUM,NNT,1-IFIRST,NDMLIST(1),
     *              DMLIST(1),DMPHI(1))
            END IF
c$$$               if(rank.eq.0) print*, 'GPU:',DMPHI(1:NPNUM)
c$$$            CALL GPUPOT(rank,NDMLIST(1),NDMLIST(1),NNT,BODY(IFIRST),
c$$$     *           X(1,IFIRST),DMPHI(1))
c$$$               if(rank.eq.0) print*, 'GPU2:',DMPHI(1:NPNUM)
         else
            DO J=1,NPNUM,2048
               NN = 2048
               IF (NPNUM-J+1.LT.2048) NN = NPNUM-J+1
               IF(GPUDM) THEN
                  CALL GPUPOT_DM(rank,NN,NNT,1-IFIRST,NDMLIST(J),
     *                 DMLIST(J),DMPHI(J))
               ELSE
                  CALL GPUPOT_FLOAT(rank,NN,NNT,1-IFIRST,NDMLIST(J),
     &                 DMLIST(J),DMPHI(J))
               END IF
            END DO
         end if
         call cputim(tttpotb)
         ttpotg = ttpotg + (tttpotb-tttpota)*60
#else
         call cputim(tttpota)
         DMPHI(1:NPNUM) = 0.0
#ifdef PARALLEL
         i_tot = NTOT - IFIRST + 1
         if(i_tot*NPNUM.lt.10000*icore*isize) then
#endif
            DO III = 1, NPNUM
               DMPHI_T = 0.D0
!$omp parallel do private(J,L,RIJ2)
!$omp& reduction(+:DMPHI_T)
               DO J = IFIRST,NTOT
                  L = NDMLIST(III)
                  IF (J.NE.L) THEN
                     RIJ2 = 0.0D0
                     DO  K = 1,3
                        RIJ2 = RIJ2 + (X(K,L) - X(K,J))**2
                     END DO
                     DMPHI_T = DMPHI_T + BODY(J)/SQRT(RIJ2)
                  END IF
               END DO
!$omp end parallel do
               DMPHI(III) = DMPHI_T
            END DO
#ifdef PARALLEL
         else
            DMPHI_L(1:NPNUM) = 0.0
            i_block = i_tot/isize
            if(mod(i_tot,isize).ne.0) i_block = i_block + 1
            i_start = rank*i_block + ifirst
            i_end = (rank+1)*i_block + ifirst - 1
            if(rank.eq.isize-1) i_end = NTOT
            DO III = 1, NPNUM
               DMPHI_T = 0.D0
!$omp parallel do private(J,L,RIJ2) reduction(+:DMPHI_T)
               DO J = i_start,i_end
                  L = NDMLIST(III)
                  IF (J.NE.L) THEN
                     RIJ2 = 0.0D0
                     DO  K = 1,3
                        RIJ2 = RIJ2 + (X(K,L) - X(K,J))**2
                     END DO
                     DMPHI_T = DMPHI_T + BODY(J)/SQRT(RIJ2)
                  END IF
               END DO
!$omp end parallel do
               DMPHI_L(III) = DMPHI_T
            END DO
            CALL MPI_ALLREDUCE(DMPHI_L,DMPHI,NPNUM,MPI_REAL8,MPI_SUM,
     &           MPI_COMM_WORLD,ierr)
         end if
#endif
c$$$            DO IIJ = 1, NPNUM
c$$$               J = NDMLIST(IIJ)
c$$$               DO III = 1, NPNUM
c$$$                  IF (IIJ.NE.III) THEN
c$$$                     L = NDMLIST(III)
c$$$c$$$                     IF(J.EQ.L) print*,'Warning!: J.eq.L',J,L,IIJ,III
c$$$                     RIJ2 = 0.0D0
c$$$                     DO  K = 1,3
c$$$                        RIJ2 = RIJ2 + (X(K,L) - X(K,J))**2
c$$$                     END DO
c$$$                    DMPHI(IIJ) = DMPHI(IIJ) + 0.5*DMLIST(III)/SQRT(RIJ2)
c$$$                  END IF
c$$$               END DO
c$$$            END DO
c$$$            IF (NCNUM.GT.0) THEN
c$$$               DO IIJ = 1, NPNUM
c$$$                  J = NDMLIST(IIJ)
c$$$                  DO III = 1, NCNUM
c$$$                     L = NDCLIST(III)
c$$$                     IF(J.NE.L) THEN
c$$$                        RIJ2 = 0.0D0
c$$$                        DO  K = 1,3
c$$$                           RIJ2 = RIJ2 + (X(K,L) - X(K,J))**2
c$$$                        END DO
c$$$                        DMPHI(IIJ) = DMPHI(IIJ) + 
c$$$     &                       0.5*DCLIST(III)/SQRT(RIJ2)
c$$$                     END IF
c$$$                  END DO
c$$$               END DO
c$$$            END IF
         call cputim(tttpotb)
         ttpot = ttpot + (tttpotb-tttpota)*60
*     --05/20/14 22:10-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*         if(NPNUM.GE.2) then
c$$$         print*,rank,'NPNUM',NPNUM,'DMPHI',DMPHI(1:NPNUM),'ttpot',ttpot
c$$$         call flush(6)
c$$$         stop
*         end if
*     --05/20/14 22:10-lwang-end----------------------------------------*
c$$$            if(rank.eq.0)print*, 'HOST:',NPNUM,DMPHI(1:NPNUM)
c$$$            call flush(6)
c$$$            stop
#endif
         call cputim(tttpota)
*!$omp parallel do private(J,L,VI2) reduction(+:EMDOT)         
         DO J = 1, NPNUM
            L = NDMLIST(J)
            VI2 = XDOT(1,L)**2 + XDOT(2,L)**2 + XDOT(3,L)**2
            EMDOT = EMDOT + DMLIST(J)*(0.5*VI2 - DMPHI(J))
         END DO
*!$omp end parallel do
*     --11/15/13 17:06-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         if(rank.eq.0) then
c$$$            print*,'EMDOT',EMDOT,TIME,NDMLIST(1:NPNUM),DMPHI(1),VI2
c$$$         end if
*     --11/15/13 17:06-lwang-end----------------------------------------*
         NPNUM = 0
         call cputim(tttpotb)
         ttecor = ttecor + (tttpotb-tttpota)*60
      END IF
*     --11/11/13 18:27-lwang-end----------------------------------------*

      RETURN
      
      END
