      SUBROUTINE IMBHINIT
*
*      
*     Initialization of massive black hole
*
      INCLUDE 'common6.h'
      REAL*8 XIMBH(3),VIMBH(3)
      
*     Read massive black hole mass and output time interval
      if(rank.eq.0) READ (5,*) BIMBH, XIMBH(1:3),VIMBH(1:3), DTBH

#ifdef PARALLEL
      CALL MPI_BCAST(BIMBH,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(XIMBH,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(VIMBH,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DTBH,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif      

*     Reset mass of particle N, treat N as IMBH
      BODY(N) = BIMBH/ZMBAR
      KSTAR(N) = 14
      TEV(N) = MAX(2*TCRIT,1000000.0)
      DO K = 1,3
         X(K,N) = XIMBH(K)/RBAR
         XDOT(K,N) = VIMBH(K)/VSTAR
      END DO
*     Store the name of massive black hole
      NIMBH = N

*     Printing initial massive black hole data
      IF(rank.eq.0) write(6,1) NIMBH,BIMBH,XIMBH(1:3),
     &     VIMBH(1:3)
 1    format('MBH Initialization:  NAME',I12,' M[M*]',E15.5,
     &     ' X(1:3)[PC]',3E15.5,' V(1:3)[km/s]',3E15.5)
      
      RETURN
      
      END

***** Note:------------------------------------------------------------**
      
      SUBROUTINE IMBHRESTART
*
*      
*     Initialization of massive black hole
*
      INCLUDE 'common6.h'

*     Search index of massive black hole
      IBH = 0
      DO K = N,1,-1
         IF (NAME(K).EQ.NIMBH) THEN
            IBH = K
            GO TO 10
         END IF
      END DO

      if(IBH.eq.0) then
         if(rank.eq.0) print*,'BH not found! N',NIMBH
         call flush(6)
         call abort
      end if

 10   BIMBH = BODY(IBH)*ZMBAR

*     Print massive black hole data
      IF(rank.eq.0) write(6,1) IBH,NIMBH,BIMBH,X(1:3,IBH)*RBAR,
     &     XDOT(1:3,IBH)*VSTAR
 1    format('MBH: IBH',I12,' NAME',I12,' M[M*]',E15.5,
     &     ' X(1:3)[PC]',3E15.5,' V(1:3)[km/s]',3E15.5)
      call flush(6)

      RETURN

      END
      
