*     MPI_BASE

      INCLUDE 'mpif.h'
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      INTEGER icore,isernb,iserreg,iserks
      COMMON/MPIDAT/ group,rank,ierr,isize,status,icore,isernb,iserreg,
     &               iserks
