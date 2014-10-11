      subroutine ksparmpi(operator,par_type,par_index,
     &     index2,index3,par_value)
*
*
*     Save and communicate variables inside parallel KS integration
*
#ifdef PARALLEL
      Include 'kspars.h'
      Include 'common6.h'
       COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/MODES/   EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &                BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &                RP(NTMAX),ES(NTMAX),CMM(2,NTMAX),IOSC(NTMAX),
     &                NAMEC(NTMAX)

      INTEGER LISTILEN,LISTRLEN,par_type,par_index,index2,index3
      INTEGER LISTILEN_MPI(maxpe),LISTRLEN_MPI(maxpe)
      INTEGER LIDISPLACE(maxpe),LRDISPLACE(maxpe)
      INTEGER OPERATOR
      INTEGER Ibegin,Iend,Icounter
      INTEGER I_flag,I_index,I_index2,I_index3
      INTEGER I_LAST(3)
      REAL*8 par_value,par_last
      INTEGER IKSMPI(2*KMAX+10),IKSMERGE(2*KMAX+10)
      REAL*8  RKSMPI(2*KMAX+10),RKSMERGE(2*KMAX+10)
      REAL*8  EMERGE_MPI,BE3_MPI,ETIDE_MPI,ECOLL_MPI,EGRAV_MPI,E10_MPI
      REAL*8  ECDOT_MPI,EKICK_MPI
      REAL*8  EMERGEOLD,BEOLD,ETIDEOLD,ECOLLOLD,EGRAVOLD,EOLD10
      REAL*8  ECDOTOLD,EKICKOLD
      SAVE IKSMPI,RKSMPI,LISTILEN,LISTRLEN
      SAVE EMERGEOLD,BEOLD,ETIDEOLD,ECOLLOLD,EGRAVOLD,EOLD10
      SAVE ECDOTOLD,EKICKOLD,I_LAST,par_last
      DATA I_LAST /3*-1/

*     --02/24/14 19:03-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(operator.eq.K_store) then
c$$$         print*,'Rank:',rank,'CALL KSPARMPI: TYPE',par_type,'Index ',
c$$$     &        par_index,index2,index3,'Value ',par_value,'t',time
c$$$         call flush(6)
c$$$      elseif(operator.eq.K_comm)then
c$$$         print*,'Rank:',rank,'CALL KSPARMPI COMM T:',time
c$$$         call flush(6)
c$$$      else
c$$$         print*,'Rank:',rank,'CALL KSPARMPI RESET T:',time
c$$$         call flush(6)
c$$$      end if
*     --02/24/14 19:03-lwang-end----------------------------------------*
*     save data
      IF (operator.eq.K_store) then
*     Check whether the data is dup.
         if(par_index.eq.I_Last(1).and.index2.eq.I_last(2)
     &        .and.index3.eq.I_Last(3).and.par_last.eq.par_value)
     &        return
*     Backup data
         I_LAST(1) = par_index
         I_LAST(2:3) = 0
         par_last = par_value
*     integer variable
         IF (par_type.eq.K_int) then
*     save variable index
            LISTILEN = LISTILEN + 1
            IKSMPI(LISTILEN) = par_index
*     check one dimensional array
            IF(par_index.ge.IISPLIT) then
               LISTILEN = LISTILEN + 1
               IKSMPI(LISTILEN) = index2
               I_LAST(2) = index2
            end if
*     save value
            LISTILEN = LISTILEN + 1
            IKSMPI(LISTILEN) = int(par_value)
*     real*8 variable
         ELSEIF (par_type.eq.K_real8) then
*     real*8 variable
            LISTRLEN = LISTRLEN + 1
*     --04/08/14 11:39-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            if(LISTRLEN.GE.2*KMAX+10) then
c$$$               write(130+rank,*) 'LEN',LISTRLEN,'RKS',RKSMPI(1:LISTRLEN)
c$$$               call flush(130+rank)
c$$$            end if
*     --04/08/14 11:39-lwang-end----------------------------------------*
            RKSMPI(LISTRLEN) = DFLOAT(par_index)
*     check one dimensional array
            IF(par_index.ge.IRSPLIT) then
               LISTRLEN = LISTRLEN + 1
               RKSMPI(LISTRLEN) = DFLOAT(index2)
               I_LAST(2) = index2
*     check two dimensional array
               IF(par_index.ge.IRSPLIT2) then
                  LISTRLEN = LISTRLEN + 1
                  RKSMPI(LISTRLEN) = DFLOAT(index3)
                  I_LAST(3) = index3
               END IF
            END IF
*     save value
            LISTRLEN = LISTRLEN + 1
            RKSMPI(LISTRLEN) = par_value
         END IF

*     communication
      ELSEIF (operator.eq.K_comm) then
*     Gather all integer data
         call MPI_ALLGATHER(LISTILEN,1,MPI_INTEGER,LISTILEN_MPI(1),
     &        1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_ALLGATHER(LISTRLEN,1,MPI_INTEGER,LISTRLEN_MPI(1),
     &        1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         LIDISPLACE(1) = 0
         LRDISPLACE(1) = 0
         DO K=1,isize-1
            LIDISPLACE(K+1) = LIDISPLACE(K) + LISTILEN_MPI(K)
            LRDISPLACE(K+1) = LRDISPLACE(K) + LISTRLEN_MPI(K)
         END DO
*     Integer part:
         IF(LIDISPLACE(isize).GE.isize.OR.LISTILEN_MPI(isize).GT.1) THEN
            IKSMPI(1) = LISTILEN
            CALL MPI_ALLGATHERV(IKSMPI(1),LISTILEN,MPI_INTEGER,
     *           IKSMERGE(1),LISTILEN_MPI,LIDISPLACE,MPI_INTEGER,
     &           MPI_COMM_WORLD,ierr)
*     --02/25/14 16:39-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            LASTI = LIDISPLACE(isize)+LISTILEN_MPI(isize)
c$$$            print*,'rank',rank,'T',time,'IKSMPI',IKSMERGE(1:LASTI)
c$$$            call flush(6)
*     --02/25/14 16:39-lwang-end----------------------------------------*
*     Initial the flag and index
*     I_flag: 0: read index, 1: read value, 2: read index2, 3: read index3
            I_FLAG = 0
            I_index = 0
            I_index2 = 0
            I_index3 = 0
*     For MPI processors loop
            Icounter = 0
            IEND = 0
*     LOOP
 1          IBEGIN = IEND + 2
            IEND = IKSMERGE(IEND+1) + IEND
            Icounter = Icounter + 1
            DO I=IBEGIN, IEND
               IF(I_flag.eq.0) then
                  I_index = IKSMERGE(I)
                  IF(I_index.GE.IISPLIT) then
                     I_FLAG = 2
                  ELSE
                     I_FLAG = 1
                  END IF
               ELSEIF(I_flag.eq.1) then
                  if(I_index.eq.K_DELAY) then
                     call delay_bcast(icounter-1)
                  elseif(I_index.eq.K_JCLOSE) then
                     JCLOSE = IKSMERGE(I)
                  elseif(I_index.eq.K_JCMAX) then
                     JCMAX = IKSMERGE(I)
                  elseif(I_index.eq.K_IQCOLL) then
                     IQCOLL = IKSMERGE(I)
                  elseif(I_index.eq.K_IPHASE) then
                     IPHASE = IKSMERGE(I)
                  elseif(I_index.eq.K_KSPAIR) then
                     KSPAIR = IKSMERGE(I)
                  elseif(I_index.eq.K_KICK) then
                     call kick(icounter-1,-147,0,0.0)
                  elseif(I_index.eq.K_KSTAR) then
                     KSTAR(I_index2) = IKSMERGE(I)
                  else
                     PRINT*,'Error! No Integer index ',I_index,
     *                    ' For KS MPI communication'
                     call flush(6)
                     call abort()
                  end if
                  I_flag = 0
               ELSEIF(I_flag.eq.2) then
                  I_index2 = IKSMERGE(I)
                  I_FLAG = 1
               ELSE
                  PRINT*,'Error flag for reading data in ksparmpi ('
     *                 ,I_FLAG,')!'
                  call flush(6)
                  call abort()
               END IF
*     --02/25/14 17:00-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$               print*,'rank',rank,'I',i,'FLAG',i_flag,'index',I_index,
c$$$     &              I_index2,I_index3,'counter',icounter
c$$$               call flush(6)
*     --02/25/14 17:00-lwang-end----------------------------------------*
            END DO
            IF (Icounter.LT.isize) GO TO 1
         END IF
*     REAL*8 part
         IF(LRDISPLACE(isize).GE.isize.OR.LISTRLEN_MPI(isize).GT.1) THEN
            RKSMPI(1) = DFLOAT(LISTRLEN)
            CALL MPI_ALLGATHERV(RKSMPI(1),LISTRLEN,MPI_REAL8,
     *           RKSMERGE(1),LISTRLEN_MPI,LRDISPLACE,MPI_REAL8,
     &           MPI_COMM_WORLD,ierr)
*     --02/25/14 16:39-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            LASTI = LRDISPLACE(isize)+LISTRLEN_MPI(isize)
c$$$            print*,'rank',rank,'T',time,'RKSMPI',RKSMPI(1:LISTRLEN)
c$$$            print*,'rank',rank,'T',time,'RKSMERGE',RKSMERGE(1:LASTI),
c$$$     &           'IRSPLIT',IRSPLIT
c$$$            call flush(6)
*     --02/25/14 16:39-lwang-end----------------------------------------*
*     Initial energy
            EMERGE_MPI = 0
            BE3_MPI = 0
            ETIDE_MPI = 0
            ECOLL_MPI = 0
            EGRAV_MPI = 0
            E10_MPI = 0
            ECDOT_MPI = 0
            EKICK_MPI = 0
            
*     Initial the flag and index
*     I_flag: 0: read index, 1: read value, 2: read index2, 3: read index3
            I_FLAG = 0
            I_index = 0
            I_index2 = 0
            I_index3 = 0
*     For MPI processors loop
            Icounter = 0
            IEND = 0
*     LOOP
 2          IBEGIN = IEND + 2
            IEND = INT(RKSMERGE(IEND+1)) + IEND
            Icounter = Icounter + 1
            DO I=IBEGIN, IEND
               IF(I_flag.eq.0) then
                  I_index = Int(RKSMERGE(I))
                  IF(I_index.GE.IRSPLIT) then
                     I_FLAG = 2
                     IF(I_index.GE.IRSPLIT2) I_FLAG = 3
                  ELSE
                     I_FLAG = 1
                  END IF
               ELSEIF(I_flag.eq.1) then
                  if(I_index.eq.K_EMERGE) then
                     EMERGE_MPI = EMERGE_MPI+RKSMERGE(I)
                  elseif(I_index.eq.K_BE3) then
                     BE3_MPI = BE3_MPI+RKSMERGE(I)
                  elseif(I_index.eq.K_ETIDE) then
                     ETIDE_MPI = ETIDE_MPI+RKSMERGE(I)
                  elseif(I_index.eq.K_ECOLL) then
                     ECOLL_MPI = ECOLL_MPI+RKSMERGE(I)
                  elseif(I_index.eq.K_EGRAV) then
                     EGRAV_MPI = EGRAV_MPI+RKSMERGE(I)
                  elseif(I_index.eq.K_E10) then
                     E10_MPI = E10_MPI+RKSMERGE(I)
                  elseif(I_index.eq.K_EBCH0) then
                     EBCH0 = RKSMERGE(I)
                  elseif(I_index.eq.K_ECDOT) then
                     ECDOT_MPI = ECDOT_MPI+RKSMERGE(I)
                  elseif(I_index.eq.K_EKICK) then
                     EKICK_MPI = EKICK_MPI+RKSMERGE(I)
                  elseif(I_index.eq.K_RMAX) then
                     RMAX = MAX(RMAX,RKSMERGE(I))
                  elseif(I_index.eq.K_TEV) then
                     TEV(I_index2) = RKSMERGE(I)
                  elseif(I_index.eq.K_RADIUS) then
                     RADIUS(I_index2) = RKSMERGE(I)
                  elseif(I_index.eq.K_SPIN) then
                     SPIN(I_index2) = RKSMERGE(I)
                  elseif(I_index.eq.K_STEPS) then
                     STEPS(I_index2) = RKSMERGE(I)
                  elseif(I_index.eq.K_TMDIS) then
                     TMDIS(I_index2) = RKSMERGE(I)
                  elseif(I_index.eq.K_X0DOT) then
                     X0DOT(I_index3,I_index2) = RKSMERGE(I)
                  elseif(I_index.eq.K_CM_B) then
                     CM(I_index3,I_index2) = RKSMERGE(I)
                  elseif(I_index.eq.K_CM_M) then
                     CMM(I_index3,I_index2) = RKSMERGE(I)
                  else
                     PRINT*,'Error! No Integer index ',I_index,
     *                    ' For KS MPI communication'
                     call flush(6)
                     CALL abort()
                  end if
                  I_flag = 0
               ELSEIF(I_flag.eq.2) then
                  I_index2 = int(RKSMERGE(I))
                  I_FLAG = 1
               ELSEIF(I_flag.eq.3) then
                  I_index3 = int(RKSMERGE(I))
                  I_FLAG = 2
               ELSE
                  PRINT*,'Error flag for reading data in ksparmpi ('
     *                 ,I_FLAG,')!'
                  call flush(6)
                  CALL abort()
               END IF
*     --02/25/14 17:00-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$               print*,'rank',rank,'I',i,'FLAG',i_flag,'index',I_index,
c$$$     &              I_index2,I_index3,'counter',icounter
c$$$               call flush(6)
*     --02/25/14 17:00-lwang-end----------------------------------------*
            END DO
            IF (Icounter.LT.isize) GO TO 2
*     Update energy
            IF(EMERGE_MPI.NE.0) EMERGE = EMERGEOLD + EMERGE_MPI
*     --03/05/14 21:26-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            if(EMERGE_MPI.GE.5) print*,'Warning! EMERGE_MPI Error'
*     --03/05/14 21:26-lwang-end----------------------------------------*
            IF(BE3_MPI.NE.0) BE(3) = BEOLD + BE3_MPI
            IF(ETIDE_MPI.NE.0) ETIDE = ETIDEOLD + ETIDE_MPI
            IF(ECOLL_MPI.NE.0) ECOLL = ECOLLOLD + ECOLL_MPI
            IF(EGRAV_MPI.NE.0) EGRAV = EGRAVOLD + EGRAV_MPI
            IF(E10_MPI.NE.0) E(10) = EOLD10 + E10_MPI
            IF(ECDOT_MPI.NE.0) ECDOT = ECDOTOLD + ECDOT_MPI
            IF(EKICK_MPI.NE.0) EKICK = EKICKOLD + EKICK_MPI
         END IF
*     RESET
      ELSEIF (operator.eq.K_reset) then
         LISTILEN = 1
         LISTRLEN = 1
         EMERGEOLD = EMERGE
         BEOLD = BE(3)
         ETIDEOLD = ETIDE
         ECOLLOLD = ECOLL
         EGRAVOLD = EGRAV
         EOLD10 = E(10)
         ECDOTOLD = ECDOT
         EKICKOLD = EKICK
         I_LAST(1:3) = -1
      END IF
#endif

      RETURN
      
      END

**************************************************************************      
*     For delay
      subroutine delay_bcast(iproc)
*
*
*     bcast delay common variables
*
#ifdef PARALLEL      
      include 'mpi_base.h'
      COMMON/SAVEIT/ IPH0,KS0,KCH0,KS20,JC0,JCL0,PCR0
      integer iproc

      CALL MPI_BCAST(IPH0,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KS0,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KCH0,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KS20,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(JC0,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(JCL0,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(PCR0,1,MPI_REAL8,iproc,MPI_COMM_WORLD,ierr)
#endif
      RETURN

      END
      
