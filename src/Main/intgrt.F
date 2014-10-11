      SUBROUTINE INTGRT
*     label-intgrt
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      include 'tlist.h'
      include 'omp_lib.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/BLKLVL/  JHIST(0:NMAX),JHISTR(0:NMAX),JHISTU(0:NMAX)
      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
      COMMON/XPRED/ TPRED(NMAX),TRES(KMAX),ipredall
      COMMON/NEWXV/ XN(3,NMAX),XNDOT(3,NMAX)
      INTEGER JHIST,JHISTR,LSHORT(NMAX)
      INTEGER IREG(NMAX)
#ifdef SIMD
      REAL*8 FIRR_SIMD(3,NMAX),FD_SIMD(3,NMAX)
#endif
#ifdef PARALLEL
      integer inum(maxpe),ista(maxpe),NXTLST_MPI(NMAX+500),NPIECE
      integer istart_lst,iend_lst
      REAL*8 XMPI(19,NMAX),YMPI(41,NMAX)
      SAVE NPIECE
      INTEGER IMPI(LMAX,NMAX),JMPI(11), JMPILOCAL(11,maxpe)
#endif
      LOGICAL LSTEPM,FNTLOOP,IPSTART,ipredall
      DATA LSHORT(1) /0/
      DATA IQ,ICALL,LSTEPM,STEPM,IPSTART /0,2,.FALSE.,0.03125,.true./
      SAVE IQ,ICALL,LSTEPM,STEPM,IPSTART,TMIN

*     --09/26/13 22:52-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      INTEGER NGPUC
c$$$      SAVE NGPUC
c$$$      DATA NGPUC/0/
C$$$      REAL*8 TCCT
C$$$      DATA TCCT/32.0/
C$$$      SAVE TCCT
*     --09/26/13 22:52-lwang-end----------------------------------------*

*
*       Update quantized value of STEPM for large N (first time only).
      call cputim(tttbegin)
      IF (.NOT.LSTEPM.AND.NZERO.GT.1024) THEN
         K = INT((FLOAT(NZERO)/1024.0)**0.333333)
         STEPM = 0.03125D0/2**(K-1)
         LSTEPM = .TRUE.
      END IF
*
*       Open GPU for all single particles
      IF (IPSTART) THEN
         FNTLOOP = .false.

*     Initial nxtlst
         IF(TIME+TOFF.EQ.0.D0) THEN
            NXTLIMIT = 0
            NGHOSTS = 0
            DO I = IFIRST,NTOT
               IF(BODY(I).NE.0.D0) THEN
                  NXTLIMIT = NXTLIMIT + 1
                  NXTLST(NXTLIMIT) = I
               ELSE
                  NXTLST(NMAX-NGHOSTS) = I
                  NGHOSTS = NGHOSTS + 1
               END IF
            END DO
            IF(NGHOSTS.GT.0) NXTLST(NXTLIMIT+1:NXTLIMIT+NGHOSTS) 
     &           = NXTLST(NMAX-NGHOSTS+1:NMAX)
*     --07/17/14 10:55-lwang-debug--------------------------------------*
***** Note: Print NXTLST PARAMETER-------------------------------------**
c$$$            if(rank.eq.0) write(6,*) 'NXTLIMIT',NXTLIMIT,'NGHOSTS',
c$$$     &           NGHOSTS
*     --07/17/14 10:56-lwang-end----------------------------------------*
*     Sorting NXTLST
            call sort_tlist(STEP,DTK,.true.)
            NXTLEVEL = NDTMAX
*     --07/08/14 10:21-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            print*,'NDTMAX',NDTMAX,'NDTMIN',NDTMIN
c$$$            print*,'NDTK',(NDTK(K),DTK(K),STEP(NXTLST(NDTK(K+1)+1)),
c$$$     &           K=1,NDTMAX)
*     --07/08/14 10:21-lwang-end----------------------------------------*
         END IF

*     HDF5 format output initialization (Maxwell)
#ifdef H5OUTPUT
         IF(KZ(46).GT.0.and.TTOT.EQ.0.D0) THEN
            KZ46_TMP = KZ(46)
*     At the first time, full dataset need to be stored,
*     thus templately change KZ(46) to full output option
            IF(KZ(46).EQ.1) KZ(46) = 3
            IF(KZ(46).EQ.2) KZ(46) = 4
            IF(rank.eq.0) CALL custom_output(NXTLST,NXTLEN,1)
            KZ(46) = KZ46_TMP
*
            if(rank.eq.0) then
               IF(KZ(46).eq.1)
     &              write(6,*) 'HDF5 output with time interval ',
     &              DTK(KZ(47)+1)
               IF(KZ(46).eq.3)
     &              write(6,*) 'HDF5 output every irregular ',
     &              'block step'
               IF(KZ(46).eq.2)
     &              write(6,*) 'CSV output with time interval ',
     &              DTK(KZ(47)+1)
               IF(KZ(46).eq.4)
     &              write(6,*) 'CSV output every irregular ',
     &              'block step'
            end if
         END IF
#endif

#ifdef GPU
         NN = NTOT + 10
         CALL GPUNB_OPEN(NN,rank)
*         CALL GPUPOT_INIT_FLOAT(rank)
*         CALL GPUPOT_INIT(rank,NN)
#endif
#ifdef SIMD
         NN = NTOT + NBIN0 - NPAIRS + 100
         IF (QVIR.LT.0.01.AND.TIME.EQ.0.0D0) NN = NN + 500
         CALL IRR_SIMD_OPEN(NN,LMAX,rank)
#endif
#ifdef PARALLEL         
*       Get the offset number for NXTLST_MPI
         IF (MOD(NMAX,ISIZE).EQ.0) THEN
            NPIECE = NMAX / ISIZE
         ELSE
            NPIECE = NMAX / ISIZE + 1
         END IF
#endif         
         IPSTART = .false.
      END IF
      call cputim(tttiend)
      ttintb = ttintb +(tttiend-tttbegin)*60
*
*       Search for high velocities after escape or KS/chain termination.
  999 IF (IPHASE.EQ.-1.OR.IPHASE.GE.2) THEN
          CALL HIVEL(0)
      END IF
*
#ifdef SIMD
      call cputim(ttsimda)
*     When necessary, update all particle data in AVX/SSE library
      IF(.NOT.FNTLOOP.OR.IPHASE.LT.0.OR.IPHASE.GE.3) THEN
!$omp parallel do private(I)
         DO I = IFIRST,NTOT
            CALL IRR_SIMD_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &           BODY(I),T0(I))
            CALL IRR_SIMD_SET_LIST(I,LIST(1,I))
         END DO
!$omp end parallel do
      END IF
      call cputim(ttsimdb)
      ttsimdsend = ttsimdsend +(ttsimdb-ttsimda)*60
#endif
      call cputim(tttbegin)

*     Reset control & regularization indicators.
      IPHASE = 0
      IKS = 0
*
      IF (IQ.LT.0) ICALL = 0
      IQ = 0
      call cputim(tttiend)
      ttintb = ttintb +(tttiend-tttbegin)*60
*
*       Find all particles due at next block time.
    1 CONTINUE
*
*
*       Redetermine TMIN after main change to catch new small steps after chain.
      ICALL = ICALL + 1
      IF (ICALL.EQ.2) GO TO 999
*
      call cputim(tttbegin)

*     See whether need to add delay list into NXTLST
      IF(NLSTDELAY(1).GT.0) call delay_add_tlist(T0,STEP,DTK)

*     Determine NXTLEN, TMIN and NXTLEVEL
      call next_tlist(TMIN,DTK,T0)

*     --07/07/14 23:44-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      IF(TMIN.GE.0.71682739257812500) THEN
c$$$         print*,'NXTLEN',NXTLEN,'NXTL',NXTLEVEL,'DTK',DTK(NXTLEVEL),
c$$$     &        'TMIN',TMIN,'T',TIME
c$$$      END IF
c$$$      DO L = 1, NXTLEN
c$$$         J=NXTLST(L)
c$$$         IF(T0(J)+STEP(J).NE.TMIN) then
c$$$            print*,'J',J,'N',NAME(J),'TMIN',TMIN,'T0+STEP',T0(J)+STEP(J)
c$$$     &           ,'T0',T0(J),'STEP',STEP(J),'DTKM',NDTMAX,DTK(NDTMAX),
c$$$     &           'NXTL',NXTLEVEL,'NXTLEN',NXTLEN,'L',L,'NXTLIMIT',
c$$$     &           NXTLIMIT
c$$$            print*,'NDTK',(K,NDTK(K),DTK(K),K=1,NDTMAX)
c$$$            call flush(6)
c$$$            call abort()
c$$$         END IF
c$$$C         IF(J.EQ.9951) print*,'L',L,'J',J,'T',TIME
c$$$      END DO
c$$$      DO L = NXTLEN+1,NXTLIMIT
c$$$         J = NXTLST(L)
c$$$         IF(T0(J)+STEP(J).LE.TMIN) THEN
c$$$            print*,'J',J,'N',NAME(J),'TMIN',TMIN,'T0+STEP',T0(J)+STEP(J)
c$$$     &           ,'T0',T0(J),'STEP',STEP(J),'DTKM',NDTMAX,DTK(NDTMAX),
c$$$     &           'NXTL',NXTLEVEL,'NXTLEN',NXTLEN,'L',L,'LOST!'
c$$$            print*,'NDTK',(K,NDTK(K),DTK(K),K=1,NDTMAX)
c$$$            call flush(6)
c$$$            call abort()
c$$$         END IF
c$$$C         IF(J.EQ.9951) print*,'L',L,'J',J,'T',TIME
c$$$      END DO
*     --07/07/14 23:44-lwang-end----------------------------------------*

*     Determine MPI NXTLST and Regular block
      NREG = 0
      DO L = 1, NXTLEN
         J = NXTLST(L)
#ifdef PARALLEL                  
         NXTLST_MPI(MOD(L-1,isize)*NPIECE + (L-1)/ISIZE +1) = J
#endif   
         TIMENW(J) = T0(J) + STEP(J)
         IF (TIMENW(J).GE.T0R(J) + STEPR(J)) THEN
            NREG = NREG + 1
            IREG(NREG) = J
         END IF
         IF (J.GT.N) THEN
            IPAIR = J - N
            IF (LIST(1,2*IPAIR-1).GT.0) NSTEPB = NSTEPB + 1
         END IF
      END DO
*
      LSKIP = NXTLEN 
      IF (LSKIP.LE.50) THEN
*       Update short timestep list for regularization search.
          CALL SHORT(NXTLEN,NXTLST,LSHORT)
      END IF
*
*     --11/08/13 23:15-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(rank.eq.0.and.time.ge.TCCT) then
c$$$         write(120,*) time,ifirst,n,ntot
c$$$         do L=1,NTOT
c$$$            write(120,*) name(L),tev(L)-tev0(L),step(L),stepr(L)
c$$$         END DO
c$$$         TCCT = TCCT + 0.1
c$$$      end if
c$$$      if(time+toff.ge.59.9) then
c$$$         if(rank.eq.0) call mydump(1,1)
c$$$         call mpi_barrier(MPI_COMM_WORLD,ierr)
c$$$         stop
c$$$      end if
*     --11/08/13 23:15-lwang-end----------------------------------------*
*       Set new time and save block time (for regularization terminations).
      TIME = TMIN
*     --11/09/13 19:12-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$*      IF(DMOD(TIME,step(nxtlst(1))).NE.0) then
c$$$         DTT = TIME - TBLOCK
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)
c$$$      if(time.ge.69.7665.and.rank.eq.0) then
c$$$      if(time.ge.0.71387481689453125) then
c$$$         write(6,*) 'TIME',time,'LEN',nxtlen,nxtlst(1),
c$$$     &     'NREG',NREG,
c$$$     &     'n',name(nxtlst(1)),'iks',iks,
c$$$     &     'STEP I R',step(nxtlst(1)),stepr(nxtlst(1)),
c$$$     *     'T0',t0(nxtlst(1)),'DM',DMOD(time,step(nxtlst(1))),
c$$$     &     'tb',tblock
C$$$         call flush(102+rank)
c$$$         call flush(6)
c$$$         call mpi_barrier(MPI_COMM_WORLD,ierr)  
c$$$*         stop
c$$$      end if
*     --11/09/13 19:12-lwang-end----------------------------------------*
      TBLOCK = TIME
      TTOT = TIME + TOFF
*     --11/03/13 10:22-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*      if(time.ge.3.25285E+01) then
c$$$         write(100+rank,*),'t',time,'nxtlen',nxtlen,'nxtlst',nxtlst(1),
c$$$     &     'n',name(nxtlst(1)),'iks',iks,'tpr',tprev,
c$$$     &     'STEP I R',step(nxtlst(1)),stepr(nxtlst(1)),
c$$$     *     'T0',t0(nxtlst(1)),'nstepi',nstepi
c$$$         call flush(100+rank)
*      end if
c$$$      if(time.ge.69.7668) then
c$$$         call adjust
c$$$         call flush(6)
c$$$         stop
c$$$      end if
*     --11/03/13 10:22-lwang-end----------------------------------------*
*
*     Check option for advancing interstellar clouds.
      IF (KZ(13).GT.0) THEN
          CALL CLINT
      END IF
*
*       Check optional integration of cluster guiding centre.
      IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
          IF (KZ(14).EQ.3.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
              CALL GCINT
          END IF
*       Include mass loss by gas expulsion (Kroupa et al. MN 321, 699).
          IF (MPDOT.GT.0.0D0.AND.TIME + TOFF.GT.TDELAY) THEN
              MP = MP0/(1.0 + MPDOT*(TIME + TOFF - TDELAY))
          END IF
      END IF

      call cputim(tttiend)
      ttintb = ttintb +(tttiend-tttbegin)*60

*
*       Include commensurability test (may be suppressed if no problems).
*     IF (STEP(IMIN).LT.1.0E-15.OR.DMOD(TIME,STEP(IMIN)).NE.0.0D0) THEN
*         WRITE (6,1005) IMIN, NAME(IMIN), NSTEPI, TIME, STEP(IMIN), 
*    &                   TIME/STEP(IMIN)
*1005     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
*    &                        2I5,I11,F12.5,1P,E9.1,0P,F16.4)
*         CALL ABORT
*     END IF
*
*       Check for new regularization at end of block.
      IF (IKS.GT.0) THEN
          TIME = TPREV
          IPHASE = 1
          GO TO 100
      END IF
*
*       Check next adjust time before beginning a new block.
      IF (TIME.GT.TADJ) THEN
          TIME = TADJ
          IPHASE = 3
          GO TO 100
      END IF
*
*       Also check output time in case DTADJ & DELTAT not commensurate.
      IF (TIME.GT.TNEXT) THEN
          TIME = TNEXT
          CALL OUTPUT
          GO TO 1
      END IF
*
*       See whether to advance any close encounters at first new time.
      IF (TIME.GT.TPREV) THEN
         call cputim(tt5)
*     --08/27/13 16:25-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            if (time.ge.3.90625000000000000E-003) then
c$$$               do iii=1,ntot
*         if(name(iii).eq.5287) then
*         write(113+rank,*)'5287a',iii,
*     *           t0(iii),step(iii),ibarcount,time,tprev
*           call flush(6)
*         end if
c$$$                  if(name(iii).eq.23871) then
*      print*,'suba',iii,x0(1,iii),x0dot(1,iii),time,tprev
c$$$                    print*,rank,'suba',iii,t0(iii),step(iii),time,tprev,
c$$$     *                    x(1,iii),x0(1,iii)
c$$$                    call flush(6)
c$$$                  end if
c$$$               end do
c$$$            end if
*          if(time.ge.6E-002) print*,rank,time
*          if(time.ge.6.6E-02) stop
*     --09/04/13 15:10-lwang-end----------------------------------------*
c$$$      call flush(6)
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)            
c$$$         if(time.ge.0.73) then
c$$$            write(100+rank,*),'bsub time',time,'iq',iq,'tpr',
c$$$     &           tprev,'tblock',tblock
c$$$            call flush(100+rank)
c$$$         end if
c$$$         if(rank.eq.0.and.time.ge.0.3756) then
c$$$            print*,'asub iq',iq,'time',time,'tprev',tprev
c$$$         end if
*     --08/27/13 16:25-lwang-end-debug----------------------------------*
         CALL SUBINT(IQ,I10)
*     --08/27/13 16:31-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         if(time.eq.1.70898437500000000E-003) then
c$$$            print*,'j',12079,'n',name(12079),'tpred',tpred(12079),
c$$$     ^           't',time
c$$$            call flush(6)
c$$$         end if
c$$$         if(time.ge.0.73) then
c$$$            write(100+rank,*),'asub time',time,'iq',iq,'tpr'
c$$$     &           ,tprev,'tblock',tblock
c$$$            call flush(100+rank)
c$$$         end if
c$$$      if(time.ge.TCCT)
c$$$         if(rank.eq.0.and.iq.ne.0) then
c$$$            print*,'bsub iq',iq,'time',time,'tprev',tprev
c$$$         end if
c$$$         if(iq.ne.0) print*,rank,'subint ',iq,time
c$$$         call flush(6)
c$$$            if (time.ge.3.90625000000000000E-003) then
c$$$               do iii=1,ntot
c$$$                  if(name(iii).eq.23871) then
c$$$c$$$*         write(113+rank,*)'5287b',iii,
c$$$c$$$*     *           t0(iii),step(iii),ibarcount,time,tprev
c$$$c$$$*         end if
c$$$c$$$         if(time.ge.1.41601562500000000E-002.and.name(iii).eq.5287) then
c$$$                    print*,rank,'subb',iii,t0(iii),step(iii),time,tprev,
c$$$     *                    x(1,iii),x0(1,iii)
c$$$                    call flush(6)
c$$$c$$$           stop
c$$$                 end if
c$$$              end do
c$$$              end if
c$$$         if(IPHASE.EQ.8) print*,'IPHASE',IPHASE,'IQ',IQ,'T',TIME
*     --08/27/13 16:31-lwang-end-debug----------------------------------*
         call cputim(tt6)
         ttksintgrt = ttksintgrt + (tt6-tt5)*60.
         ttks = ttks + (tt6-tt5)*60.
*
         IF (IQ.LT.0) GO TO 999
      END IF

*
*       Check regularization criterion for single particles.
      call cputim(tttbegin)
      IKS = 0
      ISMIN = 0
      DSMIN = DTMIN
      IF (LSKIP.LE.50) THEN
*       Search only in prepared list of short-step particles. (R.Sp.)
         ISHORT = LSHORT(1)
         DO L = 2,ISHORT+1
            I = LSHORT(L)
*       Search for minimum timestep candidate for not ordered steps (R.Sp.)
*       Beware that members of LSHORT may be members of KS pair (R.Sp.)
c$$$          If(step(i).LT.1.0D-8) THEN
c$$$             WRITE (6,24)  TIME, I, NAME(I), NXTLEN, NSTEPI,
c$$$     &            STEP(I), STEPR(I), NTOT, NPAIRS
c$$$ 24          FORMAT (' SMALL STEP!!  T I NAME LEN #I SI SR NT NP',F10.5,
c$$$     &            3I6,I11,1P,2E10.2,2I11)
c$$$             CALL FLUSH(6)
c$$$             STOP
c$$$  END IF
            IF (STEP(I).LT.DTMIN.AND.STEP(I).LT.DSMIN.AND.I.LE.N.AND.
     &           I.GE.IFIRST) THEN
               DSMIN = STEP(I)
               ISMIN = I
            END IF
         END DO
      END IF
*
*       See whether dominant body can be regularized.
      IF(ISMIN.GT.0) THEN
         CALL SEARCH(ISMIN,IKS)
*       Include close encounter search for low-eccentric massive binaries.
         IF (IKS.EQ.0.AND.NXTLEVEL.LE.7.AND.STEP(ISMIN).LT.4.0*DTMIN) 
     &        THEN
*       Consider massive single bodies in absence of subsystems.
            IF (ISMIN.LE.N.AND.BODY(I).GT.2.0*BODYM.AND.NSUB.EQ.0) THEN
*
*       Obtain two-body elements and relative perturbation.
               JMIN = 0
               CALL ORBIT(ISMIN,JMIN,SEMI,ECC,GI)
*     
               EB = -0.5*BODY(ISMIN)*BODY(JMIN)/SEMI
               IF (EB.LT.EBH.AND.GI.LT.0.25.AND.JMIN.GE.IFIRST) THEN
                  APO = SEMI*(1.0 + ECC)
*     Check eccentricity (cf. max perturbation) and neighbour radius.
                  IF (ECC.LT.0.25.AND.APO.LT.0.02*RS(ISMIN)) THEN
                     if(rank.eq.0) PRINT*, ' KS TRY: NAM E A EB ',
     *                    NAME(ISMIN), NAME(JMIN), ECC, SEMI, EB
                     CALL FLUSH(6)
                     IKS = IKS + 1
                     ICOMP = ISMIN
                     JCOMP = JMIN
                  END IF
               END IF
            END IF
         END IF
      END IF

      call cputim(tt1)
      ttintb=ttintb +(tt1-tttbegin)*60.
*
*       Initialize counters for irregular & regular integrations.
      TPREV = TIME
      NBLOCK = NBLOCK + 1

*     Initialize counters for irregular & regular integrations.
*      NREG = 0
*
*     write(6,*)' irr ',time,nxtlen,(nxtlst(k),step(k),k=1,5)
*     call flush(6)
*       Advance the irregular step for all particles in the current block.
*       Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
      IF(KZ(33).GT.0)JHIST(NXTLEN) = JHIST(NXTLEN) + 1
*
*     --08/27/13 16:38-lwang-debug--------------------------------------*
***** Note: Before nbint-----------------------------------------------**
c$$$      do L = 1 ,NXTLEN
c$$$         J = NXTLST(L)
c$$$         IF(rank.eq.0.and.NAME(J).EQ.12833) print*,'J',J,'N',NAME(J),'t'
c$$$     &        ,time,'T0',T0(J),'STEP',STEP(J),'TN',TIMENW(J)
c$$$      end do
c$$$      J = 12195
c$$$      JJ = 2*(J-N)-1
c$$$      print*,rank,'J',J,'N',NAME(J),'X',X(1,J),'X1',X(1,JJ),
c$$$     &     'X2',X(1,JJ+1),'T',TIME,'TPR',TPRED(J),'NNB',LIST(1,JJ)
c$$$      call flush(6)
*****      GET NEIGHBOR LIST      
c$$$      if(rank.eq.0.and.nxtlen.ge.5000) then
c$$$         do j=1,nxtlen
c$$$            jj = nxtlst(j)
c$$$            write(102+NGPUC,*) jj, list(1,jj)
c$$$         end do
c$$$*         tcct = tcct + 0.1
c$$$         NGPUC = NGPUC + 1
c$$$      end if
c$$$         if(time.eq.0.73070335388183594) then
C$$$      if(time.ge.TCCT) then
C$$$c$$$         do i=1,ifirst-1
C$$$c$$$            if(name(i).eq.16107) then
C$$$c$$$            write(100+rank,*) i,name(i),step(i),t0(i),time,body(i),
C$$$c$$$     *           gamma((i+1)/2),list(1,i),name((i+1)/2+N)
C$$$c$$$*               write(100+rank,*) i,name(i),x0(1,i),x0dot(1,i),
C$$$c$$$*     *              f(1,i),fdot(1,i),time
C$$$c$$$            end if
C$$$c$$$         end do
C$$$c$$$        write(102+rank,*) 'NXTLEN :',nxtlen,'t',time,'ibar',
C$$$c$$$     &     ibarcount,'nreg',nreg,'ireg',ireg(1),'ip',iphase
c$$$        do L=1,NXTLEN
c$$$           J = NXTLST(L)
c$$$           if(name(j).eq.1) then
c$$$c$$$              if(rank.eq.0.and.name(j).eq.2570) then
c$$$c$$$                 rmmjj=0.0
c$$$c$$$                 DO Jk=2,list(1,j)+1
c$$$c$$$                    jj=list(jk,j)
c$$$c$$$                    rmjjj=body(jj)/((X(1,jj)-X(1,J))**2
c$$$c$$$     &                   +(x(2,jj)-x(2,j))**2+
c$$$c$$$     ^                   (x(3,jj)-x(3,j))**2)
c$$$c$$$                    if(rmjjj.gt.rmmjj) then
c$$$c$$$                       immj=jj
c$$$c$$$                       rmmjj=rmjjj
c$$$c$$$                    end if
c$$$c$$$                 end do
c$$$c$$$                 print*,'J',J,'N',NAME(J),'JP',immj,'N',
c$$$c$$$     &                NAME(immj),' f',rmmjj,
c$$$c$$$     ^                'fv',rmmjj*body(j)/body(immj),'FI',
c$$$c$$$     &                sqrt(FI(1,J)**2+FI(2,J)**2+FI(3,J)**2),
c$$$c$$$     &                'T',TIME,'SI',STEP(J)
c$$$c$$$                 call flush(6)
c$$$c$$$              end if
c$$$            write(102+rank,109),l,j,name(j),x0(1:3,j),x0dot(1:3,j),
c$$$     *           t0(j),t0r(j),step(j),stepr(j),
c$$$     *           f(1:3,j),fdot(1:3,j),fi(1:3,j),fidot(1:3,j),
c$$$     &             d0(1:3,j),d1(1:3,j),
c$$$     *           d2(1:3,j),d3(1:3,j),d0r(1:3,j),d1r(1:3,j),
c$$$     *           d2r(1:3,j),d3r(1:3,j),body(j),time,list(1,j),
c$$$     &           NREG,NPAIRS,NXTLEN
c$$$C$$$           if(rank.eq.0.and.d2r(1,j).ne.d2r(1,j)) STOP
c$$$            write(102+rank,*)'LIST',LIST(1:LIST(1,J)+1,J)
c$$$c$$$            CALL IRR_SIMD_SET_JP(J,X0(1,J),X0DOT(1,J),
c$$$c$$$     &           F(1,J),FDOT(1,J),BODY(J),T0(J))
c$$$c$$$            call irr_simd_set_list(j,list(1,j))
c$$$           call flush(102+rank)
c$$$*           stop
c$$$          end if
c$$$        end do
c$$$*         end if
c$$$ 109    format(/,'L',I7,'J',I7,'  N',I7,'  X0',3E25.17,'  X0D',3E25.17,
c$$$     &       '  T0',F21.17,'  T0R',F21.17,'  STEP',F20.17,
c$$$     &       '  STEPR',F20.17,'  F',3E25.17,' FD',3E25.17,
c$$$     &       '  FI',3E25.17,'  FIDOT',3E25.17,'  D0',3E25.17,
c$$$     &       '  D1',3E25.17,'  D2',3E25.17,'  D3',3E25.17,
c$$$     &       '  D0R',3E25.17,'  D1R',3E25.17,'  D2R',3E25.17,
c$$$     &       '  D3R',3E25.17,'  M',E25.17,'  T',F21.17,'  NB',I4,
c$$$     &       '  NREG',I8,'  NPAIRS',I8,'  NXTLEN',I8)
c$$$        call mpi_barrier(MPI_COMM_WORLD,ierr)
C$$$C$$$        if (TCCT.gt.2.0) stop
C$$$        TCCT = TCCT + 0.05
C$$$      end if
c$$$      if(time.ge.2.0) stop
c$$$      if(time.ge.32.91898596286773682) stop
c$$$      if (time.ge.TCCT) then
c$$$        call mpi_barrier(MPI_COMM_WORLD,ierr)        
c$$$      end if
c$$$         if (time.ge.20.852127075195313) then
c$$$            print*,rank,time,'nxtlst ',(nxtlst(k),k=1,nxtlen)
c$$$         end if
c$$$      if(time.ge.54.192613631486893) stop
*     --08/27/13 16:38-lwang-end-debug----------------------------------*

*     Do full prediction when regular block step come
      ipredall = .false.
      IF(NREG.GT.0) THEN
         call cputim(tttpre1)
*       Predict all particles (except TPRED=TIME) in C++ on host. 
*         CALL CXVPRED(IFIRST,NTOT,TIME,T0,X0,X0DOT,F,FDOT,X,XDOT,TPRED)
         call xbpredall
         call cputim(tttpre2)
         ttpre = ttpre + (tttpre2-tttpre1)*60
         ttpreall = ttpreall + (tttpre2-tttpre1)*60
      END IF

*     Predict chain variables and perturbers at new block-time.
      IF (NCH.GT.0) THEN
          CALL JPRED_int(ICH,TIME)
          CALL XCPRED(2)
      END IF

#ifdef PARALLEL
      IF(NXTLEN.LE.isernb)THEN
#endif
*
         call cputim(ttnb1)
*#ifndef SIMD
*     Predict x and xdot
*         call xbpred(NREG,1,NXTLEN,NXTLST)
*#endif
*         call cputim(ttnb2)
*         ttpre = ttpre + (ttnb2-ttnb1)*60.
#ifdef SIMD
*     --04/18/14 23:07-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         DO L = 1,NXTLEN
c$$$            I = NXTLST(L)
c$$$            call IRR_SIMD_SET_JP(I,X0(1,I),x0dot(1,I),f(1,i),fdot(1,i),
c$$$     &           body(i),t0(i))
c$$$         END DO
*     --04/18/14 23:07-lwang-end----------------------------------------*
*     Use AVX/SSE to calculate irregular force and fdot
         call IRR_SIMD_FIRR_VEC(TIME,NXTLEN,NXTLST,FIRR_SIMD,FD_SIMD)
         call cputim(ttsimdc)
         ttsimdcalc = ttsimdcalc + (ttsimdc-ttnb1)*60.
*     --04/18/14 22:27-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         DO L = 1,NXTLEN
c$$$            I = NXTLST(L)
c$$$            IF(TIME.EQ.30.04376220703125000.and.I.eq.7660) THEN
c$$$               NNNNB = LIST(1,I)+1
c$$$               DO LL = 2, NNNNB
c$$$                  K = LIST(LL,I)
c$$$                  write(102+rank,*),'K',K,'N',NAME(K),'X0',X0(1:3,K),
c$$$     &           'XD',X0DOT(1:3,K),'F',F(1:3,K),'FD',FDOT(1:3,K),'T',
c$$$     ^                 TIME,'I',I,'N',NAME(I)
c$$$               END DO
c$$$            write(102+rank,*) 'L',L,'I',I,'N',NAME(I),'X',X0(1:3,I),
c$$$     &           'XD',X0DOT(1:3,I),'F',F(1:3,I),'FD',FDOT(1:3,I),
c$$$     &           'FIRR',FIRR_SIMD(1:3,L),'FD',FD_SIMD(1:3,L)
c$$$     &           ,'T',TIME,'NXT',NXTLEN,'NNB',NNNNB,'IFIRST',IFIRST
c$$$            call flush(102+rank)
c$$$         end if
c$$$         end do
*     --04/18/14 22:27-lwang-end----------------------------------------*
!$omp parallel do if(NXTLEN.GE.ITHREAD) private(L,I)
         DO L = 1,NXTLEN
            I = NXTLST(L)
            CALL NBINT_COR(I,FIRR_SIMD(1,L),FD_SIMD(1,L))
         END DO
!$omp end parallel do

#else         
*     Neighbor Integration
!$omp parallel do if(NXTLEN.GE.ITHREAD) private(L,I)
         DO L = 1,NXTLEN
            I = NXTLST(L)
*     --04/04/14 15:24-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            IF(name(i).eq.8098) then
c$$$               nnb=list(1,i)+1
c$$$               write(102+rank,*) 'i',i,'n',name(i),'x',x(1,i),
c$$$     &              'xd',xdot(1,i),'nnb',nnb,'t',time
c$$$               do lj=2,nnb
c$$$                  j=list(lj,i)
c$$$                  write(102+rank,*) 'j',j,'n',name(j),
c$$$     ^                 'x',x(1,j),'xd',xdot(1,j),'t',time
c$$$                  if(j.ge.n) then
c$$$                     j1=2*(j-n)-1
c$$$                     if(list(1,j1).gt.0) then
c$$$                        write(102+rank,*) 'j1',j1,'n',name(j1),
c$$$     ^                       'x',x(1,j1),x(1,j1+1),'xd',xdot(1,j1),
c$$$     &                       xdot(1,j1+1),'t',time,
c$$$     &                       'u',u(1,j-n),'ud',udot(1,j-n),
c$$$     &                       'fu',fu(1,j-n),'fud',fudot(1,j-n),
c$$$     &                       'fudot2',fudot(1,j-n),
c$$$     &                       'fudot2',fudot2(1,j-n),'u0',u0(1,j-n)
c$$$                     end if
c$$$                  end if
c$$$               end do
c$$$               call flush(102+rank)
c$$$            end if
*     --04/04/14 15:24-lwang-end----------------------------------------*
            CALL NBINT(I)
         END DO
!$omp end parallel do
#endif
*     --04/02/14 15:33-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$        do L=1,nxtlen
c$$$           J = nxtlst(L)
c$$$            write(102+rank,*),'NBj',j,'n',name(j),
c$$$     &          'xn',xn(1:3,j),'xnd',xndot(1:3,j),
c$$$     *           's',step(j),'fi',fi(1:3,j),'fid',fidot(1:3,j),
c$$$     *           'd2',d2(1:3,j),'d3',d3(1:3,j),'t',time
c$$$           call flush(102+rank)
c$$$        end do
*     --04/02/14 15:33-lwang-end----------------------------------------*

!$omp parallel do if(NXTLEN.GE.ITHREAD) private(L,I,DTR)
         DO L = 1,NXTLEN
            I = NXTLST(L)
            X0(1:3,I) = XN(1:3,I)
c$$$            X(1:3,I) = XN(1:3,I)
            X0DOT(1:3,I) = XNDOT(1:3,I)
c$$$            XDOT(1:3,I) = XNDOT(1:3,I)
            D0(1:3,I) = FI(1:3,I)
            D1(1:3,I) = FIDOT(1:3,I)
*     Save new block step and update T0 & next time
            T0(I) = TIME
            TIMENW(I) = T0(I) + STEP(I)
*     Set non-zero indicator for new regular force.
            IF (T0R(I) + STEPR(I).GT.TIME) THEN
*     Extrapolate regular force & first derivatives to obtain F & FDOT.
              DTR = TIME - T0R(I)
              F(1,I) = 0.5*(FRDOT(1,I)*DTR + FR(1,I) + FI(1,I))
              F(2,I) = 0.5*(FRDOT(2,I)*DTR + FR(2,I) + FI(2,I))
              F(3,I) = 0.5*(FRDOT(3,I)*DTR + FR(3,I) + FI(3,I))
              FDOT(1,I) = ONE6*(FRDOT(1,I) + FIDOT(1,I))
              FDOT(2,I) = ONE6*(FRDOT(2,I) + FIDOT(2,I))
              FDOT(3,I) = ONE6*(FRDOT(3,I) + FIDOT(3,I))
           END IF
         END DO
!$omp end parallel do

         call cputim(tt3)
         ttirr = ttirr + (tt3-ttnb1)*60.
*
#ifdef PARALLEL
*       start PARALLEL section for nbint
        ELSE
*
           call cputim(tt998)
           call mpi_barrier(MPI_COMM_WORLD,ierr)
           call cputim(tt999)
           ttbar = ttbar + (tt999-tt998)*60
           ibarcount=ibarcount+1
           ttbarnb = ttbarnb + (tt999-tt998)*60
*           print*, 'Nbint Barrier 1: rank iphase ttbarnb ttbar',
*     &          ' dt nxtlen bcount',rank,ixxx,ttbarnb,ttbar,tttdet,
*     &          nxtlen,ibarcount
*
           istart_lst = rank * NPIECE + 1
           ishift = MOD(NXTLEN,isize)
           idivide = NXTLEN / isize
           if (rank.LT.ishift) then
              iend_lst = idivide + istart_lst
           else
              iend_lst = idivide + istart_lst - 1
           end if
           
           nl = NXTLEN
*     
           inl = nl/isize
           jsize = isize*inl
           idiff = nl - jsize
           irun = 0
*     
           do ix = 1,isize
              inum(ix)=inl
              if(ix.le.idiff)inum(ix) = inum(ix) + 1
              ista(ix) = irun+1
              if(ista(ix).gt.nl)inum(ix) = 0
              irun = irun + inum(ix)
           end do
*     
           istart = ista(rank+1)
           iend = ista(rank+1) + inum(rank+1) - 1
*
*     if(time.lt.0.1d0)then
*     print*,' rank ',rank,' NXTLEN ',NXTLEN,TIME,
*    *    ' istart,iend=',istart,iend
*     end if
*
c$$$      call cputim(tt998)
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)
c$$$      call cputim(tt999)
c$$$      tttdet=(tt999-tt998)*60
c$$$      ttbar = ttbar + (tt999-tt998)*60
c$$$      ibarcount=ibarcount+1
c$$$      ttbarnb = ttbarnb + (tt999-tt998)*60
*      print*, 'Nbint Barrier 2: rank iphase ttbarnb ttbar',
*     &     ' dt nxtlen bcount',rank,ixxx,ttbarnb,ttbar,tttdet,
*     &     nxtlen,ibarcount
*
*     if(time.lt.0.1d0)then
*     print*,' rank ',rank,' NXT in DO L ',L,TIME,
*    *    ' istart,iend=',istart,iend
*     end if

*     --11/09/13 19:09-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$           print*,'TEST ',time,rank,istart_lst,iend_lst,
c$$$     &          nxtlst_mpi(istart_lst:iend_lst)
c$$$           call flush(6)
c$$$           call mpi_barrier(MPI_COMM_WORLD,ierr)
*     --11/09/13 19:09-lwang-end----------------------------------------*
           call cputim(ttnb1)
*#ifndef SIMD
*     Predict x and xdot
*           call xbpred(NREG,istart_lst,iend_lst,NXTLST_MPI)
*#endif
*           call cputim(ttnb2)
*           ttpre = ttpre + (ttnb2-ttnb1)*60.

#ifdef SIMD
*     Use AVX/SSE to calculate irregular force and fdot
           NPART_LEN = iend_lst - istart_lst + 1
           call IRR_SIMD_FIRR_VEC(TIME,NPART_LEN,NXTLST_MPI(istart_lst),
     &          FIRR_SIMD,FD_SIMD)
           call cputim(ttsimdc)
           ttsimdcalc = ttsimdcalc + (ttsimdc-ttnb1)*60.
           
           LL = istart
           I_STEP_LEN = ithread*icore
           DO L = istart_lst,iend_lst,I_STEP_LEN
              if(iend_lst-L+1.GE.I_STEP_LEN) THEN
                 I_LOOP_LEN = I_STEP_LEN - 1
              else
                 I_LOOP_LEN = iend_lst - L
              end if
!$omp parallel do private(LK,I,IC_K)
              DO LK = 0, I_LOOP_LEN
                 I = NXTLST_MPI(L + LK)
                 IC_K = L + LK - istart_lst + 1
                 CALL NBINT_COR(I,FIRR_SIMD(1,IC_K),FD_SIMD(1,IC_K))
                 XMPI(1:3,LL+LK) = XN(1:3,I)
                 XMPI(4:6,LL+LK) = XNDOT(1:3,I)
                 XMPI(7:9,LL+LK) = D2(1:3,I)
                 XMPI(10:12,LL+LK) = D3(1:3,I)
                 XMPI(13:15,LL+LK) = FI(1:3,I)
                 XMPI(16:18,LL+LK) = FIDOT(1:3,I)
                 XMPI(19,LL+LK) = STEP(I)
              END DO
!$omp end parallel do
              LL = LL + I_STEP_LEN
           END DO
*     --04/18/14 21:46-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         if(rank.eq.0) then
c$$$            print*,'NXTLEN',NPART_LEN,
c$$$     &           'FIRR_SIMD',FIRR_SIMD(1:3,1:NPART_LEN+1),
c$$$     &       'T', TIME
c$$$            call flush(6)
c$$$         end if
c$$$c$$$         call mpi_barrier(MPI_COMM_WORLD,ierr)
c$$$         stop
*     --04/18/14 21:46-lwang-end----------------------------------------*
#else 
*     Neighbor Integration
           LL = istart
           I_STEP_LEN = ithread*icore
           DO L = istart_lst,iend_lst,I_STEP_LEN
              if(iend_lst-L+1.GE.I_STEP_LEN) THEN
                 I_LOOP_LEN = I_STEP_LEN - 1
              else
                 I_LOOP_LEN = iend_lst - L
              end if
!$omp parallel do private(LK,I)
              DO LK = 0, I_LOOP_LEN
                 I = NXTLST_MPI(L + LK)
                 CALL NBINT(I)
                 XMPI(1:3,LL+LK) = XN(1:3,I)
                 XMPI(4:6,LL+LK) = XNDOT(1:3,I)
                 XMPI(7:9,LL+LK) = D2(1:3,I)
                 XMPI(10:12,LL+LK) = D3(1:3,I)
                 XMPI(13:15,LL+LK) = FI(1:3,I)
                 XMPI(16:18,LL+LK) = FIDOT(1:3,I)
                 XMPI(19,LL+LK) = STEP(I)
*     --04/09/14 20:56-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$                 write(102+rank,*),'NB I',I,'N',Name(i),'FI',FI(1,I),
c$$$     &                'FD',FIDOT(1,I),'XN',XN(1,I),'X',X(1,I),'T',time
c$$$                 IF(XN(1,I).NE.XN(1,I)) STOP
c$$$                 call flush(102)
*     --04/09/14 20:56-lwang-end----------------------------------------*
              END DO
!$omp end parallel do
              LL = LL + I_STEP_LEN
           END DO
#endif           
*
           call cputim(tt3)
           ttirr = ttirr + (tt3-ttnb1)*60.
*
*        Distribute variables into private vectors again T3E (R.Sp.)
*
           isend = rank + 1
           if(isend.eq.isize)isend = 0
           irecv = rank - 1
           if(irecv.eq.-1)irecv = isize - 1
*
           do ir = 0,isize-2
*     
              irank = rank - ir
              if(irank.lt.0)irank=irank+isize
*     
              istsen=ista(irank+1)
              icnt = inum(irank+1)
*     
              if(irank.eq.0)irank=isize
              istrec = ista(irank)
              icnt2 = inum(irank)
*
*     if(time.lt.0.1d0.and.icnt.gt.0)then
*     print*,' NXT: rank t',rank,time,' ir ',ir,' send ',istsen,
*    *    ' thru ',istsen+icnt-1,' to ',isend,' cnt ',icnt,
*    *    ' istart,iend=',istart,iend
*     end if
*     if(time.lt.0.1d0.and.icnt2.gt.0)then
*     print*,' NXT: rank t',rank,time,' ir ',ir,' recv ',istrec,
*    *    ' thru ',istrec+icnt2-1,' fr ',irecv,' cnt2 ',icnt2,
*    *    ' istart,iend=',istart,iend
*     end if
*
#ifdef PUREMPI
              call cputim(tta)
              CALL MPI_SENDRECV(XMPI(1,istsen),19*icnt,MPI_REAL8,isend,
     *             rank,XMPI(1,istrec),19*icnt2,MPI_REAL8,irecv,irecv,
     *             MPI_COMM_WORLD,status,ierr)
              call cputim(ttb)
              call mpi_barrier(MPI_COMM_WORLD,ierr)
              call cputim(tt999)
              ttbar = ttbar + (tt999-ttb)*60
              ttbarnb = ttbarnb + (tt999-ttb)*60
              ibarcount=ibarcount+1
*     print*, 'Nbint Barrier 3: rank iphase ttbarnb ttbar',
*     &     ' dt nxtlen bcount',rank,ixxx,ttbarnb,ttbar,tttdet,
*     &     nxtlen,ibarcount
              xtsub1 = xtsub1 + dble(19*8*(icnt+icnt2))
              ttsub = ttsub + (ttb-tta)*60.
#endif
*
#ifdef SHMEM
              call barrier()
              call shmem_get(XMPI(1,istrec),XMPI(1,istrec),19*icnt2,
     &             irecv)
#endif
*
           end do
*
*          call mpi_barrier(MPI_COMM_WORLD,ierr)
*
           call cputim(ttmov3)
           Loffset = 0
           DO LK = 0, isize-1
              istart = LK*NPIECE + 1
              if (LK.LT.ishift) then
                 iend = istart + idivide
              else
                 iend = istart + idivide - 1
              end if
!$omp parallel do if(idivide.GE.ITHREAD) private(L,LL,I,DTR)
              DO LL = istart, iend
                 I = NXTLST_MPI(LL)
                 L = LL - Loffset
                 X0(1:3,I) = XMPI(1:3,L)
*              X(1:3,I) = XMPI(1:3,L)
                 X0DOT(1:3,I) = XMPI(4:6,L)
*              XDOT(1:3,I) = XMPI(4:6,L)
                 D2(1:3,I) = XMPI(7:9,L)
                 D3(1:3,I) = XMPI(10:12,L)
                 FI(1:3,I) = XMPI(13:15,L)
                 D0(1:3,I) = FI(1:3,I)
                 FIDOT(1:3,I) = XMPI(16:18,L)
                 D1(1:3,I) = FIDOT(1:3,I)
                 STEP(I) = XMPI(19,L)
*     Save new block step and update T0 & next time
                 T0(I) = TIME
                 TIMENW(I) = T0(I) + STEP(I)
*     Set non-zero indicator for new regular force.
                 IF (T0R(I) + STEPR(I).GT.TIME) THEN
*     Extrapolate regular force & first derivatives to obtain F & FDOT.
                    DTR = TIME - T0R(I)
                    F(1,I) = 0.5*(FRDOT(1,I)*DTR + FR(1,I) + FI(1,I))
                    F(2,I) = 0.5*(FRDOT(2,I)*DTR + FR(2,I) + FI(2,I))
                    F(3,I) = 0.5*(FRDOT(3,I)*DTR + FR(3,I) + FI(3,I))
                    FDOT(1,I) = ONE6*(FRDOT(1,I) + FIDOT(1,I))
                    FDOT(2,I) = ONE6*(FRDOT(2,I) + FIDOT(2,I))
                    FDOT(3,I) = ONE6*(FRDOT(3,I) + FIDOT(3,I))
                 END IF
              END DO
!$omp end parallel do
              Loffset = Loffset + (LK+1)*NPIECE - iend
           END DO
*
           call cputim(tt4)
           ttmov = ttmov + (tt4-ttmov3)*60.
*     
        END IF
*          End PARALLEL section for nbint
#endif
*

c$$$        call cputim(ttt32)
c$$$**!$omp parallel do if(NXTLEN.GE.ithread)
c$$$**!$omp& default(shared) private(L,I,DTR,K)
c$$$        DO L = 1,NXTLEN
c$$$           I = NXTLST(L)
c$$$*       Save new block step and update T0 & next time
c$$$           T0(I) = TIME
c$$$           TIMENW(I) = T0(I) + STEP(I)
c$$$           
c$$$*     Set non-zero indicator for new regular force.
c$$$           IF (T0R(I) + STEPR(I).GT.TIME) THEN
c$$$*     *!$omp critical 
c$$$*              NREG = NREG + 1
c$$$*              IREG(NREG) = I
c$$$**!$omp end critical
c$$$*          ELSE
c$$$*       Extrapolate regular force & first derivatives to obtain F & FDOT.
c$$$              DTR = TIME - T0R(I)
c$$$              F(1,I) = 0.5*(FRDOT(1,I)*DTR + FR(1,I) + FI(1,I))
c$$$              F(2,I) = 0.5*(FRDOT(2,I)*DTR + FR(2,I) + FI(2,I))
c$$$              F(3,I) = 0.5*(FRDOT(3,I)*DTR + FR(3,I) + FI(3,I))
c$$$              FDOT(1,I) = ONE6*(FRDOT(1,I) + FIDOT(1,I))
c$$$              FDOT(2,I) = ONE6*(FRDOT(2,I) + FIDOT(2,I))
c$$$              FDOT(3,I) = ONE6*(FRDOT(3,I) + FIDOT(3,I))
c$$$* Higher order extrapolation?
c$$$*                 F(K,I) = FI(K,I) + FR(K,I) + DTR*(FRDOT(K,I)
c$$$*    *                + DTR*(D2R(K,I)/2.D0 + DTR*D3R(K,I)/6.D0))
c$$$*                 FDOT(K,I) = FIDOT(K,I) + FRDOT(K,I)
c$$$*    *                + DTR*(D2R(K,I) + DTR*D3R(K,I)/2.D0)
c$$$*                 F(K,I) = F(K,I)/2.D0
c$$$*                 FDOT(K,I) = FDOT(K,I)/6.D0
c$$$           END IF
c$$$          
c$$$*     X0(K,I) = XN(K,I)
c$$$*     X0DOT(K,I) = XNDOT(K,I)
c$$$        END DO
c$$$**!$omp end parallel do
*
*     --09/26/13 12:46-lwang-debug--------------------------------------*
***** Note: After nbint------------------------------------------------**
c$$$         if(time.ge.0.73070144653320313) then
c$$$            j = 15967
c$$$            kk = 1
c$$$c$$$         do L=1,nxtlen
c$$$c$$$            J = nxtlst(L)
c$$$c$$$            if(name(j).eq.16107) then
c$$$            write(6,*),rank,'Anbint',j,'n',name(j),'x0',x0(kk,j),
c$$$     *        'x0dot',x0dot(kk,j),'t0',t0(j),
c$$$     &           'step',step(j),'stepr',stepr(j),
c$$$     *        'f',f(kk,j),'fdot',fdot(kk,j),
c$$$     *       'fi',fi(kk,j),'fidot',fidot(kk,j),
c$$$     *        'd0',d0(kk,j),'d1',d1(kk,j),
c$$$     *        'd2',d2(kk,j),'d3',d3(kk,j),
c$$$     *        'd0r',d0r(kk,j),'d1r',d1r(kk,j),
c$$$     *        'd2r',d2r(kk,j),'d3r',d3r(kk,j),
c$$$     *        'body',body(j),'time',time,'nnb',list(1,j)
c$$$*            call flush(100+rank)
c$$$            call flush(6)
c$$$c$$$            end if
c$$$c$$$         end do
c$$$         call mpi_barrier(MPI_COMM_WORLD,ierr)
c$$$c$$$*        if (TCCT.gt.0.111) stop
c$$$c$$$*        TCCT = TCCT + 0.001
c$$$      end if
c$$$      if(time.ge.0.73070335388183594) stop
*     --09/26/13 12:46-lwang-end----------------------------------------*

*
*       See whether any KS candidates are in the same block.
        IF (IKS.GT.0) THEN
*       Accept same time, otherwise reduce STEP(ICOMP) and/or delay.
           IF (T0(JCOMP).EQ.T0(ICOMP)) THEN
*             call cputim(ttnb1)
*             call xbpredall
*             call cputim(ttnb2)
*             ttpre = ttpre + (ttnb2-ttnb1)*60.
              I = ICOMP
              ICOMP = MIN(ICOMP,JCOMP)
              JCOMP = MAX(I,JCOMP)
              call jpred_int(icomp,time)
              call jpred_int(jcomp,time)
*     --09/25/13 21:44-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$              if(time.ge.0.3) then
c$$$                 print*,rank,'icomp',icomp,name(icomp),
c$$$     &            x(1,icomp),xdot(1,icomp),x0(1,icomp),x0dot(i,icomp)
c$$$              end if
*     --09/25/13 21:44-lwang-end----------------------------------------*
           ELSE IF (T0(JCOMP) + STEP(JCOMP).LT.T0(ICOMP)) THEN
              STEP(ICOMP) = 0.5D0*STEP(ICOMP)
              TIMENW(ICOMP) = T0(ICOMP) + STEP(ICOMP)
              IKS = 0
           ELSE
              IKS = 0
           END IF
        END IF
*
      NSTEPI = NSTEPI + NXTLEN

*
c$$$      call cputim(ttt33)
c$$$      ttnewt = ttnewt + (ttt33 - ttt32)*60.

***** --Regular force calculation--------------------------------------**
      IF(NREG.GT.0)THEN
      
#ifdef GPU
*       Send all single particles to GPU memory
         call cputim(tt51)
         NN = NTOT - IFIRST + 1
         CALL GPUNB_SEND(NN,BODY(IFIRST),X(1,IFIRST),XDOT(1,IFIRST))
         call cputim(tt52)
*     --09/26/13 16:58-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(rank.eq.0) print*,'NREG',NXTLEN,NREG,TIME,NGPUC,tt52-tt51
c$$$      NGPUC = NGPUC + 1
c$$$      if(rank.eq.0) then
c$$$         do i=1,nreg
c$$$            L=ireg(i)
c$$$            write(100,*) NREG, Name(L), STEPR(L),TIME
c$$$         end do
c$$$      end if
*     --09/26/13 16:58-lwang-end----------------------------------------*
         ttgrcomm = ttgrcomm + (tt52-tt51)*60.0
#endif

*
*       Block-Step Level Diagnostics (R.Sp. 29.Apr. 1993)
         IF(KZ(33).GT.1)JHISTR(NREG) = JHISTR(NREG) + 1
*

#ifdef PARALLEL
         IF(NREG.LE.iserreg)THEN
#endif
*
            call cputim(tt1)
#ifdef GPU
            CALL CALC_REG_ON_GPU(IREG,1,NREG)
            call cputim(tt2)
            ttgrcalc = ttgrcalc + (tt2-tt1)*60.0
#else
*
            DO L = 1,NREG
               I = IREG(L)
               CALL REGINT(I,LIST(1,I))
            END DO
#endif            
*
!$omp parallel do private(L,I)
            DO L = 1,NREG
               I = IREG(L)
               F(1,I) = 0.5D0*(FI(1,I) + FR(1,I))
               F(2,I) = 0.5D0*(FI(2,I) + FR(2,I))
               F(3,I) = 0.5D0*(FI(3,I) + FR(3,I))
               FDOT(1,I) = ONE6*(FIDOT(1,I) + FRDOT(1,I))
               FDOT(2,I) = ONE6*(FIDOT(2,I) + FRDOT(2,I))
               FDOT(3,I) = ONE6*(FIDOT(3,I) + FRDOT(3,I))               
               X0(1:3,I) = XN(1:3,I)
c$$$               X(1:3,I) = XN(1:3,I)
               X0DOT(1:3,I) = XNDOT(1:3,I)
c$$$               XDOT(1:3,I) = XNDOT(1:3,I)
*     Copy neighbor list
*               NNB = IMPI(1,L) + 1 
*               LIST(1:NNB,I) = IMPI(1:NNB,L)
               IF(LIST(1,I).GT.0) then
!$omp critical                   
                  RSMIN = MIN(RSMIN,RS(I))
!$omp end critical
               end if
#ifdef SIMD
*     Update neighbor list for AVX/SSE library
               CALL IRR_SIMD_SET_LIST(I,LIST(1,I))
#endif
*     --04/09/14 18:41-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$               J = I
c$$$               write(102+rank,119),j,name(j),x0(1:3,j),x0dot(1:3,j),
c$$$     *              t0(j),t0r(j),step(j),stepr(j),
c$$$     *              f(1:3,j),fdot(1:3,j),fi(1:3,j),fidot(1:3,j),
c$$$     &              d0(1:3,j),d1(1:3,j),
c$$$     *              d2(1:3,j),d3(1:3,j),d0r(1:3,j),d1r(1:3,j),
c$$$     *              d2r(1:3,j),d3r(1:3,j),body(j),time,list(1,j),
c$$$     &              NREG,NPAIRS
c$$$               call flush(102+rank)
c$$$c$$$  end if
c$$$ 119        format(/,'AJ',I7,'  N',I7,'  X0',3E25.17,'  X0D',3E25.17,
c$$$     &           '  T0',F21.17,'  T0R',F21.17,'  STEP',F20.17,
c$$$     &           '  STEPR',F20.17,'  F',3E25.17,' FD',3E25.17,
c$$$     &           '  FI',3E25.17,'  FIDOT',3E25.17,'  D0',3E25.17,
c$$$     &           '  D1',3E25.17,'  D2',3E25.17,'  D3',3E25.17,
c$$$     &           '  D0R',3E25.17,'  D1R',3E25.17,'  D2R',3E25.17,
c$$$     &           '  D3R',3E25.17,'  M',E25.17,'  T',F21.17,'  NB',I4,
c$$$     &           '  NREG',I8,'  NPAIRS',I8)
*     --04/09/14 18:41-lwang-end----------------------------------------*
            END DO
!$omp end parallel do
*     
            NBSUM = 0
            call cputim(tt2)
            ttreg = ttreg + (tt2-tt1)*60.

#ifdef PARALLEL
*       Start PARALLEL section for regint
         ELSE
*
            call cputim(tt998)
            call mpi_barrier(MPI_COMM_WORLD,ierr)
            call cputim(tt999)
            ibarcount=ibarcount+1
            ttbar = ttbar + (tt999-tt998)*60
            ttbarreg = ttbarreg +(tt999-tt998)*60

*
            nl = NREG
*     
            inl = nl/isize
            jsize = isize*inl
            idiff = nl - jsize
            irun = 0
*     
            do ix = 1,isize
               inum(ix)=inl
               if(ix.le.idiff)inum(ix) = inum(ix) + 1
               ista(ix) = irun+1
               if(ista(ix).gt.nl)inum(ix) = 0
               irun = irun + inum(ix)
            end do
*     
            istart = ista(rank+1)
            iend = ista(rank+1) + inum(rank+1) - 1
c$$$      if (ixxxx.le.10000000.and.time.gt.1.33.and.rank.eq.0) then
c$$$         PRINT*,' Start reg block rank ',rank,
c$$$     *        ' TIME=',TIME,' NREG=',NREG
c$$$         PRINT*,' Start reg block rank ',rank,
c$$$     *        ' istart,iend=',istart,iend
c$$$         call flush(6)
c$$$      end if
*
            NC5 = NBFULL
            NC6 = NBVOID
            NC11 = NBSMIN
            NC15 = NBDIS2
            NC19 = NBDIS
            NC20 = NLSMIN
            NC30 = NFAST
            NCX = NICONV
            NCY = NBFAST
            NBSUM = NBFLUX
            NC12 = NBPRED
*

            call cputim(tt1)
#ifdef GPU
            CALL CALC_REG_ON_GPU(IREG,istart,iend)
            call cputim(tt2)
            ttgrcalc2 = ttgrcalc2 + (tt2-tt1)*60.0
#else
*
            DO L = istart,iend
               I = IREG(L)
*
               NNB = LIST(1,I) + 1
               IMPI(1:NNB,L) = LIST(1:NNB,I)
*
               CALL REGINT(I,IMPI(1,L))
            END DO
            call cputim(tt2)
#endif
            ttreg = ttreg + (tt2-tt1)*60.0

!$omp parallel do private(L,I,NNB)
            DO L = istart,iend
               I = IREG(L)
#ifdef GPU
               NNB = LIST(1,I) + 1
               IMPI(1:NNB,L) = LIST(1:NNB,I)
#endif               
               YMPI(37,L) = STEPR(I)
               YMPI(38,L) = STEP(I)
               YMPI(39,L) = TIMENW(I)
               YMPI(40,L) = RS(I)
               YMPI(41,L) = T0R(I)
*
               YMPI(1:3,L) = XN(1:3,I)
               YMPI(4:6,L) = XNDOT(1:3,I)
               YMPI(7:9,L) = FI(1:3,I)
               YMPI(10:12,L) = FR(1:3,I)
               YMPI(13:15,L) = FIDOT(1:3,I)
               YMPI(16:18,L) = FRDOT(1:3,I)
               YMPI(19:21,L) = D1R(1:3,I)
               YMPI(22:24,L) = D2R(1:3,I)
               YMPI(25:27,L) = D3R(1:3,I)
*     Save corrections of force polynomials from routine fpcorr
               IF (KZ(38).GT.0.OR.I.GT.N) THEN
                  YMPI(28:30,L) = D1(1:3,I)
                  YMPI(31:33,L) = D2(1:3,I)
                  YMPI(34:36,L) = D3(1:3,I)
               END IF
            END DO
!$omp end parallel do      

            JMPI(1) = NBFULL - NC5
            JMPI(2) = NBVOID - NC6
            JMPI(3) = NBSMIN - NC11
            JMPI(4) = NBDIS2 - NC15
            JMPI(5) = NBDIS - NC19
            JMPI(6) = NLSMIN - NC20
            JMPI(7) = NFAST - NC30
            JMPI(8) = NICONV - NCX
            JMPI(9) = NBFAST - NCY
            JMPI(10) = NBFLUX - NBSUM
            JMPI(11) = NBPRED - NC12
            call cputim(tt3)
            ttmov = ttmov +(tt3-tt2)*60.
*
*
*        Distribute variables into private vectors again T3E (R.Sp.)
            isend = rank + 1
            if(isend.eq.isize)isend = 0
            irecv = rank - 1
            if(irecv.eq.-1)irecv = isize - 1
*
            do ir = 0,isize-2
*     
               irank = rank - ir
               if(irank.lt.0)irank=irank+isize
*     
               istsen=ista(irank+1)
               icnt = inum(irank+1)
*     
               if(irank.eq.0)irank=isize
               istrec = ista(irank)
               icnt2 = inum(irank)
*
c$$$      if (ixxxx.le.10000000.and.time.gt.1.33) then
c$$$         print*,' INTGRT-R: rank t',rank,time,' ir ',ir,' send ',istsen,
c$$$     *        ' thru ',istsen+icnt-1,' to ',isend,' cnt ',icnt
c$$$         print*,' INTGRT-R: rank t',rank,time,' ir ',ir,' recv ',istrec,
c$$$     *        ' thru ',istrec+icnt2-1,' fr ',irecv,' cnt2 ',icnt2
c$$$         call flush(6)
c$$$      end if
*
#ifdef PUREMPI
               call cputim(tta)
               CALL MPI_SENDRECV(YMPI(1,istsen),41*icnt,MPI_REAL8,isend,
     *              rank,YMPI(1,istrec),41*icnt2,MPI_REAL8,irecv,irecv,
     *              MPI_COMM_WORLD,status,ierr)
*     
c$$$  call mpi_barrier(MPI_COMM_WORLD,ierr)
               CALL MPI_SENDRECV(IMPI(1,istsen),lmax*icnt,MPI_INTEGER,
     *              isend,rank,IMPI(1,istrec),lmax*icnt2,MPI_INTEGER,
     *              irecv,irecv,MPI_COMM_WORLD,status,ierr)
               call cputim(ttb)
               call mpi_barrier(MPI_COMM_WORLD,ierr)
               call cputim(tt999)
               ibarcount=ibarcount+1
               ttbar = ttbar + (tt999-ttb)*60
               ttbarreg = ttbarreg +(tt999-ttb)*60
               xtsub2 = xtsub2 + dble((41*8+lmax*4)*(icnt+icnt2))
               ttsub2 = ttsub2 + (ttb-tta)*60.
#endif
*
#ifdef SHMEM
*
               call barrier()
               call shmem_get(YMPI(1,istrec),YMPI(1,istrec),41*icnt2,
     *              irecv)
               call shmem_get(IMPI(1,istrec),IMPI(1,istrec),lmax*icnt2,
     *              irecv)
#endif
*
            END DO
*
            CALL MPI_ALLGATHER(JMPI(1),11,MPI_INTEGER,
     *           JMPILOCAL(1,1),11,MPI_INTEGER,MPI_COMM_WORLD,ierr)
*     
            call cputim(tt998)
            call mpi_barrier(MPI_COMM_WORLD,ierr)
            call cputim(tt999)
            ibarcount=ibarcount+1
            ttbar = ttbar + (tt999-tt998)*60

!$omp parallel do if(NREG.GE.ITHREAD) private(L,I,NNB)
            DO L = 1,NREG
               I = IREG(L)
               X0(1:3,I) = YMPI(1:3,L)
*               X(1:3,I) = YMPI(1:3,L)
               X0DOT(1:3,I) = YMPI(4:6,L)
*               XDOT(1:3,I) = YMPI(4:6,L)
               FI(1:3,I) = YMPI(7:9,L)
               FR(1:3,I) = YMPI(10:12,L)
               FIDOT(1:3,I) = YMPI(13:15,L)
               FRDOT(1:3,I) = YMPI(16:18,L)
               D0(1:3,I) = FI(1:3,I)
               D0R(1:3,I) = FR(1:3,I)
               D1R(1:3,I) = YMPI(19:21,L)
               D2R(1:3,I) = YMPI(22:24,L)
               D3R(1:3,I) = YMPI(25:27,L)
*     Save corrections of force polynomials from routine fpcorr
               IF (KZ(38).GT.0.OR.I.GT.N) THEN
                  D1(1:3,I) = YMPI(28:30,L)
                  D2(1:3,I) = YMPI(31:33,L)
                  D3(1:3,I) = YMPI(34:36,L)
               END IF
               STEPR(I) = YMPI(37,L)
               STEP(I) = YMPI(38,L)
               TIMENW(I) = YMPI(39,L)
               RS(I) = YMPI(40,L)
               T0R(I) = YMPI(41,L)
               F(1,I) = 0.5D0*(FI(1,I) + FR(1,I))
               F(2,I) = 0.5D0*(FI(2,I) + FR(2,I))
               F(3,I) = 0.5D0*(FI(3,I) + FR(3,I))               
               FDOT(1,I) = ONE6*(FIDOT(1,I) + FRDOT(1,I))
               FDOT(2,I) = ONE6*(FIDOT(2,I) + FRDOT(2,I))
               FDOT(3,I) = ONE6*(FIDOT(3,I) + FRDOT(3,I))
*     Copy neighbor list
               NNB = IMPI(1,L) + 1 
               LIST(1:NNB,I) = IMPI(1:NNB,L)
               IF(LIST(1,I).GT.0) THEN
!$omp critical                  
                  RSMIN = MIN(RSMIN,RS(I))
!$omp end critical
               END IF
#ifdef SIMD               
               CALL IRR_SIMD_SET_LIST(I,LIST(1,I))
#endif
            END DO
!$omp end parallel do       
*
            NBFULL = NC5
            NBVOID = NC6
            NBSMIN = NC11
            NBDIS2 = NC15
*     NRCONV = NC18
            NBDIS = NC19
            NLSMIN = NC20
            NFAST = NC30
            NICONV = NCX
            NBFAST = NCY
            NBFLUX = NBSUM
            NBPRED = NC12
            DO J=1,isize
               NBFULL = NBFULL + JMPILOCAL(1,J)
               NBVOID = NBVOID + JMPILOCAL(2,J)
               NBSMIN = NBSMIN + JMPILOCAL(3,J)
               NBDIS2 = NBDIS2 + JMPILOCAL(4,J)
               NBDIS = NBDIS + JMPILOCAL(5,J)
               NLSMIN = NLSMIN + JMPILOCAL(6,J)
               NFAST = NFAST + JMPILOCAL(7,J)
               NICONV = NICONV + JMPILOCAL(8,J)
               NBFAST = NBFAST + JMPILOCAL(9,J)
               NBFLUX = NBFLUX + JMPILOCAL(10,J)
               NBPRED = NBPRED + JMPILOCAL(11,J)
            END DO
            call cputim(tt3)
            ttmov = ttmov + (tt3-tt999)*60.
*
*      call mpi_barrier(MPI_COMM_WORLD,ierr)
*      call cputim(tt999)
*      ibarcount=ibarcount+1
*      ttbar = ttbar + (tt999-tt3)*60
         END IF
*         End PARALLEL section for regint
#endif
*
c$$$         call cputim(tt333) 
         NSTEPR = NSTEPR + NREG
         NBLCKR = NBLCKR + 1
c$$$*
c$$$c$$$**!$omp parallel do if(NREG.GE.ithread) private(L,I,K)      
c$$$c$$$         DO L = 1,NREG
c$$$c$$$            I = IREG(L)
c$$$c$$$*     
c$$$c$$$            LIST(1:LMAX,I) = IMPI(1:LMAX,L)
c$$$c$$$*
c$$$c$$$*       Check minimum neighbor sphere since last output
c$$$c$$$**!$omp critical          
c$$$c$$$            IF(LIST(1,I).GT.0)RSMIN = MIN(RSMIN,RS(I))
c$$$c$$$**!$omp end critical              
c$$$c$$$*
c$$$c$$$         END DO
c$$$c$$$**!$omp end parallel do      
c$$$         call cputim(tt334)
c$$$         ttmov = ttmov +(tt334-tt333)*60
*
*          OPEN(98,STATUS='OLD',ERR=123)
*         print*,' last reg block t=',time,' length=',nreg
*         print*,' first 10 =',(name(ireg(l)),l=1,min(nreg,10))
*         call flush(6)
*         CLOSE(98)
*123      CONTINUE
*
      END IF

*     Send corrected active particles to GPUIRR library.
#ifdef SIMD
      call cputim(ttirrs1)
!$omp parallel do private(I, L)
      DO L = 1,NXTLEN
         I = NXTLST(L)
         CALL IRR_SIMD_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &        BODY(I),T0(I))
      END DO
!$omp end parallel do
      call cputim(ttirrs2)
      ttsimdsend = ttsimdsend + (ttirrs2-ttirrs1)*60.
#endif

*     Resorting the NXTLST
      call cputim(tt335)
      call sort_tlist(STEP,DTK,.false.)
      call cputim(tt336)
      ttintb =ttintb + (tt336-tt335)*60
      ttnewt = ttnewt + (tt336 - tt335)*60.      
*     --07/08/14 15:38-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      DO L = 1, NDTMAX
c$$$         IF(NDTK(L).LT.NDTK(L+1)) then
c$$$            print*,'Error: L',L,'NDTK(L)',NDTK(L),'NDTK(L+1)',NDTK(L+1)
c$$$            call flush(6)
c$$$            call abort()
c$$$         END IF
c$$$      END DO
c$$$      IF(NGHOSTS.GT.0) THEN
c$$$         DO L = 1, NGHOSTS
c$$$            J = NXTLST(NXTLIMIT+L)
c$$$            IF(STEP(J).LE.DTK(1)) THEN
c$$$               print*,'Error: Ghost particle with normal step! J',J,
c$$$     &              'L',L,'STEP',STEP(J),'T',TIME,'NGHOSTS',
c$$$     &              NXTLST(NXTLIMIT+1:NXTLIMIT+NGHOSTS)
c$$$               call flush(6)
c$$$               call abort()
c$$$            END IF
c$$$         END DO
c$$$      END IF
c$$$      J = 751
c$$$      K=K_STEP(STEP(J),DTK)
c$$$      DO L = NDTK(K+1)+1, NDTK(K)
c$$$         IF(NXTLST(L).EQ.J) THEN
c$$$            print*,'J',J,'N',NAME(J),'STEP',STEP(J),'T',TIME,'NXTLEVEL',
c$$$     &           NXTLEVEL,'L',L,'DTK',DTK(K),'K',K,'IPHASE',IPHASE
c$$$            GO TO 1022
c$$$         end if
c$$$      END DO
c$$$      print*,'J',J,'N',NAME(J),'NOT FIND STEP',STEP(J),'K',K,'NXTLEVEL',
c$$$     &     NXTLEVEL,'IPHASE',IPHASE
c$$$      call flush(6)
c$$$C      call abort()
c$$$ 1022 continue
*     --07/08/14 15:38-lwang-end----------------------------------------*
*
#ifdef H5OUTPUT
*     HDF5 format output (Maxwell)
      IF(KZ(46).GT.0) THEN
         IF(rank.eq.0) CALL custom_output(NXTLST,NXTLEN,1)
      END IF
#endif

*     Output of Massive black hole data
      IF (KZ(24).EQ.1) then
         call BHPLOT
      END IF
c$$$*     --10/27/13 10:33-lwang-add----------------------------------------*
c$$$***** if ks need terminate, do prediction for icomp and jcomp-----------*
c$$$      if (iks.gt.0) then
c$$$         call xvpred(ICOMP,0)
c$$$         call xvpred(JCOMP,0)
c$$$      end if
c$$$*     --10/27/13 10:33-lwang-end----------------------------------------*
*
**!$omp parallel do if (nreg.ge.ithread) private(L,I)      
*       DO L =1,NREG
*          I = IREG(L)
*          DO 82 K = 1,3
*             X0(K,I) = XN(K,I)
*             X0DOT(K,I) = XNDOT(K,I)
* 82       CONTINUE
*       END DO
**!$omp end parallel do      
*     Needed to update only active particles, avoid sending all at reg block
*     Take care is this really efficient?
*         CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
*    &                                            BODY(I),T0(I))
      
*
*       Exit on KS termination, new multiple regularization or merger.
      IF (IQ.GT.0) THEN
          NBPREV = 0
          IF (IQ.GE.4.AND.IQ.NE.7) THEN
              CALL DELAY(IQ,-1)
          ELSE
*       Ensure correct KS index (KSPAIR may denote second termination).
              KSPAIR = KVEC(I10)
              IPHASE = IQ
          END IF
          GO TO 100
      END IF
*
*       Perform optional check on high-velocity particles at major times.
      IF (LISTV(1).GT.0) THEN
          IF (DMOD(TIME,STEPM).EQ.0.0D0) THEN
              CALL SHRINK(TMIN)
              IF (LISTV(1).GT.0) THEN
                  CALL HIVEL(-1)
              END IF
          END IF
      END IF
*
*       Check optional mass loss time.
*     --09/27/13 10:39-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      IF (NREG.GT.0.or.(TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0))
c$$$     &     THEN
c$$$         if(TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0) then
c$$$            NGPUC = NGPUC + 1
c$$$         end if
c$$$         if(rank.eq.0) then
c$$$            write(140,*) NREG,NGPUC,TIME
c$$$            call flush(140)
c$$$         end if
c$$$      end if
c$$$      if(time.ge.3.5373535156250000) then
c$$$      ii=3173
c$$$*      print*,rank,'bm',ii,name(ii),f(1,ii),fdot(1,ii),body(ii),time
c$$$      print*,rank,'bm',ii,name(ii),d3(1,ii),d2r(1,ii),d3r(1,ii),time
c$$$      call flush(6)
c$$$      end if
*     --09/27/13 10:39-lwang-end----------------------------------------*
      IF (KZ(19).GT.0) THEN
*       Delay until time commensurate with 1000-year step (new polynomials).
          IF (NREG.NE.0.AND.TIME.GT.TMDOT
     *        .AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
c$$$          IF (TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
             call cputim(tttmdota)
             IF (KZ(19).GE.3) THEN
                CALL MDOT
             ELSE
                CALL MLOSS
             END IF
             call cputim(tttmdotb)
             ttmdot = ttmdot +(tttmdotb-tttmdota)*60
             IF (IPHASE.LT.0.OR.IPHASE.EQ.1.OR.IPHASE.EQ.2) GO TO 999
          END IF
      END IF
*
*     --09/27/13 10:39-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if (time.ge.3.5373535156250000) then
c$$$      ii=3173
c$$$*      print*,rank,'bm',ii,name(ii),f(1,ii),fdot(1,ii),body(ii),time
c$$$      print*,rank,'bm',ii,name(ii),d3(1,ii),d2r(1,ii),d3r(1,ii),time
c$$$      call flush(6)
c$$$      end if
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)      
*     --09/27/13 10:39-lwang-end----------------------------------------*
*       Advance counters and check timer & optional COMMON save (NSUB = 0).
      NTIMER = NTIMER + NXTLEN
      IF (NTIMER.LT.NMAX) GO TO 1
      NTIMER = 0
      NSTEPS = NSTEPS + NMAX
*
      IF (NSTEPS.GE.100*NMAX.AND.NSUB.EQ.0) THEN
         NSTEPS = 0
         IF (KZ(1).GT.1) CALL MYDUMP(1,1)
      END IF
*     
C*     Check option for general binary search.
C      IF (KZ(4).NE.0.AND.TIME - TLASTS.GT.DELTAS) THEN
C         CALL EVOLVE(0,0)
C      END IF
*     
*     Include facility for termination of run (create dummy file STOP).
      IF(rank.EQ.0)THEN
         OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
         IF (IO.EQ.0) THEN
            CLOSE (99)
            IF (NSUB.EQ.0.and.rank.eq.0) WRITE (6,90)
 90         FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
            CPU = 0.0
         END IF
      END IF
*     
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      TCOMP = (TCOMP-TTOTA)*60.
*     
      IF (TCOMP.LT.CPU) GO TO 1
*     
*     Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
*     Specify zero step to enforce termination.
         STEPS(1:NSUB) = 0.0D0
         NTIMER = NMAX
         GO TO 1
      END IF
*     
*     Terminate run with optional COMMON save.
      IF (KZ(1).NE.0) THEN
         CPUTOT = CPUTOT + TCOMP - CPU0
         CALL MYDUMP(1,1)
         if(rank.eq.0)
     &        WRITE (6,98)  TOFF, TIME, TIME+TOFF, TCOMP, CPUTOT/60.0, 
     &        ERRTOT, DETOT
 98      FORMAT (//,9X,'COMMON SAVED AT TOFF/TIME/TTOT =',1P,E16.8,
     &        '  TCOMP =',E16.8,'  CPUTOT =',E16.8,
     &                  '  ERRTOT =',F10.6,'  DETOT =',F10.6)
      END IF
*
*     Determine time interval and step numbers per time unit
      TIMINT = TIME + TOFF - TINIT
*
#ifdef PARALLEL
      IF(rank.EQ.0)THEN
#endif
         WRITE (6,195)  rank,TIMINT,NSTEPI-NIR,NSTEPB-NIB,NSTEPR-NRGL,
     &        NSTEPU-NKS
 195     FORMAT (//,I9,' INTEGRATION INTERVAL =',F8.2,3X,' NIRR=',I11,
     &        ' NIRRB=',I11,' NREG=',I11,' NKS=',I11)
         WRITE (6,196)  (NSTEPI-NIR)/TIMINT,(NSTEPB-NIB)/TIMINT,
     &        (NSTEPR-NRGL)/TIMINT,(NSTEPU-NKS)/TIMINT
 196     FORMAT (//,9X,' PER TIME UNIT: NIRR=',1P,D12.5,' NIRRB=',
     &        D12.5,' NREG=',D12.5,' NKS=',D12.5)
#ifdef PARALLEL
      END IF
#endif
#ifdef GPU
      CALL GPUNB_CLOSE
*      CALL GPUPOT_CLOSE
*      CALL GPUPOT_CLOSE_FLOAT
#endif
#ifdef SIMD
      CALL IRR_SIMD_CLOSE(rank)
#endif
#ifdef PARALLEL
      call cputim(tt998)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call cputim(tt999)
      ibarcount=ibarcount+1
      ttbar = ttbar + (tt999-tt998)*60
      CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
#endif
      STOP
*
 100  CONTINUE
*
*     Set current global time.
      TTOT = TIME + TOFF

      RETURN
*
      END
