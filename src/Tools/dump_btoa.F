      PROGRAM DUMP_BTOA
*
*
*     Small routine transfor fort.1 or fort.2 from binary format to ANSI format
*
      IMPLICIT REAL*8  (A-H,O-Z)
      INCLUDE 'params.h'
      INCLUDE 'tlist.h'
#ifdef TT
      INCLUDE 'tt.h'
#endif

      INTEGER IUNIT, IRORW
      CHARACTER*10 FILE_I
*     IUNIT = file unit
*     IRORW = 1 (read binary, write ANSI) 2 (read ANSI, write binary)
*     Notice!!: Should be careful the params.h is used 

*     NO (MCL)
*     NP (NTMAX)
*     NQ (MMAX)
*     NS (MMAX)
      
      PARAMETER  (NA=84,NB=168,NC=265,NE=12,NM=20,NL=99,
     &     ND1=52,ND2=288)
      INTEGER*4  C_NAMES,C_COUNTS,C_STARSI,C_RAND2
      REAL*8  C_PARAMS,C_STARSR,C_PLPOT,C_GALAXY

      INTEGER I_NMAX,I_KMAX,I_LMAX,I_MMAX,I_MLD,I_MLR,I_MLV,I_MCL
      INTEGER I_NCMAX,I_NTMAX,I_maxpe,I_ithread
      INTEGER NPARTMP(10)
      
      COMMON/NAMES/  NTOT,NPAIRS,NTTOT,C_NAMES(NA)
      COMMON/COUNTS/ C_COUNTS(NB)
      COMMON/PARAMS/ C_PARAMS(NC) 
      COMMON/STARS/  C_STARSR(ND1),C_STARSI(ND2),LISTR(MLR),LISTD(MLD),
     &               LISTV(MLV)
      COMMON/PLPOT/  C_PLPOT(NE)
      COMMON/BLOCKS/ TPREV,TBLOCK,DTK(64),KVEC(2*KMAX)
      COMMON/RAND2/  C_RAND2(NL)
      COMMON/GALAXY/ C_GALAXY(NM)
      COMMON/CLOUDS/ XCL(3,MCL),XDOTCL(3,MCL),BODYCL(MCL),RCL2(MCL),
     &               CLM(MCL),CLMDOT(MCL),CLDOT,VCL,SIGMA,RB2,PCL2,
     &               TCL,STEPCL,NCL,NEWCL
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      COMMON/RCHE/   CMRCH(13,MMAX),NAMER(2,MMAX),KSTARR(3,MMAX)
      COMMON/BINARY/ ZM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &               HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &               NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      COMMON/NBODY/  X(3,NMAX),XDOT(3,NMAX),X0(3,NMAX),X0DOT(3,NMAX),
     &               F(3,NMAX),FDOT(3,NMAX),BODY(NMAX),RS(NMAX),
     &               FI(3,NMAX),D1(3,NMAX),D2(3,NMAX),D3(3,NMAX),
     &               FR(3,NMAX),D1R(3,NMAX),D2R(3,NMAX),D3R(3,NMAX),
     &               STEP(NMAX),T0(NMAX),STEPR(NMAX),T0R(NMAX),
     &               TIMENW(NMAX),RADIUS(NMAX),TEV(NMAX),TEV0(NMAX),
     &               BODY0(NMAX),EPOCH(NMAX),SPIN(NMAX),XSTAR(NMAX),
     &               ZLMSTY(NMAX),FIDOT(3,NMAX),D0(3,NMAX),
     &               FRDOT(3,NMAX),D0R(3,NMAX),KSTAR(NMAX),NAME(NMAX),
     &               LIST(LMAX,NMAX)
      COMMON/PAIRS/  U(4,KMAX),U0(4,KMAX),UDOT(4,KMAX),FU(4,KMAX),
     &               FUDOT(4,KMAX),FUDOT2(4,KMAX),FUDOT3(4,KMAX),
     &               H(KMAX),HDOT(KMAX),HDOT2(KMAX),HDOT3(KMAX),
     &               HDOT4(KMAX),DTAU(KMAX),TDOT2(KMAX),TDOT3(KMAX),
     &               R(KMAX),R0(KMAX),GAMMA(KMAX),SF(7,KMAX),H0(KMAX),
     &               FP0(4,KMAX),FD0(4,KMAX),KBLIST(10*KMAX),
     &               KSLOW(KMAX),TBLIST

*     READ parameters:
      WRITE(6,*) 'Please input 2 arguments:'
      WRITE(6,*) 'Unit Number # of binary data (will read fort.#):'
      READ(5,*) IUNIT
      WRITE(6,*) 'Transfer direction (1: bianry->ANSI; 2: ANSI->binary)'
      READ(5,*) IRORW

*     Reading binary and write ANSI
      IF(IRORW.EQ.1) GO TO 1
*     Reading ANSI and write binary
      IF(IRORW.EQ.2) GO TO 2

*     FILE name
 20   format('fort.',I1,I4)

*     Reading binary
 1    WRITE(FILE_I,20) IUNIT
      write(6,*) 'Reading file: ',FILE_I,'; Writing file: dump.dat'
      OPEN (UNIT=IUNIT,STATUS='OLD',ACTION='READ',FORM='UNFORMATTED',
     &     FILE=FILE_I)

      READ (IUNIT) NPARTMP
*     Safe check

      IF(NPARTMP(3).NE.INT(LMAX)) then
         write(6,*) 'Error! LMAX mismatch (Data:',NPARTMP(3),
     &        'Current:',LMAX
         call flush(6)
         call abort()
      END IF
      
      IF(NPARTMP(4).NE.INT(MMAX)) then
         write(6,*) 'Error! MMAX mismatch (Data:',NPARTMP(4),
     &        'Current:', MMAX
         call flush(6)
         call abort()
      END IF
      
      IF(NPARTMP(5).NE.INT(MLD)) then
         write(6,*) 'Error! MLD mismatch (Data:',NPARTMP(5),
     &        'Current:', MLD
         call flush(6)
         call abort()
      END IF
      
      IF(NPARTMP(6).NE.INT(MLR)) then
         write(6,*) 'Error! MLR mismatch (Data:',NPARTMP(6),
     &        'Current:', MLR
         call flush(6)
         call abort()
      END IF

      IF(NPARTMP(7).NE.INT(MLV)) then
         write(6,*) 'Error! MLV mismatch (Data:',NPARTMP(7),
     &        'Current:', MLV
         call flush(6)
         call abort()
      END IF
      
      IF(NPARTMP(8).NE.INT(MCL)) then
         write(6,*) 'Error! MCL mismatch (Data:',NPARTMP(8),
     &        'Current:',MCL
         call flush(6)
         call abort()
      END IF
      
      IF(NPARTMP(9).NE.INT(NCMAX)) then
         write(6,*) 'Error! NCMAX mismatch (Data:',NPARTMP(9),
     &        'Current:',NCMAX
         call flush(6)
         call abort()
      END IF
      
      IF(NPARTMP(10).NE.INT(NTMAX)) then
         write(6,*) 'Error! NTMAX mismatch (Data:',NPARTMP(10),
     &        'Current:' ,NTMAX
         call flush(6)
         call abort()
      END IF
      
      READ (IUNIT) ntot,npairs,nttot,
*     NAMES,COUNTS,PARAMS
     &     C_NAMES,C_COUNTS,C_PARAMS,
*     STARS
     &     C_STARSR,C_STARSI,LISTR,LISTD,LISTV,
*     PLPOT
     &     C_PLPOT,
*     BLOCKS
     &     TPREV,TBLOCK,DTK,
*     RAND2,GALAXY
     &     C_RAND2,C_GALAXY,
*     CLOUDS
     &     XCL,XDOTCL,BODYCL,
     &     RCL2,CLM,CLMDOT,CLDOT,VCL,SIGMA,RB2,PCL2,
     &     TCL,STEPCL,NCL,NEWCL,
*     MODES
     &     EB0,ZJ0,ECRIT,
     &     AR,BR,EOSC,EDEC,TOSC,
     &     RP,ES,CM,IOSC,NAMEC,
*     RCHE
     &     CMRCH,NAMER,KSTARR,
*     BINARY
     &     ZM,XREL,VREL,HM,UM,UMDOT,TMDIS,
     &     NAMEM,NAMEG,KSTARM,IFLAG
*
      if (ntot.gt.nmax) then
         write (*,*) "DANGER NTOT ",ntot," > NMAX ",nmax,"!"
         call flush(6)
         stop
      end if

      if (npairs.gt.kmax) then
         write (*,*) "DANGER NPAIRS ",npairs," > KMAX ",kmax,"!"
         call flush(6)
         stop
      end if 

      NTSAVE = NTOT
      IF (NTTOT.GT.0) THEN
         NTOT = NTTOT
      END IF
      read (IUNIT) ((x(k,i),k=1,3),i=1,ntot),((x0(k,i),k=1,3),i=1,ntot),
     *     ((x0dot(k,i),k=1,3),i=1,ntot),((xdot(k,i),k=1,3),i=1,ntot),
     *     ((f(k,i),k=1,3),i=1,ntot),((fdot(k,i),k=1,3),i=1,ntot),
     *     (body(i),i=1,ntot),(rs(i),i=1,ntot),
     *     ((fi(k,i),k=1,3),i=1,ntot),((d1(k,i),k=1,3),i=1,ntot),
     *     ((d2(k,i),k=1,3),i=1,ntot),((d3(k,i),k=1,3),i=1,ntot),
     *     ((fr(k,i),k=1,3),i=1,ntot),((d1r(k,i),k=1,3),i=1,ntot),
     *     ((d2r(k,i),k=1,3),i=1,ntot),((d3r(k,i),k=1,3),i=1,ntot),
     *     (step(i),i=1,ntot),(t0(i),i=1,ntot),(stepr(i),i=1,ntot),
     *     (t0r(i),i=1,ntot),(timenw(i),i=1,ntot),(radius(i),i=1,ntot),
     *     (tev(i),i=1,ntot),
     *     (tev0(i),i=1,ntot),(body0(i),i=1,ntot),(epoch(i),i=1,ntot),
     *     (spin(i),i=1,ntot),(xstar(i),i=1,ntot),(zlmsty(i),i=1,ntot),
     *     ((fidot(k,i),k=1,3),i=1,ntot),((d0(k,i),k=1,3),i=1,ntot),
     *     ((frdot(k,i),k=1,3),i=1,ntot),((d0r(k,i),k=1,3),i=1,ntot),
     *     (kstar(i),i=1,ntot),(name(i),i=1,ntot)
*
      read (IUNIT) ((u(k,i),k=1,4),i=1,npairs),((u0(k,i),k=1,4),i=1,
     *     npairs),((udot(k,i),k=1,4),i=1,npairs),((fu(k,i),k=1,4),i=1,
     *     npairs),((fudot(k,i),k=1,4),i=1,npairs),((fudot2(k,i),k=1,4),
     *     i=1,npairs),((fudot3(k,i),k=1,4),i=1,npairs),(h(i),i=1,
     *     npairs),(hdot(i),i=1,npairs),(hdot2(i),i=1,npairs),
     *     (hdot3(i),i=1,npairs),(hdot4(i),i=1,npairs),(dtau(i),
     *     i=1,npairs),(tdot2(i),i=1,npairs),(tdot3(i),i=1,npairs),
     *     (r(i),i=1,npairs),(r0(i),i=1,npairs),(gamma(i),i=1,npairs),
     *     ((sf(k,i),k=1,7),i=1,npairs),(h0(i),i=1,npairs),((fp0(k,i),
     *     k=1,4),i=1,npairs),((fd0(k,i),k=1,4),i=1,npairs),
     *     (kblist(i),i=1,10*npairs),(kslow(i),i=1,npairs),tblist
*
      read (IUNIT) (list(1,i),(list(k,i),k=2,list(1,i)+2),i=1,ntot)

      read (IUNIT) NXTLIMIT,NGHOSTS,NXTLST(1:NXTLIMIT+NGHOSTS),NXTLEN,
     *       NDTK(1:64),NDTMIN,NDTMAX,NXTLEVEL,NLSTDELAY(1),
     *       (NLSTDELAY(K),K=2,NLSTDELAY(1)+1)

#ifdef TT
*** FlorentR
      read (IUNIT) ttunit, nbtt, ttmode
      read (IUNIT) (((ttens(k,i,kk),k=1,3),i=1,3),kk=1,nbtt),
     *    ((tteff(k,i),k=1,3),i=1,3),((dtteff(k,i),k=1,3),i=1,3),
     *    (tttime(i),i=1,nbtt)
#endif

*     --06/21/14 10:41-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
*      print*,'TPREV,TBLOCK,DTK,KEV,',TPREV,TBLOCK,DTK(64),NAMEM(1)
*     --06/21/14 10:41-lwang-end----------------------------------------*

*     Write ANSI
      OPEN (UNIT=10,STATUS='NEW',ACTION='WRITE',FORM='FORMATTED',
     &     FILE='dump.dat')

*     Write all params.h
      write(10,*) NMAX,KMAX,LMAX,MMAX,MLD,MLR,MLV,MCL,NCMAX,NTMAX,
     &     maxpe,ithread

*     Write as mydump.F
      write(10,*) ntot,npairs,nttot,
     &     C_NAMES,C_COUNTS,C_PARAMS,C_PLPOT,C_RAND2,C_GALAXY,
*     BLOCKS
     &     TPREV,TBLOCK,DTK

*     STARS
      write(10,*) C_STARSR,C_STARSI,LISTR,LISTD,LISTV

*     CLOUDS
      write(10,*) XCL,XDOTCL,BODYCL,
     &     RCL2,CLM,CLMDOT,CLDOT,VCL,SIGMA,RB2,PCL2,
     &     TCL,STEPCL,NCL,NEWCL

*     MODES
      write(10,*) EB0,ZJ0,ECRIT,AR,BR,EOSC,EDEC,TOSC,
     &     RP,ES,CM,IOSC,NAMEC

*     RCHE
      write(10,*) CMRCH,NAMER,KSTARR

*     BINARY
      write(10,*) ZM,XREL,VREL,HM,UM,UMDOT,TMDIS,
     &     NAMEM,NAMEG,KSTARM,IFLAG
 
      IF (NTTOT.GT.0) THEN
         NTOT = NTTOT
      END IF
*     NBODY
      write(10,*) x(1:3,1:ntot),x0(1:3,1:ntot),
     *     xdot(1:3,1:ntot),x0dot(1:3,1:ntot)
      write(10,*) f(1:3,1:ntot),fdot(1:3,1:ntot),body(1:ntot),
     *     rs(1:ntot)
      write(10,*) fi(1:3,1:ntot),d1(1:3,1:ntot),d2(1:3,1:ntot),
     *     d3(1:3,1:ntot),fr(1:3,1:ntot),d1r(1:3,1:ntot),
     *     d2r(1:3,1:ntot),d3r(1:3,1:ntot)
      write(10,*) step(1:ntot),t0(1:ntot),stepr(1:ntot),t0r(1:ntot),
     *     timenw(1:ntot)
      write(10,*) radius(1:ntot),tev(1:ntot),tev0(1:ntot),
     *     body0(1:ntot),epoch(1:ntot),spin(1:ntot),
     *     xstar(1:ntot),zlmsty(1:ntot)
      write(10,*) fidot(1:3,1:ntot),d0(1:3,1:ntot),
     *     frdot(1:3,1:ntot),d0r(1:3,1:ntot),
     *     kstar(1:ntot)
*     PAIRS
      write(10,*) u(1:4,1:npairs),u0(1:4,1:npairs),udot(1:4,1:npairs),
     *     fu(1:4,1:npairs),fudot(1:4,1:npairs),fudot2(1:4,1:npairs),
     *     fudot3(1:4,1:npairs)
      write(10,*) h(1:npairs),hdot(1:npairs),hdot2(1:npairs),
     *     hdot3(1:npairs),hdot4(1:npairs)
      write(10,*) dtau(1:npairs),tdot2(1:npairs),tdot3(1:npairs)
      write(10,*) r(1:npairs),r0(1:npairs),gamma(1:npairs),
     *     sf(1:7,1:npairs)
      write(10,*) h0(1:npairs),fp0(1:4,1:npairs),
     *     fd0(1:4,1:npairs),tblist
      write(10,*) kblist(1:10*npairs),kslow(1:npairs),
     *     name(1:ntot)
      write(10,*) NXTLIMIT,NGHOSTS,NXTLST(1:NXTLIMIT+NGHOSTS),NXTLEN,
     *     NDTK(1:64),NDTMIN,NDTMAX,NXTLEVEL,
     *     (NLSTDELAY(K),K=1,NLSTDELAY(1)+1)

*     LIST
      DO i=1,ntot
         write(10,*) list(1:list(1,i)+2,i)
      END DO

#ifdef TT
*** FlorentR
      write (10,*) ttunit, nbtt, ttmode
      write (10,*) (((ttens(k,i,kk),k=1,3),i=1,3),kk=1,nbtt),
     *      ((tteff(k,i),k=1,3),i=1,3),((dtteff(k,i),k=1,3),i=1,3),
     *      (tttime(i),i=1,nbtt)
*** FRenaud
#endif

      END FILE 10
      CLOSE (UNIT=10)

      GO TO 999

*     Reading ANSI
 2    OPEN (UNIT=10,STATUS='OLD',ACTION='READ',FORM='FORMATTED',
     &     FILE='dump.dat')

*     Reading previous params.h
      read(10,*) I_NMAX,I_KMAX,I_LMAX,I_MMAX,I_MLD,I_MLR,I_MLV,I_MCL,
     &     I_NCMAX,I_NTMAX,I_maxpe,I_ithread

      write(6,*) 'Previous params.h:'
      write(6,10) I_NMAX,I_KMAX,I_LMAX,I_MMAX,I_MLD,I_MLR,I_MLV,I_MCL,
     &     I_NCMAX,I_NTMAX,I_maxpe,I_ithread

      write(6,*) 'Current params.h:'
      write(6,10) NMAX,KMAX,LMAX,MMAX,MLD,MLR,MLV,MCL,NCMAX,NTMAX,
     &     maxpe,ithread

 10   format('NMAX=',I10,' KMAX=',I8,' LMAX=',I6,' MMAX=',
     &     I6,' MLD=',I4,' MLR=',I4,' MLV=',I4,' MCL=',I4,
     &     ' NCMAX=',I4,' NTMAX=',I5,' maxpe=',I7,' ithread=',
     &     I4)

*     Reading dump file
      read(10,*) ntot,npairs,nttot,
     &     C_NAMES,C_COUNTS,C_PARAMS,C_PLPOT,C_RAND2,C_GALAXY,
*     BLOCKS
     &     TPREV,TBLOCK,DTK

*     Safe check
      if(ntot.gt.NMAX) then
         write(6,*) 'Error: NTOT (',ntot,')> NMAX (',nmax,')!'
         call flush(6)
         stop
      end if

      if(nttot.gt.NMAX) then
         write(6,*) 'Error: NTTOT (',nttot,')> NMAX (',nmax,')!'
         call flush(6)
         stop
      end if

      if(npairs.gt.KMAX) then
         write(6,*) 'Error: NPAIRS (',npairs,')> KMAX (',KMAX,')!'
         call flush(6)
         stop
      end if

      if(I_LMAX.GT.LMAX) then
         write(6,*) 'Error: Old LMAX (',I_LMAX,')> New LMAX (',LMAX,')!'
         call flush(6)
         stop
      end if
         
      if(I_MMAX.GT.MMAX) then
         write(6,*) 'Error: Old MMAX (',I_MMAX,')> New MMAX (',MMAX,')!'
         call flush(6)
         stop
      end if

      if(I_MLD.GT.MLD) then
         write(6,*) 'Error: Old MLD (',I_MLD,')> New MLD (',MLD,')!'
         call flush(6)
         stop
      end if
      
      if(I_MLR.GT.MLR) then
         write(6,*) 'Error: Old MLR (',I_MLR,')> New MLR (',MLR,')!'
         call flush(6)
         stop
      end if

      if(I_MLV.GT.MLV) then
         write(6,*) 'Error: Old MLV (',I_MLV,')> New MLV (',MLV,')!'
         call flush(6)
         stop
      end if

      if(I_MCL.GT.MCL) then
         write(6,*) 'Error: Old MCL (',I_MCL,')> New MCL (',MCL,')!'
         call flush(6)
         stop
      end if

      if(I_NTMAX.GT.NTMAX) then
         write(6,*) 'Error: Old NTMAX (',I_NTMAX,')> New NTMAX (',NTMAX
     &        ,')!'
         call flush(6)
         stop
      end if

*     STARS
      read(10,*) C_STARSR,C_STARSI,LISTR(1:I_MLR),LISTD(1:I_MLD),
     &     LISTV(1:I_MLV)
      
*     CLOUDS
      read(10,*) XCL(1:3,1:I_MCL),XDOTCL(1:3,1:I_MCL),BODYCL(1:I_MCL),
     &     RCL2(1:I_MCL),CLM(1:I_MCL),CLMDOT(1:I_MCL),CLDOT,VCL,SIGMA,
     &     RB2,PCL2,TCL,STEPCL,NCL,NEWCL

*     MODES
      read(10,*) EB0(1:I_NTMAX),ZJ0(1:I_NTMAX),ECRIT(1:I_NTMAX),
     &     AR(1:I_NTMAX),BR(1:I_NTMAX),EOSC(1:4,1:I_NTMAX),
     &     EDEC(1:I_NTMAX),TOSC(1:I_NTMAX),RP(1:I_NTMAX),ES(1:I_NTMAX),
     &     CM(1:2,1:I_NTMAX),IOSC(1:I_NTMAX),NAMEC(1:I_NTMAX)
      
*     RCHE
      read(10,*) CMRCH(1:13,1:I_MMAX),NAMER(1:2,1:I_MMAX),
     &     KSTARR(1:3,1:I_MMAX)

*     BINARY
      read(10,*) ZM(1:4,1:I_MMAX),XREL(1:3,1:I_MMAX),VREL(1:3,1:I_MMAX),
     &     HM(1:I_MMAX),UM(1:4,1:I_MMAX),UMDOT(1:4,1:I_MMAX),
     &     TMDIS(1:I_MMAX),NAMEM(1:I_MMAX),NAMEG(1:I_MMAX),
     &     KSTARM(1:I_MMAX),IFLAG(1:I_MMAX)
      
      IF (NTTOT.GT.0) THEN
         NTOT = NTTOT
      END IF

*     NBODY
      read(10,*) x(1:3,1:ntot),x0(1:3,1:ntot),
     *     xdot(1:3,1:ntot),x0dot(1:3,1:ntot)
      read(10,*) f(1:3,1:ntot),fdot(1:3,1:ntot),body(1:ntot),
     *     rs(1:ntot)
      read(10,*) fi(1:3,1:ntot),d1(1:3,1:ntot),d2(1:3,1:ntot),
     *     d3(1:3,1:ntot),fr(1:3,1:ntot),d1r(1:3,1:ntot),
     *     d2r(1:3,1:ntot),d3r(1:3,1:ntot)
      read(10,*) step(1:ntot),t0(1:ntot),stepr(1:ntot),t0r(1:ntot),
     *     timenw(1:ntot)
      read(10,*) radius(1:ntot),tev(1:ntot),tev0(1:ntot),
     *     body0(1:ntot),epoch(1:ntot),spin(1:ntot),
     *     xstar(1:ntot),zlmsty(1:ntot)
      read(10,*) fidot(1:3,1:ntot),d0(1:3,1:ntot),
     *     frdot(1:3,1:ntot),d0r(1:3,1:ntot),
     *     kstar(1:ntot)
*     PAIRS
      read(10,*) u(1:4,1:npairs),u0(1:4,1:npairs),udot(1:4,1:npairs),
     *     fu(1:4,1:npairs),fudot(1:4,1:npairs),fudot2(1:4,1:npairs),
     *     fudot3(1:4,1:npairs)
      read(10,*) h(1:npairs),hdot(1:npairs),hdot2(1:npairs),
     *     hdot3(1:npairs),hdot4(1:npairs)
      read(10,*) dtau(1:npairs),tdot2(1:npairs),tdot3(1:npairs)
      read(10,*) r(1:npairs),r0(1:npairs),gamma(1:npairs),
     *     sf(1:7,1:npairs)
      read(10,*) h0(1:npairs),fp0(1:4,1:npairs),
     *     fd0(1:4,1:npairs),tblist
      read(10,*) kblist(1:10*npairs),kslow(1:npairs),
     *     name(1:ntot)
      read(10,*) NXTLIMIT,NGHOSTS,NXTLST(1:NXTLIMIT+NGHOSTS),NXTLEN,
     *     NDTK(1:64),NDTMIN,NDTMAX,NXTLEVEL,
     *     (NLSTDELAY(K),K=1,NLSTDELAY(1)+1)
*     LIST
      DO i=1,ntot
         read(10,*) list(1,i),list(2:list(1,i)+1,i)
      END DO

#ifdef TT
*** FlorentR
      read (10,*) ttunit, nbtt, ttmode
      read (10,*) (((ttens(k,i,kk),k=1,3),i=1,3),kk=1,nbtt),
     *     ((tteff(k,i),k=1,3),i=1,3),((dtteff(k,i),k=1,3),i=1,3),
     *     (tttime(i),i=1,nbtt)
#endif
      
*     Write new binary data
      WRITE(FILE_I,20) IUNIT
      OPEN (UNIT=IUNIT,STATUS='NEW',ACTION='WRITE',FORM='UNFORMATTED',
     &     FILE=FILE_I)
      write(6,*) 'Reading file: dump.dat;  Writing file:',FILE_I

*     INITIAL KVEC
      DO I=1,KMAX
         KVEC(2*I-1) = I
         KVEC(2*I) = I
      END DO
*
      NPARTMP(1) = NMAX
      NPARTMP(2) = KMAX
      NPARTMP(3) = LMAX
      NPARTMP(4) = MMAX
      NPARTMP(5) = MLD
      NPARTMP(6) = MLR
      NPARTMP(7) = MLV
      NPARTMP(8) = MCL
      NPARTMP(9) = NCMAX
      NPARTMP(10) = NTMAX

      write(IUNIT) NPARTMP
      
      write(IUNIT) ntot,npairs,nttot,
*     NAMES,COUNTS,PARAMS
     &     C_NAMES,C_COUNTS,C_PARAMS,
*     STARS
     &     C_STARSR,C_STARSI,LISTR,LISTD,LISTV,
*     PLPOT
     &     C_PLPOT,
*     BLOCKS
     &     TPREV,TBLOCK,DTK,
*     RAND2,GALAXY
     &     C_RAND2,C_GALAXY,
*     CLOUDS
     &     XCL,XDOTCL,BODYCL,
     &     RCL2,CLM,CLMDOT,CLDOT,VCL,SIGMA,RB2,PCL2,
     &     TCL,STEPCL,NCL,NEWCL,
*     MODES
     &     EB0,ZJ0,ECRIT,
     &     AR,BR,EOSC,EDEC,TOSC,
     &     RP,ES,CM,IOSC,NAMEC,
*     RCHE
     &     CMRCH,NAMER,KSTARR,
*     BINARY
     &     ZM,XREL,VREL,HM,UM,UMDOT,TMDIS,
     &     NAMEM,NAMEG,KSTARM,IFLAG
*     NBODY
      WRITE(IUNIT) ((x(k,i),k=1,3),i=1,ntot),((x0(k,i),k=1,3),i=1,ntot),
     *     ((x0dot(k,i),k=1,3),i=1,ntot),((xdot(k,i),k=1,3),i=1,ntot),
     *     ((f(k,i),k=1,3),i=1,ntot),((fdot(k,i),k=1,3),i=1,ntot),
     *     (body(i),i=1,ntot),(rs(i),i=1,ntot),
     *     ((fi(k,i),k=1,3),i=1,ntot),((d1(k,i),k=1,3),i=1,ntot),
     *     ((d2(k,i),k=1,3),i=1,ntot),((d3(k,i),k=1,3),i=1,ntot),
     *     ((fr(k,i),k=1,3),i=1,ntot),((d1r(k,i),k=1,3),i=1,ntot),
     *     ((d2r(k,i),k=1,3),i=1,ntot),((d3r(k,i),k=1,3),i=1,ntot),
     *     (step(i),i=1,ntot),(t0(i),i=1,ntot),(stepr(i),i=1,ntot),
     *     (t0r(i),i=1,ntot),(timenw(i),i=1,ntot),(radius(i),i=1,ntot),
     *     (tev(i),i=1,ntot),
     *     (tev0(i),i=1,ntot),(body0(i),i=1,ntot),(epoch(i),i=1,ntot),
     *     (spin(i),i=1,ntot),(xstar(i),i=1,ntot),(zlmsty(i),i=1,ntot),
     *     ((fidot(k,i),k=1,3),i=1,ntot),((d0(k,i),k=1,3),i=1,ntot),
     *     ((frdot(k,i),k=1,3),i=1,ntot),((d0r(k,i),k=1,3),i=1,ntot),
     *     (kstar(i),i=1,ntot),(name(i),i=1,ntot)
*     PAIRS
      write (IUNIT) ((u(k,i),k=1,4),i=1,npairs),((u0(k,i),k=1,4),i=1,
     *     npairs),((udot(k,i),k=1,4),i=1,npairs),((fu(k,i),k=1,4),i=1,
     *     npairs),((fudot(k,i),k=1,4),i=1,npairs),((fudot2(k,i),k=1,4),
     *     i=1,npairs),((fudot3(k,i),k=1,4),i=1,npairs),(h(i),i=1,
     *     npairs),(hdot(i),i=1,npairs),(hdot2(i),i=1,npairs),  
     *     (hdot3(i),i=1,npairs),(hdot4(i),i=1,npairs),(dtau(i),
     *     i=1,npairs),(tdot2(i),i=1,npairs),(tdot3(i),i=1,npairs),
     *     (r(i),i=1,npairs),(r0(i),i=1,npairs),(gamma(i),i=1,npairs),
     *     ((sf(k,i),k=1,7),i=1,npairs),(h0(i),i=1,npairs),((fp0(k,i), 
     *     k=1,4),i=1,npairs),((fd0(k,i),k=1,4),i=1,npairs),
     *     (kblist(i),i=1,10*npairs),(kslow(i),i=1,npairs),tblist
*
      write(IUNIT) ((list(k,i),k=1,list(1,i)+2),i=1,ntot)

      write(IUNIT) NXTLIMIT,NGHOSTS,NXTLST(1:NXTLIMIT+NGHOSTS),NXTLEN,
     *     NDTK(1:64),NDTMIN,NDTMAX,NXTLEVEL,
     *     (NLSTDELAY(K),K=1,NLSTDELAY(1)+1)

#ifdef TT
*** FlorentR
      write (IUNIT) ttunit, nbtt, ttmode
      write (IUNIT) (((ttens(k,i,kk),k=1,3),i=1,3),kk=1,nbtt),
     *     ((tteff(k,i),k=1,3),i=1,3),((dtteff(k,i),k=1,3),i=1,3),
     *     (tttime(i),i=1,nbtt)
*** FRenaud
#endif

      END FILE IUNIT
      CLOSE (UNIT=IUNIT)

 999  END
