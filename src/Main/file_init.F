      SUBROUTINE FILE_INIT(ISTART)
*
*
*       Opening of Files with proper names.
*       ------------------------------------
*       ISTART=0: Open all other files after reading input data.
*       ISTART=1: Open only unit 1 for mydump.
*       ISTART=2: Open only unit 2 for mydump.
*
      INCLUDE 'common6.h'
*
*
      CHARACTER*10 FILE(100)
      CHARACTER*4 FRAG(100)
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
      FRAG(1) = 'fort'
      FRAG(2) = 'fort'
      FRAG(3) = 'conf'
C      FRAG(4) = 'bdat'
      FRAG(7) = 'lagr'
      FRAG(8) = 'bdat'
*      FRAG(9) = 'bdat'
*
      JFMIN = 1
      JFMAX = 9
      IF (ISTART.EQ.1) THEN
          JFMIN = 1
          JFMAX = 1
      END IF
      IF (ISTART.EQ.2) THEN
          JFMIN = 2
          JFMAX = 2
      END IF
*       Initialize the file names for read and write.
      DO 100 JF = JFMIN,JFMAX
          WRITE (FILE(JF),112) FRAG(JF),JF
 100  CONTINUE
 112  FORMAT(A4,'.',I1,T4)
*      FILE(3)  = 'conf.3'
      FILE(4)  = 'degen.4   '
*      FILE(5)  = 'tail.5'
      FILE(10) = 'dat.10    '
      FILE(11) = 'esc.11    '
      FILE(12) = 'hiarch.12 '
      FILE(13) = 'coll.13   '
      FILE(14) = 'shrink.14 '
      FILE(15) = 'mix.15    '
      FILE(16) = 'hirect.16 '
      FILE(17) = 'binev.17  '
      FILE(18) = 'pbin.18   '
*      FILE(19) = 'bwdat.19  '
      FILE(20) = 'symb.20   '
      FILE(22) = 'rocdeg.22 '
*      FILE(23) = 'ibeigen.23'
      FILE(24) = 'coal.24   '
      FILE(25) = 'sediag.25 '
      FILE(29) = 'highv.29  '
      FILE(30) = 'global.30 '
      FILE(31) = 'lagr1.31  '
      FILE(32) = 'lagr2.32  '
      FILE(33) = 'ns.33     '
      FILE(34) = 'bh.34     '
      FILE(35) = 'event.35  '
      FILE(38) = 'sediag.38 '
      FILE(39) = 'hbin.39   '
*      FILE(40) = 'data.h5par'
      FILE(41) = 'nbflow.41 '
*      FILE(42) = 'ibcoll.42 '
      FILE(43) = 'sediag.43 '
      FILE(44) = 'hinc.44   '
      FILE(45) = 'mbh.45    '
      FILE(46) = 'mbhnb.46  '
*      FILE(52) = 'itid3.52  '
      FILE(54) = 'hypcep.54 '
      FILE(55) = 'hypcec.55 '
      FILE(71) = 'cirdiag.71'
      FILE(73) = 'histab.73 '
      FILE(75) = 'cirdiag.75'
      FILE(77) = 'kscrit.77'
*      FILE(82) = 'bev.82    '
*      FILE(83) = 'sev.83    '
      FILE(81) = 'chstab.81 '
      FILE(84) = 'merger.84 '
      FILE(85) = 'roche.85  '
      FILE(89) = 'quastab.89'
*      FILE(87) = 'hdat.87   '
      FILE(91) = 'bs.91     '
      FILE(95) = 'wdcirc.95 '
      FILE(96) = 'cirdiag.96'
#ifdef ENSEMBLE
      DO 101 JF = JFMIN,JFMAX
      if (isize.lt.10) then
          WRITE (FILE(JF),113) FRAG(JF),rank,JF
 113  FORMAT(A4,'000',I1,'.',I1)
      else if (isize.ge.10.and.isize.lt.100) then
          WRITE (FILE(JF),115) FRAG(JF),rank,JF
 115  FORMAT(A4,'00',I2,'.',I1)
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(JF),117) FRAG(JF),rank,JF
 117  FORMAT(A4,'0',I3,'.',I1)
      else if (isize.ge.1000) then
          WRITE (FILE(JF),119) FRAG(JF),rank,JF
 119  FORMAT(A4,I4,'.',I1)
      end if
 101  CONTINUE
*
      if (isize.lt.10) then
          WRITE (FILE(10),213) rank
          WRITE (FILE(11),313) rank
          WRITE (FILE(12),413) rank
          WRITE (FILE(13),513) rank
          WRITE (FILE(15),613) rank
          WRITE (FILE(82),713) rank
          WRITE (FILE(83),813) rank
 213  FORMAT('dat000',I1,'.10')
 313  FORMAT('esc000',I1,'.11')
 413  FORMAT('hia000',I1,'.12')
 513  FORMAT('hid000',I1,'.13')
 613  FORMAT('per000',I1,'.15')
 713  FORMAT('bev000',I1,'.13')
 813  FORMAT('sev000',I1,'.15')
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(10),215) rank
          WRITE (FILE(11),315) rank
          WRITE (FILE(12),415) rank
          WRITE (FILE(13),515) rank
          WRITE (FILE(15),615) rank
          WRITE (FILE(82),715) rank
          WRITE (FILE(83),815) rank
 215  FORMAT('dat00',I2,'.10')
 315  FORMAT('esc00',I2,'.11')
 415  FORMAT('hia00',I2,'.12')
 515  FORMAT('hid00',I2,'.13')
 615  FORMAT('per00',I2,'.15')
 715  FORMAT('bev00',I2,'.13')
 815  FORMAT('sev00',I2,'.15')
      else if (isize.ge.100.and.isize.lt.1000) then
          WRITE (FILE(10),217) rank
          WRITE (FILE(11),317) rank
          WRITE (FILE(12),417) rank
          WRITE (FILE(13),517) rank
          WRITE (FILE(15),617) rank
          WRITE (FILE(82),717) rank
          WRITE (FILE(83),817) rank
 217  FORMAT('dat0',I3,'.10')
 317  FORMAT('esc0',I3,'.11')
 417  FORMAT('hia0',I3,'.12')
 517  FORMAT('hid0',I3,'.13')
 617  FORMAT('per0',I3,'.15')
 717  FORMAT('bev0',I3,'.13')
 817  FORMAT('sev0',I3,'.15')
      else if (isize.ge.1000) then
          WRITE (FILE(10),219) rank
          WRITE (FILE(11),319) rank
          WRITE (FILE(12),419) rank
          WRITE (FILE(13),519) rank
          WRITE (FILE(15),619) rank
          WRITE (FILE(82),719) rank
          WRITE (FILE(83),819) rank
 219  FORMAT('dat',I4,'.10')
 319  FORMAT('esc',I4,'.11')
 419  FORMAT('hia',I4,'.12')
 519  FORMAT('hid',I4,'.13')
 619  FORMAT('per',I4,'.15')
 719  FORMAT('bev',I4,'.13')
 819  FORMAT('sev',I4,'.15')
      end if
#endif
#ifdef PARALLEL
      if(rank.eq.0) then
#endif
      IF (KZ(1).GT.0.AND.ISTART.EQ.1)
     &OPEN (UNIT=1,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(1))
*
      IF (KZ(2).GT.0.AND.ISTART.EQ.2)
     &OPEN (UNIT=2,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(2))
*
      IF (ISTART.GT.0) RETURN
*
C      IF (KZ(3).GT.0)
C     &OPEN (UNIT=3,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FILE(3),
C     &             ACCESS='APPEND')
      IF (KZ(9).GE.3)
     &OPEN (UNIT=4,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(4),
     &             ACCESS='APPEND')
      IF (KZ(7).GE.2)
     &OPEN (UNIT=7,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(7),
     &             ACCESS='APPEND')
      IF (KZ(8).GT.0.OR.NBIN0.GT.0)
     &OPEN (UNIT=8,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(8),
     &             ACCESS='APPEND')
C      IF (KZ(8).GE.2.OR.NBIN0.GT.0)
C     &OPEN (UNIT=9,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(9),
C     &             ACCESS='APPEND')
      IF (KZ(22).GT.0)
     &OPEN (UNIT=10,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(10))
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4)
     &OPEN (UNIT=11,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(11),
     &             ACCESS='APPEND')
      IF (KZ(18).EQ.1.OR.KZ(18).EQ.3)
     &OPEN (UNIT=12,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(12),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=13,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(13),
     &             ACCESS='APPEND')

      OPEN (UNIT=14,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(14),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=15,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(15),
     &             ACCESS='APPEND')
      IF (KZ(27).EQ.2.or.KZ(34).GT.0.or.KZ(28).GT.0)
     &OPEN (UNIT=16,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(16),
     &             ACCESS='APPEND')
      IF (KZ(9).GE.2)
     &OPEN (UNIT=17,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(17),
     &             ACCESS='APPEND')
      IF (KZ(8).GT.0)
     &OPEN (UNIT=18,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(18),
     &             ACCESS='APPEND')
*      IF (KZ(8).GE.2)
*     &OPEN (UNIT=19,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(19),
*     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=20,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(20),
     &             ACCESS='APPEND')
      IF (KZ(34).GT.0)
     &OPEN (UNIT=22,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(22),
     &             ACCESS='APPEND')
C      IF ((KZ(8).EQ.1).or.(KZ(8).GE.3).AND.KZ(42).EQ.6)
C     &OPEN (UNIT=23,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(23),
C     &             ACCESS='APPEND')
C      IF (KZ(19).GE.3) 
C     &OPEN (UNIT=24,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(24),
C     &             ACCESS='APPEND')
      IF (KZ(27).GT.0) 
     &OPEN (UNIT=25,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(25),
     &             ACCESS='APPEND')

      OPEN (UNIT=29,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(29),
     &             ACCESS='APPEND')

      OPEN (UNIT=30,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(30),
     &             ACCESS='APPEND')
      IF (KZ(7).GE.5)
     &OPEN (UNIT=31,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(31),
     &             ACCESS='APPEND')
      IF (KZ(7).GE.5)
     &OPEN (UNIT=32,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(32),
     &             ACCESS='APPEND')
      IF (KZ(19).GT.0.OR.KZ(27).GT.0) 
     &OPEN (UNIT=35,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(35),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=38,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(38),
     &             ACCESS='APPEND')
      IF (KZ(9).EQ.1.or.KZ(9).EQ.3)
     &OPEN (UNIT=39,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(39),
     &             ACCESS='APPEND')

      OPEN (UNIT=41,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(41),
     &             ACCESS='APPEND')
C      IF (KZ(8).EQ.1.or.KZ(8).GE.3)
C     &OPEN (UNIT=42,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(42),
C     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=43,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(43),
     &             ACCESS='APPEND')
      IF (KZ(27).GT.0)
     &OPEN (UNIT=44,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(44),
     &             ACCESS='APPEND')
      IF (KZ(24).EQ.1)
     &OPEN (UNIT=45,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(45),
     &             ACCESS='APPEND')
      IF (KZ(24).EQ.1)
     &OPEN (UNIT=46,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(46),
     &             ACCESS='APPEND')
C      IF (KZ(14).EQ.3)
C     &OPEN (UNIT=52,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(44),
C     &             ACCESS='APPEND')
      IF (KZ(14).EQ.3)
     &OPEN (UNIT=54,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(54),
     &             ACCESS='APPEND')
      IF (KZ(14).EQ.3)
     &OPEN (UNIT=55,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(55),
     &             ACCESS='APPEND')
      IF (KZ(27).GT.0)
     &OPEN (UNIT=71,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(71),
     &             ACCESS='APPEND')
      IF (KZ(15).GT.0)
     &OPEN (UNIT=73,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(73),
     &             ACCESS='APPEND')
      IF (KZ(27).EQ.2)
     &OPEN (UNIT=75,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(75),
     &             ACCESS='APPEND')
      IF (KZ(16).GT.2.AND.(KZ(30).GT.0.OR.KZ(30).EQ.-1))
     &OPEN (UNIT=77,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(77),
     &             ACCESS='APPEND')
C      IF (KZ(12).GT.0)
C     &OPEN (UNIT=82,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(82),
C     &             ACCESS='APPEND')
C      IF (KZ(12).GT.0)
C     &OPEN (UNIT=83,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(83),
C     &             ACCESS='APPEND')
      IF (KZ(15).GT.0.AND.(KZ(30).GT.0.or.KZ(30).EQ.-1))
     &OPEN (UNIT=81,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(81),
     &             ACCESS='APPEND')
      IF (KZ(8).GE.2)
     &OPEN (UNIT=84,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(84),
     &             ACCESS='APPEND')
      IF (KZ(34).GT.0)
     &OPEN (UNIT=85,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(85),
     &             ACCESS='APPEND')
      IF (KZ(15).GE.3)
     &OPEN (UNIT=89,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(89),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3) 
     &OPEN (UNIT=91,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(91),
     &             ACCESS='APPEND')
      IF (KZ(27).GT.0) 
     &OPEN (UNIT=95,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(95),
     &             ACCESS='APPEND')
      IF (KZ(27).GT.0) 
     &OPEN (UNIT=96,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(96),
     &             ACCESS='APPEND')
#ifdef PARALLEL
      end if
#endif
*
      RETURN
      END
