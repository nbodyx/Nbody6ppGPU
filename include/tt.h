*** FlorentR - New block used for tt treatment
      PARAMETER (NBTTMAX=100000)
      COMMON/TT/     TTENS(3,3,NBTTMAX),TTEFF(3,3),DTTEFF(3,3),
     &               TTTIME(NBTTMAX), TTUNIT, NBTT, TTMODE,
     &               IKEYPOT(11)
*** FRenaud
      INTEGER TTMODE,NBTT,IKEYPOT
      REAL*8 TTENS,TTEFF,DTTEFF,TTTIME,TTUNIT