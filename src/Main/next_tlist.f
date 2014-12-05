      SUBROUTINE next_tlist(TMIN,DTK,T0)
*
*      
*     Determine next nxtlen
*     -----------------------------------------
*
*     TMIN:     Next integrating time
*
      include 'params.h'
      include 'tlist.h'
      REAL*8 TMIN,DTK(64),T0(NMAX)
*
*     New Block Time
      TMIN = T0(NXTLST(1)) + DTK(NDTMAX)
*     Do not miss larger step level
      IF(DMOD(TMIN,DTK(NXTLEVEL)).EQ.0) THEN
 1       IF(NXTLEVEL.GT.NDTMIN) THEN
            IF(DMOD(TMIN,DTK(NXTLEVEL-1)).EQ.0) THEN
               NXTLEVEL = NXTLEVEL - 1
               GO TO 1
            END IF
         END IF
*     Check whether need go to smaller step level
      ELSE
 2       IF(NXTLEVEL.LT.NDTMAX) THEN
            NXTLEVEL = NXTLEVEL + 1
            IF(DMOD(TMIN,DTK(NXTLEVEL)).NE.0) GO TO 2
         ELSE
            write(6,*) 'Error: Smallest time step reached! TMIN =',TMIN,
     &           'NDTMAX',NDTMAX,'DTK(NDTMAX)',DTK(NDTMAX),'NXTLST(1)',
     &           NXTLST(1),'T0(1)',T0(NXTLST(1))
            call flush(6)
            call abort()
         END IF
      END IF
*     Next block list range
      NXTLEN = NDTK(NXTLEVEL)
      
      RETURN

      END
