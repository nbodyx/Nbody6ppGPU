      SUBROUTINE shrink_tlist
*
*
*     Shrink NDTMAX and NDTMIN if necessary
*     
      include 'params.h'
      include 'tlist.h'

*     Check whether need to modify NDTMAX and NDTMIN
 50   IF(NDTK(NDTMAX).EQ.0.AND.NDTMAX.GT.NDTMIN) THEN
         NDTMAX = NDTMAX - 1
         GO TO 50
      END IF

 51   IF(NDTK(NDTMIN).EQ.NDTK(NDTMIN+1)
     &     .AND.NDTK(NDTMIN).EQ.NXTLIMIT) THEN
         NDTMIN = NDTMIN + 1
         GO TO 51
      END IF
      
      RETURN
      
      END
      
