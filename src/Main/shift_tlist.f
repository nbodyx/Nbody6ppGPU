      SUBROUTINE shift_tlist(I,DT,DL,DTK)
*     
*     
*     Shift particle in block list due to LIST_SF
*     -----------------------------------
*     
*     I: particle index
*     DT: particle irregular step
*     DL: shifting step level

      include 'params.h'
      include 'tlist.h'
      REAL*8 DT,DTK(64)
      INTEGER I,DL
      INTEGER L,J,NL

*     Current step level
      L = K_STEP(DT,DTK)
*     escape for ghost case
      IF(L.LT.0) RETURN
*     Particle position in NXTLST
      J = NDTK(L+1)+1
 1    IF(NXTLST(J).NE.I) THEN
         IF(J.LT.NDTK(L)) THEN
            J = J + 1
            GO TO 1
         ELSE 
            IF(NLSTDELAY(1).GT.0) THEN
               J = 2
*     If find in NLSTDELAY, just return
 2             IF (NLSTDELAY(J).EQ.I) RETURN
               J = J + 1
               IF (J.LE.NLSTDELAY(1)+1) GO TO 2
            END IF
            write(6,*) 'Error: Index ',I,' not found in step level ',
     &           L,'!'
            call flush(6)
            call abort()
         END IF
      END IF
*     New step level
      NL = L + DL
*     See whether the shift is reasonable
      IF(NL.GT.64) then
         write(6,*) 'Error!: Too small Step: I',I,
     &        'Step level',NL
         call flush(6)
         call abort()
      END IF
*     See whether need to decrease NDTMIN or increase NDTMAX
      IF(NL.LT.NDTMIN.AND.NL.GE.1) NDTMIN = NL
      IF(NL.GT.NDTMAX) NDTMAX = NL
*     If exceed maximum step, treat as special particle like ghost and move index outside nxtlimit 
      IF(NL.LT.1) THEN
         NL = 1
         NXTLIMIT = NXTLIMIT - 1
      END IF
*     Shifting to larger step level (smaller time step)
      IF(DL.GT.0) THEN
*     Shift first index of current level to position J
         NXTLST(J) = NXTLST(NDTK(L+1)+1)
         DO LK = L+1, NL-1
*     Increase NDTK by one
            NDTK(LK) = NDTK(LK) + 1
*     Shift first index to the end position of level LK
            NXTLST(NDTK(LK)) = NXTLST(NDTK(LK+1)+1)
C            NXTK(NDTK(LK)) = LK
         END DO
*     Increase target level NDTK by one
         NDTK(NL) = NDTK(NL) + 1
*     Store index I in the last position of new level NL
         NXTLST(NDTK(NL)) = I
C         NXTK(NDTK(NL)) = NL
*     Shifting to smaller step level (larger time step)
      ELSE
*     Shift last index of current level to position J
         NXTLST(J) = NXTLST(NDTK(L))
         DO LK = L, NL+1, -1
*     Shift last index to the first position of level LK
            NXTLST(NDTK(LK)) = NXTLST(NDTK(LK-1))
C            NXTK(NDTK(LK)) = LK-1
*     Reduce NDTK(LK) by one
            NDTK(LK) = NDTK(LK) - 1
         END DO
*     Store index I in the last position of new level NL
         NXTLST(NDTK(NL)) = I
C         NXTK(NDTK(NL)) = NL
      END IF

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
