      SUBROUTINE UPDATE(IPAIR)
*
*
*       List modifications after KS termination.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      LOGICAL FUPDATE
*
*     Get index of ipair c.m. and components
      IPAIRCM = N + IPAIR
      I2_IPAIR = 2*IPAIR
      I1_IPAIR = I2_IPAIR - 1

*     Set largest index of any regularized member.
      ILAST = NPAIRS + 1
      ILASTCM = NTOT + 1
      I2_ILAST = 2*ILAST
      I1_ILAST = I2_ILAST - 1
*
*
*       Rename lists containing single regularized components.
      DO 10 J = 1,NTOT
*       Skip renaming if first neighbour exceeds regularized components.
***** Note:Christoph pointed out a problem where the neighbour list of a
*****     zero mass ghost particle grew too much due to KS terminations.
*****     This can be dealt with in several ways, the simplest being to
*****     skip the updating procedure in the two big loops of update.f
*****     where one c.m. particle is replaced by two KS components.
c$$$          IF (LIST(2,J).GT.ILAST) GO TO 50
c$$$          IF (LIST(2,J).GT.ILAST.OR.BODY(J).EQ.0.0D0) GO TO 50
         IF (BODY(J).EQ.0.0D0) GO TO 10
         NNB = LIST(1,J) + 1
         IF (NNB.EQ.1) GO TO 10

*     Indicator and index of ipair and last pair
         ID_IPAIR = 0
         ID_IPAIRCM = 0
         INUM_IPAIR = 0
         IPAIRLIMIT = 0
         IPCMLIMIT = 0
         ID_ILAST = 0
         ID_ILASTCM = 0
         ILASTLIMIT = 0
         INUM_ILAST = 0
         ILCMLIMIT = 0
         FUPDATE = .false.

         DO L = 2, NNB
            IF(IPAIRLIMIT.EQ.0.AND.LIST(L,J).GT.I2_IPAIR)
     &           IPAIRLIMIT = L
            IF(LIST(L,J).GT.I2_ILAST) THEN
               ILASTLIMIT = L
               GO TO 20
            END IF
            IF(LIST(L,J).EQ.I1_IPAIR) THEN
               ID_IPAIR = L
               INUM_IPAIR = 1
            END IF
            IF(LIST(L,J).EQ.I2_IPAIR) THEN
               IF(INUM_IPAIR.EQ.1) THEN
                  INUM_IPAIR = 2
               ELSE
                  ID_IPAIR = L
                  INUM_IPAIR = 1
               END IF
            END IF
            IF(LIST(L,J).EQ.I1_ILAST) THEN
               ID_ILAST = L
               INUM_ILAST = 1
            END IF
            IF(LIST(L,J).EQ.I2_ILAST) THEN
               IF(INUM_ILAST.EQ.1) THEN
                  INUM_ILAST = 2
               ELSE
                  ID_ILAST = L
                  INUM_ILAST = 1
               END IF
            END IF
         END DO

 20      DO L = NNB, 2, -1
            IF(ILCMLIMIT.EQ.0.AND.LIST(L,J).LT.ILASTCM)
     &           ILCMLIMIT = L
            IF(LIST(L,J).LT.IPAIRCM) THEN
               IPCMLIMIT = L
               GO TO 30
            END IF
            IF(LIST(L,J).EQ.ILASTCM) ID_ILASTCM = L
            IF(LIST(L,J).EQ.IPAIRCM) THEN
               ID_IPAIRCM = L
            END IF
         END DO

*     Check the case ipair member is in list
 30      IF(ID_IPAIR.GT.0) THEN
*     Safe check
            IF(ID_IPAIRCM.GT.0) THEN
               write(6,*) 'Error! List including both Ipair member',
     &              IPAIR,'p',ID_IPAIR,' and its c.m.',ILASTCM,'p',
     &              ID_IPAIRCM,', Name(CM)',
     &              NAME(ILASTCM),' TIME',TIME,'r',rank
               CALL FLUSH(6)
               call abort()
            end if
*     Case both ipair and ilast members are in list
            IF(ID_ILAST.GT.0) THEN
*     Safe check
               IF(ID_ILASTCM.GT.0) THEN
                  write(6,*) 'Error! List including both Last pair ',
     &                 'member',ILAST,'p',ID_ILAST,' and its c.m.'
     &                 ,ILASTCM,'p',ID_ILASTCM,
     &                 ', Name(CM)',NAME(ILASTCM),' TIME',TIME,'r',rank
                  CALL FLUSH(6)
                  CALL abort()
               END IF
*     If ipair is ilast or both ipair and ilast exist with same member
*     number in list, don't need to shift other data
               IF(INUM_IPAIR.EQ.INUM_ILAST.AND.INUM_IPAIR.eq.1) THEN
                  IF(LIST(ID_IPAIR,J).EQ.I1_IPAIR) THEN
                     IF(LIST(ID_ILAST,J).EQ.I2_ILAST) THEN
                        LIST(ID_IPAIR,J) = I2_IPAIR
                        LIST(ID_ILAST,J) = I1_ILAST
                        FUPDATE = .true.
                     END IF
                  ELSE IF(LIST(ID_IPAIR,J).EQ.I2_IPAIR) THEN
                     IF(LIST(ID_ILAST,J).EQ.I1_ILAST) THEN
                        LIST(ID_IPAIR,J) = I1_IPAIR
                        LIST(ID_ILAST,J) = I2_ILAST
                        FUPDATE = .true.                     
                     END IF
                  END IF
*     Case INUM_IPAIR = 1 & INUM_ILAST = 2
               ELSE IF(INUM_IPAIR.LT.INUM_ILAST) THEN
*     Make clear which member of IPAIR in the new position
                  IF(LIST(ID_IPAIR,J).EQ.I1_IPAIR) THEN
                     LIST(ID_ILAST+1,J) = I1_ILAST
                  ELSE
                     LIST(ID_ILAST+1,J) = I2_ILAST
                  END IF
*     shift all index from IPAIR + 1 to ILAST -1 to right by one
                  LIST(ID_IPAIR+2:ID_ILAST,J) =
     &                 LIST(ID_IPAIR+1:ID_ILAST-1,J)
*     Notice now I1/2_IPAIR is new index of I1/2_ILAST, reset ID_IPAIR(+0/1) value
                  LIST(ID_IPAIR,J) = I1_IPAIR
                  LIST(ID_IPAIR+1,J) = I2_IPAIR
                  FUPDATE = .true.
*     Case INUM_IPAIR = 2 & INUM_ILAST = 1
               ELSE IF(INUM_IPAIR.GT.INUM_ILAST) THEN
*     Make clear which member of ILAST
                  IF(LIST(ID_ILAST,J).EQ.I1_ILAST) THEN
                     LIST(ID_IPAIR,J) = I1_IPAIR
                  ELSE
                     LIST(ID_IPAIR,J) = I2_IPAIR
                  END IF
*     shift all index from ipair + 2 to ilast -1 to left by one
                  LIST(ID_IPAIR+1:ID_ILAST-2,J) =
     &                 LIST(ID_IPAIR+2:ID_ILAST-1,J)
*     Notice now I1/2_ILAST is new index of I1/2_IPAIR, reset ID_ILAST(+0/1) value
                  LIST(ID_ILAST-1,J) = I1_ILAST
                  LIST(ID_ILAST,J) = I2_ILAST
                  FUPDATE = .true.                  
               END IF
*     Case only ipair member is in the list
            ELSE
               IF(ILASTLIMIT.EQ.0) ILASTLIMIT = NNB + 1
               IF(INUM_IPAIR.EQ.1) THEN
                  ITMP = 0
                  IF(LIST(ID_IPAIR,J).EQ.I1_IPAIR) THEN
                     ITMP = I1_ILAST
                  ELSE
                     ITMP = I2_ILAST
                  END IF
*     shift all index from ipair+inum_ipair to ilastlimit to left by 1
                  LIST(ID_IPAIR:ILASTLIMIT-2,J) =
     &                 LIST(ID_IPAIR+1:ILASTLIMIT-1,J)
                  LIST(ILASTLIMIT-1,J) = ITMP
                  FUPDATE = .true.
*     shift all index from ipair+inum_ipair to ilastlimit to left by 2
               ELSE
                  LIST(ID_IPAIR:ILASTLIMIT-3,J) =
     &                 LIST(ID_IPAIR+2:ILASTLIMIT-1,J)
                  LIST(ILASTLIMIT-2,J) = I1_ILAST
                  LIST(ILASTLIMIT-1,J) = I2_ILAST
                  FUPDATE = .true.                  
               END IF
            END IF
*     Case only ilast member is in the list
         ELSE IF(ID_ILAST.GT.0) THEN
*     Safe check
            IF(ID_ILASTCM.GT.0) THEN
               write(6,*) 'Error! List including both Last pair ',
     &              'member',ILAST,'p',ID_ILAST,' and its c.m.'
     &              ,ILASTCM,'p',ID_ILASTCM,
     &              ', Name(CM)',NAME(ILASTCM),' TIME',TIME,'r',rank
               CALL FLUSH(6)
               CALL abort()
            END IF
*     IF(IPAIRLIMIT.EQ.0) IPAIRLIMIT = 2
            IF(INUM_ILAST.EQ.1) THEN
               ITMP = 0
               IF(LIST(ID_ILAST,J).EQ.I1_ILAST) THEN
                  ITMP = I1_IPAIR
               ELSE
                  ITMP = I2_IPAIR
               END IF
               LIST(IPAIRLIMIT+1:ID_ILAST,J) =
     &              LIST(IPAIRLIMIT:ID_ILAST-1,J)
               LIST(IPAIRLIMIT,J) = ITMP
               FUPDATE = .true.               
            ELSE
               LIST(IPAIRLIMIT+2:ID_ILAST+1,J) =
     &              LIST(IPAIRLIMIT:ID_ILAST-1,J)
               LIST(IPAIRLIMIT,J) = I1_IPAIR
               LIST(IPAIRLIMIT+1,J) = I2_IPAIR
               FUPDATE = .true.               
            END IF
         END IF
*     Check the ipair c.m. is in the list
         IF(ID_IPAIRCM.GT.0) THEN
*     Safe check
            IF(ID_IPAIR.GT.0) THEN
               write(6,*) 'Error! List including both Ipair member',
     &              IPAIR,'p',ID_IPAIR,' and its c.m.',ILASTCM,'p',
     &              ID_IPAIRCM,', Name(CM)',
     &              NAME(ILASTCM),' TIME',TIME,'r',rank
               CALL FLUSH(6)
               CALL abort()
            END IF
*     Case the ilast c.m. is also in the list
            IF(ID_ILASTCM.GT.0) THEN
               LIST(ILASTLIMIT+2:ID_ILASTCM+1,J) =
     &              LIST(ILASTLIMIT:ID_ILASTCM-1,J)
               LIST(ILASTLIMIT,J) = I1_ILAST
               LIST(ILASTLIMIT+1,J) = I2_ILAST
*     Increase neighbor list by 1
               LIST(1,J) = LIST(1,J) + 1
               IF(LIST(1,J).GE.NNBMAX) THEN
                  write(6,*),'Warning! J',J,'LIST',LIST(1,J),
     &                 ' GE NNBMAX',NNBMAX,'T',TIME
                  call flush(6)
               END IF
               FUPDATE = .true.               
*     Case ipair c.m. and ilast member is in the list
            ELSE
*     Shift index from ID_IPAIRCM+2 to NNB to right by 1
               LIST(ID_IPAIRCM+2:NNB+1,J) =
     &              LIST(ID_IPAIRCM+1:NNB,J)
*     Shift index from ILASTLIMIT+2 to IPAIRCM-1 to right by 2
               LIST(ILASTLIMIT+2:ID_IPAIRCM+1,J) =
     &              LIST(ILASTLIMIT:ID_IPAIRCM-1,J)
*     Reset ILASTLIMIT(+0/1) to new IPAIR index (Original ILAST index)
               LIST(ILASTLIMIT,J) = I1_ILAST
               LIST(ILASTLIMIT+1,J) = I2_ILAST
               LIST(1,J) = LIST(1,J) + 1
               IF(LIST(1,J).GE.NNBMAX) THEN
                  write(6,*),'Warning! J',J,'LIST',LIST(1,J),
     &                 ' GE NNBMAX',NNBMAX,'T',TIME
                  call flush(6)
               END IF
               FUPDATE = .true.               
            END IF   
*     Case only ilast c.m. is in the list
         ELSE IF(ID_ILASTCM.GT.0) THEN
*     right shift by 1 from LPCMLIMIT
            IF(IPCMLIMIT.EQ.0) IPCMLIMIT = 1
            LIST(IPCMLIMIT+2:ID_ILASTCM,J) =
     &           LIST(IPCMLIMIT+1:ID_ILASTCM-1,J)
*     Notice the ilast c.m. now has new position of IPAIRCM
            LIST(IPCMLIMIT+1,J) = IPAIRCM
            FUPDATE = .true.
         END IF

#ifdef SIMD
*     Update neighbor list in AVX/SSE library
         IF(J.GE.IFIRST.AND.FUPDATE.AND.TIME.NE.0.0D0) THEN
            CALL IRR_SIMD_SET_LIST(J,LIST(1,J))
         END IF
#endif         
*     --04/15/14 15:37-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         IF(FUPDATE) THEN
c$$$            NNB = LIST(1,J)
c$$$            write(120+rank,*)'IPAIR',IPAIR,'IPAIRCM',IPAIRCM,
c$$$     &           'ILAST',ILAST,'ILASTCM',ILASTCM,
c$$$     &           'NAME(IL)',NAME(I1_ILAST),NAME(I2_ILAST),
c$$$     &           'NAME(IP)',NAME(I1_IPAIR),NAME(I2_IPAIR),
c$$$     &           'NAME(IPCM)',NAME(IPAIRCM),'NAME(ILCM)',
c$$$     &           NAME(ILASTCM),'NTOT',NTOT
c$$$            write(120+rank,*)'NEW',LIST(1:NNB+1,J),FUPDATE
c$$$            write(120+rank,*)'INUM_ILAST',INUM_ILAST,'INUM_IPAIR',
c$$$     &           INUM_IPAIR,'ID_ILAST',ID_ILAST,'ID_IPAIR',ID_IPAIR,
c$$$     &           'ID_ILASTCM',ID_ILASTCM,'ID_IPAIRCM',ID_IPAIRCM
c$$$            write(120+rank,*)'LIMIT: IPAIR',IPAIRLIMIT,'ILAST',
c$$$     &           ILASTLIMIT,'IPCM',IPCMLIMIT,'ILCM',ILCMLIMIT,'T',TIME
c$$$            call flush(120+rank)
c$$$         END IF
*     --04/15/14 15:37-lwang-end----------------------------------------*
 10   CONTINUE
*         
*       Modify the list of previously regularized binaries.
      NNB = LISTR(1) - 1
      L = 0
   91 L = L + 2
   92 IF (L.GT.NNB + 1) GO TO 96
      J = LISTR(L)
      K = LISTR(L+1)
*       First check the current two-body separation of any old pairs.
      RJK2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 +
     &                              (X(3,J) - X(3,K))**2
*       Remove pair if RJK > 4*RMIN when special procedure not needed.
      IF (RJK2.LT.16.0*RMIN**2) GO TO 91
      DO 94 K = L,NNB
          LISTR(K) = LISTR(K+2)
   94 CONTINUE
      NNB = NNB - 2
      GO TO 92
*
*       Add ICOMP & JCOMP to LISTR (maximum of MLR/2 - 1 pairs).
   96 IF (NNB.GT.MLR - 4) THEN
*       Note that NNB is one less than the actual membership.
          DO 98 K = 2,NNB
              LISTR(K) = LISTR(K+2)
   98     CONTINUE
          NNB = NNB - 2
*       Removal of the oldest KS pair.
      END IF
      LISTR(NNB+3) = ICOMP
      LISTR(NNB+4) = JCOMP
      LISTR(1) = NNB + 3
*
*       Copy flag index of disrupted pair (set in KSTERM).
      IFLAG = JLIST(3)
*       Add primordial pairs to LISTD (skip new KS pairs or primordials).
      IF (IFLAG.EQ.0.OR.IABS(JLIST(1) - JLIST(2)).EQ.1) GO TO 110
*
*       Check list of disrupted component names.
      NNB = LISTD(1) - 1
      KCOMP = 0
      DO 100 K = 2,NNB+1,2
          IF (LISTD(K).EQ.JLIST(1).AND.LISTD(K+1).EQ.JLIST(2)) KCOMP = 1
  100 CONTINUE
*
*       Include both components unless already members.
      IF (KCOMP.EQ.0) THEN
          IF (NNB.GT.MLD - 4) THEN
              DO 102 K = 2,NNB
                 LISTD(K) = LISTD(K+2)
  102         CONTINUE
              NNB = NNB - 2
          END IF
*       Add most recent names at the end (limit is MLD/2 - 1 pairs).
          LISTD(NNB+3) = JLIST(1)
          LISTD(NNB+4) = JLIST(2)
          LISTD(1) = NNB + 3
      END IF
      IF (rank.eq.0.and.IFLAG.NE.-1) THEN
          WRITE (6,104)  IPAIR, IFLAG, JLIST(1), JLIST(2)
  104     FORMAT (' LISTD INCONSISTENCY!!  IPAIR IFLAG NAMES ',2I5,2I8)
      END IF
*
*       Update list of high velocity particles containing c.m. members.
 110  L = 2
 125  IF (LISTV(L).EQ.IPAIRCM) THEN
*     Remove old c.m.
         LISTV(L) = LISTV(NNB+1)
         LISTV(1) = LISTV(1) - 1
         GO TO 130
      END IF
*     Replace ilast c.m. to new position (ipair c.m.)
 130  IF(LISTV(L).EQ.ILASTCM) LISTV(L) = IPAIRCM
      IF(L.LT.LISTV(1)+1) THEN
         L = L + 1
         GO TO 125
      END IF
*
      RETURN
*
      END
