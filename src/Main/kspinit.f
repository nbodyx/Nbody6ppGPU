      subroutine kspinit
*     
*     
*     Generate perturber list for primordial binaries during initialization
*     ------------------------------------------------------

      include 'common6.h'
      REAL*8  UI(4),VI(4)
*     

      TIME0 = TIME
!$omp parallel do 
!$omp& private(ipair,icm,i1,i2,eb,semi,tk,vi,ui,tp,imod,K)
      DO IPAIR = 1, NPAIRS
         ICM = IPAIR + N
         I2 = 2*IPAIR
         I1 = I2 - 1

*     Form perturber list.
         CALL KSLIST(IPAIR)
*     
*     Transform any unperturbed hard binary to apocentre and set time-step.
         EB = H(IPAIR)*BODY(I1)*BODY(I2)/BODY(ICM)
         SEMI = -0.5*BODY(ICM)/H(IPAIR)
         IF (LIST(1,I1).EQ.0.AND.EB.LT.EBH) THEN
            TK = TWOPI*SEMI*SQRT(SEMI/BODY(ICM))
*     Note TIME is not commensurate after KSPERI (cf. CHTERM & STEPS).
            DO 55 K = 1,4
               UI(K) = U(K,IPAIR)
               VI(K) = UDOT(K,IPAIR)
 55         CONTINUE
*     Determine pericentre time (TP < 0 if TDOT2 < 0) and add TK/2.
            CALL TPERI(SEMI,UI,VI,BODY(ICM),TP)
*     Note: apocentre to apocentre gives almost zero step.
            STEP(I1) = 0.5*MIN(TK,STEP(ICM)) - TP
*     Transform KS variables to peri and by pi/2 to apocentre (skip apo).
            IF (ABS(TDOT2(IPAIR)).GT.1.0E-12.OR.R(IPAIR).LT.SEMI) THEN
               TIME = TIME0
               CALL KSPERI(IPAIR)
               CALL KSAPO(IPAIR)
*     Reset TIME to quantized value (small > 0 or < 0 possible initially).
            ELSE IF (TDOT2(IPAIR).GT.0.0) THEN
               TDOT2(IPAIR) = -1.0E-20
            END IF
         END IF
*     
*     Estimate an appropriate KS slow-down index for G < GMIN.
         IMOD = 1
         IF (LIST(1,I1).EQ.0.AND.SEMI.GT.0.0) THEN
            TK = TWOPI*SEMI*SQRT(SEMI/BODY(ICM))
            IF (KZ(26).GT.0.AND.STEP(ICM).GT.TK) THEN
               IMOD = 1 + INT(LOG(STEP(ICM)/TK)/0.69)
               IMOD = MIN(IMOD,5)
            END IF
         END IF

*     Obtain polynomials for perturbed KS motion (standard case & merger).
         TIME = TIME0
         CALL KSPOLY(IPAIR,IMOD)

         LIST(2,I2) = -1
*     Diagnostics
c         if(rank.eq.0) then
c            write(6,60) IPAIR, NAME(ICM), NAME(I1),NAME(I2),STEP(ICM),
c     &           STEP(I1),LIST(1,ICM),LIST(1,I1),GAMMA(IPAIR),IMOD
c 60         format('KS PAIR: IPAIR',I10,' NAME(ICM,I1,I2)',3I10,
c     &           ' STEP(ICM,I1)',1P,2E17.5,0P,' NNB',I5,' NP',I5,
c     &           ' GAMMA',1P,E17.5,0P,'  IMOD',I3)
c         end if
      END DO
!$omp end parallel do      
      TIME = TIME0

      RETURN

      END
