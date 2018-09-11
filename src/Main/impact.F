      SUBROUTINE IMPACT(I)
*
*
*       Multiple collision or merger search.
*       ------------------------------------
*
      Include 'kspars.h'
      INCLUDE 'common6.h'
*     Unsafe for parallel!!, modified value: tmdis 
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
*     safe for parallel, no value change, used value: names
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      CHARACTER*8  WHICH1
      REAL*8  XX(3,3),VV(3,3),MFAC
      INTEGER LISTQ(100)
      LOGICAL flag_delay,MCF
*     For diagnostic, safe for parallel
      SAVE LISTQ,QCHECK
      DATA IZARE,LISTQ(1),QCHECK /0,0,0.0D0/
*
*
*     Initial flag_delay to .false.
      flag_delay = .false.

*     Indicate whether merger check is already done
      MCF =.false.
      
*       Set index of KS pair & both components of c.m. body #I.
      IPAIR = I - N
      call jpred(I,TIME,TIME)
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      NTTRY = NTTRY + 1
      PERT1 = 0.0
      PERT2 = 0.0
      JCOMP = IFIRST
      NP = 0
      KS2 = 0
      RMAX2 = 1.0
      TTOT = TIME + TOFF
      RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                               (X(3,I) - RDENS(3))**2
*
*       Search c.m. neighbours if binary has at most two perturbers.
      J1 = I1
      IF (LIST(1,J1).LE.2) J1 = I
      NNB2 = LIST(1,J1) + 1
*     
*       Set semi-major axis, eccentricity & apocentre of inner binary.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      A0 = SEMI
      ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
      APO = ABS(SEMI)*(1.0 + ECC)

*     Find the dominant body (JCOMP) and nearest perturber (JMAX).
      RCRIT3 = 10.0*APO**3/BODY(I)/GMIN
      RCRIT2 = 10.0D4*RMIN*RMIN
      DO 10 L = 2,NNB2
          J = LIST(L,J1)
          call jpred(J,TIME,TIME)
          RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                  (X(3,I) - X(3,J))**2
          RIJ = SQRT(RIJ2)
          RIJ3 = RIJ2*RIJ
          IF(J.GT.N) THEN
             JPAIR = J-N
             RIJ = MAX(0.0,RIJ-R(JPAIR))
             RIJ3 = RIJ*RIJ*RIJ
          END IF
*     J1 = I is neighbor list
          IF (J1.EQ.I.AND.
     &         (RIJ3.GT.BODY(J)*RCRIT3.AND.RIJ2.GT.RCRIT2)) GO TO 10
          NP = NP + 1
          JLIST(NP) = J
          PERT = BODY(J)/(RIJ2*SQRT(RIJ2))
          IF (PERT.GT.PERT2) THEN 
              IF (PERT.GT.PERT1) THEN
                  RJMIN2 = RIJ2
                  JMAX = JCOMP
                  JCOMP = J
                  PERT2 = PERT1
                  PERT1 = PERT
              ELSE
                  JMAX = J
                  PERT2 = PERT
                  RMAX2 = RIJ2
              END IF
          END IF
   10 CONTINUE
*
*       Include safety check on rare case of no perturbers inside 100*RMIN.
C      IF (NP.LT.2.AND.NNB2.GE.3) THEN
C         RPCRIT2 = 100*RPCRIT2
C         GO TO 5
C      END IF
      IF (NP.EQ.0) GO TO 100
*
#ifdef DECISION_CHECK
      write(6,*) 'IMPACT: I=',I, ' N=',N,' NAME=',NAME(I),
     &     ' NP=',NP,' JCOMP=',JCOMP,' NAME(JCOMP)=',NAME(JCOMP),
     &     ' JMAX=',JMAX,' NAME(JMAX)=',NAME(JMAX),
     &     ' NCH=',NCH,' Time=',TIME,' NNB=',NNB2
#endif

      CALL JPRED(JCOMP,TIME,TIME)
      RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
*
*     Mass ratio factor
      MFAC = MAX(sqrt(BODY(JCOMP)/BODY(I)),1.0)
      
*       Specify larger perturbation for optional chain regularization.
      IF ((KZ(30).GT.0.OR.KZ(30).EQ.-1).AND.NCH.EQ.0) THEN
          GSTAR = 100.0*GMIN
          KCHAIN = 1
      ELSE
          GSTAR = GMIN
          KCHAIN = 0
C*       Specify indicator -1 for allowing TRIPLE & QUAD but not CHAIN.
C          IF (KZ(30).EQ.-2) KCHAIN = -1
      END IF
*
*       Only accept inward motion or small secondary perturbation.
      PERT3 = 2.0*R(IPAIR)**3*PERT2/BODY(I)
*     --03/17/14 21:11-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(i.eq.5099.and.tblock.ge.4.55627441406250000E-002) then
c$$$      print*,rank,'impact 101 RDOT',RDOT,'PERT3',PERT3,'GSTAR',GSTAR,'t'
c$$$     *     ,time
c$$$      call flush(6)
c$$$      end if
*     --03/17/14 21:11-lwang-end----------------------------------------*
#ifdef DECISION_CHECK
      IF (RDOT.GT.0.0.OR.PERT3.GT.100.0*GSTAR)
     &     write(6,*) 'IMPACT: Outward motion ',
     &     'or secondary perturbation too large ',
     &     'RDOT(',RDOT,')>0 || PERT3(',
     &     PERT3,')>100 GSTAR(',GSTAR,') L.117 quit'
#endif
      IF (RDOT.GT.0.0.OR.PERT3.GT.100.0*GSTAR) GO TO 100
*
*       Include impact parameter test to distinguish different cases.
      A2 = (XDOT(1,I) - XDOT(1,JCOMP))**2 + 
     &     (XDOT(2,I) - XDOT(2,JCOMP))**2 +
     &     (XDOT(3,I) - XDOT(3,JCOMP))**2
      RIJ = SQRT(RJMIN2)
      A3 = 2.0/RIJ - A2/(BODY(I) + BODY(JCOMP))
      SEMI1 = 1.0/A3
      A4 = RDOT**2/(SEMI1*(BODY(I) + BODY(JCOMP)))
      ECC1 = SQRT((1.0D0 - RIJ/SEMI1)**2 + A4)
      PMIN = SEMI1*(1.0D0 - ECC1)
*
*
*       Quit on hyperbolic orbit with large impact parameter.
#ifdef DECISION_CHECK
      IF (ECC1.GT.1.0.AND.PMIN.GT.50.0*SEMI)
     &     write(6,*) 'IMPACT: hyperbolic orbit ',
     &     'and large impact parameter ',
     &     'ECC1(',ECC1,'>1.0 && PMIN(',
     &     PMIN,')>50 SEMI(',SEMI,') L.143 quit'
#endif
      IF (ECC1.GT.1.0.AND.PMIN.GT.50.0*SEMI) GO TO 100
*
*       Include rectification for non-circular binaries with KSTAR = 10 & 12.
      IF (KZ(27).EQ.2.AND.KSTAR(I).GE.10.AND.KSTAR(I).LE.18) THEN
          IF (ECC.GT.0.1.AND.MOD(KSTAR(I),2).EQ.0) THEN
              RM = MAX(RADIUS(I1),RADIUS(I2))
              ICIRC = -1
              CALL INDUCE(IPAIR,JCOMP,EMAX,EMIN,ICIRC,TC,ANGLE,TG,EDAV)
              if(rank.eq.0)
     &        WRITE (6,15)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &                      KSTAR(I), LIST(1,I1), ECC, SEMI, RM, PMIN,
     &                      GAMMA(IPAIR), TC, 360.0*ANGLE/TWOPI, EMAX
   15         FORMAT (' NON-CIRCULAR    NM K* NP E A R* PM G TC IN EX ',
     &                                  2I7,4I4,F7.3,1P,5E9.1,0P,2F7.2)
*       Circularize the orbit instantaneously for short TC.
              IF (TC.LT.-100.0) THEN
*       Set temporary unperturbed orbit (avoids new KSPOLY in DEFORM).
                  NP = LIST(1,I1)
                  LIST(1,I1) = 0
                  DT1 = STEP(I1)
                  TIME0 = TIME
                  CALL KSRECT(IPAIR)
                  QP = SEMI*(1.0 - ECC)
                  ERR = ABS(QP - R(IPAIR))/QP
*       Deform orbit to circular eccentricity after determining apocentre.
                  IF (R(IPAIR).GT.SEMI.OR.ERR.GT.1.0D-04) THEN
*       Reduce eccentric anomaly by pi for inward motion.
                      IF (TDOT2(IPAIR).LT.0.0D0) THEN
                          CALL KSAPO(IPAIR)
                      END IF
                  END IF
*       Predict to pericentre and transform by pi to exact apocentre.
                  CALL KSPERI(IPAIR)
                  CALL KSAPO(IPAIR)
                  ECCM = 0.002
                  CALL DEFORM(IPAIR,ECC,ECCM)
                  LIST(1,I1) = NP
*        Resolv X & XDOT and initialize KS polynomial at apocentre time.
                  IF (NP.EQ.0) DT1 = 0.0
                  TIME = TIME0 - DT1
                  CALL RESOLV(IPAIR,1)
                  CALL KSPOLY(IPAIR,1)
              ELSE
                  KSTAR(I) = 0
*     ks MPI communication
                  call ksparmpi(K_store,K_int,K_KSTAR,I,0,KSTAR(I))
              END IF
          END IF
      END IF
*
*       Form binding energy of inner & outer binary.
      EB = BODY(I1)*BODY(I2)*H(IPAIR)/BODY(I)
      IF(ABS(EB).LT.1.0D-10) EB = -1.0D-10
      EB1 = -0.5*BODY(JCOMP)*BODY(I)/SEMI1
*
*       Obtain the total perturbing force acting on body #I & JCOMP.
      CALL FPERT(I,JCOMP,NP,PERT)
*
*       Choose maximum of dominant scalar & total vectorial perturbation.
      PERT = PERT*RJMIN2/(BODY(I) + BODY(JCOMP))
      PERT4 = 2.0*RJMIN2*RIJ*PERT2/(BODY(I) + BODY(JCOMP))
      PERTM = MAX(PERT4,PERT)
*
*       Estimate relative perturbation at apocentre from actual value.
      RA = SEMI1*(1.0 + ECC1)
      IF (SEMI1.LT.0.0) RA = RIJ
      GI = PERT*(RA/RIJ)**3
*       Copy coordinates and velocities to local variables.
      DO 42 K = 1,3
          XX(K,1) = X(K,I1)
          XX(K,2) = X(K,I2)
          XX(K,3) = X(K,JCOMP)
          VV(K,1) = XDOT(K,I1)
          VV(K,2) = XDOT(K,I2)
          VV(K,3) = XDOT(K,JCOMP)
  42  CONTINUE
*
*     Delay for large distance ratio outside RMIN with strong perturbation.
#ifdef DECISION_CHECK
      IF ((RIJ.GT.MFAC*RMIN.AND.PERT.GT.GMAX).OR.RIJ.GT.5*MFAC*RMIN)
     &     write(6,*) 'IMPACT: RIJ(',RIJ,')>MFAC(',MFAC,
     &     ')*RMIN(',RMIN,') &&',
     &     'PERT(',PERT,')>GMAX(',GMAX,') || RIJ>5 RMIN L.247 quit'
#endif
      IF ((RIJ.GT.MFAC*RMIN.AND.PERT.GT.GMAX).OR.RIJ.GT.5*MFAC*RMIN)
     &     GO TO 100

*     Determine the inclination (in radians).
      CALL INCLIN(XX,VV,X(1,I),XDOT(1,I),ANGLE)
      
*       Use combined semi-major axis for binary-binary collision.
      IF (JCOMP.GT.N) THEN
          JPAIR = JCOMP - N
          SEMI2 = -0.5D0*BODY(JCOMP)/H(JPAIR)
*     Get secondary binary eccentricity
          ECC2 = (1.0D0 - R(JPAIR)/SEMI2)**2 +
     &         TDOT2(JPAIR)**2/(BODY(JCOMP)*SEMI2)
          ECC_2 = SQRT(ECC2)
          J1 = 2*JPAIR - 1
          EB2 = -0.5*BODY(J1)*BODY(J1+1)/SEMI2
*       Define SEMI0 as smallest binary in case IPAIR denotes widest pair.
          SEMI0 = MIN(ABS(SEMI),ABS(SEMI2))
          SEMIX = MAX(SEMI,SEMI2)
          ECCX = MAX(ECC,ECC_2)
          APO = APO + MAX(ABS(SEMI2),R(JPAIR))
          SEMIS = SEMI + SEMI2
*       Do not allow negative or soft cross section
#ifdef DECISION_CHECK
          IF (SEMIS.GT.2.0*RMIN)
     &         write(6,*) 'IMPACT: B-B SEMIS(',SEMIS,')>2.0*RMIN(',
     &         RMIN,') L.233 quit'
#endif
          IF (SEMIS.GT.2.0*RMIN) GO TO 100
*       Consider merger for PMIN > SEMI and large semi-major axis ratio.
#ifdef DECISION_CHECK
          IF (PMIN.GT.SEMIS.AND.SEMI2.GT.20.0*SEMI0)
     &         write(6,*) 'IMPACT: B-B PMIN(',PMIN,')>SEMIS(',SEMIS,
     &         ') && SEMI2(',SEMI2,')>20.0 SEMI0 (',SEMI0,') L 244'
#endif          
          IF (PMIN.GT.SEMIS.AND.SEMI2.GT.20.0*SEMI0) GO TO 30
*      END IF

*     try stability check with maximum ecc and semi of inners
          QST = QSTAB(ECCX,ECC1,ANGLE,BODY(I1),BODY(I2),BODY(JCOMP))
          RP = PMIN/SEMIX
#ifdef DECISION_CHECK
          IF (QST.LT.RP) write(6,*) 'IMPACT: B-B QST(',QST,')<RP(',RP,
     &         ') && GI(',GI,')<=0.02 L416'
#endif
*     The GI should be small and consistent with the later check
          IF (QST.LT.RP.AND.GI.LE.0.02.AND.PERT.LE.GMAX) GO TO 30
*     Consider merger for PMIN > SEMI and large semi-major axis ratio.
*          IF (PMIN.GT.SEMI.AND.SEMI2.GT.20.0*SEMI0
*     &         .AND.1.0/SEMI.GE.0.5/RMIN) GO TO 30
*       Do not allow negative or soft cross section
*       Add PMIN > 2*RMIN to avoid chain with large PMIN
*          IF (1.0/SEMI.LT.0.5/RMIN.AND.PMIN.GT.2*RMIN) GO TO 100
#ifdef DECISION_CHECK
          IF (PMIN.GT.4.0*SEMIX.AND.
     &         (ECC1.GT.0.9.AND.ECC1.LT.1.0))
     &         write(6,*) 'IMPACT: L348'
#endif      
          IF ((PMIN.GT.4.0*SEMIX.AND.ECC1.GT.0.9.AND.ECC1.LT.1.0)
     &         .AND..not.MCF) GO TO 30
          IF (PMIN.GT.0.1*RMIN) THEN
#ifdef DECISION_CHECK
             IF (PMIN.GT.A0 + SEMI2) 
     &            write(6,*) 'IMPACT: L340'
#endif
             IF ((PMIN.GT.A0 + SEMI2).AND..not.MCF) GO TO 30
          END IF
      END IF

C*     Check separation in case of chain regularization.
C 16   IF (KCHAIN.GT.0) THEN
*       Form effective gravitational radius (combine triple & quad).
 16      EBT = EB + EB1
          ZMM = BODY(I1)*BODY(I2) + BODY(I)*BODY(JCOMP)
*       Set length of chain for decision-making (also used at termination).
          RSUM = R(IPAIR) + RIJ
          RI = R(IPAIR)
          IF (JCOMP.GT.N) THEN
              EBT = EBT + EB2
              ZMM = ZMM + BODY(J1)*BODY(J1+1)
              RSUM = RSUM + R(JPAIR)
              RI = MAX(R(JPAIR),RI)
          END IF
          RGRAV = ZMM/ABS(EBT)
*       Employ initial size as upper limit in case of weakly bound system.
          RGRAV = MIN(RGRAV,RMIN)
*       Save initial energy in binaries for routine SETSYS.
          EBCH0 = EBT - EB1
*     ks MPI communication EBCH0
          call ksparmpi(K_store,K_real8,K_EBCH0,0,0,EBCH0)
*     Use RIJ instead of RSUM in 3*RGRAV test (increases initial RIJ).
*     When there are two wide binaries with big mass ratio, add PMIN GT
*     2*RMIN, switch on Chain
*         IF ((RIJ.GT.MAX(3.0*RGRAV,RMIN).OR.RSUM.GT.2.0*RMIN).AND.
*      &         PMIN.GT.2.0*RMIN) GO TO 30
C          GI = 2.0*BODY(JCOMP)*(RI/RIJ)**3/BODY(I)
#ifdef DECISION_CHECK
          IF (RIJ.GT.MAX(3.0*RGRAV,RMIN).OR.RSUM.GT.2.0*RMIN)
     &         write(6,*) 'IMPACT: RIJ(',RIJ,')>MAX(3.0*RGRAV(',RGRAV,
     &         '),RMIN(',RMIN,') || RSUM(',RSUM,')>2.0 RMIN(',RMIN,
     &         ') && GI(',GI,')<0.05 L.286'
#endif
          IF (((RIJ.GT.MAX(3.0*RGRAV,RMIN).OR.RSUM.GT.2.0*RMIN))
     &         .AND.GI.LE.0.05.AND..not.MCF) GO TO 30
*       Enforce KS orbit using MERGE for high eccentricity if PMIN > 10*RI.
#ifdef DECISION_CHECK
          IF (ECC1.GT.0.99.AND.PMIN.GT.10.0*RI.AND.
     &        PERTM.LT.GMAX)
     &         write(6,*) 'IMPACT: ECC1(',ECC1,')>0.99 && PMIN(',
     &         PMIN,')>10 RI(',RI,') && PERTM(',PERTM,')>GMAX(',
     &         GMAX,') L 295'
#endif
          IF ((ECC1.GT.0.99.AND.PMIN.GT.10.0*RI.AND.
     &        PERTM.LT.GMAX).AND..not.MCF) GO TO 40
          IF (KZ(27).GT.0.AND.JCOMP.GT.N) THEN
              IF (SEMI0.LT.SEMI2) J1 = I1
              RT = 4.0*MAX(RADIUS(J1),RADIUS(J1+1))
*       Do not allow large distance ratio for nearly synchronous binary.
#ifdef DECISION_CHECK
              IF (SEMI0.GT.RT.AND.RI.GT.25.0*SEMI0)
     &             write(6,*) 'IMPACT: SEMI0(',SEMI0,')> RT(',RT,
     &             ') && RI(',RI,') > 25 SEMI0 L 300'
#endif              
              IF ((SEMI0.GT.RT.AND.RI.GT.25.0*SEMI0)
     &             .AND..not.MCF) GO TO 30
              IF (MIN(SEMI0,SEMI2).LT.0.05*RIJ) THEN
#ifdef DECISION_CHECK
                 write(6,*) 'IMPACT: MAX(SEMI0(',SEMI0,'),SEMI2(',
     &                SEMI2,')<0.1 RIJ(',RIJ,') L 306'
#endif
                 IF ((MAX(SEMI0,SEMI2).LT.0.1*RIJ) 
     &                .AND..not.MCF) GO TO 30                 
              END IF
          END IF
C      END IF
          
*     Avoid perturbed case when chain cannot be used
#ifdef DECISION_CHECK
      IF (KCHAIN.EQ.0.AND.GI.GT.1e-6) THEN
         write(6,*) 'KCHAIN.EQ.0.AND.GI(',GI,').GT.1e-6 L387 quit'
      END IF
#endif
      IF (KCHAIN.EQ.0.AND.GI.GT.1e-6) GO TO 100
*
*       Include special case of strong interraction and large ECC1.
      IF (ECC1.GT.0.9.AND.GAMMA(IPAIR).GT.0.01) THEN
          IF (APO.LT.0.01*RMIN.AND.PMIN.LT.2.5*APO) GO TO 17
      END IF
*
*       Adopt triple, quad or chain regularization for strong interactions.
*     IF ((APO.GT.0.01*RMIN.OR.JCOMP.GT.N).AND.PMIN.GT.APO) GO TO 30
****** Not good for a strong perturbed case     
C#ifdef DECISION_CHECK
C      IF ((APO.GT.0.01*RMIN.OR.JCOMP.GT.N).AND.PMIN.GT.1.5*APO)
C     &     write(6,*) 'IMPACT: L322'
C#endif      
C      IF ((APO.GT.0.01*RMIN.OR.JCOMP.GT.N).AND.PMIN.GT.1.5*APO) GO TO 30
*     IF (APO.LE.0.01*RMIN.AND.PMIN.GT.2.0*APO) GO TO 30

*     Remove case with large RIJ
#ifdef DECISION_CHECK
      IF (((RIJ.GT.MFAC*RMIN.AND.SEMI1.GT.0.0).OR.RIJ.GT.2*MFAC*RMIN)
     &     .AND.PERT.GT.GMAX)
     &     write(6,*) 'IMPACT: RIJ(',RIJ,')> MFAC(',MFAC,
     &     ')*RMIN(',RMIN,') && SEMI1(',SEMI1,
     &     ')>0 || RIJ>2.0 MFAC(',MFAC,')* RMIN && PERT(',PERT,
     &     ')>GMAX(',GMAX,') L.303 quit'
#endif
      IF (((RIJ.GT.MFAC*RMIN.AND.SEMI1.GT.0.0).OR.RIJ.GT.2*MFAC*RMIN)
     &     .AND.PERT.GT.GMAX)
     &     GO TO 100

C#ifdef DECISION_CHECK
C      IF (PERTM.GT.100.0*GSTAR)
C     &     write(6,*) 'IMPACT: PERTM(',PERTM,')>100*GSTAR(',GSTAR,
C     & ') L334'
C#endif
C      IF ((PERTM.GT.100.0*GSTAR).AND..not.MCF) GO TO 30
*
*       Check almost stable triples (factor 1.2 is experimental).
 17   IF (JCOMP.LE.N.AND.PMIN.GT.2.5*SEMI) THEN
          CALL HISTAB(IPAIR,JCOMP,PMIN,RSTAB)
C          RA = SEMI1*(1.0 + ECC1)
C          IF (SEMI1.LT.0.0) RA = RIJ
C          GI = PERT*(RA/RIJ)**3
*     using stablility criterion
          QST = QSTAB(ECC,ECC1,ANGLE,BODY(I1),BODY(I2),BODY(JCOMP))
          RP = PMIN/SEMI
#ifdef DECISION_CHECK
          IF (QST.LT.RP) write(6,*) 'IMPACT: QST(',QST,')<RP(',RP,
     &         ') && GI(',GI,')<=0.02 L416'
#endif
*     The GI should be small and consistent with the later check
          IF (QST.LT.RP.AND.GI.LE.0.02.AND.PERT.LE.GMAX.AND.
     &         .not.MCF) GO TO 40
*       Use estimated apocentre perturbation for decision-making.
C          IF (PMIN.GT.1.2*RSTAB) THEN
C#ifdef DECISION_CHECK
C              IF (GI.LT.0.05) write(6,*) 'IMPACT: GI<0.05 L361'
C#endif 
C              IF (GI.LT.0.05.AND..not.MCF) GO TO 30
C*       Choose chain for critical case of highly eccentric outer orbit.
C#ifdef DECISION_CHECK
C              IF (ECC1.LT.0.95) write(6,*) 'IMPACT: ECC1(',
C     &             ECC1,')<0.95 L.324'
C#endif
C              IF (ECC1.LT.0.95) GO TO 100
C          ELSE IF (PMIN.GT.0.7*RSTAB) THEN
CC*       Treat marginally stable triple according to external perturbation.
CC#ifdef DECISION_CHECK
CC              IF (GI.LT.0.05) write(6,*) 'IMPACT: GI<0.05 L373'
CC#endif 
CC              IF (GI.LT.0.05.AND..not.MCF) GO TO 30
C#ifdef DECISION_CHECK
C              IF (GI.LT.1.0.OR.ECC1.LT.0.9)
C     &             write(6,*), 'IMPACT: GI(',GI,')<1.0 || ECC1(',
C     &             ECC1,')<0.9 L.334'
C#endif              
C              IF (GI.LT.1.0.OR.ECC1.LT.0.9) GO TO 100
C          END IF
C#ifdef DECISION_CHECK
C          IF (PMIN.GT.0.6*RSTAB.AND.PMIN.LT.0.9*RSTAB)
C     &         write(6,*), 'PMIN(',PMIN,')>0.6 RSTAB(',RSTAB,
C     &         ') && PMIN<0.9 RSTAB L.341'
C#endif          
C          IF (PMIN.GT.0.6*RSTAB.AND.PMIN.LT.0.9*RSTAB) GO TO 100
C 
*     Delay for large distance ratio outside 0.5*RMIN.
#ifdef DECISION_CHECK
          IF (RIJ.GT.MAX(10.0*APO,0.5*MFAC*RMIN).AND.PERT.GT.GMAX)
     &         write(6,*) 'IMPACT: RIJ(',RIJ,')>MAX(10 APO(',
     &         APO,'), 2.0 MFAC(',MFAC,
     &         ')*RMIN(',RMIN,')) && PERT(',PERT,
     &         ')> GMAX(',GMAX,') L.348 quit'
#endif
          IF (RIJ.GT.MAX(10.0*APO,0.5*MFAC*RMIN).AND.PERT.GT.GMAX)
     &         GO TO 100
C#ifdef DECISION_CHECK
C          IF (RIJ.GT.10.0*APO)
C     &         write(6,*) 'IMPACT: RIJ(',
C     &         RIJ,')>10 APO(',APO,') L.353'
C#endif
C          IF (RIJ.GT.10.0*APO) GO TO 100
C#ifdef DECISION_CHECK
C          IF (PMIN.GT.2.5*APO) 
C     &         write(6,*) 'IMPACT: PMIN(',PMIN,')>2.5 APO(',
C     &         APO,') L 401' 
C#endif             
C** Suppress this condition since it do not give proper criterion for stable system  (L.Wang)
*     Instead, qstab is used
C          IF (PMIN.GT.2.5*APO) GO TO 40
      END IF
C#ifdef DECISION_CHECK
C      IF (PMIN.GT.3.0*SEMI.AND.JCOMP.LE.N)
C     &     write(6,*) 'IMPACT: L 408'
C#endif
C      IF (PMIN.GT.3.0*SEMI.AND.JCOMP.LE.N) GO TO 40
*
C      IF (JCOMP.GT.N) THEN
C#ifdef DECISION_CHECK
C          IF (RIJ.GT.10.0*APO)
C     &         write(6,*) 'IMPACT: JCOMP>N && RIJ(',
C     &         RIJ,')>10 APO(',APO,') L.365'
C#endif
C          IF (RIJ.GT.10.0*APO) GO TO 100
C      END IF
*       Skip chain if merged binary or chain c.m. (denoted by NAME <= 0).
#ifdef DECISION_CHECK
      IF (NAME(I).LT.0)
     &     write(6,*) 'IMPACT: NAME(I)(',NAME(I),
     &     ')<0 quit'
#endif
      IF (NAME(I).LT.0) GO TO 100

*     Terminate merger if exist
#ifdef DECISION_CHECK
      IF (NAME(JCOMP).LT.0)
     &     write(6,*) 'IMPACT: NAME(JCOMP)(',NAME(JCOMP),
     &     ')<0 term merge quit'
#endif
      IF (NAME(JCOMP).LT.0) then
         JPAIR = JCOMP - N
         J1 = 2*JPAIR - 1
         JCOMP = I
         KSPAIR = JPAIR
         IPHASE= 9
         call DELAY(0,0)
         call ksparmpi(K_store,K_int,K_DELAY,0,0,0.0)
         flag_delay=.true.
         GO TO 100
      endif

*     quit in the case of chain
      IF (NAME(I).EQ.0.OR.NAME(JCOMP).EQ.0) GO TO 100
*
*       Compare with existing subsystem of same type (if any).
      IF (NSUB.GT.0.AND.KCHAIN.LE.0) THEN
          PERIM = R(IPAIR) + RIJ
          IF (JCOMP.GT.N) PERIM = PERIM + R(JPAIR)
          IGO = 0
          CALL PERMIT(PERIM,IGO)
#ifdef DECISION_CHECK
          IF (IGO.GT.0)
     &         write(6,*) 'IMPACT: NSUB(',NSUB,')>0 && KCHAIN(',
     &         KCHAIN,')<0 && IGO(',IGO,')>0 L.389 quit'
#endif
          IF (IGO.GT.0) GO TO 100
      END IF
*
      
*       Skip all multiple regs on zero option (mergers done by #15 > 0).
#ifdef DECISION_CHECK
      IF (KZ(30).EQ.0) write(6,*) 'IMPACT: KZ(30)==0 L.397 quit'
#endif
      IF (KZ(30).EQ.0) GO TO 100
C*       Allow CHAIN only (#30 = -1) or TRIPLE & QUAD only (#30 < -1).
C      IF (KZ(30).EQ.-2.AND.KCHAIN.EQ.0) GO TO 100
*
      WHICH1 = ' TRIPLE '
      IF (JCOMP.GT.N) WHICH1 = ' QUAD   '
      IF (KCHAIN.GT.0) WHICH1 = ' CHAIN  '
*
      IF (rank.eq.0.and.H(IPAIR).GT.0.0) THEN
          WRITE (6,18)  I, JCOMP, ECC, ECC1, SEMI1, RIJ, GAMMA(IPAIR)
   18     FORMAT (' HYP CHAIN    I J E E1 A1 RIJ G  ',
     &                           2I6,2F7.3,1P,3E9.1)
      END IF
*
      IF (KZ(15).GT.1.OR.KZ(30).GT.1) THEN
          RI = SQRT((X(1,I) - RDENS(1))**2 +
     &              (X(2,I) - RDENS(2))**2 +
     &              (X(3,I) - RDENS(3))**2)
          VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
          PD = TWOPI*SEMI*SQRT(DABS(SEMI)/(BODY(I)+BODY(JCOMP)))
     &                             *TSTAR*365.24D6
          PD1 = TWOPI*SEMI1*SQRT(DABS(SEMI1)/BODY(I))*TSTAR*365.24D6
          if(rank.eq.0)
     &    WRITE (6,20)  WHICH1, TTOT, I, JCOMP, IPAIR, NAME(I1),
     &         NAME(I2), NAME(JCOMP), NAME(I), KSTAR(I1), KSTAR(I2),
     &         KSTAR(JCOMP), KSTAR(I), BODY(I1),BODY(I1)*ZMBAR,
     &         BODY(I2), BODY(I2)*ZMBAR, BODY(JCOMP),BODY(JCOMP)*ZMBAR,
     &         BODY(I)+BODY(JCOMP),(BODY(I)+BODY(JCOMP))*ZMBAR,
     &         R(IPAIR),R(IPAIR)*SU,H(IPAIR),SEMI,ECC,EB,PD,
     &         SEMI1,ECC1,EB1,PD1,PERT4, RIJ, RIJ*SU, PMIN, EB1/EB, 
     &         LIST(1,I1), RADIUS(I1)*SU,RADIUS(I2)*SU,
     &         RADIUS(JCOMP)*SU,RI,VI
   20     FORMAT (/,' NEW',A8,1P,' TIME',E12.5,' INCM',I9,
     &         ' I3',I9,' INPAIR',I9,' NAME(I1~3,INCM)',4I10,
     &         ' K*(I1~3,IN)',4I4,' M(I1)[NB,M*]',1P,E10.2,0P,F10.3,
     &         ' M(I2)[NB,M*]',1P,E10.2,0P,F10.3,
     &         ' M(I3)[NB,M*]',1P,E10.2,0P,F10.3,
     &         ' M(Icm)[NB,M*]',1P,E10.2,0P,F10.3,
     &         ' R12[NB,R*]',1P,2E10.2,' H[NB]',E10.2,
     &         ' SEMI0[NB]',E12.4,' ECC0',0P,F12.5,1P,' EB0[NB]',E10.2,
     &         ' P0[DAY]',E10.2,
     &         ' SEMI1[NB]',E12.4,' ECC1',0P,F12.5,1P,' EB1[NB]',E10.2,
     &         ' P1[DAY]',E10.2,
     &         ' GAMMA4',E10.3,'  RIN3[NB,R*]',2E10.2,
     &         ' PERIM[NB]',E8.1,'  EB1/EB0',
     &         E8.1,'  NP',I4,' RS(I1~3)[R*]',
     &         3E10.2,'  RI[NB]', E10.2,' VI[NB]',2E10.2)
          call flush(6)
      END IF
*
*     --01/03/14 13:46-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      call abort()
c$$$      print*,rank,'IMPACT I1',i1,'n',name(i1),'x',x(1,i1),
c$$$     &     'xdot',xdot(1,i1),'fdot',fdot(1,i1),'time',time
c$$$      call flush(6)
c$$$      call mpi_barrier(MPI_COMM_WORLD,ierr)
*     --01/03/14 13:46-lwang-end----------------------------------------*

      JCMAX = 0
*       Include any close single or c.m. perturber (cf. routine SETSYS).
      IF (JMAX.GT.0.AND.JMAX.NE.JCOMP) THEN
         IF(SQRT(RMAX2).LT.MIN(2.0D0*RSUM,RMIN)
     &        .AND.NAME(JMAX).GT.0) THEN
            IF (JCOMP.GT.N.AND.JMAX.GT.N) THEN
               JCMAX = 0
            ELSE
               if(rank.eq.0) THEN
                  WRITE (6,21)  NAME(JCOMP), NAME(JMAX), RSUM, 
     &                 SQRT(RMAX2)
 21               FORMAT (' B+2 CHAIN    NAM RSUM RMX ',2I7,1P,2E10.2)
               END IF
               CALL JPRED(JMAX,time,time)
               JCMAX = JMAX
*     ks MPI communication JCMAX
               call ksparmpi(K_store,K_int,K_JCMAX,0,0,JCMAX)
            END IF
         END IF
      END IF
*
*       Save global index of intruder for TRIPLE or CHAIN.
      JCLOSE = JCOMP
*
*       Check B-B interaction for switch of IPAIR & JPAIR or inert binary.
      IF (KCHAIN.GT.0.AND.JCOMP.GT.N) THEN
          K1 = 2*JPAIR - 1
          if(rank.eq.0)
     &    WRITE (6,22)  NAME(I1), NAME(I2), NAME(K1), NAME(K1+1),
     &                  KSTAR(I), KSTAR(JCOMP), ECC, ECC1, A0, SEMI2,
     &                  RIJ, SEMI1, PMIN
   22     FORMAT (' CHAIN B-B    NAME(I1-I4)',4I10,'  K*(INCM)',I3,
     &         '  K*(OCM)',I3,'  ECC0',F6.3,'  ECC1',F6.3,'  SEMI0',
     &         1P,E10.2'  SEMI2',E10.2,'  RIN3',E10.2,'  SEMI1',E10.2,
     &         '  PERIM',E10.2,0P)
          RT = 4.0*MAX(RADIUS(I1),RADIUS(I2))
          IF (SEMI0.LT.4.0*RT.AND.LIST(1,J1).EQ.0.OR.
     &        MIN(SEMI0,SEMI2).LT.0.01*RIJ) THEN
*       Ensure that widest binary comes first (more similar to triple).
              IF (SEMI0.LT.SEMI2) THEN
                  KPAIR = JPAIR
                  JPAIR = IPAIR
                  IPAIR = KPAIR
                  JCLOSE = N + JPAIR
              END IF
*       Check reduction of c.m. index (JPAIR becomes JPAIR - 1 if > IPAIR).
              IF (JPAIR.GT.IPAIR) JCLOSE = JCLOSE - 1
              IF (KZ(26).LT.2) THEN
*       Replace unperturbed near-synchronous binary by inert body in CHAIN.
                  JCOMP = 0
                  if(rank.eq.0)
     &            WRITE (6,25)  SEMI0, RIJ, R(JPAIR), GAMMA(JPAIR)
   25             FORMAT (' INERT BINARY    A RIJ R G ',1P,4E10.2)
              END IF
          ELSE
              JCLOSE = 0
          END IF
      END IF
*
*       Set phase indicator for calling TRIPLE or QUAD from MAIN.
      IPHASE = 4
      KSPAIR = IPAIR
*
*       Include the case of two interacting KS pairs.
      IF (JCOMP.GT.N) THEN
          IPHASE = 5
*       Switch pair indices and rename JCOMP if JPAIR has smaller size.
          IF (STEP(J1).LT.STEP(I1).AND.LIST(1,I1).GT.0) THEN
              KSPAIR = JPAIR
              JCOMP = I
              KS2 = IPAIR
          ELSE
              KS2 = JPAIR
          END IF
          IF (KZ(27).LE.0.AND.JPAIR.GT.IPAIR) THEN
              IF (JCLOSE.GT.0) JCLOSE = JCLOSE - 1
          END IF
*       Terminate smallest pair first and reduce second index if higher.
*         CALL KSTERM
C          IF (KS2.GT.KSPAIR) KS2 = KS2 - 1
*       Check whether KS2 is the last pair, if so, set the index as first pair since they will be exchanged
          IF (KS2.EQ.NPAIRS) KS2 = KSPAIR
          
      END IF
*
*       See whether chain regularization indicator should be switched on.
      IF (KCHAIN.GT.0) THEN
          IPHASE = 8
      END IF
*
*       Save KS indices and delay initialization until end of block step.
      CALL DELAY(KCHAIN,KS2)
*     ks MPI delay
      call ksparmpi(K_store,K_int,K_DELAY,0,0,0.0)
      flag_delay=.true.
c$$$*
c$$$*       Prepare procedure for chain between hierarchy and single body (9/99).
c$$$      IF (NAME(I).LT.0.AND.NAME(I).GE.-NZERO.AND.JCOMP.LE.N) THEN
c$$$*       Indentify merged ghost particle JG.
c$$$          CALL FINDJ(I1,JG,IM)
c$$$          if(rank.eq.0)
c$$$     &    WRITE (6,28)  NAME(I), NAME(JCOMP), NAME(JG),ECC1, PMIN, RIJ
c$$$   28     FORMAT (' HI CHAIN    NAM E1 PM RIJ ',I7,2I6,F7.3,1P,2E10.2)
c$$$          JJ = JCOMP
c$$$*       Terminate the merger in the usual way.
c$$$          KSPAIR = IPAIR
c$$$          IPHASE = 7
c$$$          CALL RESET
c$$$          ZMU = BODY(2*NPAIRS-1)*BODY(2*NPAIRS)/BODY(NTOT)
c$$$          EBCH0 = EBCH0 + ZMU*H(NPAIRS)
c$$$*       Specify chain indicator and define the two single particles.
c$$$          IPHASE = 8
c$$$          JCMAX = JG
c$$$          JCLOSE = JJ
c$$$          KSPAIR = NPAIRS
c$$$*       Set relevant variables in DELAY before terminating inner binary.
c$$$          CALL DELAY(KCHAIN,KS2)
c$$$          CALL DELAY(IPHASE,-1)
c$$$*       Initialize new chain of the 4 members JMAX, JCLOSE & KS components.
c$$$          ISUB = 0
c$$$          CALL CHAIN(ISUB)
c$$$*       Note that IPHASE = -1 now and INTGRT goes back to the beginning.
c$$$      ELSE IF (NAME(I).LT.-NZERO.OR.NAME(JCOMP).LT.0.OR.
c$$$     &        (NAME(I).LT.0.AND.JCOMP.GT.N)) THEN
c$$$*       Continue until KS termination on MERGE2 or merger with JCOMP > N.
c$$$          IPHASE = 0
c$$$      END IF
*
      GO TO 100
*
*       Begin check for merger of stable hierarchical configuration.
 30   MCF = .true.
#ifdef DECISION_CHECK
      write(6,*) 'IMPACT: Merger check'
#endif
      RA = SEMI1*(1.0 + ECC1)
      IF (SEMI1.LT.0.0) RA = RIJ
*
*       Identify formation of wide quadruple before merger is accepted.
      IF (JCOMP.GT.N.AND.ECC1.LT.1.0.AND.SEMI1.LT.0.1*RSCALE) THEN
          NNB = LISTQ(1) - 1
          K = 0
*       See whether current list contains first inner/outer component.
          NAM1 = NAME(2*JPAIR-1)
          DO 32 L = 2,NNB+2
              IF (NAM1.EQ.LISTQ(L)) K = K + 1
   32     CONTINUE
*       Generate diagnostics of first five outer orbits every half period.
          IF (K.LE.5.AND.TIME.GT.QCHECK.AND.KZ(15).GE.3) THEN
              ZMB = BODY(I) + BODY(JCOMP)
              RI = SQRT(RI2)
              TK = SEMI1*SQRT(SEMI1/ZMB)
              QCHECK = TIME + MIN(0.5*TWOPI*TK,0.1*TCR)
              TK = DAYS*TK
*       Check the stability criterion.
              QST = QSTAB(ECC,ECC1,ANGLE,BODY(I1),BODY(I2),BODY(JCOMP))
              PCR = QST*SEMIX
C              PCR = stability(BODY(I1),BODY(I2),BODY(JCOMP),ECC,ECC1,
C     &                                                  0.0D0)*SEMIX
              if(rank.eq.0)
     &        WRITE (89,33)  TTOT, NAME(2*IPAIR-1), NAM1, K, RI,
     &                       ECC1, EB, EB2, EB1, TK, PMIN, PCR
   33         FORMAT (' QUAD  Time[NB] NAME(I1) NAME(J1) LQ RI[NB] ',
     &             'ECC1 EB[NB] EB2[NB] EB1[NB] P1[days] PERIM[NB] ',
     &             'PCR[NB] ',1P,E20.11,0P,2I12,I4,F12.5,F8.4,1P,
     &             3E12.3,3E14.5)
              CALL FLUSH(89)
*       Remove two oldest members if list is too big.
              IF (NNB.GT.96) THEN
                  DO 34 K = 2,NNB
                      LISTQ(K) = LISTQ(K+2)
   34             CONTINUE
                  NNB = NNB - 2
              END IF
*       Add current names (inner & outer) at end and update membership.
              LISTQ(NNB+3) = NAME(2*IPAIR-1)
              LISTQ(NNB+4) = NAME(2*JPAIR-1)
              LISTQ(1) = NNB + 3
          END IF
      END IF
*
*       Allow temporary merger of inner part of extremely eccentric orbit.
      RFAC = 10.0*RMIN
      IF (ECC1.GT.0.99.AND.RA.GT.RFAC) THEN
          IF (RIJ.LT.0.1*SEMI1) RFAC = RA
      END IF
*
*       Increase apocentre tolerance to local scale factor for EB1 < EBS.
      EBS = 0.25*EBH/SQRT(1.0 + SQRT(RI2)/RSCALE)
      IF (EB1.LT.EBS) THEN
          H2 = (RC**2 + RI2)/FLOAT(NC+10)**0.66667
          RH = 6.0*SQRT(H2/CMSEP2)
          RFAC = MAX(RFAC,RH)
*       Extend maximum apocentre for massive systems (less perturbers).
          IF (BODY(I) + BODY(JCOMP).GT.10.0*BODYM) RFAC = 2.0*RFAC
      END IF
 
*       Skip merger for hyperbolic & soft binding energy or large apocentre.
*     ZF = 1.0
*     IF (BODY(I)*SMU.LT.0.4) ZF = 1.5
      IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
*       Use orbital velocity condition instead of binding energy for planets.
#ifdef DECISION_CHECK
          IF (SEMI1.GT.2.0*RMIN)
     &        write(6,*) 'IMPACT: L769'
#endif
          IF (SEMI1.GT.2.0*RMIN) GO TO 16
      ELSE IF (EB.GT.EBH.OR.EB1.GT.EBS.OR.
     &         (RA.GT.RFAC.AND.GI.GT.0.1)) THEN
C*         quit for hyperbolic case if perturber is merger or chain
C#ifdef DECISION_CHECK
C         IF(EB1.GT.0.0.AND.NAME(JCOMP).LE.0)
C     &        write(6,*) 'IMPACT: EB1(',EB1,')>0 && NAME(JCOMP)<0 L773'
C#endif
C         IF(EB1.GT.0.0.AND.NAME(JCOMP).LE.0) GO TO 16
#ifdef DECISION_CHECK
          write(6,*) 'IMPACT: EB(',EB,')>EBH(',EBH,')||EB1(',EB1,
     &     ')>EBS(',EBS,')||(RA(',RA,')>RFAC(',RFAC,') && ',
     &        'GI (',GI,')>0.1 L774'
#endif
          GO TO 16
       END IF
*
C*       Estimate the relative apocentre perturbations on body #I & JCOMP.
C      IF (ECC1.LT.0.95) THEN
C          PERT = PERT*(RA/RIJ)**3
C      ELSE
C          PERT = PERT*(ABS(SEMI1)/RIJ)**3
C      END IF
C      PERTA = PERT4*(RA/RIJ)**3
C*
*       Check tidal capture option (synchronous or evolving binary orbit).
      IF (KZ(27).GT.0) THEN
*       Skip merger if outer component would suffer tidal dissipation.
***       IF (SEMI1*(1.0 - ECC1).LT.4.0*RADIUS(JCOMP)) GO TO 100
*       Do not allow merger if Roche overflow or mass loss during next orbit.
          TK = TWOPI*SEMI1*SQRT(SEMI1/(BODY(I) + BODY(JCOMP)))
          TM = MIN(TEV(I1),TEV(I2),TEV(JCOMP),TEV(I))
          IF (KZ(34).GT.0.AND.TM - TIME.LT.STEPX) THEN
              RT = 5.0*MAX(RADIUS(I1),RADIUS(I2))
              IF (A0.LT.RT.OR.KSTAR(I).GT.0) THEN
                  CALL TRFLOW(IPAIR,DTR)
                  IF ((MOD(KSTAR(I),2).EQ.1.AND.DTR.LT.STEPX).OR.
     &                 DTR.LT.TK) GO TO 100
              END IF
              IF (JCOMP.GT.N.AND.KSTAR(JCOMP).GT.0) THEN
                  CALL TRFLOW(JPAIR,DTR)
                  IF(MOD(KSTAR(JCOMP),2).EQ.1.OR.DTR.LT.STEPX) GOTO 100
              END IF
          END IF
*
*       Ensure SLEEP for circularizing binary with TCIRC > 1.0E+06.
          QPERI = A0*(1.0 - ECC)
          RM = MAX(RADIUS(I1),RADIUS(I2))
          IF (KSTAR(I).EQ.-2.AND.QPERI.GT.10.0*RM) THEN
              ICIRC = -1
              CALL INDUCE(IPAIR,JCOMP,EMAX,EMIN,ICIRC,TC,ANGLE,TG,EDAV)
              IF (TC.GT.1.0D+06) THEN
                  if(rank.eq.0)
     &            WRITE (6,35)  NAME(I1), KSTAR(I1), KSTAR(I2), ECC,
     &                          EMAX, QPERI/RM, EDAV, a0, PMIN, TC
   35             FORMAT (' IMPACT SLEEP    NM K* E EX QP/R ED A PM TC',
     &                                      I7,2I4,2F8.4,F6.1,1P,4E10.2)
                  KSTAR(I) = 0
*     ks MPI communication
                  call ksparmpi(K_store,K_int,K_KSTAR,I,0,KSTAR(I))
                  NSLP = NSLP + 1
                  II = -I
                  CALL SPIRAL(II)
              END IF
          END IF
*
*       Delay merger for recently updated standard binary and short TCIRC.
          DT = MIN(TEV(I1),TEV(I2)) - TIME
          IF (KSTAR(I).EQ.0.AND.NAME(I).GT.0.AND.DT.LT.TK) THEN
              ICIRC = -1
              CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
*             if(rank.eq.0)
*    &        WRITE (6,36) NAME(I1), ECC, TTOT, RADIUS(I1)*SU, QPERI, TC
*  36         FORMAT (' TCIRC    NAM E T R* QP TC ',
*    &                           I6,F7.3,F8.3,F7.1,1P,2E10.2)
*       Beware possible termination by routine HMDOT using QPERI < 3*RADIUS.
              IF (TC.LT.2000.0.AND.ECC.GT.0.002) GO TO 100
          END IF
          IF (KZ(19).GE.3) THEN
              IF (MIN(TEV(I1),TEV(I2)).LT.TIME + TK) GO TO 100
          END IF
*       Skip chaotic binary (KSTAR = -2 is treated below).
          IF (KSTAR(I).EQ.-1.OR.KSTAR(JCOMP).EQ.-1) GO TO 100
      END IF
*
CCC Not sure whether this is good, suppress for test (L.Wang)      
C*       Ensure consistency of estimated perturbations with termination.
C      PERT = PERT + PERTA
C      IF (NP.LE.3) THEN
C          NP = LIST(1,I)
C*       Copy neighbour list (routine FPERT skips body #JCOMP).
C          DO 38 L = 1,NP
C              JLIST(L) = LIST(L+1,I)
C   38     CONTINUE
C      END IF
C*
C*       Allow highly eccentric outer orbit (estimated PERT may be excessive).
C      IF (ECC1.GT.0.98.AND.RIJ.LT.0.1*SEMI1) THEN
C          PERT = PERT4
C          GO TO 40
C      END IF
C*
C*       Evaluate the actual perturbation.
C      CALL FPERT(I,JCOMP,NP,PERT2)
C      PERT2 = PERT2*RJMIN2/(BODY(I) + BODY(JCOMP))
C      GI = PERT2*(RA/RIJ)**3
C      IF (PERT4.GT.GMAX.OR.GI.GT.0.05.OR.RIJ.GT.6*RMIN) GO TO 100
C
 40   MCF = .true.
*       Skip merger if an outer binary is fairly perturbed or not hard.
      IF (JCOMP.GT.N) THEN
          IF (GAMMA(JPAIR).GT.1.0E-03.OR.EB2.GT.EBH) THEN
#ifdef DECISION_CHECK
             write(6,*) 'IMPACT: GAMMA(JPAIR)(',GAMMA(JPAIR),')>1e-3',
     &            ' || EB2(',EB2,')>EBH(',EBH,') L756'
#endif
             GO TO 16
          END IF
       END IF
*
#ifdef DECISION_CHECK
      IF (PERT.GT.GMAX.OR.GI.GT.0.02) 
     &     write(6,*) 'IMPACT: PERT(',PERT,') > GMAX(',GMAX,
     &     ') || GI (',GI,') > 0.02; ECC1(',ECC1,
     &     ') SEMI1(',SEMI1,') L824'
#endif
      IF (PERT.GT.GMAX.OR.GI.GT.0.02) GO TO 16

#ifdef DECISION_CHECK
      IF (GI.GT.0.1.AND.RIJ.GT.RMIN)
     &     write(6,*) 'IMPACT: GI(',GI,')>0.1 && RIJ(',RIJ,')>RMIN(',
     &     RMIN,') L934'
#endif
      IF (GI.GT.0.1.AND.RIJ.GT.RMIN) GO TO 16

#ifdef DECISION_CHECK
      IF (GI.GT.0.01.AND.RIJ.GT.5.0*RMIN)
     &     write(6,*) 'IMPACT: GI(',GI,')>0.01 && RIJ(',RIJ,')>5RMIN(',
     &     RMIN,') L942'
#endif
      IF (GI.GT.0.01.AND.RIJ.GT.5.0*RMIN) GO TO 16

#ifdef DECISION_CHECK
      IF (GI.GT.0.25) write(6,*) 'IMPACT: GI(',GI,')>0.25 L947'
#endif
      IF (GI.GT.0.25) GO TO 16
      
*       Switch to direct integration for planetary systems if GI > 1D-04.
      IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
#ifdef DECISION_CHECK
          IF (GI.GT.1.0D-04) 
     &        write(6,*) 'IMPACT: L831'
#endif
          IF (GI.GT.1.0D-04) GO TO 16
      END IF

      IF (NMERGE.EQ.MMAX - 1) THEN
          IF (NWARN.LT.1000) THEN
              NWARN = NWARN + 1
              if(rank.eq.0)
     &        WRITE (6,41)  NMERGE
   41         FORMAT (5X,'WARNING!    MERGER LIMIT    NMERGE =',I4)
          END IF
#ifdef DECISION_CHECK
          write(6,*) 'IMPACT: L846'
#endif
          GO TO 16
      END IF
*
*       Do not allow merger in the inner region of perturbed eccentric orbit.
      IF (RIJ.LT.SEMI1.AND.LIST(1,I1).GT.0) THEN
*       Note: moved down from label 30 with 0.1*SEMI1 replacing 2*PMIN.
          IF (ECC1.GT.0.95.AND.RIJ.LT.0.1*SEMI1) THEN
#ifdef DECISION_CHECK
             write(6,*) 'IMPACT: ECC1 (',ECC1,')>0.95 && RIJ(',
     &            RIJ,')<0.1*SEMI1(',SEMI1,') L856 quit'
#endif
              GO TO 100
          END IF
      END IF
*
*     -----------------------------------------------------------------------
*       Form coefficients for stability test (Valtonen, Vistas Ast 32, 1988).
*     AM = (2.65 + ECC)*(1.0 + BODY(JCOMP)/BODY(I))**0.3333
*     FM = (2.0*BODY(JCOMP) - BODY(I))/(3.0*BODY(I))
*       Note that routine NEWTEV in MERGE2 replaces suppressed part.
*     IF (KZ(19).GE.3) THEN
*         TM = MIN(TEV(I1),TEV(I2),TEV(JCOMP),TEV(I))
*         IF (MIN(NAME(I),NAME(JCOMP)).LT.0.AND.TM-TIME.LT.0.2) THEN
*             GO TO 100
*         END IF
*     END IF
*
*       Expand natural logarithm for small arguments.
*     IF (ABS(FM).LT.0.67) THEN
*         BM = FM*(1.0 - (0.5 - ONE3*FM)*FM)
*     ELSE
*         BM = LOG(1.0D0 + FM)
*     END IF
*
*       Adopt mass dependent criterion of Harrington (A.J. 82, 753) & Bailyn.
*     PCRIT = AM*(1.0 + 0.7*BM)*SEMI
*     -----------------------------------------------------------------------
*
*       Form hierarchical stability ratio (Eggleton & Kiseleva 1995).
*     QL = BODY(I)/BODY(JCOMP)
*     Q1 = MAX(BODY(I2)/BODY(I1),BODY(I1)/BODY(I2))
*     Q3 = QL**0.33333
*     Q13 = Q1**0.33333
*     AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*     EK = AR*SEMI*(1.0D0 + ECC)
*
*       Choose the most dominant triple in case of two binaries.
      IF (JCOMP.GT.N) THEN
*       Adopt 10% fudge factor with linear dependence on smallest ratio.
          YFAC = 1.0 + 0.1*MIN(SEMI2/SEMI,SEMI/SEMI2)
      ELSE
          YFAC = 1.0
      END IF
*
*       Prepare inclination evaluation for triple or widest inner binary.
      IF (JCOMP.GT.N) THEN
*       Ensure widest inner binary (swap is OK for termination or ZARE).
          IF (SEMI.LT.SEMI2) THEN
              ECC2 = (1.0 - R(JPAIR)/SEMI2)**2 +
     &                             TDOT2(JPAIR)**2/(BODY(JCOMP)*SEMI2)
              ECC = SQRT(ECC2)
              KPAIR = IPAIR
              IPAIR = JPAIR
              JPAIR = KPAIR
              I1 = 2*IPAIR - 1
              I2 = I1 + 1
              JJ = I
              I = JCOMP
              JCOMP = JJ
              SEMIZ = SEMI2
              SEMI2 = SEMI
              SEMI = SEMIZ
          END IF
      END IF
*
*       Resolve binary (even perturbed KS not always done).
      IF (ECC1.LT.1.0) THEN
          CALL RESOLV(IPAIR,1)
      END IF
*
*       Employ the basic stability criterion for fast check (ECC1 < 1).
*     IF (ECC1.LT.1.0) THEN
*         Q = BODY(JCOMP)/BODY(I)
*         XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
*         PCRIT = 2.8*XFAC**0.4*SEMI*(1.0 - 0.6*ANGLE/TWOPI)
*         IF (PCRIT.GT.2.0*PMIN) GO TO 100
*     END IF
*
*       Evaluate the general stability function (Mardling MNRAS, 2008).
      IF (ECC1.LT.1.0.AND.YFAC.LT.1.02) THEN
          BJ = BODY(JCOMP)
          EOUT = ECC1
*       Increase tolerance near sensitive stability boundary (RM 10/2008).
          IF (EOUT.GT.0.8) THEN
              DE = 0.5*(1.0 - EOUT)
              DE = MIN(DE,0.01D0)
*       Evaluate outer eccentricity derivative due to dominant perturber.
              IF (JMAX.NE.JCOMP.AND.PERT.GT.GMIN) THEN
                  CALL EDOT(I,JCOMP,JMAX,SEMI1,ECC1,ECCDOT)
*       Include a small effect of positive derivative over 10 orbits.
                  IF (ECCDOT.GT.0.0) THEN
                      TK1 = TWOPI*SEMI1*
     &                      SQRT(SEMI1/(BODY(I) + BODY(JCOMP)))
                      DE = DE - 10.0*MIN(ECCDOT*TK1,0.001D0)
*                     WRITE (6,443)  ECC1, DE, ECCDOT, TK1, PERT
* 443                 FORMAT (' EDOT!!    E1 DE ED TK G ',
*    &                                    2F9.5,1P,3E9.1)
                  END IF
              END IF
*       Allow extra tolerance after 1000 tries. (suppress to avoid
*       inconsistent of different processors in parallel)
c$$$              IF (NMTRY.GE.1000) DE = MIN(1.0D0 - EOUT,0.02D0)
CC              EOUT = EOUT - DE
              EOUT = MIN(EOUT - DE,0.9999D0)
              PMIN = SEMI1*(1.0 - EOUT)
          END IF
C          NST = NSTAB(SEMI,SEMI1,ECC,EOUT,ANGLE,BODY(I1),BODY(I2),BJ)
          QST = QSTAB(ECC,EOUT,ANGLE,BODY(I1),BODY(I2),BJ)
          RP = PMIN/SEMI
*     Using the new qstab (L.Wang)
          IF (QST.LT.RP) THEN
*     --10/08/14 21:40-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
C          print*,'NST SEMI',SEMI,'SEMI1',SEMI1,'ECC',ECC,'EOUT',EOUT,
C     &         'ANGLE',ANGLE,'BODY1',BODY(I1),'BODY(I2)',BJ
*     --10/08/14 21:40-lwang-end-debug----------------------------------*
c          IF (NST.EQ.0) THEN
c              PCRIT = 0.98*PMIN*(1.0 - PERT)
c              PCR = stability(BODY(I1),BODY(I2),BODY(JCOMP),ECC,EOUT,
c     &                                                          ANGLE)
c              PCR = PCR*SEMI
*       Reduce the limit slightly, also allowing for perturbation.
              PCRIT = 0.99*PMIN*(1.0 - PERT)
              PCR = QST*SEMI
*       Specify reduced peri if old criterion < PCRIT/2 (avoids switching).
              IF (PCR.LT.0.5*PCRIT) THEN
                  PCRIT = 0.75*PCRIT
              END IF
              IF (PCRIT.LT.YFAC*PCR.AND.PERT.LT.0.02.AND.
     &            NMTRY.LT.10) THEN
                  ALPH = 360.0*ANGLE/TWOPI
                  FAIL = PMIN*(1-PERT) - YFAC*PCR
                  WRITE (6,43)  TTOT, ALPH, ECC, ECC1, PMIN, FAIL, PERT
   43             FORMAT (' NEWSTAB    T INC EI EO PMIN FAIL PERT ',
     &                                 F7.1,F7.2,2F8.4,1P,3E10.2)
              END IF
c*     if(rank.eq.0)then
c*     WRITE (57,444)  BODY(I1),(X(K,I1),K=1,3),(XDOT(K,I1),K=1,3)
c*     WRITE (57,444)  BODY(I2),(X(K,I2),K=1,3),(XDOT(K,I2),K=1,3)
c*     WRITE (57,444)BODY(JCOMP),(X(K,JCOMP),K=1,3),(XDOT(K,JCOMP),K=1,3)
c* 444 FORMAT (' ',1P,E14.6,1P,3E18.10,3E14.6)
c*     end if
c              IF (NMTRY.GE.1000) THEN
c*                 if(rank.eq.0)
c*    &            WRITE (6,44)  TTOT, NAME(I1), ECC1, EOUT, PCRIT, PMIN
c*  44             FORMAT (' MARGINAL STAB    T NM E1 EOUT PCR PMIN ',
c*    &                                       F7.1,I7,2F8.4,1P,2E10.2)
c                  NMTRY = 0
c              END IF
*       Increase the limit slightly, also allowing for perturbation.
          ELSE
              NMTRY = NMTRY + 1
              PCRIT = 1.01*PMIN
          END IF
      ELSE
*     Not good, used original way
C          PCRIT = PMIN
          PCR = stability(BODY(I1),BODY(I2),BODY(JCOMP),ECC,ECC1,ANGLE)
          PCRIT = PCR*SEMI
      END IF
*
*       Check whether the main perturber dominates the outer component.
      IF (JMAX.NE.JCOMP) THEN
          RIJ2 = (X(1,JMAX) - X(1,JCOMP))**2 +
     &           (X(2,JMAX) - X(2,JCOMP))**2 +
     &           (X(3,JMAX) - X(3,JCOMP))**2
          FMAX = (BODY(JMAX) + BODY(JCOMP))/RIJ2
#ifdef DECISION_CHECK
          IF (FMAX.GT.(BODY(I) + BODY(JCOMP))/RJMIN2)
     &         write(6,*) 'IMPACT: L1016'
#endif
          IF (FMAX.GT.(BODY(I) + BODY(JCOMP))/RJMIN2) GO TO 16
      END IF
*
*       Determine time-scale for stability (absolute or approximate).
      PM1 = PMIN*(1.0 - 2.0*PERT)
      CALL TSTAB(I,ECC1,SEMI1,PM1,YFAC,ITERM)
#ifdef DECISION_CHECK
      IF (ITERM.GT.0) 
     &     write(6,*) 'IMPACT: ITERM(',ITERM,')>0 L1026'
#endif
      IF (ITERM.GT.0) GO TO 100
*
*       Check perturbed stability condition.
#ifdef DECISION_CHECK
      IF (PMIN*(1.0 - PERT).LT.YFAC*PCRIT) 
     &     write(6,*) 'IMPACT: PMIN(',PMIN,')*(1-PERT(',
     &     PERT,')<YFAC(',YFAC,')*PCRIT(',PCRIT,') L1033'
#endif
      IF (PMIN*(1.0 - PERT).LT.YFAC*PCRIT) GO TO 16
*
*       Extend active Roche case up to end of look-up time (bug fix 01/09).
      IF (KSTAR(I).GE.11.AND.MOD(KSTAR(I),2).NE.0) THEN
          DT = TEV(I) - TIME
#ifdef DECISION_CHECK
          IF (DT.LT.10.0*STEPX) 
     &     write(6,*) 'IMPACT: L1042'
#endif
          IF (DT.LT.10.0*STEPX) GO TO 100
*     KS MPI communication for tmdis
          TMDIS(NMERGE+1) = MIN(TIME + DT, TMDIS(NMERGE+1))
          call ksparmpi(K_store,K_real8,K_TMDIS,NMERGE+1,0,
     &         TMDIS(NMERGE+1))
      END IF
*
*       Obtain growth time and inclination for significant eccentricity.
      IF (SEMI1.GT.0.0.AND.ECC.GT.0.1) THEN
          ICIRC = -1
          CALL INDUCE(IPAIR,JCOMP,EMAX,EMIN,ICIRC,TC,ANGLE,TG,EDAV)
*       Prevent termination for TC < 2000 in HMDOT but allow small EMAX & DE.
          IF (KZ(27).EQ.2.AND.TC.LT.2000.0) THEN
              DE = ABS(EMAX - ECC)
              DT = MIN(TEV(I1),TEV(I2),TEV(I)) - TIME
*       Enforce an update next block-step for small remaining times.
              IF (DT.LT.0.1.AND.MOD(KSTAR(I),2).EQ.0) THEN
                  TEV(I1) = TIME + 0.1
                  TEV(I2) = TIME + 0.1
*     ks MPI COMMUNICATION tev
                  call ksparmpi(K_store,K_real8,K_TEV,i1,0,TEV(I1))
                  call ksparmpi(K_store,K_real8,K_TEV,i2,0,TEV(I2))
                  if(rank.eq.0)
     &            WRITE (6,46)  TIME+TOFF, NAME(I1), KSTAR(I), ECC, TC
   46             FORMAT (' ENFORCED TEV UPDATE    T NM K* E TC ',
     &                                             F9.2,I8,I4,F8.4,F7.1)
              END IF
*       Accept KSTAR = 10, coasting Roche or small eccentricity.
              IF (MOD(KSTAR(I),2).EQ.0) THEN
                  CALL TRFLOW(IPAIR,DTR)
                  IF (DTR.LT.STEP(I)) THEN
                      TEV(I) = TIME + DTR
*     ks MPI COMMUNICATION tev
                      call ksparmpi(K_store,K_real8,K_TEV,I,0,TEV(I))
                  END IF
              ELSE
#ifdef DECISION_CHECK
                 IF (EMAX.GT.0.8.OR.DE.GT.0.2.OR.DT.LT.0.1)
     &                write(6,*) 'IMPACT: L1082'
#endif
                  IF (EMAX.GT.0.8.OR.DE.GT.0.2.OR.DT.LT.0.1) GO TO 100
              END IF
          END IF
      END IF
*
*       Check circularization and dissipation time (exclude Roche stages).
      IF (KZ(27).EQ.2.AND.KSTAR(I).LT.10.AND.ECC1.LT.1.0.AND.
     &    NAME(I).GT.0) THEN
          ECC0 = ECC
          CALL DECIDE(IPAIR,SEMI,ECC,EMAX,EMIN,TG,EDAV,IQ)
#ifdef DECISION_CHECK
          IF (IQ.GT.0.OR.IPHASE.LT.0) 
     &     write(6,*) 'IMPACT: L1096'
#endif
          IF (IQ.GT.0.OR.IPHASE.LT.0) GO TO 100
          TK1 = TWOPI*SEMI1*SQRT(SEMI1/(BODY(I) + BODY(JCOMP)))
*         IF (TMDIS(NMERGE+1) - TIME.LT.TK1) GO TO 100
*         EK = EK*(1.0 + ECC)/(1.0 + ECC0)
          PCRIT = PCRIT*(1.0 + ECC)/(1.0 + ECC0)
      END IF
*
*       Perform safety check on radii for case of no circularization.
      IF (KZ(27).EQ.0) THEN
          IF (SEMI*(1.0 - ECC).LT.2.0*MAX(RADIUS(I1),RADIUS(I2))) THEN
              IF (KZ(19).EQ.0) THEN
                  GO TO 100
              END IF
          END IF
      END IF
      
*
*       Include rare case of circularizing outer binary.
      IF (KZ(27).EQ.2.AND.JCOMP.GT.N.AND.KSTAR(JCOMP).EQ.-2) THEN
          ECC2 = (1.0 - R(JPAIR)/SEMI2)**2 +
     &                               TDOT2(JPAIR)**2/(BODY(JCOMP)*SEMI2)
          ECCJ = SQRT(ECC2)
*       See whether to reduce the look-up time TMDIS (no skip here).
          CALL DECIDE(JPAIR,SEMI2,ECCJ,EMAXJ,EMIN,TG,EDAV,IQ)
          IF (IPHASE.LT.0) GO TO 100
      END IF
*
*       Check Zare exchange stability criterion and create diagnostics.
      IF (SEMI1.GT.0.0) THEN
          CALL ZARE(I1,I2,SP,JCOMP)
          Q = BODY(JCOMP)/BODY(I)
*       Note inclination is determined by routine INCLIN for ECC < 0.1.
          IF (SP.LT.1.0.AND.ANGLE.LT.0.17) THEN
              IZARE = IZARE + 1
              IF (IZARE.LT.200) THEN
              if(rank.eq.0)then
              WRITE (6,48)  TTOT, Q, ECC, ECC1, SEMI, PMIN, PCRIT,
     &                      YFAC, SP
              WRITE (6,47) I,JCOMP,N,I1,I2,RIJ,SEMI1
   47         FORMAT (' I JCOMP N I1 I2 RIJ A1   ',5I6,1P,2E10.2)
              CALL FLUSH(7)
              WRITE (6,48)  TTOT, Q, ECC, ECC1, SEMI, PMIN, PCRIT,
     &                      YFAC, SP
   48         FORMAT (' ZARE TEST  Time[NB] M(I3)/M(INCM) ECC0 ECC1 ',
     &             'SEMI PERIM  PCR YF SP ',
     &             F8.2,F5.1,2F7.3,1P,3E9.1,0P,2F6.2)
              end if
              END IF
          END IF
          if(rank.eq.0)
     &    WRITE (73,49)  TTOT, Q, ECC, ECC1, SEMI, PMIN, PCRIT,
     &                   TG, SP, 360.0*ANGLE/TWOPI, KSTAR(I)
   49     FORMAT ('STAB   Time[NB] M(I3)/M(INCM) ECC0 ECC1 SEMI0[NB] ',
     &         'PERIM[NB] PCR[NB] TG[NB] SP INA[deg] K* ',
     &         1P,E19.10,0P,F12.5,2F10.5,1P,4E15.6,0P,F12.4,F7.1,I4)
          CALL FLUSH(73)
*         IF (KSTAR(I1).GE.10) TEV(I1) = 1.0E+10
*         IF (KSTAR(I2).GE.10) TEV(I2) = 1.0E+10
      END IF
*
*       Specify the final critical pericentre using the fudge factor.
      PCRIT = YFAC*PCRIT
*
      IF (NAME(I).GT.N.AND.NAME(JCOMP).GT.N.AND.ECC1.GT.1.0) GO TO 100
*       Skip if #JCOMP is a chain c.m. but allow bound double hierarchy.
      IF (NAME(JCOMP).EQ.0) GO TO 100
      IF (ECC1.GT.1.0.AND.MIN(NAME(I),NAME(JCOMP)).LT.0) GO TO 100
      DO 55 ISUB = 1,NSUB
          IF (NAME(JCOMP).EQ.NAMES(1,ISUB)) GO TO 100
   55 CONTINUE
*       Do not allow the formation of a SEPTUPLET.
*     IF ((NAME(I).LT.-2*NZERO.AND.JCOMP.GT.N).OR.
*    &     NAME(JCOMP).LT.-2*NZERO) GO TO 100
*
*       Include diagnostics for double hierarchy or optional standard case.
      IF (NAME(I).LT.0.OR.NAME(JCOMP).LT.0) THEN
          IF (KZ(15).GT.1) THEN
              WHICH1 = ' MERGE2 '
          RI = SQRT((X(1,I) - RDENS(1))**2 +
     &              (X(2,I) - RDENS(2))**2 +
     &              (X(3,I) - RDENS(3))**2)
          VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
          PD = TWOPI*SEMI*SQRT(DABS(SEMI)/(BODY(I)+BODY(JCOMP)))
     &                             *TSTAR*365.24D6
          PD1 = TWOPI*SEMI1*SQRT(DABS(SEMI1)/BODY(I))*TSTAR*365.24D6
          if(rank.eq.0)
     &    WRITE (6,20)  WHICH1, TTOT, I, JCOMP, IPAIR, NAME(I1),
     &         NAME(I2), NAME(JCOMP), NAME(I), KSTAR(I1), KSTAR(I2),
     &         KSTAR(JCOMP), KSTAR(I), BODY(I1),BODY(I1)*ZMBAR,
     &         BODY(I2), BODY(I2)*ZMBAR, BODY(JCOMP),BODY(JCOMP)*ZMBAR,
     &         BODY(I)+BODY(JCOMP),(BODY(I)+BODY(JCOMP))*ZMBAR,
     &         R(IPAIR),R(IPAIR)*SU,H(IPAIR),SEMI,ECC,EB,PD,
     &         SEMI1,ECC1,EB1,PD1,PERT4, RIJ, RIJ*SU, PMIN, EB1/EB, 
     &         LIST(1,I1), RADIUS(I1)*SU,RADIUS(I2)*SU,
     &         RADIUS(JCOMP)*SU,RI,VI
          END IF
*       Note rare case of two hierarchies merging and identify ghost names.
          IF (NAME(I).LT.0.AND.NAME(JCOMP).LT.0) THEN
              CALL FINDJ(I1,JI,IM)
              J1 = 2*JPAIR - 1
              CALL FINDJ(J1,JJ,JM)
              if(rank.eq.0)
     &        WRITE (6,60)  NAME(I1), NAME(JI), NAME(I1+1), NAME(J1),
     &                      NAME(JJ), NAME(J1+1), ECC, ECC1, SEMI,
     &                      SEMI1, PMIN, PCRIT
   60         FORMAT (' HI MERGE    NAM E E1 A A1 PM PC ',
     &                              6I7,2F7.3,1P,4E10.2)
          END IF
      ELSE IF (KZ(15).GT.1) THEN
          WHICH1 = ' MERGER '
          IF (JCOMP.GT.N) WHICH1 = ' QUAD   '
          RI = SQRT((X(1,I) - RDENS(1))**2 +
     &              (X(2,I) - RDENS(2))**2 +
     &              (X(3,I) - RDENS(3))**2)
          VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
          PD = TWOPI*SEMI*SQRT(DABS(SEMI)/(BODY(I)+BODY(JCOMP)))
     &                             *TSTAR*365.24D6
          PD1 = TWOPI*SEMI1*SQRT(DABS(SEMI1)/BODY(I))*TSTAR*365.24D6
          if(rank.eq.0)
     &    WRITE (6,20)  WHICH1, TTOT, I, JCOMP, IPAIR, NAME(I1),
     &         NAME(I2), NAME(JCOMP), NAME(I), KSTAR(I1), KSTAR(I2),
     &         KSTAR(JCOMP), KSTAR(I), BODY(I1),BODY(I1)*ZMBAR,
     &         BODY(I2), BODY(I2)*ZMBAR, BODY(JCOMP),BODY(JCOMP)*ZMBAR,
     &         BODY(I)+BODY(JCOMP),(BODY(I)+BODY(JCOMP))*ZMBAR,
     &         R(IPAIR),R(IPAIR)*SU,H(IPAIR),SEMI,ECC,EB,PD,
     &         SEMI1,ECC1,EB1,PD1,PERT4, RIJ, RIJ*SU, PMIN, EB1/EB, 
     &         LIST(1,I1), RADIUS(I1)*SU,RADIUS(I2)*SU,
     &         RADIUS(JCOMP)*SU,RI,VI
      END IF
*
*       Check for diagnostic output of quadruples.
      IF (SEMI1.GT.0.0.AND.JCOMP.GT.N.AND.KZ(15).GE.3) THEN
          ZMB = BODY(I) + BODY(JCOMP)
          TK = DAYS*SEMI1*SQRT(SEMI1/ZMB)
          if(rank.eq.0)
     &    WRITE (89,65)  TTOT, NAME(2*IPAIR-1), NAME(2*JPAIR-1),
     &                   LISTQ(1), SQRT(RI2), ECC1, EB, EB2, EB1,
     &                   TK, PMIN, PCRIT
   65     FORMAT (' QUAD#   T NAM LQ RI E1 EB EB2 EB1 P1 PM PC ',
     &                      F8.1,2I6,I4,F6.2,F8.4,1P,3E12.3,3E9.1)
      END IF
*
*       Generate a diagnostic file of stable hierarchies (suppressed).
C      IF (ECC1.LT.-1.0) THEN
C          RI = SQRT(RI2)/RC
C          if(rank.eq.0)
C     &    WRITE (80,70)  TPHYS, RI, NAME(JCOMP), QL, Q1, ECC, ECC1,
C     &                   SEMI, SEMI1, PCRIT/PMIN, 360.0*ANGLE/TWOPI,EMAX
C   70     FORMAT (F8.1,F5.1,I6,2F6.2,2F6.3,1P,2E10.2,0P,F5.2,F6.1,F6.3)
C          CALL FLUSH(80)
C      END IF
*
*       Copy pair index and set indicator for calling MERGE from MAIN.
      KSPAIR = IPAIR
      IPHASE = 6
*****   Cannot call delay with -2 here in parallel ks
c$$$      CALL DELAY(IPHASE,-2)
*
*       Save KS indices and delay merger until end of block step.
      CALL DELAY(KS2,KS2)
*     ks MPI delay
      if(.not.flag_delay) call ksparmpi(K_store,K_int,K_DELAY,0,0,0.0)
      if(.not.flag_delay) call ksparmpi(K_store,K_int,K_DELAY,0,0,0.0)
*
 100  IF (IPHASE.NE.8) JCMAX = 0
*
      RETURN
*
      END
