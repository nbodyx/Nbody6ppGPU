      SUBROUTINE NBINT_COR(I,FIRR,FD)
*
*
*       Irregular corrector for SIMD version
*       ----------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/NEWXV/ XN(3,NMAX),XNDOT(3,NMAX)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FREG(3),FD(3),FDUM(3)
*      LOGICAL RPRED(NMAX)
*
*       Obtain irregular force & first derivative.
*      if (.not.NOPRED(I)) RPRED(I)=.false.
*      if (NOPRED(I).and.RPRED(I)) print*,rank,K,'not predict',time
*       NOPRED(I) = .true.
      call jpred_int(I,TIME)
      XI(1:3) = X(1:3,I)
      XIDOT(1:3) = XDOT(1:3,I)
*          PHII(I) = 0.D0
      NNB0 = LIST(1,I)
!$omp critical      
      NIRRF = NIRRF + NNB0
      NBPRED = NBPRED + NNB0
!$omp end critical      
*
*       Assume small mass at centre for special case of no neighbours.
      IF (NNB0.EQ.0) THEN
          RI2 = XI(1)**2 + XI(2)**2 + XI(3)**2
          FIJ = 0.01*BODYM/(RI2*SQRT(RI2))
          RDOT = 3.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) +
     &                                 XI(3)*XIDOT(3))/RI2
          DO 10 K = 1,3
              FIRR(K) = -FIJ*XI(K)
              FD(K) = -(XIDOT(K) - RDOT*XI(K))*FIJ
   10     CONTINUE
          IF (I.GT.N) IPAIR = I - N
          GO TO 70
      END IF
*
*       Choose force loop for single particle or regularized c.m. body.
      IF (I.LE.N) GO TO 20
*     Set KS pair index.
      IPAIR = I - N
*       Adopt c.m. approximation for small total perturbation.
      I1 = 2*IPAIR - 1
      IF (LIST(1,I1).GT.0) THEN
*     --03/03/14 20:51-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(time.ge.14.968627929687500.and.name(i).eq.18207) then
c$$$         print*,'bcmf',rank,'I',I,'n',name(i),'fi',fi(1,i),
c$$$     &        'fin',firr(1),'s',step(i),'fd',fidot(1,i),
c$$$     &        'fdn',fd(1),'nb',list(1,i),'pn',list(1,2*(i-n)-1)
c$$$         call flush(6)
c$$$      end if
*     --03/03/14 20:51-lwang-end----------------------------------------*
         call cmfirr_cor(I,I1,FIRR,FD)
*     --03/30/13 22:34-lwang-debug--------------------------------------*
c$$$      if (time.ge.0.3.and.name(i).eq.6193) then
c$$$         write(111+rank,*),'anb',I,name(i),FIRR,FD,TIME,nnb0,list(1,i1)
c$$$*      if (rank.eq.0.and.nstepr.gt.4314800) then
c$$$*         write (6,*) name(i), nnb,time, freg(1)
c$$$* 118     format ('name, nnb, freg, time', 2I6, 1F2.6, 1E10.2)
c$$$         call flush(111+rank)
c$$$      end if
c$$$*         write(111+rank,*) 'nb',I,i1,FIRR,FD,TIME,nnb0,list(1,i1)
*     --10/15/13 19:05-lwang-end----------------------------------------*
          GO TO 70
      END IF
*
 20   IF (NPAIRS.GT.0) then
         call cmfirr_ucor(I,NNB0,LIST(1,I),FIRR,FD)
      END IF
*     --03/07/14 20:37-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if (time.eq.4.2424011230468750.and.name(i).eq.2602) then
c$$$         write(6,*),rank,'AC I',I,'N',name(i),'FI',FIRR,'FD',FD,'T',TIME
c$$$     &        ,'NB',nnb0
c$$$         call flush(6)
c$$$      end if
*     --03/07/14 20:37-lwang-end----------------------------------------*
*     
*       Include treatment for regularized clump.
      IF (NCH.GT.0) THEN
*       Distinguish between chain c.m. and any other particle.
*       Note: in NBODY6++ CHFIRR and FCHAIN are called with IR=1 since
*             chain prediction and perturber list are updated in integrator.
          IF (NAME(I).EQ.0) THEN
              CALL CHFIRR(I,1,XI,XIDOT,FIRR,FD)
          ELSE
*       See if chain perturber list contains body #I.
              NP1 = LISTC(1) + 1
              DO 65 L = 2,NP1
                  J = LISTC(L)
                  IF (J.GT.I) GO TO 70
                  IF (J.EQ.I) THEN
                      CALL FCHAIN(I,1,XI,XIDOT,FIRR,FD)
                      GO TO 70
                  END IF
   65         CONTINUE
          END IF
      END IF
*
*       Check option for external tidal field.
*       Use predicted force of previous step for compatibility with regint.
   70 DT = TIME - T0(I)
      IF (KZ(14).GT.0) THEN
          DO 75 K = 1,3
              FREG(K) = FR(K,I) + DT*FRDOT(K,I)
   75     CONTINUE
          CALL XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDUM,0)
      END IF
*
*       Include the corrector and set new F, FDOT, D1, D2 & D3.
      DTSQ = DT**2
      DT6 = 6.0/(DT*DTSQ)
      DT2 = 2.0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DT13 = ONE3*DT
*
*     --03/03/14 20:51-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(time.eq.0.96801757812500000.and.rank.eq.0) then
c$$$         print*,rank,'I',I,'n',name(i),'fi',fi(1,i),
c$$$     &        'fin',firr(1),'s',step(i),'fd',fidot(1,i),
c$$$     &        'fdn',fd(1),'nb',list(1,i)
c$$$         call flush(6)
c$$$      end if
*     --03/03/14 20:51-lwang-end----------------------------------------*
      DO 80 K = 1,3
         DF = FI(K,I) - FIRR(K)
         FID = FIDOT(K,I)
         SUM = FID + FD(K)
         AT3 = 2.0*DF + DT*SUM
         BT2 = -3.0*DF - DT*(SUM + FID)
*       Use here new variables for consistency in parallel execution (R.Sp.)
          XN(K,I) = XI(K) + (0.6*AT3 + BT2)*DTSQ12
          XNDOT(K,I) = XIDOT(K) + (0.75*AT3 + BT2)*DT13
*
          FI(K,I) = FIRR(K)
          FIDOT(K,I) = FD(K)
*       Use total force for irregular step (cf. Makino & Aarseth PASJ, 1992).
          FDUM(K) = FIRR(K) + FR(K,I)
*
          D2(K,I) = (3.0*AT3 + BT2)*DT2
          D3(K,I) = AT3*DT6
*       NOTE: These are real derivatives!
   80 continue
*
      TTMP = TSTEP(FDUM,FD,D2(1,I),D3(1,I),ETAI)
      DT0 = TTMP
*     --03/03/14 20:29-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(rank.eq.0.and.time.ge.4.1667480468750000.and.name(i).eq.1155)
c$$$     &     then
c$$$         print*,rank,'nbint i',i,'n',name(i),'firr',firr(1),'fd',fd(1),
c$$$     &        'step',step(i),'dt0',dt0,'t',time,'x',xi(1),'xd',xidot(1),
c$$$     *        'nb',list(1,i),'step',step(i)
c$$$         call flush(6)
c$$$      end if
      if(dt0.le.0.1*step(i)) then
         write(6,81) rank,I,name(i),DT0/STEP(I),dt0,step(i),stepr(i),
     &        FI(1,i),FIDOT(1,i),D2(1,i),D3(1,i),time,t0(i),t0r(i),
     &        LIST(1,I),LIST(2,I),NAME(LIST(2,I))
         call flush(6)
 81      format(I3,' Warning!: Irregular step jumping! I',I7,' N',I7,
     &        ' ratio',E10.3,' dt0',F20.17,' step',F20.17,
     &        ' stepr',F20.17,' FI',E12.5,
     &        ' FD',E12.5,' D2',E12.5,' D3',E12.5,' t',F21.17,
     &        ' t0',F21.17,' t0r',F21.17,
     &        ' NB',I4,' LIST1',I7,' N1',I7)
c$$$         if(name(i).eq.1155) stop
c$$$         print*,'N',N,'BODY',BODY(I)
c$$$         do k =2,LIST(1,I)+1
c$$$            kk =LIST(K,I)
c$$$            print*, 'I',KK,'N',NAME(KK),'STEP',STEP(KK),'M',BODY(KK)
c$$$         end do
      end if
*     --03/03/14 20:29-lwang-end----------------------------------------*
*
*       Suggestion due to Winston Sweatman
*     DVV = (XDOT(1,I)-X0DOT(1,I))**2 + (XDOT(2,I)-X0DOT(2,I))**2 +
*    &     (XDOT(3,I)-X0DOT(3,I))**2
*     FFD = FIRR(1)**2 + FIRR(2)**2 + FIRR(3)**2
*     ETAIW = ETAI
*     TTMPW = ETAIW*DVV*BODY(I)/FFD
*
*     PRINT*,' irr I=',I,' TTMP,TTMPW,RATIO=',
*    &  TTMP,TTMPW,TTMP/TTMPW
*
*     IF(TTMP.GT.TTMPW)THEN
*     IGT = IGT + 1
*     ELSE
*     ILE = ILE + 1
*     END IF
*     IF(MOD(IGT+ILE,100).EQ.0)PRINT*,' irr IGT,ILE=',IGT,ILE
*
*     TTMP = MAX(TTMPW,TTMP)
*     DT0 = TTMP
*
*     IF (I.GT.N) THEN
*       Check for hierarchical configuration but exclude small perturbations.
*         IF (H(IPAIR).LT.-ECLOSE.AND.KZ(36).GT.0) THEN
*             IF (GAMMA(IPAIR).GT.1.0E-04) THEN
*                 CALL KEPLER(I,TTMP)
*                 DT0 = TTMP
*             END IF
*         END IF
*     END IF
*
*       Include convergence test for large step (cf. Makino, Ap.J. 369, 200).
      IF (TTMP.GT.STEPJ.AND.N.GT.1000) THEN
         DV2 = 0.0
         F2 = 0.0
*       Use only low order predicted value here.
         DO 85 K = 1,3
            DV2 = DV2 + (XIDOT(K) - XNDOT(K,I))**2
            F2 = F2 + FIRR(K)**2
   85    CONTINUE
*       Employ Makino criterion to avoid over-shooting (cf. Book, 2.16).
         DTJ = STEP(I)*(1.0D-06*STEP(I)**2*F2/DV2)**0.1
         TTMP = MIN(TTMP,DTJ)
      END IF
*
*       Select discrete value (increased by 2, decreased by 2 or unchanged).
      IF (TTMP.GT.2.0*STEP(I)) THEN
          IF (DMOD(TIME,2.0*STEP(I)).EQ.0.0D0) THEN
              TTMP = MIN(2.0*STEP(I),SMAX)
          ELSE
              TTMP = STEP(I)
          END IF
      ELSE IF (TTMP.LT.STEP(I)) THEN
          TTMP = 0.5*STEP(I)
            IF (TTMP.GT.DT0) THEN
                TTMP = 0.5*TTMP
          END IF
      ELSE
          TTMP = STEP(I)
      END IF
*
      STEP(I) = TTMP
*
*     --03/30/13 22:34-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if (time.ge.0.3.and.name(i).eq.6193) then
c$$$      j=i
c$$$         write(105+rank,*),'nb',j,name(j),'x0',(x0(kk,j),kk=1,3),
c$$$     *     'x0dot',(x0dot(kk,j),kk=1,3),'xn',(xn(kk,j),kk=1,3),
c$$$     *     'xndot',(xndot(kk,j),kk=1,3),
c$$$     *     't0',t0(j),'step',step(j),
c$$$     *        'stepr',stepr(j),
c$$$     *        'f',(f(kk,j),kk=1,3),'fdot',(fdot(kk,j),kk=1,3),
c$$$     *       'fi',(fi(kk,j),kk=1,3),'fidot',(fidot(kk,j),kk=1,3),
c$$$     *        'd0',(d0(kk,j),kk=1,3),'d1',(d1(kk,j),kk=1,3),
c$$$     *        'd2',(d2(kk,j),kk=1,3),'d3',(d3(kk,j),kk=1,3),
c$$$     *        'd0r',(d0r(kk,j),kk=1,3),'d1r',(d1r(kk,j),kk=1,3),
c$$$     *        'd2r',(d2r(kk,j),kk=1,3),'d3r',(d3r(kk,j),kk=1,3),
c$$$     *        'body',body(j),'time',time,'list',list(1,j)
c$$$*      if (rank.eq.0.and.nstepr.gt.4314800) then
c$$$*         write (6,*) name(i), nnb,time, freg(1)
c$$$* 118     format ('name, nnb, freg, time', 2I6, 1F2.6, 1E10.2)
c$$$         call flush(105+rank)
c$$$      end if
*     --03/30/13 22:34-lwang-end----------------------------------------*
      
      RETURN
*
      END
