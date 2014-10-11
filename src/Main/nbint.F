      SUBROUTINE NBINT(I)
*
*
*       Irregular integration.
*       ----------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/NEWXV/ XN(3,NMAX),XNDOT(3,NMAX)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FREG(3),FD(3),FDUM(3),DV(3)
*      LOGICAL RPRED(NMAX)
*
*       Obtain irregular force & first derivative.
*      if (.not.NOPRED(I)) RPRED(I)=.false.
*      if (NOPRED(I).and.RPRED(I)) print*,rank,K,'not predict',time
*       NOPRED(I) = .true.
      call jpred_int(I,TIME)
      DO 5 K = 1,3
          XI(K) = X(K,I)
          XIDOT(K) = XDOT(K,I)
          FIRR(K) = 0.0D0
          FD(K) = 0.0D0
    5 CONTINUE
*          PHII(I) = 0.D0
          NNB0 = LIST(1,I)
!$omp critical      
          NIRRF = NIRRF + NNB0
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
#ifdef GPU
*     --10/15/13 19:05-lwang-cm-body------------------------------------*
***** Note: Make calculation of irr for pertuber cm consistent with regint**
      DO  L = 2,NNB0+1
         K = LIST(L,I)
         call jpred_int(K,TIME)
         A1 = X(1,K) - XI(1)
         A2 = X(2,K) - XI(2)
         A3 = X(3,K) - XI(3)
         DV(1) = XDOT(1,K) - XIDOT(1)
         DV(2) = XDOT(2,K) - XIDOT(2)
         DV(3) = XDOT(3,K) - XIDOT(3)
         RIJ2 = A1*A1 + A2*A2 + A3*A3
*     
         DR2I = 1.0/RIJ2
         DR3I = BODY(K)*DR2I*SQRT(DR2I)
         DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
         FIRR(1) = FIRR(1) + A1*DR3I
         FIRR(2) = FIRR(2) + A2*DR3I
         FIRR(3) = FIRR(3) + A3*DR3I
*     PHII(I) = PHII(I) - DR3I*RIJ2
         FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
         FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
         FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
*     --03/30/13 22:34-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         IF(FIRR(1).NE.FIRR(1)) then
c$$$            print*,'XI',XI(1), 'X',X(1,K),'XD',XDOT(1,K)
c$$$     &           ,'RIJ2',RIJ2,'K',K,'N',NAME(K),
c$$$     &           'I',I,'N',NAME(I)
c$$$            call flush(6)
c$$$            STOP
c$$$         end if
*         if(time.eq.3.53173828125000000.and.name(i).eq.1260) then
c$$$         if(rank.eq.0) then
c$$$         if(time.ge.4.1667480468750000.and.name(i).eq.1155) then
c$$$            if(name(k).eq.4174) then
c$$$            write(6,*),rank,'I',I,'K',K,'N',name(K),'FI',FIRR,'FD',FD,
c$$$     &           'RIJ2',RIJ2,'XDI',XIDOT,'XDJ',XDOT(1:3,K),'MJ',body(k),
c$$$     &              'T',TIME,'T0J',T0(K)
c$$$            jj=2*(k-n)-1
c$$$            print*,'x1',x(1:3,jj),'x2',x(1:3,jj+1),
c$$$     &           't',time,'xi',xi(1:3),'m',body(jj),body(jj+1),body(i)
c$$$            do JJ=2,LIST(1,K)
c$$$               JK= LIST(JJ,K)
c$$$               if (name(JK).eq.1417) print*,'1 NB',LIST(1,K),' JJ',JJ
c$$$            end do
c$$$            end if
c$$$            call flush(6)
c$$$         end if
c$$$         end if
c$$$*     --03/30/13 22:34-lwang-end----------------------------------------*
      END DO
#endif
      
*       Choose force loop for single particle or regularized c.m. body.
      IF (I.LE.N) GO TO 20
*
*       Set KS pair index.
      IPAIR = I - N
*       Adopt c.m. approximation for small total perturbation.
      I1 = 2*IPAIR - 1
      IF (LIST(1,I1).GT.0) THEN
#ifdef GPU
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
#else         
*       Obtain irregular force on perturbed c.m. body (including any chain).
         CALL CMFIRR(I,IPAIR,XI,XIDOT,FIRR,FD)
#endif          
          GO TO 70
      END IF
*
#ifdef GPU
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
#else      
*       Set neighbour number & list index of the last single particle.
   20 NNB1 = NNB0 + 1
      NNB2 = NNB1
   25 IF (LIST(NNB2,I).LE.N) GO TO 30
      NNB2 = NNB2 - 1
      IF (NNB2.GT.1) GO TO 25
*       Include special case of only c.m. neighbours.
      GO TO 40
*
*       Sum over single particles (unperturbed case included).
 30   CONTINUE
**!$omp parallel do if(NNB2.GE.5*ITHREAD) default(shared)
**!$omp& private(L,K,A1,A2,A3,DV,RIJ2,DR2I,DR3I,DRDV)
**!$omp& reduction(+:FIRR,FD,PHIT)
      DO 35 L = 2,NNB2
          K = LIST(L,I)
          call jpred_int(K,TIME)
*          call jpred(K,time,time)
*         if (.not.NOPRED(K)) RPRED(K)=.false.
*         if (NOPRED(K).and.RPRED(K)) print*,rank,K,'not predict',time
*           NOPRED(K) = .true.
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          DV(1) = XDOT(1,K) - XIDOT(1)
          DV(2) = XDOT(2,K) - XIDOT(2)
          DV(3) = XDOT(3,K) - XIDOT(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(K)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
*          PHII(I) = PHII(I) - DR3I*RIJ2
          FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
*     --10/14/13 23:46-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          if (time.ge.3.9062000000000000E-003.and.name(i).eq.1689) then
c$$$             print*,k,name(k),rij2,ifirst
c$$$             call flush(6)
c$$$          end if
*     --10/14/13 23:46-lwang-end----------------------------------------*

   35 CONTINUE
**!$omp end parallel do
*
*       See whether any c.m. neighbours should be included.
      IF (NNB2.EQ.NNB1) GO TO 60
*
   40 NNB3 = NNB2 + 1
*       Set index for distinguishing c.m. or resolved components.
      KDUM = 0
*
*       Sum over regularized c.m. neighbours.
**!$omp parallel do IF(NNB1-NNB3.GE.5*ITHREAD) default(shared)
**!$omp& private(L,K,A1,A2,A3,DV,RIJ2,DR2I,DR3I,DRDV,J,KDUM)
**!$omp& reduction(+:FIRR,FD)
      DO 50 L = NNB3,NNB1
          K = LIST(L,I)
          call jpred_int(K,TIME)
*          call jpred(k,time,time)
*          if (.not.NOPRED(K)) RPRED(K)=.false.
*          if (NOPRED(K).and.RPRED(K)) print*,rank,K,'not predict',time
*          NOPRED(K) = .true.
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          DV(1) = XDOT(1,K) - XIDOT(1)
          DV(2) = XDOT(2,K) - XIDOT(2)
          DV(3) = XDOT(3,K) - XIDOT(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
*       See whether c.m. approximation applies (ignore unperturbed case).
          J = K - N
          KDUM = 2*J - 1
          IF (RIJ2.GT.CMSEP2*R(J)**2.OR.LIST(1,KDUM).EQ.0) GO TO 48
*
          K = KDUM
*       Sum over individual components of pair #J.
   45     A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          DV(1) = XDOT(1,K) - XIDOT(1)
          DV(2) = XDOT(2,K) - XIDOT(2)
          DV(3) = XDOT(3,K) - XIDOT(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
*       Adopt c.m. approximation outside the effective perturber sphere.
   48     DR2I = 1.0/RIJ2
          DR3I = BODY(K)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
*     --10/14/13 23:46-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          if (time.ge.3.9062000000000000E-003.and.name(i).eq.1689) then
c$$$             print*,k,name(k),j+n,name(j+n),rij2,ifirst,x(1,k),x(1,j+n)
c$$$     *            ,time,tprev
c$$$             call flush(6)
c$$$          end if
*     --10/14/13 23:46-lwang-end----------------------------------------*
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
          FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
          IF (K.EQ.KDUM) THEN
              K = K + 1
              GO TO 45
          END IF
   50 CONTINUE
**!$omp end parallel do
 60   CONTINUE
#endif

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
