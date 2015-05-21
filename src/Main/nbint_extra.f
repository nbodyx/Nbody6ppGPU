      SUBROUTINE NBINT_EXTRA(I,FIRR,FD)
*
*
*     C.M. correction and extra force to irregular force
*     ------------------------
      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
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
*
      RETURN
*
      END
