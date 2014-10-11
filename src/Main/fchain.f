      SUBROUTINE FCHAIN(I,IR,XI,XIDOT,FIRR,FD)
*
*
*       Irregular force & derivative due to chain.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
*     Safe for parallel
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  XI(3),XIDOT(3),DX(3),DV(3),FIRR(3),FD(3),XIS(3),VIS(3)
*
*
*       Use c.m. values for correction of perturbed KS (call from KCPERT).
      I2 = 0
      IF (I.LT.IFIRST) THEN
          IPAIR = KVEC(I)
          IF (I.EQ.2*IPAIR) I2 = 1
          ICM = N + IPAIR
          call jpred_int(icm,time)
*       Save local variables for individual chain contributions.
          DO 1 K = 1,3
              XIS(K) = XI(K)
              VIS(K) = XIDOT(K)
              XI(K) = X(K,ICM)
              XIDOT(K) = XDOT(K,ICM)
    1     CONTINUE
      END IF
*
*       Evaluate terms for the original chain c.m. interaction.
      DR2 = 0.0
      DRDV = 0.0
      call jpred_int(ICH,TIME)
      DO 5 K = 1,3
          DX(K) = X(K,ICH) - XI(K)
          DV(K) = XDOT(K,ICH) - XIDOT(K)
          DR2 = DR2 + DX(K)**2
          DRDV = DRDV + DX(K)*DV(K)
    5 CONTINUE
      DR2I = 1.0/DR2
      DR3I = BODY(ICH)*DR2I*SQRT(DR2I)
      DRDV = 3.0*DRDV*DR2I
*
*       Subtract force and first derivative from current values.
      DO 10 K = 1,3
          FIRR(K) = FIRR(K) - DX(K)*DR3I
          FD(K) = FD(K) - (DV(K) - DX(K)*DRDV)*DR3I
   10 CONTINUE
*
*     --03/14/14 15:45-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(I.eq.1605) then
c$$$         write(125+rank,*)'fchain xi', xi(1),'xid',xidot(1),
c$$$     &        'xch',x(1,ich),'xdc',xdot(1,ich),'fd',fd(1),'t',time
c$$$         call flush(125+rank)
c$$$      end if
*     --03/14/14 15:45-lwang-end----------------------------------------*
*       Restore XI & XIDOT for KS components (perturbed case).
      IF (I.LT.IFIRST) THEN
          DO 12 K = 1,3
              XI(K) = XIS(K)
              XIDOT(K) = VIS(K)
   12     CONTINUE
      END IF
*
*       Resolve chain coordinates & velocities using the predicted c.m.
      IF (IR.EQ.0.AND.I2.EQ.0) THEN
          CALL XCPRED(0)
      END IF
*
*       Obtain contributions from all members of the chain.
      DO 20 J = 1,NCH
          DR2 = 0.0
          DRDV = 0.0
          DO 15 L = 1,3
              DX(L) = XC(L,J) - XI(L)
              DV(L) = UC(L,J) - XIDOT(L)
              DR2 = DR2 + DX(L)**2
              DRDV = DRDV + DX(L)*DV(L)
   15     CONTINUE
*     --03/14/14 15:45-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      if(I.eq.1605) then
c$$$         write(125+rank,*)'fchain j',j,'xc', xc(1,j),'uc',uc(1,j),'bc',
c$$$     &        bodyc(j),'t',time
c$$$         call flush(125+rank)
c$$$      end if
*     --03/14/14 15:45-lwang-end----------------------------------------*
          DR2I = 1.0/DR2
          DR3I = BODYC(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*DRDV*DR2I
*
*       Add force & first derivative.
          DO 18 L = 1,3
              FIRR(L) = FIRR(L) + DX(L)*DR3I
              FD(L) = FD(L) + (DV(L) - DX(L)*DRDV)*DR3I
   18     CONTINUE
   20 CONTINUE
*
      RETURN
*
      END
