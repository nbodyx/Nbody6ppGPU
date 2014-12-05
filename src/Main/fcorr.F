      SUBROUTINE FCORR(I,DM,KW)
*
*
*       Total force corrections due to masss loss.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      REAL*8  A(6),VOLD(3)
      LOGICAL IKICK
*      REAL*8  XBACK(3,LMAX),XDBACK(3,LMAX)
*
*
*       Save the velocity components and square velocity.
      call xbpredall
*      call jpred(I,TIME,TIME)
      VI2 = 0.0
      DO 1 K = 1,3
          VOLD(K) = XDOT(K,I)
          VI2 = VI2 + XDOT(K,I)**2
    1 CONTINUE
*
*       Include velocity kick in case of new WD, NS, BH or massless SN.
      IF (KW.GE.10.AND.KW.LE.15) THEN
          IKICK = .TRUE.
*       Distinguish between single star (first time only) and binary.
          IF (I.LE.N.AND.KW.NE.KSTAR(I)) THEN
              CALL KICK(I,1,KW,DM)
          ELSE IF (I.GT.N) THEN
              IPAIR = I - N
              CALL KICK(IPAIR,0,KW,DM)
          END IF
      ELSE
          IKICK = .FALSE.
      END IF
*
*       Define consistent c.m. variables for KS mass loss (exclude Roche).
      VFAC = 0.0
      IF (I.GT.N.AND.KSTAR(I).LE.10) THEN
          I2 = 2*(I - N)
          I1 = I2 - 1
          VF2 = 0.0
          DV2 = 0.0
          BODYI = BODY(I)
          IF (BODY(I).EQ.0.0D0) BODYI = BODY(I1) + BODY(I2)
          DO 8 K = 1,3
              X(K,I) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/BODYI
              XDOT(K,I) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/
     &                                                           BODYI
              X0(K,I) = X(K,I)
              X0DOT(K,I) = XDOT(K,I)
              VF2 = VF2 + XDOT(K,I)**2
              DV2 = DV2 + (XDOT(K,I) - VOLD(K))**2
    8     CONTINUE
          VFAC = SQRT(VF2/VI2)
      END IF
*
*       Correct potential energy, all forces & first derivatives.
      POTJ = 0.0D0
c$$$*       Backup all predicted x and xdot to xback and xdback
c$$$      NNB = LIST(1,I)
c$$$      DO II = 1, NNB
c$$$         IL = LIST(II+1,I)
c$$$         CALL JPRED(IL,TIME,TIME)
c$$$         XBACK(1,II) = X0(1,IL)
c$$$         XBACK(2,II) = X0(2,IL)
c$$$         XBACK(3,II) = X0(3,IL)
c$$$         X0(1,IL) = X(1,IL)
c$$$         X0(2,IL) = X(2,IL)         
c$$$         X0(3,IL) = X(3,IL)
c$$$         XDBACK(1,II) = X0DOT(1,IL)
c$$$         XDBACK(2,II) = X0DOT(2,IL)
c$$$         XDBACK(3,II) = X0DOT(3,IL)
c$$$         X0DOT(1,IL) = XDOT(1,IL)
c$$$         X0DOT(2,IL) = XDOT(2,IL)
c$$$         X0DOT(3,IL) = XDOT(3,IL)
c$$$      END DO

!$omp parallel do 
!$omp& private(J,NNB1,RIJ2,RIJDOT,DX,DXDOT,RIJ,
!$omp& RDVDOT,A3,A4,A5,A6,A7,A,K)
!$omp& reduction(+:POTJ)
      DO 40 J = IFIRST,NTOT
#ifdef SIMD
          IF (J.EQ.I) GO TO 35
#else
          IF (J.EQ.I) GO TO 40
#endif
          RIJ2 = 0.0D0
          RIJDOT = 0.0D0
          RDVDOT = 0.0D0
*
          DO 10 K = 1,3
              A(K) = X(K,I) - X(K,J)
              A(K+3) = VOLD(K) - XDOT(K,J)
              RIJ2 = RIJ2 + A(K)**2
              RIJDOT = RIJDOT + A(K)*A(K+3)
              RDVDOT = RDVDOT + A(K)*(XDOT(K,I) - VOLD(K))
   10     CONTINUE
*
          RIJ = SQRT(RIJ2)
          POTJ = POTJ + BODY(J)/RIJ
          A3 = 1.0/(RIJ2*RIJ)
          A4 = BODY(I)*A3
          A5 = DM*A3
          A6 = 3.0*RIJDOT/RIJ2
          A7 = 3.0*RDVDOT/RIJ2
*
          DO 15 K = 1,3
              A(K+3) = (A(K+3) - A6*A(K))*A5
              IF (IKICK) THEN
*       Include FDOT corrections due to increased velocity.
                  A(K+3) = A(K+3) + (XDOT(K,I) - VOLD(K))*A4
                  A(K+3) = A(K+3) - A7*A(K)*A4
              END IF
   15     CONTINUE
*
*       Use neighbour list of #J to distinguish irregular & regular terms.
          NNB1 = LIST(1,J) + 1
          DO 25 L = 2,NNB1
              IF (LIST(L,J).EQ.I) THEN
                  DO 20 K = 1,3
                      F(K,J) = F(K,J) - 0.5*A(K)*A5
                      FI(K,J) = FI(K,J) - A(K)*A5
                      FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
                      D1(K,J) = D1(K,J) - A(K+3)
                      FIDOT(K,J) = FIDOT(K,J) - A(K+3)
   20             CONTINUE
c$$$                  IF(IKICK.AND.T0(J)+0.5D0*STEP(J).GE.TBLOCK) THEN
c$$$                     STEP(J) = 0.5D0*STEP(J)
c$$$                     IF(T0(J)+0.5D0*STEP(J).GE.TBLOCK) THEN
c$$$                        STEP(J) = 0.5D0*STEP(J)
c$$$                     END IF
c$$$                     TIMENW(J) = T0(J) + STEP(J)
c$$$                  END IF

*     --03/19/14 22:48-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$          if(name(J).eq.7) then
c$$$             print*,rank,'J',J,'F',F(1,J),'FD',FDOT(1,J),
c$$$     &            'FI',FI(1,J),'FID',D1(1,J),'FC',0.5*A(1)*A5,'FDC',
c$$$     &            -ONE6*A(4),
c$$$     &            't',time,'XDOT',XDOT(1,J)
c$$$             call flush(6)
c$$$          end if
*     --03/19/14 22:48-lwang-end----------------------------------------*
#ifdef SIMD
                  GO TO 35
#else
                  GO TO 40
#endif
              END IF
   25     CONTINUE      
          DO 30 K = 1,3
              F(K,J) = F(K,J) - 0.5*A(K)*A5
              FR(K,J) = FR(K,J) - A(K)*A5
              FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
              D1R(K,J) = D1R(K,J) - A(K+3)
              FRDOT(K,J) = FRDOT(K,J) - A(K+3)
   30     CONTINUE
#ifdef SIMD
*     Reset AVX/SSE particle data
 35       call IRR_SIMD_SET_JP(J,X0(1,J),X0DOT(1,J),F(1,J),FDOT(1,J),
     &         BODY(J),T0(J))
#endif
   40 CONTINUE
!$omp end parallel do

c$$$*       Revert X0 and X0dot to original values
c$$$      DO II = 1, NNB
c$$$         IL = LIST(II+1,I)
c$$$         X0(1,IL) = XBACK(1,II)
c$$$         X0(2,IL) = XBACK(2,II)
c$$$         X0(3,IL) = XBACK(3,II)
c$$$         X0DOT(1,IL) = XDBACK(1,II)
c$$$         X0DOT(2,IL) = XDBACK(2,II)
c$$$         X0DOT(3,IL) = XDBACK(3,II)
c$$$      END DO
*
*       Update the potential and kinetic energy loss.
      EMDOT = EMDOT - DM*POTJ + 0.5*DM*VI2
*
*       Modify energy loss further for c.m. body (exclude Roche cases).
      IF (I.GT.N.AND.KSTAR(I).LE.10) THEN
          ECDOT = ECDOT - 0.5*BODY(I)*VI2*(VFAC**2 - 1.0)
      END IF
*
*       See whether linearized tidal terms should be included.
      IF (KZ(14).GT.0.AND.KZ(14).LT.3) THEN
          EMDOT = EMDOT - 0.5*DM*(TIDAL(1)*X(1,I)**2 +
     &                            TIDAL(3)*X(3,I)**2)
      END IF
*
*       Check optional Plummer potential.
      IF (KZ(14).EQ.4.OR.KZ(14).EQ.3) THEN
          RI2 = AP2
          DO 50 K = 1,3
              RI2 = RI2 + X(K,I)**2
  50      CONTINUE
          EMDOT = EMDOT - DM*MP/SQRT(RI2)
      END IF
*
*       Accumulate energy loss for conservation check (not used).
      E(12) = EMDOT
*
      RETURN
*
      END
