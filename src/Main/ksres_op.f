      SUBROUTINE KSRES_OP(IPAIR,I1,I2,Q,RDOT,KCASE)
*
*
*       Coordinates & velocities of KS pair for output
*       ------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  Q(3),RDOT(3),UI(4),V(4),A1(3,4)
*
*
*       Note: KDUM = 1 for prediction of U & UDOT and = 2 for termination.
      I1 = 2*IPAIR - 1
      I2 = 2*IPAIR
*
*       Ensure appropriate prediction (second pair in binary collision).
      IF (LIST(1,I1).EQ.0.OR.T0(I1).EQ.TIME) THEN
*       Set current values of regularized coordinates & velocities.
          DO 1 K = 1,4
              UI(K) = U0(K,IPAIR)
              V(K) = UDOT(K,IPAIR)
    1     CONTINUE
*
*       Copy binding energy for routine ENERGY.
          HT = H(IPAIR)
          GO TO 4
      END IF
*
*       Predict U & UDOT to order FUDOT3 using third-order time inversion.
      A2 = 1.0/R(IPAIR)
      A3 = A2*(TIME - T0(I1))
*
*       See whether the time interval should be modified by KSLOW procedure.
      IF (KSLOW(IPAIR).GT.1) THEN
          IMOD = KSLOW(IPAIR)
          A3 = A3/FLOAT(ISLOW(IMOD))
      END IF
*
*       Expand regularized interval to third order.
      A4 = 3.0D0*TDOT2(IPAIR)**2*A2 - TDOT3(IPAIR)
      DTU = ((ONE6*A4*A3 - 0.5D0*TDOT2(IPAIR))*A2*A3 + 1.0)*A3
*       Apply safety test near small pericentre or for unperturbed motion.
      IF (DTU.GT.DTAU(IPAIR)) DTU = 0.8*DTAU(IPAIR)
      IF (DTU.LT.0.0) DTU = 0.0
*
      DTU1 = 0.2D0*DTU
      DTU2 = DTU/24.0D0
      DO 3 K = 1,4
          UI(K) = ((((FUDOT3(K,IPAIR)*DTU1 + FUDOT2(K,IPAIR))*DTU2 +
     &                          FUDOT(K,IPAIR))*DTU + FU(K,IPAIR))*DTU +
     &                                  UDOT(K,IPAIR))*DTU + U0(K,IPAIR)
          V(K) = (((FUDOT3(K,IPAIR)*DTU2 + ONE6*FUDOT2(K,IPAIR))*DTU +
     &              3.0D0*FUDOT(K,IPAIR))*DTU + 2.0D0*FU(K,IPAIR))*DTU +
     &                                                     UDOT(K,IPAIR)
    3 CONTINUE
*
*       Predict current binding energy per unit mass for routine ADJUST.
      HT = (((HDOT4(IPAIR)*DTU2 + ONE6*HDOT3(IPAIR))*DTU +
     &             0.5D0*HDOT2(IPAIR))*DTU + HDOT(IPAIR))*DTU + H(IPAIR)
*
*       Form current transformation matrix and two-body separation.
    4 CALL MATRIX(UI,A1)
      RI = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
*
*       Obtain relative coordinates & velocities.
      DO 7 J = 1,3
          Q(J) = 0.0D0
          RDOT(J) = 0.0D0
          DO 6 K = 1,4
              Q(J) = Q(J) + A1(J,K)*UI(K)
              RDOT(J) = RDOT(J) + 2.0D0*A1(J,K)*V(K)/RI
    6     CONTINUE
    7 CONTINUE
*
      IF(KCASE.EQ.1) THEN
         I = N + IPAIR

*     Set coordinates & velocities (XDOT for KDUM = 1) of KS components.
         DO 15 K = 1,3
            X(K,I1) = X(K,I) + BODY(I2)*Q(K)/BODY(I)
            X(K,I2) = X(K,I) - BODY(I1)*Q(K)/BODY(I)
            XDOT(K,I1) = XDOT(K,I) + BODY(I2)*RDOT(K)/BODY(I)
            XDOT(K,I2) = XDOT(K,I) - BODY(I1)*RDOT(K)/BODY(I)
 15      CONTINUE
      END IF
*
      RETURN
*
      END
