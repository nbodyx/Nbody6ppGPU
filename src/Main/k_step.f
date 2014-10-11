      function k_step(dt,dtk)
*
*     
*     Get level of step
*
      INTEGER K_STEP
      REAL*8 DT,DTK(64),DT1

*     Check ghost particle
      IF (DT.GT.DTK(1)) THEN
*     Indicator ghost particle
         K_STEP = -1
         RETURN
      END IF
*
      K_STEP = 1

*       Compare predicted step with discrete values decreasing by 1/32.
      DT1 = DTK(1)
    1 IF (DT1.GT.DT) THEN
          DT1 = 0.03125D0*DT1
          K_STEP = K_STEP + 5
          GO TO 1
      END IF
*
*       Increase by 2 until predicted step is exceeded (then correct level).
    4 IF (DT1.LT.DT) THEN
          DT1 = 2.0D0*DT1
          K_STEP = K_STEP - 1
          GO TO 4
      END IF

      RETURN
*
      END

