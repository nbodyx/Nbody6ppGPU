      SUBROUTINE TTFORCE(XI,XIDOT,FM,FD,DT)
*
*
*       Compute the galactic force and its time derivative (mode B)
*       -----------------------------------------------------------
*
*** FlorentR - new subroutine

      INCLUDE 'common6.h'
      INCLUDE 'tt.h'
      INCLUDE 'galaxy.h'

      INTEGER TTTDEP
      REAL*8 XI(3),XIDOT(3),FM(3),FD(3), DT, TTTMP

      REAL*8 TTX(75), TTDX(4), TT2DX(4), TT12DX(4)
      REAL*8 TTPHI(25)

      SAVE TTDX

      TTDX(1) = 4d-5 * SQRT(XI(1)**2+XI(2)**2+XI(3)**2)
      IF(TTDX(1) .EQ. 0) TTDX(1) = 1.0
      TTDX(2) = TTDX(1)
      TTDX(3) = TTDX(1)

* update the timestep only for the guiding centre. Else, use the stored
* value
      IF(DT .NE. 0.0) THEN
        TTDX(4) = DT
      ENDIF

      DO K=1,4
        TT2DX(K) = 2.0 * TTDX(K)
        TT12DX(K) = 12.0 * TTDX(K)
      ENDDO
      
* Build the stencil using the value of TTDX
* TTX contains the position of the 25 points, relative to XI
* The position XI is number 1
* TTX(1:25)  = x positions
* TTX(26:50) = y positions
* TTX(51:75) = z positions

      TTX = (/ 0, 
     &        0, 0,-1, 1, 0, 0,
     &        0, 0,-1, 1, 0, 0,-1, 1,-2, 2,-1, 1, 0, 0,-1, 1, 0, 0,
     &        0, 
     &        0,-1, 0, 0, 1, 0,
     &        0,-1, 0, 0, 1,-2,-1,-1, 0, 0, 1, 1, 2,-1, 0, 0, 1, 0,
     &        0,
     &        -1, 0, 0, 0, 0, 1,
     &        -2,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2 /)
        
      DO I=1,25
        TTX(I) = TTX(I)*TTDX(1)
        TTX(I+25) = TTX(I+25)*TTDX(2)
        TTX(I+50) = TTX(I+50)*TTDX(3)
      END DO

* Get the galactic potential around the position XI
      DO I=1,25
         CALL TTGALAXY(XI(1)+TTX(I), XI(2)+TTX(I+25), XI(3)+TTX(I+50),
     &        TG, ZMBAR, RBAR, TSTAR, VSTAR, TTTMP, TTTDEP, IKEYPOT)
         TTPHI(I) = TTTMP
      END DO

* Get the galactic force at XI:

      FM(1) = (TTPHI(17)-8.*TTPHI(5)+8.*TTPHI(4)-TTPHI(16))/TT12DX(1)
      FM(2) = (TTPHI(20)-8.*TTPHI(6)+8.*TTPHI(3)-TTPHI(13))/TT12DX(2)
      FM(3) = (TTPHI(25)-8.*TTPHI(7)+8.*TTPHI(2)-TTPHI(8))/TT12DX(3)

* Get the tidal tensor at XI:

      TTENS(1,1,1) = (TTPHI(1)-TTPHI(17)-TTPHI(16)+TTPHI(1)) 
     &        / (TT2DX(1)*TT2DX(1))
      TTENS(1,2,1) = (TTPHI(18)-TTPHI(19)-TTPHI(14)+TTPHI(15))
     &        / (TT2DX(1)*TT2DX(2))
      TTENS(1,3,1) = (TTPHI(22)-TTPHI(23)-TTPHI(10)+TTPHI(11))
     &        / (TT2DX(1)*TT2DX(3))
      TTENS(2,2,1) = (TTPHI(1)-TTPHI(20)-TTPHI(13)+TTPHI(1))
     &        / (TT2DX(2)*TT2DX(2))
      TTENS(2,3,1) = (TTPHI(21)-TTPHI(24)-TTPHI(9)+TTPHI(12))
     &        / (TT2DX(2)*TT2DX(3))
      TTENS(3,3,1) = (TTPHI(1)-TTPHI(25)-TTPHI(8)+TTPHI(1))
     &        / (TT2DX(3)*TT2DX(3))
      TTENS(2,1,1) = TTENS(1,2,1)
      TTENS(3,1,1) = TTENS(1,3,1)
      TTENS(3,2,1) = TTENS(2,3,1)


      FD(1) = TTENS(1,1,1) * XIDOT(1) + TTENS(1,2,1) * XIDOT(2) 
     &   + TTENS(1,3,1) * XIDOT(3)
      FD(2) = TTENS(2,1,1) * XIDOT(1) + TTENS(2,2,1) * XIDOT(2)
     &   + TTENS(2,3,1) * XIDOT(3)
      FD(3) = TTENS(3,1,1) * XIDOT(1) + TTENS(3,2,1) * XIDOT(2)
     &   + TTENS(3,3,1) * XIDOT(3)

* if the potential is time dependent
      IF(TTTDEP .GT. 0 .AND. TTDX(4) .NE. 0.0) THEN

* compute partial F / partial t at t+dt and t-dt.
        DO I=2,7
           CALL TTGALAXY(XI(1)+TTX(I), XI(2)+TTX(I+25), XI(3)+TTX(I+50),
     &          TG-TTDX(4), ZMBAR, RBAR, TSTAR, VSTAR, TTTMP, TTTDEP, 
     &          IKEYPOT)
           TTPHI(I) = TTTMP

           CALL TTGALAXY(XI(1)+TTX(I), XI(2)+TTX(I+25), XI(3)+TTX(I+50),
     &        TG+TTDX(4), ZMBAR, RBAR, TSTAR, VSTAR, TTTMP, TTTDEP,
     &          IKEYPOT)
           TTPHI(I+6) = TTTMP
        END DO

* sum T*v and partial F / partial t
        FD(1) = FD(1)+(-TTPHI(4)+TTPHI(5)+TTPHI(10)-TTPHI(11)) 
     &        / (TT2DX(4)*TT2DX(1))
        FD(2) = FD(2)+(-TTPHI(3)+TTPHI(6)+TTPHI(9)-TTPHI(12))
     &        / (TT2DX(4)*TT2DX(2))
        FD(3) = FD(3)+(-TTPHI(2)+TTPHI(7)+TTPHI(8)-TTPHI(13))
     &        / (TT2DX(4)*TT2DX(3))
      ENDIF

      RETURN
      END

* Positions of the TTX points
*              
*                     Z ^
*                       |    ^  Y
*                       |   /
*                       |  /
*                       | / 
*                       +--------> X
*
*
*                      +--------+--------+--------+--------+  
*                     /        /        /        /        / 
*                    /        /        /        /        /
*                   +--------+--------+--------+--------+
*                  /        /        /        /        /
*                 /        /        /        /        /
*                +--------+------- 25 ------+--------+    z = +2 * TTDX
*               /        /        /        /        /
*              /        /        /        /        /
*             +--------+--------+--------+--------+
*            /        /        /        /        /
*           /        /        /        /        /
*          +--------+--------+--------+--------+
*          
*          
*                      +--------+--------+--------+--------+
*                     /        /        /        /        /
*                    /        /        /        /        /
*                   +--------+------- 24 ------+--------+
*                  /        /        /        /        /
*                 /        /        /        /        /
*                +------- 22 ----- 7 ------ 23 ------+    z = +1 * TTDX
*               /        /        /        /        /
*              /        /        /        /        /
*             +--------+------- 21 ------+--------+
*            /        /        /        /        /
*           /        /        /        /        /
*          +--------+--------+--------+--------+
*          
*          
*                      +--------+------- 20 ------+--------+
*                     /        /        /        /        /
*                    /        /        /        /        /
*                   +------- 18----- 6 ------ 19 ------+
*                  /        /        /        /        /
*                 /        /        /        /        /
*                16 ----- 4 ----- 1 ------- 5 ----- 17     z = 0
*               /        /        /        /        /
*              /        /        /        /        /
*             +------ 14 ------ 3 ------ 15 ------+
*            /        /        /        /        /
*           /        /        /        /        /
*          +--------+------ 13 -------+--------+
*          
*          
*                      +--------+--------+--------+--------+
*                     /        /        /        /        /
*                    /        /        /        /        /
*                   +--------+------ 12 -------+--------+
*                  /        /        /        /        /
*                 /        /        /        /        /
*                +------ 10 ------ 2 ------ 11 ------+     z = -1 * TTDX
*               /        /        /        /        /
*              /        /        /        /        /
*             +--------+------- 9 -------+--------+
*            /        /        /        /        /
*           /        /        /        /        /
*          +--------+--------+--------+--------+
*          
*          
*                      +--------+--------+--------+--------+
*                     /        /        /        /        /
*                    /        /        /        /        /
*                   +--------+--------+--------+--------+
*                  /        /        /        /        /
*                 /        /        /        /        /
*                +--------+------- 8 -------+--------+     z = -2 * TTDX
*               /        /        /        /        /
*              /        /        /        /        /
*             +--------+--------+--------+--------+
*            /        /        /        /        /
*           /        /        /        /        /
*          +--------+--------+--------+--------+
*
