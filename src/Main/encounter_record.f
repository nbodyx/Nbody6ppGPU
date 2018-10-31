      SUBROUTINE ENCOUNTER_RECORD(NXTLEN,NXTLST)
*     label-encounter
*
*     close encounter record
*
*
      INCLUDE 'common6.h'
      include 'omp_lib.h'
      
      COMMON/ENCOUNTER/ IMINR_PRE(NMAX) 
      REAL*8 DX(3), DV(3), DXV
      INTEGER NXTLEN, NXTLST(NMAX)
      INTEGER(4) IUNIT

!$OMP PARALLEL PRIVATE(II,I,J,DX,DV,DXV,IUNIT)
      IUNIT = omp_get_thread_num()
      IUNIT = IUNIT + 100
!$OMP DO
      DO II = 1, NXTLEN
         I = NXTLST(II)
         IF (LIST(1,I).GT.0) THEN
            J = IMINR(I)
            CALL JPRED_int(J,TIME)

            IF(J.LE.0) call abort
            
            DX(1:3) = X(1:3,J) - X0(1:3,I)
            DV(1:3) = XDOT(1:3,J) - X0DOT(1:3,I)
            
            DXV = DX(1)*DV(1) + DX(2)*DV(2) + DX(3)*DV(3)

            IF (IMINR_PRE(I).EQ.J.AND.DXV.GT.0) THEN
               if(rank.eq.0) then
                  call HYPERBOLIC_RECORD(NAME(I),NAME(J),
     &                 BODY(I),BODY(J),X0(1,I),X0DOT(1,I),
     &                 DX,DV,DXV,IUNIT,TIME)
               end if
               IMINR_PRE(I) = -2
            ELSE
               IF(DXV.LT.0) IMINR_PRE(I) = J
            END IF 
         END IF 
      END DO
!$OMP END DO
!$OMP END PARALLEL

      RETURN

      END

*-----------------------------------------------

      SUBROUTINE ENCOUNTER_INIT

      include 'params.h'
      include 'omp_lib.h'

      COMMON/ENCOUNTER/ IMINR_PRE(NMAX)
      INTEGER IMINR_PRE

      IMINR_PRE(1:NMAX) = -10

      RETURN

      END

*-----------------------------------------------


      SUBROUTINE HYPERBOLIC_RECORD(I,J,M1,M2,X1,V1,DX,DV,DXV,OUNIT,TIME)

      implicit none
      INTEGER I,J,OUNIT
      REAL*8 M1,M2,X1(3),V1(3),DX(3),DV(3),DXV,TIME
      REAL*8 E,EBIN,XR,XR2,VR2,M12,SEMI,ECC,PERI,VINF,VP,B
      
      XR2 = DX(1)*DX(1) + DX(2)*DX(2) + DX(3)*DX(3)
      XR = SQRT(XR2)
      VR2 = DV(1)*DV(1) + DV(2)*DV(2) + DV(3)*DV(3)
      M12 = M1 + M2
      
      E = VR2/M12 - 2.0/XR
      SEMI = -1.0/E
      EBIN = M1*M2*E/2.0;

      ECC = SQRT((1-XR/SEMI)**2 + DXV*DXV/SEMI/M12)
      PERI = SEMI*(1-ECC)

      VINF = SQRT(ABS(M12/SEMI))
      VP = VINF*SQRT((ECC+1)/ABS(ECC-1))

      B = -SEMI*SQRT(ABS(ECC*ECC-1))

      WRITE(OUNIT,1) TIME,I,J,M1,M2,X1(1:3),V1(1:3),
     &     SEMI,ECC,PERI,VINF,VP,B,EBIN
 1    FORMAT(1P,E26.14,0P,2I8,1P,15E26.14,0P)

      RETURN

      END
