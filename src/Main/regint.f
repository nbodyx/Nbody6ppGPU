      SUBROUTINE REGINT(IOFF,NI,IREG,IFIRST,N,X,XDOT,BODY,
     &     DTR,BODYM,M_FLAG,NB_FLAG,NNBOPT,RS,FREG,FDR,LISTGP,POT)
*
*
*       Regular integration.
*       --------------------
*     M_FLAG: 1: mass weighted neighbor radius
*     N_FLAG: 1: reduce neighbor radius by (NNBOPT/NNB)**0.333 if NNB overflow
*     POT: positive potential value

      PARAMETER (maxthr=1024)
      include 'params.h'
      REAL*8 X(3,NMAX),XDOT(3,NMAX),BODY(NMAX),RS(NMAX),DTR(NMAX)
      REAL*8 FREG(3,maxthr),FDR(3,maxthr),POT(NMAX),BODYM
      INTEGER LISTGP(LMAX,maxthr),IREG(NMAX)
      INTEGER IOFF,NI,IFIRST,N,NNBOPT,M_FLAG,NB_FLAG
*     local variables
      REAL*8 RSM,RS2,A1,A2,A3,DV(3),RIJ2,DR2I,DR3I,DRDV,DRDP
      REAL*8 RI,DP(3),RIJP
      INTEGER I,J
      
!$omp parallel do 
!$omp& private(II,I,RS2,NNB,J,A1,A2,A3,DV,RIJ2,DR2I,DR3I,
!$omp&  DRDV,DRDP,DP,RIJP,RSM,RI)
      DO II = 1, NI
         I = IREG(II+IOFF-1)
         
*     Set neighbor radius limit
 2       RS2 = RS(I)*RS(I)
*     Mass weighted neighbor radius
         IF(M_FLAG.EQ.1) RS2 = RS2/BODYM

*       Start count at 2 and subtract 1 at the end to avoid ILIST(NNB+1).
         NNB = 1
*
*       Initialize scalars for forces & derivatives.
         FREG(1:3,II) = 0.0D0
         FDR(1:3,II) = 0.0

         DO 1 J = IFIRST,N
            IF (J.EQ.I) GO TO 1
            A1 = X(1,J) - X(1,I)
            A2 = X(2,J) - X(2,I)
            A3 = X(3,J) - X(3,I)
*     Predicted coordinates avoids spurious force differences.
            DV(1) = XDOT(1,J) - XDOT(1,I)
            DV(2) = XDOT(2,J) - XDOT(2,I)
            DV(3) = XDOT(3,J) - XDOT(3,I)
*     
            RIJ2 = A1*A1 + A2*A2 + A3*A3
*     Velocity crterion
            DP(1) = A1 + DV(1)*DTR(I)
            DP(2) = A2 + DV(2)*DTR(I)
            DP(3) = A3 + DV(3)*DTR(I)
            RIJP = DP(1)*DP(1) + DP(2)*DP(2) + DP(3)*DP(3)

            DR2I = 1.0/RIJ2
            DR3I = BODY(J)*DR2I*SQRT(DR2I)
            DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
            DRDP = 3.0*DRDV*DR2I
*     
*     IF mass weighted option is switched on, use RIJ2*MASS for comparing
            RSM = MIN(RIJ2,RIJP)
            IF(M_FLAG.EQ.1) RSM = RSM/BODY(J)
*     Use distance and predicted distance criterion for neighbor list
            IF (RSM.GT.RS2) THEN
               FREG(1,II) = FREG(1,II) + A1*DR3I
               FREG(2,II) = FREG(2,II) + A2*DR3I
               FREG(3,II) = FREG(3,II) + A3*DR3I
               FDR(1,II) = FDR(1,II) + (DV(1) - A1*DRDP)*DR3I
               FDR(2,II) = FDR(2,II) + (DV(2) - A2*DRDP)*DR3I
               FDR(3,II) = FDR(3,II) + (DV(3) - A3*DRDP)*DR3I
            ELSE
*     Increase neighbour counter and obtain current irregular force.
               NNB = NNB + 1
               IF(NNB.LE.LMAX-3) LISTGP(NNB,II) = J
            END IF
*     Obtain potential.
            POT(I) = POT(I) + DR3I*RIJ2
 1       CONTINUE
*     Check neighbor list overflow
         IF(NNB.GT.LMAX-3) THEN
*     Dianostics saved in unit 41
            RI = sqrt(x(1,I)**2 + x(2,I)**2 + x(3,I)**2)
            WRITE (41,3) I, NNB, RS(I), RI
 3          FORMAT (' OVERFLOW! I',I8,'  NBNEW',I6,'  RNB[NB]',F10.4,
     &           '  RI[NB]',F10.4)
*     Reduce the neighbor radius
            IF(NB_FLAG.EQ.1) THEN
               RS(I) = FLOAT(NNBOPT)/FLOAT(NNB)**0.3333*RS(I)
            ELSE
               RS(I) = 0.9*RS(I)
            END IF
*     Retry calculation
            GO TO 2
         END IF
         LISTGP(1,II) = NNB - 1
      END DO
!$omp end parallel do

      RETURN
*
      END
