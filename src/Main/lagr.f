      SUBROUTINE LAGR(C)
*
*
*       Lagrangian radii.
*       -----------------
*
      parameter (NLENS=18)
      INCLUDE 'common6.h'
      INCLUDE 'omp_lib.h'
      INTEGER NPARTC(NLENS), NCORE
      REAL*8  R2(NMAX),RSNGL(NMAX),RBIN(NMAX),C(3)
      REAL*8  FLAGR(NLENS),RLAGR(NLENS),AVMASS(NLENS),AVMRC
      REAL*8  RSLAGR(NLENS),RBLAGR(NLENS)
      REAL*8  SIGR2(NLENS),SIGT2(NLENS),SIG2(NLENS)
      REAL*8  SIGR2C,SIGT2C,SIG2C
      REAL*8  VROT(NLENS),VROTC
      REAL*8  VRAVE(NLENS),VTAVE(3,NLENS),VAVE(3,NLENS)
      REAL*8  VTAVEV(NLENS),VAVEV(NLENS)
      REAL*8  VRI(NMAX), VTI(3,NMAX),VTITMP(3),VITMP(3)
      REAL*8  VR_CORE,VT_CORE(3),V_CORE(3),V_COREV,VT_COREV
      INTEGER ISLIST(NMAX),IBLIST(NMAX)
*
*     Lagrangian radii fraction of total mass
      DATA FLAGR/0.001D0,0.003D0,0.005D0,0.01D0,0.03D0,0.05D0,0.1D0,
     &     0.2D0,0.3D0,0.4D0,0.5D0,0.6D0,0.7D0,0.8D0,0.9D0,0.95D0,
     &     0.99D0,1.0D0/
*
*     half mass radius index
      DATA IRLAGRH /11/
      
      IF(KZ(7).EQ.2.OR.KZ(7).EQ.4) THEN
*     Get initial total mass
         ZMASS0 = ZSMASS0 + ZBMASS0
      ELSE
         ZMASS0 =ZMASS
      END IF

*     Particle number counts
      NP = 0
      NSNGL = 0
      NBIN  = 0
      IF (KZ(8).GT.0) THEN
*     Need to exclude massive black hole mass
*     Set square radii of resolved binaries
         IF (KZ(24).EQ.1) THEN
            NPP = 0
            NCMB = -100000
            DO I = 1,IFIRST-1
               IF(NAME(I).NE.NIMBH) THEN
                  NPP = NPP + 1
                  R2(NPP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                 (X(3,I) - C(3))**2
                  JLIST(NPP) = I
               ELSE
                  NCMB = KVEC(I) + N
               END IF
            END DO
*
*     Set square radii of single stalrs
            DO I = IFIRST,N
               IF(NAME(I).NE.NIMBH) THEN
                  NSNGL = NSNGL + 1
                  R2(NPP+NSNGL) = (X(1,I) - C(1))**2 + 
     &                 (X(2,I) - C(2))**2 + (X(3,I) - C(3))**2
                  RSNGL(NSNGL) = R2(NPP+NSNGL)
                  JLIST(NPP+NSNGL) = I
                  ISLIST(NSNGL) = I
               END IF
            END DO
*     Set square radii of c.m.
            DO I = N+1, NTOT
               IF(NAME(I).NE.NCMB) THEN
                  NBIN = NBIN + 1
                  RBIN(NBIN) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                 (X(3,I) - C(3))**2
                  IBLIST(NBIN) = I
               END IF
            END DO
            NP = NPP + NSNGL
         ELSE
*     Set square radii of resolved binaries
!$omp parallel do private(I)
            DO I = 1,IFIRST-1
               R2(I) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &              (X(3,I) - C(3))**2
               JLIST(I) = I
            END DO
!$omp end parallel do

*
*     Set square radii of single stalrs
!$omp parallel do private(I)
            DO I = IFIRST,N
               R2(I) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &              (X(3,I) - C(3))**2
               RSNGL(I-IFIRST+1) = R2(I)
               JLIST(I) = I
               ISLIST(I-IFIRST+1) = I
            END DO
!$omp end parallel do

*     Set square radii of c.m.
!$omp parallel do private(I)
            DO I = N+1, NTOT
               RBIN(I-N) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &              (X(3,I) - C(3))**2
               IBLIST(I-N) = I
            END DO
!$omp end parallel do
            NP = N
            NBIN = NTOT - N
            NSNGL = N - IFIRST + 1
         END IF

*     Sort square distances of all particles with respect to the centre C.
         CALL SORT1(NP,R2,JLIST)
*     Sort square distances of all singles with respect to the centre C.
         CALL SORT1(NSNGL,RSNGL,ISLIST)
*     Sort square distances of c.m. with respect to the centre C.
         CALL SORT1(NBIN,RBIN,IBLIST)
*
      ELSE
*     Only consider singles and resolved K.S.
*     Set square radii
*     exclude massive black hole
         IF (KZ(24).EQ.1) THEN
            DO I = 1,N
               IF(NAME(I).NE.NIMBH) THEN
                  NP = NP + 1
                  R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                 (X(3,I) - C(3))**2
                  JLIST(NP) = I
               END IF
            END DO
         ELSE
!$omp parallel do private(I)
            DO I = 1,N
               R2(I) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &              (X(3,I) - C(3))**2
               JLIST(I) = I
            END DO
!$omp end parallel do
            NP = N
         END IF
*       Sort square distances with respect to the centre C.
         CALL SORT1(NP,R2,JLIST)
      END IF

      IF(KZ(7).EQ.1.OR.KZ(7).EQ.2.OR.KZ(7).EQ.4) THEN
*     Determine half-mass radius by using current total mass.
         ZM = 0.0
         ZMH = 0.5D0*ZMASS
*     If massive black hole exist, exclude it from the rscale determination
         IF (KZ(24).EQ.1) ZMH = ZMH - BIMBH/ZMBAR
         I = 0
 1       I = I + 1
         IM = JLIST(I)
         ZM = ZM + BODY(IM)
         IF (ZM.LT.ZMH) GO TO 1

*     Replace approximate half-mass radius by actual value.
         RSCALE = SQRT(R2(I))
      END IF

*     escape calculation if only rscale is needed
      IF (KZ(7).EQ.1) GO TO 100
*     
*  Determine the Lagrangian radii for specified mass fractions.
*     RLAGR = Lagrangian radius
*     AVMASS = average mass of a spherical shell with radius R2(I)
*     AVMRC = average mass inside of core radius RC
*     NPARTC = particle counter within a shell
*     NCORE = particle counter for the core
*     RC = Core radius (calculated in core.f)
*     VAVE = mass weighted average velocity within a shell
*     VTAVE = mass weighted average tangential velocity within a shell
*     VRAVE = mass weighted average radial velocity within a shell
*     SIG2 = mass weighted velocity dispersion square within a shell
*     SIGR2 = mass weighted radius velocity dispersion square within a shell
*     SIGT2 = mass weighted tangential velocity dispersion square within a shell
*     VROT = mass weighted average rotational velocity projected in x-y plane within a shell
      VR_CORE = 0.0D0
      VT_CORE(1:3) = 0.0D0
      V_CORE(1:3) = 0.0D0
      SIGR2C = 0.0D0
      SIGT2C = 0.0D0
      SIG2C =0.0D0
      VROTC = 0.0D0
      AVMRC = 0.D0
      RC2 = RC*RC
      NCORE = 0
      I = 0
      ZM = 0.D0
*
      DO 15 J = 1,NLENS
         IF(KZ(7).EQ.2.OR.KZ(7).EQ.3.OR.J.EQ.1) THEN
            VROT(J) = 0.0D0
            VRAVE(J) = 0.0D0
            VAVE(1:3,J) = 0.0D0
            VTAVE(1:3,J) = 0.0D0
            AVMASS(J) = 0.0D0
            NPARTC(J) = 0
         ELSE
            VROT(J) = VROT(J-1)
            VRAVE(J) = VRAVE(J-1)
            VAVE(1:3,J) = VAVE(1:3,J-1)
            VTAVE(1:3,J) = VTAVE(1:3,J-1)
            AVMASS(J) = AVMASS(J-1)
            NPARTC(J) = NPARTC(J-1)
         END IF
*
 20      I = I + 1
*     If reach last particle finish the loop
         IF(I.GT.NP) THEN
            IF(NPARTC(J).NE.0) THEN
               II = I - 1
 31            IM = JLIST(II)
               IF(BODY(IM).EQ.0.0D0) THEN
                  II = II - 1
                  GO TO 31
               ELSE
                  RLAGR(J) = SQRT(R2(II))
                  I = I - 1
                  GO TO 15
               END IF
            ELSE
               RLAGR(J) = RLAGR(J-1)
               VTAVEV(J) = VTAVEV(J-1)
               VAVEV(J) = VAVEV(J-1)
               I = I - 1
               GO TO 15
            END IF
         END IF
*     Sorted list index
         IM = JLIST(I)
*     escape ghost particle, but add it into particle count
         IF(BODY(IM).EQ.0.0D0) THEN
            NPARTC(J) = NPARTC(J) + 1
            GO TO 20
         END IF
*     Cumulated mass
         ZM = ZM + BODY(IM)
*     Cumulated velocity
         DO K = 1,3
            VAVE(K,J) = VAVE(K,J) + BODY(IM)*XDOT(K,IM)
         END DO
*     Radial velocity can be calculated by velocity .dot. radial direction vector
         VRI(I) = 0.D0
         DO K = 1,3
            VRI(I) = VRI(I) + XDOT(K,IM)*(X(K,IM)-C(K))/DSQRT(R2(I))
         END DO
*     Cumulate radial velocity
         VRAVE(J) = VRAVE(J) + BODY(IM)*VRI(I)
*     Radial velocity square
         VR2I = VRI(I)*VRI(I)
*     Tangential velocity
         DO K = 1,3
            VTI(K,I) = XDOT(K,IM) - VRI(I)*(X(K,IM)-C(K))/DSQRT(R2(I))
*     Cumulate tangential velocity
            VTAVE(K,J) = VTAVE(K,J) + BODY(IM)*VTI(K,I)
         END DO
*     X-Y plane projected position square
         RR12 = X(1,IM)**2 + X(2,IM)**2
*     X-Y plane radial velocity value * sqrt(RR12) 
         XR12 = 0.D0
         DO K = 1,2
            XR12 = XR12 + XDOT(K,IM)*(X(K,IM)-C(K))
         END DO
*     Rotational velocity
         VROT1 = XDOT(1,IM) - XR12/RR12*X(1,IM)
         VROT2 = XDOT(2,IM) - XR12/RR12*X(2,IM)
*     Rotational direction sign
         XSIGN = VROT1*X(2,IM)/DSQRT(RR12) - VROT2*X(1,IM)/DSQRT(RR12)
         VROTM = DSQRT(VROT1**2+VROT2**2)
*     Cumulate rotational velocity
         IF(XSIGN.GT.0.D0) THEN
            VROT(J) = VROT(J) + BODY(IM)*VROTM
         ELSE
            VROT(J) = VROT(J) - BODY(IM)*VROTM
         END IF
*     Cumulate mass and particle number
         AVMASS(J) = AVMASS(J) + BODY(IM)
         NPARTC(J) = NPARTC(J) + 1
*
*     Independent determination of mass in core radius.
         IF (R2(I).LT.RC2) THEN
            VR_CORE = VR_CORE + BODY(IM)*VRI(I)
            DO K = 1,3
               V_CORE(K) = V_CORE(K) + BODY(IM)*XDOT(K,IM)
               VT_CORE(K) = VT_CORE(K) + BODY(IM)*VTI(K,I)
            END DO
            IF(XSIGN.GT.0.D0) THEN
               VROTC = VROTC + BODY(IM)*VROTM
            ELSE
               VROTC = VROTC - BODY(IM)*VROTM
            END IF
            AVMRC = AVMRC + BODY(IM)
            NCORE = NCORE + 1
         END IF
*
*     Check whether mass within Langrangian radius is complete.
         IF (I.LT.NP.AND.ZM.LT.FLAGR(J)*ZMASS0) GO TO 20

*     Get average within a shell
         RLAGR(J) = SQRT(R2(I))
 15   CONTINUE

*     Get Half mass radius
      IF(KZ(7).NE.1.AND.KZ(7).NE.2.AND.KZ(7).NE.4) 
     &     RSCALE = RLAGR(IRLAGRH)

*     Get average within core radius
      VR_CORE = VR_CORE/AVMRC
      V_COREV = 0.0D0
      VT_COREV = 0.0D0
      DO K = 1,3
         VT_CORE(K) = VT_CORE(K)/AVMRC
         VT_COREV = VT_COREV + VT_CORE(K)*VT_CORE(K)
         V_CORE(K) = V_CORE(K)/AVMRC
         V_COREV = V_COREV + V_CORE(K)*V_CORE(K)
      END DO
      VROTC = VROTC/AVMRC

*     Get velocity dispersion
      NSTART = 1
      NSHELL = 0
      DO 16 J = 1,NLENS
         NSTART = NSHELL + 1
         IF(KZ(7).EQ.2.OR.KZ(7).EQ.3.OR.J.EQ.1) THEN
            SIGR2(J) = 0.0D0
            SIGT2(J) = 0.0D0
            SIG2(J) = 0.0D0
            NSHELL = NSHELL + NPARTC(J)
         ELSE
            SIGR2(J) = SIGR2(J-1)
            SIGT2(J) = SIGT2(J-1)
            SIG2(J) = SIG2(J-1)
            NSHELL = NPARTC(J)
         END IF

         DO IK = NSTART, NSHELL
            IM = JLIST(IK)
*     Radial direction
            VR2I = VRI(IK) - VRAVE(J)
            VR2I = VR2I*VR2I
*     Tangential direction and original direction
            VTI2 = 0.0D0
            VI2 = 0.0D0
            DO K = 1,3
               VTITMP(K) = VTI(K,IK) - VTAVE(K,J)
               VTI2 = VTI2 + VTITMP(K)*VTITMP(K)
               VITMP(K) = XDOT(K,IM) - VAVE(K,J)
               VI2 = VI2 + VITMP(K)*VITMP(K)
            END DO
*     Cumulate radial velocity dispersion
            SIGR2(J) = SIGR2(J) + BODY(IM)*VR2I
*     Cumulate tangential velocity dispersion
            SIGT2(J) = SIGT2(J) + BODY(IM)*VTI2/2.D0
*     Cumulate velocity dispersion
            SIG2(J) = SIG2(J) + BODY(IM)*VI2/3.D0
*     Only for core radius
            IF (R2(IK).LT.RC2) THEN
               VR2I = VRI(IK) - VR_CORE
               VR2I = VR2I*VR2I
               SIGR2C = SIGR2C + BODY(IM)*VR2I

               VTI2 = 0.0D0
               VI2 = 0.0D0
               DO K = 1,3
                  VTITMP(K) = VTI(K,IK) - VT_CORE(K)
                  VTI2 = VTI2 + VTITMP(K)*VTITMP(K)
                  VITMP(K) = XDOT(K,IM) - V_CORE(K)
                  VI2 = VI2 + VITMP(K)*VITMP(K)
               END DO
               SIGT2C = SIGT2C + BODY(IM)*VTI2/2.D0
               SIG2C = SIG2C + BODY(IM)*VI2/3.D0
            END IF
         END DO
 16   CONTINUE

*     Average velocity dispersion for core region
      SIGR2C = SIGR2C/AVMRC
      SIGT2C = SIGT2C/AVMRC
      SIG2C = SIG2C/AVMRC
      AVMRC = AVMRC/NCORE
*
*     Final average
      DO J = 1, NLENS
         IF(NPARTC(J).GT.0) THEN
            VRAVE(J) = VRAVE(J)/AVMASS(J)
            VTAVEV(J) = 0.0D0
            VAVEV(J) = 0.0D0
            DO K = 1,3
               VTAVE(K,J) = VTAVE(K,J)/AVMASS(J)
               VTAVEV(J) = VTAVEV(J) + VTAVE(K,J)*VTAVE(K,J)
               VAVE(K,J) = VAVE(K,J)/AVMASS(J)
               VAVEV(J) = VAVEV(J) + VAVE(K,J)*VAVE(K,J)
            END DO
            VTAVEV(J) = SQRT(VTAVEV(J))
            VAVEV(J) = SQRT(VAVEV(J))
            VROT(J) = VROT(J)/AVMASS(J)

            SIGR2(J) = SIGR2(J)/AVMASS(J)
            SIGT2(J) = SIGT2(J)/AVMASS(J)
            SIG2(J) = SIG2(J)/AVMASS(J)
            AVMASS(J) = AVMASS(J)/NPARTC(J)
         END IF
      END DO

      IF (KZ(8).GT.0) THEN
*  Determine the Lagrangian radii for singles and binaries separately.
         ZMS = 0.0D0
         I = 0
*
         DO 17 J = 1, NLENS
 21         I = I + 1
            IF(I.GT.NSNGL) THEN
               RSLAGR(J) = RSLAGR(J-1)
               GO TO 17
            END IF
            IM = ISLIST(I)
            IF(BODY(IM).EQ.0.0D0) GO TO 21
            ZMS = ZMS + BODY(IM)
            IF (I.LT.NSNGL.AND.ZMS.LT.FLAGR(J)*ZSMASS0) GO TO 21
            RSLAGR(J) = SQRT(RSNGL(I))
 17      CONTINUE
*
         ZMB = 0.0D0
         I = 0
*     
         DO 18 J = 1, NLENS
 22         I = I + 1
            IF(I.GT.NBIN) THEN
               RBLAGR(J) = RBLAGR(J-1)
               GO TO 18
            END IF
            IM = IBLIST(I)
            IF(BODY(IM).EQ.0.0D0) GO TO 22
            ZMB = ZMB + BODY(IM)
            IF (I.LT.NBIN.AND.ZMB.LT.FLAGR(J)*ZBMASS0) GO TO 22
            RBLAGR(J) = SQRT(RBIN(I))
 18      CONTINUE
*
      END IF
*
*     Write on diagnostics of RLAGR, AVMASS, NPARTC.
      if(rank.eq.0)then
         IF (KZ(7).GE.2) THEN
            WRITE (6,40) (FLAGR(K),K=1,NLENS)
 40         FORMAT (/,11X,'TIME   M/MT:',1P,18(1X,D9.2),7X,'<RC')
            WRITE (6,41) TTOT, (RLAGR(K),K=1,NLENS),RC
 41         FORMAT (3X,D12.4,' RLAGR: ',1P,19(1X,D9.2))
*
            IF (KZ(8).GT.0 .OR. NBIN0.GT.0) THEN
               WRITE (6,401) TTOT, (RSLAGR(K),K=1,NLENS)
 401           FORMAT (3X,D12.4,' RSLAGR: ',1P,18(1X,D9.2))
               WRITE (6,402) TTOT, (RBLAGR(K),K=1,NLENS)
 402           FORMAT (3X,D12.4,' RBLAGR: ',1P,18(1X,D9.2))
            END IF
*
            WRITE (6,42) TTOT, (AVMASS(K),K=1,NLENS),AVMRC
 42         FORMAT (3X,D12.4,' AVMASS:',1P,19(1X,D9.2))
            WRITE (6,43) TTOT, (NPARTC(K),K=1,NLENS),NCORE
 43         FORMAT (3X,D12.4,' NPARTC:',19I10)
         END IF
      end if
*
*     Write all data on unit 7.
      if(rank.eq.0)then
         IF (KZ(7).GE.2) THEN
            IF(kstart.eq.1.and.time.eq.0.0) then
               write (7,50) 2,3+NLENS,3+2*NLENS,3+3*NLENS,4+4*NLENS,
     &              5+5*NLENS,6+6*NLENS,7+7*NLENS,8+8*NLENS,9+9*NLENS,
     &              10+10*NLENS,11+11*NLENS,12+12*NLENS,13+13*NLENS,
     &              14+14*NLENS
 50            format('## Time, [Column number | Label]: ',
     &              I4,' R_{lagr} ',I4,' R_{lagr,s} ',
     &              I4,' R_{lagr,b} ',I4,' <M> ',I4,' N_{shell} ',
     &              I4,' <V_x> ',I4,' <V_y> ',I4,' <V_z> ',I4,'<V> ',
     &              I4,' <V_r> ',I4,' <V_t> ',I4,' \Sigma^2 ',
     &              I4,' \Sigma_r^2 ',I4,' \Sigma_t^2 ',
     &              I4,' <V_{rot.}> ',
     &              '(All units are NB)')
               write (7,51) 'TIME',
*     1-5               
     &              FLAGR(1:NLENS),'<RC',
     &              FLAGR(1:NLENS),
     &              FLAGR(1:NLENS),
     &              FLAGR(1:NLENS),'<RC',
     &              FLAGR(1:NLENS),'<RC',
*     6-10
     &              FLAGR(1:NLENS),'<RC',
     &              FLAGR(1:NLENS),'<RC',
     &              FLAGR(1:NLENS),'<RC',
     &              FLAGR(1:NLENS),'<RC',
     &              FLAGR(1:NLENS),'<RC',
*     11-15
     &              FLAGR(1:NLENS),'<RC',
     &              FLAGR(1:NLENS),'<RC',
     &              FLAGR(1:NLENS),'<RC',
     &              FLAGR(1:NLENS),'<RC',
     &              FLAGR(1:NLENS),'<RC'
 51            format(A4,1P,1X,18E10.2,1X,A4,2(1X,18E10.2),
     &              12(1X,18E10.2,1X,A4))
            end if

            write (7,60) TTOT,
     &           RLAGR(1:NLENS),RC,
     &           RSLAGR(1:NLENS),RBLAGR(1:NLENS),
     &           AVMASS(1:NLENS),AVMRC,
     &           NPARTC(1:NLENS),NCORE,
     &           VAVE(1,1:NLENS),V_CORE(1),
     &           VAVE(2,1:NLENS),V_CORE(2),
     &           VAVE(3,1:NLENS),V_CORE(3),
     &           VAVEV(1:NLENS),V_COREV,
     &           VRAVE(1:NLENS),VR_CORE,
     &           VTAVEV(1:NLENS),VT_COREV,
     &           SIG2(1:NLENS),SIG2C,
     &           SIGR2(1:NLENS),SIGR2C,
     &           SIGT2(1:NLENS),SIGT2C,
     &           VROT(1:NLENS),VROTC
 60         format(75E26.17,19I12,190E26.17)
            CALL FLUSH(7)
         END IF
      end if
*
 100  RETURN
*
      END


