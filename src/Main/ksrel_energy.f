      subroutine ksrel_energy(IPAIR,M1,M2,XI,VI,DT,FRELP,FRELD,KZ45,
     &     TIME)
*     
*     Relativistic energy correction
*     ------------------------------------------------------------
*     Caluculations of the relativistic energy (used in adjust.F).
*     ------------------------------------------------------------
*     
*     This calculation approximates the relativistic energy associated
*     with the binary using a semi-classical radiative approach:
*     E = int(F dot v)dt ~ F dot v *timestep
*     + (dF/dt dot v + F dot dv/dt)*timestep^2
*     and should provide a self-consistent relativistic energy for
*     high values of c
      include 'params.h'
      INCLUDE 'postnewton.h'
      REAL*8 XI(6),VI(6),M1,M2,FRELP(3),FRELD(3),DT,TIME
      INTEGER IPAIR,KZ45
      REAL*8 KSV(3),KSX(3),MASAS,MOR,MU
      REAL*8 NX,NY,NZ,NDV
      REAL*8 RP,CL,CL2,CL4,Cl5,CL6,RIJ2,RIJ,PI2
      REAL*8 A25POWER(3),A25DPOWER(3),AT(3),AN(3)
C      REAL*8 A2_5,B2_5,A2_5D,B2_5D
      REAL*8 DE25N,DEALL,DDEALL
C      REAL*8 AS(3),AP(3)
C      REAL*8 POWER1S,POWER1P,POWER2S,POWER2P
      REAL*8 V12,V22,N12V1,N12V2,V1V2,V1_V22
C      REAL*8 EKIN,ENEWT

      PI2 = 9.869604401089359
*     Set speed of light CL and its powers.
      CL = CLIGHTPN
      CL2 = CL*CL
      CL4 = CL2*CL2
      CL5 = CL4*CL
      CL6 = CL2*CL4
      
      DO K = 1,3
         KSX(K) = XI(K)-XI(K+3)       
         KSV(K) = VI(K)-VI(K+3)
      END DO

      MASAS = M1 + M2
      MU = M1*M2/MASAS
C      ETA = MU/MASAS
      
      RIJ2 = KSX(1)**2+KSX(2)**2+KSX(3)**2
      RIJ = SQRT(RIJ2)

      NX = KSX(1)/RIJ
      NY = KSX(2)/RIJ
      NZ = KSX(3)/RIJ

      MOR = MASAS/RIJ

      AN(1) = - MOR*NX/RIJ
      AN(2) = - MOR*NY/RIJ
      AN(3) = - MOR*NZ/RIJ
      
      AT(1) = FRELP(1) + AN(1)
      AT(2) = FRELP(2) + AN(2)
      AT(3) = FRELP(3) + AN(3)
      V1_V22 = KSV(1)*KSV(1)+KSV(2)*KSV(2)+KSV(3)*KSV(3)

      RP = (KSX(1)*KSV(1)+KSX(2)*KSV(2)+KSX(3)*KSV(3))/RIJ
      NDV = NX*KSV(1) + NY*KSV(2) + NZ*KSV(3)

      V12 = VI(1)*VI(1)+VI(2)*VI(2)+VI(3)*VI(3)
      V22 = VI(4)*VI(4)+VI(5)*VI(5)+VI(6)*VI(6)
      N12V1 = NX*VI(1)+NY*VI(2)+NZ*VI(3)
      N12V2 = NX*VI(4)+NY*VI(5)+NZ*VI(6)
      V1V2 = VI(1)*VI(4)+VI(2)*VI(5)+VI(3)*VI(6)
*     
**    Suppressed 
C      A2_5  = PN2_5(1)
C      A2_5D = PN2_5(2)
C      B2_5  = PN2_5(3)
C      B2_5D = PN2_5(4)

*     PN1,2,2.5,3 energy correction
      DIFFS(4) = (M1*M1*M2/(2.*RIJ2)+3.*M1*V12*V12/8.+M1*M2/RIJ
     &     *(-0.25*(N12V1)
     &     *(N12V2)+1.5*V12
     &     -1.75*(V1V2)))/CL2
     &     +(M2*M2*M1/(2.*RIJ2)+3.*M2*V22*V22/8.+M1*M2/RIJ
     &     *(-0.25*(N12V1)
     &     *(N12V2)+1.5*V22
     &     -1.75*(V1V2)))/CL2
     &     +(-M1*M1*M1*M2/(2./RIJ2*RIJ)-19.*M1*M1*M2*M2
     &     /(8.*RIJ2*RIJ)+5.*M1*V12**3/16.+M1*M2/RIJ*(3./8.*
     &     N12V1**3*N12V2+3./16.*N12V1**2*N12V2**2-9./8.*N12V1
     &     *N12V2*V12-13./8.*N12V2**2*V12+21./8.*V12**2+13./8.*
     &     N12V1**2*V1V2+3./4.*N12V1*N12V2*V1V2-55./8.*V12*V1V2
     &     +17./8.*V1V2**2+31./16.*V12*V22)+M1*M1*M2/RIJ2*
     &     (29./4.*N12V1**2-13./4.*N12V1*N12V2+0.5*N12V2**2-
     &     1.5*V12+7./4.*V22))/CL4
     &     +(-M2*M2*M2*M1/(2./RIJ2*RIJ)-19.*M2*M2*M1*M1
     &     /(8.*RIJ2*RIJ)+5.*M2*V22**3/16.+M1*M2/RIJ*(3./8.*
     &     N12V2**3*N12V1+3./16.*N12V2**2*N12V1**2-9./8.*N12V2
     &     *N12V1*V22-13./8.*N12V1**2*V22+21./8.*V22**2+13./8.*
     &     N12V2**2*V1V2+3./4.*N12V2*N12V1*V1V2-55./8.*V22*V1V2
     &     +17./8.*V1V2**2+31./16.*V22*V12)+M1*M2*M2/RIJ2*
     &     (29./4.*N12V2**2-13./4.*N12V2*N12V1+0.5*N12V1**2-
     &     1.5*V22+7./4.*V12))/CL4
     &     +(35.*M1*V12**4/128.
     &     +M1*M2/RIJ*(-5./16.*N12V1**5*N12V2-5./16.*N12V1**4
     &     *N12V2**2-5./32.*N12V1**3*N12V2**3
     &     +19./16.*N12V1**3*N12V2*V12+15./16.*N12V1**2
     &     *N12V2**2*V12+3./4.*N12V1*N12V2**3*V12
     &     +19./16.*N12V2**4*V12-21./16.*N12V1*N12V2*V12**2
     &     -2.*N12V2**2*V12**2
     &     +55./16.*V12**3-19./16.*N12V1**4*V1V2-N12V1**3
     &     *N12V2*V1V2
     &     -15./32.*N12V1**2*N12V2**2*V1V2+45./16.*N12V1**2
     &     *V12*V1V2
     &     +1.25*N12V1*N12V2*V12*V1V2+11./4.*N12V2**2*V12
     &     *V1V2-139./16.*V12**2*V1V2
     &     -0.75*N12V1**2*V1V2**2+5./16.*N12V1*N12V2*V1V2**2
     &     +41./8.*V12*V1V2**2+1./16.*V1V2**3
     &     -45./16.*N12V1**2*V12*V22-23./32.*N12V1*N12V2
     &     *V12*V22+79./16.*V12**2*V22-161./32.*V12*V1V2*V22)
     &     +M1*M1*M2/RIJ2*(-49./8.*N12V1**4+75./8.*N12V1**3
     &     *N12V2-187./8.*N12V1**2*N12V2**2
     &     +247./24.*N12V1*N12V2**3+49./8.*N12V1**2*V12
     &     +81./8.*N12V1*N12V2*V12
     &     -21./4.*N12V2**2*V12+11./2.*V12**2-7.5*N12V1**2
     &     *V1V2-1.5*N12V1*N12V2*V1V2
     &     +21./4.*N12V2**2*V1V2-27.*V12*V1V2+55./2.*V1V2**2
     &     +49./4.*N12V1**2*V22
     &     -27./2.*N12V1*N12V2*V22+3./4.*N12V2**2*V22+55./4.*
     &     V12*V22-28.*V1V2*V22+135./16.*V22**2)
     &     +3.*M1**4*M2/(8.*RIJ2*RIJ2)+M1**3*M2**2/RIJ2/RIJ2*
     &     9707./420. +M1**2*M2**2/RIJ2/RIJ*(547./12.*N12V1**2
     &     -3115./48.*N12V1*N12V2-123./64.*N12V1*NDV*PI2-575./18.
     &     *V12+41./64.*PI2*(VI(1)*KSV(1)+VI(2)*KSV(2)+VI(3)*KSV(3))
     &     +4429./144.*V1V2)
     &     +M1**3*M2/RIJ2/RIJ*(-44627./840.*N12V1**2+32027.
     &     /840.*N12V1*N12V2+1.5*N12V2**2+24187./2520.*V12
     &     -27967./2520.*V1V2+1.25*V22))/CL6
     &     +(35.*M2*V22**4/128.
     &     +M1*M2/RIJ*(-5./16.*N12V2**5*N12V1-5./16.*N12V2**4
     &     *N12V1**2-5./32.*N12V2**3*N12V1**3
     &     +19./16.*N12V2**3*N12V1*V22+15./16.*N12V2**2
     &     *N12V1**2*V22+3./4.*N12V2*N12V1**3*V22
     &     +19./16.*N12V1**4*V22-21./16.*N12V2*N12V1*V22**2
     &     -2.*N12V1**2*V22**2
     &     +55./16.*V22**3-19./16.*N12V2**4*V1V2-N12V2**3
     &     *N12V1*V1V2
     &     -15./32.*N12V2**2*N12V1**2*V1V2+45./16.*N12V2**2
     &     *V22*V1V2
     &     +1.25*N12V2*N12V1*V22*V1V2+11./4.*N12V1**2*V22
     &     *V1V2-139./16.*V22**2*V1V2
     &     -0.75*N12V2**2*V1V2**2+5./16.*N12V2*N12V1*V1V2**2
     &     +41./8.*V22*V1V2**2+1./16.*V1V2**3
     &     -45./16.*N12V2**2*V22*V12-23./32.*N12V2*N12V1
     &     *V22*V12+79./16.*V22**2*V12-161./32.*V22*V1V2*V12)
     &     +M1*M2*M2/RIJ2*(-49./8.*N12V2**4+75./8.*N12V2**3
     &     *N12V1-187./8.*N12V2**2*N12V1**2
     &     +247./24.*N12V2*N12V1**3+49./8.*N12V2**2*V22
     &     +81./8.*N12V2*N12V1*V22
     &     -21./4.*N12V1**2*V22+11./2.*V22**2-7.5*N12V2**2
     &     *V1V2-1.5*N12V2*N12V1*V1V2
     &     +21./4.*N12V1**2*V1V2-27.*V22*V1V2+55./2.*V1V2**2
     &     +49./4.*N12V2**2*V12
     &     -27./2.*N12V2*N12V1*V12+3./4.*N12V1**2*V12+55./4.*
     &     V22*V12-28.*V1V2*V12+135./16.*V12**2)
     &     +3.*M2**4*M1/(8.*RIJ2*RIJ2)+M2**3*M1**2/RIJ2/RIJ2*
     &     9707./420. +M2**2*M1**2/RIJ2/RIJ*(547./12.*N12V2**2
     &     -3115./48.*N12V2*N12V1-123./64.*N12V2*NDV*PI2-575./18.
     &     *V22-41./64.*PI2*(VI(4)*KSV(1)+VI(5)*KSV(2)+VI(6)*KSV(3)) !-41/64 sign changed by asimetry
     &     +4429./144.*V1V2)
     &     +M2**3*M1/RIJ2/RIJ*(-44627./840.*N12V2**2+32027.
     &     /840.*N12V2*N12V1+1.5*N12V1**2+24187./2520.*V22
     &     -27967./2520.*V1V2+1.25*V12))/CL6
*******
     &     +4./5.*M1*M1*M2/(CL5*RIJ2)*N12V1*(V1_V22
     &     -2.*(M1-M2)/RIJ)
     &     -4./5.*M2*M2*M1/(CL5*RIJ2)*N12V2*(V1_V22
     &     -2.*(M2-M1)/RIJ)


*     Use PN energy formula correction
*     PN radiation energy
      DE25N = 8./15./CL5*(M1*M2)**2/(RIJ2*RIJ2)*
     &     (12*V1_V22-11*(NX*KSV(1)+NY*KSV(2)+NZ*KSV(3))**2)
      
      IF(RFLAG(IPAIR).EQ.0) THEN
         DE25(1,IPAIR) = DE25N
      END IF

      DIFFS(5) = DIFFS(5) + 0.5*(DE25N+DE25(1,IPAIR))*DT

      DE25(1,IPAIR) = DE25N


      IF(KZ45.EQ.2) THEN
*     Correction by integrating F dot v dt with 4th order accuracy
*     ----------------------
         DO 173 MM = 1,3
c$$$  AP(MM) = M2*FRELP(MM)/MASAS
c$$$  AS(MM) = -M1*FRELP(MM)/MASAS
C     AS(MM) = -M1/M2*AP(MM)
            A25POWER(MM) = FRELP(MM)
            A25DPOWER(MM) = FRELD(MM)
 173     CONTINUE
*     
*     Calculte the power radiated at each order. (suppressed)
c$$$      POWER1P = (AP(1)*VI(1)+AP(2)*VI(2)+
c$$$     &     AP(3)*VI(3))*M1
c$$$      POWER1S = (AS(1)*VI(4)+AS(2)*VI(5)+
c$$$     &     AS(3)*VI(6))*M2
c$$$      DIFFS(1) = POWER1P
c$$$      DIFFS(2) = POWER1S
c$$$      DIFFS(3) = (POWER1P+POWER1S)

*     DE and DEDOT
         DEALL = MU*(A25POWER(1)*KSV(1)
     &        + A25POWER(2)*KSV(2)
     &        + A25POWER(3)*KSV(3))
         DDEALL = MU*((A25DPOWER(1)*KSV(1)
     &        + A25DPOWER(2)*KSV(2)
     &        + A25DPOWER(3)*KSV(3))
     &        + A25POWER(1)*AT(1)
     &        + A25POWER(2)*AT(2)
     &        + A25POWER(3)*AT(3))
         
*     Check whether there is break of PN treatment
*     Notice here DE25(2,IPAIR) is used to save the previous time, not energy!
         IF(RFLAG(IPAIR).EQ.1.AND.
     &        (TIME-DE25(2,IPAIR))/DT.GT.1.1) THEN
C            print*,'TIME,T0',TIME,DE25(2,IPAIR),'DT',DT
            RFLAG(IPAIR) = 0
         END IF

         DE25(2,IPAIR) = TIME

         IF(RFLAG(IPAIR).EQ.0) THEN
            PN1ENERGY(1,IPAIR) = DEALL
            PN1ENERGY(2,IPAIR) = DDEALL
            RFLAG(IPAIR) = 1
         END IF
C         ELSE
         PREVRELENERGY = PREVRELENERGY 
     &        - 0.5*(DEALL+PN1ENERGY(1,IPAIR))*DT
     &        - 1.0/12.0*(PN1ENERGY(2,IPAIR)-DDEALL)*DT*DT

         PN1ENERGY(1,IPAIR) = DEALL
         PN1ENERGY(2,IPAIR) = DDEALL

C         END IF
*     
      ELSE
         IF (RFLAG(IPAIR).EQ.0) THEN
            RFLAG(IPAIR) = 1
            PN1ENERGY(2,IPAIR) = DIFFS(4)
         END IF
         PN1ENERGY(1,IPAIR) = PN1ENERGY(2,IPAIR)
         PN1ENERGY(2,IPAIR) = DIFFS(4)
         PREVRELENERGY = PREVRELENERGY + PN1ENERGY(2,IPAIR)
     &        -PN1ENERGY(1,IPAIR) + DIFFS(5)
      END IF


C     E_CORRECT = DIFFS(3)*DT
C     IF (E_CORRECT.GT.0.) THEN
C     PREVRELENERGY = PREVRELENERGY + E_CORRECT
C     RELENERGY(IPAIR) = RELENERGY(IPAIR) + DIFFS(5) - DIFFS(4)
C     ELSE
C     PREVRELENERGY = PREVRELENERGY - E_CORRECT
C     RELENERGY(IPAIR) = RELENERGY(IPAIR) - E_CORRECT
C     END IF
*     
*     RELTSTEP is defined in ksint.F
C     RELENERGY(IPAIR) = RELENERGY(IPAIR)
C     &     + (POWER1P + POWER1S)*DT
C     PREVRELPOTENERGY = PREVRELPOTENERGY + 
C     &   (POWER1P + POWER1S)*DT
C     write (*,*) "RELTSTEP AND STEP",RELTSTEP(IPAIR),STEP(IPAIR)
*     
*     Second order correction
*     -----------------------


C     RELENERGY(IPAIR) = RELENERGY(IPAIR)+
C     &     (POWER2P+POWER2S)*DT*DT

*     =================================================================
*     Calculate additional influence of perturbers on energy correction
*     =================================================================
*     
*     Sum over perturberlist
*     Acc of pair particles given as AP(3) and AS(3) above
C     NNB2 = LIST(1,I1) + 1
C     DO 101 J = 2, NNB2
C     K = LIST(J,I1)
C     A1 = X(1,K) - XI(1)
C     A2 = X(2,K) - XI(2)
C     A3 = X(3,K) - XI(3)
C     RIJ2 = A1*A1 + A2*A2 + A3*A3
C     RIJ = SQRT(RIJ2)
C     RELPOTENERGY(IPAIR) = RELPOTENERGY(IPAIR)+
C     &    M1*BODY(K)*(AP(1)*A1+AP(2)*A2+AP(3)*A3)*
C     &    DT*DT/(2.*RIJ2*RIJ)
C     *     And for 2nd body in pair
C     A1 = X(1,K) - XI(4)
C     A2 = X(2,K) - XI(5)
C     A3 = X(3,K) - XI(6)
C     RIJ2 = A1*A1 + A2*A2 + A3*A3
C     RIJ = SQRT(RIJ2)
C     RELPOTENERGY(IPAIR) = RELPOTENERGY(IPAIR)+
C     &    M2*BODY(K)*(AS(1)*A1+AS(2)*A2+AS(3)*A3)*
C     &    DT*DT/(2.*RIJ2*RIJ)
C     101  CONTINUE
*     ===============================================================
*     
*     
*     RELTSTEP is defined in ksint.F
C     WRITE (*,*) 'debugg', (FRELP(1)*KSV(1)+FRELP(2)*KSV(2)+
C     &     FRELP(3)*KSV(3))*RELTSTEP(IPAIR)*MU
C     RELENERGY(IPAIR) = RELENERGY(IPAIR)
C     &    + (FRELP(1)*KSV(1)+FRELP(2)*KSV(2)+
C     &     FRELP(3)*KSV(3))*DT*MU

C     &     + (POWER2P + POWER2S)*RELTSTEP(IPAIR)*RELTSTEP(IPAIR)

C     RELSUM = RELSUM + (POWER1P+POWER1S)*RELTSTEP(IPAIR)
C     RELSUM = RELSUM + (A25POWER(1)*KSV(1)+A25POWER(2)*KSV(2)+
C     &     A25POWER(3)*KSV(3))*RELTSTEP(IPAIR)*MASAS
C     RELSUM = RELSUM + (FRELP(1)*KSV(1)+FRELP(2)*KSV(2)+
C     &     FRELP(3)*KSV(3))*RELTDT*MU
C     RELENERGY(IPAIR) = RELENERGY(IPAIR) + 
C     & 0.5*(FRELD(1)*KSV(1)+FRELD(2)*KSV(2)+
C     & FRELD(3)*KSV(3)+
C     & FRELP(1)*AT(1)+FRELP(2)*AT(2)+FRELP(3)*AT(3))*
C     & DT*DT*MU
*     WRITE (*,*) 'debug ', POWER1P, POWER1S, POWER2P, POWER2S
*     write (*,*) "2nd order ",RELENERGY2_5(IPAIR)
*     
*     CALL RENG (IPAIR,I1,KSX,KSV)
*     
*     -------------------------------------------
*     End calculation of the relativistic energy.
*     -------------------------------------------
*     


C     WRITE (*,*) 'Relativistic'
*     checking derivative and other tests***
*     DIFFS(1) = (C2_5(1)-DIFFS(2))/STEP(1)
*     DIFFS(2) = C2_5(1)
*     * CHECK DERIVATIVES **
*     DIFFS(1) = (DIFFS(1)-C2D(3))/C2D(3)
C     DIFFS(2) = SQRT(L(1)*L(1)+L(2)*L(2)+L(3)*L(3))*CL/(BODY(1)+
C     &   BODY(2))**2
C     DIFFS(3) = SQRT(SPN(1,1)**2+SPN(2,1)**2+SPN(3,1)**2)
C     DIFFS(4) = SQRT(SPN(1,2)**2+SPN(2,2)**2+SPN(3,2)**2)
C     DIFFS(1) = SQRT((L(1)+SPN(1,1)*BODY(1)*BODY(1)/CL+SPN(1,2)*
C     &   BODY(2)*BODY(2)/CL)**2+(L(2)+SPN(2,1)*BODY(1)*BODY(1)/CL+
C     &   SPN(2,2)*BODY(2)*BODY(2)/CL)**2+(L(3)+SPN(3,1)*BODY(1)*
C     &   BODY(1)/CL+SPN(3,2)*BODY(2)*BODY(2)/CL)**2)*CL/(BODY(1)+
C     &   BODY(2))**2
C     IF ((A3D.GT.1.0D30).OR.(A3D.LT.-1.0D30)) THEN
*     DIFFS(1) = (DIFFS(1)-C2_5D(1))/C2_5D(1)
C     ELSE
C     DIFFS(1) = 0.
C     END IF
*     /checking derivative**
C     DIFFS(1) = SQRT((SU3(1)*SU3(1)+SU3(2)*SU3(2)+SU3(3)*SU3(3))/
C     &  (KSS(1)*KSS(1)+KSS(2)*KSS(2)+KSS(3)*KSS(3)))*STEP(1)

*     * PN strength survey
C     ANEWTON = MOR/RIJ
C     DIFFS(1) = SQRT((A1*NX+B1*KSV(1))**2+(A1*NY+B1*KSV(2))**2+
C     & (A1*NZ+B1*KSV(3))**2)/CL2
C     DIFFS(2) = SQRT((A2*NX+B2*KSV(1))**2+(A2*NY+B2*KSV(2))**2+
C     & (A2*NZ+B2*KSV(3))**2)/CL4
C     DIFFS(3) = SQRT((A2_5*NX+B2_5*KSV(1))**2+(A2_5*NX+B2_5*KSV(1))**2+
C     & (A2_5*NZ+B2_5*KSV(3))**2)/CL5
C     DIFFS(4) = SQRT((A3*NX+B3*KSV(1))**2+(A3*NY+B3*KSV(2))**2+
C     & (A3*NZ+B3*KSV(3))**2)/CL6
C     DIFFS(5) = SQRT((A3_5*NX+B3_5*KSV(1))**2+(A3_5*NY+B3_5*KSV(2))**2+
C     & (A3_5*NZ+B3_5*KSV(3))**2)/CL7
C     DIFFS(6) = SQRT(C1_5(1)**2+C1_5(2)**2+C1_5(3)**2)/ANEWTON/CL2
C     DIFFS(7) = SQRT(C2(1)**2+C2(2)**2+C2(3)**2)/ANEWTON/CL4
C     DIFFS(8) = SQRT(C2_5(1)**2+C2_5(2)**2+C2_5(3)**2)/ANEWTON/CL4
C     DIFFS(9) = SQRT(SQM(1)**2+SQM(2)**2+SQM(3)**2)/ANEWTON/CL4
C     DIFFS(5) = (SPN(1,1)*LU(1)+SPN(2,1)*LU(2)+SPN(3,1)*LU(3))/
C     & SQRT(SPN(1,1)**2+SPN(2,1)**2+SPN(3,1)**2)
C     DIFFS(6) = (SPN(1,2)*LU(1)+SPN(2,2)*LU(2)+SPN(3,2)*LU(3))/
C     & SQRT(SPN(1,2)**2+SPN(2,2)**2+SPN(3,2)**2)
C     DIFFS(7) = (SPN(1,1)*SPN(1,2)+SPN(2,1)*SPN(2,2)+
C     &  SPN(3,1)*SPN(3,2))/DIFFS(3)/DIFFS(4)
C     DIFFS(5) = 180.*ACOS(DIFFS(5))/3.14
C     DIFFS(6) = 180.*ACOS(DIFFS(6))/3.14
C     DIFFS(7) = 180.*ACOS(DIFFS(7))/3.14

*     If the actual total spin is lower than the prediction,
*     it is definitely time to merge!
C     IF (DIFFS(1).LT.AFIN) THEN
C     WRITE (*,*) 'SPIN CRITERION FULFILLED, Merge in next cycle!'
*     Set some flag here to merge in next iteration...
C     END IF
      
C     DIFFS(8) = AFIN

*     *plot complete L,S1,S2 information
C     DIFFS(1) = L(1)+SPN(1,1)*BODY(1)**2/CL+SPN(1,2)*BODY(2)**2/CL
C     DIFFS(2) = L(2)+SPN(2,1)*BODY(1)**2/CL+SPN(2,2)*BODY(2)**2/CL
C     DIFFS(3) = L(3)+SPN(3,1)*BODY(1)**2/CL+SPN(3,2)*BODY(2)**2/CL
C     LABS = DIFFS(1)*DIFFS(1)+DIFFS(2)*DIFFS(2)+DIFFS(3)*DIFFS(3)
C     DIFFS(1) = DIFFS(1)/SQRT(LABS)
C     DIFFS(2) = DIFFS(2)/SQRT(LABS)
C     DIFFS(3) = DIFFS(3)/SQRT(LABS)
C     DIFFS(4) = SPN(1,1)*BODY(1)**2/CL
C     DIFFS(5) = SPN(2,1)*BODY(1)**2/CL
C     DIFFS(6) = SPN(3,1)*BODY(1)**2/CL
C     DIFFS(7) = SPN(1,2)*BODY(2)**2/CL
C     DIFFS(8) = SPN(2,2)*BODY(2)**2/CL
C     DIFFS(9) = SPN(3,2)*BODY(2)**2/CL
C     DIFFS(4) = LU(1)
C     DIFFS(5) = LU(2)
C     DIFFS(6) = LU(3)
**********************

      
C     VWHOLE = MAX(SQRT(V12),SQRT(V22))

c$$$  EKIN = M1*V12/2.0 + M2*V22/2.0
c$$$  EPOTEN = -M1*M2/RIJ
***   add spin energy
C     ESPIN = M1*M1*(SQRT(SPN(1,SPNI1)*SPN(1,SPNI1)+
C     &       SPN(2,SPNI1)*SPN(2,SPNI1)+SPN(3,SPNI1)*SPN(3,SPNI1)))+M2*
C     &       M2*(SQRT(SPN(1,SPNI2)*SPN(1,SPNI2)+
C     &       SPN(2,SPNI2)*SPN(2,SPNI2)+SPN(3,SPNI2)*SPN(3,SPNI2)))

c$$$  ENEWT = EKIN + EPOTEN
C     WRITE (*,*) 'ENERGYCHECK',TIME+TOFF, RELENERGY(IPAIR), ENEWT
C     IF (I1.EQ.1) WRITE (*,*) TIME+TOFF, DIFFS(1), DIFFS(6),
C     &   DIFFS(7), DIFFS(8)
c$$$  KSA1(1) =  M2*A(1)/MASAS
c$$$  KSA1(2) =  M2*A(2)/MASAS
c$$$  KSA1(3) =  M2*A(3)/MASAS
c$$$  
c$$$  KSA2(1) =  -M1*A(1)/MASAS
c$$$  KSA2(2) =  -M1*A(2)/MASAS
c$$$  KSA2(3) =  -M1*A(3)/MASAS
c$$$  
c$$$  POWER1 = (KSA1(1)*VI(1)+KSA1(2)*VI(2)+KSA1(3)*VI(3))
c$$$  &     *M1
c$$$  POWER2 = (KSA2(1)*VI(4)+KSA2(2)*VI(5)+KSA2(3)*VI(6))
c$$$  &     *M2

      RETURN

      END
