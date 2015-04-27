      subroutine ksrel_output(time,xi,vi,n1,n2,m1,m2,dt,np,hi,gi_pn,
     &     gi_nopn)
*
*     Output for PN events
*     ------------------------------------
*
      include 'params.h'
      include 'postnewton.h'
      REAL*8 time,XI(6),VI(6),M1,M2,DT,HI,GI_PN,GI_NOPN
      INTEGER N1,N2,NP
      REAL*8 XX12(3),VV12(3),V1_2,V2_2,V1_22,MU,M12,ETA,RIJ,ETEMP
      REAL*8 LANG(3),L2,LPN1,L2PN,ARAD,ER2,POTMIN,POTMINPN
      REAL*8 SEMI,SEMIPN,ECC2,ECC,ECCPN,GE,TZ,CL2,PI,A_EIN
      
      CL2 = CLIGHTPN*CLIGHTPN
      PI = 3.141592654

      XX12(1) = XI(1)-XI(4)
      XX12(2) = XI(2)-XI(5)
      XX12(3) = XI(3)-XI(6)
      VV12(1) = VI(1)-VI(4)
      VV12(2) = VI(2)-VI(5)
      VV12(3) = VI(3)-VI(6)         
      V1_2 = VI(1)**2+VI(2)**2+VI(3)**2
      V2_2 = VI(4)**2+VI(5)**2+VI(6)**2
      V1_22 = VV12(1)**2+VV12(2)**2+VV12(3)**2
      M12 = M1+M2
      MU = M1*M2/M12
      ETA = MU/M12
      RIJ = SQRT(XX12(1)**2+XX12(2)**2+XX12(3)**2)

      ETEMP = (0.5*V1_22 - M12/RIJ)
C      ETEMP = HI
      LANG(1) = (XX12(2)*VV12(3)-XX12(3)*VV12(2))
      LANG(2) = (XX12(3)*VV12(1)-XX12(1)*VV12(3))
      LANG(3) = (XX12(1)*VV12(2)-XX12(2)*VV12(1)) 
      L2 = (LANG(1)**2+LANG(2)**2+LANG(3)**2) !/(MU)**2
      LPN1 = 0.5*(1.0-3.0*ETA)*V1_22 + (3.0+ETA)*M12/RIJ
      L2PN = L2*(1.0+LPN1/CL2)**2
*
C     now follows an expression from memmesheimer gopakumar 2004 for a and e
      ARAD=M12/(-2.0*ETEMP)*(1.0-2.0*ETEMP/(4.0*CL2)*(-7.0+ETA)
     &     +(-2.0*ETEMP)**2/(16.0*CL2*CL2)*(1.0+10.0*ETA+ETA*ETA-1.0/
     &     (2.0*ETEMP*L2)*(-68.0+44.0*ETA)))
      ER2 = 1.0+2.0*ETEMP*L2-2.0*ETEMP/(4.0*CL2)*
     &     (24.0-4.0*ETA+5.0*(-3.0+ETA)*(-2.0*ETEMP*L2))
     &     + (-2.0*ETEMP)**2/(8.0*CL2*CL2)*
     &     (52.0+2.0*ETA+2.0*ETA*ETA
     &     -(80.0-55.0*ETA+4.0*ETA*ETA)*(-2.0*ETEMP*L2)
     &     +8.0/(2.0*ETEMP*L2)*(-17.0+11.0*ETA))
*
C     IF ((ETEMP+DIFFS(4)/MU).LT.0.) THEN
      POTMIN = -(M12*M12/(2.*L2))
      POTMINPN = -(M12*M12/(2.*L2PN))
      IF ((ETEMP.LT.0.).AND.((ETEMP.GE.POTMIN)
     &     .OR.(ETEMP.GE.POTMINPN))) THEN ! adding restriction over ETEMP after AND
         SEMI = -M12/(2.*ETEMP)
         SEMIPN = -M12/(2.*(ETEMP+DIFFS(4)/MU))
         ECC2 = 1.+(2.*ETEMP*L2)/(M12*M12)
         ECC = SQRT(ECC2)
         ECCPN = SQRT(1.+(2.*(ETEMP+DIFFS(4)/MU)*L2PN)/
     &        (M12*M12))
*     Einstein shift
         A_EIN = 6.0*PI*M12/(SEMI*CL2*(1-ECC2))
*     Merging timescale estimation
         GE = (1.0 - ECC2)**3.5
     &        /(1.0 + (73.0/24.0 + 37.0*ECC2/96.0)*ECC2)
         TZ = 5.0/64.0*CLIGHTPN**5*GE*SEMI**4/(M1*M2*M12)
         WRITE (50,*) TIME,TZ,V1_22/CLIGHTPN,A_EIN,GI_PN,GI_NOPN,DT,NP, 
     &        N1,N2,M1,M2,SEMIPN,ECCPN,SEMI,ECC,
     &        XI(1:6),VI(1:6),XX12(1:3),VV12(1:3),
     &        HI*MU,PREVRELENERGY,DIFFS(4),DIFFS(5)
      END IF

      call flush(50)

      RETURN

      END
      
