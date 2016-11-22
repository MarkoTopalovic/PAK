C==========================================================================
C==========================================================================
CS                      MAKSIMOVICEV MATERIJALNI MODEL
CE                      MAKSIMOVIC MATERIAL MODEL 
C==========================================================================
C==========================================================================
C     Poslednja izmena: 
C     Opis: 
C           
C==========================================================================
CE    SUBROUTINE D3M45
CE               TI3445
C
      SUBROUTINE D3M45(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC) 
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    PROGRAM FOR DEFINITION OF LOCATIONS AT INTEGRATION PIONT LEVEL
C
      include 'paka.inc'
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG 
C
      DIMENSION TAU(6),DEF(6) 
C
      IF(IDEBUG.GT.0) PRINT *, ' D3M45'
C
      LFUN=MREPER(1) 
      LNTA=MREPER(2)
      MATE=MREPER(4)
C
      LTAU=LPOCG
      LDEFT=LTAU + 6
      LDEFPP=LDEFT + 6
      LEMP=LDEFPP + 6
      LXT=LEMP + 1
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 6
      LDEFP1=LDEFT1 + 6
      LEMP1=LDEFP1 + 6
      LXTDT=LEMP1 + 1
C
      CALL TI3445(PLAST(LTAU),PLAST(LDEFT),PLAST(LDEFPP),
     1            PLAST(LEMP),PLAST(LXT),
     &            PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LEMP1),PLAS1(LXTDT), 
     &            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC)
C
      RETURN
      END
C  ========================================================================
C
      SUBROUTINE TI3445(TAUT,DEFT,DEFPP,EMP,XT,
     &                  TAU1,DEF1,DEFP1,EMP1,XTDT,
     &                  FUN,NTA,MATE,TAU,DEF,IRAC,IBTC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS    PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
CS    MAKSIMOVICEVOG MODEL 
CE    PROGRAM FOR STRESS INTEGRATION FOR MAKSIMOVIC MODEL
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     &               DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /ITERBR/ ITER
      COMMON /CONMAT/ AE,EP,DVT
      COMMON /MATERb/ korz(100,100,3),evg(100,100,3)
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAUT(6),DEFT(6),DEFPP(6),TAU(6),DEF(6),TAU1(6),DEF1(6), 
     &          DEFP1(6)
      DIMENSION FUN(11,*),NTA(*),DSIG(6),DEPS(6),DFDS(6),DGDS(6),ALAM(6)
     &         ,DDEFE(6)
      dimension di1ds(6),dj2dds(6),dfds1(6),dfds2(6),dfds3(6),dfds12(6)
     &                  ,dj3dds(6),dgds1(6),dgds2(6),dgds3(6),dgds12(6)
     &                  ,dfds34(6),dtdj3dT(6),dtdj2dT(6),dthds(6)
     &                  ,dphids(6)
C
      IF(IDEBUG.EQ.1) PRINT *, 'TI3445'
C
      IF(IRAC.EQ.2) RETURN
c==========================================================================
CE    BASIC KONSTANTS
      TOLL = 1.D-6                ! tolerance 
      tnull= 1.d-4
      MAXT = 500                  ! max. no. of iterations
      pi   = 4.d0*atan(1.d0)
      atriq=dsqrt(3.d0)
c==========================================================================
c     Material constants
      E    = FUN(1,MAT)           ! young's modulus
      ANI  = FUN(2,MAT)           ! poisson's ratio 
c
      fib  = FUN(3,MAT)           ! basic angle
      dfi  = FUN(4,MAT)           ! angle diference
      pn   =-FUN(5,MAT)           ! parameter
C==========================================================================
      fib  = fib*pi/180.d0        ! deg. to rad.
      dfi  = dfi*pi/180.d0        ! deg. to rad.
C==========================================================================
c
c     elasticity matrix for 3d
      call mel3el(elast,e,ani)    ! formiranje matrice elasticnosti
c
c     {de}
      call jedna1(deps,def,6)     ! deps(i)=def(i)
      call zbirmm(deps,deft,6)    ! deps(i)=deps(i)-deft(i)
c
c     deltaSigma_E 
      call clear(dsig,6)                 ! {dSigma}=0
      call mnozm1(dsig,elast,deps,6,1,6) !dsig(i,j)=elast(i,k)*deps(k,j)
c
c     Sigma_E
      call zbir2b(tau,taut,dsig,6)       ! tau(i)=taut(i)+dsig(i)
c
      call clear(ddefp,6)                ! {deP}=0
      if(kor.eq.1.and.iter.eq.0) goto 400
c
c     I1
      call ainvI1(ai1,tau)
c
c     Tension cutoff
      aT=0.d0
      if(ai1.gt.aT) then
        do I=1,6
            tau(I)=0.d0
        enddo
      goto 400
      endif 
c
c     I2
      call ainvI2(ai2,tau)
c
c     I3
      call ainvI3(ai3,tau)
c
c     J2D  
      call ainvJ2d(aj2d,tau)
c
c     Sqrt(J2D)
      aj2d=dsqrt(aj2d*aj2d)
      aj2dq=dsqrt(aj2d)
c
c     J3D
      call ainvJ3d(aj3d,ai1,ai2,ai3)
c
c     J3D/J2D^3/2  
      aj2d3d=aj3d/(aj2dq**3)
      if(aj2d.le.tnull) then
        aj2d3d=2.d0*atriq/9.
        aj2d=tnull
      endif
c
c     Lode's angle argument 
      alode=3.d0*atriq*aj2d3d/2.
c
      if(alode.gt. 1.d0) alode= 0.98d0
      if(alode.lt.-1.d0) alode=-0.98d0
c
c     Lode's angle (Theta)
      theta=1.d0/3.*dasin(-alode)
c
c     Increment of plastic strain is set to zerro
      call clear(DDEFP,6)
c
c     SigmaM
      sigmam=ai1/3.d0
c
c     fiM
      phim=fib+dfi/2.d0
c
c     pAV
      pav=pn/3.*(3.d0-dsin(phim))/(1.d0-(dsin(phim))**2)
c
c     fi
      phi=fib+dfi/(1.d0+sigmam/pav)
c
c     Maksimovic failure curve
      call funMAK(Fmake,ai1,aj2dq,theta,phi)
c
c     Yielding check 
c     Fmak>0
      if(Fmake.gt.toll) goto 100
c     Fmak<0
      if(Fmake.le.toll) goto 400
c
c     Not satisfied any condition
      stop 'Maksimovic passed all conditions!'
c==========================================================================
c         PLASTIC DOMAIN
c==========================================================================
  100     continue
c--------------------------------------------------------------------------  
c ******* dF/dI1 *********************************
          dfdi1= dsin(phi)/3.
c
c ******* dF/dJ2D ********************************
          dfdj2d= 1.d0/(2.d0*aj2dq)*(dcos(theta)
     &           -1.d0/atriq*dsin(theta)*dsin(phi))
c     
c ******* dF/dTheta ******************************
          dfdth=-aj2dq*(dcos(theta)*dsin(phi)/atriq+dsin(theta))
c
c ******* dF/dPhi ******************************
          dfdph= (ai1/3.-aj2dq/atriq*dsin(theta))*dcos(phi)
c
c--------------------------------------------------------------------------         
c ******* {dTheta/dSigma}T ***********************
c--------------------------------------------------------------------------
c         pomocna
          akor=4.d0/3.-(9.d0*aj3d**2)/(aj2d**3)
          akoren=dsqrt(akor)
c
c ******* dTheta/dJ2D ****************************
          dtdj2d=3.d0*aj3d/(2.d0*akoren*aj2dq**5)
          if(akor.le.toll) dtdj2d= 0.d0
c  
c ******* dTheta/dJ3D ****************************
          dtdj3d=-1.d0/(akoren*aj2dq**3)
          if(akor.le.toll) dtdj3d=0.d0
c
c         Sigma*
          tauz=(tau(1)*tau(2)+tau(2)*tau(3)+tau(3)*tau(1)
     &         -tau(4)**2-tau(5)**2-tau(6)**2)/3.
c
c         Sigma_m
          tauem=(tau(1)+tau(2)+tau(3))/3.
c
c ******* {dI1/dSigma}T **************************
          di1ds(1)=1.d0
          di1ds(2)=1.d0
          di1ds(3)=1.d0
          di1ds(4)=0.d0
          di1ds(5)=0.d0
          di1ds(6)=0.d0
c
c ******* {dJ2D/dSigma}T *************************
          dj2dds(1)=(2.d0*tau(1)      -tau(2)      -tau(3))/3.
          dj2dds(2)=(    -tau(1) +2.d0*tau(2)      -tau(3))/3.
          dj2dds(3)=(    -tau(1)      -tau(2) +2.d0*tau(3))/3.
          dj2dds(4)= 2.d0*tau(4)
          dj2dds(5)= 2.d0*tau(5)
          dj2dds(6)= 2.d0*tau(6)
c
c ******* {dJ3D/dSigma}T *************************
          dj3dds(1)= tauz-tau(5)**2+tau(2)*tau(3)+tauem*(tau(2)+tau(3))
     &              +2.d0*tauem**2
          dj3dds(2)= tauz-tau(6)**2+tau(1)*tau(3)+tauem*(tau(1)+tau(3))
     &              +2.d0*tauem**2
          dj3dds(3)= tauz-tau(4)**2+tau(1)*tau(2)+tauem*(tau(1)+tau(2))
     &              +2.d0*tauem**2
          dj3dds(4)=-2.d0*tau(3)*tau(4)+2.d0*tau(5)*tau(6)
     &              -2.d0*tau(4)*tauem
          dj3dds(5)=-2.d0*tau(1)*tau(5)+2.d0*tau(6)*tau(4)
     &              -2.d0*tau(5)*tauem
          dj3dds(6)=-2.d0*tau(2)*tau(6)+2.d0*tau(4)*tau(5)
     &              -2.d0*tau(6)*tauem
c
c ******* dTheta/dSigma (zbir) *******************
          call jednak(dtdj2dT,dj2dds,dtdj2d,6)
          call jednak(dtdj3dT,dj3dds,dtdj3d,6)
          call zbir2b(dthds,dtdj3dT,dtdj2dT,6)
c
c ******* {dphi/dSigma}T *************************
          dphids(1)=-3.d0*dfi*pav/(ai1+3.d0*pav)**2
          dphids(2)=-3.d0*dfi*pav/(ai1+3.d0*pav)**2
          dphids(3)=-3.d0*dfi*pav/(ai1+3.d0*pav)**2
          dphids(4)= 0.d0
          dphids(5)= 0.d0
          dphids(6)= 0.d0
c
c -------------------------------------------------------------------------  
c         {dF/dSigma}T
          call jednak(dfds1,di1ds,dfdi1,6)
          call jednak(dfds2,dj2dds,dfdj2d,6)
          call jednak(dfds3,dthds,dfdth,6)
          call jednak(dfds4,dphids,dfdph,6)
          call clear(dfds3,6)
c
          call zbir2b(dfds12,dfds1,dfds2,6)
          call zbir2b(dfds34,dfds3,dfds4,6)
          call zbir2b(dfds,dfds12,dfds34,6)
c
c--------------------------------------------------------------------------          
c         {dGmc/dSigma}T
          call jedna1(dgds,dfds,6)
c 
c -------------------------------------------------------------------------  
c         lambda
          call clear(alam,6)
          call mnozt1(alam,dfds,elast,6,6)      ! {dFmc/dSigma}T*[Ce]
          fcg=dot(alam,dgds,6)    ! {dFmc/dSigma}T*[Ce]*{dGmc/dSigma}
          fce=dot(alam,deps,6)    ! {dFmc/dSigma}T*[Ce]*{de}
          dlam=fce/fcg
c
C==========================================================================
c         Inicijalizacija za bisekcije
          I=0
          dlam1 = 0.d0
          dlam2 = 0.1*dlam
          Fmak1  = Fmake
          Fmax  = 1.d3*Fmake
          call jednak(ddefp,dgds,dlam2,6)  ! {deP}=dL{dGmc/dSigma}
          call oduz2b(ddefe,deps,ddefp,6)  ! {deE}={de}-{deP}
          call clear(dsig,6)               ! {dSigma}=0
          call mnozi1(dsig,elast,ddefe,6,6)! {dSigma}=[Ce]{deE}
          call zbir2b(tau,taut,dsig,6)
          call ainvI1(ai1,tau)
          call ainvJ2d(aj2d,tau)
          aj2dq=dsqrt(aj2d)
          call funMAK(Fmak2,ai1,aj2dq,theta,phi)
          call clear(ddefp,6)
          write(3,*)'dlam,dlam1,dlam2,Fmak,Fmak1,Fmak2',
     &               dlam,dlam1,dlam2,Fmak,Fmak1,Fmak2
c            
C==========================================================================
CD        LOCAL LOOP
          write(3,*)' I,      dlam,      dlam1,      dlam2,  
     &     Fmak,      Fmak1,      Fmak2,'
  110     I = I + 1
c--------------------------------------------------------------------------          
            Fp=(Fmak2-Fmak1)/(dlam2-dlam1)
            dlam=dlam2-Fmak2/Fp
            dlam1=dlam2
            dlam2=dlam
            call jednak(ddefp,dgds,dlam,6)    ! {deP}=dL{dGmc/dSigma}
            call oduz2b(ddefe,deps,ddefp,6)   ! {deE}={de}-{deP}
            call clear(dsig,6)                ! {dSigma}=0
            call mnozi1(dsig,elast,ddefe,6,6) ! {dSigma}=[Ce]{deE}
c
c           {Sigma_t+dt}={Sigma_t}+{dSigma}
            call zbir2b(tau,taut,dsig,6)
c
c           I1
            call ainvI1(ai1,tau)
c
c           J2D
            call ainvJ2d(aj2d,tau)
c
c           sqrt(J2D)
            aj2dq=dsqrt(aj2d)
c
c           SigmaM
            sigmam=ai1/3.d0
c
c           fiM
            phim=fib+dfi/2.d0
c
c           pAV
            pav=pn/3.*(3.d0-dsin(phim))/(1.d0-(dsin(phim))**2)
c
c           fi
            phi=fib+dfi/(1.d0+sigmam/pav)
C
c           Maksimovic yield curve
            call funMAK(Fmak,ai1,aj2dq,theta,phi)
            Fmak1=Fmak2
            Fmak2=Fmak
c
            write(3,1001) I,dlam,dlam1,dlam2,Fmak,Fmak1,Fmak2
c
            if(I.gt.maxt) then
                stop 'Max. num. of bisection in Maksimovic model!'
            endif
c
            if(Fmak.gt.Fmax) then
                write(3,*)'Fmak,Fmax ',Fmak,Fmax
                stop 'Solution diverges in Maksimovic model!'
            endif
c
            if(dabs(Fmak).gt.toll) goto 110
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      goto 400 
C==========================================================================
CE    UPDATES FOR NEXT STEP
  400 CONTINUE
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      demp=(ddefp(1)+ddefp(2)+ddefp(3))/3.
      emp1=emp+demp
      evp1=3.d0*emp1
C
      CALL ZBIR2B(DEFP1,DEFPP,DDEFP,6)
C
c========================================================================
c     Corection of values from previous step when convergence is reatched
      call jedna1(def1,def,6)
      call jedna1(tau1,tau,6)
c
      return
C==========================================================================
 1001 FORMAT(I3,6E12.4)
      end
C==========================================================================
C==========================================================================
      SUBROUTINE ainvI1(ai1,tau)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     first stress invariant I1
C
      dimension tau(6)
      ai1=tau(1)+tau(2)+tau(3)
      return
      end
C==========================================================================
      SUBROUTINE ainvI2(ai2,tau)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     second stress invariant I2
C
      dimension tau(6)
      ai2=-tau(1)*tau(2)-tau(2)*tau(3)-tau(3)*tau(1)
     &    +tau(4)**2+tau(5)**2+tau(6)**2
      return
      end
C==========================================================================
      SUBROUTINE ainvI3(ai3,tau)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     third stress invariant I2
C
      dimension tau(6)
      ai3= tau(1)*tau(2)*tau(3)
     &    -tau(1)*tau(5)**2-tau(2)*tau(6)**2-tau(3)*tau(4)**2
     &    +2.d0*tau(4)*tau(5)*tau(6)
      return
      end
C==========================================================================
      SUBROUTINE ainvJ2d(aj2d,tau)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     second deviatoric stress invariant J2D
C
      dimension tau(6)
      aj2d=1.d0/6.*((tau(1)-tau(2))**2  +
     &              (tau(2)-tau(3))**2  +
     &              (tau(3)-tau(1))**2) +
     &               tau(4)**2+tau(5)**2+tau(6)**2
      return
      end
C==========================================================================
      SUBROUTINE ainvJ3d(aj3d,ai1,ai2,ai3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     third deviatoric stress invariant J3D
C
      aj3d=ai3+1.d0/3.*ai1*ai2+2.d0/27.*ai1**3
      return
      end
C==========================================================================
      SUBROUTINE funMAK(F,ai1,aj2dq,theta,phi)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     Maksimovic function
C
      atriq=dsqrt(3.d0)
      F=ai1*dsin(phi)/3.+aj2dq*(dcos(theta)-dsin(theta)*dsin(phi)
     &    /atriq)
      return
      end
C==========================================================================