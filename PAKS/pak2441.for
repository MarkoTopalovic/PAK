C=======================================================================
C=======================================================================
CS    Drucker-Prager materijalni model sa kapom 2D ELEMENT
CE    Drucker-Prager cap podel 2D ELEMENT
C=======================================================================
C=======================================================================
CE    SUBROUTINE D2M41
CE               TI2441
C
      SUBROUTINE D2M41(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC)
c      USE PLAST3D
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
      DIMENSION TAU(4),DEF(4)
C
      IF(IDEBUG.GT.0) PRINT *, ' D2M41'
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MATE=MREPER(4)
C
      LTAU=LPOCG
      LDEFT=LTAU + 4*IDVA
      LDEFPP=LDEFT + 4*IDVA
      LEMP=LDEFPP + 4*IDVA
      LXT=LEMP + 1*IDVA
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 4*IDVA
      LDEFP1=LDEFT1 + 4*IDVA
      LEMP1=LDEFP1 + 4*IDVA
      LXTDT=LEMP1 + 1*IDVA
C
      CALL TI2441(A(LTAU),A(LDEFT),A(LDEFPP),
     1            A(LEMP),A(LXT),
     1            A(LTAU1),A(LDEFT1),A(LDEFP1),
     1            A(LEMP1),A(LXTDT), 
     1            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC)
C
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI2441(TAUT,DEFT,DEFPP,EMP,XT,
     1                  TAU1,DEF1,DEFP1,EMP1,XTDT,
     1                  FUN,NTA,MATE,TAU,DEF,IRAC,IBTC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS    PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
CS    DRUCKER-PRAGER MODEL 
CE    PROGRAM FOR STRESS INTEGRATION FOR DRUCKER-PRAGER MODEL
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /TAUD2/ TAUD(4),DEFDPR(4),DEFDS(4),DDEFP(4),
     1               DETAU(4),DDEF(4)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /ITERBR/ ITER
      COMMON /CONMAT/ AE,EP,DVT
      COMMON /MATERb/ korz(100,100,3),evg(100,100,3)
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAUT(4),DEFT(4),DEFPP(4),TAU(4),DEF(4),TAU1(4),DEF1(4), 
     +          DEFP1(4)
      DIMENSION FUN(11,*),NTA(*),DSIG(4),DEPS(4),DFDS(4),DGDS(4),ALAM(4)
     +         ,CP(4,4),POM(4), CEP(4,4),DDEFE(4),dsigp(4)
     +         ,DFdDS(4),DFcDS(4),DGdDS(4),DGcDS(4),ALAMd(4),ALAMc(4)
     +         ,al1p(4),al2p(4),ddefpd(4),ddefpc(4),taue(4)
c
      IF(IDEBUG.EQ.1) PRINT *, 'TI2441'
c
      IF(IRAC.EQ.2) RETURN
c     basic constants
      toll=  1.d-8
      tol =  1.d-6
      tole=  1.d-12
      maxt=  500
c====================================================================
c
c     material constants
      E       = FUN(1,MAT) 
      ani     = FUN(2,MAT)
c 
      ak      = FUN(3,MAT)
      alf     = FUN(6,MAT)
      T       = FUN(7,MAT)
      x0      =-FUN(8,MAT)
c
      W       =-FUN(9,MAT)
      D       =-FUN(10,MAT)
c====================================================================
c
c     elasticity matrix for 3d
      call MEL2el(elast,e,ani)
c      call wrr6(elast,36,'elas')
c
c     de
c      call wrr6(def,4,'def ')
      call jedna1(deps,def,4)     ! deps(i)=def(i)
c      call wrr6(deft,4,'deft')
      call zbirmm(deps,deft,4)    ! deps(i)=deps(i)-deft(i) 
c      call wrr6(deps,4,'deps')
c
c     dsigEL
      call clear(dsig,4)
      call mnozm1(dsig,elast,deps,6,1,6) ! dsig(i,j)=elast(i,k)*deps(k,j)
c
c      call wrr6(dsig,4,'dsig')
c      stop 'stop'
c
c     sigEL
      call zbir2b(tau,tauT,dsig,4)       ! tau(i)=taut(i)+dsig(i)
c
c     I1
      ai1e=tau(1)+tau(2)+tau(4)          ! I1=sigma1+sigma2+sigma3
c
c     J2de
      aj2de=1.d0/6.d0*((tau(1)-tau(2))*(tau(1)-tau(2))+
     +                 (tau(2)-tau(4))*(tau(2)-tau(4))+
     +                 (tau(4)-tau(1))*(tau(4)-tau(1)))+
     +                  tau(3)*tau(3)
c
c     sqrt(j2de)
      aj2de=dsqrt(aj2de*aj2de)
      aj2dqe=dsqrt(aj2de)
c
c     tension cutoff limit
      if(ai1e.gt.T) then
          do I=1,4
              tau(I)=T/3.d0
          enddo
          tau(3)=0
          goto 400
      endif 
c
c     increment of plastic strain is zerro in elastic domain
      ddefp(1)=0.d0
      ddefp(2)=0.d0
      ddefp(3)=0.d0
      ddefp(4)=0.d0
c
      demp=0.d0
c
c     Drucker-Prager failure surface
      fdpe=alf*ai1e+aj2dqe-ak
      fce=xt-ai1e
	Iobl=0
c
      ai1=ai1e
	call jedna1(taue,tau,4)
c
c     domain identification
      if(fdpe.gt.toll) goto 100
      goto 400
      stop 'Drucker-Prager passed all conditions!!!'
c========================================================================
c         elastic-plastic domain
c
c         dF/dSgima
  100     iobl=1
          dfds(1)=alf+( 2.d0*tau(1)-tau(2)-tau(4))/(6.d0*aj2dqe)
          dfds(2)=alf+(-tau(1)+2.d0*tau(2)-tau(4))/(6.d0*aj2dqe) 
          dfds(3)=tau(3)/aj2dqe
          dfds(4)=alf+(-tau(1)-tau(2)+2.d0*tau(4))/(6.d0*aj2dqe)
c
c         dG/dSigma = dF/dSigma
          call jedna1(dgds,dfds,4)          ! dgds(i)=dfds(i)
c
c         labdaA
          call clear(alam,4)
          call mnozt1(alam,dfds,elast,6,6)  ! alam(i)=dfds(j)*elast(j,i)
          fcg=dot(alam,dgds,4)
          dlam=dot(alam,deps,4)
          dlam=dlam/fcg
c
c         initialization
          call jednak(ddefp,dgds,dlam,4)
          call oduz2b(ddefe,deps,ddefp,4)
          call clear (dsig,4)
          call mnozi1(dsig,elast,ddefe,6,6)
c
c         sigma
          call zbir2b(tau,taut,dsig,4)      ! tau(i)=taut(i)+dsig(i)
c
c         i1
          ai1=tau(1)+tau(2)+tau(4)        
c
c         j2d
          aj2d=1.d0/6.d0*((tau(1)-tau(2))**2 +
     &                    (tau(2)-tau(4))**2 +
     &                    (tau(4)-tau(1))**2)+
     &                     tau(3)**2
c         sqrt(j2d)
          aj2d=dsqrt(aj2d*aj2d)
          aj2dq=dsqrt(aj2d)
          fdpm=alf*ai1+aj2dq-ak

c
          I=0
          dlamp=0.d0
          dlamm=dlam
          dx=0.001*dlamm
          fdp=fdpm
          fdpp=fdpe
c          dlam=0.d0
          af=2.d0
          ib=0
          jp=2
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c          write(3,*)' I,      dlam,      dlamp,      dlamm,  
c     &      Fdp,       Fdpp,       Fdpm,'
  110     I = I + 1
              call bisectg (dlam,dlamm,dlamp,dx,fdp,fdpm,fdpp,af,ib,jp)
              call jednak(ddefp,dgds,dlam,4)
              call oduz2b(ddefe,deps,ddefp,4)
              call clear (dsig,4)
              call mnozi1(dsig,elast,ddefe,6,6)
c
c             sigma
              call zbir2b(tau,taut,dsig,4)      ! tau(i)=taut(i)+dsig(i)

c             I1
              ai1=tau(1)+tau(2)+tau(4)        
c
c             j2d
              aj2d=1.d0/6.d0*((tau(1)-tau(2))**2 +
     &                        (tau(2)-tau(4))**2 +
     &                        (tau(4)-tau(1))**2)+
     &                         tau(3)**2
c             sqrt(j2d)
              aj2d=dsqrt(aj2d*aj2d)
              aj2dQ=dsqrt(aj2d)
c
              demp=(ddefp(1)+ddefp(2)+ddefp(4))/3.d0
              emp1=emp+demp
              evp1=3.d0*emp1
c
c              if(evp1.gt.(0.999*w)) then
c                  xtdt=x0-dlog(1.d0-evp1/w)/d
c              else
c                  xtdt=ai1e
c              endif
c
c              if(xtdt.gt.x0) xtdt=x0

c             Drucker-Prager failure surface
              fdp=alf*ai1+aj2dq-ak
c             write(3,1001) I,dlam,dlamp,dlamm,Fdp,Fdpp,Fdpm
c
              if(I.gt.maxt) then                                      
                  stop 'Max. num. of iteration in D-P model I'
              endif
          if(dabs(fdp).gt.tol) goto 110
c
      goto 400 
c========================================================================
c
c     updates for next step
  400 continue
c      write(3,*) 'kor,Iobl',kor,Iobl
c
      call zbir2b(defp1,defpp,ddefp,4)
c      write(3,*) 'Iobl,fc,fd',Iobl,fc,fd
c========================================================================
c     Corection of values from previous step when convergence is
c     reatched
      call jedna1(def1,def,4)
      call jedna1(tau1,tau,4)
c     call wrr6(def1,6,'def1')
      return   
      end
c=======================================================================
C=======================================================================