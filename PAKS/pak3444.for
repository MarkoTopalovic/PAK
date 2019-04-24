C =========================================================================
C =========================================================================
CS               GENERALIZOVANI HOEK-BROWN MATERIJALNI MODEL
CE               GENERALIZED HOEK-BROWN MATERIAL MODEL
C =========================================================================
C     Poslednja izmena: 16.12.2013.
C     Opis: Izmenjen deo oko granicnih vrednosti Theta
C           Dodat uslov za dlam=0 u plasticnosti
C
C =========================================================================
CE    SUBROUTINE D3M44
CE               TI3444 
C
      SUBROUTINE D3M44(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC)
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
      IF(IDEBUG.GT.0) PRINT *, ' D3M44'
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
      CALL TI3444(PLAST(LTAU),PLAST(LDEFT),PLAST(LDEFPP),
     1            PLAST(LEMP),PLAST(LXT),
     &            PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LEMP1),PLAS1(LXTDT), 
     &            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC)
C
      RETURN
      END
C
C  ========================================================================
C
      SUBROUTINE TI3444(TAUT,DEFT,DEFPP,EMP,XT,
     &                  TAU1,DEF1,DEFP1,EMP1,XTDT,
     &                  FUN,NTA,MATE,TAU,DEF,IRAC,IBTC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS    PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
CS    HOEK-BROWN MODEL 
CE    PROGRAM FOR STRESS INTEGRATION FOR HOEK-BROWN MODEL
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
     &         ,DDEFE(6),dsiga(6),dsigb(6),dsigc(6)
      dimension di1ds(6),dj2dds(6),dthds(6),dj3dds(6),dgds1(6)
     &         ,dgds2(6),dtds1(6),dtds2(6),dfds1(6),dfds2(6)
     &         ,taua(6),taub(6),tauc(6),ddefea(6),ddefeb(6),ddefec(6)
     &         ,ddefpa(6),ddefpb(6),ddefpc(6),tulz(6)
C
      IF(IDEBUG.EQ.1) PRINT *, 'TI3444'
C
      IF(IRAC.EQ.2) RETURN
c
CE    BASIC KONSTANTS
      TOLL =  1.D-6               ! tolerance
      toldl=  1.d-12
      tolq=  1.d-10
      MAXT =  500                 ! max. no. of iterations
      atriq=  dsqrt(3.d0)
      pi   =  4.d0*atan(1.d0)
c =========================================================================
c     Material constants
      E    = FUN(1,MAT)           ! young's modulus
      ANI  = FUN(2,MAT)           ! poisson's ratio 
c
      smc  = FUN(3,MAT)           ! material constant
      em   = FUN(4,MAT)           ! material constant 
      emd  = FUN(5,MAT)           ! material constant
      es   = FUN(6,MAT)           ! material constant
      ahb  = FUN(7,MAT)           ! material constant
c =========================================================================
c     ELASTIC DOMAIN
c =========================================================================
c     elasticity matrix for 3d
      call mel3el(elast,e,ani)    ! formiranje matrice elasticnosti
c
c     {de}
      call jedna1(deps,def,6)     ! deps(i)=def(i)
      call zbirmm(deps,deft,6)    ! deps(i)=deps(i)-deft(i) 
c
c     deltaSigma_E
      call clear(dsig,6)
      call mnozi1(dsig,elast,deps,6,6) !dsig(i)=elast(i,k)*deps(k)
c
c     Sigma_E
      call zbir2b(tau,taut,dsig,6)       ! tau(i)=taut(i)+dsig(i)
      call jedna1(tulz,tau,6)
c
c     Prva invarijanta napona I1
      ai1=tau(1)+tau(2)+tau(3)           ! i1=sigma1+sigma2+sigma3
      ai1e=ai1
c
c     Tension cutoff
      sigt=es*smc/em
      ai1t=3.0d0*sigt
      if(ai1.ge.ai1t) then
c      write(3,*)'tension cutoff 1, sigt, ai1t, ai1', 
c     &                             sigt, ai1t, ai1
c      CALL WRR6 (TAU,6,'TAU  ')
        do i=1,6
          IF(I.LE.3)THEN
            tau(i)=0.98d0*sigt
          ELSE
            tau(i)=0.d0
          ENDIF
        enddo
        go to 400
      endif 
c
c     Druga invarijanta napona I2
      ai2= tau(1)*tau(2)+tau(2)*tau(3)+tau(3)*tau(1)
     &    -tau(4)**2-tau(5)**2-tau(6)**2
c
c     Treca invarijanta napona I3
      ai3= tau(1)*tau(2)*tau(3)
     &    -tau(1)*tau(5)**2-tau(2)*tau(6)**2-tau(3)*tau(4)**2
     &    +2.d0*tau(4)*tau(5)*tau(6)
c
c     Druga invarijanta devijatora napona J2D
      aj2d=1.d0/6.d0*((tau(1)-tau(2))**2  +
     &                (tau(2)-tau(3))**2  +
     &                (tau(3)-tau(1))**2) +
     &                 tau(4)**2+tau(5)**2+tau(6)**2
c
c     Treca invarijanta devijatora napona J3D
      aj3d=ai3-1.d0/3.d0*ai1*ai2+2.d0/27.d0*ai1**3
c
c     Sqrt(J2D)
      AJ2D=DSQRT(AJ2D*AJ2D)
      aj2dq=dsqrt(aj2d)
      aj2dqe=aj2dq
c
c     J3D/J2D^3/2  
      if(dabs(aj2dq).lt.tolq) then
        aj2d3d=0.d0
      else
        aj2d3d=aj3d/(aj2dq**3)
      endif
c
c     Lode's angle argument 
      alode=-3.d0*atriq/2.d0*aj2d3d
      if(alode.gt. 3.98d0) then
           alode= 0.98d0
      endif
      if(alode.lt.-3.98d0) then
           alode=-0.98d0
      endif
      if(alode.gt. 0.98d0) alode= 1.98d0-alode
      if(alode.lt.-0.98d0) alode=-1.98d0-alode
      if(alode.gt. 1.98d0) alode=-1.98d0+alode
      if(alode.lt.-1.98d0) alode= 1.98d0+alode
      if(alode.gt. 2.98d0) alode= 3.98d0-alode
      if(alode.lt.-2.98d0) alode=-3.98d0-alode
c
      if(alode.lt.0.0d0) alode=-alode
c     Lode's angle (Theta)
      theta=1.d0/3.d0*dasin(alode)
c      zero=0.d0
c      if(theta.lt.zero) theta = 0.d0
      thetad=theta*180.d0/pi
c
c     Increment of plastic strain is zerro in elastic domain
      call clear(DDEFP,6)
c
      demp=0.d0
c
c     Generalised Hoek-Brown yield surface
      call funHB(Fhb,ai1,aj2dq,smc,em,es,ahb,theta)
c      write(3,*)'Fhb,ai1,aj2dq,aj3d,theta',Fhb,ai1,aj2dq,aj3d,theta
c
      Fhbe=Fhb
c
C =======================================================================
c     Yielding check  
c     Fph>0
      if(Fhbe.gt.toll) goto 100
c     Fph<0
      if(Fhbe.le.toll) goto 400
c
c     U slucaju prolaska svih uslova
      stop 'ERROR! Generalized Hoek-Brown passed all conditions!!!'
c =========================================================================
c     PLASTIC DOMAIN
c ==========================================================================
  100   continue
c ------------------------------------------------      
c ***** {dG/dSigma}T *****************************
c ------------------------------------------------
c ***** dG/dI1 ***********************************
        dgdi1=emd/3.d0*(smc**(1.d0/ahb-1.d0))
c
c ***** dG/dJ2D **********************************
        dgdj2d =1.d0/aj2dq*(2.d0**(1./ahb-1.)*dcos(theta)*(dcos(theta)*
     &           aj2dq)**(1./ahb-1.)/ahb+(dcos(theta)-dsin(theta)
     &           /atriq)*emd*smc**(1./ahb-1.)/2.d0)
        if(aj2dq.lt.tolq) dgdj2d=0.d0
c
c ***** {dI1/dSigma}T ****************************
        di1ds(1)=1.d0
        di1ds(2)=1.d0
        di1ds(3)=1.d0
        di1ds(4)=0.d0
        di1ds(5)=0.d0
        di1ds(6)=0.d0
c
c ***** {dJ2D/dSigma}T ***************************
        dj2dds(1)=(2.d0*tau(1)      -tau(2)      -tau(3))/3.d0
        dj2dds(2)=(    -tau(1) +2.d0*tau(2)      -tau(3))/3.d0
        dj2dds(3)=(    -tau(1)      -tau(2) +2.d0*tau(3))/3.d0
        dj2dds(4)=2.d0*tau(4)
        dj2dds(5)=2.d0*tau(5)
        dj2dds(6)=2.d0*tau(6)
c
c ***** {dG/dSigma}T (suma)***********************
        call jednak(dgds1,di1ds,dgdi1,6)
        call jednak(dgds2,dj2dds,dgdj2d,6)
        call zbir2b(dgds,dgds1,dgds2,6)
c
c -------------------------------------------------         
c ***** {dF/dSigma}T *****************************
c -------------------------------------------------
c ***** dF/dI1 ***********************************
        dfdi1=(em*smc**(1.d0/ahb-1.d0))/3.d0
c
c ***** dF/dJ2D **********************************
!        if(aj2dq.le.toll) stop 'aj2dq<0!, 2'
        dfdj2d =2.d0**(1.d0/ahb-1.d0)*dcos(theta)*(dcos(theta)*
     &          aj2dq)**(1.d0/ahb-1.d0)/ahb/aj2dq+
     &         (dcos(theta)-dsin(theta)/atriq)*em*smc**(1.d0/ahb-1.d0)/
     &          2.d0/aj2dq
        if(aj2dq.le.toll) dfdj2d=0.d0
c
c ***** {dF/dSigma}T (sum)************************
        call jednak(dfds1,di1ds,dfdi1,6)
        call jednak(dfds2,dj2dds,dfdj2d,6)
c
        call zbir2b(dfds,dfds1,dfds2,6)
c
c --------------------------------------------------------------------------  
c ***** dLambda **********************************
        call clear(alam,6)
        call mnozt1(alam,dfds,elast,6,6)      ! {dF/dSigma}T*[Ce]
        call clear(fcg,1)
        fcg=dot(alam,dgds,6)    ! {dF/dSigma}T*[Ce]*{dG/dSigma}
        call clear(fce,1)
        fce=dot(alam,deps,6)    ! {dF/dSigma}T*[Ce]*{de}
        dlam=fce/fcg
c
        dlamp=dlam
        if(dabs(dlam).lt.toldl) then
            dlam=toldl
c            goto 400
        endif
c
c -------------------------------------------------------------------------  
cr      Inicijalizacija za bisekcije
c -------------------------------------------------------------------------  
       I=0
       dlama = 0.d0
       dlamb = 0.25d0*dlam
       dlamc = 0.50d0*dlam
       Fmax  = 1.0d3*Fhbe
       call jednak(ddefpa,dgds,dlama,6)   ! {deP}=dL{dGmc/dSigma}
       call jednak(ddefpb,dgds,dlamb,6)   ! {deP}=dL{dGmc/dSigma}
       call jednak(ddefpc,dgds,dlamc,6)   ! {deP}=dL{dGmc/dSigma}
       call oduz2b(ddefea,deps,ddefpa,6)  ! {deE}={de}-{deP}
       call oduz2b(ddefeb,deps,ddefpb,6)  ! {deE}={de}-{deP}
       call oduz2b(ddefec,deps,ddefpc,6)  ! {deE}={de}-{deP}
       call clear(dsiga,6)                ! {dSigma}=0
       call clear(dsigb,6)                ! {dSigma}=0
       call clear(dsigc,6)                ! {dSigma}=0
       call mnozi1(dsiga,elast,ddefea,6,6)! {dSigma}=[Ce]{deE}
       call mnozi1(dsigb,elast,ddefeb,6,6)! {dSigma}=[Ce]{deE}
       call mnozi1(dsigc,elast,ddefec,6,6)! {dSigma}=[Ce]{deE}
       call zbir2b(taua,taut,dsiga,6)
       call zbir2b(taub,taut,dsigb,6)
       call zbir2b(tauc,taut,dsigc,6)
       call ainvI1(ai1a,taua)
       call ainvI1(ai1b,taub)
       call ainvI1(ai1c,tauc)
       call ainvJ2d(aj2da,taua)
       call ainvJ2d(aj2db,taub) 
       call ainvJ2d(aj2dc,tauc)
       aj2dqa = dsqrt(aj2da)
       aj2dqb=dsqrt(aj2db)
       aj2dqc=dsqrt(aj2dc) 
       call funHB(Fhba,ai1a,aj2dqa,smc,em,es,ahb,theta)
       call funHB(Fhbb,ai1b,aj2dqb,smc,em,es,ahb,theta) 
       call funHB(Fhbc,ai1c,aj2dqc,smc,em,es,ahb,theta)
       dlam = 0.0d0 
C =================================================
CD      LOCAL ITERATIONS                            
c        write(3,*)' I,      dlam,      dlamp,      dlamm,  
c     &      Fhb,       Fhbp,       Fhbm,kor'
  110 I = I + 1
c -------------------------------------------------   
c      if (kor.eq.2.and.iter.eq.1) then
c          write(*,*) 'korak 2, iteracija 1'
c      endif
c -------------------------------------------------          
            dlam=dlamb+(Fhbb/Fhba*(Fhba/Fhbc*(Fhbb/Fhbc-Fhba/Fhbc)*
     &          (dlamc-dlamb)-(1.d0-Fhbb/Fhbc)*(dlamb-dlama)))/
     &          ((Fhba/Fhbc-1.d0)*(Fhbb/Fhbc-1.d0)*(Fhbb/Fhba-1.d0))
c
            if(dabs(dlam).gt.10) then
                write(*,*) 'dlam= ',dlam
                write(*,*) 'dlama,dlamb,dlamc ',dlama,dlamb,dlamc
                write(*,*) 'Fhba,Fhbb,Fhbc ',Fhba,Fhbb,Fhbc
              stop
            endif
c
            dlama=dlamb    
            dlamb=dlamc
            dlamc=dlam
c
            call jednak(ddefp,dgds,dlam,6)   ! {deP}=dL{dGmc/dSigma}
            call oduz2b(ddefe,deps,ddefp,6)  ! {deE}={de}-{deP}
            call clear(dsig,6)
            call mnozi1(dsig,elast,ddefe,6,6)! {dSigma}=[Ce]{deE}
c
c           {Sigma_t+dt}={Sigma_t}+{dSigma}
            call zbir2b(tau,taut,dsig,6)
c
c           Prva invarijanta napona
            ai1=tau(1)+tau(2)+tau(3) 
c
            if(ai1.ge.ai1t) then
c            write(3,*)'tension cutoff 2'
              do i=1,6
                IF(I.LE.3)THEN
                  tau(i) = sigt
                ELSE
                  tau(i) = 0.d0
                ENDIF
              enddo
              goto 400
	      endif 
c           J2D
            aj2d=1.d0/6.d0*((tau(1)-tau(2))**2 +
     &                      (tau(2)-tau(3))**2 +
     &                      (tau(3)-tau(1))**2)+
     &                       tau(4)**2+tau(5)**2+tau(6)**2
c            if(dabs(aj2d).lt.toll) aj2d=toll
c
c           sqrt(J2D)
            aj2dq=dsqrt(aj2d)
c
c           Generalised Hoek-Brown yield curve
            call funHB(Fhb,ai1,aj2dq,smc,em,es,ahb,theta)
            Fhba=Fhbb
            Fhbb=Fhbc
            Fhbc=Fhb
c 
            if(I.gt.maxt) then
            stop 'Max. num. of bisection in Gen. Hoek-Brown model!'
            endif
c
c            if(kor.eq.2.and.ITER.eq.1) then
c              write(3,1000) kor,i,theta
c              write(3,1001) Fhbe,ai1e,aj2dqe  
c              write(3,1002) Fhb,ai1,aj2dq  
c            endif
 1000 FORMAT('kor,i,theta     ',2i18,f18.4)
 1001 FORMAT('Fhbe,ai1e,aj2dqe',3f18.4)
 1002 FORMAT('Fhb,ai1,aj2dq   ',3f18.4)
c            
            if(Fhb.gt.Fmax) then
               write(3,*)'Fhb,Fhbe,Fmax,ai1e,aj2dqe',
     &                    Fhb,Fhbe,Fmax,ai1e,aj2dqe
                CALL WRR6 (tulz,6,'tulz ')
                CALL WRR6 (TAU,6,'TAU  ')
                CALL WRR6 (dsig,6,'dsig ')
                write(3,*)'E, ANI, NE',E, ANI, NE
                write(3,*)'ai1, aj2d, dlam',ai1, aj2d, dlam
               stop 'Solution diverges in Gen. Hoek-Brown model!'
            endif          
c
!            tollam=dlam/dlamp
            tollam=dlamb-dlamc
c
      if(dabs(Fhb).gt.toll.and.dabs(tollam).gt.toldl) goto 110
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      demp=(ddefp(1)+ddefp(2)+ddefp(3))/3.d0
      emp1=emp+demp
      evp1=3*emp1 
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      goto 400 
C ========================================================================= 
! 1001 FORMAT(I3,6E12.4,I3)
CE    UPDATES FOR NEXT STEP
  400 CONTINUE
C 
C =======================================================================
c     Remaining payload (Racunanje rastojanja od loma)
      ass=aj2dq
      atriq=  dsqrt(3.d0) 
      asm= 1./(24*dcos(theta)*dcos(theta))*((-3.*dcos(theta)+atriq*
     &     dsin(theta))*em*smc+atriq*dsqrt(smc*(-16.*dcos(theta)*
     &     dcos(theta)*ai1*em+48.*dcos(theta)*dcos(theta)*es*smc+
     &     (3.*dcos(theta)*dcos(theta)-2.*atriq*dcos(theta)*
     &     dsin(theta)+dsin(theta)*dsin(theta))*em*em*smc)))
      sfactor=(1.d0-ass/asm)*100
      if(sfactor.gt.100.d0) sfactor=100.d0
      if(sfactor.lt.0.d0) sfactor=0.d0
      xtdt=sfactor
C
c =======================================================================
c     Corection of values from previous step when convergence is reatched
      CALL ZBIR2B(DEFP1,DEFPP,DDEFP,6)
      call jedna1(def1,def,6)
      call jedna1(tau1,tau,6) 
      return   
      end
C ======================================================================
      SUBROUTINE funHB(F,ai1,aj2dq,smc,em,es,ahb,theta)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     Hoek-Brown yield function
C
      atriq=  dsqrt(3.d0) 
c
      F =1.d0/3.d0*ai1*em*smc**(1./ahb-1.)-es*smc**(1./ahb)+
     &     2**(1./ahb)*(aj2dq*dcos(theta))**(1./ahb)+
     &     em*aj2dq*smc**(1./ahb-1.)*(dcos(theta)-dsin(theta)
     &     /atriq)
      return
      end
C ======================================================================= 