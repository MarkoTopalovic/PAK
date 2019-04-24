C=======================================================================
C=======================================================================
CS                      CAM-CLAY MATERIJALNI MODEL
CE                      CAM-CLAY MATERIAL MODEL
C=======================================================================
C=======================================================================
C     Poslednja izmena: 15.06.2016. 
C     Opis: Udaljenje od uslova loma
C           
C           
C=======================================================================
C    SUBROUTINE D3M9
C               TAUI35
C               TEQBI3
C               PRILA3 
C               DEVEQ3
C
      SUBROUTINE D3M9(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC,lpoc0)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
C
CE    PROGRAM FOR DEFINITION OF LOCATIONS AT INTEGRATION PIONT LEVEL
C
      include 'paka.inc'
C 
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C 
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
C
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /OBNOVA/ IRUSI
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAU(6),DEF(6)
C
      IF(IDEBUG.GT.0) PRINT *, ' D3M9'
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MATE=MREPER(4)
C
      LTAU=LPOCG
      LDEFT=LTAU + 6
      LDEFPP=LDEFT + 6
      LPOT=LDEFPP + 6
      LEMP=LPOT + 1
      LOCR=LEMP + 1
      LXT=LOCR+1
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 6
      LDEFP1=LDEFT1 + 6
      LPOT1=LDEFP1 + 6
      LEMP1=LPOT1 + 1
      LOCR1=LEMP1 + 1
      LXTDT=LOCR1+1
C
      IF(IRUSI.EQ.1) THEN
      LTAU0=LPOC0
      LDEFT0=LTAU0 + 6
      LDEFP0=LDEFT0 + 6
      LPOT0=LDEFP0 + 6
      LEMP0=LPOT0 + 1
      LOCR0=LEMP0 + 1
      LXT0=LOCR0+1
      ELSE
      LTAU0=1
      LDEFT0=LTAU0 + 6
      ALLOCATE (PLAS0(12))
      ENDIF

C
      CALL TI3409(PLAST(LTAU),PLAST(LDEFT),PLAST(LDEFPP),PLAST(LPOT),
     &            PLAST(LEMP),PLAST(LOCR),PLAST(LXT),
     &            PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),PLAS1(LPOT1),
     &            PLAS1(LEMP1),PLAS1(LOCR1),PLAS1(LXTDT),
     &            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC,
     &            PLAS0(LTAU0),PLAS0(LDEFT0))
C
      IF(IRUSI.NE.1) DEALLOCATE(PLAS0)
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI3409(TAUT,DEFT,DEFPP,P0T,
     &                  EMP,OCR,XT,
     &                  TAU1,DEF1,DEFP1,P0TDT,
     &                  EMP1,OCR1,XTDT,
     &                  FUN,NTA,MATE,TAU,DEF,IRAC,IBTC,
     &                  tau0,def0)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE    PROGRAM FOR STRESS INTEGRATION FOR CAM-CLAY CAP MODEL
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
     &          DEFP1(6),tau0(6),def0(6),TAUD0(6),TAUDP(6)
      DIMENSION FUN(9,MATE),NTA(*),
     &          di1ds(6),dj2dds(6),deps(6),dsig(6),dfds1(6),dfds2(6),
     &          dfds(6),dqds(6),dgds(6),alam(6),ddefe(6),ajed(6),
     &          dgds1(6),dgds2(6),dtaup(6),taup(6),
     &          ddefpa(6),ddefea(6),dsiga(6),taua(6),
     &          ddefpb(6),ddefeb(6),dsigb(6),taub(6),
     &          ddefpc(6),ddefec(6),dsigc(6),tauc(6),dfdsn(6)
C
      IF(IDEBUG.EQ.1) PRINT *, 'TI3409'
C
      IF(IRAC.EQ.2) RETURN
c=======================================================================
CE    BASIC KONSTANTS
      TOLL = 1.D-8                 ! tolerance 
      TOLI = 1.D-6
      TOLDL= 1.d-10
      MAXT = 500                  ! max. no. of iterations
      pi   = 4.d0*atan(1.d0)
      advaq=dsqrt(2.0D0)
      atriq=dsqrt(3.d0)
c=======================================================================
CE    MATERIAL CONSTANTS
      E    =  FUN(1,MAT)
      ANI  =  FUN(2,MAT)
      at   =  FUN(3,MAT)
      iel  =  INT(FUN(4,MAT))
C
      aem  =  FUN(5,MAT)
      alamb=  FUN(6,MAT)
      aka  =  FUN(7,MAT)
      p0   = -FUN(8,MAT)
      e0   =  FUN(9,MAT)
C=======================================================================
c
c     Elasticity modulus using porosity
      if(iel.eq.1.and.kor.ne.1) then
        evtt=deft(1)+deft(2)+deft(3)
        ett=(1.0d0+e0)*dexp(evtt)-1.0d0
        taumt=(taut(1)+taut(2)+taut(2))/3.
        call eTan(E,ani,ett,aka,taumt)
      endif
c 
      emin=E*1.d-04   
      if(e.lt.emin) then 
        stop 'E < 1e-4*E0'
      endif
c
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
c     Stress correction
      taut(1)=taut(1)-at      
      taut(2)=taut(2)-at      
      taut(3)=taut(3)-at      
      taut(4)=taut(4)
      taut(5)=taut(5)     
      taut(6)=taut(6) 
c      
c     Sigma_E
      call zbir2b(tau,taut,dsig,6)       ! tau(i)=taut(i)+dsig(i)
c
      call clear(ddefp,6)                ! {deP}=0
c
      demp=0.d0
c
c     p0 
      IF(KOR.EQ.1)then
        p0t=p0
      ENDIF
      p0tdt=p0t
c      
c     I1
      call ainvI1(ai1,tau)
      taum=ai1/3.
c
c     J2D
      call ainvJ2d(aj2d,tau)
      aj2d=dsqrt(aj2d*aj2d)
c
c     q 
      aq=dsqrt(3.d0*aj2d)
c
      if(dabs(ai1).lt.toll.and.aq.lt.toll) goto 400
c
c     ks
      aks=alamb-aka
c
c     eV
      evt=deft(1)+deft(2)+deft(3)
      evtdt=def(1)+def(2)+def(3)
c
      et=(1.0d0+e0)*dexp(evt)-1.0d0
      etdt=(1.0d0+e0)*dexp(evtdt)-1.0d0
c
c     bv
      abv=aks/(3.d0*(1.d0+etdt))
c
c     Cam-Clay failure curve
!      call funCC(Fcce,ai1,aq,aem,p0t)
!      Fcc=Fcce
c
      p0tc=aq**2/(taum*aem**2)+taum 
      OCR=p0t/p0tc
c      
c      if(kor.eq.3.and.iter.eq.53) then
c        write(*,*) 'kor = ', kor
c      endif
c
      if(ai1.gt.0.d0) then
        do I=1,6
           tau(I)=0.d0
        enddo
      goto 400
      endif 
c
c     cccccc   OVAJ USLOV PREISPITATI!!! cccccccccccccccccc
      ip0t=0                                           !ccc
      if(ip0t.eq.1) then                               !ccc
        p0t2=p0t/2.                                    !ccc
            if(OCR.lt.0.6d0.and.taum.lt.p0t2) then     !ccc
                p0t=p0tc*0.6d0                         !ccc
            p0tdt=p0t                                  !ccc
        endif                                          !ccc
      endif                                            !ccc
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Cam-Clay failure curve
      call funCC(Fcce,ai1,aq,aem,p0t) 
      Fcc=Fcce

c
c     Yielding check
c     Fcc>0
      if(Fcce.gt.toll) goto 100
c     Fcc<0
      if(Fcce.le.toll) goto 400 
c
      stop 'Not satisfied any condition in Cam-Clay model!'
c=======================================================================
c         PLASTIC DOMAIN
c=======================================================================
  100     continue
c         Set OCR=1
          OCR=1.
c-----------------------------------------------------------------------
c         Derivative dF/dSgima
c-----------------------------------------------------------------------
c ********* dF/dI1 *****************************************************
          dfdi1=aem**2*(2.d0*ai1-3.d0*p0t)/(9.d0*p0t**2)
c
c ********* dF/dq ******************************************************
          dfdq=2.d0*aq/p0t**2
c
c ********* dF/dp0 *****************************************************
          dfdp0=-(2.d0*ai1**2*aem**2-3.d0*ai1*aem**2*p0t+18.d0*aq**2)/
     &           (9.d0*p0t**3)
c
c ********* {dI1/dSigma}T **********************************************
          di1ds(1)=1.d0
          di1ds(2)=1.d0
          di1ds(3)=1.d0
          di1ds(4)=0.d0
          di1ds(5)=0.d0
          di1ds(6)=0.d0
c
c ********* dq/dJ2D ****************************************************
          dqdj2d=0.d0
          if(dabs(aj2d).gt.toll) then
              dqdj2d=(dsqrt(3.d0))/(2.d0*dsqrt(aj2d))
          endif
c     
c ********* {dJ2D/dSigma}T *********************************************
          dj2dds(1)=(2.d0*tau(1)      -tau(2)      -tau(3))/3.
          dj2dds(2)=(    -tau(1) +2.d0*tau(2)      -tau(3))/3.
          dj2dds(3)=(    -tau(1)      -tau(2) +2.d0*tau(3))/3.
          dj2dds(4)= 2.d0*tau(4)
          dj2dds(5)= 2.d0*tau(5)
          dj2dds(6)= 2.d0*tau(6)
c
c ********* dq/dSigma **************************************************
          call jednak(dqds,dj2dds,dqdj2d,6)
c
c ----------------------------------------------------------------------
c         {dF/dSigma}T
          call jednak(dfds1,dqds,dfdq,6)
          call jednak(dfds2,di1ds,dfdi1,6)
c
          call zbir2b(dfds,dfds1,dfds2,6)
c          
c ----------------------------------------------------------------------
c         {dG/dSigma}T={dF/dSigma}T
          call jedna1(dgds,dfds,6)
c ----------------------------------------------------------------------
c         mT
          ajed(1)=1.d0
          ajed(2)=1.d0
          ajed(3)=1.d0
          ajed(4)=0.d0
          ajed(5)=0.d0
          ajed(6)=0.d0
c ----------------------------------------------------------------------
c         lambda
          call clear(alam,6)
          call mnozt1(alam,dfds,elast,6,6)      ! {dFmc/dSigma}T*[Ce]
          fcg1=dot(alam,dgds,6)   ! {dFmc/dSigma}T*[Ce]*{dGmc/dSigma}
c
          fcg2=dot(ajed,dgds,6)
          fcg2=p0t*(1.d0+etdt)/aks*dfdp0*fcg2
!          fcg2=0.d0
c
          fcg=fcg1+fcg2
c
          fce=dot(alam,deps,6)    ! {dFmc/dSigma}T*[Ce]*{de}
c
          dlam=fce/fcg
		  dlamp=dlam
          if(dabs(dlam).lt.toldl) dlam=toldl
c          
C ======================================================================
          call jednak(ddefp,dgds,dlam,6)   ! {deP}=dL{dGmc/dSigma}
          call oduz2b(ddefe,deps,ddefp,6)  ! {deE}={de}-{deP}
          call clear(dsig,6)               ! {dSigma}=0
          call mnozi1(dsig,elast,ddefe,6,6)! {dSigma}=[Ce]{deE}
          call zbir2b(tau,taut,dsig,6)
          call ainvI1(ai1,tau)
          call ainvJ2d(aj2d,tau) 
          aj2d=dsqrt(aj2d*aj2d)
          aq=dsqrt(3.d0*aj2d)
          taum=ai1/3.
          demp=(ddefp(1)+ddefp(2)+ddefp(3))/3.
          dempp=demp
          p0tdt=p0t*dexp(-demp/abv)
          call funCC(Fccp,ai1,aq,aem,p0tdt)
c
C ======================================================================
!          if(Fccp.gt.0.d0)then
!            write(3,*) 'Fccp, dlam = ',Fcc,dlam
!            stop ' Unable to find solution!, Fccp>0' 
!          endif 
C ======================================================================
c         Inicijalizacija za bisekcije
          I=0
          dlama=0.d0
          dlamb=0.2d0*dlam
          dlamc=0.4d0*dlam
          Fmax  = 1.d5*Fcce
          call jednak(ddefpa,dgds,dlama,6)  ! {deP}=dL{dGmc/dSigma}
          call jednak(ddefpb,dgds,dlamb,6)  ! {deP}=dL{dGmc/dSigma}
          call jednak(ddefpc,dgds,dlamc,6)  ! {deP}=dL{dGmc/dSigma}
          call oduz2b(ddefea,deps,ddefpa,6)  ! {deE}={de}-{deP}
          call oduz2b(ddefeb,deps,ddefpb,6)  ! {deE}={de}-{deP}
          call oduz2b(ddefec,deps,ddefpc,6)  ! {deE}={de}-{deP}
          call clear(dsiga,6)               ! {dSigma}=0
          call clear(dsigb,6)               ! {dSigma}=0
          call clear(dsigc,6)               ! {dSigma}=0
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
          aqa=dsqrt(3.d0*aj2da)
          aqb=dsqrt(3.d0*aj2db)
          aqc=dsqrt(3.d0*aj2dc)
          dempa=(ddefpa(1)+ddefpa(2)+ddefpa(3))/3.
          dempb=(ddefpb(1)+ddefpb(2)+ddefpb(3))/3.
          dempc=(ddefpc(1)+ddefpc(2)+ddefpc(3))/3.
          p0tdta=p0t*dexp(-dempa/abv)
          p0tdtb=p0t*dexp(-dempb/abv)
          p0tdtc=p0t*dexp(-dempc/abv)
          call funCC(Fcca,ai1a,aqa,aem,p0tdta)
          call funCC(Fccb,ai1b,aqb,aem,p0tdtb)
          call funCC(Fccc,ai1c,aqc,aem,p0tdtc)
		  dlam=0.d0
c            
C=======================================================================
CD        LOCAL LOOP
!          write(3,*)' I,      dlam,      dlama,      dlamb,      dlamc,
!     &      Fcc,       Fcca,       Fccb,       Fccc,'
  110     I = I + 1
c-----------------------------------------------------------------------
c
c ********* Bretn's method **********
c
            dlam=dlamb+(Fccb/Fcca*(Fcca/Fccc*(Fccb/Fccc-Fcca/Fccc)*
     &          (dlamc-dlamb)-(1.d0-Fccb/Fccc)*(dlamb-dlama)))/
     &          ((Fcca/Fccc-1.d0)*(Fccb/Fccc-1.d0)*(Fccb/Fcca-1.d0))
c
            dlama=dlamb    
            dlamb=dlamc
            dlamc=dlam
c
            call jednak(ddefp,dgds,dlam,6)       ! {deP}=dL{dGmc/dSigma}
            call oduz2b(ddefe,deps,ddefp,6)      ! {deE}={de}-{deP}
            call clear(dsig,6)                   ! {dSigma}=0
            call mnozi1(dsig,elast,ddefe,6,6)    ! {dSigma}=[Ce]{deE}
c
c           {Sigma_t+dt}={Sigma_t}+{dSigma}
            call zbir2b(tau,taut,dsig,6)
c
c           I1
            call ainvI1(ai1,tau)
c
            taum=ai1/3.
c
c           J2D
            call ainvJ2d(aj2d,tau)
            aj2d=dsqrt(aj2d*aj2d)
c
c           q
            aq=dsqrt(3.d0*aj2d)
c 
            if(dabs(ai1).lt.toll.and.aq.lt.toll) goto 400
c
            demp=(ddefp(1)+ddefp(2)+ddefp(3))/3.
c            
            tollam=dlam/dlamp
c
            p0tdt=p0t*dexp(-demp/abv)
c
c           write(3,1002) I,kor,iter,taum,aq
c
c           Cam-Clay failure curve
            call funCC(Fcc,ai1,aq,aem,p0tdt)
c
            Fcca=Fccb
            Fccb=Fccc
            Fccc=Fcc
c
c           write(3,*)'Fcce, Fcc', Fcce, Fcc
c            write(3,1001) I,dlam,dlama,dlamb,dlamc,Fcc,Fcca,Fccb,Fccc
c
            if(I.gt.maxt) then 
               stop 'Max. num. of bisection in Cam-Clay model!'
            endif
c
            if(Fcc.gt.Fmax) then
               write(3,*)'Fcc,Fmax ',Fcc,Fmax
               stop 'Solution diverges in Cam-Clay model!'
            endif
c 
            if(dabs(Fcc).gt.toll.and.dabs(tollam).gt.toll) goto 110 
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      goto 400 
C=======================================================================
c
CE    UPDATES FOR NEXT STEP
  400 CONTINUE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
c=======================================================================
c     Return to original values (Sigma_ij*=Sigma_ij+Sigma0*delta_ij)
      tau(1)=tau(1)+at      
      tau(2)=tau(2)+at      
      tau(3)=tau(3)+at      
      tau(4)=tau(4)
      tau(5)=tau(5)     
      tau(6)=tau(6) 
c
C========================================================================
c     Remaining payload (Racunanje rastojanja od loma)
      ass=aq
      atriq=  dsqrt(3.d0) 
      asm= aem*dsqrt(taum*(p0-taum))
c
      sfactor=(1.d0-ass/asm)*100
      if(sfactor.gt.100.d0) sfactor=100.d0
      if(sfactor.lt.0.d0) sfactor=0.d0
c
      xtdt=sfactor
      xt  =sfactor
c      write(3,*) 'xtdt', xtdt
c========================================================================
c     Corection of values from previous step when convergence is reatched
      call zbir2b(defp1,defpp,ddefp,6)
      call jedna1(def1,def,6)
      call jedna1(tau1,tau,6)
      emp1=emp+demp
      evp1=3.d0*emp1
c
      return
C=======================================================================
 1001 FORMAT(I3,8E12.4)
 1002 FORMAT(3I3,2f12.8)
 1003 FORMAT(6E12.4)
      end
C=======================================================================
C=======================================================================
      SUBROUTINE funCC(F,ai1,aq,aem,p0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     Cam-Clay function
C
      F=(aq/p0)**2-aem**2*ai1/(3.d0*p0)*(1.d0-ai1/(3.d0*p0))
      return
      end
C=======================================================================
C=======================================================================
      SUBROUTINE eTan(EE,ani,et,aka,taum)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     tE 
C
      EE=-3.d0*(1.d0-2.d0*ani)*(1.d0+et)/aka*taum
      return
      end
C=======================================================================