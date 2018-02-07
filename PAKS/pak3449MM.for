C       Micunov materijalni model
C      Oznacen brojem 49 za sada iskopiran Rakicev DP 41
      
C
C=======================================================================
C=======================================================================
CS    GENERALISANI MODEL SA KAPOM  3/D ELEMENT
CE    GENERALIZED CAP MODEL 3D ELEMENT
C=======================================================================
C=======================================================================
CE    SUBROUTINE D3M49
CE               TI3441
C
      SUBROUTINE D3M49(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC)
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
      IF(IDEBUG.GT.0) PRINT *, ' D3M41'
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
      CALL TI3449(PLAST(LTAU),PLAST(LDEFT),PLAST(LDEFPP),
     1            PLAST(LEMP),PLAST(LXT),
     1            PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LEMP1),PLAS1(LXTDT), 
     1            A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC)
C
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI3449(TAUT,DEFT,DEFPP,EMP,XT,
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
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1               DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8 
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /ITERBR/ ITER
      COMMON /CONMAT/ AE,EP,DVT
      COMMON /MATERb/ korz(100,100,3),evg(100,100,3)
      COMMON /CDEBUG/ IDEBUG
      common /ak/ akapa(100000,10),aw(100000,10)
      common /comeplast/  ceplast(100000,8,6)
      common /prolaz/ iprolaz(100000,8)
C
      DIMENSION TAUT(6),DEFT(6),DEFPP(6),TAU(6),DEF(6),TAU1(6),DEF1(6), 
     +          DEFP1(6)
      DIMENSION FUN(11,*),NTA(*),DSIG(6),DEPS(6),DFDS(6),DGDS(6),ALAM(6)
     +         ,CP(6,6),POM(6), CEP(6,6),DDEFE(6),TAUDE(6)
     +         ,DFdDS(6),DFcDS(6),DGdDS(6),DGcDS(6),ALAMd(6),ALAMc(6)
     +         ,al1p(6),al2p(6),ddefpd(6),ddefpc(6),taue(6)
     1         ,stress(6),e_new(6),e_elas_n(6),e_elas_n1(6),kroneker(6)!MM
     1         ,a(17),astress0(6),astress(6),eplas(6),eplas0(6)
     1         ,eplasStari(6),d_eplas(6),Replas(6),a_mu(6),statev(100)
C
      parameter (one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0,zero=0.0d0)!MM
      data newton,toler,temp0,coef,yield0/8 ,1.d-6,273.,0.0d0,25.0d0/!MM
      IF(IDEBUG.EQ.1) PRINT *, 'TI3449'
C
      IF(IRAC.EQ.2) RETURN
CE    BASIC KONSTANTS
      TOL =   1.D-8
      TOLL=  -1.D-5
      MAXT=   500.D0
      ntens = 6 !MM
      ndi = 3 !MM
      tol1 = 1.d-7 !MM
	tol2 = 1.d-7 !MM
      tol2 = 1.0d-6!MM
	tolk = 1.d-5 !MM
      do k2=1, ndi
         kroneker(k2)     = one
         kroneker(k2+ndi) = zero
      end do
C====================================================================
C
CE.   MATERIAL CONSTANTS 
      E       = FUN(1,MAT) 
      ANI     = FUN(2,MAT) 
C
      AK      = FUN(3,MAT)
      ALF     = FUN(6,MAT)
      T       = FUN(7,MAT)
      X0      =-FUN(8,MAT)
C
      W       =-FUN(9,MAT)
      D       =-FUN(10,MAT)

      
C====================================================================
C
      CM      = E/(1-2.D0*ANI)
      G       = 0.5D0*E/(1.D0+ ANI)
      DUM     = 0.01D0
      XTMAX   =-DLOG(DUM)/D+X0
C
      IF(KOR.EQ.1) THEN
       XT=X0
      ENDIF
      
      
              
      
!C -----------------------------------------------------------
!c
!c     paneerselvam phd   (page 165 tens constants)
!c
      
      	  
      a3 = -1.81e5
      a4 = -0.906e5

      a6 = 0.
      a7 = 0.
      a8 = 0.
      a9 = 0.
      x = 1.03e-3
      a_l = 2          
      alfamm = 14.2  !17.1   
      h = 50.
      beta = 8.23e4     !1.34e4
      a_m  = one
      gama = -5.e-3 
      alpha_t  = 4.e-5
	  
	!stari
	!a1 = 70.3e2
	!a2 = 0.036e5 !7.5e4
	!a5 = -0.6046e5
	!novi
	a1 =  0.036e5	  
      a2 = -0.6046e5
	a5 = 70.3e2
!c -----------------------------------------------------------
      a(1)=a1
      a(2)=a2
      a(3)=a3
      a(4)=a4
      a(5)=a5
      a(6)=a6
      a(7)=a7
      a(8)=a8
      a(9)=a9
      a(10)=x
      a(11)=a_l
      a(12)=alfamm
      a(13)=h
      a(14)=beta
      a(15)=m
      a(16)=gama
      a(17)=alpha_t
       
!  KONSTANTE
!     NDI: Broj direkthih komponenti napona u datom trenutku!     NDI=3
!     NSHR: Broj smicajnih komponenti napona u datom trenutku!      NSHR=3
!     NTENS: velicina niza napona ili deformacija (NDI + NSHR)!     NTENS=6
!     NOEL: Broj elementa
      noel=1 !zameniti sa realnom brojem elementa iz paka
!     NPT: Broj integracione ta?ke
      npt=1
!  UILAZNE VELICINE U UMAT
!     STRAN(NTENS): Niz koji sadrzi ukupne deformacije na pocetku inkrenenta
!     DSTRAN(NTENS): Niz inkremenata deformacija
!  IZLAZNE VELICINE    
!     STRESS(NTENS): TENZOR NAPONA OVO JE I ULAZ I IZLAZ
!     NA POCETKU DOBIJEMO TENZOR NAPONA NA POCETKU INKREMENTA
!     PA U UMATU TREBA DA IZRACUNAMO TENZOR NAPONA NA KRAJU INKREMENTA    
!     DDSDDE(NTENS,NTENS)   PARCIJALNI IZVOD INKREMENTA NAPONA PO INKREMENTU DEFORMAICJE
!     OVO JE CISTO IZLAZNA VELICINA KOJA SE RACUNA U UMATU
!     STRESS JE MNOGO VAZNIJI OD DDSDDE
!     DDSDDE SE KORISTI U ABAQUSU ZA PROVERU KONVERGENCIJE       
!  OSTALE ULAZNE ILI IZLAZNE VELICINE NE KORISTIMO
    
!  UNUTRASNJE PROMENLJIVE
!     STATEV NIZ UNUTRASNJIH PROMENLJIVIH
!     NSTATV BROJ UNUTRASNJIH PROMENLJIVIH
!     OVO VISE NE KORISTIMO SADA SVE CUVAMO U COMMON STRUKTURAMA          

!c********************************************************************
!c
!c               predictor corector algorithm
!c
!c********************************************************************     

!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!c***************      1) predictor phase      ***********************  
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      
      

      
      
      
      
C====================================================================
C
CE    ELASTICITY MATRIX FOR 3D
      CALL MEL3EL(ELAST,E,ANI)
      IF(IRAC.EQ.2) RETURN
C     
CE    MEAN STRAIN emE braon knjiga 2.2.21
      EMT=(DEF(1)+DEF(2)+DEF(3))/3.0D0
C
CE    TOTAL DEVIATORIC STRAINS e'=e-em braon knjiga 2.2.23
      DO I=1,6
         IF(I.LE.3) THEN
            DEFDPR(I)=DEF(I)-EMT
         ELSE
            DEFDPR(I)=0.5D0*DEF(I) !smicuce komponente e = gama/2
         ENDIF 
      ENDDO

      !  UILAZNE VELICINE U UMAT
!     STRAN(NTENS): Niz koji sadrzi ukupne deformacije na pocetku inkrenenta
!     DSTRAN(NTENS): Niz inkremenata deformacija
  ! u paku stran je DEF    

      
CE    TRIAL ELASTIC DEVIATORIC STRAINS e''=e'-emp' braon knjiga 2.2.23
      DO I=1,6
         IF(I.LE.3) THEN
            DEFDS(I)=DEFDPR(I)-DEFPP(I)+EMP
         ELSE
            DEFDS(I)=DEFDPR(I)-DEFPP(I)
         ENDIF
      ENDDO
C
CE    TRIAL ELASTIC MEAN STRAIN em''=em-emp
      EMS=EMT-EMP
! ucitava iz prethodnog koraka
	  
       
      do k1=1,ntens
	stress(k1)=0
	!eplas0(k1) = ceplast(noel,npt,k1) 
		
	!e_new(k1)     = stran(k1)+dstran(k1)-eplas0(k1)
	!e_elas_n(k1)  = stran(k1)+dstran(k1)-eplas0(k1)
	!e_elas_n1(k1) = stran(k1)+dstran(k1)-eplas0(k1)
          
      e_new(k1)     = DEFDS(k1)+kroneker(k1)*EMS!MM
	e_elas_n(k1)  = DEFDS(k1)+kroneker(k1)*EMS!MM
	e_elas_n1(k1) = DEFDS(k1)+kroneker(k1)*EMS!MM
      enddo
      a_kapa0  = akapa(noel,npt)  ! 6.20  ! ucitava iz prethodnog koraka
	  if (a_kapa0.lt.tolk) a_kapa0=tolk
! ucitava iz prethodnog koraka
C
      call stressMM(a,ntens,e_elas_n,astress0)

CE    ELASTIC SOLUTION SE=2G*e'' braon knjiga 2.2.24
!     e'  DEFDS devijatorska deformacija     
      DO I=1,6
          TAUD(I)=2.0D0*G*DEFDS(I)
          TAUD(I)=astress0(I)
      ENDDO
CD
      DO I=1,6
          TAUDE(I)=TAUD(I)
      ENDDO
CE    TRIAL ELASTIC MEAN STRESS sigmaE=cm*em'' braon knjiga 2.2.17
      SMTE=EMS*CM
C=================================================================
C     znaci sto je u abaqusu s_dev to je u paku TAUDE
CE    SECOND INVARIANT OF DEVIATORIC STRESSES J2d braon knjiga 1.2.44
!     ! a_j2 = ( s_dev(1)**2+s_dev(2)**2+s_dev(3)**2)/two+
!     !1          s_dev(4)**2+s_dev(5)**2+ s_dev(6)**2
      
      AJ2DE=0.5*(TAUDE(1)*TAUDE(1)+TAUDE(2)*TAUDE(2)+TAUDE(3)*TAUDE(3))+
     1           TAUDE(4)*TAUDE(4)+TAUDE(5)*TAUDE(5)+TAUDE(6)*TAUDE(6)
c     CALL WRR6 (TAUDE,6,'TAUDE')
      
      
      a_j2 = AJ2DE
      
      call loadingf(f10,astress0,a,ntens,ndi,s_dev,a_j2,a_mu) ! 6.21
		f0   = abs(f10) - h*a_kapa0 
    
      call stressMM(a,ntens,e_elas_n1,astress)    
                   
          !	call loadingf(f1,astress,a,ntens,ndi,s_dev,a_j2,a_mu) ! 6.21
	!	f   = abs(f1) - h*a_kapa0 
                                      
	    call noviddsdde(a,ddsdde,ntens,e_elas_n)
	    ! call    hyperconstitutive (a,ddsdde,ntens,e_elas_n1)
          
          do k1=1,ntens        
		    eplas(k1)  = eplas0(k1)
			stress(k1) = astress0(k1)
          enddo
	    dkapa = 0
      
          
          !c------------------ compute  loading surface f-----------------------
          call loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu) ! 6.21
        
          f = abs(f1) - h*a_kapa0 
		a_kapa = a_kapa0
	if (f.le.zero) then 
		goto 500   
      endif
      
      dkapau = 0
	call loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu)
      f  = abs(f1) - h*a_kapa				
	a_kxl = x+a_kapa0**a_l	
      
!c------------------  end of elastic predictor ----------------------
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!c***************      1) corector phase      ***********************  
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
          
      
      
      do kewton = 1,newton
            dkapa0  =  (((f/beta)**1)/a_kxl)
			a_kapa = a_kapa0 + dkapa0
			a_kxl = x+a_kapa**a_l
			dkapa  =  (((f/beta)**1)/a_kxl)
			
			do k1 = 1,ntens
			    eplasStari(k1)   = eplas(k1)
			    d_eplas(k1)   =  (dkapa)*a_mu(k1)*dtime
			    eplas(k1)   = eplas0(k1) + d_eplas(k1)
			    Replas(k1) = eplas(k1) - eplasStari(k1)
			enddo
				deplas_int = (Replas(1)**2+Replas(2)**2+Replas(3)**2+
     1    two*Replas(4)**2+two*Replas(5)**2+two*Replas(6)**2)**0.5
	 
	        f  = abs(f1) - h*a_kapa
	        a_kxl = x+a_kapa**a_l
			skonvergencija = abs(a_kapa - a_kapa0)
			if ((skonvergencija.lt.tol1).and.(deplas_int.lt.tol2)) then
			!write(6,*) 'radi'
			goto 33
			endif
			
      enddo
	
		!dkapa = a_kapa - a_kapa0
 33        continue

          do k1 = 1,ntens
		    !eplas(k1)   =  dkapa*a_mu(k1) !6.14
		    eplas(k1)   = eplas(k1) + (a_kapa-a_kapa0)*a_mu(k1)*dtime
		    !eplas(k1)   = eplas(k1) + dkapa*a_mu(k1)*dtime
		    !eplas(k1)   = eplas(k1) + dkapa*a_mu(k1)*dtime
		    !e_elas_n1(k1) = stran(k1)+dstran(k1)-eplas(k1)
		    ceplast(noel,npt,k1) = eplas(k1)
              akapa(noel,npt) = a_kapa  ! 6.20  ! ucitava iz prethodnog koraka
		enddo
		
		call stressMM(a,ntens,e_elas_n1,astress)

	    do k1=1,ntens        
			stress(k1) =astress(k1)
              TAU(k1) =astress(k1)
          enddo
	
	    call noviddsdde(a,ddsdde,ntens,e_elas_n1)
	   	!invarijante	
	

      
      AJ2DEQ  =DSQRT(AJ2DE)
C
CE    FIRST STRESS INVARIANT I1E=3*sigmamE braon knjiga 2.2.15 1.2.21
      AI1E    =3.0D0*SMTE
C
C ================================================================     

C     DRUCKER-PRAGER LINE
      FDP=ALF*AI1E+AJ2DEQ-AK
C      
C     CAP LINE
      FC=XT-AI1E
C
      SMTDT=SMTE
      DELEMP=0.d0
CE    ELASTIC SOLUTION
!      IF(FDP.LE.TOL) GO TO 400
      goto 400
C==================================================================
C
CE    PLASTIC STRAIN
C
CE    DRUCKER-PRAGER YIELDING FDP>0 & FC<0
      IF(FDP.GT.TOL)GOTO 100    
C      
      STOP ' NONDEFINED AREA - DRUCKER-PRAGER '
C===================================================================
C===================================================================
C
CE    DRUCKER-PRAGER YIELDING FDP>0 & FC<0
  100 CONTINUE
      IYIELD=1
c     WRITE(3,*)'OBLAST 1'
      IF(DABS(ALF).LT.TOL) THEN
CE        INCREMENT OF MEAN PLASTIC STRAIN Demp
        AJ2DQ=AK
C         sqrt(J2D)=sqrt(J2DE)-Dlam*G
CE        Dlambda
        DLAM=(AJ2DEQ-AJ2DQ)/G
          DELEMP=0.D0
          ELSE
CE      INCREMENT OF MEAN PLASTIC STRAIN Demp
        DELEMP=(ALF*AI1E+AJ2DEQ-AK)/(3.0D0*ALF*CM+G/ALF)
CE      Dlambda.
        DLAM=DELEMP/ALF
C       sqrt(J2D)=sqrt(J2DE)-Dlam*G (6.2.29)
        AJ2DQ=AJ2DEQ-DLAM*G
      ENDIF
C
CE    MEAN STRESS sigmam=sigmamE-cm*Demp
      SMTDT=SMTE-CM*DELEMP
C
CE    FIRST STRESS INVARIANT
      AI1=3.0D0*SMTDT
C
CE    MEAN PLASTIC STRAIN emp=emp+Demp
      EMP1=EMP+DELEMP
C
C====================================================================
CE    DEVIATORIC STRESS S
      DEN=1.0D0+DLAM*G/AJ2DQ
      DO I=1,6
         TAUD(I)=TAUDE(I)/DEN
      ENDDO
C
CE    ADD INCREMENT OF DEVIATORIC PLASTIC STRAIN ep'=ep'+De'p 
      COEF=DLAM/AJ2DQ/2.0D0
      DO I=1,6
         DEFP1(I)=DEFPP(I)+COEF*TAUD(I)
      ENDDO
C
      AJ2D=AJ2DQ*AJ2DQ
      GO TO 400
C====================================================================
C====================================================================
C
CE    UPDATES FOR NEXT STEP
  400 CONTINUE
CE    ADD INCREMENT OF MEAN PLASTIC STRAIN ep=ep+Demp
      DO I=1,6
         IF(I.LE.3) THEN
            DEFP1(I)=DEFP1(I)+DELEMP
         ELSE
            DEFP1(I)=DEFP1(I)
         ENDIF
      ENDDO
C
CE    STRESS CALCULATION sigma=S+sigmam
      DO I=1,6
         IF(I.LE.3) THEN
            TAU(I)=TAUD(I)+SMTDT
         ELSE
            TAU(I)=TAUD(I)
         ENDIF
      ENDDO
  
C========================================================================
  500 CONTINUE
      DO I=1,6
          DEF1(I)=DEF(I) 
          TAU1(I)=TAU(I)
      ENDDO
!     cuva za naredni korak
      dw = 0
	do k1=1,ntens
	    dw=dw+stress(k1)*d_eplas(k1)
		eplas(k1)   = eplas(k1) + dkapa*a_mu(k1)*dtime
          ceplast(noel,npt,k1)  =  eplas(k1)
      enddo 
		akapa(noel,npt) = a_kapa
		w = aw(noel,npt)
		w=w+dw
		aw(noel,npt)=w
      
      RETURN
      END
      
      
      
      
      
      subroutine loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu)
    
!     ulazne promenljive   
!         real a(17), alfa,gama       
!         integer ntens,ndi
!         real  stress(ntens)
!     ulazne promenljive  

!     izlazne promenljive
!         real  f1
!         real a_j2
!         real a_mu(ntens)
!         real s_dev(ntens)
!     izlazne promenljive     

      dimension stress(ntens),s_dev(ntens),a(17),a_mu(ntens),
     1   kroneker(ntens),st(ntens)
      
      parameter (one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0) 
      alfa =  a(12)
      gama = a(16)
      s_i1=0
	s_i2=0
	s_i3=0
      do k2=1, ndi
         kroneker(k2)     = one
         kroneker(k2+ndi) = zero
	   st(k2) = stress(k2) ! unutrasnja inverzija znaka  za racunanje funkcije tecenja
!         OVO JE JAKO BITNO DA LI JE TACNO RADI ZA SAMO PRVA 3 NA ZA 4,5,6????
      end do 
            
!c     s_napon invarijante  i funkcija f 

!     s_i1 = stress(1)+stress(2)+stress(3)
	s_i1 = st(1)+st(2)+st(3)
	  
    !  s_i2 = (stress(1)**2+stress(2)**2+stress(3)**2)/two+
    ! 1       stress(4)**2+stress(5)**2+stress(6)**2
    ! 1       -(s_i1**2)/two
	 s_i2 = st(1)*st(2)+st(2)*st(3)+st(1)*st(3) !plava knjiga 1.4.6
     1 		-st(4)**2-st(5)**2-st(6)**2
!         Zasto nije isto kao u braon knjizi 1.2.21????
		
!      s_i3 = stress(1)*stress(2)*stress(3)+
!     1 two*stress(4)*stress(5)*stress(6)-stress(1)*stress(6)**2- 
!     1 stress(2)*stress(5)**2-stress(3)*stress(4)**2
       s_i3 = st(1)*st(2)*st(3)+
     1 two*st(4)*st(5)*st(6)-st(1)*st(6)**2- 
     1 st(2)*st(5)**2-st(3)*st(4)**2
!     plava knjiga 1.4.9. braon knjiga 1.2.21
      do  k1=1,ndi
          !s_dev(k1)     = stress(k1) - s_i1/three
          !s_dev(k1+ndi) = stress(k1+ndi)
		s_dev(k1)     = st(k1) - s_i1/three
          s_dev(k1+ndi) = st(k1+ndi)
      enddo
  
      a_j2 = ( s_dev(1)**2+s_dev(2)**2+s_dev(3)**2)/two+
     1          s_dev(4)**2+s_dev(5)**2+ s_dev(6)**2 
     
      f1 = s_i1*s_i2 + alfa*s_i3 
	   
      do k1=1,ntens
	    if (a_j2.gt.0) then  !5.21
              a_mu(k1) = gama*kroneker(k1)+(0.5*s_dev(k1)/(a_j2**0.5))
		endif
      enddo   
      return
      end
      
!c---------------------------------subroutine loadingf           
           
      subroutine noviddsdde(a,ddsdde,ntens,e_elas)

!     ulazne promenljive   
!         real a(17), a1,a2,a3,a4,a5,a6,a7,a8,a9      
!         integer ntens
!         real  e_elas(ntens)
!     ulazne promenljive  

!     izlazne promenljive
!         real ddsdde(ntens,ntens)
!     izlazne promenljive

      dimension ddsdde(ntens,ntens), a(17), e_elas(ntens)

      a1=a(1)
      a2=a(2)
      a3=a(3)
      a4=a(4)
      a5=a(5)
      a6=a(6)
      a7=a(7)
      a8=a(8)
      a9=a(9)
      do k1=1,ntens
	  do k2=2,ntens
	    ddsdde(k1,k2)=0
	  enddo
      enddo

       ddsdde(1,1)= 2*a1 + 4*a5 + 4*e_elas(1)*a2 + 
	1   4*e_elas(1)*a4 + 6*a3*(2*e_elas(1) + 
	1 2*e_elas(2) + 2*e_elas(3)) + 2*a4*(e_elas(1) + 
	1 e_elas(2) + e_elas(3))
 
       ddsdde(1,2)= 4*a5 + 2*e_elas(1)*a4 + 2*e_elas(2)*a4 + 
	1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
      
       ddsdde(1,3)= 4*a5 + 2*e_elas(1)*a4 + 2*e_elas(3)*a4 + 
	1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(1,4)= e_elas(4)*a2 + e_elas(4)*a4 + e_elas(4)*a2 + 
	1   e_elas(4)*a4
       ddsdde(1,5)= e_elas(5)*a2 + e_elas(5)*a4 + e_elas(5)*a2 + 
	1   e_elas(5)*a4
       ddsdde(1,6)= E23*a4 + e_elas(6)*a4

       ddsdde(2,1)= 4*a5 + 2*e_elas(1)*a4 + 2*e_elas(2)*a4 + 
	1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(2,2)= 2*a1 + 4*a5 + 4*e_elas(2)*a2 + 4*e_elas(2)*a4 + 
	1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3)) + 
	1   2*a4*(e_elas(1) + e_elas(2) + e_elas(3))
     
       ddsdde(2,3)= 4*a5 + 2*e_elas(2)*a4 + 2*e_elas(3)*a4 + 
	1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(2,4)= e_elas(4)*a2 + e_elas(4)*a4 + e_elas(4)*a2 + 
	1   e_elas(4)*a4
       ddsdde(2,5)= e_elas(5)*a4 + e_elas(5)*a4
       ddsdde(2,6)= E23*a2 + E23*a4 + e_elas(6)*a2 + e_elas(6)*a4
             
       ddsdde(3,1)= 4*a5 + 2*e_elas(1)*a4 + 2*e_elas(3)*a4 + 
	1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(3,2)= 4*a5 + 2*e_elas(2)*a4 + 2*e_elas(3)*a4 + 
	1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(3,3)= 2*a1 + 4*a5 + 4*e_elas(3)*a2 + 4*e_elas(3)*a4 + 
	1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3)) + 
	1   2*a4*(e_elas(1) + e_elas(2) + e_elas(3))
     
       ddsdde(3,4)= e_elas(4)*a4 + e_elas(4)*a4
       ddsdde(3,5)= e_elas(5)*a2 + e_elas(5)*a4 + e_elas(5)*a2 + 
	1   e_elas(5)*a4
       ddsdde(3,6)= E23*a2 + E23*a4 + e_elas(6)*a2 + e_elas(6)*a4
       
       ddsdde(4,1)= 2*e_elas(4)*a2 + 2*e_elas(4)*a4
       ddsdde(4,2)= 2*e_elas(4)*a2 + 2*e_elas(4)*a4
       ddsdde(4,3)= 2*e_elas(4)*a4
       ddsdde(4,4)= a1 + a2*(e_elas(1) + e_elas(2)) + a4*(e_elas(1) + 
	1   e_elas(2) + e_elas(3))
       ddsdde(4,5)= e_elas(6)*a2
       ddsdde(4,6)= e_elas(5)*a2
       
       ddsdde(5,1)= 2*e_elas(5)*a2 + 2*e_elas(5)*a4
       ddsdde(5,2)= 2*e_elas(5)*a4
       ddsdde(5,3)= 2*e_elas(5)*a2 + 2*e_elas(5)*a4
       ddsdde(5,4)= e_elas(6)*a2
       ddsdde(5,5)= a1 + a2*(e_elas(1) + e_elas(3)) + a4*(e_elas(1) + 
	1   e_elas(2) + e_elas(3))
       ddsdde(5,6)= e_elas(4)*a2
       
       ddsdde(6,1)= 2*e_elas(6)*a4
       ddsdde(6,2)= 2*e_elas(6)*a2 + 2*e_elas(6)*a4
       ddsdde(6,3)= 2*e_elas(6)*a2 + 2*e_elas(6)*a4
       ddsdde(6,4)= e_elas(5)*a2
       ddsdde(6,5)= e_elas(4)*a2
       ddsdde(6,6)= a1 + a2*(e_elas(2) + e_elas(3)) + a4*(e_elas(1) +  
	1   e_elas(2) + e_elas(3))
       return
      end
           
!c----------------------subroutine noviddsdde 
      
      subroutine hyperconstitutive(a,ddsdde,ntens,e_elas)
      
!     ulazne promenljive   
!         real a(17), a1,a2,a3,a4,a5,a6,a7,a8,a9      
!         integer ntens
!         real  e_elas(ntens)
!     ulazne promenljive  

!     izlazne promenljive
!         real ddsdde(ntens,ntens)
!     izlazne promenljive

      dimension ddsdde(ntens,ntens), a(17), e_elas(ntens)

      a1=a(1)
      a2=a(2)
      a3=a(3)
      a4=a(4)
      a5=a(5)
      a6=a(6)
      a7=a(7)
      a8=a(8)
      a9=a(9)
      do k1=1,ntens
	  do k2=2,ntens
	    ddsdde(k1,k2)=0
        enddo
      enddo

       ddsdde(1,1)= 2*a1 + a2 + 2*e_elas(1)*a4 + 2*e_elas(1)*a5 + 
     1 3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3)) + a4*(e_elas(1) 
     1 + e_elas(2) + e_elas(3))
 
       ddsdde(1,2)= 2*a1 + e_elas(1)*a4 + e_elas(2)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
      
       ddsdde(1,3)= 2*a1 + e_elas(1)*a4 + e_elas(3)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(1,4)= 2*e_elas(4)*a4 + 2*e_elas(4)*a5
       ddsdde(1,5)= 2*e_elas(5)*a4 + 2*e_elas(5)*a5
       ddsdde(1,6)= 2*e_elas(6)*a4

       ddsdde(2,1)= 2*a1 + e_elas(1)*a4 + e_elas(2)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(2,2)= 2*a1 + a2 + 2*e_elas(2)*a4 + 2*e_elas(2)*a5 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3)) + a4*(e_elas(1) 
     1  + e_elas(2) + e_elas(3))
     
       ddsdde(2,3)= 2*a1 + e_elas(2)*a4 + e_elas(3)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(2,4)= 2*e_elas(4)*a4 + 2*e_elas(4)*a5
       ddsdde(2,5)= 2*e_elas(5)*a4
       ddsdde(2,6)= 2*e_elas(6)*a4 + 2*e_elas(6)*a5
             
       ddsdde(3,1)= 2*a1 + e_elas(1)*a4 + e_elas(3)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(3,2)= 2*a1 + e_elas(2)*a4 + e_elas(3)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(3,3)= 2*a1 + a2 + 2*e_elas(3)*a4 + 2*e_elas(3)*a5 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3)) + a4*(e_elas(1) 
     1 + e_elas(2) + e_elas(3))
     
       ddsdde(3,4)= 2*e_elas(4)*a4
       ddsdde(3,5)= 2*e_elas(5)*a4 + 2*e_elas(5)*a5
       ddsdde(3,6)= 2*e_elas(6)*a4 + 2*e_elas(6)*a5
       
       ddsdde(4,1)= 2*e_elas(4)*a4 + 2*e_elas(4)*a5
       ddsdde(4,2)= 2*e_elas(4)*a4 + 2*e_elas(4)*a5
       ddsdde(4,3)= 2*e_elas(4)*a4
       ddsdde(4,4)= 2*a2 + a5*(2*e_elas(1) + 2*e_elas(2)) + 
     1  2*a4*(e_elas(1) + e_elas(2) + e_elas(3))
       ddsdde(4,5)= 2*e_elas(6)*a5
       ddsdde(4,6)= 2*e_elas(5)*a5
       
       ddsdde(5,1)= 2*e_elas(5)*a4 + 2*e_elas(5)*a5
       ddsdde(5,2)= 2*e_elas(5)*a4
       ddsdde(5,3)= 2*e_elas(5)*a4 + 2*e_elas(5)*a5
       ddsdde(5,4)= 2*e_elas(6)*a5
       ddsdde(5,5)= 2*a2 + a5*(2*e_elas(1) + 2*e_elas(3)) + 
     1  2*a4*(e_elas(1) + e_elas(2) + e_elas(3))
       ddsdde(5,6)= 2*e_elas(4)*a5
       
       ddsdde(6,1)= 2*e_elas(6)*a4
       ddsdde(6,2)= 2*e_elas(6)*a4 + 2*e_elas(6)*a5
       ddsdde(6,3)= 2*e_elas(6)*a4 + 2*e_elas(6)*a5
       ddsdde(6,4)= 2*e_elas(5)*a5
       ddsdde(6,5)= 2*e_elas(4)*a5
       ddsdde(6,6)= 2*a2 + a5*(2*e_elas(2) + 2*e_elas(3)) + 
     1  2*a4*(e_elas(1) + e_elas(2) + e_elas(3))    
       return
      end         

!c----------------------subroutine hyperconstitutive     
      
      subroutine stressMM(a,ntens,e_elas,stress)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!     ulazne promenljive   
!         real a(17), a1,a2,a3,a4,a5,a6,a7,a8,a9      
!         integer ntens
!         real  e_elas(ntens)
!     ulazne promenljive  

!     izlazne promenljive
!         real stress(ntens,ntens)
!     izlazne promenljive
      

      dimension stress(ntens), a(17), e_elas(ntens)

      a1=a(1)
      a2=a(2)
      a3=a(3)
      a4=a(4)
      a5=a(5)
      a6=a(6)
      a7=a(7)
      a8=a(8)
      a9=a(9)
      do k1=1,ntens
	    stress(k1)=0 
      enddo
      
!invarijante
	e_i1n = e_elas(1)+e_elas(2)+e_elas(3)
		
      e_i2n = e_elas(1)*e_elas(2)+e_elas(2)*e_elas(3)+
	1 e_elas(3)*e_elas(1)-e_elas(4)*e_elas(4)-
	1 e_elas(5)*e_elas(5)-e_elas(6)*e_elas(6)
	
	e_i3n = e_elas(1)*e_elas(2)*e_elas(3)-
	1 e_elas(1)*e_elas(6)*e_elas(6)-
	1 e_elas(2)*e_elas(5)*e_elas(5)-
	1 e_elas(3)*e_elas(4)*e_elas(4)+
	1 2*e_elas(4)*e_elas(5)*e_elas(6)
!invarijante
          
                    
   	stress(1)= (2*a5*e_i1n+3*a3*e_i1n*e_i1n+a4*e_i2n)+
	1 (a1+a4*e_i1n)*e_elas(1)+a2*(e_elas(1)*e_elas(1)
	1 +e_elas(4)*e_elas(4)+e_elas(5)*e_elas(5))
	
	stress(2)= (2*a5*e_i1n+3*a3*e_i1n*e_i1n+a4*e_i2n)+
	1 (a1+a4*e_i1n)*e_elas(2)+a2*(e_elas(4)*e_elas(4)
	1 +e_elas(2)*e_elas(2)+e_elas(6)*e_elas(6))
	
	 stress(3)= (2*a5*e_i1n+3*a3*e_i1n*e_i1n+a4*e_i2n)+
	1 (a1+a4*e_i1n)*e_elas(3)+a2*(e_elas(5)*e_elas(5)
	1 +e_elas(6)*e_elas(6)+e_elas(3)*e_elas(3))
	
	
      stress(4)= (a1+a4*e_i1n)*e_elas(4)
	1 +a2*(e_elas(1)*e_elas(4)
	1 +e_elas(4)*e_elas(2)+e_elas(5)*e_elas(6))
	
      stress(5)= (a1+a4*e_i1n)*e_elas(5)
	1 +a2*(e_elas(1)*e_elas(5)
	1 +e_elas(4)*e_elas(6)+e_elas(5)*e_elas(3))
	
      stress(6)= (a1+a4*e_i1n)*e_elas(6)
	1 +a2*(e_elas(4)*e_elas(5)
	1 +e_elas(2)*e_elas(6)+e_elas(6)*e_elas(3))   
      
      return
      end

!c------------------------subroutine stressMM