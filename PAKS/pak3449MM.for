C       Micunov materijalni model
C      Oznacen brojem 49 nastao kopiranjem i modifikovanjem DP 41
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
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAU(6),DEF(6)
C
      IF(IDEBUG.GT.0) PRINT *, ' D3M49'
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      MATE=MREPER(4)
C
      LTAU=LPOCG
      LDEFT=LTAU + 6
      LDEFPP=LDEFT + 6
      LEMP=LDEFPP + 6
      LKAPA=LEMP + 1
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 6
      LDEFP1=LDEFT1 + 6
      LEMP1=LDEFP1 + 6
      LXTDT=LEMP1 + 1
C
      CALL TI3449(LDEFPP,PLASTMM(LDEFPP),PLASTMM(LKAPA),
     1            PLAS1(LTAU1),TAU,DEF,IRAC)
C
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI3449(LDEFPP,DEFPP,a_kapaP,TAU1,TAU,DEF,IRAC)

      IMPLICIT NONE
      
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERBR/ ITER
      COMMON /CDEBUG/ IDEBUG
      
      !pak
      INTEGER LNKDT,LDTDT,LVDT,NDT,KOR
      DOUBLE PRECISION DT,VREME
      DOUBLE PRECISION ELAST,XJ,ALFA,TEMP0,DET
      INTEGER NLM,KK
      INTEGER ITER,IDEBUG,IRAC,LDEFPP
      DOUBLE PRECISION DEFPP,DEF,DEFE,TAU,TAU1
      !pak
      
      !mm
      INTEGER kewton,newton,NDI,NTENS,K1,K2,I
      DOUBLE PRECISION s_dev,stress
      DOUBLE PRECISION one,two,three,six,zero,kroneker
      DOUBLE PRECISION TOL2,TOLK
      DOUBLE PRECISION a,a1,a2,a3,a4,a5,a6,a7,a8,a9,beta,gama,h
      DOUBLE PRECISION x,a_l,alfamm,a_m,a_kapa0,a_kapaP
      DOUBLE PRECISION a_kxl,a_kapa,a_j2,eplasstari,skonvergencija
      DOUBLE PRECISION deplas_int,Replas,d_eplas,dkapa0,dkapa,a_mu
      DOUBLE PRECISION eplas,eplas0,f,f1,ALFMMM,dtime
      !mm
C
      DIMENSION DEF(6),DEFE(6),DEFPP(6),TAU(6),TAU1(6),a(17) 
     1 ,stress(6),s_dev(6),kroneker(6),eplas(6),eplas0(6)!MM
     1 ,eplasStari(6),d_eplas(6),Replas(6),a_mu(6)
C
      
      one=1.0d0
      two=2.0d0
      three=3.0d0
      six=6.0d0
      zero=0.0d0!MM
      
      IF(IDEBUG.EQ.1) PRINT *, 'TI3449'
      IF(IRAC.EQ.2) RETURN
      
      kewton = 1
      newton = 50000
 
CE    BASIC KONSTANTS
      ntens = 6 !MM velicina niza napona ili deformacija
      ndi = 3 !MM Broj direkthih komponenti napona u datom trenutku!
      tol2 = 1.0d-6!MM
	tolk = 1.d-4 !MM
      do k2=1, ndi
         kroneker(k2)     = one
         kroneker(k2+ndi) = zero
      end do
C====================================================================
C
CE.   MATERIAL CONSTANTS 
!      E       = FUN(1,MAT) STARO CITANJE IZ DAT-A ZA DP 
!c     paneerselvam phd   (page 165 tens constants)
!c
      x = 1.03d-3
      a_l = 2          
      alfamm = 14.2d0  !17.1   
      h = 50.d0
      beta = 8.23d4     !1.34e4
      a_m  = one

	!stari
	!a1 = 70.3e2
	!a2 = 0.036e5 !7.5e4
	!a5 = -0.6046e5
	!novi
	a1 =  0.036d5	  
      a2 = -0.6046d5
      a3 = -1.81d5
      a4 = -0.906d5
	a5 = 70.3d2
      a6 = 0.
      a7 = 0.
      a8 = 0.
      a9 = 0.
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
      a(12)=14.2d0 !a(12)=alfamm
      a(13)=50.d0 !a(13)=h
      a(14)=8.23d4 !a(14)=beta
      a(15)=one !a(15)=m
      a(16)=-5.d-3 !a(16)=gama !
      a(17)=4.d-5 !a(17)=alpha_t !alpha_t  = 4.d-5      

!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!c***************      1) predictor phase      ***********************  
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      dtime=DT
      
      a_kapa0 = a_kapaP
      if (a_kapa0.lt.tolk) a_kapa0=tolk
      a_kapa = a_kapa0
      
      DO  I=1,6
      eplas(I)   = DEFPP(I)
      eplas0(I)   = DEFPP(I)
      END DO   
      
CE    TRIAL ELASTIC STRAINS
      DO I=1,6 
            DEFE(I)=DEF(I)-DEFPP(I)
      ENDDO
    
      call noviddsdde(ELAST,a,ntens,ndi,DEFE)
      call stressMM(a,ntens,ndi,DEFE,TAU)    
      call loadingf(f1,TAU,a,ntens,ndi,s_dev,a_j2,a_mu) ! 6.21
        
      f = abs(f1) - h*a_kapa0 

      !if (LDEFPP.eq.13) then
      !write(*,*) 'nn', TAU(2), DEF(2)
      !endif
	if (f.gt.zero) then 
!c------------------  end of elastic predictor ----------------------
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!c***************      2) corector phase      ************************  
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      do kewton = 1,newton

          call loadingf(f1,TAU,a,ntens,ndi,s_dev,a_j2,a_mu)
          f = abs(f1) - h*a_kapa
          a_kxl = x+a_kapa0**a_l
          dkapa0  =  (((f/beta)**1)/a_kxl)*dtime
          a_kapa = a_kapa0 + dkapa0
          a_kxl = x+a_kapa**a_l
          dkapa  =  (((f/beta)**1)/a_kxl)*dtime
      ! compute residual S          
          skonvergencija = abs(a_kapa - a_kapa0)
          
          do k1 = 1,ntens  
          eplasStari(k1)   = eplas(k1)
          d_eplas(k1)   =  (dkapa)*a_mu(k1)
          eplas(k1)   = eplas0(k1) + d_eplas(k1)
          Replas(k1) = eplas(k1) - eplasStari(k1)
          enddo
      ! compute residual R          
          deplas_int = (Replas(1)**2+Replas(2)**2+Replas(3)**2+
     1    two*Replas(4)**2+two*Replas(5)**2+two*Replas(6)**2)**0.5
     
       ! update for the next iteration kapa, epsilon, sigma   
          a_kapa0 = a_kapa
          DO I=1,6 
            DEFE(I)=DEF(I)-eplas(I)
          ENDDO
          call stressMM(a,ntens,ndi,DEFE,TAU)
!      check convergence          
       if ((skonvergencija.lt.tolk).and.(deplas_int.lt.tol2)) then
          goto 52
       endif
            
      enddo  !do kewton = 1,newton  
!c------------------  end of plastic corrector ----------------------  
          
      endif !if (f.gt.zero) then 

!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!c***************      3) save phase      ****************************  
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
CE    UPDATES FOR NEXT STEP

 52   continue      
      DO  I=1,6
      DEFPP(I)=eplas(I)
      TAU1(I)=TAU(I)
      END DO
      a_kapaP = a_kapa
      
      RETURN
      END
C
C  =====================================================================
C
      subroutine loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu)
      IMPLICIT NONE
      
!     ulazne promenljive   
      DOUBLE PRECISION a(17), alfa,gama       
      integer k1,k2,ntens,ndi
      DOUBLE PRECISION  stress(ntens),st(ntens)
!     ulazne promenljive

!     unutrasnje promenljive
      DOUBLE PRECISION one,zero,two,three,kroneker,s_i1,s_i2,s_i3
!     unutrasnje promenljive
      
!     izlazne promenljive
      DOUBLE PRECISION f1,a_j2,a_mu(ntens),s_dev(ntens) 
!     izlazne promenljive     

      dimension kroneker(ntens)
      
      f1 = 0
      alfa =  a(12)
      gama = a(16)
      s_i1=0
	s_i2=0
	s_i3=0
      one=1.0d0
      zero=0.0d0
      two=2.0d0
      three=3.0d0
      do k2=1, ndi
         kroneker(k2)     = one
         kroneker(k2+ndi) = zero
      end do
      do k2=1, 6
	   st(k2) = stress(k2)
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
      
      RETURN
      END
C
C  =====================================================================
C   
      subroutine noviddsdde(ddsdde,a,ntens,ndi,e_elas)
      IMPLICIT NONE
      
!     ulazne promenljive   
      DOUBLE PRECISION a(17),a1,a2,a3,a4,a5,a6,a7,a8,a9,e_elas(ntens)
      integer k1,k2,ntens,ndi 
!     ulazne promenljive  

!     izlazne promenljive
      DOUBLE PRECISION ddsdde(ntens,ntens)
!     izlazne promenljive

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
       ddsdde(1,6)= e_elas(6)*a4 + e_elas(6)*a4
 
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
       ddsdde(2,6)= e_elas(6)*a2 + e_elas(6)*a4 
     1  + e_elas(6)*a2 + e_elas(6)*a4
             
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
       ddsdde(3,6)= e_elas(6)*a2 + e_elas(6)*a4 
     1  + e_elas(6)*a2 + e_elas(6)*a4
       
       ddsdde(4,1)= 2*e_elas(4)*a2 + 2*e_elas(4)*a4
       ddsdde(4,2)= 2*e_elas(4)*a2 + 2*e_elas(4)*a4
       ddsdde(4,3)= 2*e_elas(4)*a4
       ddsdde(4,4)= a1 + a2*(e_elas(1) + e_elas(2)) 
	1 + a4*(e_elas(1) + e_elas(2) + e_elas(3))
       ddsdde(4,5)= e_elas(6)*a2
       ddsdde(4,6)= e_elas(5)*a2
       
       ddsdde(5,1)= 2*e_elas(5)*a2 + 2*e_elas(5)*a4
       ddsdde(5,2)= 2*e_elas(5)*a4
       ddsdde(5,3)= 2*e_elas(5)*a2 + 2*e_elas(5)*a4
       ddsdde(5,4)= e_elas(6)*a2
       ddsdde(5,5)= a1 + a2*(e_elas(1) + e_elas(3))
	1 + a4*(e_elas(1) + e_elas(2) + e_elas(3))
       ddsdde(5,6)= e_elas(4)*a2
       
       ddsdde(6,1)= 2*e_elas(6)*a4
       ddsdde(6,2)= 2*e_elas(6)*a2 + 2*e_elas(6)*a4
       ddsdde(6,3)= 2*e_elas(6)*a2 + 2*e_elas(6)*a4
       ddsdde(6,4)= e_elas(5)*a2
       ddsdde(6,5)= e_elas(4)*a2
       ddsdde(6,6)= a1 + a2*(e_elas(2) + e_elas(3))
	1 + a4*(e_elas(1) + e_elas(2) + e_elas(3))
      
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE stressMM(a,ntens,ndi,e_elas,stress)  
      IMPLICIT NONE
      
!     ulazne promenljive   
      DOUBLE PRECISION a(17),a1,a2,a3,a4,a5,a6,a7,a8,a9,e_elas(ntens)
      integer k1,ntens,ndi 
!     ulazne promenljive
         
!     unutrasnje invarijante         
      DOUBLE PRECISION e_i1n, e_i2n, e_i3n
!     unutrasnje invarijante
         
!     izlazne promenljive
      DOUBLE PRECISION stress(ntens)
!     izlazne promenljive

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
      
      RETURN
      END
C
C  =====================================================================
C