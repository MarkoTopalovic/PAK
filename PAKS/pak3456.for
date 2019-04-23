C=======================================================================
C
CS    CONCRETE MODEL WITH DAMAGE (22.07.2015.)
CE
      SUBROUTINE D3M56(TAU,DEF,IRAC,LPOCG,LPOC1,AU)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'paka.inc'
C
CE    PROGRAM FOR DEFINITION OF LOCATIONS AT INTEGRATION POINT LEVEL
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /REPERM/ MREPER(4)
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAU(6),DEF(6)
      DIMENSION AU(*)
      REAL AU
C
      IF(IDEBUG.GT.0) PRINT *, ' D3M56'
C  
C     REPERI ZA MATERIJALE U PAK-U
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)  
C           
C     REPERI ZA VELICINE U TRENUTKU T       
C
      LTAU = LPOCG
      LDEF = LTAU + 6
      LKAPA = LDEF + 6
      LDEFP = LKAPA + 2
      LUDL = LDEFP + 6
      LDEE = LUDL + 2
      LUPT = LDEE + 1      
C    
C     REPERI ZA VELICINE U TRENUTKU T+dT
C
      LTAU1 = LPOC1
      LDEF1 = LTAU1 + 6
      LKAPA1 = LDEF1 + 6
      LDEFP1 = LKAPA1 + 2
      LUDL1 = LDEFP1 + 6
      LDEE1 = LUDL1 + 2
      LUPT1 = LDEE1 + 1
C
C     POZIV SUBROUTINE TI3455
C      
      CALL TI3456(PLAS1(LDEFP1),PLAST(LDEFP),PLAS1(LKAPA1),
     * PLAST(LKAPA),TAU,PLAS1(LTAU1),PLAST(LTAU),DEF,
     * PLAS1(LDEF1),PLAST(LDEF),A(LFUN),IRAC,DT,KOR,MATE,
     * PLAST(LDEE),PLAS1(LDEE1),PLAS1(LUDL1),PLAS1(LUDL1+1),
     * PLAST(LUDL),PLAST(LUDL+1),PLAST(LUPT),PLAS1(LUPT1),AU(LCEL))
C
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI3456(DEFPT1,DEFPT,AKAPATT,AKAPAT,
     * TAU,TAU1,TAUT,DEF,DEF1,DEFT,PROPS,IRAC,DTIME,
     * KORAK,MATE,DEEEET,DEEEE1,UDL1T,UDL1C,UDLT,UDLC,
     * XT,XTDT,MCVEL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION DEF(6),DEFT(6),DEF1(6),TAU(6),TAUT(6),TAUDD(6),
     * DEFDS(6),TAUEF(6),DT(6,6),DEFPT(6),DEFPTT(6),TAUPR(3),
     * TAUDDGL(6),DDEF(6),QPTRT(3,3),PROPS(20,MATE),TAUEFGL(6),
     * TSG(6,6),TAUEF1(6),TAUDDPR(3),QPTRTD(3,3),AJED(6),TAUDDGL1(6),
     * TAUEFG(6),AKAPAT(2),AKAPATT(2),hTEN(2,3),HATEN(2),ANABFI(3),
     * QRESIZV(2),AKAPA(2),QRES(2),DEKAPA(2),DEDEF(6),TAU1(6),
     * DEDEF1(6),DDEKAPA(2),TSS(6,6),QIZV(2),CEPD(6,6),TAUDD1(6),
     * DEFPT1(6),TAUEFG1(6)
      DIMENSION MCVEL(*)
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DETT,NLM,KK
      COMMON /CDEBUG/ IDEBUG
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP  
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8      
      COMMON /CEPMAT/ INDCEP
      COMMON /ITERBR/ ITER  
      COMMON /ZAPREM/ ZAPRE      
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
C
      IF(IDEBUG.GT.0) PRINT *, ' TI3456'
C
C     INDIKATOR STAMPANJA PORUKA IMESS (=0-NE STAMPA; =1-STAMPA)
      IMESS=0
      NMM=NLM
      DEGG=0.95D0
      IF(ICVEL.EQ.1) NMM=MCVEL(NLM)
C
      TOLL=1.d-8
      TOLLL=-1.d-8
      DEEEE=0.d0
C
      E     = PROPS(1,MAT) !Young modulus
      V     = PROPS(2,MAT) !Poissont coefficient
      fcpri = PROPS(3,MAT) !maximal compressive stress
      ftpri = PROPS(4,MAT) !maximal tensile stress
      ACONC = PROPS(5,MAT) !ac
      DCCC  = PROPS(6,MAT) !dc
      ACONT = PROPS(7,MAT) !at
      DTTT  = PROPS(8,MAT) !dt
!       AALFFP = 0.23D0
      Ggc   = PROPS(12,MAT)
      Ggt   = PROPS(13,MAT)
      alch  = PROPS(14,MAT)
      AVRTEMP = PROPS(15,MAT)
!       GAMMA = 3.D0
      ESS0  = PROPS(16,MAT)
      AALFFP = PROPS(17,MAT)  
      AALFF = PROPS(18,MAT)  
      GAMMA = PROPS(19,MAT) 
      ADCR = PROPS(20,MAT)      
!       AKSI  = 0.D0
C      AKSI = PROPS(12,MAT)
!       AALFF = 0.12D0
C
      IF(ftpri.lt.1.d-8)THEN
      ftpri=0.1d0*fcpri
      PROPS(4,MAT)=0.1d0*fcpri
      ENDIF
      IF(aconc.lt.1.d-8)THEN
      PROPS(5,MAT)=5.d0
      aconc=5.d0
      ENDIF  
      IF(acont.lt.1.d-8)THEN
      PROPS(7,MAT)=.5d0
      acont=.5d0
      ENDIF          
      IF(DCCC.lt.1.d-8)THEN
      PROPS(6,MAT)=0.4d0
      DCCC=0.4d0
      ENDIF    
      IF(DTTT.lt.1.d-8)THEN
      PROPS(8,MAT)=0.51d0
      DTTT=0.51d0
      ENDIF   
      IF(AALFFP.lt.1.d-8)THEN
      PROPS(17,MAT)=0.23D0
      AALFFP=0.23D0
      ENDIF  
      IF(AALFF.lt.1.d-8)THEN
      PROPS(18,MAT)=0.12D0      
      AALFF=0.12D0
      ENDIF  
      IF(GAMMA.lt.1.d-8)THEN
      PROPS(19,MAT)=3.D0
      GAMMA=3.D0
      ENDIF   
      IF(ADCR.lt.1.d-8)THEN
      PROPS(20,MAT)=1.d0
      ADCR=1.D0
      ENDIF 
C
C     STRENGTH REDUCTION FACTOR	SAMO ZA PRIMER PRITUNI-10P0.DAT  
      SRF=1.	  
      IF(KORAK.GT.100) THEN
	     SRF=1.+(KORAK-100)/100
		 ENRED=((1.-AALFF)*fcpri-GAMMA*ftpri)/SRF
         AALFF=AALFF/SRF
		 fcpri=fcpri/SRF
		 ftpri=ftpri/SRF
		 Ggc=Ggc/SRF
		 Ggt=Ggt/SRF
		 GAMMA=((1.-AALFF)*fcpri-ENRED)/ftpri
	  ENDIF
C
      if(alch.lt.1.d-8)then
      alch=ZAPRE**(1.d0/3.d0)
      endif
C      
      gc=Ggc/alch
      gt=Ggt/alch
C
      CBCONC=(DLOG(1.D0-DCCC))/(DLOG((1.D0+ACONC)/2.D0/ACONC))        
      CBCONT=(DLOG(1.D0-DTTT))/(DLOG((1.D0+ACONT)-
     1        DSQRT(1.D0+ACONT**2))-DLOG(2.D0*ACONT))
      f0c   = fcpri*4.D0*ACONC/((1.D0+ACONC)**2)
      f0t   = ftpri   
C
      IF(KORAK.EQ.1) THEN
      AKAPAT(1)=0.D0
      AKAPAT(2)=0.D0
      CALL CLEAR (DEFPT,6)
      ENDIF     
C      
      fiikt1=1.D0+ACONT*(2.D0+ACONT)*AKAPAT(1)
      fiikc1=1.D0+ACONC*(2.D0+ACONC)*AKAPAT(2)
      ft1=(((1.D0/ACONT)*(1.D0+ACONT-
     1    DSQRT(fiikt1)))**(1.D0-CBCONT))*f0t*DSQRT(fiikt1)
      fc1=(((1.D0/ACONC)*(1.D0+ACONC-
     1    DSQRT(fiikc1)))**(1.D0-CBCONC))*f0c*DSQRT(fiikc1)
!       BBETTH = AKSI*AALFFP*ft0 
      GE=E/2.D0/(1.D0+V)
      AK0=E/3.D0/(1.D0-2.D0*V)
C      
      CALL CLEAR(ELAST,36)
      CALL MEL3EL(ELAST,E,V)       
C     
      IF(IRAC.EQ.2) RETURN           
C
      CALL CLEAR(TAUEFGL,6)
      CALL CLEAR(DEDEF1,6)  
      CALL CLEAR(DDEF,6)
C      
      DO 5 I=4,6
   5  DEFPT(I)=.5D0*DEFPT(I)   
C     Trial strain (Probna deformacija)
      DO 10 I=1,3
   10    DDEF(I)=DEF(I)-DEFPT(I)
      DO 15 I=4,6
   15    DDEF(I)=0.5D0*DEF(I)-DEFPT(I)
C     Probni efektivni napon TAUEF(Trial effective stress TAUEF)
      CALL CLEAR(TAUEF,6)
      CALL MNOZI1(TAUEF,ELAST,DDEF,6,6)
C
C     Glavni pravci
      CALL CLEAR(TAUPR,3)
      CALL CLEAR(QPTRT,9)  
      CALL CLEAR(TSG,36)
      CALL CLEAR(TSS,36)
      CALL CLEAR(ANABFI,3)
      CALL CLEAR(AKAPA,2)
      CALL CLEAR(TAUDDGL1,6)      
      CALL GLAVNTR(TAUEF,TAUPR)
      CALL GLAPR3TR(TAUEF,QPTRT,TAUPR)
C     Transformation matrices      
      CALL TRANSE(TSG,QPTRT)
      CALL TRANSS(TSS,QPTRT)
C      
      TAUEFGL(1)=TAUPR(1)
      TAUEFGL(2)=TAUPR(2)
      TAUEFGL(3)=TAUPR(3)  
      TAUEFGL(4)=0.d0
      TAUEFGL(5)=0.d0
      TAUEFGL(6)=0.d0
!       write(3,*)'TAUEFGL',TAUEFGL       
      DO 20 I=1,3
      J=I+3
      AJED(J)=0.D0
   20 AJED(I)=1.D0
C 
      AITR1=TAUEFGL(1)+TAUEFGL(2)+TAUEFGL(3) 
CS    SREDNJA UKUPNA I PLASTICNA DEFORMACIJA     
      EMT = (DDEF(1)+DDEF(2)+DDEF(3))/3.D0
!       EMP = (DEFPT(1)+DEFPT(2)+DEFPT(3))/3.D0
C     Trial deviatoric strain  
      CALL CLEAR(DEFDS,6)
      CALL CLEAR(TAUDD,6)    
      CALL CLEAR(TAUDD1,6)
      CALL CLEAR(TAUEFG,6)
      DO 30 I=1,3
   30    DEFDS(I)=DDEF(I)-(EMT)
      DO 40 I=4,6
   40    DEFDS(I)=DDEF(I)
C
      DO 50 I=1,6
   50    TAUDD(I)=DEFDS(I)*2.D0*GE
C   
      CALL JEDNA1(TAUDD1,TAUDD,6)
C 
      CALL CLEAR(TAUDDGL,6)
      CALL GLAVNTR(TAUDD,TAUDDPR)
      CALL GLAPR3TR(TAUDD,QPTRTD,TAUDDPR)
      TAUDDGL(1)=TAUDDPR(1)
      TAUDDGL(2)=TAUDDPR(2)
      TAUDDGL(3)=TAUDDPR(3) 
      TAUDDGL(4)=0.d0
      TAUDDGL(5)=0.d0
      TAUDDGL(6)=0.d0       
!       write(3,*)'TAUDDGL',TAUDDGL      
C
! !       IF(DABS(TAUEFGL(1)).GT.DABS(TAUEFGL(3))) THEN
! ! !       IF(AITR1.GT.0.D0) THEN
      IF(TAUEFGL(1).GT.TOLLL) THEN
        BBETT = (1.D0-AALFF)*fc1/ft1-(1.D0+AALFF)
      ELSE
	BBETT = GAMMA
      ENDIF
C
C     Ekvivalentni napon
      SEFF=0.d0
      SEFF=STRS_DOT(TAUDDGL,TAUDDGL)
      SEQV=DSQRT(3.D0/2.D0*SEFF)
      SEF=DSQRT(SEFF)
C      
      CALL JEDNA1(AKAPA,AKAPAT,2)
      CALL JEDNA1(TAUDDGL1,TAUDDGL,6)
      CALL JEDNA1(TAUEFG,TAUEFGL,6)      
      IF((DABS(TAUEFGL(1))+DABS(TAUEFGL(2))+
     1    DABS(TAUEFGL(3))).GT.0.D0) THEN
      Ar=0.5D0+0.5D0*AITR1/(DABS(TAUEFGL(1))+
     1   DABS(TAUEFGL(2))+DABS(TAUEFGL(3)))
      ELSE
      Ar=.5D0
      ENDIF
!       if(iter.eq.0) go to 301      
C     Calculate trial yield stress
!       write(3,*)'AALFF,AITR1,SEQV,BBETT,TAUEFGL,fc1',
!      1           AALFF,AITR1,SEQV,BBETT,TAUEFGL,fc1
      CALL FUNBET(FYD,AALFF,AITR1,SEQV,BBETT,TAUEFGL,fc1)
!       write(3,*)'FYDsmo ovde',FYD
      l=0
      ALAMBDA=0.D0
!       ADCR=0.D0
      DEEEE=DEEEET
      IF(FYD.LT.1.D-8) GOTO 301
C
      CALL CLEAR(DEKAPA,2)  
!       IF(DEEEET.GE.0.94999999999D0) THEN
!       DEEEE=0.95D0
!       goto 301  
!       ENDIF      
!       l=0
 300  CONTINUE
!        IF(AKAPA(1).LT.-1.d-8) STOP 'KAPA_T < 0'
!        IF(AKAPA(2).LT.-1.d-8) STOP 'KAPA_C < 0'   
       IF(AKAPA(1).LT.-1.d-8) THEN
       DEEEE=DEGG
       goto 301
       ENDIF
       IF(AKAPA(2).LT.-1.d-8) THEN
       DEEEE=DEGG
       goto 301
       ENDIF       
C
!       ALAMBDA=(AALFF*AITR1+SEQV+BBETT*TAUEFGL(1)-(1.D0-AALFF)*fc1)/
!      1 (9.D0*AK0*AALFFP*AALFF+DSQRT(6.D0)*GE*SEF1/DSQRT(BBETTH**2.D0+
!      1 SEF1**2.D0)+BBETT*(2.D0*GE*TAUDDGL(1)/
!      1 DSQRT(BBETTH**2.D0+SEF1**2.D0)+3.D0*AK0*AALFFP))
C
      fiikt1=1.D0+ACONT*(2.D0+ACONT)*AKAPA(1)
      fiikc1=1.D0+ACONC*(2.D0+ACONC)*AKAPA(2)
      ft1=(((1.D0/ACONT)*(1.D0+ACONT-
     1    DSQRT(fiikt1)))**(1.D0-CBCONT))*f0t*DSQRT(fiikt1)
      fc1=(((1.D0/ACONC)*(1.D0+ACONC-
     1    DSQRT(fiikc1)))**(1.D0-CBCONC))*f0c*DSQRT(fiikc1)
!       write(3,*)'TAUEFGL',TAUEFGL
!       write(3,*)'ft1,fc1,fiikt1,fiikc1',ft1,fc1,fiikt1,fiikc1
! !       IF(DABS(TAUEFGL(1)).GT.DABS(TAUEFGL(3))) THEN
! ! !       IF(AITR1.GT.0.D0) THEN
      IF(TAUEFG(1).GT.TOLLL) THEN
! ! !       IF(DABS(TAUEFGL(1)).GT.DABS(TAUEFGL(3))) THEN
	BBETT = (1.D0-AALFF)*fc1/ft1-(1.D0+AALFF)
      ELSE
	BBETT=GAMMA
      ENDIF
C
! !       write(3,*)'AALFF,AITR1,SEQV,BBETT,TAUEFGL(1),fc1,AK0,AALFFP,
! !      1GE,TAUDDGL(1),SEF',AALFF,AITR1,SEQV,BBETT,TAUEFGL(1),fc1,
! !      1AK0,AALFFP,GE,TAUDDGL(1),SEF
C
      ALAMBDA=(AALFF*AITR1+SEQV+BBETT*TAUEFGL(1)-
     1 (1.D0-AALFF)*fc1)/(9.D0*AK0*AALFFP*AALFF+
     1 DSQRT(6.D0)*GE+BBETT*(2.D0*GE*TAUEFGL(1)/
     1            SEF+3.D0*AK0*AALFFP-
     1            GE*AITR1/3.D0/SEF))
!       write(3,*)'ALAMBDA',ALAMBDA    
      IF(ALAMBDA.LT.0.D0.AND.IMESS.EQ.1) 
     1   write(*,*)'UPOZORENJE ALAMBDA<0',
     1               ' MAT=',MAT,' NMODM=',NMODM,' NLM=',NMM
C
!       ADCR=1.d0   
C     
      SEF1=SEF-2.D0*GE*ALAMBDA
C
      DO 70 I=1,6
  70 	TAUDDGL1(I)=TAUDDGL(I)*SEF1/SEF
C
      CALL CLEAR(ANABFI,3)
      DO 95 J=1,3
  95  ANABFI(J)=TAUDDGL1(J)/SEF1+AALFFP*AJED(J)
C  
      DO 80 I=1,6
	TAUEFG(I)=TAUEFGL(I)-ALAMBDA*(2.D0*GE*TAUEFGL(I)/
     1            SEF+3.D0*AK0*AALFFP*AJED(I)-
     1            GE*AITR1/3.D0/SEF*AJED(I))
  80  CONTINUE
      IF(ADCR.GT.1.D0) goto 301  
      IF(DEEEET.GE.DEGG) THEN
      DEEEE=DEGG
      goto 301  
      ENDIF
C
      DETTT=1.D0-(((1.D0+ACONT-
     1    DSQRT(fiikt1))/ACONT)**CBCONT)
      DECCC=1.D0-(((1.D0+ACONC-
     1    DSQRT(fiikc1))/ACONC)**CBCONC)
C
!       AITR11=TAUEFG(1)+TAUEFG(2)+TAUEFG(3) 
!       IF((DABS(TAUEFG(1))+DABS(TAUEFG(2))+
!      1    DABS(TAUEFG(3))).GT.0.D0) THEN      
!       Ar=.5D0+.5D0*AITR11/(DABS(TAUEFG(1))+
!      1   DABS(TAUEFG(2))+DABS(TAUEFG(3)))
!       ELSE
!       Ar=.5D0
!       ENDIF     
!       write(3,*)'Ar,AITR11',Ar,AITR11   
      DO 90 I=1,2
         DO 90 J=1,3      
  90  hTEN(I,J)=0.D0
      hTEN(1,1)=Ar*ft1*(1.D0-DETTT)/gt
      hTEN(2,3)=-(1.D0-Ar)*fc1*(1.D0-DECCC)/gc    
! ! !       write(3,*)'ft1,fc1,gt,gc',ft1,fc1,gt,gc  
C
! ! !       write(3,*)'ANABFI',ANABFI
      CALL CLEAR(HATEN,2)
      CALL MNOZI1(HATEN,hTEN,ANABFI,2,3)
C
!       write(3,*)'AKAPA,AKAPAT,ALAMBDA',AKAPA,AKAPAT,ALAMBDA
!       write(3,*)'HATEN',HATEN      
      DO 100 I=1,2
 100    QRES(I)=-AKAPA(I)+AKAPAT(I)+ALAMBDA*HATEN(I)
!       write(3,*)'QRES',QRES  
C
      QRESN2=STRS_DOT2(QRES,QRES)
      QRESN=DSQRT(QRESN2)
      DO 122 I=1,2
        IF(AKAPA(I).GT.ADCR) THEN 
        IF(IMESS.EQ.1)
     1  write(*,*)'POTPUNA DEGRADACIJA - K(',I,')>',ADCR,'=',AKAPA(I),
     1               ' MAT=',MAT,' NMODM=',NMODM,' NLM=',NMM
        write(3,*)'POTPUNA DEGRADACIJA - K(',I,')>',ADCR,'=',AKAPA(I),      
     1               ' MAT=',MAT,' NMODM=',NMODM,' NLM=',NMM
        AKAPA(I)=ADCR     
        DEEEE=DEGG    
        GOTO 301  
        ENDIF
  122 CONTINUE 
      IF(QRESN.LT.1.d-8) GOTO 301
C
      l=l+1
! !       write (3,*)'l=',l    
!       CALL QNUMIZV(QRESIZV,AKAPA,ALAMBDA,HATEN,ACONC,
!      1             ACONT,CBCONC,CBCONT,f0t,f0c,
!      1             gt,gc,TAUEFGL,TAUDDGL,AALFF,AITR1,GE,SEQV,
!      1             AK0,AALFFP,SEF,QRES)
      call clear(QIZV,2)        
      CALL QIZVOD(QIZV,AKAPA,ALAMBDA,ACONC,
     1            ACONT,CBCONC,CBCONT,f0t,f0c,
     1            ANABFI,gt,gc,Ar,TAUDDGL,AALFFP,GE,SEF,AK0,
     1            TAUEFGL,AALFF,BBETT,AITR1,SEQV) 
!       write(3,*)'QRESIZV',QRESIZV
! !       write(3,*)'QIZV',QIZV      
C
      CALL CLEAR (DDEKAPA,2)
      DO 105 I=1,2
 105     DDEKAPA(I)=-QRES(I)/QIZV(I)
C 
      DO 106 I=1,2
      IF(DDEKAPA(I).GT.0.1D0) THEN
      DDEKAPA(I)=0.1D0
      ENDIF
 106  CONTINUE     
C  
      DO 110 I=1,2
 110    DEKAPA(I)=DEKAPA(I)+DDEKAPA(I) 
      DO 127 I=1,2
        IF(DEKAPA(I).LT.-1.d-8) THEN
        IF(IMESS.EQ.1)
     1  write(*,*)'DEKAPA(',I,')<0=',DEKAPA(I), 
     1               ' MAT=',MAT,' NMODM=',NMODM,' NLM=',NMM
        ENDIF
  127 CONTINUE   
! ! !       write(*,*)'DEKAPA,AKAPAT,AKAPA',DEKAPA,AKAPAT,AKAPA 
      CALL CLEAR(AKAPA,2)      
      DO 120 I=1,2
 120  AKAPA(I)=AKAPAT(I)+DEKAPA(I)
! !       write (3,*)'NLM',NLM
! !       write (3,*)'QRES,QIZV',QRES,QIZV      
! !       write (3,*)'AKAPA,DEKAPA,AKAPAT',AKAPA,DEKAPA,AKAPAT    
      
      DO 121 I=1,2
        IF(AKAPA(I).GT.ADCR) THEN 
        IF(IMESS.EQ.1)
     1  write(*,*)'POTPUNA DEGRADACIJA - K(',I,')>',ADCR,'=',AKAPA(I),
     1               ' MAT=',MAT,' NMODM=',NMODM,' NLM=',NMM
        write(3,*)'POTPUNA DEGRADACIJA - K(',I,')>',ADCR,'=',AKAPA(I),  
     1               ' MAT=',MAT,' NMODM=',NMODM,' NLM=',NMM
        AKAPA(I)=ADCR
        DEEEE=DEGG
        GOTO 301     
        ENDIF
        IF(AKAPA(I).LT.-1.d-8) THEN 
        IF(IMESS.EQ.1)
     1  write(*,*)'NEGATIVNO - K(',I,')<0=',AKAPA(I),'SMANJITI KORAK',
     1               ' MAT=',MAT,' NMODM=',NMODM,' NLM=',NMM
        write(3,*)'NEGATIVNO - K(',I,')<0=',AKAPA(I),'SMANJITI KORAK',   
     1               ' MAT=',MAT,' NMODM=',NMODM,' NLM=',NMM
        IF(IMESS.EQ.1)
     1  write(*,*)'AKAPA=',AKAPA    
        DEEEE=DEGG
        GOTO 301     
        ENDIF        
  121 CONTINUE        
      if(l.gt.100) then 
        IF(IMESS.EQ.1)
     1  write(*,*)'QRESN',QRESN,
     1         ' MAT=',MAT,' NMODM=',NMODM,' NLM=',NMM
!       stop 'Prekoracio broj lokalnih iteracija>100'
!       write(*,*)'Prekoracio broj lokalnih iteracija>100'
      DEGG=DEGG-0.01D0
      ADCR=ADCR-0.01D0
      endif 
      GOTO 300
C      
 301  CONTINUE 
C
      CALL CLEAR(DEDEF,6)
      DO 130 I=1,3
         J=I+3
         DEDEF(I)=ALAMBDA*ANABFI(I)
 130     DEDEF(J)=0.D0        
!       write(3,*)'DEDEF',DEDEF
      CALL CLEAR(DEDEF1,6) 
      CALL MNOZI2(DEDEF1,TSS,DEDEF,6,6)
!       write(3,*)'DEDEF1',DEDEF1      
C
      ESS=ESS0+(1.D0-ESS0)*Ar
C
      DETTT=1.D0-(((1.D0+ACONT-
     1    DSQRT(fiikt1))/ACONT)**CBCONT)
      DECCC=1.D0-(((1.D0+ACONC-
     1    DSQRT(fiikc1))/ACONC)**CBCONC)
C
      INDGPR=0
      IF(INDGPR.EQ.1)THEN
      DO 654 I=1,6
      IF(TAUEFG(I).GT.0.d0) THEN
      TAUEFG1(I)=TAUEFG(I)*(1.D0-ESS*DETTT)
      ELSE
      TAUEFG1(I)=TAUEFG(I)*(1.D0-DECCC)
      ENDIF
  654 CONTINUE   
      ENDIF
C      
      CALL CLEAR(TAUEF1,6)
      CALL CLEAR(DEFPTT,6)
      IF(INDGPR.EQ.1)THEN
      CALL MNOZI2(TAUEF1,TSG,TAUEFG1,6,6)
      ELSE
      CALL MNOZI2(TAUEF1,TSG,TAUEFG,6,6)
      ENDIF
!       write(3,*)'TAUEF1',TAUEF1      
      CALL JEDNA1(DEF1,DEF,6)
      DO 150 I=1,6
  150    DEFPTT(I)=DEFPT(I)+DEDEF1(I)
      DO 160 I=4,6
  160    DEFPTT(I)=2.D0*DEFPTT(I)
!       write(3,*)'DEFPTT',DEFPTT
C
      ESS=ESS0+(1.D0-ESS0)*Ar
!       write(3,*)'ESS0,ESS,Ar',ESS0,ESS,Ar
!       ESS=1.D0
C
      IF(DEEEE.LT.DEGG) THEN
      DEEEE=1.D0-(1.D0-DECCC)*(1.D0-ESS*DETTT)
      ENDIF
!       DEEEE=1.D0-(1.D0-(1.d0-Ar)*DECCC)*(1.D0-Ar*DETTT)
!       DEEEE=0.D0
!       write(3,*)'DE,DEC,DET',DEEEE,DECCC,DETTT 
C
!       write(3,*)'ELAST',ELAST
      IF(ISKNP.NE.2.AND.INDCEP.EQ.0) THEN
      IF(DEEEE.GT.TOLL) THEN
      CALL CLEAR(CEPD,36)          
      CALL CEPDAMAGE(CEPD,AKAPA,ALAMBDA,ACONC,
     1                ACONT,CBCONC,CBCONT,f0t,f0c,
     1                ANABFI,gt,gc,Ar,TAUDDGL,AALFFP,GE,SEF,AK0,
     1                TAUEFGL,AALFF,BBETT,AITR1,ELAST,DEEEE,
     1                TAUEF1,TSS,TSG,TAUDD1)
      CALL CLEAR (ELAST,36)
      CALL JEDNA1 (ELAST,CEPD,36)
!       write(3,*)'CEPD',CEPD
      ENDIF
      ENDIF
!       DO 170 I=1,6
!          DO 170 J=1,6
!   170    ELAST(I,J)=ELAST(I,J)*(1.D0-DEEEE)
!       write(3,*)'ELAST'
!       DO 358 I=1,6
!   358    write(3,"(6E12.3)")(ELAST(I,J),J=1,6)    
C
      CALL CLEAR(TAU,6)
      CALL CLEAR(TAU1,6)  
      CALL CLEAR(DEFPT1,6)    
      CALL CLEAR(AKAPATT,2)     
       
      DO 180 I=1,6
  180    TAU(I)=(1.D0-DEEEE)*TAUEF1(I)
!       DO 190 I=1,6
!   190    TAU1(I)=TAU(I)  
      CALL JEDNA1(TAU1,TAU,6)
      CALL JEDNA1(DEFPT1,DEFPTT,6)      
!       write(3,*)'TAU1',TAU1  
      CALL JEDNA1(AKAPATT,AKAPA,2)
!       write(84,*)'AKAPATT',AKAPATT
      DEEEE1=DEEEE
C
!       write(555,*)'ft1,f0t,ACONT,CBCONT',ft1,f0t,ACONT,CBCONT
      UDL1C=fc1/f0c*(1.D0-DECCC)*100.d0
      UDL1T=ft1/f0t*(1.D0-DETTT)*100.d0
!       write(3,*)'UDL1C,UDL1T',UDL1C,UDL1T
      UDLC=UDL1C
      UDLT=UDL1T
C
      AITR111=TAUEFG(1)+TAUEFG(2)+TAUEFG(3) 
      asm=(1.d0-AALFF)*fc1-AALFF*AITR111-BBETT*TAUEFG(1)      
      aj2d=1.d0/6.d0*((TAUEF1(1)-TAUEF1(2))**2 +
     &                (TAUEF1(2)-TAUEF1(3))**2 +
     &                (TAUEF1(3)-TAUEF1(1))**2)+
     &                 TAUEF1(4)**2+TAUEF1(5)**2+TAUEF1(6)**2
      aj2dq=dsqrt(3.d0*aj2d)    
      ass=aj2dq
      sfactor=(1.d0-ass/asm)*100  
      if(sfactor.gt.100.d0) sfactor=100.d0
      if(sfactor.lt.0.d0) sfactor=0.d0
      xtdt=sfactor      
C      
      RETURN
      END 
c     ******************************************************************  
      SUBROUTINE FUNBET(FYD,AALFF,AITR1,SEQV0,BBETT,TAUDDGL,CCCC)   
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)          
C
      COMMON /CDEBUG/ IDEBUG
C
c     ********************* A   R   R   A   Y   S **********************
      DIMENSION TAUDDGL(6)
C
      IF(IDEBUG.GT.0) PRINT *, ' FUNBET'     
C
      FYD=0.D0
      FYD=(AALFF*AITR1+SEQV0+
     1    BBETT*TAUDDGL(1))/(1.D0-AALFF)-CCCC
C                                                    
      RETURN
      END
c     ******************************************************************        
C 
      FUNCTION STRS_DOT2(Q,S)    
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      
      DIMENSION Q(2),S(2)
C
      D=0.D0 
C      
      DO 10 I=1,2
 10       D=D+Q(I)*S(I)
C 
      STRS_DOT2=D    
C
      RETURN
      END  
c     ******************************************************************   
C
      SUBROUTINE QIZVOD(QIZV,AKAPA0,ALAMBDA1,ACONC,
     1                ACONT,CBCONC,CBCONT,f0t,f0c,
     1                ANABFI1,gt,gc,Ar1,TAUDDGL,AALFFP,GE,SEF,AK0,
     1                TAUEFGL,AALFF,BBETT,AITR1,SEQV)
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)          
C
      COMMON /CDEBUG/ IDEBUG
C
c     ********************* A   R   R   A   Y   S **********************
      DIMENSION AKAPA0(2),hTEN1(2,3),HATEN(2),ANABFI1(3),
     1         QIZV(2),hTENPOKAPA(2,3),hTENPOSIGT(3),hTENPOSIGC(3),
     1         HPOKAPA(2),QPOKAPA(2),QTPOSIG(3),QCPOSIG(3),SIGPOAL(3),
     1         TAUDDGL(6),AJED(6),TAUEFGL(6),QSIGKAP(2),ALPOKAPA(2),
     1         AJED13(3),AJED31(3),BLABLA1(3),BLABLA2(3),
     1         AJED66(6,6),AJED16(6),AJED6(6),EFPOSIG(3),QIZV2(2)
C
      IF(IDEBUG.GT.0) PRINT *, ' QIZVOD'    
C
      CALL CLEAR(QPOKAPA,2)
      CALL CLEAR(QSIGKAP,2)   
      CALL CLEAR(HATEN,2)   
      CALL CLEAR(ALPOKAPA,2) 
      CALL CLEAR(HPOKAPA,2)          
C
      TOLLL=0.d0
C
      AJED13(1)=1.D0
      AJED13(2)=0.D0
      AJED13(3)=0.D0
      AJED31(1)=0.D0
      AJED31(2)=0.D0
      AJED31(3)=1.D0      
C
      fiikt0=1.D0+ACONT*(2.D0+ACONT)*AKAPA0(1)
      fiikc0=1.D0+ACONC*(2.D0+ACONC)*AKAPA0(2)
      DETTTI=1.D0-(((1.D0+ACONT-
     1    DSQRT(fiikt0))/ACONT)**CBCONT)
      DECCCI=1.D0-(((1.D0+ACONC-
     1    DSQRT(fiikc0))/ACONC)**CBCONC)     
      ft11=f0t*(SQRT(fiikt0))*((1.D0/ACONT)*(1.D0+ACONT
     1    -SQRT(fiikt0)))**(1.D0-CBCONT)
      fc11=f0c*(SQRT(fiikc0))*((1.D0/ACONC)*(1.D0+ACONC
     1    -SQRT(fiikc0)))**(1.D0-CBCONC)
C
      DO 10 I=1,2
         DO 10 J=1,3      
  10  hTEN1(I,J)=0.D0
      hTEN1(1,1)=Ar1*ft11*(1.D0-DETTTI)/gt
      hTEN1(2,3)=-(1.D0-Ar1)*fc11*(1.D0-DECCCI)/gc
C      
      ftepokapa=f0t*((1.D0-CBCONT)*(((1.D0/ACONT*(1.D0+ACONT-
     1          DSQRT(fiikt0)))**(-CBCONT))*(-1.D0)*(2.D0+ACONT)
     1          /2.D0)+(((1.D0+ACONT-
     1          DSQRT(fiikt0))/ACONT)**(1.D0-CBCONT))*ACONT*
     1          (2.D0+ACONT)/2.D0/DSQRT(fiikt0))       
      fcepokapa=f0c*((1.D0-CBCONC)*(((1.D0/ACONC*(1.D0+ACONC-
     1          DSQRT(fiikc0)))**(-CBCONC))*(-1.D0)*(2.D0+ACONC)
     1          /2.D0)+(((1.D0+ACONC-
     1          DSQRT(fiikc0))/ACONC)**(1.D0-CBCONC))*ACONC*
     1          (2.D0+ACONC)/2.D0/DSQRT(fiikc0))    
C
      ftpokapa=f0t*(2.D0+ACONT)*((1.D0+ACONT)/2.D0/DSQRT(fiikt0)-1.d0)     
      fcpokapa=f0c*(2.D0+ACONC)*((1.D0+ACONC)/2.D0/DSQRT(fiikc0)-1.d0)
C
      DO 20 I=1,2
         DO 20 J=1,3      
  20  hTENPOKAPA(I,J)=0.D0
      hTENPOKAPA(1,1)=Ar1/gt*ftpokapa
      hTENPOKAPA(2,3)=-(1.D0-Ar1)/gc*fcpokapa    
!       write(3,*)'hTENPOKAPA',hTENPOKAPA(1,1),hTENPOKAPA(2,3)
C
      IF((DABS(TAUEFGL(1))+DABS(TAUEFGL(2))+DABS(TAUEFGL(3))).
     1    GT.0.D0) THEN
      IF(DABS(TAUEFGL(1)).GT.0.D0) THEN
      arposig1=((DABS(TAUEFGL(1))+DABS(TAUEFGL(2))+DABS(TAUEFGL(3)))-
     1          AITR1*TAUEFGL(1)/DABS(TAUEFGL(1)))/2.D0/
     1 ((DABS(TAUEFGL(1))+DABS(TAUEFGL(2))+DABS(TAUEFGL(3)))**2.D0)
      else
      arposig1=0.d0
      endif
      IF(DABS(TAUEFGL(2)).GT.0.D0) THEN  
      arposig2=((DABS(TAUEFGL(1))+DABS(TAUEFGL(2))+DABS(TAUEFGL(3)))-
     1          AITR1*TAUEFGL(2)/DABS(TAUEFGL(2)))/2.D0/
     1 ((DABS(TAUEFGL(1))+DABS(TAUEFGL(2))+DABS(TAUEFGL(3)))**2.D0)
      else
      arposig2=0.d0
      endif     
      IF(DABS(TAUEFGL(3)).GT.0.D0) THEN
      arposig3=((DABS(TAUEFGL(1))+DABS(TAUEFGL(2))+DABS(TAUEFGL(3)))-
     1          AITR1*TAUEFGL(3)/DABS(TAUEFGL(3)))/2.D0/
     1 ((DABS(TAUEFGL(1))+DABS(TAUEFGL(2))+DABS(TAUEFGL(3)))**2.D0)
      else
      arposig3=0.d0
      endif
      ELSE
      arposig1=0.D0
      arposig2=0.D0
      arposig3=0.D0
      ENDIF
! ! !       write(3,*)'arposig1,arposig2,arposig3',arposig1,arposig2,arposig3  
C
      hTENPOSIGT(1)=ft11*(1.D0-DETTTI)/gt*arposig1*ANABFI1(1)
      hTENPOSIGT(2)=ft11*(1.D0-DETTTI)/gt*arposig2*ANABFI1(1)
      hTENPOSIGT(3)=ft11*(1.D0-DETTTI)/gt*arposig3*ANABFI1(1)
      hTENPOSIGC(1)=fc11*(1.D0-DECCCI)/gc*arposig1*ANABFI1(3)
      hTENPOSIGC(2)=fc11*(1.D0-DECCCI)/gc*arposig2*ANABFI1(3)
      hTENPOSIGC(3)=fc11*(1.D0-DECCCI)/gc*arposig3*ANABFI1(3)      
! ! !       write(3,*)'hTENPOSIGT',hTENPOSIGT
! ! !       write(3,*)'hTENPOSIGC',hTENPOSIGC  
!       write(3,*)'ALAMBDA1',ALAMBDA1 
!       write(3,*)'TAUDDGL',TAUDDGL
!       write(3,*)'SEF',SEF      
C     
      IF(TAUEFGL(1).GT.TOLLL) THEN
!       IF(DABS(TAUEFGL(1)).GT.DABS(TAUEFGL(3))) THEN
!       IF(AITR1.GT.0.D0) THEN
      BETPOKAPT=(1.d0-AALFF)*(-fc11*ftepokapa)/(ft11**2.D0)
      BETPOKAPC=(1.d0-AALFF)*fcepokapa/ft11 
      ELSE
      BETPOKAPT=0.d0
      BETPOKAPC=0.d0
      ENDIF
C
      CALL CLEAR(HATEN,2) 
      CALL MNOZI1(HATEN,hTEN1,ANABFI1,2,3)  
      CALL CLEAR(HPOKAPA,2)     
      CALL MNOZI1(HPOKAPA,hTENPOKAPA,ANABFI1,2,3)     
C     
      CALL CLEAR(QPOKAPA,2) 
      DO 30 I=1,2
  30    QPOKAPA(I)=ALAMBDA1*HPOKAPA(I)-1.D0
!       DO 30 I=1,2
!   30    QPOKAPA(I)=-1.D0 ! jer je lambda*dH/dkapa<<1.D0   
C     
      CALL CLEAR(QTPOSIG,3) 
      CALL CLEAR(QCPOSIG,3)       
      DO 40 I=1,3
  40    QTPOSIG(I)=ALAMBDA1*(hTENPOSIGT(I)+hTEN1(1,1)/SEF
     1  *(AJED13(I)-TAUDDGL(1)*TAUDDGL(I)/(SEF**2)))
      DO 50 I=1,3
  50    QCPOSIG(I)=ALAMBDA1*(hTENPOSIGC(I)+hTEN1(2,3)/SEF
     1  *(AJED31(I)-TAUDDGL(3)*TAUDDGL(I)/(SEF**2)))  
! ! !       DO 40 I=1,3
! ! !   40    QTPOSIG(I)=ALAMBDA1*(hTEN1(1,1)/SEF
! ! !      1  *(AJED13(I)-TAUDDGL(1)*TAUDDGL(I)/(SEF**2)))
! ! !       DO 50 I=1,3
! ! !   50    QCPOSIG(I)=ALAMBDA1*(hTEN1(2,3)/SEF
! ! !      1  *(AJED31(I)-TAUDDGL(3)*TAUDDGL(I)/(SEF**2))) 
C
!       write(3,*)'QTPOSIGizQ' 
!       write(3,"(3E12.3)")(QTPOSIG(J),J=1,3)    
!       write(3,*)'QCPOSIGizQ' 
!       write(3,"(3E12.3)")(QCPOSIG(J),J=1,3)  
C
!       DO 40 I=1,3
!   40    QTPOSIG(I)=ALAMBDA1*hTENPOSIGT(I)
!       DO 50 I=1,3
!   50    QCPOSIG(I)=ALAMBDA1*hTENPOSIGC(I) 
C  
      CALL CLEAR(SIGPOAL,3)
      DO 60 I=1,3
        SIGPOAL(I)=2.D0*GE*TAUEFGL(I)/
     1            SEF+3.D0*AK0*AALFFP-
     1            GE*AITR1/3.D0/SEF
  60  CONTINUE
C
      CALL CLEAR(QSIGKAP,2)
      QSIGKAP(1)=DOT(QTPOSIG,SIGPOAL,3)
      QSIGKAP(2)=DOT(QCPOSIG,SIGPOAL,3)
C      
C     PROVERITI OVAJ IZVOD
C
! ! ! !       ALPOKAPT=(BETPOKAPT*TAUEFGL(1)*(9.D0*AK0*AALFFP*AALFF+
! ! ! !      1 DSQRT(6.D0)*GE+BBETT*(2.D0*GE*TAUDDGL(1)/SEF+
! ! ! !      1 3.D0*AK0*AALFFP))-(AALFF*AITR1+SEQV+BBETT*TAUEFGL(1)-
! ! ! !      1 (1.D0-AALFF)*fc11)*(BETPOKAPT*(2.D0*GE*TAUDDGL(1)/SEF+
! ! ! !      1 3.D0*AK0*AALFFP)))/((9.D0*AK0*AALFFP*AALFF+
! ! ! !      1 DSQRT(6.D0)*GE+BBETT*(2.D0*GE*TAUDDGL(1)/SEF+
! ! ! !      1 3.D0*AK0*AALFFP))**2.D0)
! ! ! ! C
! ! ! !       ALPOKAPC=((BETPOKAPC*TAUEFGL(1)-(1.D0-AALFF)*fcepokapa)*
! ! ! !      1 (9.D0*AK0*AALFFP*AALFF+
! ! ! !      1 DSQRT(6.D0)*GE+BBETT*(2.D0*GE*TAUDDGL(1)/SEF+
! ! ! !      1 3.D0*AK0*AALFFP))-((AALFF*AITR1+SEQV+BBETT*TAUEFGL(1)-
! ! ! !      1 (1.D0-AALFF)*fc11)*(BETPOKAPC*(2.D0*GE*TAUDDGL(1)/SEF+
! ! ! !      1 3.D0*AK0*AALFFP))))/((9.D0*AK0*AALFFP*AALFF+
! ! ! !      1 DSQRT(6.D0)*GE+BBETT*(2.D0*GE*TAUDDGL(1)/SEF+
! ! ! !      1 3.D0*AK0*AALFFP))**2.D0)  
      ALPOKAPT=((BETPOKAPT*TAUEFGL(1)*(9.D0*AK0*AALFFP*AALFF+
     1 DSQRT(6.D0)*GE+BBETT*(2.D0*GE*TAUEFGL(1)/
     1            SEF+3.D0*AK0*AALFFP-
     1 GE*AITR1/3.D0/SEF)))-(AALFF*AITR1+SEQV+BBETT*TAUEFGL(1)-
     1 (1.D0-AALFF)*fc11)*(BETPOKAPT*(2.D0*GE*TAUEFGL(1)/
     1            SEF+3.D0*AK0*AALFFP-
     1 GE*AITR1/3.D0/SEF)))/((9.D0*AK0*AALFFP*AALFF+
     1 DSQRT(6.D0)*GE+BBETT*(2.D0*GE*TAUEFGL(1)/
     1            SEF+3.D0*AK0*AALFFP-
     1 GE*AITR1/3.D0/SEF))**2.D0)
C
      ALPOKAPC=((BETPOKAPC*TAUEFGL(1)-(1.D0-AALFF)*fcepokapa)*
     1 (9.D0*AK0*AALFFP*AALFF+
     1 DSQRT(6.D0)*GE+BBETT*(2.D0*GE*TAUEFGL(1)/
     1            SEF+3.D0*AK0*AALFFP-
     1 GE*AITR1/3.D0/SEF))-((AALFF*AITR1+SEQV+BBETT*TAUEFGL(1)-
     1 (1.D0-AALFF)*fc11)*(BETPOKAPC*(2.D0*GE*TAUEFGL(1)/
     1            SEF+3.D0*AK0*AALFFP-
     1 GE*AITR1/3.D0/SEF))))/((9.D0*AK0*AALFFP*AALFF+
     1 DSQRT(6.D0)*GE+BBETT*(2.D0*GE*TAUEFGL(1)/
     1            SEF+3.D0*AK0*AALFFP-
     1 GE*AITR1/3.D0/SEF))**2.D0)  
! ! !       EFPOKAPT=-TAUEFGL(1)*fc11*ftpokapa/(ft11**2.D0)
! ! !       EFPOKAPC=(TAUEFGL(3)/ft11-1.D0)*fcpokapa
! ! !       
! ! !       DO 61 I=1,3
! ! !    61   EFPOSIG(I)=(AALFF+DSQRT(3.D0/2.D0)*
! ! !      1             TAUDDGL(I)/SEF+BBETT*AJED13(I))    
! ! ! C
! ! !       EFSIGLA=DOT(EFPOSIG,SIGPOAL,3) 
! ! ! C
! ! !       ALPOKAPT=-EFPOKAPT/EFSIGLA    
! ! !       ALPOKAPC=-EFPOKAPC/EFSIGLA        
C      
      ALPOKAPA(1)=ALPOKAPT
      ALPOKAPA(2)=ALPOKAPC  
C        
      CALL CLEAR(QIZV2,2)
      CALL CLEAR(QIZV,2)
      DO 65 I=1,2
  65    QIZV2(I)=(QSIGKAP(I)+HATEN(I))*ALPOKAPA(I)
      DO 70 I=1,2
  70    QIZV(I)=QPOKAPA(I)+(QSIGKAP(I)+HATEN(I))*ALPOKAPA(I)
! !       write(3,*)'QPOKAPA',QPOKAPA  
! !       write(3,*)'QSIGKAP',QSIGKAP
! !       write(3,*)'HATEN',HATEN
! !       write(3,*)'ALPOKAPA',ALPOKAPA
! !       write(3,*)'QIZV2',QIZV2
C
      RETURN
      END  
c     ******************************************************************   
C
      SUBROUTINE CEPDAMAGE(CEPD,AKAPA0,ALAMBDA1,ACONC,
     1                ACONT,CBCONC,CBCONT,f0t,f0c,
     1                ANABFI1,gt,gc,Ar1,TAUDDGLL,AALFFP,GE,SEF,AK0,
     1                TAUEFGLL,AALFF,BBETT,AITR1,ELASTT,DEEEE,
     1                TAUEFFF,TSSS,TSGG,TAUDDDD)
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)          
C
      COMMON /CDEBUG/ IDEBUG
C
c     ********************* A   R   R   A   Y   S **********************
      DIMENSION AKAPA0(2),hTEN1(2,3),HATEN(2),ANABFI1(3),ELASTT(6,6),
     1         CEPD(6,6),hTENPOKAPA(2,3),hTENPOSIGT(3),hTENPOSIGC(3),
     1         HPOKAPA(2),QPOKAPA(2),QTPOSIG(6),QCPOSIG(6),SIGPOAL(3),
     1         TAUDDGLL(6),AJED(6),TAUEFGLL(6),QSIGKAP(2),ALPOKAPA(2),
     1         AJED66(6,6),AJED16(6),AJED6(6),SIGDSIG(6,6),CEPEF(6,6),
     1         CINVEL(6,6),TAUEFFF(6),DSIG(6),DKAP(2),AKAPSIG(2,6),
     1        ANALAMP(6,6),ANABFI(3),ALAMPOSIG(6),TSSS(6,6),TSGG(6,6),
     1         AKAP11(2),AKAPTPOSIG(6),AKAPCPOSIG(6),QPOSIG(2,6),
     1         EFPOKAP(2),EFPOSIG(6),ALAMPOSIGL(6),FPKAPQPSIG(6),
     1         AKAP12(2),LJA(6),MJA(6),LJAA(6),MJAA(6),CEPD1(6,6),
     1         ANABFII(6),AKAP112(2,6),AKAP222(2,6),ANALAMPGL(6,6),
     1         AKAPSIGG(2,6),EFPOSIG2(6,6),AJED13(3),AJED31(3),
     1         TAUDDDD(6),QPOSI(2,6),QTPOSI(6),QCPOSI(6),CEPD2(6,6),
     1         QQTPOSIG(6),QQCPOSIG(6),CEPEF3CL(6,6),AMATR(6,6),
     1         CEPEF3(6,6),AKAP2221(2)
C
      IF(IDEBUG.GT.0) PRINT *, ' CEPDAMAGE'    
C
      CALL CLEAR(QPOKAPA,2)
      CALL CLEAR(QSIGKAP,2)   
      CALL CLEAR(HATEN,2)   
      CALL CLEAR(ALPOKAPA,2)   
      CALL CLEAR(HPOKAPA,2) 
      CALL CLEAR(FPKAPQPSIG,6)
      CALL CLEAR(ALAMPOSIG,6)
      CALL CLEAR(AKAP222,12)
C
      DO 210 I=1,6
         DO 210 J=1,6
      IF(I.EQ.J) THEN
      AJED66(I,J)=1.D0
      ELSE
      AJED66(I,J)=0.D0
      ENDIF
  210 CONTINUE    
      DO 211 I=1,6
      IF(I.EQ.1) THEN
      AJED16(I)=1.D0
      ELSE
      AJED16(I)=0.D0
      ENDIF
  211 CONTINUE  
      DO 212 I=1,6
      IF(I.LE.3) THEN
      AJED6(I)=1.D0
      ELSE
      AJED6(I)=0.D0
      ENDIF
  212 CONTINUE    
      AJED13(1)=1.D0
      AJED13(2)=0.D0
      AJED13(3)=0.D0
      AJED31(1)=0.D0
      AJED31(2)=0.D0
      AJED31(3)=1.D0 
      CALL CLEAR(AMATR,36)
      AMATR(1,1)=2.D0/3.D0
      AMATR(2,2)=2.D0/3.D0
      AMATR(3,3)=2.D0/3.D0
      AMATR(1,2)=-1.D0/3.D0
      AMATR(1,3)=-1.D0/3.D0
      AMATR(2,3)=-1.D0/3.D0
      AMATR(2,1)=-1.D0/3.D0
      AMATR(3,2)=-1.D0/3.D0
      AMATR(3,1)=-1.D0/3.D0
      AMATR(4,4)=1.D0
      AMATR(5,5)=1.D0
      AMATR(6,6)=1.D0
C
      fiikt0=1.D0+ACONT*(2.D0+ACONT)*AKAPA0(1)
      fiikc0=1.D0+ACONC*(2.D0+ACONC)*AKAPA0(2)
      DETTTI=1.D0-(((1.D0+ACONT-
     1    DSQRT(fiikt0))/ACONT)**CBCONT)
      DECCCI=1.D0-(((1.D0+ACONC-
     1    DSQRT(fiikc0))/ACONC)**CBCONC)     
      ft11=f0t*(SQRT(fiikt0))*((1.D0/ACONT)*(1.D0+ACONT
     1    -SQRT(fiikt0)))**(1.D0-CBCONT)
      fc11=f0c*(SQRT(fiikc0))*((1.D0/ACONC)*(1.D0+ACONC
     1    -SQRT(fiikc0)))**(1.D0-CBCONC)
C
      DO 10 I=1,2
         DO 10 J=1,3      
  10  hTEN1(I,J)=0.D0
      hTEN1(1,1)=Ar1*ft11*(1.D0-DETTTI)/gt
      hTEN1(2,3)=-(1.D0-Ar1)*fc11*(1.D0-DECCCI)/gc
C      
      ftepokapa=f0t*((1.D0-CBCONT)*(((1.D0/ACONT*(1.D0+ACONT-
     1          DSQRT(fiikt0)))**(-CBCONT))*(-1.D0)*(2.D0+ACONT)
     1          /2.D0/DSQRT(fiikt0))*DSQRT(fiikt0)+(((1.D0+ACONT-
     1          DSQRT(fiikt0))/ACONT)**(1.D0-CBCONT))*ACONT*
     1          (2.D0+ACONT)/2.D0/DSQRT(fiikt0))
      fcepokapa=f0c*((1.D0-CBCONC)*(((1.D0/ACONC*(1.D0+ACONC-
     1          DSQRT(fiikc0)))**(-CBCONC))*(-1.D0)*(2.D0+ACONC)
     1          /2.D0/DSQRT(fiikc0))*DSQRT(fiikc0)+(((1.D0+ACONC-
     1          DSQRT(fiikc0))/ACONC)**(1.D0-CBCONC))*ACONC*
     1          (2.D0+ACONC)/2.D0/DSQRT(fiikc0))    
C
      ftpokapa=f0t*((1.D0+ACONT)*(2.D0+ACONT)/2.D0/DSQRT(fiikt0)-
     1         (2.D0+ACONT))
      fcpokapa=f0c*((1.D0+ACONC)*(2.D0+ACONC)/2.D0/DSQRT(fiikc0)-
     1         (2.D0+ACONC)) 
C
      DO 20 I=1,2
         DO 20 J=1,3      
  20  hTENPOKAPA(I,J)=0.D0
      hTENPOKAPA(1,1)=Ar1/gt*ftpokapa
      hTENPOKAPA(2,3)=-(1.D0-Ar1)/gc*fcpokapa    
C
C     PROVERITI + -
C
      arposig11=(SIGN(1.D0,TAUEFGLL(1))*(DABS(TAUEFGLL(1))+
     1 DABS(TAUEFGLL(2))+DABS(TAUEFGLL(3)))-SIGN(1.D0,TAUEFGLL(1))*
     1 (TAUEFGLL(1)+TAUEFGLL(2)+TAUEFGLL(3)))/((DABS(TAUEFGLL(1))+
     1 DABS(TAUEFGLL(2))+DABS(TAUEFGLL(3)))**2)/2.D0
      arposig23=(SIGN(1.D0,TAUEFGLL(3))*(DABS(TAUEFGLL(1))+
     1 DABS(TAUEFGLL(2))+DABS(TAUEFGLL(3)))-SIGN(1.D0,TAUEFGLL(3))*
     1 (TAUEFGLL(1)+TAUEFGLL(2)+TAUEFGLL(3)))/((DABS(TAUEFGLL(1))+
     1 DABS(TAUEFGLL(2))+DABS(TAUEFGLL(3)))**2)/2.D0  
C
      hTENPOSIGT(1)=ft11*(1.D0-DETTTI)/gt*arposig11*ANABFI1(1)
      hTENPOSIGT(2)=ft11*(1.D0-DETTTI)/gt*arposig11*ANABFI1(1)
      hTENPOSIGT(3)=ft11*(1.D0-DETTTI)/gt*arposig11*ANABFI1(1)
      hTENPOSIGC(1)=arposig23*fc11*(1.D0-DECCCI)/gc*ANABFI1(3)
      hTENPOSIGC(2)=arposig23*fc11*(1.D0-DECCCI)/gc*ANABFI1(3)
      hTENPOSIGC(3)=arposig23*fc11*(1.D0-DECCCI)/gc*ANABFI1(3)
!       write(3,*)'hTENPOSIGT',hTENPOSIGT 
!       write(3,*)'hTENPOSIGC',hTENPOSIGC      
C     
!       BETPOKAPT=(1.d0-AALFF)*(-fc11*ftepokapa)/ft11**2.D0
!       BETPOKAPC=(1.d0-AALFF)*(ft11*fcepokapa)/ft11**2.D0      
C
      CALL MNOZI1(HATEN,hTEN1,ANABFI1,2,3)  
      CALL MNOZI1(HPOKAPA,hTENPOKAPA,ANABFI1,2,3)  
C     
      DO 30 I=1,2
  30    QPOKAPA(I)=ALAMBDA1*HPOKAPA(I)-1.D0
!       DO 30 I=1,2
!   30    QPOKAPA(I)=-1.D0 ! jer je lambda*dH/dkapa<<1.D0   
C     
!       DO 40 I=1,3
!         J=I+3
!         QTPOSIG(J)=0.D0
!   40    QTPOSIG(I)=ALAMBDA1*hTENPOSIGT(I)
!       DO 50 I=1,3
!         J=I+3
!         QCPOSIG(J)=0.D0      
!   50    QCPOSIG(I)=ALAMBDA1*hTENPOSIGC(I)  
      DO 40 I=1,3
        J=I+3
        QTPOSIG(J)=0.D0      
  40    QTPOSIG(I)=ALAMBDA1*(hTEN1(1,1)/SEF
     1  *(AJED13(I)-TAUDDGLL(1)*TAUDDGLL(I)/(SEF**2)))
      DO 50 I=1,3
        J=I+3
        QCPOSIG(J)=0.D0          
  50    QCPOSIG(I)=ALAMBDA1*(hTEN1(2,3)/SEF
     1  *(AJED31(I)-TAUDDGLL(3)*TAUDDGLL(I)/(SEF**2)))    
C
      DO 51 I=1,2
         DO 51 J=1,6
      IF(I.EQ.1)THEN
      QPOSIG(I,J)=QTPOSIG(J)
      ELSE
      QPOSIG(I,J)=QCPOSIG(J)
      ENDIF
  51  CONTINUE  
C
      DO 53 I=1,2
         DO 53 J=1,6
      IF(I.EQ.1)THEN
      QQTPOSIG(J)=QTPOSIG(J)/QPOKAPA(I)
      ELSE
      QQCPOSIG(J)=QCPOSIG(J)/QPOKAPA(I)
      ENDIF
  53  CONTINUE  
!       write(3,*)'QTPOSIG' 
!       write(3,"(6E12.3)")(QTPOSIG(J),J=1,6)    
!       write(3,*)'QCPOSIG' 
!       write(3,"(6E12.3)")(QCPOSIG(J),J=1,6)  
      CALL CLEAR (QTPOSI,6)
      CALL CLEAR (QCPOSI,6)      
      CALL MNOZI2(QTPOSI,TSGG,QQTPOSIG,6,6)
      CALL MNOZI2(QCPOSI,TSGG,QQCPOSIG,6,6)   
!       write(3,*)'QTPOSI' 
!       write(3,"(6E12.3)")(QTPOSI(J),J=1,6)    
!       write(3,*)'QCPOSI' 
!       write(3,"(6E12.3)")(QCPOSI(J),J=1,6)       
C
      CALL CLEAR (QPOSI,12)
      DO 52 I=1,2
         DO 52 J=1,6
      IF(I.EQ.1)THEN
      QPOSI(I,J)=QTPOSI(J)
      ELSE
      QPOSI(I,J)=QCPOSI(J)
      ENDIF
  52  CONTINUE
C      
      EFPOKAPT=-TAUEFGLL(1)*fc11*ftpokapa/(ft11**2.D0)
      EFPOKAPC=(TAUEFGLL(3)/ft11-1.D0)*fcpokapa
!       EFPOKAPT=-TAUEFGL(1)*(1.D0-AALFF)*fc11*ftpokapa/(ft11**2.D0)
!       EFPOKAPC=-(1.D0-AALFF)*fcpokapa
      EFPOKAP(1)=EFPOKAPT/QPOKAPA(1)
      EFPOKAP(2)=EFPOKAPC/QPOKAPA(2)
C
      CALL MNOZM1(FPKAPQPSIG,EFPOKAP,QPOSIG,1,6,2)
!       write(3,*)'FPKAPQPSIG',FPKAPQPSIG      
C
!       CALL MNOZM1(EFPKAPH,EFPOKAP,HATEN,1,1,2)   
      EFPKAPH=DOT(EFPOKAP,HATEN,2)
!       write(3,*)'EFPKAPH',EFPKAPH
C
      CALL CLEAR(EFPOSIG,6)
!       write(3,*)'TAUDDGL,AJED6,AJED16',TAUDDGL,AJED6,AJED16
      DO 60 I=1,6
   60   EFPOSIG(I)=AALFF*AJED6(I)+DSQRT(3.D0/2.D0)*
     1             TAUDDGLL(I)/SEF+BBETT*AJED16(I)
C     
!       DO 61 I=1,3
!    61   EFPOSIG(I)=(AALFF+DSQRT(3.D0/2.D0)*
!      1             TAUDDGL(I)/SEF+BBETT*AJED13(I))/(1.D0-AALFF)        
C
!       write(3,*)'EFPOSIG',EFPOSIG
      DO 61 I=1,6
   61   ALAMPOSIGL(I)=(EFPOSIG(I)-FPKAPQPSIG(I))/EFPKAPH  
C
!       write(3,*)'ALAMPOSIGL' 
!       write(3,"(6E12.3)")(ALAMPOSIGL(J),J=1,6) 
      CALL CLEAR (ALAMPOSIG,6)
      CALL MNOZI2(ALAMPOSIG,TSGG,ALAMPOSIGL,6,6)   
C
      DKAPT=CBCONT*((1.D0/ACONT*(1.D0+ACONT-DSQRT(fiikt0)))**
     1     (CBCONT-1.D0))*(2.D0+ACONT)/2.D0/DSQRT(fiikt0)
      DKAPC=CBCONC*((1.D0/ACONC*(1.D0+ACONC-DSQRT(fiikc0)))**
     1     (CBCONC-1.D0))*(2.D0+ACONC)/2.D0/DSQRT(fiikc0)
      DKAP(1)=DKAPT
      DKAP(2)=DKAPC
C
!       CALL MNOZM1(AKAP112,QPOSIG,TSS,2,6,6)
      DO 58 I=1,2
      AKAP2221(I)=HATEN(I)/QPOKAPA(I)
  58  CONTINUE 
C  
      CALL MNOZM1(AKAP222,AKAP2221,ALAMPOSIG,2,6,1)
!       CALL MNOZI2(AKAP11,TSS,QTPOSIG,6,6)
!       CALL MNOZI2(AKAP12,TSS,QCPOSIG,6,6)    
!       DO 54 I=1,6
!       AKAPTPOSIG(I)=AKAP11(I)+HATEN(1)*ALAMPOSIG(I)
!       AKAPCPOSIG(I)=AKAP12(I)+HATEN(2)*ALAMPOSIG(I) 
!   54  CONTINUE
! ! ! ! !       write(3,*)'QPOSI'
! ! ! ! !       DO 355 I=1,2
! ! ! ! !   355    write(3,"(6E12.3)")(QPOSI(I,J),J=1,6)  
! ! ! ! !       write(3,*)'AKAP222'
! ! ! ! !       DO 354 I=1,2
! ! ! ! !   354    write(3,"(6E12.3)")(AKAP222(I,J),J=1,6)   
C  
      DO 54 I=1,2
         DO 54 J=1,6
      AKAPSIG(I,J)=-QPOSI(I,J)-AKAP222(I,J)
  54  CONTINUE
C
!       write(3,*)'TAUDD' 
!       write(3,"(6E12.3)")(TAUDD(J),J=1,6)   
      DO 56 I=1,6
      ANABFI(I)=TAUDDDD(I)/SEF+AALFFP*AJED6(I)
  56  CONTINUE
!       write(3,*)'ANABFI' 
!       write(3,"(6E12.3)")(ANABFI(J),J=1,6)   
      CALL CLEAR(ANALAMP,36)
      CALL CLEAR(CEPEF,36)      
C
!       write(3,*)'ALAMPOSIG'
!       write(3,"(6E12.3)")(ALAMPOSIG(J),J=1,6)       
      CALL MNOZM1(ANALAMP,ANABFI,ALAMPOSIG,6,6,1)   
! ! ! !       write(3,*)'DKAP' 
! ! ! !       write(3,"(2E12.3)")(DKAP(J),J=1,2)  
! ! ! !       write(3,*)'AKAPSIG'
! ! ! !       DO 356 I=1,2
! ! ! !   356    write(3,"(6E12.3)")(AKAPSIG(I,J),J=1,6)   
      CALL CLEAR(DSIG,6)
      CALL MNOZM1(DSIG,DKAP,AKAPSIG,1,6,2)
! ! ! !       write(3,*)'DSIG' 
! ! ! !       write(3,"(6E12.3)")(DSIG(J),J=1,6)    
! ! ! !       write(3,*)'TAUEFFF' 
! ! ! !       write(3,"(6E12.3)")(TAUEFFF(J),J=1,6)
      CALL CLEAR(SIGDSIG,36)
      CALL MNOZM1(SIGDSIG,TAUEFFF,DSIG,6,6,1)
! ! ! C
! ! !       write(3,*)'ELAST'
! ! !       DO 369 I=1,6
! ! !   369    write(3,"(6E12.3)")(ELASTT(I,J),J=1,6)  
      CALL JEDNA1(CINVEL,ELASTT,36)   
      CALL MINV(CINVEL,6,DUM,LJA,MJA)    
C
!       write(3,*)'ANALAMP'
!       DO 317 I=1,6
!   317    write(3,"(6E12.3)")(ANALAMP(I,J),J=1,6) 
!       DO 65 I=1,6
!          DO 65 J=1,6
!   65    CEPEF(I,J)=CINVEL(I,J)+ANALAMP(I,J)+ALAMBDA1*1.D0/SEF*
!      1  (AJED66(I,J)-1.D0/SEF/SEF*TAUDDDD(I)*TAUDDDD(J))
!       write(3,*)'TAUDDDD',TAUDDDD
      DO 64 I=1,6
         DO 64 J=1,6
  64    CEPEF3CL(I,J)=ALAMBDA1/SEF*
     1  (AJED66(I,J)-TAUDDDD(I)*TAUDDDD(J)/(SEF**2))
      CALL CLEAR(CEPEF3,36)
      CALL MNOZM1(CEPEF3,CEPEF3CL,AMATR,6,6,6)
!       write(3,*)'CEPEF3CL'
!       DO 328 I=1,6
!   328    write(3,"(6E12.3)")(CEPEF3CL(I,J),J=1,6)   
!        write(3,*)'CEPEF3'
!       DO 329 I=1,6
!   329   write(3,"(6E12.3)")(CEPEF3(I,J),J=1,6)    
C 
      CALL CLEAR(CEPEF,36) 
      DO 65 I=1,6
         DO 65 J=1,6
  65    CEPEF(I,J)=CINVEL(I,J)
!       write(3,*)'CINVEL'
!       DO 368 I=1,6
!   368    write(3,"(6E12.3)")(CINVEL(I,J),J=1,6)    
!       write(3,*)'CEPEF'
!       DO 378 I=1,6
!   378    write(3,"(6E12.3)")(CEPEF(I,J),J=1,6)      
C     
      CALL MINV(CEPEF,6,DUMM,LJAA,MJAA)  
! ! !       write(3,*)'CEPEFINV'
! ! !       DO 328 I=1,6
! ! !   328    write(3,"(6E12.3)")(CEPEF(I,J),J=1,6)        
! ! ! !       write(3,*)'SIGDSIG'
! ! ! !       DO 357 I=1,6
! ! ! !   357    write(3,"(6E12.3)")(SIGDSIG(I,J),J=1,6)       
!       DO 70 I=1,6
!          DO 70 J=1,6
!   70    CEPD1(I,J)=(1.D0-DEEEE)*AJED66(I,J)-SIGDSIG(I,J)
      CALL CLEAR(CEPD1,36) 
      DO 70 I=1,6
         DO 70 J=1,6
  70     CEPD1(I,J)=(1.D0-DEEEE)*AJED66(I,J)
C
      CALL CLEAR(CEPD2,36)
      CALL CLEAR(CEPD,36)      
      CALL MNOZM1(CEPD2,CEPD1,CEPEF,6,6,6)
!       CALL JEDNA1(CEPD,CEPD2,36)
! ! !       write(3,*)'CEPD2'
! ! !       DO 307 I=1,6
! ! !   307    write(3,"(6E12.3)")(CEPD2(I,J),J=1,6)          
      DO 80 I=1,6
         DO 80 J=1,6
      IF(I.NE.J) THEN
      CEPD(I,J)=(CEPD2(I,J)+CEPD2(J,I))/2.D0
      ELSE
      CEPD(I,J)=CEPD2(I,J)
      ENDIF
  80  CONTINUE    
! ! !       write(3,*)'CEPD'
! ! !       DO 358 I=1,6
! ! !   358    write(3,"(6E12.3)")(CEPD(I,J),J=1,6)         
C
      RETURN
      END 
C
c     ******************************************************************   
      SUBROUTINE QNUMIZV(QRESIZV,AKAPA0,ALAMBDA1,HATEN0,ACONC,
     1                   ACONT,CBCONC,CBCONT,f0t,f0c,
     1                   gt,gc,TAUEFGL,TAUDDGL,AALFF,AITR1,GE,SEQV,
     1                   AK0,AALFFP,SEF,QRES22)
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)          
C
      COMMON /CDEBUG/ IDEBUG
C
c     ********************* A   R   R   A   Y   S **********************
      DIMENSION AKAPA0(2),AKAPA1(2),hTEN1(2,3),HATEN1(2),ANABFI(3),
     1          QRESIZV(2),QRES0(2),QRES1(2),HATEN0(2),TAUEFGL(6),
     1          TAUDDGL(6),TAUEFG(6),AJED(6),TAUDDGL1(6),QRES22(2)
C
      IF(IDEBUG.GT.0) PRINT *, ' QNUMIZV'    
C
      TOLEPS=1.d-8
      TOLLL=0.d0
C
      DO 5 I=1,3
      J=I+3
      AJED(J)=0.D0
   5  AJED(I)=1.D0
C
      DO 10 I=1,2
  10    AKAPA1(I)=AKAPA0(I)+TOLEPS
C
      fiikt0=1.D0+ACONT*(2.D0+ACONT)*AKAPA1(1)
      fiikc0=1.D0+ACONC*(2.D0+ACONC)*AKAPA1(2)
      DETTTI=1.D0-(((1.D0+ACONT-
     1    DSQRT(fiikt0))/ACONT)**CBCONT)
      DECCCI=1.D0-(((1.D0+ACONC-
     1    DSQRT(fiikc0))/ACONC)**CBCONC)     
      ft11=f0t*(SQRT(fiikt0))*((1.D0/ACONT)*(1.D0+ACONT
     1    -SQRT(fiikt0)))**(1.D0-CBCONT)
      fc11=f0c*(SQRT(fiikc0))*((1.D0/ACONC)*(1.D0+ACONC
     1    -SQRT(fiikc0)))**(1.D0-CBCONC)
C
      IF(AITR1.GT.0.D0) THEN
! ! !       IF(TAUEFGL(1).GT.TOLLL) THEN
! ! !       IF(DABS(TAUEFGL(1)).GT.DABS(TAUEFGL(3))) THEN      
	BBETT = (1.D0-AALFF)*fc11/ft11-(1.D0+AALFF)
      ELSE
	BBETT=0.d0
      ENDIF
C
      ALAMBDA2=(AALFF*AITR1+SEQV+BBETT*TAUEFGL(1)-
     1 (1.D0-AALFF)*fc11)/(9.D0*AK0*AALFFP*AALFF+
     1 DSQRT(6.D0)*GE+BBETT*(2.D0*GE*TAUDDGL(1)/SEF+
     1 3.D0*AK0*AALFFP))
C     
      SEF1=SEF-2.D0*GE*ALAMBDA2
C
      DO 70 I=1,6
  70 	 TAUDDGL1(I)=TAUDDGL(I)*SEF1/SEF
C  
      DO 80 I=1,6
	TAUEFG(I)=TAUEFGL(I)-ALAMBDA2*(2.D0*GE*TAUDDGL(I)/
     1            SEF+3.D0*AK0*AALFFP*AJED(I))
  80  CONTINUE
      AITR11=TAUEFG(1)+TAUEFG(2)+TAUEFG(3) 
      IF((DABS(TAUEFG(1))+DABS(TAUEFG(2))+
     1    DABS(TAUEFG(3))).GT.0.D0) THEN      
      Ar1=0.5D0+0.5D0*AITR11/(DABS(TAUEFG(1))+
     1   DABS(TAUEFG(2))+DABS(TAUEFG(3)))
      ELSE
      Ar1=1.D0
      ENDIF  
      DO 95 J=1,3
  95  ANABFI(J)=TAUDDGL1(J)/SEF1+AALFFP*AJED(J)      
C  
C
      DO 20 I=1,2
         DO 20 J=1,3      
  20  hTEN1(I,J)=0.D0
      hTEN1(1,1)=Ar1*ft11*(1.D0-DETTTI)/gt
      hTEN1(2,3)=-(1.D0-Ar1)*fc11*(1.D0-DECCCI)/gc
C
!       write(3,*)'ANABFI1izv',ANABFI1
!       write(3,*)'hTEN1izv',hTEN1
      CALL CLEAR(HATEN1,6)
      CALL MNOZI1(HATEN1,hTEN1,ANABFI,2,3)
C
!       write(3,*)'HATEN1izv',HATEN1
      DO 30 I=1,2
  30  QRES0(I)=ALAMBDA1*HATEN0(I)
C
      DO 40 I=1,2
  40  QRES1(I)=-AKAPA1(I)+AKAPA0(I)+ALAMBDA2*HATEN1(I)
!       write(3,*)'QRES1,QRES0izv',QRES1,QRES0
C
      DO 50 I=1,2
  50    QRESIZV(I)=(QRES1(I)-QRES22(I))/TOLEPS
C
      RETURN
      END
c     ******************************************************************   
C=======================================================================
      SUBROUTINE GLAVNTR(S,PRINC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.   PROGRAM
CS.      ZA RACUNANJE GLAVNIH NAPONA ZA 3/D ELEMENT
CE.   PROGRAM
CE.      TO CALCULATE PRINCIPAL STRESSES FOR 3/D ELEMENT
C . 
C ......................................................................
C 
      COMMON /DEFNAP/ NAPDEF
      COMMON /CDEBUG/ IDEBUG
      DIMENSION S(*),PRINC(3)
C OVU TOLERANCIJU PROVERITI 1.D-08 -STARA TOLERANCIJA
C      DATA TOL/1.D-08/
      DATA TOL/1.D-08/
C
      IF(IDEBUG.GT.0) PRINT *, ' GLAVNTR'
      IZILE=0
C
C   INARIJANTE TENZORA
C
C      ONE=0.9999999999999999D0
      ONE=1.D0
      CC=(S(1)+S(2)+S(3))/3.D0
      DO 20 I=1,3                                                              
   20 PRINC(I)=S(I)-CC                                                        
      C2=(PRINC(1)*PRINC(1)+PRINC(2)*PRINC(2)+PRINC(3)*PRINC(3))*.5D0+
     1    S(4)*S(4)+S(5)*S(5)+S(6)*S(6)             
c      C2=-PRINC(1)*PRINC(2)-PRINC(2)*PRINC(3)-PRINC(3)*PRINC(1) 
c     1   +S(4)*S(4)+S(5)*S(5)+S(6)*S(6) 
      C3=PRINC(1)*(PRINC(2)*PRINC(3)-S(5)*S(5))+
     1   S(4)*(S(5)*S(6)-S(4)*PRINC(3))+        
     1   S(6)*(S(4)*S(5)-PRINC(2)*S(6))                                          
      DUM=S(1)*(S(2)*S(3)-S(5)*S(5))+
     1    S(4)*(S(5)*S(6)-S(4)*S(3))+        
     1    S(6)*(S(4)*S(5)-S(2)*S(6))                                          
C
      SQ2 =DSQRT(2.D0)
      PI23=2.D0*DACOS(0.5D0)
C      DUM=DUM*1.D-8
C      IF(DABS(DUM).LT.1.D-9) DUM=1.D-9
      DUM=DUM*TOL
      IF(DABS(DUM).LT.0.1D0*TOL) DUM=0.1D0*TOL
      DUMM=DUM
C      WRITE(3,*) 'DUMM,C2,CC',DUMM,C2,CC
      IF(C2.LT.DABS(DUM))THEN
         DO 330 I=1,3
            PRINC(I)=CC
C ZILE 
C           ZBOG DSQRT(PRINC(I))
            IF(NAPDEF.EQ.1.AND.PRINC(I).LT.DUMM) PRINC(I)=DABS(PRINC(I))
C ZILE
  330    CONTINUE
C ZILE
          IF(IZILE.EQ.1) THEN
             PRINC(1)=PRINC(1)+DUMM
             PRINC(2)=PRINC(2)-DUMM
          ENDIF
C ZILE
         RETURN
      ENDIF
      T = DSQRT(C2/1.5D0)                                                           
      A = C3*SQ2/(T*T*T)
C
C     IF ( A .LT. -1.D0 ) A=-ONE
C     IF ( A .GT.  1.D0 ) A= ONE
      IF(DABS(A).GT.ONE) A=DSIGN(ONE,A)
      A=DACOS(A)/3.D0                                                             
      T=T*SQ2
C
      PRINC(1)=T*DCOS(A)                                                       
      PRINC(2)=T*DCOS(A-PI23)                                                
      PRINC(3)=T*DCOS(A+PI23)                                                
C
      DO 240 I=1,3                                                              
         PRINC(I)=PRINC(I)+CC                                                    
C ZILE 
C        ZBOG DSQRT(PRINC(I))
         IF(NAPDEF.EQ.1.AND.PRINC(I).LT.DUMM) PRINC(I)=DABS(PRINC(I))
C ZILE
  240 CONTINUE
C  S O R T I R A N J E
      DO 220 I=1,2                                                             
      DO 215 J=I+1,3                                                             
         IF (PRINC(I).GT.PRINC(J)) GO TO 215                                     
         CC=PRINC(I)                                                              
         PRINC(I)=PRINC(J)                                                       
         PRINC(J)=CC                                                              
  215 CONTINUE
  220 CONTINUE                                                                  
C
C ZILE
       IF(IZILE.EQ.1) THEN
          IF(DABS(PRINC(1)-PRINC(2)).LT.DUMM) THEN
             PRINC(1)=PRINC(1)+DUMM
             PRINC(2)=PRINC(2)-DUMM
          ENDIF
          IF(DABS(PRINC(2)-PRINC(3)).LT.DUMM) THEN
             PRINC(2)=PRINC(2)+DUMM
             PRINC(3)=PRINC(3)-DUMM
          ENDIF
       ENDIF
C ZILE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE GLAPR3TR(S,V,PRINC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.      RACUNANJE GLAVNIH PRAVACA NAPONA ZA 3/D ELEMENT
CE.      CALCULATE PRINCIPAL STRESS DIRECTIONS FOR 3/D ELEMENT
C . 
C ......................................................................
C 
      COMMON /CDEBUG/ IDEBUG
      DIMENSION S(*),V(3,*),VP(3),A(6),PRINC(3)
C OVU TOLERANCIJU PROVERITI 1.D-06 -STARA TOLERANCIJA
C ZA SAVIJANJE
C      DATA TOL/1.D-06/
C ZA ZATEZANJE
C      DATA TOL/1.D-08/
      DATA TOL/1.D-06/
C
      IF(IDEBUG.GT.0) PRINT *, ' GLAPR3TR'
CE  CHECK FOR ALL SAME ROOTS
      IF(DABS(PRINC(1)-PRINC(2)).LT.TOL.AND.
     1   DABS(PRINC(2)-PRINC(3)).LT.TOL) THEN
        CALL CLEAR(V,9)
        V(1,1)=1.D0
        V(2,2)=1.D0
        V(3,3)=1.D0
        RETURN
      ENDIF
CE  CHECK FOR TWO SAME ROOTS
      IP=2
      IDB=0
      IF(DABS(PRINC(1)-PRINC(2)).LT.TOL)THEN
        IP=3
        IS=1
        IT=2
        IDB=1
      ELSEIF(DABS(PRINC(2)-PRINC(3)).LT.TOL)THEN
        IP=1
        IS=2
        IT=3
        IDB=1
      ELSEIF(DABS(PRINC(1)-PRINC(3)).LT.TOL)THEN
        IP=2
        IS=3
        IT=1
        IDB=1
      ENDIF
C
    1 CONTINUE
C
        DO 10 I=1,3
   10   A(I) = S(I)-PRINC(IP)
        DO 20 I=4,6
   20   A(I) = S(I)
C
      VP(1)=A(4)*A(5)-A(2)*A(6)
      VP(2)=A(4)*A(6)-A(1)*A(5)
      VP(3)=A(1)*A(2)-A(4)*A(4)
      AI2 =DOT(VP,VP,3)
      CALL JEDNA1(V(1,IP),VP,3)
      VP(1)=A(5)*A(6)-A(3)*A(4)
      VP(2)=A(1)*A(3)-A(6)*A(6)
      VP(3)=A(4)*A(6)-A(1)*A(5)
      AI22=DOT(VP,VP,3)
      IF(AI22.GT.AI2) THEN
        AI2=AI22
        CALL JEDNA1(V(1,IP),VP,3)
      ENDIF
      VP(1)=A(2)*A(3)-A(5)*A(5)
      VP(2)=A(5)*A(6)-A(3)*A(4)
      VP(3)=A(4)*A(5)-A(2)*A(6)
      AI22=DOT(VP,VP,3)
      IF(AI22.GT.AI2) THEN
        CALL JEDNA1(V(1,IP),VP,3)
      ENDIF
      CALL JEDV(V(1,IP),V(2,IP),V(3,IP))
C
      IF(IDB.EQ.0) THEN
        IF(IP.EQ.2) THEN
          IP=3
          GO TO 1
        ELSEIF(IP.EQ.3) THEN
          IP=2
          IS=3
          IT=1
        ENDIF
      ELSEIF(IDB.EQ.1) THEN
        VP(1)=1.D0
        VP(2)=0.D0
        VP(3)=VP(2)
        AI2 =DOT(V(1,IP),VP,3)
        IF(DABS(1.D0-DABS(AI2)).LE.TOL) THEN
C        IF(DABS(1.D0-AI2).LE.TOL) THEN
          VP(2)=VP(1)
          VP(1)=VP(3)
        ENDIF
        CALL AXBV(V(1,IP),VP,V(1,IS))
        CALL JEDV(V(1,IS),V(2,IS),V(3,IS))
      ENDIF
C     GLAVNI PRAVCI PO KOLONAMA
      CALL AXBV(V(1,IP),V(1,IS),V(1,IT))
C     TRANSPONOVANJE MATRICE V DA GLAVNI PRAVCI BUDU PO VRSTAMA
      DO 200 I=1,2
      DO 200 J=I+1,3
      DUM=V(I,J)
      V(I,J)=V(J,I)
  200 V(J,I)=DUM
C     DA PROJEKCIJA PRVOG GLAVNOG PRAVCA NA X OSU BUDE POZITIVNA 
C      IF(V(1,1).LT.0.D0) THEN
C         DO 300 I=1,3
C            V(1,I)=-V(1,I)
C            V(3,I)=-V(3,I)
C  300    CONTINUE
C      ENDIF
      RETURN
      END
C=======================================================================      