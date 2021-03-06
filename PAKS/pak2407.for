C=======================================================================
C
CS   GENERALISANI MODEL SA KAPOM  3/D ELEMENT
CE   GENERALIZED CAP MODEL 3D ELEMENT
C
CE    SUBROUTINE D2M14
C               TAUI35
C               TEQBI3
C               PRILA3
C               DEVEQ3
C
      SUBROUTINE D2M7(TAU,DEF,IRAC,LPOCG,LPOC1,IBTC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
CS     NA NIVOU INTEGRACIONE TACKE
CE     PROGRAM FOR DEFINITION OF LOCATIONS AT INTEGRATION PIONT LEVEL
      include 'paka.inc'
      
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
C
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAU(4),DEF(4)
C
      IF(IDEBUG.GT.0) PRINT *, ' D2M14'
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
      LELT=LXT + 1*IDVA
      LEMP0=LELT + 1*IDVA
      LIPL=LEMP0 + 1*IDVA
C
      LTAU1=LPOC1
      LDEFT1=LTAU1 + 4*IDVA
      LDEFP1=LDEFT1 + 4*IDVA
      LEMP1=LDEFP1 + 4*IDVA
      LXTDT=LEMP1 + 1*IDVA
      LELTDT=LXTDT + 1*IDVA
      LEMP01=LELTDT + 1*IDVA
      LIPL1=LEMP01 + 1*IDVA
C
      CALL TI207(A(LTAU),A(LDEFT),A(LDEFPP),A(LEMP),A(LXT),A(LELT),
     1       A(LEMP0),A(LIPL),
     1       A(LTAU1),A(LDEFT1),A(LDEFP1),A(LEMP1),A(LXTDT),A(LELTDT),
     1       A(LEMP01),A(LIPL1), 
     1       A(LFUN),A(LNTA),MATE,TAU,DEF,IRAC,IBTC)
C
      RETURN
      END
C
C  =====================================================================
C
      SUBROUTINE TI207(TAUT,DEFT,DEFPP,EMP,XT,ELT,EMP2,PL,
     1                  TAU1,DEF1,DEFP1,EMP1,XTDT,ELTDT,EMP0,PL1,
     1                  FUN,NTA,MATE,TAU,DEF,IRAC,IBTC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CS     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA 
CS     GENERALISANI MODEL SA KAPOM
CE     PROGRAM FOR STRESS INTEGRATION FOR GENERALIZED CAP MODEL
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
C
      COMMON /TAUD2/ TAUD(4),DEFDPR(4),DEFDS(4),DDEFP(4),
     1              DETAU(4),DDEF(4)
      COMMON/MAT2D/E,ANI,ET,TEQY0
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
      COMMON/PLASTI/LPLAST,LPLAS1,LSIGMA
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
C
      COMMON/ITERBR/ITER
C
      COMMON /CONMAT/ AE,EP,DVT
      COMMON /MATERb/ korz(100,100,3),evg(100,100,3)
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION TAUT(*),DEFT(*),DEFPP(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1          DEFP1(*)
      DIMENSION FUN(11,*),NTA(*),FUN2(30),elas(4,4)
      IF(IDEBUG.EQ.1) PRINT *, 'TI214'
C
      IF(IRAC.EQ.2) RETURN
C     INDIKATOR KONTROLNE STAMPE
      IST=0
CS     OSNOVNE KONSTANTE
CE     BASIC KONSTANTS
C
      TOLFY=1.0D-8
C
      IPL =PL
      IPL1=PL1
      DVT=2.0D0/3.0D0
      SQ2=DSQRT(2.0D0)
      DJP=DSQRT(1.5D0)
      MAXT = 100
      IPLSS=IETYP
      TOL = 1.D-9
      TOLL = 1.D-8
      TOLDP = 1.D-5
      TOLYLD = 1.D-8
      SQRT2 = SQRT(2.D0)
      AF =1.0D0
      AF3 = AF*AF*AF
      EMP2=EMP0
C
CE.     MATERIAL CONSTANTS
C
      E     = FUN(1,MAT)
      ANI   = FUN(2,MAT)
C
      A     = FUN( 3,MAT)
      B1    = -FUN( 4,MAT)
      C     = FUN( 5,MAT)
      TETA  = -FUN( 6,MAT)
      AT  =   -FUN( 7,MAT)
      AI1A0  = FUN( 8,MAT)
C
      W     =  FUN( 9,MAT)
      D     = -FUN(10,MAT)
      R     = FUN(11,MAT)
C
      B1C = B1*C
C
      izadji=0
C     pocetak slabljenja po 10%
      if(korz(7,mat,1).gt.0.and.korz(7,mat,1).le.kor) then
         e=e*(1.-float(kor-korz(7,mat,1)+1)/10.)
      endif
c     potpuno iskovanje
      if(korz(7,mat,2).gt.0.and.korz(7,mat,2).le.kor) then
         e=e/1.e9
         izadji=1
      endif
c     odumiranje i zamena materijala
      if(korz(7,mat,3).gt.0.and.korz(7,mat,3).le.kor) then
         e=evg(7,mat,1)
         ani=evg(7,mat,2)
         do i=1,4
            def(i)=def(i)-deft(i)
         enddo
         izadji=1
      endif
c
      IF(KOR.GT.1) SM0=0.0D0
      FUN2(1)=E
      FUN2(2)=ANI
      IDUM=MAT
      MAT=1
C
      CALL MEL01(FUN2)
      MAT=IDUM

      if(izadji.eq.1) then
         k4=3
         DO 80 I=1,K4
         DO 80 J=1,K4
   80       ELAS(I,J)=ELAST(I,J)
         CALL CLEAR(TAU,4) 
         CALL MNOZI1(TAU,ELAS,DEF,K4,K4)
         go to 500
      endif
C
CS     INICIJALIZACIJA OSNOVNIH VELICINA
CE     INITIAL VALUES OF VARIABLES
C
      IF(KOR.EQ.1) THEN
      EPVT = 0.D0
      EPPV = 0.D0
      AI1 = 0.D0
      AJ2D = 0.D0
      EL = 0.D0
      ELT = 0.D0
      DKV=0.0D0
C
       AI1A = AI1A0
       DELEMP=0.0D0
C
      XP = 0.D0
      AI1AS = AI1A
      ELS = ELT
      EL =0.0D0
      XT=AI1A0
      I = 0 
C
CE     BISECTION LOOP
C
  481 I = I + 1
C
CE     SOLUTION FOR EL (NEWTON ITERATIONS)
C
      F = EL-AI1A-R*(A-C*DEXP(-B1*EL)+TETA*EL) 
C
      FP = 1.- R*(C*B1*DEXP(-B1*EL)+TETA) 
C
      DEL = F/FP
      EL = EL - DEL
C
      IF  (I.GT.MAXT) THEN                                                  
          WRITE (6,1010) I,XM,XP,X,F,FM,FP
          STOP 'DP RESENJE X'                                                                  
      ENDIF                                                                     
C                                                                               
      IF (DABS(DEL).GT.TOLL) GO TO 481
C
      AI1A0=AI1A
      XT = AI1A0
      XTDT = AI1A0
      ELT = EL
      ELTDT = EL
      XSQ = AI1A*AI1A
      DAI1A = D*AI1A
      DEPV =-W*(1.-DEXP(-D*AI1A))
      EMP0=DEPV/3.
        ENDIF
C
      CM = E/(1-2.*ANI)
      G = 0.5*E/(1.+ ANI)
       AE2 = 1./G
       AE = 0.5*AE2
       DDL=1.0D-7
       AM =1.0/CM
C
c     odavde uzeto
      AE=(1.0D0+ANI)/E
      DEFPP(3)=2.0*DEFPP(3)
      DO 16 I=1,4
      TAU(I)=0.0D0
      DDEF(I)=DEF(I)-DEFPP(I)
   16 CONTINUE
C   PROVERITI
      DEFPP(3)=0.50*DEFPP(3)
C
      KS=4
      IF(IPLSS.EQ.1)KS=4
C
      DO 11 I=1,KS
      DO 11 J=1,KS
   11 TAU(I)=TAU(I)+ELAST(I,J)*DDEF(J)             
C
      IF(IRAC.EQ.2.OR.ITER.EQ.0) GO TO 500
C
      SMTE=(TAU(1)+TAU(2)+TAU(4))/3.0D0
C
      IF(3.0*SMTE.GT.AT.OR.PL1.LT.0.0D0) THEN
      DO 111 I=1,KS
  111 TAU(I)=0.0D0
      PL1=-1.0D0
       GO TO 500
      ENDIF 
C
CS     ODREDIVANJE DEVIJATORA UKUPNE DEFORMACIJE
CE    TOTAL STRAIN DEVIATOR 
C
      EMT = (DEF(1)+DEF(2)+DEF(4))/3.0D0
C
        DEFDPR(1)=DEF(1)-EMT
        DEFDPR(2)=DEF(2)-EMT
        DEFDPR(4)=DEF(4)-EMT
        DEFDPR(3)=0.5D0*DEF(3)
C
      DEFDS(1)=DEFDPR(1)-DEFPP(1)+EMP
      DEFDS(2)=DEFDPR(2)-DEFPP(2)+EMP
      DEFDS(4)=DEFDPR(4)-DEFPP(4)+EMP
      DEFDS(3)=DEFDPR(3)-DEFPP(3)
      IF(IST.EQ.1) call wrr6(defds,4,'defs')
C
C      IF(IBTC.EQ.1) 
C
CS     ODREDIVANJE POJEDINIH CLANOVA FUNKCIJE EFEKTIVNOG NAPONA
CE    EFFECTIVE STRESS FUNCTION DEFINOTION
C
      DKV = DEFDS(1)*DEFDS(1)+DEFDS(2)*DEFDS(2)+DEFDS(4)*DEFDS(4)
     1      +2.0D0*DEFDS(3)*DEFDS(3)
      DD=DSQRT(DKV)
      DSQ = DKV
C
CS     PROVERA ELASTICNOG RESENJA
CS     ODREDIVANJE NAPONA KOJI ODGOVARA ELASTICNOM RESENJU
CE     PLASTIC YIELDING CHECK
CE     ELASTIC SOLUTION
C
      DO 32 I=1,4
   32 TAUD(I)=DEFDS(I)/AE
C      IF(IBTC.EQ.1) 
C      WRITE(3,4487)(TAUD(I),I=1,4)
C
      AJ2DE=0.5*(TAUD(1)*TAUD(1)+TAUD(2)*TAUD(2)+TAUD(4)*TAUD(4))+
     1      TAUD(3)*TAUD(3)
      AJ2DQ=DSQRT(AJ2DE)
C
      EMS = EMT - EMP
      SMTE =EMS/AM
      SMTDT=SMTE
      IF(3.0*SMTE.GT.AT.OR.PL1.LT.0.0D0) THEN
      DO 117 I=1,KS
  117 TAU(I)=0.0D0
      PL1=-1.0D0
       GO TO 500
      ENDIF 
C
      AI1A=XT
      AI1 = 3.0*SMTE
      AI1L = AI1-ELT
      XML = AI1A - ELT
      RSQ = R*R
C
      FC = AI1L*AI1L + RSQ*AJ2DE - XML*XML
      FDP = AJ2DQ - (A - C*DEXP(-B1*AI1))-TETA*AI1
C
      IF((FC.LE.TOLYLD .AND. FDP.LE.TOLYLD).OR.
     * (FDP.LE.TOLYLD.AND.AI1.GT.ELT)) THEN
          IELAST = 0
          GO TO 500
      ENDIF
      DTDT = SQRT(DKV)
      DTDT2 = DTDT*SQRT2
      DSQ2 = DTDT/SQRT2
C
CE     PLASTIC DEFORMATION
C
      IELAST = 1
      IF (FDP.GT.TOLYLD .AND. (AI1-ELT).GT.TOLYLD) GO TO 100
      IF (FC.GT.TOLYLD .AND. ((AI1-ELT).GT.TOLYLD)) THEN
          IELAST = 0
          GO TO 400
      ENDIF
      IF (FC.GT.TOLYLD .AND.(AI1-ELT).LT.TOLYLD) GO TO 200
CE     DRUCKER-PRAGER YIELDING
C
  100 XP = 0.D0
      IYIELD=1
      PL1=1.0D0
      XM = 0.0D0
      FM = 0.D0
      FP = 0.0D0
      JP = 2
      DX = -0.1*EMT
      X = 0.0D0
      IB = 0
      I = 0 
      AF=5.0D0
C
CE     BISECTION LOOP
C
   10 I = I + 1
C
      SMTDT = CM*(EMS-X)
      AI1 = 3.*SMTDT
      DEX = DEXP(-B1*AI1)
      XL = -X/(B1C*DEX+TETA)
      AJ2DQ = (DTDT2 - XL)/AE2
      BCT=B1C*DEX+TETA
      DXLE=-1.0/BCT-3.0*X*B1*B1C*DEX*CM/(BCT*BCT)
      F1P=0.5*DXLE/AE-3.0*(B1C*DEX-TETA)*CM
      F = AJ2DQ - (A - C*DEX)-TETA*AI1
      DX=-F/F1P
      X=X-DX
      IF  (I.GT.MAXT) THEN                                                  
C                                                                               
CE         NUMBER OF TRIALS FOR DEP EXCEEDS MAXIT - STOP                      
          WRITE (6,*) 'NLM,MAT',NLM,MAT
          WRITE (6,1010) I,XM,XP,X,F,FM,FP
          STOP 'DP BISEKCIJE'                                                                 
      ENDIF                                                                     
1010  FORMAT (' NO SOLUTION IN BISECTING FOR F=0. '/
     1 ' IBIS =',I5/
     3 ' X-MINUS =',D14.6/' X-PLUS=',D14.6/' XL=',D14.6/
     4        ' F =',D14.6/' F-MINUS=',D14.6/' F-PLUS=',D14.6)
C                                                                               
      IF (DABS(DX).GT.TOLL) GO TO 10
C
      XTDT = XT
      DELEMP = X
      EMP1 = EMP + DELEMP
      AI1A = -DLOG(1.+3.*EMP1/W)/D + AI1A0
C
      XP = 0.D0
      AI1AS = AI1A
      ELS = ELT
      EL =0.0D0
      DX = 0.02*AI1A
      X = DX
      IB = 0
      I = 0 
      IF(KOR.EQ.1) THEN
       EPVT=3.*EMP0
      ENDIF 
C
CE     BISECTION LOOP
C
  281 I = I + 1
C
CE     SOLUTION FOR EL (NEWTON ITERATIONS)
      RP=0.0D0
C
      F = EL-AI1A-R*(A-C*DEXP(-B1*EL)+TETA*EL) 
C
      FP = 1.- R*(C*B1*DEXP(-B1*EL)+TETA) 
C
      DEL = F/FP
      EL = EL - DEL
C
      IF  (I.GT.MAXT) THEN                                                  
          WRITE (6,1010) I,XM,XP,X,F,FM,FP
          STOP 'DP RESENJE X'                                                                  
      ENDIF                                                                     
C                                                                               
      IF (DABS(DEL).GT.TOLL) GO TO 281
C                                                                               
      ELTDT = EL
      XTDT = AI1A
      COEF = 0.5*XL/AJ2DQ
      DEN = AE + COEF
      DO 615 I=1,4
  615 TAUD(I)=DEFDS(I)/DEN
      IF(IST.EQ.1) call wrr6(taud,4,'taud')
      DO 710 I=1,4
      DEFP1(I)=DEFPP(I)+COEF*TAUD(I)
  710 CONTINUE
      IF(IST.EQ.1) call wrr6(DEFP1,4,'DEFP')
      GO TO 400
C
CE     CAP YIELDING
C
  200 XP = ELT
      XM = 0.0D0
      FM = 0.D0
      FP = FC
      JP = 2
      PL1=2.0D0
      IYIELD=2
      IMAXE=0
      AI1AS = AI1A
      ELS = ELT
      EL =ELT
      DX = 0.02*AI1A
       DEL = 0.5*AI1A0
      X = DX
      IB = 0
      I = 0 
      EPVT=3.*EMP
C
CE     BISECTION LOOP
C
  210 I = I + 1
C
CE     SOLUTION FOR EL (NEWTON ITERATIONS)
C
      RP=0.0D0
      IF(IMAXE.EQ.0) THEN
       AJ2L = A-C*DEXP(-B1*EL)+TETA*EL
       AI1A = EL - R*AJ2L
      ENDIF
C
      XSQ = AI1A*AI1A
      AI1AP=AI1A-XT
      XSQ = AI1AP*AI1AP
      DAI1AP = D*AI1AP
      DEPV =-W*(DEXP(-D*XT)-DEXP(-D*AI1A))
C
      EEMP1=EPVT+DEPV
      DDEPV =-W*D*DEXP(-D*AI1A)
C
      DELEMP = DEPV/3.
      SMTDT = CM*(EMS-DELEMP)
      AI1 = 3.*SMTDT
      XL = DELEMP/(2.*(AI1-EL))
      RSQ = R*R
      DEN = AE + RSQ*XL
      AJ2D = 0.5*DSQ/(DEN*DEN)
      AJ2DQ = DSQRT(AJ2D)
      AI1L = AI1 - EL
      XML = AI1A - EL
C
      X =AI1A
      DXL = 1.-R*(B1C*DEXP(-B1*EL)+TETA)-RP*(A-C*DEXP(-B1*EL)+TETA*EL)        
      DXDEPV = -W*D*DEXP(-D*AI1A)
CC
      DDEPV = W*D*DEXP(-D*AI1A)*(D*(AI1A-XT)-1.)
C
      DXDEPM = DXDEPV*DXL/3.
C
      DAI1L = -3.*CM*DXDEPM
C
      DJ2DL = - RSQ*DSQ/(DEN*DEN*DEN)
C
      DLCL = 0.5*(DXDEPM*AI1L - DELEMP*(DAI1L-1.))/(AI1L*AI1L)
C
      DJ2L = DJ2DL*DLCL
      F = AI1L*AI1L + RSQ*AJ2D - XML*XML
C
      FP = 2.*AI1L*(DAI1L-1.) +2.*R*AJ2D*RP + RSQ*DJ2L-
     1      2.*XML*(DXL - 1.)
C
      DEL = F/FP
C
      EL = EL - DEL
C
C     CALL BISECD (EL,XM,XP,DEL,F,FM,FP,AF,JP,IB)
C
      IF  (I.GT.MAXT) THEN                                                  
C                                                                               
CE         NUMBER OF TRIALS FOR DEP EXCEEDS MAXIT - STOP                      
CC                                                                               
          WRITE (6,1010) I,F,FP
          STOP 'KAPA'                                                                  
      ENDIF                                                                     
C                                                                               
C                                                                               
      IF (DABS(DEL).GT.TOLL) GO TO 210
C
      ELTDT = EL
      XTDT = AI1A
      EMP1 = EMP + DELEMP
      DO 625 I=1,4
  625 TAUD(I)=DEFDS(I)/DEN
      COEF = XL*RSQ
      DO 720 I=1,4
      DDEFP(I)=COEF*TAUD(I)
  720 CONTINUE
      DO 722 I=1,4
      DEFP1(I)=DEFPP(I)+DDEFP(I)
  722 CONTINUE
      GO TO 400
C
C     VERTEX YIELDING
C
  
C 300 XP = 0.D0
      XP = 0.D0
      XM = 0.0D0
      FM = 0.D0
      FP = FDP
      JP = 2
      IYIELD=3
      AI1AS = AI1A
      ELS = ELT
      DX = 0.01*AI1A
      X = DX
      IB = 0
      I = 0 
      EPVT=3.*EMP
      IF(KOR.EQ.1) EPVT=3.*EMP0
C
C     BISECTION LOOP
C
  310 I = I + 1
C
C     SOLUTION FOR EL (NEWTON ITERATIONS)
C
      EL = ELS + X
      AI1A = EL - R*(A-C*DEXP(-B1*EL)+TETA*EL)
      DEPV =-W*D*DEXP(-D*AI1A)*(AI1A-XT)
      DELEMP = DEPV/3.
      SMTDT = CM*(EMS-DELEMP)
      AI1 = 3.*SMTDT
      XML = AI1A - EL
      XL = -R*DSQ2/XML - AE
      F = 3.0*CM*(EMS-DELEMP)-EL
C
      CALL BISECD (X,XM,XP,DX,F,FM,FP,AF,JP,IB)
C
      IF  (I.GT.MAXT) THEN                                                  
C                                                                               
CE         NUMBER OF TRIALS FOR DEP EXCEEDS MAXIT - STOP                      
C                                                                               
          WRITE (6,1010) I,XM,XP,X,F,FM,FP
          STOP                                                                  
      ENDIF                                                                     
C                                                                               
      IF (IB.EQ.0) GO TO 310                                                   
C                                                                               
      IF (ABS(DX).GT.TOLL) GO TO 310
C
      ELTDT = EL
      XTDT = AI1A
      DEN = AE + XL
      EMP1 = EMP + DELEMP
      IF(KOR.EQ.1) EMP1=EMP0+DELEMP
      DO 635 I=1,4
  635 TAUD(I)=DEFDS(I)/DEN
      DO 730 I=1,4
      DEFP1(I)=DEFPP(I)+XL*TAUD(I)
  730 CONTINUE
C
      GO TO 400
C
C     UPDATES FOR NEXT STEP
C
  400 CONTINUE
      DEFP1(1)=DEFP1(1)+DELEMP
      DEFP1(2)=DEFP1(2)+DELEMP
      DEFP1(4)=DEFP1(4)+DELEMP
C
C     ODREDIVANJE PLASTICNE DEFORMACIJE
CE     STRESS CALCULATION
C
      DO 260 I=1,4
      TAU(I)=TAUD(I)+SMTDT
       IF(I.EQ.3) THEN
        TAU(I)=TAUD(I)
       ENDIF
  260 CONTINUE
C
C    ODREDIVANJE ELASTOPLASTICNE MATRICE
C
C      CALL MEL214(TAUD,DEFDS,R,RP,XL,CM,AJ2D,DDEPV,AI1L,
C     1           FP,XML,DEN,TETA,B1,C,SMTDT,DELEMP,IYIELD,IBTC)
C
C      ENDIF
C
C     KORIGOVANJE VELICINA IZ PRETHODNOG KORAKA KAD JE POSTIGNUTA
C     KONVERGENCIJA
CE    CORECTION OF VALUES FROM PREVIOUS STEP WHEN CONVERGENCE IS
CE    REATCHED
  500 CONTINUE
      DO 290 I=1,4
      DEF1(I)=DEF(I)
  290 TAU1(I)=TAU(I)
      RETURN
      END
C======================================================================
C
      SUBROUTINE MEL214(TAUD,DEFDS,R,RP,XL,CM,AJ2D,DDEPV,AI1L,
     1                 FP,XML,DEN,TETA,BB1,C1,SMTDT,DELEMP,IYIELD,IBTC)
CC      SUBROUTINE MEL39(TAU,DEFDS,SMTDT,P0TDT,DLAM,EAL,AMM,AER,AM,
CC     1           RAZ,DKV,DELEMP,ICSS)
C
C----------------------------------------------------------------------
CS    PROGRAM ZA FORMIRANJE MATRICE C<E> ILI C<EP> ZA CAM-CLAY MODEL
CE    PROGRAM FOR CONSTITUTIVE MATRIX C<E> OR C<EP>
C----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /CDEBUG/ IDEBUG
C
C
C
C
      DIMENSION C(6,6),EPSD(6),TAUD(6),DEFDS(6)
      IF(IDEBUG.EQ.1) PRINT *, 'MEL09 '
C
C
      DVT=2.0/3.0
      TR=1.0/3.0
      SQ2=DSQRT(2.0D0)
C
      DO 1 I=1,3
      EPSD(I)=DEFDS(I)
    1 EPSD(I+3)=2.0*DEFDS(I+3)
C
      D2=EPSD(1)**2+EPSD(2)**2+EPSD(3)**2+
     *(EPSD(4)**2+EPSD(5)**2+EPSD(6)**2)/2.
C
       DD=DSQRT(D2)
       AA=DEN
       RSQ=R*R
       IF(IYIELD.EQ.1) THEN
        XC=XL*CM
        B1C=BB1*CM
        B1CE=BB1*C1*DEXP(-3.0*BB1*SMTDT)
        B1= 1.0/(DELEMP+3.0*TETA)-3.0*B1C
        B2= 0.5*SQ2*CM/DD         
        B3= 0.5*B1C+3.0*CM*(B1CE+TETA)
        B4= (0.5*SQ2*XC/(AJ2D*DD) + (1.0+0.5*XC/AJ2D)*B1*B2/B3)/AJ2D
        B5= 3.0*CM*(B1CE+TETA-0.5*B1C*XL)
        B6=0.5*(1.0+0.5*XC/AJ2D)*(B1*B5+3.0*B1C*XL)/(AJ2D*DEN) 
        BP=B4
        BP6=B6
        B711=B2*CM/B3
        B715=CM*(1.0-B5)
       ENDIF
       IF(IYIELD.EQ.2) THEN
        B7 = FP*DDEPV/3.0D0
        B8 = 0.5*(B7*(1.0D0+3.0D0*XL*CM)+XL)/XML
        B9 = AA/DEN
        B10= -2.0*R*AJ2D*(R*B8+2.0*XL*RP)/DEN
        B12=-2.0*AI1L*(3.0*CM*B7+1.0)+2.0*R*AJ2D*RP-RSQ*B10-
     1      2.0*XML*(FP-1.0)
        B11=-RSQ*B9/B12
        B13=(RSQ*B8+2.0*R*XL*RP)*B11
        B14=1.5*CM*XL/AI1L
        B15=-6.0*AI1L*CM/B12
        B16=R*(R*(B8*B15-B14)+2.0*XL*RP)/DEN
        BP= B13
        BP6=B16
        B711=B7*B11*CM
        B715=CM*(1.0-B7*B15)
       ENDIF
CC
      DO 4 I=1,4
      DO 4 J=1,4
      ELAST(I,J)=0.0D0
    4 C(I,J)=0.0D0
C
C   MATRICA __Cik' NADVUCENO
C
      DO 5 I=1,4
      DO 5 J=1,4
      C(I,J)=-BP*TAUD(I)*EPSD(J)
      IF(I.EQ.J)C(I,J)=C(I,J)+AA
    5 CONTINUE
C
C    MATRICA C'ij
      DO 6 I=1,3
      DO 6 J=1,3
      DO 16 K=1,3
       IF(J.EQ.K) THEN
       TT=DVT
       ELSE
       TT=-TR
       ENDIF
      ELAST(I,J)=ELAST(I,J)+C(I,K)*TT
   16 ELAST(I+3,J)=ELAST(I+3,J)+C(I+3,K)*TT
    6 CONTINUE
C
      DO 40 I=1,3
      DO 40 J=1,3
      ELAST(I,J)=ELAST(I,J)-BP6*TAUD(I)
   40 ELAST(I+3,J)=ELAST(I+3,J)-B16*TAUD(I+3)
      DO 7 I=1,4
      DO 7 J=4,4
      ELAST(I,J)=C(I,J)/2.0
    7 CONTINUE
C
      DO 8 I=1,3
      DO 8 J=1,3
      DO 26 K=1,3
       IF(J.EQ.K) THEN
       TT=DVT
       ELSE
       TT=-TR
       ENDIF
   26  ELAST(I,J)=ELAST(I,J)-B711*EPSD(K)*TT
       ELAST(I,J+3)=ELAST(I,J+3)-B711*EPSD(J+3)*0.5
       ELAST(I,J)=ELAST(I,J)+B715
    8 CONTINUE
C
C
      DO 30 I=1,4
      DO 30 J=I,4
30    ELAST(I,J)=(ELAST(I,J)+ELAST(J,I))/2.
C
      DO 31 I=1,4
      DO 31 J=1,I
31    ELAST(I,J)=ELAST(J,I)
      RETURN
      END
C
C   =================================================================
C
      SUBROUTINE BISECD (X,XM,XP,DX,F,FM,FP,AF,JP,IB)
C
CE    SUBROUTINE TO SOLVE A NONLINEAR EQUATION BY BISECTION
CE
CE    NOTATION :
CE
CE    X  = X-ARGUMENT (POSITIVE)
CE    XM = X-MINUS
CE    XP = X-PLUS
CE    F  = FUNCTION
CE    FM = F-MINUS
CE    FP = F-PLUS
CE    JP = INDICATOR FOR SIGN OF FIRST DERIVATIVE
CE         1 = FIRST DERIVATIVE .GT.0
CE         2 = FIRST DERIVATIVE .LT.0
CE    DX = ICREMENT OF X
CE    IB = BISECTION FLAG
CE         0 = SEARCH FOR XM AND XP
CE         1 = XM AND XP ARE OBTAINED (SEARCH FOR SOLUTION BY BISECTION)
CE    AF = ACCELERATION FACTOR FOR DX (EQ.3.D0)
CE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, 'BISECD'
      IF (IB.EQ.1) GO TO 100
C
      DX = AF*DX
      IF (JP.EQ.2) GO TO 10
      IF  (F.GE.0.D0) THEN
          FP = F
          XP = X
          X = 0.5*(XP+XM)
          IB = 1
          RETURN
      ELSE
          XM = X
          X = X + DX
          FM = F
      ENDIF
      RETURN
C
   10 IF  (F.LE.0.D0) THEN
          FM = F
          XM = X
          X = 0.5*(XM+XP)
          IB = 1
          RETURN
      ELSE
          XP = X
          FP = F
          X = X + DX
      ENDIF
      RETURN
C
C
C     BISECTION
C
  100 IF (JP.EQ.2) GO TO 150
C
      IF (F.GE.0.D0) GO TO 120
      XM1 = X
      X = X + F/(FM-F)*(X-XM)
      IF (X.GE.XP) X = 0.5*(XM1+XP)
      XM = XM1
      FM = F
      DX = X - XM
      RETURN
  120 XP1 = X
      X = X - F/(FP-F)*(XP-X)
      IF (X.LE.XM) X = 0.5*(XM+XP1)
      XP = XP1
      FP = F
      DX = XP - X
      RETURN
C
  150 IF (F.GE.0.D0) GO TO 170
      XM1 = X
      X = 0.5*(XM1+XP)
      XM = XM1
      FM = F
      DX = XM - X
      RETURN
  170 XP1 = X
      X = 0.5*(XM+XP1)
      XP = XP1
      FP = F
      DX = X - XP1
      RETURN
C
      END


