C=======================================================================
C
C   TERMO-PLASTICNOST 3/D ELEMENT  -  MESOVITO OJACANJE    (22.11.1993)
C
C    SUBROUTINE D3M14
C               TI314
C               MEL3T
C               CEP314
C
C=======================================================================
      SUBROUTINE D3M14(TAU,DEF,TGT,IRAC,LPOCG,LPOC1)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
C      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
C
      LTAU  =LPOCG
      LDEFT =LTAU   + 6
      LDEFPT=LDEFT  + 6
      LALFAT=LDEFPT + 6
      LTEQT =LALFAT + 6
      LDQPT =LTEQT  + 1
      LIPL  =LDQPT  + 1
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 6
      LDEFP1=LDEFT1 + 6
      LALFA1=LDEFP1 + 6
      LTEQT1=LALFA1 + 6
      LDQPT1=LTEQT1 + 1
      LIPL1 =LDQPT1 + 1
C
      CALL TI314(PLAST(LIPL),PLAST(LDEFPT),
     1            PLAST(LALFAT),PLAST(LTEQT),PLAST(LDQPT),
     1            PLAS1(LIPL1),PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LALFA1),PLAS1(LTEQT1),PLAS1(LDQPT1),
     1            A(LFUN),MATE,TAU,DEF,TGT,A(LTEM),A(LNTA),IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI314( PL ,DEFPT,ALFAT,TEMT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEM, DEFQP,
     1                   FUN,MATE,TAU,DEF,TGT,TREF,NTFUN,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     POTPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA
C     TERMO-ELASTOPLASTICAN MATERIJAL SA MESOVITIM OJACANJEM
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /MATIZO/ E,V
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /CEPMAT/ INDCEP
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1          DEFP(*),ALFAT(*),ALFA1(*),SHET(6),SEHET(6),ETHERM(3)
      DIMENSION FUN(4,MATE*4,*),TREF(*),NTFUN(*)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
CE  INITIAL DATA
C
      IPL =PL
C
CE  E,V,ALFA,TEQY0,CY,AN,TEMP0,EM
C
      MAT4=(MAT-1)*4
      MATE4=MATE*4
      DO 6 J=1,4
        NFE=MAT4+J
        CALL BTAB(FUN,NTFUN,NFE,MATE4,TGT,NL,IND,4)
        IF (IND.EQ.2) GO TO 300
        IF (IND.EQ.1) THEN
          EVA=FUN(2,NFE,1)
        ELSE
          AMU=TGT-FUN(1,NFE,NL)
          DEN=FUN(1,NFE,NL+1)-FUN(1,NFE,NL)
          EVA=((FUN(2,NFE,NL+1)-FUN(2,NFE,NL))/DEN)*AMU+FUN(2,NFE,NL)
        END IF
        IF (J.EQ.1) E   = EVA
        IF (J.EQ.2) V = EVA
        IF (J.EQ.3) THEN
          DO 7 K=1,3
    7       ALFA(K) = EVA
          DO 71 K=4,6
   71      ALFA(K) = 0.D0
        END IF
        IF (J.EQ.4) THEN
          TEQY0=EVA
          IF (IND.EQ.1) THEN
            CY=FUN(3,NFE,1)
            AN=FUN(4,NFE,1)
          ELSE
            CY=((FUN(3,NFE,NL+1)-FUN(3,NFE,NL))/DEN)*AMU+FUN(3,NFE,NL)
            AN=((FUN(4,NFE,NL+1)-FUN(4,NFE,NL))/DEN)*AMU+FUN(4,NFE,NL)
          END IF
        END IF
    6 CONTINUE
      TEMP0 =TREF(MAT)
      EM    =FUN(3,MAT4+1,1)
C
CE    AUXILIARY CONSTANTS
C
      DVT   =2.D0/3.
      G2    =E/(1.+V)
      CM    =E/(1.-2.*V)
      ANCY  =AN*CY
      EM1   =1.-EM
      AN1   =AN-1.
C
CE    YIELD STRESS
C
      TEQY=TEQY0+CY*(EM*DEFQPT)**AN
      CALL MEL3T
      IF (IRAC.EQ.2) RETURN
C
CE    THERMAL STRAIN
C
      TEM=TGT
      IF(KOR.GT.1.AND.ITER.EQ.0) TGT=TGT-TEMT+TEMP0
      CALL STERM3(ETHERM,TGT)
      ETH = ETHERM(1)
C     ZA ITER=0, TERMICKE SILE ZA POMERANJA OD TERMICKIH DEFORMACIJA
      IF(ITER.EQ.0) GO TO 900
C
CE    TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 10 I=4,6
   10   DEFPT(I)=0.5*DEFPT(I)
C
CE    MEAN STRAIN, AND ESEKUNDUM
C
      EMT = (DEF(1)+DEF(2)+DEF(3))/3.
      DO 20 I=1,3
   20   DEFDS(I)=DEF(I)-EMT
      DO 25 I=4,6
   25   DEFDS(I)=0.5*DEF(I)
      IF (IATYP.NE.4) THEN
        DO 26 I=1,6
   26     DEFDS(I)=DEFDS(I)-DEFPT(I)
      END IF
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      DO 40 I=1,6
        TAUD(I) =G2*DEFDS(I)
        SEHET(I)=TAUD(I)-ALFAT(I)
   40   SHET(I) =SEHET(I)
      TEQ=DSQRT(1.5*TENDOT(SHET))
      TEQSEH=TEQ
C
CE   2)  CHECK FOR YIELDING
C
      IF ((TEQ-TEQY)/TEQY.LT.1.D-5) THEN
        DEFQP=DEFQPT
        DO 600 I=1,6
  600     DEFP(I)=DEFPT(I)
        GO TO 500
      END IF
C
CE   3)  SOLUTION IS ELASTO-PLASTIC.  OBTAIN ZERO OF THE ESF.
C
      PL1=1.0D0
C
CE    BISECTION
C
      DEFQP = DEFQPT
      IF (DEFQP.LT.1.D-8) DEFQP=1.D-8
      EALFA=ANCY*DEFQP**AN1
      IF (EALFA.LT.1.D-8) EALFA=1.D-8
C
      AF    = 3.D0
      KB    = 0
      KT    = 0
      DQMIN = 0.D0
      DQTOL = 1.D-10
      DEPL  = 0.D0
      FHETL = TEQY-TEQSEH
      DDEP  = -0.1*FHETL/EALFA
      DDEFQP= DDEP
C
  250   KT  = KT+1
        KB1 = KB
C
        IF(KT.GT.ITMAX) THEN
          IF (ISRPS.EQ.0) WRITE(IZLAZ,2015)
          IF (ISRPS.EQ.1) WRITE(IZLAZ,6015)
          WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
          STOP
        END IF
C
        DEFQP = DEFQPT+DDEFQP
        TEQ=TEQY0+CY*(EM*DEFQP)**AN
        EALFA=ANCY*DEFQP**AN1
        CHET=DVT*EM1*EALFA
        CHETP=AN1*CHET/DEFQP
        CG2=G2+CHET
C
        FHET=TEQ+1.5*CG2*DDEFQP-TEQSEH
C
        CALL BISECB(DDEFQP,DEPL,DEPD,DDEP,
     &              FHET,FHETL,FHETD,AF,KB,DQMIN,DQTOL)
      IF (KB.EQ.-1) GO TO 255
      IF (KB1.EQ.0) GO TO 250
      IF (DABS(DDEP).GT.EPSIL.AND.
     1    (DABS(DDEP)/(DEPL+DEPD)).GT.EPSIL) GO TO 250
C
CE   4)   PLASTIC STRAIN, BACK STRESS, DEVIATORIC STRESS
C
  255 FHETP=EALFA*EM**AN+1.5*CG2+1.5*CHETP*DDEFQP
      DLAM=1.5*DDEFQP/TEQ
      DEL=1.+CG2*DLAM
      DO 165 I=1,6
        SHET(I) =SEHET(I)/DEL
        DDEFPS  =DLAM*SHET(I)
        DEFP(I) =DEFPT(I)+DDEFPS
        ALFA1(I)=ALFAT(I)+CHET*DDEFPS
  165   TAUD(I) =SHET(I)+ALFA1(I)
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF(ISKNP.NE.2.AND.INDCEP.EQ.0)
     1     CALL CEP314(SEHET,DLAM,TEQ,DDEFQP,EM,G2,CM,FHETP,TEQSEH,
     1                 EALFA,CHET,CHETP,DEL,AN)
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
      TAUM=CM*(EMT-ETH)
      DO 200 I=1,3
  200   TAU(I)=TAUD(I)+TAUM
      DO 205 I=4,6
        TAU(I)=TAUD(I)
  205   DEFP(I)=2.*DEFP(I)
C
CE  UPDATE FROM PREVIOUS STEP
C
  900 IF(ITER.EQ.0) THEN
         ETH=ETH*(ELAST(1,1)+ELAST(1,2)+ELAST(1,3))
         TAU(1)=TAU(1)-ETH
         TAU(2)=TAU(2)-ETH
         TAU(3)=TAU(3)-ETH
      ENDIF
      DO 290 I=1,6
        DEF1(I)=DEF(I)
  290   TAU1(I)=TAU(I)
      RETURN
  300 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) NFE,TGT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) NFE,TGT
      STOP
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2015 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI314')
 2005 FORMAT(///' ARGUMENT VAN OPSEGA ZADATE KRIVE U TI314'/
     1' TEMPERATURSKA FUNKCIJA BROJ =',I5/
     2' ARGUMENT TEMPERATURA =',1PD12.4)
C-----------------------------------------------------------------------
 6015 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI314')
 6005 FORMAT(///' ARGUMENT IS OUT OF RANGE IN TI314'/
     1' TEMPERATURE FUNCTION  =',I5/
     2' ARGUMENT TEMPERATURE  =',1PD12.4)
C-----------------------------------------------------------------------
       END
C======================================================================
      SUBROUTINE CEP314(SEHET,DLAM,TEQ,DDEFQP,EM,G2,CM,FHETP,TEQSEH,
     1                 EALFA,CHET,CHETP,DEL,AN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ELAST )
CE     ELASTO-PLASTIC  CEP MATRIX
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION SEHET(*),CP(3,3)
C
      TRIPO =1.5D0
      DVA   =2.D0
      TR    =1.D0/3.
      AMULT =G2/DEL
C
      AP=TRIPO*G2/FHETP/TEQSEH
      BP=TRIPO*(1.-EALFA*DDEFQP*EM**AN/TEQ)/TEQ
      DP=AMULT*(BP-DLAM*DLAM*CHETP)/DEL*AP
      AA=AMULT*(1.+CHET*DLAM)
C
      DO 25 I=1,3
        DO 20 J=I,3
          CP(I,J)=-DP*SEHET(I)*SEHET(J)
   20   CONTINUE
        CP(I,I)=CP(I,I)+AA
   25 CONTINUE
      DO 30 I=1,3
        DO 30 J=4,6
          ELAST(I,J)=-DP*SEHET(I)*SEHET(J)
   30 CONTINUE
      DO 40 I=4,6
        DO 40 J=I,6
          ELAST(I,J)=-DP*SEHET(I)*SEHET(J)
   40 CONTINUE
      DO 45 I=4,6
        ELAST(I,I)=ELAST(I,I)+0.5*AA
   45 CONTINUE
C
      ELAST(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,3)+CM)
      ELAST(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,3)+CM)
      ELAST(1,3)=TR*(DVA*CP(1,3)-CP(1,1)-CP(1,2)+CM)
      ELAST(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,3)+CM)
      ELAST(2,3)=TR*(DVA*CP(2,3)-CP(1,2)-CP(2,2)+CM)
      ELAST(3,3)=TR*(DVA*CP(3,3)-CP(1,3)-CP(2,3)+CM)
C
      DO 50 I=1,6
        DO 50 J=I,6
          ELAST(J,I)=ELAST(I,J)
   50 CONTINUE
C
      RETURN
      END
C======================================================================
      SUBROUTINE MEL3T
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
CE     FORM ( ELAST ) MATRIX
C
      COMMON /MATIZO/ E,V
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
C
CE     NULL ( ELAST )
C
      DO 15 I=1,6
        DO 15 J=I,6
   15     ELAST(I,J)=0.0D0
C
      ELAST(1,1)=E*(1.-V)/(1.+V)/(1.-2.*V)
      ELAST(2,2)=ELAST(1,1)
      ELAST(3,3)=ELAST(1,1)
      ELAST(1,2)=ELAST(1,1)*V/(1.-V)
      ELAST(1,3)=ELAST(1,2)
      ELAST(2,3)=ELAST(1,2)
      ELAST(4,4)=ELAST(1,1)*(1.-2.*V)/(1.-V)/2.
      ELAST(5,5)=ELAST(4,4)
      ELAST(6,6)=ELAST(4,4)
      DO 50 I=1,6
        DO 50 J=I,6
   50     ELAST(J,I)=ELAST(I,J)
      RETURN
      END
C=======================================================================
C
C   PUZANJE 3/D ELEMENT                          (18.04.1994)
C
C    SUBROUTINE D3M15
C               TI315
C               CEC315
C
C=======================================================================
      SUBROUTINE D3M15(TAU,DEF,TGT,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
C
      LTAU  =LPOCG
      LDEFT =LTAU   + 6*IDVA
      LDEFCT=LDEFT  + 6*IDVA
      LOPLUT=LDEFCT + 6*IDVA
      LOMINT=LOPLUT + 6*IDVA
      LTEFT =LOMINT + 6*IDVA
      LDQCT =LTEFT  + 1*IDVA
      LIOR  =LDQCT  + 1*IDVA
      LPTAUT=LIOR   + 1*IDVA
      LA2CT =LPTAUT + 1*IDVA
      LTEQT =LA2CT  + 1*IDVA
C
      LTAU1 =LPOC1
      LDEF1 =LTAU1  + 6*IDVA
      LDEFC1=LDEF1  + 6*IDVA
      LOPLU1=LDEFC1 + 6*IDVA
      LOMIN1=LOPLU1 + 6*IDVA
      LTEF1 =LOMIN1 + 6*IDVA
      LDQC1 =LTEF1  + 1*IDVA
      LIOR1 =LDQC1  + 1*IDVA
      LPTAU1=LIOR1  + 1*IDVA
      LA2C1 =LPTAU1 + 1*IDVA
      LTEQT1=LA2C1  + 1*IDVA
C
      CALL TI315(A(LIOR),A(LDEFCT),A(LOPLUT),A(LOMINT),A(LTEFT),
     &            A(LDQCT),A(LTEQT),A(LTEQT1),
     &            A(LIOR1),A(LTAU1),A(LDEF1),A(LDEFC1),A(LOPLU1),
     &            A(LOMIN1),A(LTEF1),A(LDQC1),
     &            A(LA2CT),A(LA2C1) ,A(LPTAUT),A(LPTAU1),
     &            A(LFUN),MATE,TAU,DEF,TGT,A(LTEM),A(LNTA),IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI315(ORI ,DEFCT,OPLUT,OMINT,TEFT ,DEFQCT,TEMT ,TEM  ,
     1                 ORI1,TAU1 ,DEF1 ,DEFC ,OPLUS,OMINS ,TEF  ,DEFQC,
     &                 A2CT,A2C  ,PTAUT,PTAU ,
     1                 FUN ,MATE ,TAU  ,DEF  ,TGT  ,TREF  ,NTFUN,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     POTPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA
C     TERMO-ELASTICAN MATERIJAL SA PUZANJEM
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /SRPSKI/ ISRPS
      COMMON /MATIZO/ E,V
      COMMON /ITERBR/ ITER
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      DIMENSION DEFCT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),DEFC(*),
     1          OPLUT(*),OMINT(*),OPLUS(*),OMINS(*),ETHERM(3)
      DIMENSION FUN(2,MATE*3,*),TREF(*),NTFUN(*)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
CE  INITIAL DATA
C
      IOR =ORI
      IOR1=IOR
C
CE  E,V,ALFA,TEMP0
C
      MAT3=(MAT-1)*3
      MATE3=MATE*3
      DO 6 J=1,3
        NFE=MAT3+J
        CALL BTAB(FUN,NTFUN,NFE,MATE3,TGT,NL,IND,2)
        IF (IND.EQ.2) GO TO 300
        IF (IND.EQ.1) THEN
          EVA=FUN(2,NFE,1)
        ELSE
          AMU=TGT-FUN(1,NFE,NL)
          DEN=FUN(1,NFE,NL+1)-FUN(1,NFE,NL)
          EVA=((FUN(2,NFE,NL+1)-FUN(2,NFE,NL))/DEN)*AMU+FUN(2,NFE,NL)
        END IF
        IF (J.EQ.1) E   = EVA
        IF (J.EQ.2) V = EVA
        IF (J.EQ.3) THEN
          DO 7 K=1,3
    7       ALFA(K) = EVA
          DO 71 K=4,6
   71      ALFA(K) = 0.D0
        END IF
    6 CONTINUE
      TEMP0 =TREF(MAT)
C
CE    AUXILIARY CONSTANTS
C
      G2INV =(1.+V)/E
      CM    =E/(1.-2.*V)
      TRIPO =3.D0/2.
C
CE    ELASTIC CONSTUTIVE MATRIX
C
      CALL MEL3T
      IF (IRAC.EQ.2) RETURN
C
CE    THERMAL STRAIN
C
      CALL STERM3(ETHERM,TGT)
      ETH = ETHERM(1)
C
CE    TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 10 I=4,6
   10   DEFCT(I)=0.5*DEFCT(I)
C
CE    MEAN STRAIN, AND ESEKUNDUM
C
      EMT = (DEF(1)+DEF(2)+DEF(3))/3.
      DO 20 I=1,3
   20   DEFDS(I)=DEF(I)-EMT
      DO 25 I=4,6
   25   DEFDS(I)=0.5*DEF(I)
      IF (IATYP.NE.4) THEN
        DO 26 I=1,6
   26     DEFDS(I)=DEFDS(I)-DEFCT(I)
      END IF
      D2=1.5*TENDOT(DEFDS)
C
CE    ELASTIC DEVIATORIC STRESS SOLUTION (TAUD)
C
      DO 40 I=1,6
   40   TAUD(I)=DEFDS(I)/G2INV
      TEF=DSQRT(1.5*TENDOT(TAUD))
C
      LPU=0
      IF (TEF.LT.1.D-12) THEN
        LPU=-1
        DEFQC=DEFQCT
        DO 600 I=1,6
  600     DEFC(I)=DEFCT(I)
        GO TO 500
      END IF
C
CE   1)    OBTAIN ZERO OF THE ESF (BISECTION)
C
      AF    = 3.D0
      IB    = 0
      IT    = 0
      DQTOL = 1.D-10
      DCMIN = 1.D-10
      TEFL  = 1.D-10
      IF (ITER.EQ.0) TEFL=TEFT
      IF (TEFL.LT.1.D-10) TEFL=1.D-10
      TEFP  = TEFL
      DTEF  = 0.1*TEFT
      IF (DTEF.LT.1.D-8) DTEF=0.1*TEF
C
      PTAU=PTAUT+DT
      IF (ICLAW.EQ.1) THEN
        GAMP   = TRIPO*(ECTAU(TEFL,DEFQCT)-DEFQCT)/(DT*TEFL)
      ELSE
        GAMP   = TRIPO*ECDOT(TEFL,PTAU,TGT)/TEFL
      END IF
      GAMDT  = DT*GAMP
      AGAMA  = G2INV+GAMDT
      A2C    = A2CT
      FL     = D2-(AGAMA*TEFL)**2
C
      IF (ITER.EQ.0) GO TO 50
      TEF    = TEFL+DTEF
C
  100 IT  = IT+1
      IB1 = IB
C
      IF(IT.GT.ITMAX) THEN
        IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
        IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
        WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
        STOP
      END IF
C
CE    FINDING THE PSEUDO TIME (PTAU)
C
        IF (DT.LE.1.D-3.OR.ICLAW.EQ.1) THEN
          PTAU=PTAUT+DT
          GO TO 195
        END IF
        JB    = 0
        JT    = 0
        PTL    = 1.D-10
        DPT   = 0.1*(PTAUT+DT)
        GL    = DEFQCT+DT*ECDOT(TEF,PTL,TGT)-EC(TEF,PTL,TGT)
        PTAU  = PTL+DPT
C
  200   JT    = JT+1
        JB1   = JB
C
        IF(JT.GT.ITMAX) THEN
          IF (ISRPS.EQ.0) WRITE(IZLAZ,2010)
          IF (ISRPS.EQ.1) WRITE(IZLAZ,6010)
          WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
          STOP
        END IF
C
        G    = DEFQCT+DT*ECDOT(TEF,PTAU,TGT)-EC(TEF,PTAU,TGT)
C
        CALL BISECB(PTAU,PTL,PD,DPT,G,GL,GD,AF,JB,DCMIN,DQTOL)
        IF (JB.EQ.-1) GO TO 195
        IF (JB1.EQ.0) GO TO 200
      IF (DABS(DPT).GT.EPSIL.AND.
     1      (DABS(DPT)/(PTL+PD)).GT.EPSIL) GO TO 200
C
  195 IF (ICLAW.EQ.1) THEN
        GAMAC  = TRIPO*(ECTAU(TEF,DEFQCT)-DEFQCT)/(DT*TEF)
      ELSE
        GAMAC  = TRIPO*ECDOT(TEF,PTAU,TGT)/TEF
      END IF
      GAMDT  = DT*GAMAC
      AGAMA  = G2INV+GAMDT
C
      DDD    = TEF-TEFP
      IF (DABS(DDD).GE.1.D-7) THEN
        A2C  = (GAMAC-GAMP)/DDD
        GAMP = GAMAC
        TEFP = TEF
      END IF
C
      F = D2-(AGAMA*TEF)**2
C
      CALL BISECB(TEF,TEFL,TEFD,DTEF,F,FL,FD,AF,IB,DCMIN,DQTOL)
      IF (IB.EQ.-1) GO TO 50
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DTEF).GT.EPSIL.AND.
     1    (DABS(DTEF)/(TEFL+TEFD)).GT.EPSIL) GO TO 100
C
CE   2)   DEVIATORIC STRESS AND CREEP STRAIN
C
   50 DO 165 I=1,6
        TAUD(I) =DEFDS(I)/AGAMA
  165   DEFC(I) =DEFCT(I)+GAMDT*TAUD(I)
C
CE     E L A S T I C  -  C R E E P   M A T R I X   CEC
C
      IF (ISKNP.NE.2) CALL CEC315(AGAMA,A2C,TEF,CM)
C
CE   3)    CALCULATE STRESS
C
  500 CONTINUE
      TAUM=CM*(EMT-ETH)
      DO 400 I=1,3
  400   TAU(I)=TAUD(I)+TAUM
      DO 405 I=4,6
        TAU(I)=TAUD(I)
  405   DEFC(I)=2.*DEFC(I)
C
CE     THE MODIFIED EFFECTIVE CREEP STRAIN (DEFQC)
C
      IF (LPU.NE.-1) THEN
        DO 410 I=1,6
          OPLUS(I)=OPLUT(I)
  410     OMINS(I)=OMINT(I)
        CALL ORNL3(TAU,DEFC,DEFQC,OPLUS,OMINS,IOR1)
        ORI1=IOR1
      END IF
C
CE  UPDATE FROM PREVIOUS STEP
C
      DO 290 I=1,6
        DEF1(I)=DEF(I)
  290   TAU1(I)=TAU(I)
      RETURN
  300 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) NFE,TGT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) NFE,TGT
      STOP
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI315')
 2005 FORMAT(///' ARGUMENT VAN OPSEGA ZADATE KRIVE U TI315'/
     1' TEMPERATURSKA FUNKCIJA BROJ =',I5/
     2' ARGUMENT TEMPERATURA =',1PD12.4)
 2010 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI315 ',
     &        '( PSEUDO-TIME )')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI315')
 6005 FORMAT(///' ARGUMENT IS OUT OF RANGE IN TI315'/
     1' TEMPERATURE FUNCTION  =',I5/
     2' ARGUMENT TEMPERATURE  =',1PD12.4)
 6010 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI315 ',
     &        '( PSEUDO-TIME )')
C-----------------------------------------------------------------------
       END
C======================================================================
      SUBROUTINE CEC315(AGAMA,A2C,TEF,CM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEC ( ELAST )
CE     ELASTO-CREEP  CEC MATRIX
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      DIMENSION CP(3,3)
C
      DVA   =2.D0
      TR    =1.D0/3.
      A2CDT =A2C*DT
C
      AA=1./AGAMA
      IF (TEF.LT.1.D-8) TEF=1.D-8
      DP=3*A2CDT/(2*(AGAMA**3)*TEF*(A2CDT*TEF+AGAMA))
C
      DO 25 I=1,6
        DO 20 J=I,6
          EIJ=-DP*DEFDS(I)*DEFDS(J)
          IF (J.LT.4) THEN
            CP(I,J)=EIJ
          ELSE
            ELAST(I,J)=EIJ
          END IF
   20   CONTINUE
        IF (I.LT.4) THEN
          CP(I,I)=CP(I,I)+AA
        ELSE
          ELAST(I,I)=ELAST(I,I)+.5*AA
        END IF
   25 CONTINUE
C
      ELAST(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,3)+CM)
      ELAST(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,3)+CM)
      ELAST(1,3)=TR*(DVA*CP(1,3)-CP(1,1)-CP(1,2)+CM)
      ELAST(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,3)+CM)
      ELAST(2,3)=TR*(DVA*CP(2,3)-CP(1,2)-CP(2,2)+CM)
      ELAST(3,3)=TR*(DVA*CP(3,3)-CP(1,3)-CP(2,3)+CM)
C
      DO 50 I=1,6
        DO 50 J=I,6
          ELAST(J,I)=ELAST(I,J)
   50 CONTINUE
C
      RETURN
      END
C=======================================================================
C
C   TERMO-ELASTO-PLASTICNOST SA PUZANJEM -  3/D ELEMENT     (17.04.1994)
C               (IZOTROPAN SA MESOVITIM OJACANJEM)
C
C    SUBROUTINE D3M16
C               TI316
C               EPC316
C
C=======================================================================
      SUBROUTINE D3M16(TAU,DEF,TGT,IRAC,LPOCG,LPOC1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
C
      LFUN=MREPER(1)
      LNTA=MREPER(2)
      LTEM=MREPER(3)
      MATE=MREPER(4)
C
      LTAU  =LPOCG
      LDEFT =LTAU   + 6*IDVA
      LDEFPT=LDEFT  + 6*IDVA
      LALFAT=LDEFPT + 6*IDVA
      LDEFCT=LALFAT + 6*IDVA
      LOPLUT=LDEFCT + 6*IDVA
      LOMINT=LOPLUT + 6*IDVA
      LTEQT =LOMINT + 6*IDVA
      LTEFT =LTEQT  + 1*IDVA
      LDQPT =LTEFT  + 1*IDVA
      LDQCT =LDQPT  + 1*IDVA
      LIPL  =LDQCT  + 1*IDVA
      LIOR  =LIPL   + 1*IDVA
      LPTAUT=LIOR   + 1*IDVA
      LA2CT =LPTAUT + 1*IDVA
C
      LTAU1 =LPOC1
      LDEF1 =LTAU1  + 6*IDVA
      LDEFP1=LDEF1  + 6*IDVA
      LALFA1=LDEFP1 + 6*IDVA
      LDEFC1=LALFA1 + 6*IDVA
      LOPLU1=LDEFC1 + 6*IDVA
      LOMIN1=LOPLU1 + 6*IDVA
      LTEQ1 =LOMIN1 + 6*IDVA
      LTEF1 =LTEQ1  + 1*IDVA
      LDQP1 =LTEF1  + 1*IDVA
      LDQC1 =LDQP1  + 1*IDVA
      LIPL1 =LDQC1  + 1*IDVA
      LIOR1 =LIPL1  + 1*IDVA
      LPTAU1=LIOR1  + 1*IDVA
      LA2C1 =LPTAU1 + 1*IDVA
C
      CALL TI316(A(LIOR) ,A(LDEFCT),A(LOPLUT),A(LOMINT),A(LDQCT),
     &            A(LIPL) ,A(LDEFPT),A(LALFAT),A(LDQPT) ,A(LTEQT),
     &            A(LIOR1),A(LDEFC1),A(LOPLU1),A(LOMIN1),A(LDQC1),
     &            A(LIPL1),A(LDEFP1),A(LALFA1),A(LDQP1) ,A(LTEQ1),
     &            A(LTAU1),A(LDEF1) ,A(LTEFT) ,A(LTEF1) ,
     &            A(LA2CT),A(LA2C1) ,A(LPTAUT),A(LPTAU1),
     &            A(LFUN),MATE,TAU,DEF,TGT,A(LTEM),A(LNTA),IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TI316( ORI ,DEFCT,OPLUT,OMINT ,DEFQCT,
     &                  PL  ,DEFPT,ALFAT,DEFQPT,TEMT  ,
     &                  ORI1,DEFC ,OPLUS,OMINS ,DEFQC ,
     &                  PL1 ,DEFP ,ALFA1,DEFQP ,TEM   ,
     &                  TAU1,DEF1 ,TEFT ,TEF   ,
     &                  A2CT,A2C  ,PTAUT,PTAU  ,
     &                  FUN ,MATE ,TAU  ,DEF   ,TGT   ,TREF,NTFUN,IRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     POTPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA
C     TERMO-ELASTO-PLASTICAN MATERIJAL SA PUZANJEM
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     &                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     &                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /CREEPI/ ICELAW,ICLAW,ISG1,ITI1,ITH1,ISG2,ITI2,ITH2
      COMMON /SRPSKI/ ISRPS
      COMMON /MATIZO/ E,V
      COMMON /ITERBR/ ITER
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      DIMENSION TAU(*),DEF(*),TAU1(*),DEF1(*),ETHERM(3),
     &          DEFCT(*),DEFC(*),OPLUT(*),OMINT(*),OPLUS(*),OMINS(*),
     &          DEFPT(*),DEFP(*),ALFAT(*),ALFA1(*),SHET(6) ,GLITL(6),
     &          AHET(6)
      DIMENSION FUN(4,MATE*4,*),TREF(*),NTFUN(*)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
CE  INITIAL DATA
C
      IOR =ORI
      IOR1=IOR
      IPL =PL
C
CE  E,V,ALFA,TEQY0,CY,AN,TEMP0,EM
C
      MAT4=(MAT-1)*4
      MATE4=MATE*4
      DO 6 J=1,4
        NFE=MAT4+J
        CALL BTAB(FUN,NTFUN,NFE,MATE4,TGT,NL,IND,4)
        IF (IND.EQ.2) GO TO 300
        IF (IND.EQ.1) THEN
          EVA=FUN(2,NFE,1)
        ELSE
          AMU=TGT-FUN(1,NFE,NL)
          DEN=FUN(1,NFE,NL+1)-FUN(1,NFE,NL)
          EVA=((FUN(2,NFE,NL+1)-FUN(2,NFE,NL))/DEN)*AMU+FUN(2,NFE,NL)
        END IF
        IF (J.EQ.1) E   = EVA
        IF (J.EQ.2) V = EVA
        IF (J.EQ.3) THEN
          DO 7 K=1,3
    7       ALFA(K) = EVA
          DO 71 K=4,6
   71      ALFA(K) = 0.D0
        END IF
        IF (J.EQ.4) THEN
          TEQY0=EVA
          IF (IND.EQ.1) THEN
            CY=FUN(3,NFE,1)
            AN=FUN(4,NFE,1)
          ELSE
            CY=((FUN(3,NFE,NL+1)-FUN(3,NFE,NL))/DEN)*AMU+FUN(3,NFE,NL)
            AN=((FUN(4,NFE,NL+1)-FUN(4,NFE,NL))/DEN)*AMU+FUN(4,NFE,NL)
          END IF
        END IF
    6 CONTINUE
      TEMP0 =TREF(MAT)
      EM    =FUN(3,MAT4+1,1)
C
CE    AUXILIARY CONSTANTS
C
      G2INV =(1.+V)/E
      CM    =E/(1.-2.*V)
      ANCY  =AN*CY
      EM1   =1.-EM
      AN1   =AN-1.
      DVT   =2.D0/3.
      TRIPO =3.D0/2.
C
CE    YIELD STRESS
C
      TEQY=TEQY0+CY*(EM*DEFQPT)**AN
      CALL MEL3T
      IF (IRAC.EQ.2) RETURN
C
CE    THERMAL STRAIN
C
      CALL STERM3(ETHERM,TGT)
      ETH = ETHERM(1)
C
CE    TRANSFORM ENGENEER. SHEAR STRAIN INTO TENSORIAL
C
      DO 10 I=4,6
        DEFCT(I)=0.5*DEFCT(I)
   10   DEFPT(I)=0.5*DEFPT(I)
C
CE    MEAN STRAIN, AND ESEKUNDUM
C
      EMT = (DEF(1)+DEF(2)+DEF(3))/3.
      DO 20 I=1,3
   20   DEFDS(I)=DEF(I)-EMT
      DO 25 I=4,6
   25   DEFDS(I)=0.5*DEF(I)
      IF (IATYP.NE.4) THEN
        DO 26 I=1,6
   26     DEFDS(I)=DEFDS(I)-DEFPT(I)-DEFCT(I)
      END IF
      D0=TRIPO*TENDOT(DEFDS)
C
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD)
C
      AGAMA = G2INV
      DO 40 I=1,6
        TAUD(I) =DEFDS(I)/AGAMA
        SHET(I) =TAUD(I)-ALFAT(I)
   40   GLITL(I)=SHET(I)*AGAMA
      TEQ=DSQRT(TRIPO*TENDOT(SHET))
      TEQE=TEQ
      TEF=DSQRT(TRIPO*TENDOT(TAUD))
C
CE   2)  CHECK FOR YIELDING
C
      LPU=0
      IF ((TEQ-TEQY)/TEQY.LT.1.D-5.OR.ITER.EQ.0) THEN
        IPL1=0
        DEFQP=DEFQPT
        IF (TEF.LT.1.D-10) THEN
          LPU=-1
          DEFQC=DEFQCT
          DO 601 I=1,6
  601       DEFC(I)=DEFCT(I)
          GO TO 500
        END IF
      ELSE
        PL1 =1.0D0
        IPL1=1
      END IF
C
CE   3)  SOLUTION IS ELASTO-PLASTIC WITH CREEP.
C
CE         OBTAIN ZERO OF THE ESF (BISECTION)
C
      AF    = 3.D0
      IB    = 0
      IT    = 0
      DQTOL = 1.D-10
      DCMIN = 1.D-10
      DQMIN = 0.D0
      TEFL  = 1.D-10
      IF (ITER.EQ.0) TEFL=TEFT
      IF (TEFL.LT.1.D-10) TEFL=1.D-10
      TEFP  = TEFL
      DTEF  = 0.1*TEFT
      IF (DTEF.LT.1.D-8) DTEF=0.1*TEF
C
      PTAU=PTAUT+DT
      IF (ICLAW.EQ.1) THEN
        GAMP   = TRIPO*(ECTAU(TEFL,DEFQCT)-DEFQCT)/(DT*TEFL)
      ELSE
        GAMP   = TRIPO*ECDOT(TEFL,PTAU,TGT)/TEFL
      END IF
      GAMDT  = DT*GAMP
      AGAMA  = G2INV+GAMDT
      AGLHET = AGAMA
      D2     = D0
      A2C    = A2CT
C
      IF (ITER.EQ.0) GO TO 50
C
      FL     = D2-(AGLHET*TEFL)**2
      TEF    = TEFL+DTEF
C
  100 IT  = IT+1
      IB1 = IB
C
      IF(IT.GT.ITMAX) THEN
        IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
        IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
        WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
        STOP
      END IF
C
CE    FINDING THE PSEUDO TIME (PTAU)
C
        IF (DT.LE.1.D-3.OR.ICLAW.EQ.1) THEN
          PTAU=PTAUT+DT
          GO TO 195
        END IF
        JB    = 0
        JT    = 0
        PTL    = 1.D-10
        DPT   = 0.1*(PTAUT+DT)
        GL    = DEFQCT+DT*ECDOT(TEF,PTL,TGT)-EC(TEF,PTL,TGT)
        PTAU  = PTL+DPT
C
  200   JT    = JT+1
        JB1   = JB
C
        IF(JT.GT.ITMAX) THEN
          IF (ISRPS.EQ.0) WRITE(IZLAZ,2010)
          IF (ISRPS.EQ.1) WRITE(IZLAZ,6010)
          WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
          STOP
        END IF
C
        G    = DEFQCT+DT*ECDOT(TEF,PTAU,TGT)-EC(TEF,PTAU,TGT)
C
        CALL BISECB(PTAU,PTL,PD,DPT,G,GL,GD,AF,JB,DCMIN,DQTOL)
        IF (JB.EQ.-1) GO TO 195
        IF (JB1.EQ.0) GO TO 200
      IF (DABS(DPT).GT.EPSIL.AND.
     1      (DABS(DPT)/(PTL+PD)).GT.EPSIL) GO TO 200
C
  195 IF (ICLAW.EQ.1) THEN
        GAMAC  = TRIPO*(ECTAU(TEF,DEFQCT)-DEFQCT)/(DT*TEF)
      ELSE
        GAMAC  = TRIPO*ECDOT(TEF,PTAU,TGT)/TEF
      END IF
      GAMDT  = DT*GAMAC
      AGAMA  = G2INV+GAMDT
      AGLHET = AGAMA
C
      DDD    = TEF-TEFP
      IF (DABS(DDD).GE.1.D-7) THEN
        A2C  = (GAMAC-GAMP)/DDD
        GAMP = GAMAC
        TEFP = TEF
      END IF
C
CE    FINDING THE RADIUS OF THE YIELD SURFACE (TEQ)
C
        DLAM  = 0.D0
        DO 240 I=1,6
  240     GLITL(I) = DEFDS(I)-ALFAT(I)*AGAMA
        GNORM = DSQRT(TRIPO*TENDOT(GLITL))
        DEFQP = DEFQPT
        TEQY  = TEQY0+CY*(EM*DEFQP)**AN
        IF (DEFQP.LT.1.D-8) DEFQP=1.D-8
        EALFA = ANCY*DEFQP**AN1
        IF (EALFA.LT.1.D-8) EALFA=1.D-8
C
        KB    = 0
        KT    = 0
        DEPL  = 0.D0
        DDEP  = 0.1*(TEQE-TEQY)/EALFA
        FHETL = GNORM-AGAMA*TEQY
        IF (FHETL.LE.0.D0) THEN
          IPL1=0
          GO TO 260
        END IF
        IPL1=1
        DDEFQP= DDEP
C
  250   KT    = KT+1
        KB1   = KB
C
        IF(KT.GT.ITMAX) THEN
          IF (ISRPS.EQ.0) WRITE(IZLAZ,2015)
          IF (ISRPS.EQ.1) WRITE(IZLAZ,6015)
          WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
          STOP
        END IF
C
        DEFQP = DEFQPT+DDEFQP
        TEQY  = TEQY0+CY*(EM*DEFQP)**AN
        EALFA = ANCY*DEFQP**AN1
        CHET  = DVT*EM1*EALFA
        AGPOT = 1.+AGAMA*CHET
C
        FHET  = GNORM-AGAMA*TEQY-TRIPO*AGPOT*DDEFQP
C
        CALL BISECB(DDEFQP,DEPL,DEPD,DDEP,
     &              FHET,FHETL,FHETD,AF,KB,DQMIN,DQTOL)
        IF (KB.EQ.-1) GO TO 255
        IF (KB1.EQ.0) GO TO 250
      IF (DABS(DDEP).GT.EPSIL.AND.
     1      (DABS(DDEP)/(DEPL+DEPD)).GT.EPSIL) GO TO 250
C
  255   DLAM = TRIPO*DDEFQP/TEQY
        TEQ  = TEQY
        DO 251 I=1,6
  251     AHET(I) = ALFAT(I)+CHET*DEFDS(I)
        AGLHET = AGAMA+AGPOT*DLAM
        D2 = D0+3*DLAM*TDOTAB(DEFDS,AHET)+TRIPO*DLAM*DLAM*TENDOT(AHET)
C
  260 F = D2-(AGLHET*TEF)**2
C
      CALL BISECB(TEF,TEFL,TEFD,DTEF,F,FL,FD,AF,IB,DCMIN,DQTOL)
      IF (IB.EQ.-1) GO TO 50
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DTEF).GT.EPSIL.AND.
     1    (DABS(DTEF)/(TEFL+TEFD)).GT.EPSIL) GO TO 100
C
CE   4)   DEVIATORIC STRESS,PLASTIC STRAIN,BACK STRESS AND CREEP STRAIN
C
   50 DO 165 I=1,6
        TAUD(I) = DEFDS(I)/AGAMA
        IF ((IPL1.EQ.1).AND.(ITER.NE.0)) THEN
          SHET(I)  = GLITL(I)/AGLHET
          DDEFPS   = DLAM*SHET(I)
          DEFP(I)  = DEFPT(I)+DDEFPS
          ALFA1(I) = ALFAT(I)+CHET*DDEFPS
          TAUD(I)  = TAUD(I)-DDEFPS/AGAMA
        ELSE
          DEFP(I)  = DEFPT(I)
        END IF
        DEFC(I) = DEFCT(I)+GAMDT*TAUD(I)
  165 CONTINUE
C
CE     E L A S T I C  -  P L A S T I C  -  C R E E P   M A T R I X  CEPC
C
      IF (ISKNP.NE.2)
     &    CALL EPC316(AGAMA,A2C,TEF,CM,CHET,AN1,DEFQP,DLAM,AGPOT,IPL1,
     &                AHET,AGLHET,EM,AN,EALFA,DDEFQP,GNORM,GLITL,TEQY,
     &                ALFAT)
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
      IF (IPL1.NE.1) THEN
        DO 600 I=1,6
  600     DEFP(I)=DEFPT(I)
        TEQ=TEQE
      END IF
      TAUM = CM*(EMT-ETH)
      DO 400 I=1,3
  400   TAU(I) = TAUD(I)+TAUM
      DO 405 I=4,6
        TAU(I)  = TAUD(I)
        DEFP(I) = 2.*DEFP(I)
  405   DEFC(I) = 2.*DEFC(I)
C
CE     THE MODIFIED EFFECTIVE CREEP STRAIN (DEFQC)
C
      IF (LPU.NE.-1) THEN
        DO 410 I=1,6
          OPLUS(I) = OPLUT(I)
  410     OMINS(I) = OMINT(I)
        CALL ORNL3(TAU,DEFC,DEFQC,OPLUS,OMINS,IOR1)
        ORI1 = IOR1
      END IF
C
CE  UPDATE FROM PREVIOUS STEP
C
      DO 290 I=1,6
        DEF1(I) = DEF(I)
  290   TAU1(I) = TAU(I)
      RETURN
  300 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005) NFE,TGT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005) NFE,TGT
      STOP
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI316')
 2005 FORMAT(///' ARGUMENT VAN OPSEGA ZADATE KRIVE U TI316'/
     1' TEMPERATURSKA FUNKCIJA BROJ =',I5/
     2' ARGUMENT TEMPERATURA =',1PD12.4)
 2010 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI316 ',
     &        '( PSEUDO-TIME )')
 2015 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TI316 ',
     &        '( RADIJUS POVRSI TECENJA )')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI316')
 6005 FORMAT(///' ARGUMENT IS OUT OF RANGE IN TI316'/
     1' TEMPERATURE FUNCTION  =',I5/
     2' ARGUMENT TEMPERATURE  =',1PD12.4)
 6010 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI316 ',
     &        '( PSEUDO-TIME )')
 6015 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TI316 ',
     &        '( THE RADIUS OF YIELD SURFACE )')
C-----------------------------------------------------------------------
       END
C======================================================================
      SUBROUTINE EPC316(AGAMA,A2C,TEF,CM,CHET,AN1,DEFQP,DLAM,AGPOT,IPL1,
     &                  AHET,AGLHET,EM,AN,EALFA,DDEFQP,GNORM,GLITL,TEQ,
     &                  ALFAT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE EPC ( ELAST )
CE     ELASTO-PLASTIC-CREEP  EPC MATRIX
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      DIMENSION CP(3,3),AHET(*),EAHET(6),GLITL(*),ALFAT(*)
C
      DVA   =2.D0
      TR    =1.D0/3.
      A2CDT =A2C*DT
      IF (TEF.LT.1.D-8) TEF=1.D-8
C
      IF ((IPL1.NE.1).OR.(ITER.EQ.0)) THEN
        AA=1./AGAMA
        DP=3*A2CDT/(2*(AGAMA**3)*TEF*(A2CDT*TEF+AGAMA))
C
        DO 25 I=1,6
          DO 20 J=I,6
            EIJ=-DP*DEFDS(I)*DEFDS(J)
            IF (J.LT.4) THEN
              CP(I,J)=EIJ
            ELSE
              ELAST(I,J)=EIJ
            END IF
   20     CONTINUE
          IF (I.LT.4) THEN
            CP(I,I)=CP(I,I)+AA
          ELSE
            ELAST(I,I)=ELAST(I,I)+.5*AA
          END IF
   25   CONTINUE
      ELSE
C
        EMEAL=(EM**AN)*EALFA
        CHETP=CHET*AN1/DEFQP
        AMGAM=AGAMA*(EMEAL+1.5*CHETP*DDEFQP)+1.5*AGPOT
        ASPOT=1.5/GNORM/AMGAM
        BGPOT=A2CDT*(TEQ+1.5*CHET*DDEFQP+1.5*TDOTAB(GLITL,ALFAT)/GNORM)
     &              /AMGAM
        DLFAK=(1.5-DLAM*EMEAL)/TEQ
        ASIG =ASPOT*DLFAK
        BGAM =BGPOT*DLFAK
        CHDL1=1.+CHET*DLAM
        BIFAK=CHDL1*A2CDT
        DIFAK=CHETP*DLAM
        AA   =CHDL1/AGLHET
        ATEF2=AGLHET*TEF*TEF
        DO 100 I=1,6
  100     EAHET(I)=DEFDS(I)+DLAM*AHET(I)
C
        ASIGMA=DVA*BIFAK*ATEF2
        BSIGMA=DVA*AGPOT*ATEF2-3.*TDOTAB(EAHET,AHET)
        CSIGMA=DLAM*(DVA*AGAMA*ATEF2-3.*TDOTAB(EAHET,DEFDS))
        DSIGMA=ASIGMA-BGAM*BSIGMA-BGPOT*CHETP*CSIGMA+2*AGLHET*AGLHET*TEF
        IF (DSIGMA.LT.1.D-8) DSIGMA=1.D-8
        HJFAK =ASIG*BSIGMA+ASPOT*CHETP*CSIGMA
C
        DO 110 I=1,6
          BI=BIFAK*TAUD(I)
          CI=AGPOT*TAUD(I)-AHET(I)
          DI=DIFAK*(AGAMA*TAUD(I)-DEFDS(I))
          DO 120 J=I,6
            HJ=3*CHDL1*EAHET(J)-HJFAK*GLITL(J)
            SJOT=HJ/DSIGMA
            EIJ=(-BI*SJOT-CI*(ASIG*GLITL(J)-BGAM*SJOT)
     &           -DI*(ASPOT*GLITL(J)-BGPOT*SJOT))/AGLHET
            IF (J.LT.4) THEN
              CP(I,J)=EIJ
            ELSE
              ELAST(I,J)=EIJ
            END IF
  120     CONTINUE
          IF (I.LT.4) THEN
            CP(I,I)=CP(I,I)+AA
          ELSE
            ELAST(I,I)=ELAST(I,I)+.5*AA
          END IF
  110   CONTINUE
      END IF
C
      ELAST(1,1)=TR*(DVA*CP(1,1)-CP(1,2)-CP(1,3)+CM)
      ELAST(1,2)=TR*(DVA*CP(1,2)-CP(1,1)-CP(1,3)+CM)
      ELAST(1,3)=TR*(DVA*CP(1,3)-CP(1,1)-CP(1,2)+CM)
      ELAST(2,2)=TR*(DVA*CP(2,2)-CP(1,2)-CP(2,3)+CM)
      ELAST(2,3)=TR*(DVA*CP(2,3)-CP(1,2)-CP(2,2)+CM)
      ELAST(3,3)=TR*(DVA*CP(3,3)-CP(1,3)-CP(2,3)+CM)
C
      DO 50 I=1,6
        DO 50 J=I,6
          ELAST(J,I)=ELAST(I,J)
   50 CONTINUE
C
      RETURN
      END
