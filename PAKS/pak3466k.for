C=======================================================================
C
C     PLASTICNOST 3/D ELEMENT  -  MESOVITO OJACANJE    (08.08.1992)
C
C=======================================================================
C
C     ELASTO PLASTICNA ANALIZA NA OSNOVU ALGORITMA
C     IZ TABELE (4.3.1) NELINEARNA ANALIZA KONSTRUKCIJA M. ZIVKOVIC 2006
C
C=======================================================================
      SUBROUTINE D3M66D(TAU,DEF,IRAC,LPOCG,LPOC1)
      USE PLAST3D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PROGRAM ZA ODREDIVANJE LOKACIJA VELICINA KOJE SE CUVAJU
C     NA NIVOU INTEGRACIONE TACKE
C
      include 'paka.inc'
      
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /REPERM/ MREPER(4)
      COMMON /DUPLAP/ IDVA
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR      
C
      LFUN=MREPER(1)
      MATE=MREPER(4)
C
      LTAU  =LPOCG
      LDEFT =LTAU   + 6
      LDEFPT=LDEFT  + 6
      LALFAT=LDEFPT + 6
      LTEQT =LALFAT + 6
      LDQPT =LTEQT  + 1
      LIPL  =LDQPT  + 1
      LAGE  =LIPL   + 1
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 6
      LDEFP1=LDEFT1 + 6
      LALFA1=LDEFP1 + 6
      LTEQT1=LALFA1 + 6
      LDQPT1=LTEQT1 + 1
      LIPL1 =LDQPT1 + 1
      LAGE1 =LIPL1  + 1
C
      CALL TAUD366(PLAST(LIPL),PLAST(LDEFPT),
     1            PLAST(LALFAT),PLAST(LTEQT),PLAST(LDQPT),
     1         PLAS1(LIPL1),PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LALFA1),PLAS1(LTEQT1),PLAS1(LDQPT1),
     1            A(LFUN),MATE,TAU,DEF,IRAC,PLAST(LAGE),PLAS1(LAGE1),
     1            KOR)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAUD366( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,
     1                   FUN,MATE,TAU,DEF,IRAC,AGET,AGETT,KORAK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     PODPROGRAM ZA INTEGRACIJU KONSTITUTIVNIH RELACIJA ZA ELASTOPLASTIC
C     ELASTOPLASTICAN MATERIJAL SA IZOTROPNIM OJACANJEM
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TAUD3/ TAUD(6),DEFDPR(6),DEFDS(6),DDEFP(6),
     1                DETAU(6),DDEF(6)
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ITERBR/ ITER
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /VELIKE/ LCOR0,LGM0,JG,NGR,NGS,NGT,NGS4
      COMMON /CEPMAT/ INDCEP
      COMMON /VELIKD/ DETG,QP(3,3),IGLPR
      COMMON /POCNAP/ IPOCET
      COMMON /PRINCI/ PRINC(3)
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /GRADIJ/ GRAD(3,3),GRAE(3,3),GRAP(3,3)
      DIMENSION DEFPT(*),TAU(*),DEF(*),TAU1(*),DEF1(*),
     1    DEFP(*),ALFAT(*),ALFA1(*),SHET(6),GLITL(6),SHETE(6),POM(3,3)
      DIMENSION FUN(2,MATE,*)
      DIMENSION VDEF(3,3),TAU0(6),DDEFPS(6),CT(3,3,3,3),TSG(6,6)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
CE  INITIAL DATA
C
      IF(IRAC.EQ.2) RETURN
      IPL =PL
C
CE  E,ANI,TEQY0,CY,AN,EM
C
      E     =FUN(1,MAT,1)
      ANI   =FUN(2,MAT,1)
      TEQY0 =FUN(1,MAT,2)
      CY    =FUN(2,MAT,2)
      AN    =FUN(1,MAT,3)
      EM    =FUN(2,MAT,3)
      HS    =FUN(1,MAT,4)
C
      AEL=FUN(1,MAT,5)
      AELL=FUN(2,MAT,5)
      AGEC=FUN(1,MAT,6)     
      AETA=FUN(2,MAT,6)
      EPECRIT=FUN(1,MAT,7)   
C
      ISIMO=0
      IF(DABS(HS).GT.1.D-10) ISIMO=1
C     PRIVREMENO
C     INDIKATOR ZA NAPONE (0-U DEKARTOVOM SISTEMU, 1-U GLAVNIM PRAVCIMA)
      NAPGL=0
C     INDIKATOR ZA TRANSFORMACIJU NA KOSIJEVE NAPONE (0-NE,1-DA)
      NAPKO=1
C     INDIKATOR ZA POCETNE NAPONE ZA SUPERGREDNI ELEMENT (0-NE,1-DA)
      IPOCET=0
      IF(IPOCET.EQ.1) THEN
         TEQY0=TEQY0+TAU(1)
         TAU(1)=0.D0
      ENDIF
C     KORISCENJE AJOT 
      AJOT=1.D0

C
CS    MODUL SMICANJA (4.3.11) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      GE=E/2.D0/(1.D0+ANI)
C
CS    ZAPREMINSKI MODUL (4.3.8) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      DVT   =2.D0/3.D0
      C1    =(2.D0-ANI)/3.D0/(1.D0-ANI)
      CNI   =1.D0-C1
      CYAN23=DVT*CY*AN
      EM1   =1.D0-EM
      AN1   =AN-1.D0
      AN2   =AN-2.D0
      AE    =(1.D0+ANI)/E
      AEINV =1.D0/AE
      CM    =E/(1.D0-2.D0*ANI)
C
CE    YIELD STRESS
CS    NAPON TECENJA FORMULA (4.3.31) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      IF(ISIMO.EQ.0)THEN
         TEQY =TEQY0+CY*(EM*DEFQPT)**AN
      ELSE
         TEQY =TEQY0+CY*(1.D0-DEXP(-AN*EM*DEFQPT))+HS*EM*DEFQPT
      ENDIF
C
CE    ELASTIC CONSTITUTIVE MATRIX
C

      CALL MEL36(FUN,MATE)
      IF(IRAC.EQ.2) RETURN
C
C... TRANSFORM ENGENEER. PLASTIC SHEAR STRAIN INTO TENSORIAL
C
      DO 10 I=4,6
   10 DEFPT(I)=.5D0*DEFPT(I)
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM, GLITL
CS    SREDNJA UKUPNA DEFORMACIJA (4.3.7) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC    
C
      EMT = (DEF(1)+DEF(2)+DEF(3))/3.D0
      IF(IATYP.EQ.5) EMT=0.D0
      IF(IATYP.LT.4)THEN
C
CS    PROBNE ELASTICNE DEVIJATORSKE DEVIJACIJE (4.3.41) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
         DO 20 I=1,3
   20    DEFDS(I)=DEF(I)-EMT-DEFPT(I)
         DO 25 I=4,6
   25    DEFDS(I)=.5D0*DEF(I)-DEFPT(I)
      ELSEIF(IATYP.GE.4) THEN
         DO 21 I=1,3
   21    DEFDS(I)=DEF(I)-EMT
         DO 26 I=4,6
   26    DEFDS(I)=.5D0*DEF(I)
      ENDIF
C
CE    ELASTIC DEVIATORIC STRESS RADIUS  (SHET)
CS    ELASTICNI RADIJUS NAPON (4.3.51) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      IF(AGET.LT.1.d-8)THEN
      AGE=1.d0
      ELSE
      AGE=AGET
      ENDIF
      p=0
  222 CONTINUE
      p=p+1
C      call wrr6(ALFAT,6,'ALFT')
      CALL CLEAR(TAUD,6)
      DO 41 I=1,6
        TAUD(I) =AGE*DEFDS(I)*(2.D0*GE)/AJOT       
   41   SHETE(I) =TAUD(I)-ALFAT(I)
        write(3,*)'TAUDposle',TAUD    
C
CE    EFFECTIVE ELASTIC STRESS RADIUS 
CS    EFEKTIVNI ELASTICNI RADIJUS NAPON (4.3.50) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      TEQE=DSQRT(1.5D0*TENDOT(SHETE))
C
CE    CHECK FOR YIELDING
CS    PROVERA TECENJA (4.3.55) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      IF((TEQE-TEQY)/TEQY.LT.1.D-5)THEN
CE       SOLUTION IS ELASTO
         TEQ  =TEQY
         DEFQP=DEFQPT
         DO 600 I=1,6
  600    DEFP(I)=DEFPT(I)
         CALL CLEAR(DDEFPS,6)
         GO TO 500
      ENDIF
C
CE    SOLUTION IS ELASTO-PLASTIC.  OBTAIN ZERO OF THE ESF.
C
      PL1=1.D0
C
CE    ITERATIONS
C     
      DEFQP=DEFQPT
      IF(DEFQP.LE.1.D-4) DEFQP=1.D-4
      IF(ISIMO.EQ.0)THEN
         EP2=AN*CY*DEFQP**(AN-1.D0)      
      ELSE
         EP2=AN*CY*DEXP(-AN*DEFQP)+HS
      ENDIF
C
      IB = 0
      IT = 0
      AF = 5.D0
C     DEFINISANJE INKREMENTA (DDD)
      DDD= .1D0*(TEQE-TEQY)/EP2
      GG =(1.D0/2.D0/GE)*(TEQY0-DSQRT(1.5D0*TENDOT(SHETE)))
C       
CS    PROVERA TECENJA U ELASTICNOJ OBLASTI NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      IF(ISIMO.EQ.0)THEN
         AA =(1.D0/2.D0/GE)*CY*(EM**AN)
         FP=TEQE-TEQY
      ELSE
         FP =TEQE-TEQY
      ENDIF
      TAUY  = TEQY
      FM    = 0.D0
      DEPBM = 0.D0
      DEPBP = 0.D0
      DDEFQP= DDD
      ABB=1.5D0/AJOT
  100 IT=IT+1
      IB1 = IB
C
      IF(IT.GT.ITMAX) THEN
         IF (ISRPS.EQ.0) WRITE(IZLAZ,2000)
         IF (ISRPS.EQ.1) WRITE(IZLAZ,6000)
         WRITE(IZLAZ,2001)NLM,NGR,NGS,NGT
         STOP
      ENDIF
C
CS    EFEKTIVNA PLASTICNA DEFORMACIJA (4.3.57) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      DEFQP=DEFQPT+DDEFQP
C
CS    MODUL C^ (4.3.46) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      IF(ISIMO.EQ.0)THEN
         CC=(2.D0/3.D0)*CY*AN*DEFQP**(AN-1.D0)
         CHET =(1.D0-EM)*CC
         BB =ABB+1.5D0*CHET*(1.D0/2.D0/GE)
         FB=-AA*DEFQP**AN-BB*DDEFQP-GG
      ELSE
         CC=(2.D0/3.D0)*CY*AN*DEXP(-AN*DEFQP)+(2.D0/3.D0)*HS
         CHET =(1.D0-EM)*CC
         BB =1.5D0+1.5D0*CHET*(1.D0/2.D0/GE)
         FB =-(1.D0/2.D0/GE)*(CY*(1.D0-DEXP(-AN*EM*DEFQP))
     1	   -HS*EM*DEFQP)-BB*DDEFQP-GG
      ENDIF
C
      CALL BISEC (DDEFQP,DEPBM,DEPBP,DDD,FB,FM,FP,AF,IB)
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DDD).GT.EPSIL.AND.
     1    (DABS(DDD)/(DEPBM+DEPBP)).GT.EPSIL) GO TO 100
C
CE      ...   ( DEVIATORIC STRESS )
C
CS    NAPON TECENJA (4.3.38) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      IF(ISIMO.EQ.0)THEN
         TEQ =TEQY0+CY*(EM*DEFQP)**AN
      ELSE
         TEQ =TEQY0+CY*(1.D0-DEXP(-AN*EM*DEFQP))+HS*EM*DEFQP
      ENDIF
C
CS    POZITIVAN SKALAR LAMBDA (4.3.45) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      DLAM=1.5D0*DDEFQP/TEQ

CS    RADIJUS NAPONI (4.3.44) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C      
      DO 165 I=1,6
         SHET(I)=SHETE(I)/(1.+(2.*GE+CHET)*DLAM)
C
CS    DEVIJATORSKI NAPONI (4.3.43) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
  165 TAUD(I)=ALFAT(I)+(CHET*DLAM+1.D0)*SHET(I)
C
CE    DETERMINE SOLUTION 
C
C
CE    PLASTIC STRAIN ,  BACK STRESS 
C
      DO 170 I=1,6
         DDEFPS(I)=DLAM*SHET(I)
C
CS    PLASTICNE DEFORMACIJE (4.3.39) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
         DEFP(I)  =DEFPT(I)+DDEFPS(I)
C
CS    POLOZAJNI NAPONI NA KRAJU KORAKA (4.3.43) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
         ALFA1(I) =ALFAT(I)+CHET*DDEFPS(I)
  170 CONTINUE
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
C      call wrr6(GLITL,6,'GLID')
      IF(ISKNP.NE.2.AND.INDCEP.EQ.0)
     1  CALL CEP366K(SHETE,DLAM,TEQ,DEFQP,DDEFQP,EM,GE,CM,CHET,CY,AN,
     &             HS,ISIMO)
C
CE    CALCULATE STRESS
C
  500 CONTINUE
C
C
CS    SREDNJI NAPONI (4.3.6) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      IF(IATYP.LT.4)THEN
         TAUM=CM*EMT
      ELSEIF(IATYP.EQ.4) THEN
         TAUM=CM*(DETG-1.D0)/3.D0
      ELSEIF(IATYP.EQ.5) THEN
         TAUM=CM*DLOG(DETG)/3.D0
      ENDIF
C     NORMIRANA ELASTICNA DEFORMACIJA
      IF(IATYP.GE.4) THEN
         AEE=(1.D0/2.D0/GE)*AJOT
         DO 210 I=1,3
  210    DEF(I)=TAUD(I)*AEE+EMT
         DO 220 I=4,6
  220    DEF(I)=2.D0*TAUD(I)*AEE
      ENDIF
      TAUM=TAUM/AJOT
C
      APE=DEFQP/EPECRIT
      APSIE=AELL/2.D0/AGEC*(1.d0/3.d0*TEQE*TEQE/GE+3.D0*TAUM*TAUM/CM)
      AGAMMA=2.D0*AEL/AELL
      AS=1.d0-AGAMMA*APE*APSIE+AGAMMA*AGAMMA*APE*APE*
     1   (2.d0*APE-1.d0)*APSIE*APSIE
      AGE1=AS**(2.d0*APE)+AETA
      ADELTA=DABS((AGE1-AGE)/AGE)
      IF(p.gt.100) stop 'lokalne iteraacije>100'
      write(3,*)'ADELTA',ADELTA
      IF(ADELTA.GT.1.d-6) THEN
      write(3,*)'AGE racuna',AGE,AGE1
      AGE=AGE1
      GOTO 222
      ENDIF
C
CS    NAPONI NA KRAJU KORAKA (4.3.56) NELINEARNA ANALIZA KONSTRUKCIJA M.ZIVKOVIC
C
      DO 200 I=1,3
  200 TAU(I)=TAUD(I)+TAUM
      DO 205 I=4,6
         TAU(I)=TAUD(I)
  205 DEFP(I)=2.D0*DEFP(I)
C
CE    UPDATE FROM PREVIOUS STEP
C
C     OVO PROVERITI ZA MALE DEFORMACIJE
      IF(IGLPR.EQ.1) THEN
         IF(NAPKO.EQ.1) THEN
            IF(ILEDE.EQ.1) THEN
C              GLAVNE VREDNOSTI
C              INVERZNO LAMBDA
               write(3,*) 'princ=', princ(1), princ(2), princ(3)
               P1=1.D0/DSQRT(PRINC(1))
               P2=1.D0/DSQRT(PRINC(2))
               P3=1.D0/DSQRT(PRINC(3))
C              INVERZNI DESNI ELASTICNI TENZOR IZDUZENJA (Ue**-1)
               CALL DIJAD(POM,QP,QP,P1,P2,P3)
C              TENZOR ROTACIJE R
C              R = Fe * Ue**-1 
               CALL MNOZM1(VDEF,XJ,POM,3,3,3)
               write(3,*)'vdef=',((vdef(i,j),j=1,3),i=1,3)
CS             TRANSF. UNAZAD ROTIRANI KOSIJEV - KOSIJEV NAPON 
CE             TRANSFORM. BACK ROTATED COUCHY - CAUCHY STRESS
C              s = R * S * RT
               CALL PIOKOS(VDEF,TAU)
               CALL CEPMT(ELAST,CT,0)
C              Cmnop = Vmi Vnj Vok Vpl Cijkl
               CALL RRRRC(ELAST,CT,VDEF,1)
            ENDIF
         ENDIF
C        NAPONI U DEKARTOVOM SISTEMU
         CALL JEDNA1(TAU0,TAU,6)
      ENDIF
C
C     NAPON I DEFORMACIJA ZA STAMPANJE
C
      DO 290 I=1,6
         DEF1(I)=DEF(I)
         IF(IGLPR.EQ.1) THEN
            TAU1(I)=TAU0(I)
         ELSE
            TAU1(I)=TAU(I)
         ENDIF
  290 CONTINUE
C
      IF(ILEDE.EQ.1) THEN
         CALL TRANSS(TSG,QP)
C        Pg=Qs*Pd
         CALL CLEAR(DEF,6)
         CALL MNOZI1(DEF,TSG,DDEFPS,6,6)
C        call wrr6(def,6,'DEFP')
         RETURN
      ENDIF
C
C     KORIGOVANJE NORMIRANOG ELAST. DEF. TENZORA Be NA KRAJU KORAKA
C
      IF(IATYP.EQ.4) THEN
         DO 300 I=1,3
  300    DEF(I)=2.D0*DEF(I)+1.D0
      ELSEIF(IATYP.EQ.5) THEN
         CALL TRANSE(TSG,QP)
C        Eg=Qe*Ed
         CALL CLEAR(SHET,6)
         CALL MNOZI1(SHET,TSG,DEF,6,6)
C        call wrr6(shet,6,'shet')
C        OVO JE PRIBLIZNO ZA MESOVITO OJACANJE
         DO 310 I=1,3
  310    DEF(I)=DEXP(2.D0*SHET(I))
         CALL DIJADS(QP,DEF)
      ENDIF
      RETURN
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TAUI36')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TAUI36')
C-----------------------------------------------------------------------
      END
      SUBROUTINE CEP366K(SHETE,DLAM,TEQ,DEFQP,DDEFQP,EM,GE,CM,CHET,
     &              CY,AN,HS,ISIMO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ELAST )
CE     ELASTO-PLASTIC  CEP MATRIX
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION SHETE(*),SHETE0(6),CP(3,3)

      DO 49 I=1,6                
   49 SHETE0(I)=(SHETE(I)/(2.D0*GE))
     
      EMEPH = 0.D0
      IF(ISIMO.EQ.0) THEN
        EPPRIM=AN*(AN-1.)*CY*DEFQP**(AN-2.D0)      
        IF(EM.GT.1.D-10) EMEPH = EM * AN*CY*(EM*DEFQP)**(AN-1.)      
      ELSE
        EPPRIM=AN*AN*CY*DEXP(-AN*DEFQP)
        IF(EM.GT.1.D-10) EMEPH = EM * (AN*CY*DEXP(-AN*DEFQP)+HS)
      ENDIF
      EPHETP=EMEPH+(1.-EM)*EPPRIM*DDEFQP
C
      AP=((1.D0/2.D0/GE)*((2.D0/3.)*EPHETP+CHET)+1.)
     1   *DSQRT(1.5D0*TENDOT(SHETE0))
      BP=1.5D0/AP/TEQ*(1.-EMEPH*DDEFQP/TEQ)
      DP=(1.D0/2.D0/GE)+(1.+(1.D0/2.D0/GE)*CHET)*DLAM
      DP=(BP-(2.D0/3.)*(1.-EM)*DLAM*DLAM*EPPRIM/AP)/DP/DP
      AA=(1.+CHET*DLAM)/((1.D0/2.D0/GE)+(1.+(1.D0/2.D0/GE)*CHET)*DLAM)
C
      DO 22 I=1,3
        DO 20 J=I,3
   20   CP(I,J)=-DP*SHETE0(I)*SHETE0(J)
        CP(I,I)=CP(I,I)+AA
   22 CONTINUE
      DO 30 I=1,3
      DO 30 J=4,6
        ELAST(I,J)=-DP*SHETE0(I)*SHETE0(J)
   30 CONTINUE
      DO 40 I=4,6
      DO 40 J=I,6
        ELAST(I,J)=-DP*SHETE0(I)*SHETE0(J)
   40 CONTINUE
      DO 45 I=4,6
        ELAST(I,I)=ELAST(I,I)+0.5*AA
   45 CONTINUE
C
      ELAST(1,1)=(2.D0*CP(1,1)-CP(1,2)-CP(1,3)+CM)/3.D0
      ELAST(1,2)=(2.D0*CP(1,2)-CP(1,1)-CP(1,3)+CM)/3.D0
      ELAST(1,3)=(2.D0*CP(1,3)-CP(1,1)-CP(1,2)+CM)/3.D0
      ELAST(2,2)=(2.D0*CP(2,2)-CP(1,2)-CP(2,3)+CM)/3.D0
      ELAST(2,3)=(2.D0*CP(2,3)-CP(1,2)-CP(2,2)+CM)/3.D0
      ELAST(3,3)=(2.D0*CP(3,3)-CP(1,3)-CP(2,3)+CM)/3.D0
C      
      DO 50 I=1,6
      DO 50 J=I,6
        ELAST(J,I)=ELAST(I,J)
   50 CONTINUE
      RETURN
	  END