C=======================================================================
C
C   PLASTICNOST 3/D ELEMENT  -  MESOVITO OJACANJE    (08.08.1992)
C
C=======================================================================
      SUBROUTINE D3M66(TAU,DEF,IRAC,LPOCG,LPOC1)
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
C
      LTAU1 =LPOC1
      LDEFT1=LTAU1  + 6
      LDEFP1=LDEFT1 + 6
      LALFA1=LDEFP1 + 6
      LTEQT1=LALFA1 + 6
      LDQPT1=LTEQT1 + 1
      LIPL1 =LDQPT1 + 1
C
      CALL TAUI366(PLAST(LIPL),PLAST(LDEFPT),PLAST(LALFAT),
     1            PLAST(LTEQT),PLAST(LDQPT),
     1            PLAS1(LIPL1),PLAS1(LTAU1),PLAS1(LDEFT1),PLAS1(LDEFP1),
     1            PLAS1(LALFA1),PLAS1(LTEQT1),PLAS1(LDQPT1),
     1            A(LFUN),MATE,TAU,DEF,IRAC)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE TAUI366( PL ,DEFPT,ALFAT,TEQT,DEFQPT,
     1                   PL1,TAU1,DEF1,DEFP, ALFA1, TEQ, DEFQP,
     1                   FUN,MATE,TAU,DEF,IRAC)
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
     1          DEFP(*),ALFAT(*),ALFA1(*),SHET(6),GLITL(6),POM(3,3),
     1          ALFATG(6),TSS(6,6)
      DIMENSION FUN(2,MATE,*)
      DIMENSION VDEF(3,3),TAU0(6),DDEFPS(6),CT(3,3,3,3),TSG(6,6),
     1          taugl(6)
      DATA ITMAX/100/,EPSIL/1.0D-10/
C
C     INDIKATOR KONTROLNE STAMPE
      IST=0
C
CE  INITIAL DATA
C
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
C
CE    AUXILIARY CONSTANTS
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
      CALL MEL36(FUN,MATE)
C
      IF(IRAC.EQ.2) RETURN
C
CE    YIELD STRESS
C
      IF(ISIMO.EQ.0) THEN
         TEQY =TEQY0+CY*(EM*DEFQPT)**AN
      ELSE
         TEQY =TEQY0+CY*(1.D0-DEXP(-AN*EM*DEFQPT))+HS*EM*DEFQPT
      ENDIF
C     KORISCENJE AJOT 
      AJOT=1.D0
C     SA OVOM KOREKCIJOM DAJE LOSIJE REZULTATE
C      IF(ILEDE.EQ.1) AJOT=DETG
C
C... TRANSFORM ENGENEER. PLASTIC SHEAR STRAIN INTO TENSORIAL
C
      DO 10 I=4,6
   10 DEFPT(I)=.5D0*DEFPT(I)
C
C     D E V I A T O R I C   STRAIN, EPRIM, ESEKUNDUM, GLITL
C
      IF(IATYP.LT.4) THEN
         EMT = (DEF(1)+DEF(2)+DEF(3))/3.D0
      ELSEIF(IATYP.EQ.4) THEN
C         EMT = (DEF(1)+DEF(2)+DEF(3))/3.D0
         EMT=(DETG-1.D0)/3.D0
      ELSEIF(IATYP.EQ.5) THEN
C         EMT=0.D0
         EMT=DLOG(DETG)/3.D0
      ENDIF
      IF(IATYP.LT.4) THEN
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
CE   1)  ELASTIC DEVIATORIC STRESS SOLUTION  (TAUD), (GLITL), (SHET)
C
      IF(IATYP.GE.4.AND.IGLPR.EQ.1) THEN
         CALL TRANSE(TSG,QP)
         IF(IST.EQ.1) call wrr6(tsg,36,'tsg ')
         CALL TRANSS(TSS,QP)
         IF(IST.EQ.1) call wrr6(tss,36,'tss ')
C        Ag=Qs*Ad
         CALL CLEAR(ALFATG,6)
         CALL MNOZI1(ALFATG,TSS,ALFAT,6,6)
      ELSE
         CALL JEDNA1(ALFATG,ALFAT,6)
      ENDIF
      DO 41 I=1,6
         GLITL(I)=DEFDS(I)/AJOT-AE*ALFATG(I)
         TAUD(I) =DEFDS(I)*AEINV/AJOT
   41 SHET(I) =TAUD(I)-ALFATG(I)
      IF(IST.EQ.1) call wrr6(ALFAT,6,'ALFT')
      IF(IST.EQ.1) call wrr6(ALFATG,6,'ALFG')
      IF(IST.EQ.1) call wrr6(defds,6,'defs')
      IF(IST.EQ.1) call wrr6(glitl,6,'glit')
      IF(IST.EQ.1) call wrr6(taud,6,'taud')
C
CE   2)  CHECK FOR YIELDING
C
      TEQ=DSQRT(1.5D0*TENDOT(SHET))
      IF(IST.EQ.1) WRITE(3,*)'emt,teq,teqy',emt,teq,teqy
      IF(DABS(TEQY).GT.1.D-5) THEN
      IF((TEQ-TEQY)/TEQY.LT.1.D-5)THEN
         TEQ  =TEQY
         DEFQP=DEFQPT
         CALL JEDNA1(DEFP,DEFPT,6)
         CALL CLEAR(DDEFPS,6)
         GO TO 500
      ENDIF
      ENDIF
C
CE   3)  SOLUTION IS ELASTO-PLASTIC.  OBTAIN ZERO OF THE ESF.
C
      PL1=1.D0
C
CE    ITERATIONS
C     
      DEFQP=DEFQPT
      IF(DEFQP.LE.1.D-4) DEFQP=1.D-4
      IF(ISIMO.EQ.0) THEN
         EP2=AN*CY*DEFQP**AN1      
      ELSE
         EP2=AN*CY*DEXP(-AN*DEFQP)+HS
      ENDIF
C
      IF(CY.LT.1.D-6) THEN
         DUM=DSQRT(1.5D0*TENDOT(GLITL))
         DDEFQP=DVT*(DUM-AE*TEQY)
         DEFQP=DEFQPT+DDEFQP
         CC=0.D0
         EM1CC=0.D0
         GO TO 199
      ENDIF
      IB = 0
      IT = 0
      AF = 5.D0
      DDD= .1D0*(TEQ-TEQY)/EP2
      GG =AE*TEQY0-DSQRT(1.5D0*TENDOT(GLITL))
      IF(ISIMO.EQ.0) THEN
         AA =AE*CY*(EM**AN)
         FP=-AA*DEFQPT**AN-GG
      ELSE
         FP = -AE*(CY*(1.D0-DEXP(-AN*EM*DEFQPT))+HS*EM*DEFQPT)-GG
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
         WRITE(IZLAZ,2001) NLM,NGR,NGS,NGT
         STOP
      ENDIF
C
      DEFQP=DEFQPT+DDEFQP
      IF(ISIMO.EQ.0) THEN
         IF(IST.EQ.1)WRITE(3,*)'DEFQP,DEFQPT,DDEFQP',DEFQP,DEFQPT,DDEFQP
         CC=CYAN23*DEFQP**AN1
         EM1CC =EM1*CC
         BB =ABB+1.5D0*EM1CC*AE
         FB=-AA*DEFQP**AN-BB*DDEFQP-GG
      ELSE
         CC=CYAN23*DEXP(-AN*DEFQP)+DVT*HS
         EM1CC =EM1*CC
         BB =1.5D0+1.5D0*EM1CC*AE
         FB =-AE*(CY*(1.D0-DEXP(-AN*EM*DEFQP))-HS*EM*DEFQP)-BB*DDEFQP-GG
      ENDIF
C
      CALL BISEC (DDEFQP,DEPBM,DEPBP,DDD,FB,FM,FP,AF,IB)
      IF (IB1.EQ.0) GO TO 100
      IF (DABS(DDD).GT.EPSIL.AND.
     1    (DABS(DDD)/(DEPBM+DEPBP)).GT.EPSIL) GO TO 100
C
CE      ...   ( DEVIATORIC STRESS )
C
  199 IF(ISIMO.EQ.0) THEN
         TEQ =TEQY0+CY*(EM*DEFQP)**AN
      ELSE
         TEQ =TEQY0+CY*(1.D0-DEXP(-AN*EM*DEFQP))+HS*EM*DEFQP
      ENDIF
      DLAM=1.5D0*DDEFQP/TEQ
      DUM =EM1CC*DLAM
      DUH =1.D0/(AE+DLAM+DUM*AE)
      DUM =1.D0+DUM
      DO 165 I=1,6
         SHET(I)=DUH*GLITL(I)
  165 TAUD(I)=ALFATG(I)+DUM*SHET(I)
C
CE   4)  DETERMINE SOLUTION 
C
C
CE      ...   ( PLASTIC STRAIN ), ( BACK STRESS )
C
      DO 170 I=1,6
         DDEFPS(I)=DLAM*SHET(I)
  170 CONTINUE
      IF(IST.EQ.1) call wrr6(DDEFPS,6,'DEFP')
      IF(IATYP.GE.4.AND.IGLPR.EQ.1) THEN
C        Pd=QeT*Pg
         CALL CLEAR(ALFATG,6)
         CALL MNOZI2(ALFATG,TSG,DDEFPS,6,6)
      ELSE
         CALL JEDNA1(ALFATG,DDEFPS,6)
      ENDIF
      DO 171 I=1,6
         DEFP(I) =DEFPT(I)+ALFATG(I)
         ALFA1(I)=ALFAT(I)+EM1CC*ALFATG(I)
  171 CONTINUE
C
CE     E L A S T I C  -  P L A S T I C   M A T R I X   CEP
C
      IF(ISKNP.NE.2) THEN
      IF(IST.EQ.1) call wrr6(GLITL,6,'GLIT')
      IF(IATYP.GE.4.AND.IGLPR.EQ.1) THEN
         CALL JEDNA1(ALFATG,GLITL,6)
C        Ed=QeT*Eg
         CALL CLEAR(GLITL,6)
         CALL MNOZI2(GLITL,TSG,ALFATG,6,6)
      ENDIF
      IF(IST.EQ.1) call wrr6(GLITL,6,'GLID')
         IF(INDCEP.EQ.0)
     1   CALL CEP366(GLITL,DLAM,TEQ,DEFQP,DDEFQP,EM,AE,CM,CC,CY,AN,
     &              HS,ISIMO)
      ENDIF
C
CE   5)    CALCULATE STRESS
C
  500 CONTINUE
C
C     SREDNJI NAPON
         TAUM=CM*EMT
C     NORMIRANA ELASTICNA DEFORMACIJA
      IF(IATYP.GE.4) THEN
         AEE=AE*AJOT
         DO 210 I=1,3
  210    DEF(I)=TAUD(I)*AEE+EMT
         DO 220 I=4,6
  220    DEF(I)=2.D0*TAUD(I)*AEE
      ENDIF
      TAUM=TAUM/AJOT
      DO 200 I=1,3
  200 TAU(I)=TAUD(I)+TAUM
      DO 205 I=4,6
         TAU(I)=TAUD(I)
  205 DEFP(I)=2.D0*DEFP(I)
      IF(IST.EQ.1) call wrr(taud,6,'taud')
      IF(IST.EQ.1) call wrr(tau,6,'tau ')
      IF(IST.EQ.1) call wrr(defp,6,'defp')
      IF(IST.EQ.1) write(3,*) 'emt,e,ani',emt,e,ani
C
CE  UPDATE FROM PREVIOUS STEP
C
C OVO PROVERITI ZA MALE DEFORMACIJE
      CALL JEDNA1(SHET,DEF,6)
      IF(IST.EQ.1) call wrr6(SHET,6,'SHET')
      IF(IGLPR.EQ.1) THEN
C        NAPONI U GLAVNIM PRAVCIMA
         CALL JEDNA1(TAU0,TAU,6)
         call jedna1(taugl,tau,6)
C        NAPONI U DEKARTOVOM SISTEMU
C        Sd=QeT*Sg
         CALL CLEAR(TAU,6)
         CALL MNOZI2(TAU,TSG,TAU0,6,6)
C         call wrr(tau0,6,'tau0')
C         call wrr(tau,6,'tau ')
C         call wrr6(tsg,36,'tsg ')
         IF(NAPKO.EQ.1) THEN
CV            IF(IATYP.EQ.4.AND.ILEDE.EQ.0) THEN
C              GLAVNE VREDNOSTI
C              LAMBDA
CV               P1=DSQRT(PRINC(1))
CV               P2=DSQRT(PRINC(2))
CV               P3=DSQRT(PRINC(3))
C              LEVI KOSI-GRINOV DEFORMACIONI TENZOR (V)
CV               CALL DIJAD(VDEF,QP,QP,P1,P2,P3)
CS             TRANSF. ROTIRANI PIOLA KIRKOFOV - KOSIJEV NAPON 
CE             TRANSFORM. ROTATED PIOLA KIRCKOF - CAUCHY STRESS
C              s = V * S * V
CV               CALL PIOKOS(VDEF,TAU0)
CV               CALL JEDNA1(TAU,TAU0,6)
CV               CALL CEPMT(ELAST,CT,0)
C               IF(NLM.EQ.1.AND.JG.EQ.1) THEN
C                  WRITE(3,*) 'NLM,JG,ITER',NLM,JG,ITER
C                  CALL WRR3(VDEF,9,'VDEF')
C                  CALL WRR6(ELAST,36,'ELAP')
C                  CALL WRRT4(CT,1,1,2,2,3,3,'CTDP')
C                  CALL WRRT4(CT,1,2,2,3,1,3,'CTGP')
C                  CALL WRRT4(CT,2,1,3,2,3,1,'CTDP')
C               ENDIF
C              Cmnop = Vmi Vnj Vok Vpl Cijkl
CV               CALL RRRRC(ELAST,CT,VDEF,1)
C               IF(NLM.EQ.1.AND.JG.EQ.1) THEN
C                  CALL WRR6(ELAST,36,'ELAS')
C                  CALL WRRT4(CT,1,1,2,2,3,3,'CTD ')
C                  CALL WRRT4(CT,1,2,2,3,1,3,'CTG ')
C                  CALL WRRT4(CT,2,1,3,2,3,1,'CTD ')
C               ENDIF
CV            ENDIF
            IF(ILEDE.EQ.1) THEN
C              GLAVNE VREDNOSTI
C              INVERZNO LAMBDA
               P1=1.D0/DSQRT(PRINC(1))
               P2=1.D0/DSQRT(PRINC(2))
               P3=1.D0/DSQRT(PRINC(3))
C              INVERZNI DESNI ELASTICNI TENZOR IZDUZENJA (Ue**-1)
               CALL DIJAD(POM,QP,QP,P1,P2,P3)
C              TENZOR ROTACIJE R
C              R = Fe * Ue**-1 
               CALL MNOZM1(VDEF,XJ,POM,3,3,3)
CS             TRANSF. UNAZAD ROTIRANI KOSIJEV - KOSIJEV NAPON 
CE             TRANSFORM. BACK ROTATED COUCHY - CAUCHY STRESS
C              s = R * S * RT
               CALL PIOKOS(VDEF,TAU)
               CALL CEPMT(ELAST,CT,0)
C              Cmnop = Vmi Vnj Vok Vpl Cijkl
               CALL RRRRC(ELAST,CT,VDEF,1)
            ENDIF
         ENDIF
         IF(NAPGL.EQ.0) CALL JEDNA1(TAU0,TAU,6)
         CALL JEDNA1(ALFATG,DEF,6)
         CALL CLEAR(DEF,6)
C        Ed=QsT*Eg
         CALL MNOZI2(DEF,TSS,ALFATG,6,6)
C         call wrr(alfatg,4,'defa')
C         call wrr(def,4,'def ')
C         call wrr6(tss,36,'tss ')
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
      IF(IST.EQ.1) call wrr(def1,6,'def1')
      IF(IST.EQ.1) call wrr(tau1,6,'tau1')
C
      IF(ILEDE.EQ.1.OR.(ILEDE.EQ.0.AND.ICPM1.EQ.2)) THEN
         CALL JEDNA1(DEF,DDEFPS,6)
         RETURN
      ENDIF
C
C     KORIGOVANJE NORMIRANOG ELAST. DEF. TENZORA Be NA KRAJU KORAKA
C
      IF(IATYP.EQ.4) THEN
         DO 300 I=1,3
  300    DEF(I)=2.D0*DEF(I)+1.D0
      ELSEIF(IATYP.EQ.5) THEN
C         call wrr6(shet,6,'shet')
C        OVO JE PRIBLIZNO ZA MESOVITO OJACANJE
         DO 310 I=1,3
  310    DEF(I)=DEXP(2.D0*SHET(I))
         CALL DIJADS(QP,DEF)
      ENDIF
c     dijada nxn, izvod dmdb, cep
cg      IF(ISKNP.NE.2.AND.INDCEP.EQ.0.and.IATYP.GE.4.AND.IGLPR.EQ.1) 
cg     1   call formd(def,taugl)
      RETURN
C-----------------------------------------------------------------------
 2001 FORMAT( ' ELEMENT =',I6,'  IR =',I2,'  IS =',I2,'  IT =',I2)
 2000 FORMAT(/' DOSTIGNUT MAKSIMALAN BROJ BISEKCIJA U TAUI36')
C-----------------------------------------------------------------------
 6000 FORMAT(/' MAXIMUM NUMBER OF BISECTION IS REACHED IN TAUI36')
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE CEP366(GLITL,DLAM,TEQ,DEFQP,DDEFQP,EM,AE,CM,CC,CY,AN,
     &                 HS,ISIMO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICE CEP ( ELAST )
CE     ELASTO-PLASTIC  CEP MATRIX
C
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      DIMENSION GLITL(*),CP(3,3) 
C
      TRIPO =1.5D0
      DVA   =2.D0
      TR    =1.D0/3.
      DVT   =DVA*TR
      EM1   =1.-EM
      AN1   =AN-1.
      AN2   =AN-DVA
      CHET  =EM1*CC
C
      EMEPH = 0.D0
      IF(ISIMO.EQ.0) THEN
        EPPRIM=AN*AN1*CY*DEFQP**AN2      
        IF(EM.GT.1.D-10) EMEPH = EM * AN*CY*(EM*DEFQP)**AN1      
      ELSE
        EPPRIM=AN*AN*CY*DEXP(-AN*DEFQP)
        IF(EM.GT.1.D-10) EMEPH = EM * (AN*CY*DEXP(-AN*DEFQP)+HS)
      ENDIF
      EPHETP=EMEPH+EM1*EPPRIM*DDEFQP
C
      AP=(AE*(DVT*EPHETP+CHET)+1.)*DSQRT(TRIPO*TENDOT(GLITL))
      BP=TRIPO/AP/TEQ*(1.-EMEPH*DDEFQP/TEQ)
      DP=AE+(1.+AE*CHET)*DLAM
      DP=(BP-DVT*EM1*DLAM*DLAM*EPPRIM/AP)/DP/DP
      AA=(1.+CHET*DLAM)/(AE+(1.+AE*CHET)*DLAM)
C
      DO 22 I=1,3
        DO 20 J=I,3
   20   CP(I,J)=-DP*GLITL(I)*GLITL(J)
        CP(I,I)=CP(I,I)+AA
   22 CONTINUE
      DO 30 I=1,3
      DO 30 J=4,6
        ELAST(I,J)=-DP*GLITL(I)*GLITL(J)
   30 CONTINUE
      DO 40 I=4,6
      DO 40 J=I,6
        ELAST(I,J)=-DP*GLITL(I)*GLITL(J)
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
C      END
CC
CCE       FULL NEWTON ITERATIONS
CC     
C      GG =AE*TEQY0-DSQRT(1.5*TENDOT(GLITL))
C      AA =AE*CY*(EM**AN)
CC
C      DDEFQP=0
C      IT=0
C  100 IT=IT+1
CC
C      IF(IT.GT.ITMAX) THEN
C        WRITE(IZLAZ,2000)
C        STOP
C      ENDIF
CC
C      DEFQP=DEFQPT+DDEFQP
C      IF(DEFQP.LE.1.D-8) DEFQP=1.D-8
C      CC=CYAN23*DEFQP**AN1
C      EM1CC =EM1*CC
C      BB =1.5+1.5*EM1CC*AE
C      CCP=1.5*AE*EM1*DDEFQP*CYAN23*AN1*DEFQP**AN2
CC
C      DDD=-(AA*(DEFQP**AN)+BB*DDEFQP+GG)/(AA*(DEFQP**AN1)+BB+CCP)
C      DDEFQP=DDEFQP+DDD
C      IF(DABS(DDD/DEFQP).LT.EPSIL) GO TO 150
C      GO TO 100
      END
