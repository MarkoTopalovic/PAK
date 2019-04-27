C* SILE93 ----> RETURN
C=======================================================================
C
C   SUBROUTINE S993GL
C              RD93
C              SIS93
C              ELTE93
C              ALFC3
C=======================================================================
      SUBROUTINE S993GL(IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     GLAVNI UPRAVLJACKI PROGRAM ZA FORMIRANJE MATRICA I VEKTORA 
CS     3/D KONTAKT ELEMENATA
CE     MAIN PROGRAM FOR CALCULATION OF ELEMENT MATRIX FOR 
CE     3/D CONTACT ELEMENTS
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ELEMAE/ MXAE,LAE,LMXAE,LHE,LBET,LBED,LRTHE,LSKE,LLM
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /DUPLAP/ IDVA
      COMMON /SRPSKI/ ISRPS
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /KONTK3/ FMSTAT,FMDIN,EPSIL,NTSF,NTTOT,LNCVSF,LITSRF,
     &                LNELSF,LIDC,LIK,LIK1,LFSFD,LMASE,LNELAB,NWKCDY
      COMMON /RSN   / DETER,IPIVOT,IDETER
      COMMON /PENALTY/ AKN,AKS,IPENALTY
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' S993GL'
C      WRITE(3,*)'>>>>>>>>>>>>>>>>>>>>>>>IND',IND
C*2
      IPIVOT=1
C*2     
      LAU=LMAX
      CALL RD93(A(LAU),IND)
      IF(MOD(LMAX,2).EQ.0) LMAX=LMAX+1
C
      LAE=LMAX
      NCVE36=NCVE*3+6
      MXAE = 100*IDVA+NCVE36+1
      MMMM = NTTOT*3*IDVA+2*NWKCDY*IDVA+NCVE36+NCVE*NCVE+NCVE+1
      IF(NDIN.NE.0) MXAE = MAX0(MXAE,MMMM)
C... KONTROLA MEMORIJE			   
      CALL CTRREP(LMAX,LAE,MXAE,MTOT)
C
C  
      IF (IPENALTY.EQ.0) CALL SIS93(A(LAE),A(LAU),IND)
      IF (IPENALTY.EQ.1) CALL SIS93P(A(LAE),A(LAU),IND)      
C
C
      RETURN
      END
C=======================================================================
      SUBROUTINE RD93(AU,INA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     UPRAVLJACKI PROGRAM ZA UCITAVANJE ULAZNIH PODATAKA U AU
CE     MENAGEMENT ROUTINE FOR INPUT DATA IN    AU
C
      include 'paka.inc'
      
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /KONTK0/ IRESCT
      COMMON /KONTK3/ FMSTAT,FMDIN,EPSIL,NTSF,NTTOT,LNCVSF,LITSRF,
     &                LNELSF,LIDC,LIK,LIK1,LFSFD,LMASE,LNELAB,NWKCDY
C
C... IRESCT = INDIKATOR "RESTARTA KONTAKTA" (AUTOMATSKO OPTERECIVANJE)
C
      DIMENSION AU(*)
      REAL AU
C
CS     CITANJE SKALARA IZ DIREKT ACCES FILE 
CE     READ SCALARS FROM A DIRECT ACCESS FILE
C
      LSTAZA(1)=LMAX8
      READ(IELEM,REC=LMAX8)
     1 NTSF,NTTOT,NCVE,MXAU,NWKCDY,LNCVSF,LITSRF,
     1 LISNA,LNELSF,LNEL,LLMEL,LFSFD,LMASE,LNELAB,EPSIL
      LSTAZA(2)=LMAX8+1
CC      WRITE(3,*)
CC     1 'LMAX8,LAU,NTSF,NTTOT,NCVE,MXAU,LNCVSF,LITSRF'
CC      WRITE(3,*)
CC     1 LMAX8,LAU,NTSF,NTTOT,NCVE,MXAU,LNCVSF,LITSRF
CC      WRITE(3,*)
CC     1 'LISNA,LNELSF,LNEL,LLMEL,LFSFD,LMASE,LNELAB'
CC      WRITE(3,*)
CC     1 LISNA,LNELSF,LNEL,LLMEL,LFSFD,LMASE,LNELAB
      CALL READDD(AU(LFSFD),MXAU/IDVA,IELEM,LMAX8,LDUZI)
      LMAX=LAU+MXAU
      IF(MOD(LMAX,2).EQ.0) LMAX=LMAX+1
CS....  VELICINE SA KRAJU PRETHODNOG KORAKA
CE....  PREVIOUS STEP VALUES
      NPROS=(NE*3+1)/2*2/IDVA+5*NE
      LIK  =LMAX
      LIK1 =LIK +NE
      IF(MOD(LIK1,2).EQ.0) LIK1=LIK1+1
      LMAX =LIK1+2*NE+5*NE*IDVA
      IF(LMAX.GT.MTOT) CALL ERROR(1)
CS....  VELICINE ZA TEKUCI KORAK
      LSTAZA(3)=LMAX8+1
      LMA8=LMAX8
      CALL READDD(A(LIK), NPROS,IELEM,LMAX8,LDUZI)
C  ARC-LENGTH RESTART
      IF(IRESCT.GT.0)THEN
        CALL IJEDN1(A(LIK1),A(LIK),NE)
CC        LXX=LIK1+2*NE
CC        NXX=5*NE
CC        CALL CLEAR(A(LXX),NXX)
        CALL WRITDD(A(LIK),NPROS,IELEM,LMA8,LDUZI)
        IRESCT=0
      ENDIF       
      IF(INA.EQ.4)THEN
        CALL IJEDN1(A(LIK),A(LIK1),NE)
        CALL WRITDD(A(LIK),NPROS,IELEM,LMA8,LDUZI)
      ENDIF
C
      LSIGMA=LMAX
      NPROS =(NE+NTTOT)*3
      LSIGT =LSIGMA+NPROS*IDVA
      LMAX  =LSIGT +NPROS*IDVA
      LSTAZA(5)=LMAX8+1
      NPRO2=2*NPROS
      IF(INA.EQ.2.OR.INA.EQ.3)THEN
        CALL READDD(A(LSIGMA),NPRO2,IELEM,LMAX8,LDUZI)
        CALL CLEAR (A(LSIGMA),NPROS)
      ELSEIF(INA.EQ.4)THEN
        LMA8=LMAX8
        CALL READDD(A(LSIGMA),NPRO2,IELEM,LMAX8,LDUZI)
        CALL JEDNA1(A(LSIGT),A(LSIGMA),NPROS)
        CALL WRITDD(A(LSIGMA),NPRO2,IELEM,LMA8,LDUZI)
      ELSEIF(INA.EQ.5.OR.INA.EQ.7.OR.INA.EQ.8.OR.INA.EQ.9)THEN
        CALL READDD(A(LSIGMA),NPRO2,IELEM,LMAX8,LDUZI)
      ENDIF
      RETURN
      END
C======================================================================
      SUBROUTINE SIS93(AE,AU,IND)
      USE MATRICA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     GLAVNI UPRAVLJACKI PROGRAM  ZA MATRICE ELEMENATA I SISTEMA
CE     MAIN MANAGEMENT  PROGRAM  FOR ELEMENT MATRIX
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEMAE/ MXAE,LAE,LMXAE,LHE,LBET,LBED,LRTHE,LSKE,LLM
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /DUPLAP/ IDVA
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /KONTK3/ FMSTAT,FMDIN,EPSIL,NTSF,NTTOT,LNCVSF,LITSRF,
     &                LNELSF,LIDC,LIK,LIK1,LFSFD,LMASE,LNELAB,NWKCDY
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /DIM   / N9,N10,N11,N12,MAXUP
      COMMON /UPRIRI/ LUPRI
      COMMON /PROBAS/ IILS
      COMMON /ITERBR/ ITER
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION AE(*),AU(*)
      REAL AE,AU
      IF(IDEBUG.GT.0) PRINT *, ' SIS93'
C  DOPUNSKI REPERI U VEKTORU A
      LIK2 =LIK1 +NE
      LALFK=LIK2 +NE
      KORD=LCORD
      LNGLOB=LMASE +(NE+NTTOT)*IDVA
      LJDIAG=LNGLOB+NTTOT
      LXYZ=1
      LHE =LXYZ+3*(NCVE+1)*IDVA
      LTRA=LHE+27*IDVA
      LAE0=LTRA+9*IDVA
C... PROVERITI ZASTO ZEZA KOD IMPACT I ZA STA JE NEOPHODNO !!!!!
CC      IF(KOR.EQ.1.AND.ITER.EQ.0) CALL CLEAR(A(LRCTDT),JEDN)
CC      CALL CLEAR(A(LRCTDT),JEDN)
      CALL CLEAR(AE,LAE0/IDVA)
C
CS            P E T L J A    P O    E L E M N T I M A
CE            L O O P    O V E R    E L E M E N T S
C
      NCVE36=NCVE*3+6
      IF(IND.EQ.2.OR.IND.EQ.3)THEN
C       REPERI U VEKTORU ELEMENATA AE
       NWE =100
       LSKE=LAE0
       LLM =LSKE+NWE*IDVA
       MXAE=LLM +NCVE36-1
C
C       call WRR(A(LRTDT),JEDN,'U   ')
C       call WRR(A(LUPRI),JEDN,'UPRI')
       ICCMOV=0
       DO 100 NLM=1,NE
       CALL CLEAR(AE,MXAE/IDVA)
       CALL ELTE93(AE(LSKE),AU(LNCVSF),AU(LITSRF),AU(LNELSF),AU(LISNA),
     &             AE(LLM),AU(LIDC),AU(LNEL),A(KORD),A(LID),A(LRTDT),
     &             A(LRCTDT),A(LUPRI),A(LIK),A(LIK1),A(LALFK),AU(LFSFD),
     &             A(LSIGMA),AU(LNELAB),AE(LXYZ),AE(LHE),AE(LTRA),EPSIL,
     &             ICCMOV)
  100  CONTINUE
C       call WRR(A(LFTDT),JEDN,'FTDT')
C       call WRR(A(LRTDT),JEDN,'RTDT')
C       call WRR(A(LRCTDT),JEDN,'RC  ')
C       WRITE(3,*)'NWK,NWKC,NEQC,NEQ',NWK,NWKC,NEQC,NEQ
C       CALL IWRR(A(LMAXA),JEDN+1,'MAXA')
C       call WRR(A(LSK),NWK,'SKC ')
CE    STORE DISPLACEMENTS ON DISK
CS    ZAPISIVANJE POMERANJA NA DISK
C070794
      IF(IILS.NE.-1)CALL WSTAZ(A(LIPODS),LRTDT,52)
      ELSEIF(IND.EQ.5)THEN
C
CS     PRIPADAJUCE MASE CVOROVA
CE     NODES MASSES
C
C       call WRR(A(LSK),NWK,'MASE')
       CALL MASE93(ALSK,A(LMAXA),A(LID),AU(LMASE),AU(LNEL),AU(LNGLOB),
     &             NP,NE,NTTOT)
C
       LMA8=LSTAZA(2)-1
       CALL WRITDD(AU(LFSFD),MXAU/IDVA,IELEM,LMA8,LDUZI)
       RETURN
      ELSEIF(IND.EQ.7.OR.IND.EQ.8)THEN
C       REPERI U VEKTORU ELEMENATA AE
       LS  =LAE0
       LA  =LS+NCVE*NCVE*IDVA
       LC  =LA+NWKCDY*IDVA
       LB  =LC+NWKCDY*IDVA
       LLM =LB+NTTOT*3*IDVA
       LLD =LLM+NCVE36
       MXAE=LLD+NCVE-1
       MCLR=(LLM-LA)/IDVA
       CALL CLEAR(AE(LA),MCLR)
C
CS     KOREKCIJA UBRZANJA I BRZINA PRI UDARU
CE     UPDATE ACCELERATION AND VELOCITY AFFTER IMPACT
C
C   LPUU --> UBRZANJE ITERACIJA (I-1), BRZINA TRENUTAK (T)
C   LRR  --> ( LPUU ) + KOREKCIJE
       LRR =LPUU
       LSIG=LSIGMA
       LII =LIK2
       IF(IND.EQ.8)THEN
C*
         LPUU =LPUV
C010795         LPUU =LPUV
C*
         LRR =LPUV
         LSIG=LSIGMA+(NE+NTTOT)*3*IDVA
         LII =LIK
       ENDIF
C       call WRR(A(LPUU),JEDN,'PUU-')
C       call WRR(A(LRCTDT),JEDN,'RC- ')
       CALL VAUPD3(A(LRR),A(LPUU),A(LID),A(LRCTDT),AU(LMASE),AE(LA),
     &             AE(LC),AE(LB),AU(LISNA),A(LII),A(LIK1),A(LALFK),
     &             AU(LNELSF),AU(LITSRF),AU(LIDC),AU(LNEL),AU(LNCVSF),
     &             A(LSIG),AE(LLM),A(LRTDT),NP,NE,NTTOT,NTSF,DT,IND,
     &             A(LIK),AU(LJDIAG),AU(LNGLOB),AU(LNELAB),
     &             AE(LXYZ),AE(LHE),AE(LTRA),AE(LS),AE(LLD),NCVE)
C       call WRR(A(LPUU),JEDN,'PUU+')
C       call WRR(A(LRCTDT),JEDN,'RC+ ')
       CALL IJEDN1(A(LIK2),A(LIK1),NE)
C
CS     KOREKCIJA PRIRASTAJA SILA
CE     UPDATE FORCE INCREMENT
C
      ELSEIF(IND.EQ.9)THEN
        LLL=LUPRI
        IF(METOD.GT.5) LLL=N9
        CALL SILE93(AU(LIDC),A(LLL),A(LIK1),A(LALFK),NE)
      ENDIF
C
      NPROS=(NE*3+1)/2*2/IDVA+5*NE
      LMA8=LSTAZA(3)-1
      IF(IILS.EQ.-1)RETURN
      CALL WRITDD(A(LIK),NPROS,IELEM,LMA8,LDUZI)
C
      NPROS =(NE+NTTOT)*6
      LMA8=LSTAZA(5)-1
C      call WRR(A(LSIGMA),NPROS,'SIG2')
      CALL WRITDD(A(LSIGMA),NPROS,IELEM,LMA8,LDUZI)
C
      IF(IND.EQ.2)
     & WRITE(IZLAZ,*)'*** CONTACT NODES CHANCHED STATUS:',ICCMOV,'  ***'
C
      RETURN
      END
C=======================================================================
      SUBROUTINE ELTE93(SKE,NCVSF,ITSRF,NELSF,ISNA,LM,IDC,
     &                  NEL,CORD,ID,U,RC,UPRI,IK,IK1,ALFK,FSFD,SILE,
     &                  NELAB,XYZ,HE,TRA,EPSIL,ICCMOV)
      USE MATRICA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICA ELEMENATA I SISTEMA
CE     FORM ELEMENT MATRIX
C
C  IK     INDIKATOR KONTAKTA U TRENUTKU (T)
C  IK1    INDIKATOR KONTAKTA U TRENUTKU (T+DT) (U ITERACIJI)
C  NCAA   REDNI BROJ POLIGONA NA CILJNOJ POVRSINI
      include 'paka.inc'
      
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PROBAS/ IILS
C
      DIMENSION SKE(*),NCVSF(*),ITSRF(*),NELSF(NCVE,*),LM(*),NEL(NE,*),
     1          CORD(NP,*),U(*),RC(*),IDC(NE,*),ID(NP,*),FSFD(NE,*),
     2          IK(*),IK1(*),ALFK(5,*),SILE(3,*),NELAB(NCVE,*),ISNA(*),
     3          XYZ(3,*),HE(9,*),TRA(3,*),UPRI(*)
      DIMENSION DC(3),D(3),LM3(3),ILM(3),UPRI0(3)
      DIMENSION TAURSN(3,5),TAU123(3,5)
      DATA IT1/1000000/,SMALL/1.D-10/,BIG/1.D+10/,TCMAX/0.D0/
C      call IWRR(IK,NE,'IK  ')
C      call IWRR(IK1,NE,'IK1 ')
C
C  CVOR KONTAKTOR
C
      NPOL=1
      NCK=NEL(NLM,1)
C        WRITE(3,*)'*********NCK',NCK
C
C  KOEFICIJENTI TRENJA
C
      FS=FSFD(NLM,1)
      FD=FSFD(NLM,2)
C
C   JEDINICNI VEKTOR PRIRASTAJA POMERANJA
      DO 2 K=1,3
      KK=ID(NCK,K)
      UPRI0(K)=0.D0
2     IF(KK.GT.0) UPRI0(K)=UPRI(KK)
      IF(DABS(UPRI0(1)).GT.SMALL.OR.DABS(UPRI0(2)).GT.SMALL.OR.
     &   DABS(UPRI0(3)).GT.SMALL) CALL JEDV(UPRI0(1),UPRI0(2),UPRI0(3))
C      WRITE(3,*)' UPRI0 ',UPRI0
C
C
CS...... PETLJA PO CILJNIM POLIGONIMA
CE...... LOOP OVER TARGET SEGMENTS
C
C  ISRF - CILJNI SEGMENT, NNC - BROJ CVOROVA SEGMENTA
      ISRF=NEL(NLM,2)
      MAXCE=ISNA(ISRF)
      NODMAX=MAXCE+1
      NNC =NCVSF(ISRF)
      LL  =IABS(ITSRF(ISRF))-1
C  PROVERA STANJA U PRETHODNOJ ITERACIJI
      IND =IK1(NLM)/IT1
      IND1=IK(NLM)/IT1
C  KOREKCIJA KOORDINATA ZA CVOR KONTAKTOR
      ICAL=2
      CALL UPXYZ(XYZ,NCK,NELSF,CORD,U,ID,NP,MAXCE,ICAL)
      ICAL=1
      ICHEL=0
      IF(IND.GT.0) THEN
        NCAA=IK1(NLM)-IND*IT1
        NCA =NCAA+LL
C   KOREKCIJA KOORDINATA ZA CILJNE CVOROVE
        CALL UPXYZ(XYZ,NCK,NELSF(1,NCA),CORD,U,ID,NP,MAXCE,ICAL)
      NNOD=0
      XL2M=1.D15
      XL2=0.D0
       DO 4 NN=1,MAXCE
         DO 3 L=1,3
           D(L)=XYZ(L,NODMAX)-XYZ(L,NN)
    3    CONTINUE
         XL2=D(1)*D(1)+D(2)*D(2)+D(3)*D(3)
         IF(XL2.LT.XL2M)THEN
          XL2M=XL2
          NNOD=NN
         ENDIF
    4  CONTINUE
        IFLAG=1
        CALL ALFC3 (HE,TRA,XYZ,UPRI0,MAXCE,DC,R0,S0,EPSIL,PRODOR,
     &              TOLPD,NNOD,TCMAX,NPOL,NPOL,NCTC,NPOLMX,INDD,IFLAG)
        IF(INDD.EQ.0)THEN
          IND=0
          ICHEL=1
        ENDIF
        IOVR=0
      ENDIF
C      WRITE(3,*)'IND,IND1,INDD,IK1(NLM)/IT1',IND,IND1,INDD,IK1(NLM)/IT1
C
C   NA KRAJU PRETHODNE ITERACIJE NIJE BILO KONTAKTA
C
      IF(IND.EQ.0) THEN
C
C  1) NALAZENJE NAJBLIZEG CVORA   (NNOD)
C
      NNOD=0
      XL2M=1.D15
      DO 10 NC=1,NNC
       NC2 =NC+LL
C   KOREKCIJA KOORDINATA ZA CILJNE CVOROVE
       CALL UPXYZ(XYZ,NCK,NELSF(1,NC2),CORD,U,ID,NP,MAXCE,ICAL)
       XL2=0.D0
       DO 7 NN=1,MAXCE
         DO 5 L=1,3
           D(L)=XYZ(L,NODMAX)-XYZ(L,NN)
    5    CONTINUE
         XL2=D(1)*D(1)+D(2)*D(2)+D(3)*D(3)
         IF(XL2.LT.XL2M)THEN
          XL2M=XL2
          NNOD=NELSF(NN,NC2)
         ENDIF
    7  CONTINUE
   10 CONTINUE
C
C  2)     NALAZENJE POLIGONA   (NCAA)
C
      IFLAG=1
      DO 15 NCAA=1,NNC
       NCA =NCAA+LL
C   KOREKCIJA KOORDINATA ZA CILJNE CVOROVE
       CALL UPXYZ(XYZ,NCK,NELSF(1,NCA),CORD,U,ID,NP,MAXCE,ICAL)
       DO 12 NN=1,MAXCE
         IF(NELSF(NN,NCA).EQ.NNOD)THEN
           CALL ALFC3 (HE,TRA,XYZ,UPRI0,MAXCE,DC,R0,S0,EPSIL,PRODOR,
     &                 TOLPD,NN,TCMAX,NCAA,NPOL,NCTC,NPOLMX,IND,IFLAG)
           IF(IND.NE.0)GO TO 16
         ENDIF
   12  CONTINUE
   15 CONTINUE
       IF(IFLAG.EQ.2)THEN
         NCA =NPOL+LL
         CALL UPXYZ(XYZ,NCK,NELSF(1,NCA),CORD,U,ID,NP,MAXCE,ICAL)
         CALL ALFC3 (HE,TRA,XYZ,UPRI0,MAXCE,DC,R0,S0,EPSIL,PRODOR,
     &               TOLPD,NN,TCMAX,NPOL,NPOL,NCTC,NPOLMX,IND,IFLAG)
       ENDIF
   16 IOVR=IND
C.. ICHEL > 0   AKO JE PROMENJEN PODSEGMENT
      ICHEL=IND*ICHEL
C      WRITE(3,*)'nlm,NCK,DC,IND',NLM,NCK,DC,IND
      ENDIF
       LM3(1)=IDC(NLM,1)
       LM3(2)=IDC(NLM,2)
       LM3(3)=IDC(NLM,3)
   20  IF((ICHEL.NE.0.OR.IND.EQ.0).AND.IILS.NE.-1)THEN
CC        DO 22 J=1,3
CC         IJ=LM3(J)
CC         IF(IJ.GT.0)U(IJ)=0.D0
CC   22   CONTINUE
       ENDIF
C
C  3)     RAZLICITI KONTAKTNI USLOVI
C
C
C-----    (A)   NEMA KONTAKTA
C
      MDIM=MAXCE*3+6
      IF(IND.EQ.0)THEN
       IF(IILS.EQ.-1)RETURN
       MDIM=3
       SKE(1)=-1.D0
       SKE(4)=-1.D0
       SKE(6)=-1.D0
C
       CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM3,MDIM)
C
       DO 30 J=1,3
        IJ=LM3(J)
        IF(IJ.GT.0)U(IJ)=0.D0
        SILE(J,NLM)=0.D0
   30  CONTINUE
       IK (NLM)=0
       IK1(NLM)=0
       RETURN
      ENDIF
C
C-----    (B)   KONTAKT SA ILI BEZ KLIZANJA
C
      KK=-MAXCE-1
      LMILM=0
      DO 35 K=1,3
        KK=KK+MAXCE+2
        LM(KK)  =ID(NCK,K)
        ILM(K)  =IDC(NLM,K)
        LM(KK+1)=ILM(K)
        IF(LMILM.EQ.0) LMILM=ILM(K)
        K1=KK+1
          DO 32 L=1,MAXCE
          NCV=NELSF(L,NCA)
C      IF(K.EQ.1)WRITE(3,*)'NCV',NCV
   32     LM(K1+L)=ID(NCV,K)
          DO 34 J=1,NODMAX
   34     TAU123(K,J)=0.D0
   35 CONTINUE
C
C           KONTAKTNE SILE  - PRAVAC X,Y,Z U KONTAKTORU
C
      KK=-MAXCE
      DO 45 K=1,3
       KK=KK+MAXCE+2
       LMKK=LM(KK)
       TAU123(K,NODMAX)=0.D0
C       WRITE(3,*)'IOVR',IOVR
C       IF(ICHEL.EQ.0)THEN
         IF(LMKK.NE.0) TAU123(K,NODMAX)=U(LMKK)
C       ENDIF
   45 CONTINUE
C      CALL WRR(TAU123(1,NODMAX),3,'123 ')
C           TRANSFORMACIJA U LOKALNI SISTEM I UVODJENJE USLOVA KLIZANJA
C***
C*7
C   AKO JE PRETHODNI STATUS BIO KLIZANJE LOKALNO = GLOBALNO
C       WRITE(3,*)'************ IND',IND
C*290694       IF(IND.NE.2)THEN
C       IF(IND1.EQ.1)THEN
          CALL RSN123(TAURSN(1,NODMAX),TAU123(1,NODMAX),TRA,2)
C       ELSEIF(IND1.EQ.2)THEN
C*9
C        TAURSN(3,NODMAX)=TAU123(1,NODMAX)
C*** OSTALA DVA CLANA TREBA DA BUDU IZ PRETHODNOG PROLAZA !!!
C*9
C       ENDIF
C
C
C*7
C
C  PROVERA ZA OSLOBADJANJE IZ KONTAKTA
C
C      CALL WRR(TRA,9,'TRA ')
C      CALL WRR(TAURSN(1,NODMAX),3,'RSN ')
C      WRITE(3,*)'TAURSN(3,NMAX),IOVR,ICHEL',TAURSN(3,NODMAX),IOVR,ICHEL
C      IF(IOVR.NE.1.AND.TAURSN(3,NODMAX).LT.-SMALL)THEN
C      IF(ICHEL.EQ.0.AND.TAURSN(3,NODMAX).LT.-SMALL)THEN
C170794
      IF(TAURSN(3,NODMAX).LT.-SMALL.AND.PRODOR.GT.-TOLPD)THEN
        IND=0
        GO TO 20
      ENDIF
C
C           PROVERA TIPA KONTAKTA
C 
      IND=1
C... RRS = REZULTANTA U  RS RAVNI
      RRS=DSQRT(TAURSN(1,NODMAX)*TAURSN(1,NODMAX)+
     &          TAURSN(2,NODMAX)*TAURSN(2,NODMAX))
CC      IF(RRS.GE.FS*DABS(TAURSN(3,NODMAX))) IND=2
C070794      IF(RRS.GE.FS*DABS(TAURSN(3,NODMAX)).AND.DABS(FS).LT.BIG) IND=2
      RN=DABS(TAURSN(3,NODMAX))
C100794      IF(RN.GT.SMALL.AND.RRS.GE.FS*RN.AND.DABS(FS).LT.BIG) IND=2
      IF(ICHEL.EQ.0.AND.
     &   RN.GT.SMALL.AND.RRS.GE.FS*RN.AND.DABS(FS).LT.BIG) IND=2
      IF(DABS(FS).LT.SMALL) IND=2
C      CALL WRR(TAURSN(1,NODMAX),3,'*RSN')
C      WRITE(3,*)' IND',IND
C
C
C          B1)      KONTAKT BEZ KLIZANJA
C
C
      IF(IND.EQ.1) THEN
C   AKO JE PRETHODNA ITERACIJA  - SLIDING
C           TRANSFORMACIJA SILA U GLOBALNI SISTEM
       IF(IND1.EQ.2)CALL RSN123(TAURSN(1,NODMAX),TAU123(1,NODMAX),TRA,1)
C      CALL WRR(TAURSN(1,NODMAX),3,'#RSN ')
C      CALL WRR(TAU123(1,NODMAX),3,'#123 ')
C          KONTAKTNA KRUTOST 
       MDIM=MAXCE+2
      IF(IILS.NE.-1)THEN
       SKE(2) =-1.D0
       LLL=MDIM+1
       DO 55 L=1,MAXCE
       LLL=LLL+1
   55  SKE(LLL)= HE(L,1)
C      CALL WRR(HE(1,1),MAXCE,'HE  ')
C
       KK=-MAXCE-1
       DO 56 K=1,3
         KK=KK+MAXCE+2
         CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM(KK),MDIM)
   56  CONTINUE
      ENDIF
C           SILE NA CILJNIM CVOROVIMA
      KK=-MAXCE
      DO 52 K=1,3
       KK=KK+MAXCE+2
       IF(LM(KK).NE.0) THEN
        DUM=TAU123(K,NODMAX)
        SILE(K,NLM)=DUM
        DO 51 L=1,MAXCE
          NELL=NE+NELAB(L,NCA)
          DUM1=SILE(K,NELL)-HE(L,1)*DUM
          TAU123(K,L)=DUM1
C        WRITE(3,*)'K,NELL,SILE(K,NELL)',K,NELL,SILE(K,NELL)
C        WRITE(3,*)'K,L,TAU123(K,L)',K,L,TAU123(K,L)
          SILE(K,NELL)=DUM1
   51   CONTINUE
       ENDIF
   52 CONTINUE
C
      ENDIF
C
C
C             B2)      KONTAKT SA KLIZANJEM
C
C
      IF(IND.EQ.2) THEN
C          KONTAKTNA KRUTOST 
       MDIM=MAXCE+2
      IF(IILS.NE.-1)THEN
       KK=-MAXCE-1
       DO 258 K=1,3
         DUM=TRA(K,3)
         SKE(2) =-DUM
         LLL=MDIM+1
           DO 255 L=1,MAXCE
           LLL=LLL+1
  255      SKE(LLL)= HE(L,1)*DUM
         KK=KK+MAXCE+2
         KK1=KK+1
         LMDUM=LM(KK1)
         IF(LMDUM.NE.0) LM(KK1)=LMILM
C         CALL IWRR(LM(KK),MDIM,'LM**')
         CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM(KK),MDIM)
         LM(KK1)=LMDUM
  258  CONTINUE
       MDIM=2
       KK=0
       DO 260 M=1,3
       LM3(M)=0
         IF(IDC(NLM,M).NE.0.AND.IDC(NLM,M).NE.LMILM.AND.KK.LT.2)THEN
           KK=KK+1
           LM3(KK)=IDC(NLM,M)
         ENDIF
  260  CONTINUE
       SKE(1) =-1.D0
       SKE(2) = 0.D0
       SKE(3) =-1.D0
C         CALL IWRR(LM3,MDIM,'LM3*')
       CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM3,MDIM)
      ENDIF
C
C          KONTAKTNE SILE PRI KLIZANJU
       IF(RRS.GT.SMALL)THEN
         DUM=FD*DABS(TAURSN(3,NODMAX))/RRS
         DO 57 I=1,2
   57    TAURSN(I,NODMAX)=DUM*TAURSN(I,NODMAX)
       ENDIF
       CC=0.D0
       DO 96 K=1,3
   96  CC=CC+TRA(K,3)*DC(K)
       IF(DABS(CC).LT.SMALL)CC=0.D0
       DC(1)=CC
C           TRANSFORMACIJA SILA U GLOBALNI SISTEM
       CALL RSN123(TAURSN(1,NODMAX),TAU123(1,NODMAX),TRA,1)
C      CALL WRR(TAURSN(1,NODMAX),3,'@RSN ')
C      CALL WRR(TAU123(1,NODMAX),3,'@123 ')
CC      CALL WRR(TRA,9,'TRA ')
C           SILE NA CILJNIM CVOROVIMA
      KK=-MAXCE
      DO 65 K=1,3
       KK=KK+MAXCE+2
       IF(LM(KK).GT.0) THEN
        DUM=TAU123(K,NODMAX)
        SILE(K,NLM)=DUM
        DO 64 L=1,MAXCE
        NELL=NE+NELAB(L,NCA)
        DUM1=SILE(K,NELL)-HE(L,1)*DUM
        TAU123(K,L)=DUM1
        SILE(K,NELL)=DUM1
C        WRITE(3,*)'K,NELL,SILE(K,NELL)',K,NELL,SILE(K,NELL)
C        WRITE(3,*)'K,L,TAU123(K,L)',K,L,TAU123(K,L)
   64   CONTINUE
      ENDIF
   65 CONTINUE
      ENDIF
C........... END  IND=2
C
C  3)  RASPOREDJIVANJE KONTAKTNIH SILA
C
      ISLD=1
      KK=-MAXCE
      DO 80 K=1,3
       KK=KK+MAXCE+2
       KKLL=LM(KK)
       IF(KKLL.NE.0) THEN
        KKL=LM(KK-1)
        IF(KKL.NE.0) RC(KKL)=TAU123(K,NODMAX)
C        WRITE(3,*)'KKL,RC(KKL)',KKL,RC(KKL)
C100794
        IF(IND.EQ.2)   U(KKLL)=TAU123(K,NODMAX)
C.. OVO PONISTAVA REZULTAT PRETHODNIH ITERACIJA
C        IF(ICHEL.NE.0) U(KKLL)=0.D0
C100794  - ISKLJUCUJE OVO IZA !
C070794
C        IF(IND.EQ.1)THEN
C          U(KKLL)=TAU123(K,NODMAX)
C        ELSEIF(IND.EQ.2)THEN
C          IF(ISLD.EQ.1)THEN
C             U(KKLL)=TAURSN(3,NODMAX)
C             ISLD=0
C          ELSE
C             U(KKLL)=0.D0
C          ENDIF
C        ENDIF
C070794
C*8        IF(KKL.NE.0) RC(KKL)=RC(KKL)+TAU123(K,NODMAX)
        DO 70 L=1,MAXCE
          KKL=LM(KK+L)
C*8          IF(KKL.NE.0) RC(KKL)=RC(KKL)+TAU123(K,L)
          IF(KKL.NE.0) RC(KKL)=TAU123(K,L)
C        WRITE(3,*)'KKL,RC(KKL)',KKL,RC(KKL)
   70   CONTINUE
C           KONTAKTNO PRODIRANJE
       IF(IND.EQ.2.AND.KKLL.NE.LMILM) GOTO 80
        DUM=DC(K)
        IF(IND.EQ.2)DUM=DC(1)
        RC(KKLL)=DUM
C        WRITE(3,*)'KKLL,DC',KKLL,DUM
       ENDIF
   80 CONTINUE
C..  KONTAKTOR KOD GOVORI O TIPU KONTAKTA I PRVOM CVORU CILJNOG PODSEG.
      IK (NLM)=IND*IT1+NCAA
      IF(IK(NLM).NE.IK1(NLM)) ICCMOV=ICCMOV+1
      IK1(NLM)=IK(NLM)
      ALFK(1,NLM)=R0
      ALFK(2,NLM)=S0
      ALFK(3,NLM)=TRA(1,3)
      ALFK(4,NLM)=TRA(2,3)
      ALFK(5,NLM)=TRA(3,3)
CC      write(3,*)'NLM,IK1,R0,S0',NLM,IK1(NLM),R0,S0
      RETURN
      END
C=======================================================================
      SUBROUTINE ALFC3 (HE,TRA,XYZ,UPRI0,MAXCE,DC,R0,S0,EPSIL,PRODOR,
     &                  TOLPD,NN,TCMAX,NPL,NPOL,NCTC,NPOLMX,IND,IFLAG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     NALAZENJE POLOZAJA KONTAKTA
CE     FIND CONTACT LOCATION
C
C  IFLAG = 1   PRETRAZIVANJE PO POVRSINAMA (TEST A,B)
C  IFLAG = 2   PRETRAZIVANJE PO LINIJAMA (TEST C)
C  NN          CVOR POLIGONA NAJBLIZI KONTAKTORU
C  IND         INDIKATOR KONTAKTA (0-NEMA, 1-IMA)
C  DC()        KOMPONENTE PRODORA
C  TRA()       MATRICA TRANSFORMACIJE
C.........  ZA IFLAG=2   ....................................
C  NPOL        REDNI BROJ POLIGONA SA NAJVECIM (TESTC)
C  NPOLMX      UKUPAN BROJ POLIGONA NA SEGMENTU
C  NCTC        CVOR POLIGONA (NPOL) LINIJE NA KOJOJ JE PRODOR
C  TCMAX       MAKSIMALNA PROJEKCIJA
C
      DIMENSION HE(9,*),TRA(3,*),XYZ(3,*),DC(*),UPRI0(*)
      DIMENSION D(3),C1(3),C2(3),P(3),AN(3),C1P(3),PC2(3)
C      DATA RR/1.D0,-1.D0,-1.D0,1.D0/,SS/1.D0,1.D0,-1.D0,-1.D0/
C      DATA W0/0.D0/,W1/1.D0/
C
C      WRITE(3,*)'STIGO 1, IFLAG',IFLAG
      R0=0.D0
      S0=0.D0
      NODMAX=MAXCE+1
      DO 10 I=1,3
       D(I) =XYZ(I,NODMAX)-XYZ(I,NN)
       DC(I)=0.D0
   10 CONTINUE
C
C--- 2D CONTACT SURFACE
C
      IF(MAXCE.EQ.2)THEN
        MD=2
        IND=0
        DUM=0.D0
        DO 520 I=1,2
          C1(I)=XYZ(I,2)-XYZ(I,1)
          DUM=DUM+C1(I)*C1(I)
  520   CONTINUE        
        DUM=DSQRT(DUM)
C... PROBLEM DEPENDENT TOLERANCE
      TOLPD=EPSIL*DUM
        DO 525 I=1,2
  525   C1(I)=C1(I)/DUM
        R0=(C1(1)*D(1)+C1(2)*D(2))/DUM
C100794
C        WRITE(3,*)'****NN,R0',NN,R0
C... POKUSAJ NA DRUGI NACIN
        IF(DABS(R0).GT.1.D08)
     &    R0=(UPRI0(2)*D(1)-UPRI0(1)*D(2))/
     &       (UPRI0(2)*C1(1)-UPRI0(1)*C1(2))/DUM
C100794
        IF(NN.EQ.2) R0=1.+R0
C        WRITE(3,*)'>>>>NN,R0',NN,R0
C  PEGLANJE
        IF(DABS(R0).LT.EPSIL)    R0=0.D0
        IF(DABS(R0-1.).LT.EPSIL) R0=1.D0
C*1        IF(R0.LT.-EPSIL.OR.R0.GT.1.+EPSIL) RETURN
        IF(R0.LT.0.D0.OR.R0.GT.1.D0) RETURN
        IND=1
        CALL TRAHE(XYZ,HE,TRA,R0,S0,MAXCE)
C        HE(1,1)=1.-R0
C        HE(2,1)=R0
C        TRA(1,1)= W0
C        TRA(2,1)= W0
C        TRA(3,1)= W1
C        TRA(1,2)= C1(1)
C        TRA(2,2)= C1(2)
C        TRA(3,2)= W0
C        TRA(1,3)=-C1(2)
C        TRA(2,3)= C1(1)
C        TRA(3,3)= W0
CC      CALL WRR(TRA,9,'TRA ')
      ELSE
C
C--- 3D CONTACT SURFACE
C
      MD=3
      IND=1
C*      IF(IFLAG.EQ.2)THEN
C*        R0=RR(NN)+TCMAX*(RR(NCTC)-RR(NN))
C*        S0=SS(NN)+TCMAX*(SS(NCTC)-SS(NN))
C*        IFLAG=1
C*      ELSE
C
C  0) ... PROBLEM DEPENDENT TOLERANCE   (TOLPD)
C
      N1=NN+1
      N2=NN-1
      IF(NN.EQ.MAXCE) N1=1
      IF(NN.EQ.1) N2=MAXCE
      DO 20 I=1,3
       C1(I)=XYZ(I,N1)-XYZ(I,NN)
       C2(I)=XYZ(I,N2)-XYZ(I,NN)
   20 CONTINUE
      CINT2=AOBS(C1,C1)
C
      TOLPD=EPSIL*DSQRT(CINT2)
C
C  1) CVOR KONTAKTOR I CILJNI CVOR SE POKLAPAJU
C
      R0=0.D0
      S0=0.D0
      DI=DSQRT(AOBS(D,D))
C      WRITE(3,*)'DI,EPSIL,IND',DI,EPSIL,IND
      IF(DI.LT.EPSIL)THEN
        IND=-1
C        DO 15 I=1,3
C   15   XYZ(I,NODMAX)=XYZ(I,NN)
        GO TO 50
      ENDIF
C
C  2) CVOR KONTAKTOR I CILJNI CVOR SE NE POKLAPAJU
C
C
C   TEST  (C) :  D * C0    ,    C0 = CK / I CK I
C
      TESTC = AOBS( D, C1 )/CINT2
      IF(TESTC.GT.TCMAX)THEN
        TCMAX=TESTC
        NPOL =NPL
        NCTC =N1
      ENDIF
      CINT2=AOBS(C2,C2)
      TESTC = AOBS( D, C2 )/CINT2
      IF(TESTC.GT.TCMAX)THEN
        TCMAX=TESTC
        NPOL =NPL
        NCTC =N2
      ENDIF
C      WRITE(3,*)'TESTC',TESTC
C
C   NORMALA (MALO N) AN = (C1 X C2) / I C1 X C2 I
C 
      CALL AXBV( C1, C2, AN )
      CALL JEDV( AN(1), AN(2), AN(3) )
C
C   VEKTOR (MALO A)     =  - (AN * D) AN   SMESTA SE U  P
C			     
      AND = -AOBS(AN,D)      
      CALL JEDNAK(P,AN,AND,3)
C
C   VEKTOR (MALO P)  P = D + MALO A
C
      CALL ZBIRM1(P,D,3)
C
C   TEST  (A) :  (C1 X P)*AN > 0   ,   C1P = C1 X P 
C
      CALL AXBV( C1, P, C1P )
      TESTA = AOBS( AN, C1P )
C      WRITE(3,*)'TESTA',TESTA
      IF(TESTA.LT.-TOLPD)THEN
        IND=0
        RETURN
      ENDIF
C
C   TEST  (B) :  (C1 X P)*(P X C2) > 0   ,   PC2 = P X C2
C
      CALL AXBV( P, C2, PC2 )
      TESTB = AOBS( C1P, PC2 )
C      WRITE(3,*)'TESTB',TESTB
      IF(TESTB.LT.-TOLPD)THEN
        IND=0
        RETURN
      ENDIF
C
C   PROVERA DA LI PRECI NA IFLAG=2
C
      IF(DABS(TESTA).LE.TOLPD.OR.DABS(TESTB).LE.TOLPD) IND=0
      IF(IND.EQ.0)THEN
        IND=1
C        R0=RR(NN)+TCMAX*(RR(NCTC)-RR(NN))
C        S0=SS(NN)+TCMAX*(SS(NCTC)-SS(NN))
C      WRITE(3,*)'**R0,S0,IND',R0,S0,IND
C        GO TO 55
        GO TO 50
C*        IF(NPOL.EQ.NPOLMX)IFLAG=2
C*        RETURN
      ENDIF
C
C  3) NALAZENJE TACNIH KOORDINATA (R,S) PRODORA
C
   50 CALL PRODR3(XYZ,R0,S0,MAXCE,IND)
      IF(IND.EQ.-1)THEN
        IND=IABS(IND)
        R0=R0/DABS(R0)
        S0=S0/DABS(S0)
      ENDIF
C      WRITE(3,*)'R0,S0,IND',R0,S0,IND
      IF(IND.EQ.0)RETURN
C   55 CONTINUE
C
C  ... KRAJ ZA IFLAG=1
C
C*      ENDIF
C
C  4) MATRICA TRANSFORMACIJE  TRA() (LOKALNI <-> GLOBALNI SISTEM)
C
      CALL TRAHE(XYZ,HE,TRA,R0,S0,MAXCE)
C
C  ... KRAJ ZA MAXCE=2
C
      ENDIF
C
C  5) VEKTOR     (KONTAKTOR_CVOR    CILJNA_PRODORNA_TACKA_P)
C
      IF (IND.EQ.1)THEN
       DO 210 J=1,MD
        DO 200 I=1,MAXCE
  200   DC(J)=DC(J)+HE(I,1)*XYZ(J,I)
       DC(J)=-DC(J)+XYZ(J,NODMAX)
C290694       IF(DABS(DC(J)).LT.1.D-10) DC(J)=0.D0
  210  CONTINUE
C      WRITE(3,*)'DC',(DC(K),K=1,3)
C
C  6) PROVERA DA LI SE RADI O PRODORU ILI ZAZORU
C
      PRODOR = AOBS( DC, TRA(1,3) )
      IF(PRODOR.GT.TOLPD) IND=0
C      WRITE(3,*)'IND,PRODOR,TOLPD',IND,PRODOR,TOLPD
      ENDIF
      RETURN
      END      
C=======================================================================
      SUBROUTINE TRAHE(XYZ,HE,TRA,R,S,MAXCE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C  INTERPOLACIONE FUNKCIJE I MATRICA TRANSFORMACIJE
C
      DIMENSION XYZ(3,*),HE(9,*),TRA(3,*)
      DATA C0/0.D0/,C1/1.D0/,C25/.25D0/
C
C 1)  INTERPOLACIJE
C
      IF(MAXCE.EQ.4)THEN
      RP=C1+R
      RM=C1-R
      RK=C1-R*R
      SP=C1+S
      SM=C1-S
      SK=C1-S*S
      HE(1,1)=C25*RP*SP
      HE(2,1)=C25*RM*SP
      HE(3,1)=C25*RM*SM
      HE(4,1)=C25*RP*SM
      HE(1,2)=C25*SP
      HE(2,2)=-C25*SP
      HE(3,2)=-C25*SM
      HE(4,2)=C25*SM
      HE(1,3)=C25*RP
      HE(2,3)=C25*RM
      HE(3,3)=-C25*RM
      HE(4,3)=-C25*RP
C
      ELSEIF(MAXCE.EQ.3)THEN
      HE(1,1)=C1-R-S
      HE(2,1)=R
      HE(3,1)=S
      HE(1,2)=-C1
      HE(2,2)= C1
      HE(3,2)= C0
      HE(1,3)=-C1
      HE(2,3)= C0
      HE(3,3)= C1
      ELSEIF(MAXCE.EQ.2)THEN
        HE(1,1)=C1-R
        HE(2,1)=R
      ENDIF
C
C 2)  BAZNI VEKTORI
C
      IF(MAXCE.GT.2)THEN
      DO 20 I=1,3
      DO 20 J=1,2
      TRA(I,J)=C0
       DO 10 K=1,MAXCE
        TRA(I,J)=TRA(I,J)+HE(K,J+1)*XYZ(I,K)
   10  CONTINUE
   20 CONTINUE
C... NORMALA  V3
        CALL AXBV( TRA(1,1), TRA(1,2), TRA(1,3) )
        CALL JEDV( TRA(1,1), TRA(2,1), TRA(3,1) )
        CALL JEDV( TRA(1,3), TRA(2,3), TRA(3,3) )
C... DRUGI JEDINICNI VEKTOR PRAVOUGLOG SISTEMA  V2
        CALL AXBV( TRA(1,3), TRA(1,1), TRA(1,2) )
      ELSE
        TRA(1,1)= C0
        TRA(2,1)= C0
        TRA(3,1)= C1
        DO 50 I=1,2
   50   TRA(I,2)=XYZ(I,2)-XYZ(I,1)
        TRA(3,2)= C0
        CALL JEDV( TRA(1,2), TRA(2,2), TRA(3,2) )
        TRA(1,3)=-TRA(2,2)
        TRA(2,3)= TRA(1,2)
        TRA(3,3)= C0
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE PRODR3(XYZ,R0,S0,MAXCE,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C  NALAZENJE TACNIH KOORDINATA (R,S) PRODORA
C
      DIMENSION XYZ(3,*),PAR(8),DES(2)
      EQUIVALENCE (A1,PAR(1)),(A2,PAR(2)),(C1,PAR(3)),(C2,PAR(4)),
     &            (B, PAR(5)),(D, PAR(6)),(E, PAR(7)),(F, PAR(8))
      DATA C25/0.25D0/,EPSIL/1.D-8/,MAXIT/10/,SMALL/1.D-16/
      DO 5 I=1,8
       PAR(I)=0.D0
    5 CONTINUE
      IF(MAXCE.EQ.4)THEN
      DO 10 I=1,3
       AI=C25*(XYZ(I,1)+XYZ(I,2)+XYZ(I,3)+XYZ(I,4))
       BI=C25*(XYZ(I,1)-XYZ(I,2)-XYZ(I,3)+XYZ(I,4))
       CI=C25*(XYZ(I,1)+XYZ(I,2)-XYZ(I,3)-XYZ(I,4))
       DI=C25*(XYZ(I,1)-XYZ(I,2)+XYZ(I,3)-XYZ(I,4))
       DUM=XYZ(I,5)-AI
       A1 =A1+BI*DUM
       A2 =A2+CI*DUM
       C1 =C1-BI*BI
       C2 =C2-CI*CI
       B  =B +DI*DUM-BI*CI
       D  =D -BI*DI
       E  =E -CI*DI
       F  =F -DI*DI
   10 CONTINUE
C
C  I T E R A C I J E     FULL-NEWTON
C
      IT=0
      DES(1)=R0
      DES(2)=S0
   30 IT=IT+1
      R=DES(1)
      S=DES(2)
      R2=R*R
      S2=S*S
      RS=R*S
C
      FI11=C1+2.*D*S+F*S*S
      FI12=B+2.*(E*S+D*R+F*R*S)
      FI22=C2+2.*E*R+F*R*R
      DET=FI11*FI22-FI12*FI12
      IF(DABS(DET).LT.SMALL)RETURN
C
      DES(1) =-A1-C1*R-B *S-E*S2-2.*D*RS-F*R*S2
      DES(2) =-A2-B *R-C2*S-D*R2-2.*E*RS-F*R2*S
C
      DET1=DES(1)*FI22-FI12*DES(2)
      DET2=FI11*DES(2)-DES(1)*FI12
      DES(1)=DET1/DET
      DES(2)=DET2/DET
C
      ANOR=DES(1)*DES(1)+DES(2)*DES(2)
      DES(1)=R+DES(1)
      DES(2)=S+DES(2)
        IF(ANOR.LT.EPSIL.OR.IT.EQ.MAXIT) THEN
          R0=DES(1)
          S0=DES(2)
          IF((DABS(R0).GT.1.D0.OR.DABS(S0).GT.1.D0).AND.IND.GT.0) IND=0
          RETURN
        ENDIF
      GO TO 30
C
      ELSEIF(MAXCE.EQ.3)THEN
      DO 50 I=1,3
       AI= XYZ(I,1)
       BI=-XYZ(I,1)+XYZ(I,2)
       CI=-XYZ(I,1)+XYZ(I,3)
       DUM=XYZ(I,4)-AI
       A1 =A1+BI*DUM
       A2 =A2+CI*DUM
       C1 =C1-BI*BI
       C2 =C2-CI*CI
       B  =B -BI*CI
   50 CONTINUE
      DET=C1*C2-B*B
      IF(DABS(DET).LT.SMALL)RETURN
      R0=(B*A2-A1*C2)/DET
      S0=(A1*B-C1*A2)/DET
      IF((DABS(R0).GT.1.D0.OR.DABS(S0).GT.1.D0).AND.IND.GT.0) IND=0
      RETURN
C
      ENDIF
      END
C=======================================================================
      SUBROUTINE UPXYZ(XYZ,NCK,NELSF,CORD,U,ID,NP,MAXCE,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C   DEFORMISANA KONFIGURACIJA U ITERACIJI
C
      DIMENSION CORD(NP,*),U(*),ID(NP,*),NELSF(*),XYZ(3,*)
C
      GOTO (1,2)IND
    1 DO 20 N=1,MAXCE
        NN=NELSF(N)
        IF(NN.EQ.0) GO TO 20
        DO 10 L=1,3
         I=ID(NN,L)
         IF(I.GT.0)THEN
           XYZ(L,N)=CORD(NN,L)+U(I)
         ELSE
           XYZ(L,N)=CORD(NN,L)
         ENDIF
   10   CONTINUE
   20 CONTINUE
      RETURN
C
    2 NODMAX=MAXCE+1
      DO 30 L=1,3
       I=ID(NCK,L)
       IF(I.GT.0)THEN
         XYZ(L,NODMAX)=CORD(NCK,L)+U(I)
       ELSE
         XYZ(L,NODMAX)=CORD(NCK,L)
       ENDIF
   30 CONTINUE
      RETURN
      END
C=======================================================================
      FUNCTION AOBS(A,B)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C  SKALARNI PROIZVOD
C
      DIMENSION A(*),B(*)
      AOBS=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      RETURN
      END
C=======================================================================
      SUBROUTINE SILE93(IDC,U,IK,ALFK,NE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     KOREKCIJA PRIRASTAJA SILA
C  AKO JE BILO KLIZANJE U PRETHODNOJ ITERACIJI 
C  SILU TRANSFORMISE U GLOBALNI SISTEM KAO DA NEMA KLIZANJA
CE     UPDATE FORCE INCREMENT
C
      DIMENSION U(*),IDC(NE,*),IK(*),ALFK(5,*)
      DATA IT1/1000000/
      DO 100 NLM=1,NE
C      CALL WRR(ALFK(1,NLM),5,'ALFK')
C... SAMO ZA CVOROVE KOJI KLIZAJU
        IND =IK(NLM)/IT1
        IF(IND.NE.2)GO TO 100
        LM=0
        DO 10 K=1,3
   10   IF(LM.EQ.0) LM=IDC(NLM,K)
        DF=U(LM)
        DO 20 K=1,3
          LM=IDC(NLM,K)
C          WRITE(3,*)'*LM',LM
          IF(LM.EQ.0)GO TO 20
          K2=2+K
          U(LM)=DF*ALFK(K2,NLM)
   20 CONTINUE
  100 CONTINUE
      RETURN
      END
C=======================================================================
C070794
      SUBROUTINE RSN123(TAURSN,TAU123,TRA,ICAL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C   TRANSFORMACIJA SILA   RSN   NA  123
C   ICAL = 1   -   RSN ---> 123
C   ICAL = 2   -   123 ---> RSN 
C
      DIMENSION TAURSN(*),TAU123(*),TRA(3,*)        
      GO TO (1,2), ICAL
1     DO 20 J=1,3
        CC=0.D0
        DO 15 K=1,3
   15   CC=CC+TRA(J,K)*TAURSN(K)
        TAU123(J)=CC
   20 CONTINUE
      RETURN
C
2     DO 30 J=1,3
        CC=0.D0
        DO 25 K=1,3
   25   CC=CC+TRA(K,J)*TAU123(K)
        TAURSN(J)=CC
   30 CONTINUE
      RETURN
      END
