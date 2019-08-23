C====================================================================
C
C  UPRAVLJACKI PROGRAMI ZA POZIVANJE PROGRAMA ZA RACUNANJE
C  MATRICA ELEMENATA
C
C  SUBROUTUNE SO4EGL
C             KIGR1
C             INTKG
C             MELIG
C             JACKG
C=====================================================================
      SUBROUTINE S04EGL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'paka.inc'
      
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /IGRGLA/ NGRU,NGSU,NGTU,NGS,NTE,NUD,IANAL,JTAU,NWE,NT6,INEL
     *,MATG
      COMMON /IGREP/  LIGNOP,LIGNDS,LIGDST,LIGCTR,MXAUG,LMXAUG,LAG,MAXAE
     *               ,LICTRP,LIGX  ,LIGY  ,LIGZ,LAE,LLM,LIGMAT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, 'S04EGL'
CC    DIMENSION AG(*),AE(*)
CC    REAL AG,AE
C***************
C               RACUNANJE MATRICE KRUTOSTI I UNUTRASNJIH SILA
C
C***************
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX = LMAX + 1
      LAG = LMAX
C
      CALL READEG(A(LAG))
C
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX = LMAX + 1
      LAE=LMAX
      MAXAE=(18+33*NT6+2*NTE+NT6*(NT6+1)/2)*IDVA
      LMAX = LAE + MAXAE
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX = LMAX + 1
C  UCITAVANJE REPERA ZA MATERIJALE
C
      CALL UCITAM(A(LMODEL),NMODM)
C
      IF(LMAX.LT.MTOT) GO TO 70
      WRITE(IZLAZ,2009) LMAX,MTOT
 2009 FORMAT(' ',' NEDOVOLJNA DIMENZIJA U OSNOVNOM RADNOM VEKTORU A ',
     *'ZA UCITAVANJE ULAZNIH PODATAKA'/' ','POTREBNA DIMENZIJA ',
     *'LMAX = ',I10/' ','RASPOLOZIVA DIMENZIJA, MTOT=',I10)
      STOP
   70 CALL KIGR1(A(LAG),A(LAE))
      RETURN
      END
      SUBROUTINE KIGR1(AG,AE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C POTPROGRAM ZA POZIVANJE POTPROGRAMA ZA INTEGRALJENJE MATRICE
C KRUTOSTI ELEMENATA
C
      include 'paka.inc'
C
      
C
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
C
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
C
      COMMON /DUPLAP/ IDVA
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
C
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C
      COMMON /ZAPISI/ LSTAZA(5)
C
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
C
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
C
      COMMON /ITERBR/ ITER
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
C
      COMMON /MATERM/ LMODEL,LGUSM
C
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
C
      COMMON /IGRGLA/ NGRU,NGSU,NGTU,NGS,NTE,NUD,IANAL,JTAU,NWE,NT6,INEL
     *,MATG
C
      COMMON /IGREP/  LIGNOP,LIGNDS,LIGDST,LIGCTR,MXAUG,LMXAUG,LAG,MAXAE
     *               ,LICTRP,LIGX  ,LIGY  ,LIGZ,LAE,LLM,LIGMAT
      COMMON /ELEMA4/ LTSG,LBL,LBLT,LBNL,LBNLT,LH,LHI,LSKE,LMXAE
C
      DIMENSION AG(*),AE(*)
      REAL AG,AE
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, 'KIGR1'
C
      LTSG= 1
      LBL = LTSG+18*IDVA
      LBLT= LBL+3*NT6*IDVA
      LBNL= LBLT+3*NT6*IDVA
      LBNLT=LBNL+9*NT6*IDVA
      LH  = LBNLT+NT6*9*IDVA
      LHI=  LH+NTE*2*IDVA
      LSKE= LHI+9*NT6*IDVA
      LMXAE=LSKE + NT6*(NT6+1)/2*IDVA
C
      IF(ISKNP.EQ.1) KORK=1
      IF(ISKNP.EQ.2) KORK=2
C     IF(IND  .EQ.4) KORK=3
C
C INEGRACIJA MATRICE KRUTOSTI KRIVE GREDE
C
C
C
C  50 IF(ISKNP.EQ.1) GO TO 30
C
      CALL INTKG(A(LID),A(LCORD),AG(LIGNOP),AG(LIGNDS),AG(LIGMAT),
     *           AG(LIGDST),AG(LIGCTR),A(LMAXA),A(LRTDT),
     *           A(LFTDT),AG(LLM),
     *           AE(LBL),AE(LBLT),AE(LBNL),AE(LBNLT),
     1           AE(LH),AE(LSKE),AE(LTSG),KORK)
C
C     IF(IATYP.EQ.0.OR.IATYP.EQ.1) GO TO 30
      IF(IATYP.EQ.0) GO TO 30
      IF(NMODM.LT.5) GO TO 30
      NPROS=NE*NGS*(JTAU+10)*IDVA
      LMAX8=LSTAZA(4)-1
      CALL WRITDD(A(LPLAS1),NPROS/IDVA,IELEM,LMAX8,LDUZI)
   30 CONTINUE
C
C   30 IF(ISKNP.EQ.2) GO TO 40
C     CALL INTKG(A(LID),A(LCORD),AG(LIGNOP),AG(LIGNDS),AG(LIGMAT),
C    *           AG(LIGDST),AG(LIGCTR),A(LMAXA),A(LRTDT),
C    *           A (LFTDT),AG(LLM),
C    *           AE(LBL),AE(LBLT),AE(LBNL),AE(LBNLT),
C    1           AE(LH),AE(LSKE),AE(LTSG),1)
C
C ZAPISIVANJE NAPONA NA DISK IELEM
C
C     LMAX8=LSTAZA(2)
      LMAX8=LSTAZA(2)-1
CCC   CALL WRITDD(AG(LIGCTR),MXAUG/IDVA,IELEM,LMAX8,LDUZI)
C     LMAX8=LSTAZA(2)
      LMAX8=LSTAZA(2)-1
CCC   CALL READDD(AG(LIGCTR),MXAUG/IDVA,IELEM,LMAX8,LDUZI)
C
      IF(NMODM.GE.5) RETURN
      NPROS=NE*NGS*JTAU*IDVA
      LMAX8=LSTAZA(5)-1
      CALL WRITDD(A(LSIGMA),NPROS/IDVA,IELEM,LMAX8,LDUZI)
C
      RETURN
      END
      SUBROUTINE INTKG(ID,CORD,NOP,NDST,MATIG,DST,CTR,MAXA,UKUM,
     1FNEUR,LM,BL,BLT,BNL,BNLT,H,SKE,TSG,KORK)
      USE MATRICA
      USE STIFFNESS
      USE DRAKCE8
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'paka.inc'
      
C
C***********************************************************
C   PODPROGRAM ZA INTEGRACIJU MATRICE KRUTOSTI             *
C   ELEMENTA KRIVE GREDE                                   *
C***********************************************************
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /DUPLAP/ IDVA
      COMMON/ELKG/DET,WTU,R,S,T,NLM,JNGS,NCV,ND,PI2,SS,GR(3),GS(3)
     1,XJ(3,3),AK(40),BK(40),DLS(40),DLT(40),VS(40,3),VT(40,3),V1(40,3)
      COMMON /REPERM/ MREPER(4)
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ITERBR/ ITER
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /IGRGLA/ NGRU,NGSU,NGTU,NGS,NTE,NUD,IANAL,JTAU,NWE,NT6,INEL
     *,MATG
      COMMON /IGREP/  LIGNOP,LIGNDS,LIGDST,LIGCTR,MXAUG,LMXAUG,LAG,MAXAE
     *               ,LICTRP,LIGX  ,LIGY  ,LIGZ,LAE,LLM,LIGMAT
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /IGMAT / EM(3,3),IBTC,LNAP1,LPLAS,LPLA1
      COMMON /ELEMA4/ LTSG,LBL,LBLT,LBNL,LBNLT,LH,LHI,LSKE,LMXAE
      COMMON/CONSKG/C43,C23,C13,C15,C19
      COMMON /CDEBUG/ IDEBUG
      DIMENSION ID(NP,*),NOP(NE,*),MAXA(*),NDST
     1(NE,*),DST(NUD,*),CTR(NE,6,*)
      DIMENSION XG(55),WGT(55),FNEUR(*),
     1LM(NT6,*),MJA(3),NREF(11),LJA(3)
      DIMENSION CORD(NP,*),UKUM(*),MATIG(*)
C     DIMENSION AE(*),AG(1)
      DIMENSION BL(3,*),BLT(NT6,*),BNL(9,*),BNLT(NT6,*)
     1         ,H(NTE,*),SKE(*),TSG(6,3)
      DIMENSION FL(24)
C
      DATA NREF/0,1,3,6,10,15,21,28,36,45,55/
      DATA WGT/            2.D0,               1.D0,               1.D0,
     1       .555555555555556D0, .888888888888889D0, .555555555555556D0,
     2       .347854845137454D0, .652145154862546D0, .652145154862546D0,
     3       .347854845137454D0, .236926885056189D0, .478628670499366D0,
     4       .568888888888889D0, .478628670499366D0, .236926885056189D0,
     5       .171324492379170D0, .360761573048139D0, .467913934572691D0,
     6       .467913934572691D0, .360761573048139D0, .171324492379170D0,
     7       .129484966168870D0, .279705391489277D0, .381830050505119D0,
     8       .417959183673469D0, .381830050505119D0, .279705391489277D0,
     9       .129484966168870D0, .101228536290376D0, .222381034453374D0,
     9       .313706645877887D0, .362683783378362D0, .362683783378362D0,
     1       .313706645877887D0, .222381034453374D0, .101228536290376D0,
     2       .081274388361574D0, .180648160694857D0, .260610696402935D0,
     3       .312347077040003D0, .330239355001260D0, .312347077040003D0,
     4       .260610696402935D0, .180648160694857D0, .081274388361574D0,
     5       .066671344308688D0, .149451349150581D0, .219086362515982D0,
     6       .269266719309996D0, .295524224714753D0, .295524224714753D0,
     7       .269266719309996D0, .219086362515982D0, .149451349150581D0,
     8       .066671344308688D0/
      DATA XG /            0.D0,-.577350269189626D0, .577350269189626D0,
     1      -.774596669241483D0,               0.D0, .774596669241483D0,
     2      -.861136311594053D0,-.339981043584856D0, .339981043584856D0,
     3       .861136311594053D0,-.906179845938664D0,-.538469310105683D0,
     4                     0.D0, .538469310105683D0, .906179845938664D0,
     5      -.932469514203152D0,-.661209386466265D0,-.238619186083197D0,
     6       .238619186083197D0, .661209386466265D0, .932469514203152D0,
     7      -.949107912342759D0,-.741531185599394D0,-.405845151377397D0,
     8                     0.D0, .405845151377397D0, .741531185599394D0,
     9       .949107912342759D0,-.960289856497536D0,-.796666477413627D0,
     9      -.525532409916329D0,-.183434642495650D0, .183434642495650D0,
     1       .525532409916329D0, .796666477413627D0, .960289856497536D0,
     2      -.968160239507626D0,-.836031107326636D0,-.613371432700590D0,
     3      -.324253423403809D0,               0.D0, .324253423403809D0,
     4       .613371432700590D0, .836031107326636D0, .968160239507626D0,
     5      -.973906528517172D0,-.865063366688985D0,-.679409568299024D0,
     6      -.433395394129247D0,-.148874338981631D0, .148874338981631D0,
     7       .433395394129247D0, .679409568299024D0, .865063366688985D0,
     8       .973906528517172D0/
C
      IF(IDEBUG.GT.0) PRINT *, 'INTKG'
C
      LMHT1=LMHT-1
C
      IF(IATYP.EQ.0) GO TO 217
C
C KORIGOVANJE POCETNIH KOORDINATA SA UKUPNIM POMERANJIMA
C
C     DO 44 I=1,NP
C     DO 44 J=1,3
C     IF(ID(I,J).EQ.0) GO TO 44
C     K=ID(I,J)
C     GO TO (41,42,43)J
C  41 X(I)=CORD(I,1)+UKUM(K)
C     GO TO 44
C  42 Y(I)=CORD(I,2)+UKUM(K)
C     GO TO 44
C  43 Z(I)=CORD(I,3)+UKUM(K)
C  44 CONTINUE
      GO TO 208
C
  217 CONTINUE
C 217 DO 45 I=1,NP
C     DO 45 J=1,3
C     IF(ID(I,J).EQ.0) GO TO 45
C     K=ID(I,J)
C     GO TO (46,47,48)J
C  46 X(I)=CORD(I,1)
C     GO TO 45
C  47 Y(I)=CORD(I,2)
C     GO TO 45
C  48 Z(I)=CORD(I,3)
C  45 CONTINUE
C***
  208 CONTINUE
C
C
C  IZJEDNACAVANJE SA NULOM POCETNIH NAPONA
      LH1=LH-1
C***
C     PETLJA PO ELEMENTIMA
C***
      DO 10 NLM=1,NE
       IF(IATYP.EQ.0) GO TO 987
          DO 660 I=1,NT6
  660     FL(I)=0.0D0
987   CONTINUE
C-----
      NWE=NT6*(NT6+1)/2
C-----
      IF(ISKNP.EQ.2) GO TO 662
      DO 661 I=1,NWE
  661 SKE(I)=0.0D0
  662 CONTINUE
C
C PODACI O MATERIJALU
      MATG=MATIG(NLM)
      GO TO (50,50,50,50,50,50,50,50),NMODM
   50 LFUN=MREPER(1)
      CALL MELIG(A(LFUN))
C
      IF(IATYP.EQ.0.OR.ITER.EQ.0.OR.KORK.EQ.1) GO TO 209
      CALL KONSTE(NOP,ID,CORD,FNEUR,UKUM,NDST,CTR,DST,LM,0)
C
  209 IF(KORK.EQ.2.OR.KORK.EQ.3) GO TO 29
      NWE1=LSKE+NWE*IDVA-1
   29 JNGS=0
      CALL KONSTE(NOP,ID,CORD,FNEUR,UKUM,NDST,CTR,DST,LM,1)
C***
C     PETLJE PO GAUSOVIM TACKAMA
C***
      ND= 6*NCV
CCC
      DO 20 NGR=1,NGRU
      IF(NGSU.EQ.1)S=0.D0
      JGR=NREF(NGRU)+NGR
      R=XG(JGR)
      WR=WGT(JGR)
      DO 20 NGSS=1,NGSU
      IF(NGSU.GT.1) GO TO 37
      S=0.D0
      WS=1.D0
      GO TO 38
   37 JGR=NREF(NGSU)+NGSS
      S=XG(JGR)
      WS=WGT(JGR)
   38 DO 20 NGT=1,NGTU
      IF(NGTU.GT.1) GO TO 36
      T=0.D0
      WT=1.D0
      GO TO 27
   36 JGR=NREF(NGTU)+NGT
      T=XG(JGR)
      WT=WGT(JGR)
C***
C   ODREDJIVANJE MATRICE JAKOBIJANA
C***
   27 CONTINUE
      CALL JACKG(CORD,NOP,H)
C
      IF (DET.GT.1.D-10) GO TO 7
      WRITE(3,1333)NLM,NGR,NGSS,NGT
 1333 FORMAT(' ','S T O P   DET.LE.1.D-10'/' ','NLM =',I2/' ','NGR =',I2
     1/' ','NGSS=',I2/' ','NGT =',I2)
      STOP
    7 DO 9 I=1,3
      GR(I)=XJ(1,I)
    9 GS(I)=XJ(2,I)
C     IF(KORK.EQ.3) GO TO 30
C***
C     MATRICE ELEMENTA U GAUS.TA$KI
C***
      WTU=WR*WS*WT
      CALL MINV(XJ,3,DET,LJA,MJA)
      JNGS=JNGS+1
C     IF(NMODM.LT.5) GO TO 3211
C     LPLA1=LPLAST+(NLM-1)*NGS*(JTAU+10)*IDVA
C    *            +(JNGS-1)*(JTAU+10)*IDVA
C
      CALL MATEKG(ID,NOP,BL,BLT,BNL,BNLT,
     1  H,SKE,TSG,LM,UKUM,FNEUR,KORK)
C
   20 CONTINUE
C
C     TRANSFORMACIJA SILA NA GLOBALNI SISTEM
C
C     DO 600 I=1,NT
C     NN=NOP(NLM,I)
C     DO 610 K=1,3
C     TE(1,K) = V1(NN,K)
C     TE(1,K) = VS(NN,K)
C 610 TE(1,K) = VT(NN,K)
C     KL=0
C     II=(I-1)*6
C     DO 620 J=1,3
C     DO 630 K=1,3
C     FG(II+J ) = FG(II+J) + TE(K,J)*FL(II+K)
C     FG(II+J+3) = FG(II+J+3) + TE(K,J)*FL(II+K+3)
C 630 CONTINUE
C 620 CONTINUE
C 600 CONTINUE
CCC   DO 640 I=1,NT6
CCC   II=LM(I,NLM)
CCC   IF(II.EQ.0) GO TO 640
CCC   FNEUR(II)=FNEUR(II)+FL(I)
CC640 CONTINUE
C
C     IF(KORK .EQ.2.OR.KORK.EQ.3) GO TO 10
C***
C     RAZMESTANJE MATRICA I VEKTORA ELEMENATA U SISTEM
C
      IF(ISKNP.EQ.2) GO TO 10
      IF (TIPTACKANJA.EQ.1) THEN
      CALL SPAKUJ(ALSK,MAXA,SKE,LM(1,NLM),6*NTE)
      ELSE
      CALL SPAKUJMT(ALSK,MAXA,SKE,LM(1,NLM),6*NTE)
      ENDIF
C
C     GO TO 10
C***
C     KOREKCIJA KOORDINATA TACAKA ZA DEF. GLAV. PRAVACA PRESEKA
C***
C  22 CALL KONSTE(NOP,ID,X,Y,Z,CORD,FNEUR,NDST,CTR,DST,AE(LLM),0)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE MELIG(FUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C FORMIRANJE MATRICE ELEASTICNOSTI
C
      COMMON /IGRGLA/ NGRU,NGSU,NGTU,NGS,NTE,NUD,IANAL,JTAU,NWE,NT6,INEL
     *,MATG
      COMMON /IGMAT / EM(3,3),IBTC,LNAP1,LPLAS,LPLA1
C
      DIMENSION FUN(2,1)
C
      MATG=1
      EE =FUN(1,MATG)
      ANI=FUN(2,MATG)
C***
C     PODACI O MATERIJALIMA
C***
      DO 26 I=1,3
      DO 26 J=1,3
   26 EM(I,J) =0.D0
      EM(1,1) =EE
      EM(2,2) =EE/(2.0D0*(1.0D0 + ANI))
      EM(3,3) =EM(2,2)
      RETURN
      END
C     SUBROUTINE JACKG(CORD,NOP,X,Y,Z,H)
      SUBROUTINE JACKG(CORD,NOP,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ELKG/DET,WTU,R,S,T,NLM,JNGS,NCV,ND,PI2,SS,GR(3),GS(3)
     1,XJ(3,3),AK(40),BK(40),DLS(40),DLT(40),VS(40,3),VT(40,3),V1(40,3)
C
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
C
      COMMON /IGRGLA/ NGRU,NGSU,NGTU,NGS,NTE,NUD,IANAL,JTAU,NWE,NT6,INEL
     *,MATG
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
C
C
C***********************************************************************
C     FORMIRANJE INT.FUNKCIJA I NJIHOVIH IZVODA I JAKOBIJANA U GAUSOVOJ
C      TA$KI
C***********************************************************************
C
      DIMENSION XE(40,3),H(NTE,*),CORD(NP,*)
      DIMENSION NOP(NE,*)
C
C  ODREDJIVANJE POLO@AJA CVPROVA ZA RACUNANJE JAKOBIJANA
C      T.L. POCETNE KOORDINATE CORD(I,J)
C      U.L. KORIGOVANE KOORDINATE X,Y,Z
      DO 20 NC=1,NTE
      IF(NOP(NLM,NC).EQ.0) GO TO 10
      JJ=NOP(NLM,NC)
C     NCV=NCV+1
      GO TO (30,31),INEL
C ZA U.P.
   30 CONTINUE
C  30 XE(NC,1)=X(JJ)
C     XE(NC,2)=Y(JJ)
C     XE(NC,3)=Z(JJ)
C     GO TO 20
C ZA T.L.
   31 XE(NC,1)=CORD(JJ,1)
      XE(NC,2)=CORD(JJ,2)
      XE(NC,3)=CORD(JJ,3)
   20 CONTINUE
C
      GO TO(1,2,3,4,1,1),NCV
    1 CONTINUE
C***
C     ELEMENT SA DVA $VORA
C***
    2 H(1,1)=(1.-R)/2.
      H(1,2)= -0.5
      H(2,1)=(1.+R)/2.
      H(2,2)=0.5
      GO TO 10
C***
C     ELEMENT SA TRI $VORA
C***
    3 H(1,1)=R*(R-1.)/2.
      H(1,2)=R-0.5
      H(2,1)=1.-R*R
      H(2,2)=-2.*R
      H(3,1)=R*(1.+R)/2.
      H(3,2)=R+0.5
      GO TO 10
C**
C     ELEMENT SA 4 $VORA
C***
    4 DVS=9./16.
      R2=R*R
      H(1,1)=(3.*R+1.)*(3.*R-1.)*(1.-R)/16.
      H(1,2)=(-27.*R2+18.*R+1)/16.
      H(2,1)=DVS*(1.-R2)*(1.-3.*R)
      H(2,2)=DVS*(9.*R2-2.*R-3.)
      H(3,1)=DVS*(1.-R2 )*(1.+3.*R)
      H(3,2)=DVS*(-9.*R2-2.*R+3.)
      H(4,1)=(1.+R)*(1.+3.*R)*(3.*R-1.) /16.
      H(4,2)=(27.*R2+18.*R-1.)/16.
      GO TO 10
C***
C     FORMIRANJE JAKOBIJANA
C***
   10 CONTINUE
      DO 11 I=1,3
      DO 11 J=1,3
      XJ(I,J)=0.
      GO TO (12,13,14),I
   12 DO 15 K=1,NCV
   15 XJ(1,J)=XJ(1,J)+H(K,2)*(XE(K ,J)+(S*BK(K)+DLS(K))*VS(K,J)+
     1(T*AK(K)+DLT(K))*VT(K,J))
      GO TO 11
   13 DO 16 K=1,NCV
   16 XJ(2,J)=XJ(2,J)+H(K,1)*BK(K)*VS(K,J)
      GO TO 11
   14 DO 17 K=1,NCV
   17 XJ(3,J)=XJ(3,J)+H(K,1)*AK(K)*VT(K,J)
   11 CONTINUE
C***
C     DETERMINANTA JAKOBIJANA
C***
C
      DET=XJ(1,1)*(XJ(2,2)*XJ(3,3)-XJ(3,2)*XJ(2,3))+XJ(1,2)*(XJ(3,1)*XJ(
     12,3)-XJ(2,1)*XJ(3,3))+XJ(1,3)*(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))
      RETURN
      END
