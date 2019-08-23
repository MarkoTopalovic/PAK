C====================================================================
C
C  UPRAVLJACKI PROGRAMI ZA POZIVANJE PROGRAMA ZA RACUNANJE
C  MATRICA MASA ELEMENATA
C
C  SUBROUTUNE IGMAS
C             MASIG
C             INTMG
C             READMG
C             MATMAS
C             KONSTM
C=====================================================================
      SUBROUTINE IGMAS
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
C
      LAG = LMAX
C
      CALL READMG(A(LAG))
C
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX = LMAX + 1
      LAE=LMAX
C
C PROSTOR U VEKTORU AE(RADNI VEKTOR ELEMNTA)
C
      MAXAE=(20*NTE+NT6*(NT6+1)/2)*IDVA
      LMAX = LAE + MAXAE
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX = LMAX + 1
      IF(LMAX.LT.MTOT) GO TO 70
      WRITE(IZLAZ,2009) LMAX,MTOT
 2009 FORMAT(' ',' NEDOVOLJNA DIMENZIJA U OSNOVNOM RADNOM VEKTORU A ',
     *'ZA UCITAVANJE ULAZNIH PODATAKA'/' ','POTREBNA DIMENZIJA ',
     *'LMAX = ',I10/' ','RASPOLOZIVA DIMENZIJA, MTOT=',I10)
      STOP
  70  CALL MASIG(A(LAG),A(LAE))
      RETURN
      END
C=====================================================================
      SUBROUTINE MASIG(AG,AE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C POTPROGRAM ZA POZIVANJE POTPROGRAMA ZA INTEGRALJENJE MATRICE
C KRUTOSTI ELEMENATA
C
      include 'paka.inc'
      
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ITERBR/ ITER
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /IGRGLA/ NGRU,NGSU,NGTU,NGS,NTE,NUD,IANAL,JTAU,NWE,NT6,INEL
     *,MATG
      COMMON /IGREP/  LIGNOP,LIGNDS,LIGDST,LIGCTR,MXAUG,LMXAUG,LAG,MAXAE
     *               ,LICTRP,LIGX  ,LIGY  ,LIGZ,LAE,LLM,LIGMAT
C
      DIMENSION AG(*),AE(*)
      REAL AG ,AE
C
      LHE = 1
      LH  = LHE + 3*NT6*IDVA
      LSKE= LH  + NTE*2*IDVA
      LMXAE=LSKE + NT6*(NT6+1)/2*IDVA
C
      CALL INTMG(A(LCORD),AG(LIGNOP),AG(LIGNDS),AG(LIGMAT),
     *           AG(LIGDST),AG(LIGCTR),A(LMAXA),A(LSK),
     *           AG(LLM),AE(LHE),AE(LH),AE(LSKE),A(LGUSM),LH)
C
      RETURN
      END
      SUBROUTINE INTMG(CORD,NOP,NDST,MATIG,DST,CTR,MAXA,SK,
     1                 LM,HE,H,SKE,GUSM,LH)
      USE DRAKCE8
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'paka.inc'
      
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
C
      DIMENSION NOP(NE,*),MAXA(*),NDST
     1(NE,*),DST(NUD,*),CTR(NE,6,*)
      DIMENSION XG(55),WGT(55),
     1LM(NT6,*),NREF(11)
      DIMENSION CORD(NP,*),MATIG(*),GUSM(99,*)
      DIMENSION H(NTE,*),SKE(*),SK(*)
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
      LMHT1=LMHT-1
C
C  IZJEDNACAVANJE SA NULOM POCETNIH NAPONA
      LH1=LH-1
C***
C     PETLJA PO ELEMENTIMA
C***
      DO 10 NLM=1,NE
C-----
      NWE=NT6*(NT6+1)/2
C-----
C
C UNULJAVANJE MATRICE MASA ELEMEMENTA
C
      CALL CLEAR(SKE,NWE)
C
C PODACI O MATERIJALU
      MATG=MATIG(NLM)
      GUST = GUSM(NMODM,MATG)
      JNGS=0
      CALL KONSTM(NOP,CORD,NDST,CTR,DST)
C***
C     PETLJE PO GAUSOVIM TACKAMA
C***
      ND =6*NCV
      DO 20 NGR=1,NGRU
      JGR=NREF(NGRU)+NGR
      R=XG(JGR)
      WR=WGT(JGR)
      DO 20 NGSS=1,NGSU
      JGR=NREF(NGSU)+NGSS
      S=XG(JGR)
      WS=WGT(JGR)
      DO 20 NGT=1,NGTU
      JGR=NREF(NGTU)+NGT
      T=XG(JGR)
      WT=WGT(JGR)
C***
C   ODREDJIVANJE MATRICE JAKOBIJANA
C***
      CALL JACKG(CORD,NOP,H)
C
      IF (DET.GT.1.D-10) GO TO 7
      WRITE(3,1333)NLM,NGR,NGSS,NGT,DET
 1333 FORMAT(' ','S T O P   DET.LE.1.D-10'/' ','NLM =',I2/' ','NGR =',I2
     1/' ','NGSS=',I2/' ','NGT =',I2/' DET =',D12.5)
      STOP
    7 CONTINUE
C***
C     MATRICE ELEMENTA U GAUS.TA$KI
C***
      WTU=WR*WS*WT
      JNGS=JNGS+1
C
      CALL MATMAS(HE,H,SKE,LM,GUST)
C
   20 CONTINUE
      IF (TIPTACKANJA.EQ.1) THEN
      CALL SPAKUJ(SK,MAXA,SKE,LM(1,NLM),6*NTE)
      ELSE
      CALL SPAKUJMT(SK,MAXA,SKE,LM(1,NLM),6*NTE)
      ENDIF
   10 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE READMG(AG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  GLAVNI UPRAVLJACKI PROGRAM ZA UCITAVANJE ULAZNIH PODATAKA
C  ZA RACUNANJE MATRICA MASA
C
      include 'paka.inc'
      
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ITERBR/ ITER
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /IGRGLA/ NGRU,NGSU,NGTU,NGS,NTE,NUD,IANAL,JTAU,NWE,NT6,INEL
     *,MATG
      COMMON /IGREP/  LIGNOP,LIGNDS,LIGDST,LIGCTR,MXAUG,LMXAUG,LAG,MAXAE
     *               ,LICTRP,LIGX  ,LIGY  ,LIGZ,LAE,LLM,LIGMAT
      COMMON /ZAPISI/ LSTAZA(5)
C
      DIMENSION AG(*)
      REAL AG
C
C POZIVANJE PROGRAMA ZA ULAZNE PODATKE
C
      LSTAZA(1)=LMAX8
      READ(IELEM,REC=LMAX8)
     *NGRU,NGSU,NGTU,NTE,NUD,
     *MXAUG,LIGNOP,LIGNDS,LIGMAT,LIGDST,LIGCTR,LICTRP,LLM
C
      LSTAZA(2)=LMAX8+1
C
      CALL READDD(AG(LIGCTR),MXAUG/IDVA,IELEM,LMAX8,LDUZI)
      LMAX=LAG+MXAUG
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      RETURN
      END
C=========================================================================
      SUBROUTINE MATMAS(HE,H,SKE,LM,GUST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C==========================================================
C
C     FORMIRANJE MATRICA MASA ELEMENATA U GAUS.TA$KI
C==========================================================
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ITERBR/ ITER
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /IGRGLA/ NGRU,NGSU,NGTU,NGS,NTE,NUD,IANAL,JTAU,NWE,NT6,INEL
     *,MATG
      COMMON /IGREP/  LIGNOP,LIGNDS,LIGDST,LIGCTR,MXAUG,LMXAUG,LAG,MAXAE
     *               ,LICTRP,LIGX  ,LIGY  ,LIGZ,LAE,LLM,LIGMAT
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /DUPLAP/ IDVA
      COMMON/ELKG/DET,WTU,R,S,T,NLM,JNGS,NCV,ND,PI2,SS,GR(3),GS(3)
     1,XJ(3,3),AK(40),BK(40),DLS(40),DLT(40),VS(40,3),VT(40,3),V1(40,3)
C
      DIMENSION HE(3,*),HP(3,6),H(NTE,*),SKE(*),LM(NT6,*)
C
C  H - INTERPOLACIONE FUNKCIJE ( ODREDJUJU SE U JACKG)
C  HP- POMOCMA INTERPOLACIONA MATRICA (ODNOSI SE NA CVOR)
C  HE- INTERPOLACIONA MATRICA ELEMENTA (ZBIR MATRICA HP PO CVOROVIMA);
C      KORISTI SE ZA MATRICU MASA ELEMNTA
C
C      M=RO*INTEGRAL(HET * HE)DV
C
C
      ND=NCV*6
      ND3=NCV*3
C***
C     OSNOVNA PETLJA PO BROJU $VOROVA
C
      DO 10 NC=1,NCV
      JB=(NC-1)*6
      DO 5 I=1,3
      DO 5 J=1,6
    5 HP(I,J)=0.
C***
      SP=S
      TP=T
      SP=S*BK(NC)
      TP=T*AK(NC)
      DO 65 I=1,3
   65 IF(DABS(V1(NC,I)).LT.1.D-08) V1(NC,I)=0.
C
C  FORMIRANJE INTEPOLACIONIH FUNKCIJA CVOROVA
C
C      DLS = DST(NDST(NLM,NC),1)
C      DLT = DST(NDST(NLM,NC),2)
      HP(1,1)= H(NC,1)
      HP(1,5)= H(NC,1)*((TP+DLT(NC))*VT(NC,3)+(SP+DLS(NC))*VS(NC,3))
      HP(1,6)=-H(NC,1)*((TP+DLT(NC))*VT(NC,2)+(SP+DLS(NC))*VS(NC,2))
      HP(2,2)= H(NC,1)
      HP(2,4)=-H(NC,1)*((TP+DLT(NC))*VT(NC,3)+(SP+DLS(NC))*VS(NC,3))
      HP(2,6)= H(NC,1)*((TP+DLT(NC))*VT(NC,1)+(SP+DLS(NC))*VS(NC,1))
      HP(3,3)= H(NC,1)
      HP(3,4)= H(NC,1)*((TP+DLT(NC))*VT(NC,2)+(SP+DLS(NC))*VS(NC,2))
      HP(3,5)=-H(NC,1)*((TP+DLT(NC))*VT(NC,1)+(SP+DLS(NC))*VS(NC,1))
C
C FORMIRANJE INTERPOLACIONE MATRICE ELEMENTA
C
      DO 11 I=1,3
      DO 11 K=1,6
      JBK=JB + K
      HE(I,JBK) = HP(I,K)
   11 CONTINUE
C
   10 CONTINUE
C
C  MATRICA MASA ELEMENTA
C
      KK=0
      DO 20 I=1,ND
      DO 20 J=I,ND
      KK=KK+1
      IF(LM(I,NLM).EQ.0.OR.LM(J,NLM).EQ.0) GO TO 20
      XX=0.
      DO 21 K=1,3
   21 XX = XX + HE(K,I)*HE(K,J)
      SKE(KK) = SKE(KK) + XX*DET*WTU*GUST
   20 CONTINUE
C
      RETURN
      END
      SUBROUTINE KONSTM(NOP,CORD,NDST,CTR,DST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C*************************************************
C     FORMIRANJE KONSTANTNIH VEKTORA ELEMENATA   *
C     ZA MATRICU MASA                            *
C*************************************************
C
      COMMON/ELKG/DET,WTU,R,S,T,NLM,JNGS,NCV,ND,PI2,SS,GR(3),GS(3)
     1,XJ(3,3),AK(40),BK(40),DLS(40),DLT(40),VS(40,3),VT(40,3),V1(40,3)
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /IGRGLA/ NGRU,NGSU,NGTU,NGS,NTE,NUD,IANAL,JTAU,NWE,NT6,INEL
     *               ,MATG
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
C
      DIMENSION NOP(NE,*),NDST(NE,*),
     1CTR(NE,6,*),DST(NUD,*),XE(40,3)
      DIMENSION XS(3),XT(3)
      DIMENSION CORD(NP,*)
C***
C     OSNOVNA PETLJA PO CVOROVIMA ELEMENTA
C***
      NCV=0
      DO  10 NC=1,NTE
      IF(NOP(NLM,NC).EQ.0) GO TO 10
      JJ=NOP(NLM,NC)
      NCV=NCV+1
      XE(NC,1)=CORD(JJ,1)
      XE(NC,2)=CORD(JJ,2)
      XE(NC,3)=CORD(JJ,3)
      KK=0
      DO 2 I=1,3
      XS(I)=CTR(NLM,I,NC)
      XT(I)=CTR(NLM,I+3,NC)
      IF(DABS(XS(I)).LT.1.D-20.AND.DABS(XT(I)).LT.1.D-20) GO TO 2
      KK=KK+1
    2 CONTINUE
      IF(KK.NE.0) GO TO 1
      DO 4 I=1,3
      VS(NC,I)=VS(1,I)
      VT(NC,I)=VT(1,I)
      AK(NC)=AK(1)
    4 BK(NC)=BK(1)
      GO TO 3
    1 XX=XS(1)-XE(NC,1)
      YY=XS(2)-XE(NC,2)
      ZZ=XS(3)-XE(NC,3)
      BK(NC)=DSQRT(XX*XX+YY*YY+ZZ*ZZ)
      VS(NC,1)=XX/BK(NC)
      VS(NC,2)=YY/BK(NC)
      VS(NC,3)=ZZ/BK(NC)
      XX=XT(1)-XE(NC,1)
      YY=XT(2)-XE(NC,2)
      ZZ=XT(3)-XE(NC,3)
      AK(NC)=DSQRT(XX*XX+YY*YY+ZZ*ZZ)
      VT(NC,1)=XX/AK(NC)
      VT(NC,2)=YY/AK(NC)
      VT(NC,3)=ZZ/AK(NC)
    3 KK=NDST(NLM,NC)
      DLS(NC)=DST(KK,1)
      DLT(NC)=DST(KK,2)
      V1(NC,1)=VS(NC,2)*VT(NC,3)-VS(NC,3)*VT(NC,2)
      V1(NC,2)=VS(NC,3)*VT(NC,1)-VS(NC,1)*VT(NC,3)
      V1(NC,3)=VS(NC,1)*VT(NC,2)-VS(NC,2)*VT(NC,1)
   10 CONTINUE
      RETURN
      END
