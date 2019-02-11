C$DEBUG      
C=========================================================================
      SUBROUTINE RACUN(TT1,SILE,NEL,ID,NZAD,ZADVRE,NGPSIL,
     &MAXA,CORD,SKEF,SKEFN,KONT,FZAPR,VREME,TABF,TT10,UBRZ,UBRZ0,
     &AK,VECTJ,IVECT,POVSIL,GRADJN,ITFMAX,AKONST,NASLOV,
     &ICUR,VG,GG,KOJK)
      USE STIFFNESS
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /TRENUT/ TT21,H,HP,ZVHX,ZVHY,HV2,HV3,ZVXT,ZVYT,DETJ,
     1DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,B,BT,GRADH,GRADHT,TT210,
     2DEBLJ,ZHP,H2N,H2NP,GH2NT,H2NN
      COMMON /TREN1/ FPS1,FPS2,FP1,FP2,VNN,PERM,CORDY
      include 'paka.inc'
      INCLUDE 'mpif.h'
!      COMMON A(17000)
!      REAL A
      COMMON /KRITER/ IDOKL,MAXCUR
      COMMON /TACNOS/ EPSTR,MAXIT,NJRAP
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /NELDIM/ NDIMEL
      COMMON /ZADPO/ LNZADJ
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /NDESUK/ NDES,IDPRIT,IFORM
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT,GAMA
      COMMON /SCRATC/ ISCRC
      COMMON /VREPER/ NPER,NTABFT
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /PROMEN/ NJUTN,INDOPT,INDPRO
      COMMON /DODAT/ NDIMM
      COMMON /LSK1/ LSK
      COMMON /JEDANP/ INDJED,NBRF,NGL,INDTLO
      COMMON /FLUX/ KISA
      COMMON /DYNAM/ NDIN
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /STAC/ NSTAC,NSLP,IBROJ
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /VDP/ DT,NKORP,NN,NZAV
      COMMON /SILEZAP/ FZAP(2,9)
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL
      COMMON /SOLVER/ ISOLVER

      CHARACTER*80 NASLOV
      DIMENSION TT1(*),SILE(*),ZADVRE(*),CORD(3,*)
      DIMENSION TT10(*),UBRZ(*),UBRZ0(*)
      DIMENSION NEL(NDIMM,*),ID(1,*),NZAD(3,*),NGPSIL(4,*),MAXA(*)
      DIMENSION SKEF(NDES,*),AKONST(3,5,*),SKEFN(*)
      DIMENSION H2N(2,18),H2NP(2,4),ZHP(4,2),GH2NT(18,2)
      DIMENSION H2NN(2,18)
      DIMENSION POVSIL(2,*),ITFMAX(*),PR1(4)
      DIMENSION VG(3,NET,*),GG(3,NET,*)
C
      DIMENSION KONT(5,MAXLIN,*),KOJK(*)

      DIMENSION X(9),Y(9),TT21(44),TT210(44),CORDY(4)
      DIMENSION R(3,3),S(3,3),W(3,3)
      DIMENSION RK(3,3),SK(3,3),WK(3,3)
      DIMENSION H(9),ZVHX(9),ZVHY(9),HP(4)
      DIMENSION B(4,18),BT(18,4),GRADHT(4,2),GRADH(2,4)
      DIMENSION FZAPR(2,*),VECTJ(2,*),GRADJN(2,*)
      DIMENSION IVECT(*),LM2(44),ICUR(*)
      DIMENSION PERM(2,2)
      DIMENSION VREME(NPER,950)
      DIMENSION TABF(2,NTABFT,*)
      DIMENSION AM1(9,9),AK(NDES,*)
      DIMENSION AJX(8),AJY(8)
      DIMENSION AJX1(8),AJY1(8),NRED(4),NREDK(4)
      DIMENSION RS2(9)
      DIMENSION RS3(9)
      DIMENSION F36(44)
      DIMENSION AFIFI1(9,9),AFIFI(9,9)
      DIMENSION PUS(9,9)
      DATA NRED/4,2,1,3/ 
      DATA NREDK/9,3,1,7/ 
      integer ierr, myid
      integer*8 LM2

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      IF (myid.ne.0) goto 1234
C

C INDIKATOR ZA OSU TEZINE
      IF (IOSA.EQ.0) IOSA=2

C INDIKATOR ZA TRAZENJE IZOHIPSI PO FAJLOVIMA:
       INDJED=0

C ZBOG PAK-a NEQ-BROJ JEDNACINA .EQ. JEDN
       NEQ=JEDN
       ICONT=1

C IDOKL JE BROJ CVOROVA KOJI SE UZIMAJU KAO AKTIVNI NA LINIJI CURENJA 
C ICUR JE NIZ U KOME SE NALAZE CVOROVI NA LINIJI CURENJA
C
C  INICIJALIZACIJA KONVERGENCIJE ZA LINIJU CURENJA
       KONV1=1  

       MAXCUR=0
       IDOKL=0
       CALL INIDOK(NZAD,ICUR,CORD)

C UTVRDJIVANJE KOD KOLIKO ELEMENATA SE POJAVLJUJE JEDAN CVOR
      CALL MVECT(IVECT,NEL)

      DO I=1,JEDN
       TT1(I)=0.D0
       UBRZ(I)=0.D0
      ENDDO     
       
C DEFINISANJE CVOROVA ZA PROTOKE
      DO KK=1,NKONT
       DO II=1,LIN(KK)
        KOJK(KONT(2,II,KK))=KOJK(KONT(2,II,KK))+1
        KOJK(KONT(3,II,KK))=KOJK(KONT(3,II,KK))+1
       ENDDO
      ENDDO	  
C
      R(3,1)=-0.7745966692415
      R(3,2)=0.0
      R(3,3)=0.77459666924148
      S(3,1)=-0.7745966692415
      S(3,2)=0.0
      S(3,3)=0.77459666924148
      W(3,1)=0.55555555555556
      W(3,2)=0.88888888888889
      W(3,3)=0.55555555555556
C
      RK(3,1)=-1.0
      RK(3,2)=0.0
      RK(3,3)=1.0
      SK(3,1)=-1.0
      SK(3,2)=0.0
      SK(3,3)=1.0
      WK(3,1)=0.33333333333333
      WK(3,2)=1.33333333333333
      WK(3,3)=0.33333333333333
C
      R(2,1)=-0.5773502691896
      R(2,2)=0.57735026918963
      S(2,1)=-0.5773502691896
      S(2,2)=0.57735026918963
      W(2,1)=1.000000000000000
      W(2,2)=1.000000000000000

      R(1,1)=0.0
      S(1,1)=0.0
      W(1,1)=2.0
      
      PI=4.D0*DATAN(1.D0)
            
C===========================================================================
      IAXIS=3
C AKO JE PROBLEM OSNOSIMETRICAN
      IF (INDAX.EQ.1) IAXIS=4              
C
      INULTO=0
       KKORAK=1
      VVREME=0.D0
C=========================================================================
C
      DO 600 NNPER=1,NPER
      KORAK=0
   35 KORAK=KORAK+1
C
C=========================================================================   
C  RAD SA STACIONARNIM STANJEM U PRVOM KORAKU I NESTACIONARNIM U OSTALIM
C  ( GEOLOSKI PROFILI ZA CERNE)
C
      IF (KKORAK.EQ.1) NNSTAC=NSTAC
      IF (KKORAK.EQ.1.AND.NNSTAC.EQ.2) NSTAC=1
      IF (KKORAK.GT.1.AND.NNSTAC.EQ.2) NSTAC=0
C
      IF (KKORAK.GT.0) THEN
      TIME=VREME(NNPER,KORAK)
      VVREME=VVREME+TIME
      ELSE
       TIME=1.D7
       KORAK=KORAK-1
      ENDIF

      DO I=1,NPT 
       DO J=1,2
        VECTJ(J,I)=0.D0
        GRADJN(J,I)=0.D0
       ENDDO
      ENDDO

      DO I=1,JEDN
        TT10(I)=TT1(I)
        UBRZ0(I)=UBRZ(I)
      ENDDO

      ITER=0
      ITER1=0
      INDSK=1
      INDPR=0

C==========================================================================
C LABELA 100:
C==========================================================================
 100  IF (ITER.GT.MAXIT) THEN
        WRITE(IIZLAZ,3001) MAXIT
 3001 FORMAT(//
     1,'DOSTIGNUT JE MAKSIMALAN BROJ OD ',I5,' ITERACIJA BEZ POSTIZANJA 
     1KONVERGENCIJE!!'//)
       STOP
       ENDIF

	IF (NBLOCK.GT.1) THEN
       OPEN (ISCRC,FILE='ZSKLIN',FORM='UNFORMATTED',STATUS='UNKNOWN')
       REWIND ISCRC
       ELSE
       KC=NWK
      ENDIF

C===========================================================================
C ZA CVOROVE NA LINIJI PROCURIVANJA INICIJALIZACIJA ZA NESTACIONARNO
C STRUJANJE - START JE UVEK OD PRVOG NAJNIZEG CVORA PA NA VISE
C      IF (ITER.EQ.0.AND.MAXCUR.NE.0) IDOKL=1
      IF (ITER.EQ.0.AND.MAXCUR.NE.0) IDOKL=MAXCUR
C===========================================================================
      IF (MAXCUR.GT.0) THEN
        if(ISOLVER.NE.-11) CALL CLEAR(A(LSK),KC)
      ELSE
        IF (ITER.EQ.0.and.ISOLVER.NE.-11) CALL CLEAR(A(LSK),KC)
      ENDIF

      DO I=1,JEDN
       SILE(I)=0.D0
       UBRZ(I)=0.D0
      ENDDO 


C ODEDIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
      IBRGT=3
C      IF (NDIM.EQ.4) IBRGT=2
      IBRGT=2


C GLAVNA PETLJA PO ELEMENTIMA
      DO 400 NBREL=1,NET

      NMAT=NEL(NDIM+1,NBREL)
      NEMA=NEL(NDIM+2,NBREL)

      IF (NMAT.LE.0.OR.NEMA.GT.0) THEN
         WRITE(3,*) 'NBREL,NMAT,NEMA',NBREL,NMAT,NEMA
         WRITE(*,*) 'NBREL,NMAT,NEMA',NBREL,NMAT,NEMA
C         STOP
      ENDIF
      IF (NMAT.LE.0) GOTO 400
      IF (NEMA.GT.0) GOTO 400
      
       PERM(1,1)=AKONST(3,1,NMAT)
       PERM(1,2)=0.
       PERM(2,2)=AKONST(3,2,NMAT)
       PERM(2,1)=0.
      
C=========================================================================
      DO 130 KLM=1,NDIM
      X(KLM)=CORD(1,NEL(KLM,NBREL))
      Y(KLM)=CORD(2,NEL(KLM,NBREL))
C
      LM2(KLM)=ID(1,NEL(KLM,NBREL))
 130  CONTINUE
C=======================================================================
      NDIMEL=NDIM
      CALL CLEARD(TT21,NDIM)
      CALL PREB(TT21,TT1,ID,NEL,NBREL)
      CALL CLEARD(TT210,NDIM)
      CALL PREB(TT210,TT10,ID,NEL,NBREL)
 
C=======================================================================
      DO 163 K=1,NDIM
      DO 162 N=1,NDIM
        AM1(K,N)=0.D0    
        AFIFI(K,N)=0.D0    
        AFIFI1(K,N)=0.D0    
        PUS(K,N)=0.D0    
  162 CONTINUE
  163 CONTINUE

C POVRSINSKE SILE I ZAPREMINSKE SILE:
      DO 199 I=1,NDIM
      RS2(I)=0.D0
      RS3(I)=0.D0
 199  CONTINUE
C
C INTEGRACIJA U GAUSOVIM TACKAMA

       DO II=1,NDIM
         PR1(II)=TT21(II)-CORD(IOSA,NEL(II,NBREL))
        ENDDO

      NGAUS=0
      DO 180 I=1,IBRGT
      DO 170 J=1,IBRGT
      NGAUS=NGAUS+1
      CALL INTERP(R(IBRGT,I),S(IBRGT,J),0,0)

      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ*DEBLJ
      
      IF ( (((ITER.GT.0).AND.(NSTAC.EQ.1))
     1.OR.((KKORAK.GT.0).AND.(NSTAC.EQ.0))) .AND.(INDFS.EQ.1)) THEN
      
       PR=0.D0
        DO II=1,NDIM
         PR=PR+H(II)*(TT21(II)-CORD(IOSA,NEL(II,NBREL)))
        ENDDO

       IF( PR.LT.0.D0 .AND. INDFS.EQ.1 )THEN
        ANPUS=AKONST(3,3,NMAT)*ZASIC
        PERM(1,1)=ZASIC*AKONST(3,1,NMAT)
        PERM(2,2)=ZASIC*AKONST(3,2,NMAT)
       ELSE
        ANPUS=0.D0
        PERM(1,1)=0.0D0
        PERM(2,2)=0.0D0
       ENDIF
      ELSE
        PERM(1,1)=0.0D0
        PERM(2,2)=0.0D0
        ANPUS=0.D0
      ENDIF
 
      DO 165 K=1,NDIM
      DO 164 N=1,NDIM

         AM1(K,N)=AM1(K,N)+AKONST(3,3,NMAT)*H(K)*H(N)*WDT
         AFIFI(K,N)=AFIFI(K,N)+(AKONST(3,1,NMAT)*ZVHX(K)*ZVHX(N)+
     1AKONST(3,2,NMAT)*ZVHY(K)*ZVHY(N))*WDT
C
        IF (INDFS.EQ.1) THEN
         PUS(K,N)=PUS(K,N)+ANPUS*H(K)*H(N)*WDT
         AFIFI1(K,N)=AFIFI1(K,N)+(PERM(1,1)*ZVHX(K)*ZVHX(N)+
     1PERM(2,2)*ZVHY(K)*ZVHY(N))*WDT
        ENDIF
C
  164 CONTINUE
  165 CONTINUE
C      if(nbrel.eq.1) write(3,*) i,j,R(IBRGT,I),S(IBRGT,J)
C      if(nbrel.eq.1) write(3,*) perm(1,1),perm(2,2),wdt
C      if(nbrel.eq.1) call wrr(h,9,'h   ')
C      if(nbrel.eq.1) call wrr(zvhx,9,'zvhx')
C      if(nbrel.eq.1) call wrr(zvhy,9,'zvhy')
  170 CONTINUE
  180 CONTINUE

C INCIJALIZACIJA MATRICE SKEF I F36
      DO 260 I=1,NDES
      DO 258 J=1,NDES
      SKEF(I,J)=0.D0
 258  CONTINUE
      F36(I)=0.D0
 260  CONTINUE
C=========================================================================
C UNULJAVANJE VEKTORA RS2
	INDPUN=0
       DO I=1,NDIM
         RS2(I)=0.D0
       ENDDO
C KADA KISA PADA BAS NA LINIJU SLOBODNE POVRSINE:
CZ      IF ((DABS(ZASIC).GT.1.D-10).AND.(INDFS.EQ.1).AND.(MAXSIL.GT.0)) 
CZ     &THEN
CZ        IVRF=NGPSIL(4,1)
CZ        CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(IVRF),IVRF)
CZ        CALL PROCEN(RS2,PR1,FK1,NEL,NBREL,CORD) 
CZ        INDPUN=1
CZ       ENDIF
C============================================================================
C POVRSINSKE SILE:
C==============================================================================
      IF (INDPUN.EQ.1) GOTO 265

      DO 250 JBRPS=1,MAXSIL
      IF (NBREL.EQ.NGPSIL(1,JBRPS)) THEN
      NODE1=NGPSIL(2,JBRPS)
      NODE2=NGPSIL(3,JBRPS)
      N1=NEL(1,NBREL)
      N2=NEL(2,NBREL)
      N3=NEL(3,NBREL)
      N4=NEL(4,NBREL)
      IVRF=NGPSIL(4,JBRPS)
      CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(IVRF),IVRF)
      FPS1=POVSIL(1,JBRPS)*FK1
      FPS2=POVSIL(2,JBRPS)*FK1
C      
      NUMGAU=IBRGT
      DO 225 I=1,NUMGAU
C 
C KONSTANTNO R (KSI)
C
      IF ((NODE1.EQ.N1 .AND. NODE2.EQ.N4).OR.
     1(NODE1.EQ.N4 .AND. NODE2.EQ.N1)) THEN
      CALL INTERP(1.D0,S(NUMGAU,I),1,1)
       GOTO 220
      ENDIF 
 
      IF ((NODE1.EQ.N2 .AND. NODE2.EQ.N3).OR.
     1(NODE1.EQ.N3 .AND. NODE2.EQ.N2)) THEN
      CALL INTERP(-1.D0,S(NUMGAU,I),1,1)
       GOTO 220
      ENDIF 
C 
C KONSTANTNO S (ETA)
C
      IF ((NODE1.EQ.N1 .AND. NODE2.EQ.N2).OR.
     1(NODE1.EQ.N2 .AND. NODE2.EQ.N1)) THEN
      CALL INTERP(R(NUMGAU,I),1.D0,2,1)
       GOTO 220
      ENDIF 

      IF ((NODE1.EQ.N3 .AND. NODE2.EQ.N4).OR.
     1(NODE1.EQ.N4 .AND. NODE2.EQ.N3)) THEN
      CALL INTERP(R(NUMGAU,I),-1.D0,2,1)
       GOTO 220
      ENDIF 

 220    IF (KISA.EQ.1) THEN
C KADA PADA KISA ONDA FLUX DELUJE U PRAVCU -Y OSE 
C TAKO DA MORA DA SE PROJEKTUJE NA NORMALU LINIJE PO KOJOJ DELUJE
          FAK=DABS(CORD(1,NODE1)-CORD(1,NODE2))
          WDTS=FAK
        ELSE
          WDTS=W(NUMGAU,I)*DETJS*DEBLJ
        ENDIF

       DO K=1,NDIM
        DUM=TT21(K)-CORD(IOSA,NEL(K,NBREL))
      IF ( (((ITER.GT.0).AND.(NSTAC.EQ.1))
     1.OR.((KKORAK.GT.0).AND.(NSTAC.EQ.0)))
     1.AND.INDFS.EQ.1.AND.DUM.LT.0.) THEN
      ELSE
        RS2(K)=RS2(K)+H(K)*FS2*WDTS
      ENDIF
       ENDDO
 
 225  CONTINUE
 
      ENDIF
 250  CONTINUE

C      ENDIF
C KRAJ PETLJE AKO JE ITER>0
C=========================================================================

C PAKOVANJE POVRSINSKIH SILA SA DESNE STRANE
 265   DO I=1,NDIM
         F36(I)=F36(I)+RS2(I)
       ENDDO

C=========================================================================
      DO I=1,NDES
      DO J=1,NDES
       AK(I,J)=0.D0
      ENDDO
      ENDDO
C PAKOVANJE MATRICA AMUU,AMUW,AMWW,AKUU,AKUW,AKWW U MATRICE AMM,AC,AK
      DO I=1,NDIM
      DO J=1,NDIM
C sneza 10.02.2010.
c         AK(I,J)=AFIFI(I,J)-AFIFI1(I,J)
         AK(I,J)=AFIFI(I,J)
      IF ((NSTAC.EQ.0).AND.(KKORAK.GT.0)) THEN
c sneza 10.02.2010.
c         AK(I,J)=AK(I,J)+(AM1(I,J)-PUS(I,J))/TIME
         AK(I,J)=AK(I,J)+AM1(I,J)/TIME
        ENDIF
      ENDDO
      ENDDO
C============================================================================
C PAKOVANJE MATRICE AK U MATRICU SKEF
      DO 263 I=1,NDES
      DO 262 J=1,NDES
      SKEF(I,J)=AK(I,J)
  262 CONTINUE
  263 CONTINUE
C=========================================================================
C PAKOVANJE Kfifi*fi(i-1) SA DESNE STRANE
      DO I=1,NDIM
       DO J=1,NDIM
        F36(I)=F36(I)-(AFIFI(I,J)-AFIFI1(I,J))*TT21(J)
       IF ( (NSTAC.EQ.0).AND.(KKORAK.GT.0) ) THEN
         F36(I)=F36(I)-((AM1(I,J)-PUS(I,J))*(TT21(J)-TT210(J)))/TIME
       ENDIF
       ENDDO
      ENDDO
C============================================================================
C VEKTOR PROTOKA PO CVOROVIMA
      DO I=1,NDIM
       N=NEL(I,NBREL)
       JJ=ID(1,N)
       IF(JJ.GT.0) THEN

         POT=0.D0
        DO J=1,NDIM

         POT=POT+(AFIFI(I,J)-AFIFI1(I,J))*TT21(J)
C      	IF (NSTAC.NE.1) THEN
C          POT=POT+((AM1(I,J)-PUS(I,J))*(TT21(J)-TT210(J)))/TIME
C        ENDIF

        ENDDO
         UBRZ(JJ)=UBRZ(JJ)+POT
       ENDIF
      ENDDO
C
       CALL PSKEFN(SKEF,SKEFN)
       CALL SPAKDE (SILE,F36,LM2,NDES)

	IF (ISOLVER.EQ.0) THEN
        IF (MAXCUR.GT.0) THEN
          CALL SPAKUJ(A(LSK),MAXA,SKEFN,LM2,NDES)
        ELSE
cz ovaj uslov ne valja jer ne radi levu stranu za full newton      
C        IF (ITER.EQ.0.OR.NJRAP.EQ.2) 
          IF (ITER.EQ.0) 
     +                  CALL SPAKUJ(A(LSK),MAXA,SKEFN,LM2,NDES)
        ENDIF
      ELSEIF ((ISOLVER.EQ.-11).OR.(ISOLVER.EQ.1)) THEN
C      IMSL Sparse OR MUMPS
        IF(ITER.EQ.0) THEN
C         CALL SPARSE_ELEM_LEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,1)
          CALL sparseassembler_addelemmatrix(NDES,LM2,SKEF)
        ENDIF
      ENDIF
c      if(nbrel.eq.1) call wrr(skefn,ndes*(ndes+1)/2,'ske ')
c      if(nbrel.eq.167) call wrr(f36,ndes,'f36 ')
C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE

C=======================================================================
C CVOROVI NA LINIJI PROCURIVANJA       
C=======================================================================
 	 IF (MAXCUR.GT.0) THEN
        CALL CURI1(CORD,A(LSK),SILE,MAXA,TT1,ID,ICUR)
C       ELSE
C        IF(ITER.EQ.0) CALL CURI1(CORD,A(LSK),SILE,MAXA,TT1,ID,ICUR)
       ENDIF
C==========================================================================
C MUMPS alociranje matrice krutosti
1234  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ITER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ISOLVER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NBLOCK,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(MAXCUR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KKORAK,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        IF (ITER.EQ.0) THEN

         IF(ISOLVER.EQ.0) THEN
          if(myid.eq.0) then
          IF (NBLOCK.GT.1) THEN
           LLM =LRAD
           LSKE=LLM+100
           CALL SPAKUA(A(LSK),MAXA,A(LSKE),LM2,NDES,0,
     &                        A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,H,NZAD)
           CLOSE (ISCRC,STATUS='KEEP')
          ENDIF
          endif
         ELSEIF((ISOLVER.EQ.-11).OR.(ISOLVER.EQ.1)) THEN
          
          if(.not.allocated(rows)) then
            if(myid.eq.0) call sparseassembler_getnz(nonzeros)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(nonzeros,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            allocate(rows(nonzeros),STAT=istat)
            if(istat.ne.0) stop 'error allocating rows'
            allocate(columns(nonzeros),STAT=istat)
            if(istat.ne.0) stop 'error allocating columns'
            allocate(stiff(nonzeros),STAT=istat)
            if(istat.ne.0) stop 'error allocating stiff'
          endif
C           IMSL Sparse
c            CALL SPARSE_BUILD(ISPARSE_N,ISPARSE_NZ,
c     1                DSPARSE_A,ISPARSE_ROWS,ISPARSE_COLS)
          if(myid.eq.0) CALL sparseassembler_getsparse(nonzeros,
     1                 rows,columns,stiff)

          if(myid.eq.0) CALL sparseassembler_kill()

          IF(KKORAK.EQ.1) THEN
              stiff_n = JEDN
          ENDIF
         ENDIF
        ENDIF
C==========================================================================
	 LCMPC=1
        LMPC=1

	IF (MAXCUR.GT.0) THEN
         IF(ISOLVER.EQ.0) THEN
          if(myid.eq.0) CALL ZADLEV(A(LSK),MAXA,A(LNZADJ),NUMZAD,NZAD)
         ELSE
          if(myid.eq.0) CALL SPARS_ZADLEV(stiff_n,noneros,
     1             stiff,rows,columns,A(LNZADJ),NUMZAD,NZAD)
         ENDIF
         CALL RESEN(A(LSK),TT10,MAXA,JEDN,1)
        ELSE
cz proveriti i ovde uslov za full newton        
C        IF (ITER.EQ.0.OR.NJRAP.EQ.2) THEN
        IF (ITER.EQ.0) THEN
cz proveriti rad sa blokovima i zadlev
          if(myid.eq.0) then
          IF (NBLOCK.GT.1) THEN
            LLM =LRAD
            LSKE=LLM+100
            CALL SPAKUA(A(LSK),MAXA,A(LSKE),LM2,NDES,0,
     &          A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,H,NZAD)
            CLOSE (ISCRC,STATUS='KEEP')

          ENDIF
          endif
          
          IF(ISOLVER.EQ.0) THEN
            if(myid.eq.0) CALL ZADLEV(A(LSK),MAXA,A(LNZADJ),NUMZAD,NZAD)
          ELSE
            if(myid.eq.0) CALL SPARS_ZADLEV(stiff_n,noneros,
     1                        stiff,rows,columns,
     1                          A(LNZADJ),NUMZAD,NZAD)
          ENDIF
          CALL RESEN(A(LSK),TT10,MAXA,JEDN,1)
        ENDIF
        ENDIF

       NASL=0
       NASD=0
       IGL=0
       IGD=0
       TOL=1.0D-7
C=======================================================================

       DO 410 I=1,NUMZAD

       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,I)),NZAD(3,I))

       IF (NSTAC.EQ.0) THEN
C===========================================================================
C KADA JE STRUJANJE NESTACIONARNO

        IF( (NZAD(2,I).EQ.1.OR.(NZAD(2,I).EQ.3.AND.KKORAK.EQ.1))
     &                  .AND.(ID(1,NZAD(1,I)).GT.0)) THEN
         JJ=ID(1,NZAD(1,I))
          SILE(JJ)=1.D035*((ZADVRE(I)*FK1)-TT1(JJ))
       ENDIF

      IF ((KKORAK.EQ.1).AND.(NZAD(2,I).EQ.-1).AND.(INDFS.EQ.1)) 
     1SILE(ID(1,NZAD(1,I)))=1.0D35*CORD(IOSA,NZAD(1,I))

        IF( (NZAD(2,I).EQ.2).AND.(ID(1,NZAD(1,I)).GT.0).
     1AND.(CORD(IOSA,NZAD(1,I))-FK1).LE.TOL ) THEN
            JJ=ID(1,NZAD(1,I))
         SILE(JJ)=1.D035*(FK1-TT1(JJ))
        ENDIF

C===========================================================================
       ELSE
C===========================================================================
C KADA JE STRUJANJE STACIONARNO
C
C  ZA VDP TREBA OVAKO
C
      IF (NJUTN.EQ.2) THEN
C
         IF ( (ID(1,NZAD(1,I)).GT.0).AND.(ITER.EQ.0)) THEN
           SILE(ID(1,NZAD(1,I)))=1.0D35*ZADVRE(I)*FK1
         ELSE
          SILE(ID(1,NZAD(1,I)))=0.D0
         ENDIF
C
       ELSE
C
        IF (NZAD(2,I).NE.4) THEN
        IF ( (ID(1,NZAD(1,I)).GT.0) .AND. (ITER.EQ.0) ) THEN
          SILE(ID(1,NZAD(1,I)))=1.0D35*ZADVRE(I)*FK1
         ELSE 
          SILE(ID(1,NZAD(1,I)))=0.D0
         ENDIF
        ENDIF
 
        ENDIF
C===========================================================================

       ENDIF
  410  CONTINUE

      CALL RESEN(A(LSK),SILE,MAXA,JEDN,2)

           DO 440 I=1,JEDN
                TT1(I)=TT1(I)+SILE(I)
 440       CONTINUE

      WRITE(*,*)'ITER= ',ITER
      WRITE(*,*)'PERIOD= ',NNPER
      WRITE(*,*)'STEP= ',KKORAK

      CALL KONVTF(TT1,SILE,KONVV2,1,ID,ITER)

        IF (KKORAK.GT.0) THEN
          
C        CALL CKONV1(ICUR,TT1,CORD,ID,KONV1,UBRZ)

          IF (KONVV2*KONV1.EQ.0) THEN
            ITER=ITER+1
            WRITE(IIZLAZ,*)'ITER=',ITER
            GO TO 100
          ENDIF

        ENDIF

       IF ((NJUTN.NE.0).AND.(INDOPT.EQ.0)) RETURN
  
C=======================================================================
C RACUNANJE GRADIJENTA POTENCIJALA, BRZINA I ZAPREMINSKIH SILA
C=======================================================================
 	DO JJ=1,NPT
          FZAPR(1,JJ)=0.D0
          FZAPR(2,JJ)=0.D0
	ENDDO	
C
C RACUNANJE PROTOKA DUZ KONTURE UNUTAR MREZE
C
      IF (NKONT.GT.0) THEN
        CALL RKONTUR(TT1,ID,KONT,NEL,AKONST,CORD)
      ENDIF
C
C  PETLJA PO ELEMENTIMA
C
      DO 500 NBREL=1,NET
      NMAT=NEL(NDIM+1,NBREL)
      NEMA=NEL(NDIM+2,NBREL)

      IF (NMAT.LE.0) GOTO 500
      IF (NEMA.GT.0) GOTO 500
 
      CALL CLEARD(TT21,NDIM)
      CALL PREB(TT21,TT1,ID,NEL,NBREL)
      CALL CLEARD(TT210,NDIM)
      CALL PREB(TT210,TT10,ID,NEL,NBREL)
C
      DO KLM=1,NDIM
       X(KLM)=CORD(1,NEL(KLM,NBREL))
       Y(KLM)=CORD(2,NEL(KLM,NBREL))
      ENDDO

      DO KL=1,NDIM
       FZAP(1,KL)=0.D0
       FZAP(2,KL)=0.D0
      ENDDO
C INDIKATOR NJUTN-KOTESOVE INTEGRACIJE (0-NE,1-DA)
	KOTES=0
      NGAUS=0
      IBRGT=2
	IBRKT=2
	IF(KOTES.EQ.1)IBRKT=3
      DO  330 I=1,IBRKT
      DO  320 J=1,IBRKT
      NGAUS=NGAUS+1
	IF(KOTES.EQ.1) THEN
         CALL INTERP(RK(IBRKT,I),SK(IBRKT,J),0,0)
         WDT=WK(IBRKT,I)*WK(IBRKT,J)*DETJ*DEBLJ
	ELSE
         CALL INTERP(R(IBRKT,I),S(IBRKT,J),0,0)
         WDT=W(IBRKT,I)*W(IBRKT,J)*DETJ*DEBLJ
	ENDIF
C
C OKVASENOST U GAUS TACKI
        PR=0.D0
        DO II=1,NDIM
           PR=PR+H(II)*(TT21(II)-CORD(IOSA,NEL(II,NBREL)))
        ENDDO
        OKVAS=0.
        IF(PR.GT.0.D0) OKVAS=1.
C=======================================================================
C RACUNANJE GRADIJENTA POTENCIJALA:
C
      GRADFX=0.D0
      GRADFY=0.D0
C INDIKATOR ZA ISKLJUCIVANJE RACUNANJA BRZINA IZNAD FS      
      IF(PR.GT.0.D0.OR.ISIL.EQ.0) THEN
         DO II=1,NDIM
C OVDE PROVERITI BEZ CVOROVA IZNAD FS ZA GAUS INTEGRACIJU         
           GRADFX=GRADFX+ZVHX(II)*TT21(II)
           GRADFY=GRADFY+ZVHY(II)*TT21(II)
         ENDDO
C
C  GRADIJENT PRITISKA (FILTRACIONE SILE PO JEDINICI ZAPREMINE U GAUS TACKI)
C
         IF (IFIL.EQ.1)THEN
           GRADPX=-GAMA*GRADFX
           GRADPY=-GAMA*GRADFY
         ELSEIF(IFIL.EQ.2)THEN
           GRADPX=GAMA*DELTA(IOSA,1)*OKVAS
           GRADPY=GAMA*DELTA(IOSA,2)*OKVAS
         ELSEIF(IFIL.EQ.0)THEN
           GRADPX=-GAMA*(GRADFX-DELTA(IOSA,1)*OKVAS)
           GRADPY=-GAMA*(GRADFY-DELTA(IOSA,2)*OKVAS)
         ENDIF
C
         DO KL=1,NDIM
           IF(PR.LT.0.D0.AND.ISIL.EQ.1)THEN
           ELSE
             FZAP(1,KL)=FZAP(1,KL)+H(KL)*GRADPX*WDT
             FZAP(2,KL)=FZAP(2,KL)+H(KL)*GRADPY*WDT
           ENDIF
         ENDDO
C
	ENDIF
C
        IF ( (PR.LT.0.D0) .AND. (INDFS.EQ.1) )THEN
           AKX=(1.D0-ZASIC)*AKONST(3,1,NMAT)
           AKY=(1.D0-ZASIC)*AKONST(3,2,NMAT)
        ELSE
           AKX=AKONST(3,1,NMAT)
           AKY=AKONST(3,2,NMAT)
        ENDIF

C GRADIJENT U GAUS TACKI
       AJX1(NGAUS)=GRADFX
       AJY1(NGAUS)=GRADFY
       GG(1,NBREL,NGAUS)=GRADFX
       GG(2,NBREL,NGAUS)=GRADFY
C BRZINA U GAUS TACKI
       AJX(NGAUS)=-AKX*GRADFX
       AJY(NGAUS)=-AKY*GRADFY
       VG(1,NBREL,NGAUS)=-AKX*GRADFX
       VG(2,NBREL,NGAUS)=-AKY*GRADFY
C
  320  CONTINUE
  330  CONTINUE
C
	DO KK=1,4
	  NCVOR=NEL(KK,NBREL)
	  KK1=NRED(KK)
	  IF(KOTES.EQ.1) KK1=NREDK(KK)
C BRZINE U CVOROVIMA
          VECTJ(1,NCVOR)=VECTJ(1,NCVOR)+AJX(KK1)/IVECT(NCVOR)
          VECTJ(2,NCVOR)=VECTJ(2,NCVOR)+AJY(KK1)/IVECT(NCVOR)
C GRADIJENT POTENCIJALA U CVOROVIMA
	  GRADJN(1,NCVOR)=GRADJN(1,NCVOR)+AJX1(KK1)/IVECT(NCVOR) 
	  GRADJN(2,NCVOR)=GRADJN(2,NCVOR)+AJY1(KK1)/IVECT(NCVOR) 
C FILTRACIONE SILE U CVOROVIMA	
          FZAPR(1,NCVOR)=FZAPR(1,NCVOR)+FZAP(1,KK)
          FZAPR(2,NCVOR)=FZAPR(2,NCVOR)+FZAP(2,KK)
	ENDDO	
C======================================================================= 
 500  CONTINUE
C======================================================================= 
C POTPROGRAM ZA RACUNANJE PROTOKA PO KONTURI:
       CALL PRKONT1(TT1,ID,UBRZ,CORD,NZAD,IOSA,IVECT,ZADVRE,KOJK)
C========================================================================

C RACUNANJE PROTOKA DUZ KONTURE UNUTAR MREZE
C
      IF (NKONT.GT.0) THEN
        CALL RKONTV(VECTJ,KONT,NEL,CORD,VG)
      ENDIF
C
      WRITE(IIZLAZ,*)'PERIOD NUMBER ',NNPER         
      IF (KKORAK.GT.0) THEN
      CALL IZLLST(ID,TT1,KKORAK,VECTJ,VVREME,FZAPR,VG,GG,NEL,
     1            KOTES)
      CALL STAGP1(TT1,ID,NASLOV,VVREME,KKORAK,1,NPT,18,0,NET,NEL,
     1            KORAK,VECTJ,CORD,GRADJN,FZAPR,NZAD,UBRZ)
C
      CALL STAS09(TT1,CORD,ID,NPT,49,10,1,KKORAK,NZAD,NUMZAD,KONT)
      CALL STAU09(TT1,CORD,ID,NPT,49,9,1,KKORAK,NZAD,NUMZAD,KONT,TT1)
      CALL STAU09(TT1,CORD,ID,NPT,49,1,1,KKORAK,NZAD,NUMZAD,KONT,TT1)
      CALL STAU09(TT1,CORD,ID,NPT,49,11,1,KKORAK,NZAD,NUMZAD,KONT,TT1)
      CALL STAU09(VECTJ,CORD,ID,NPT,49,41,2,KKORAK,NZAD,NUMZAD,KONT,TT1)
      CALL STAU09(GRADJN,CORD,ID,NPT,49,51,2,KKORAK,NZAD,NUMZAD,KONT,
     1            TT1)
      CALL STAU09(FZAPR,CORD,ID,NPT,49,61,2,KKORAK,NZAD,NUMZAD,KONT,TT1)
      ENDIF

      KKORAK=KKORAK+1
      IF (DABS(VREME(NNPER,KORAK+1)).LT.1.D-10) THEN
        GOTO 600
      ELSE
        GOTO 35
      ENDIF
  600 CONTINUE
      IF(NSTAC.EQ.1.AND.KKORAK.EQ.2) 
     +   CALL IZBACI(TT1,CORD,ID,NPT,NEL,NET,IOSA,NDIM,IIZLAZ)
      RETURN      
      END
C=======================================================================
C=======================================================================

      SUBROUTINE IZLLST(ID,TT1,KKORAK,VECTJ,VVREME,FZAPR,VG,GG,
     1                  NEL,KOTES)
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /NDESUK/ NDES,IDPRIT,IFORM
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /PENALL/ PENALT,PRESS,IBRGT

       COMMON /SRPSKI/ ISRPS

      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT,GAMA
      COMMON /DODAT/ NDIMM
      COMMON /SILEZAP/ FZAP(2,9)
      COMMON /ICITANJE/INPT
      
      DIMENSION D4(4),TT1(*),VECTJ(2,*),ID(1,*),N4(4)
      DIMENSION FZAPR(2,*),NRED(4),NREDK(4),NEL(NDIMM,*)
      DIMENSION VG(3,NET,*),GG(3,NET,*)
      DATA NRED/4,2,1,3/ 
      DATA NREDK/9,3,1,7/ 
C
CS  POTENCIJAL U CVOROVIMA 
CE  NODE POTENTIAL 
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2001)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6001)
      WRITE (IIZLAZ,*)'STEP NUMBER  ',KKORAK

      KG = 0
   20 DO 50 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 100
      N4(I) = KG
      JJEDN=ID(1,KG)
      IF(JJEDN.EQ.0) THEN
        D4(I) = 0.
      ELSE
        D4(I) = TT1(JJEDN)
      ENDIF
   50 CONTINUE
C
  100 I = I - 1
      IF(IFORM.EQ.1) GO TO 88
      IF(INPT.EQ.1)THEN
       WRITE(IIZLAZ,5022) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      ENDIF
      GO TO 77
C
   88 IF(INPT.EQ.1)THEN
       WRITE(IIZLAZ,5023) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
      ENDIF
   77 IF(KG.LT.NPT) THEN
      GO TO 20
      ENDIF
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2028)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6028)
C
CS  GUSTINE STRUJE U PRAVCU X1     
CE  CURRENT DESNITIES IN DIRECTION X1
C

      KG = 0
  620 DO 650 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 1600
      N4(I) = KG
         D4(I) = VECTJ(1,KG)
  650 CONTINUE
C
 1600 I = I - 1
      IF(IFORM.EQ.1) GO TO 688
      IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5022) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      ENDIF
      GO TO 677
C
  688 IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5023) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
      ENDIF
  677 IF(KG.LT.NPT) THEN
      GO TO 620
      ENDIF

      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2039)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6039)
C
CS  GUSTINE STRUJE U PRAVCU X2     
CE  CURRENT DESNITIES IN DIRECTION X2
C


      KG = 0
  720 DO 750 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 1700
      N4(I) = KG
         D4(I) = VECTJ(2,KG)
  750 CONTINUE
C
 1700 I = I - 1
      IF(IFORM.EQ.1) GO TO 788
      IF(INPT.EQ.1)THEN
       WRITE(IIZLAZ,5022) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      ENDIF
      GO TO 777
C
  788 IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5023) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
      ENDIF
  777 IF(KG.LT.NPT) THEN
      GO TO 720
      ENDIF
C
CS  BRZINE
C
      NGAUS=IBRGT*IBRGT
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2059)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6059)
      DO I=1,NET
         WRITE(IIZLAZ,5000) I
      DO J=1,NDIM
         NC=NEL(J,I)
         NG=NRED(J)
	   IF(KOTES.EQ.1) NG=NREDK(J)
         IF(INPT.EQ.1) THEN
          WRITE(IIZLAZ,5027) NC,(VG(K,I,NG),K=1,2)
         ELSE
          WRITE(IIZLAZ,5007) NC,(VG(K,I,NG),K=1,2)
         ENDIF
      ENDDO
      ENDDO
C
CS  GRADIJENTI
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2069)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6069)
      DO I=1,NET
         WRITE(IIZLAZ,5000) I
      DO J=1,NDIM
         NC=NEL(J,I)
         NG=NRED(J)
	   IF(KOTES.EQ.1) NG=NREDK(J)
         IF(INPT.EQ.1) THEN
          WRITE(IIZLAZ,5027) NC,(GG(K,I,NG),K=1,2)
         ELSE
          WRITE(IIZLAZ,5007) NC,(GG(K,I,NG),K=1,2)
         ENDIF
      ENDDO
      ENDDO
C
CS  ZAPREMINSKE SILE
CE  BODY FORCES
C
C	WRITE(25,*)'C  N,    FX,         FY,      FZ'
	IF (GAMA.GT.0.D0) THEN
	   DO NN=1,NPT
             IF(INPT.EQ.1) THEN
              WRITE(25,5027) NN,(FZAPR(KK,NN),KK=1,2),0.000
             ELSE
              WRITE(25,5007) NN,(FZAPR(KK,NN),KK=1,2),0.000
             ENDIF
	   ENDDO
	ENDIF
C
CS  PROTOK KROZ KONTURU     
CE  CURRENT FLUX WITHIN CONTOUR
C
      IF (NKONT.GT.0) THEN
       QK=0.
       QM=0.
       DO L=1,NKONT
          QK=QK+QUK(L)
          QM=QM+QUM(L) 
       ENDDO
       IF(INDSC.EQ.0) THEN
        WRITE(IIZLAZ,7023) VVREME
        DO L=1,NKONT
         WRITE(IIZLAZ,5007) L,QUK(L),QUM(L),QUK(L)+QUM(L)
        ENDDO
        WRITE(IIZLAZ,5008) QK,QM,QK+QM
        WRITE(21,7023) VVREME
        DO L=1,NKONT
         WRITE(21,5007) L,QUK(L),QUM(L),QUK(L)+QUM(L)
        ENDDO 
        WRITE(21,5008) QK,QM,QK+QM
        ENDIF
      ENDIF
C
      RETURN
 5000 FORMAT(/' ELEMENT',I10)     
 5002 FORMAT(4(I5,1PD13.5))
 5003 FORMAT(4(I5,3X,F10.3))
 5005 FORMAT(F10.3,10(1PE13.5))
 5007 FORMAT(I6,3(1PE13.5))
 5008 FORMAT('UKUPNO',3(1PE13.5))
 5022 FORMAT(4(I10,1PD13.5))
 5023 FORMAT(4(I10,3X,F10.3))
 5027 FORMAT(I10,3(1PE13.5))
 2001 FORMAT(/'    U K U P N I    P O T E N C I J A L                  '
     1,'             '//
     1' CVOR  POMERANJE   CVOR  POMERANJE   CVOR  POMERANJE   CVOR POMER 
     1NJE   '/' BROJ',3(14X,'BROJ'))
 2028 FORMAT(/'     C V O R N E   B R Z I N E   U   X 1  P R A V C U '//
     1' CVOR    BRZINA    CVOR    BRZINA    CVOR    BRZINA    CVOR   BRZ       
     1INA   '/' BROJ',3(14X,'BROJ'))
 2039 FORMAT(/'     C V O R N E   B R Z I N E   U   X 2  P R A V C U '//
     1' CVOR    BRZINA    CVOR    BRZINA    CVOR    BRZINA    CVOR   BRZ       
     1INA   '/' BROJ',3(14X,'BROJ'))
 2059 FORMAT(//'     C V O R N E   B R Z I N E'/
     1' CVOR    BRZINA VX    BRZINA VY    BRZINA VZ')
 2069 FORMAT(//'     C V O R N I   G R A D I J E N T I'/
     1' CVOR GRADIJENT GX GRADIJENT GY GRADIJENT GZ')
 6001 FORMAT(/'    N O D A L    T O T A L     P O T E N C I A L        '
     1,'                 '//
     1' NODE  DISPLAC.    NODE  DISPLAC.    NODE  DISPLAC.    NODE DISPL
     1AC.   '/'  No.',3(15X,'No.'))
 6028 FORMAT(/'    N O D A L    V E L O C I T I E S   I N   D I R E C T'
     1,' I O N   X 1 '//
     1' NODE  VELOCITY    NODE  VELOCITY    NODE  VELOCITY    NODE VELOC
     1ITY   '/'  No.',3(15X,'No.'))
 6039 FORMAT(/'    N O D A L    V E L O C I T I E S   I N   D I R E C T'
     1,' I O N   X 2 '//
     1' NODE  VELOCITY    NODE  VELOCITY    NODE  VELOCITY    NODE VELOC
     1ITY   '/'  No.',3(15X,'No.'))
 6059 FORMAT(//'    N O D A L    V E L O C I T I E S'/
     1' NODE  VELOCITY VX  VELOCITY VY  VELOCITY VZ')
 6069 FORMAT(//'    N O D A L    G R A D I E N T S  '/
     1' NODE  GRADIENT GX  GRADIENT GY  GRADIENT GZ')
 7023 FORMAT(//'  F L U X   W I T H I N  C O N T O U R, TIME =',1PD13.5/
     1' KONTURA     ULAZNI       IZLAZNI      UKUPNI PROTOK-2')
      END
C======================================================================
C
C RACUNANJE PROTOKA DUZ KONTURE UNUTAR MREZE
C
C=========================================================================
      SUBROUTINE RKONTUR(TT1,ID,KONT,NEL,AKONST,CORD)
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /TRENUT/ TT21,H,HP,ZVHX,ZVHY,HV2,HV3,ZVXT,ZVYT,DETJ,
     1DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,B,BT,GRADH,GRADHT,TT210,
     2DEBLJ,ZHP,H2N,H2NP,GH2NT,H2NN
      COMMON /STAC/ NSTAC,NSLP,IBROJ
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /DODAT/ NDIMM
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL
C
      DIMENSION R(3,3),S(3,3),W(3,3)
      DIMENSION H(9),ZVHX(9),ZVHY(9),HP(4)
      DIMENSION H2N(2,18),H2NP(2,4),ZHP(4,2),GH2NT(18,2)
      DIMENSION B(4,18),BT(18,4),GRADHT(4,2),GRADH(2,4)
      DIMENSION X(9),Y(9),TT21(44),TT210(44)
      DIMENSION H2NN(2,18)
C
      DIMENSION KONT(5,MAXLIN,*),NEL(NDIMM,*),ID(1,*)
      DIMENSION AFIFI(9,9),AM1(9,9),AK(9,9),TT1(*),CORD(3,*)
      DIMENSION AKONST(3,5,*),QPR(9)
C
      R(3,1)=-0.7745966692415
      R(3,2)=0.0
      R(3,3)=0.77459666924148
      S(3,1)=-0.7745966692415
      S(3,2)=0.0
      S(3,3)=0.77459666924148
      W(3,1)=0.55555555555556
      W(3,2)=0.88888888888889
      W(3,3)=0.55555555555556


      R(2,1)=-0.57735026918963
      R(2,2)=0.5773502691896
      S(2,1)=-0.57735026918963
      S(2,2)=0.5773502691896
      W(2,1)=1.000000000000000
      W(2,2)=1.000000000000000

      R(1,1)=0.0
      S(1,1)=0.0
      W(1,1)=2.0
C
C ODEDIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
C      IBRGT=3
C      IF (NDIM.EQ.4) IBRGT=2
      IBRGT=2
C
      DO 723 KK=1,NKONT
      QUK(KK)=0.D0
      QUM(KK)=0.D0
C
	DO 134 II=1,LIN(KK)
C
      NBREL=KONT(1,II,KK)
      NMAT=NEL(NDIM+1,NBREL)
      NEMA=NEL(NDIM+2,NBREL)

      IF (NMAT.LE.0) GOTO 134
      IF (NEMA.GT.0) GOTO 134
C
      CALL CLEARD(TT21,NDIM)
      CALL PREB(TT21,TT1,ID,NEL,NBREL)
C
      DO KLM=1,NDIM
       X(KLM)=CORD(1,NEL(KLM,NBREL))
       Y(KLM)=CORD(2,NEL(KLM,NBREL))
      ENDDO

	  ICV1=KONT(2,II,KK)
	  ICV2=KONT(3,II,KK)
C
      DO 163 K=1,NDIM
      DO 162 N=1,NDIM
        AM1(K,N)=0.D0    
        AFIFI(K,N)=0.D0    
  162 CONTINUE
  163 CONTINUE
C
      DO 180 I=1,IBRGT
      DO 170 J=1,IBRGT
C
      CALL INTERP(R(IBRGT,I),S(IBRGT,J),0,0)
C  
      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ*DEBLJ
C
C OKVASENOST U GAUS TACKI
        PR=0.D0
        DO IK=1,NDIM
           PR=PR+H(IK)*(TT21(IK)-CORD(IOSA,NEL(IK,NBREL)))
        ENDDO
      IF ( PR.LT.0.D0 .AND. INDFS.EQ.1 )THEN
         AKX=(1.D0-ZASIC)*AKONST(3,1,NMAT)
         AKY=(1.D0-ZASIC)*AKONST(3,2,NMAT)
         AKS=(1.D0-ZASIC)*AKONST(3,3,NMAT)
      ELSE
         AKX=AKONST(3,1,NMAT)
         AKY=AKONST(3,2,NMAT)
         AKS=AKONST(3,3,NMAT)
      ENDIF
C
      DO 165 K=1,NDIM
      DO 164 N=1,NDIM
C
         AM1(K,N)=AM1(K,N)+AKS*H(K)*H(N)*WDT
         AFIFI(K,N)=AFIFI(K,N)+(AKX*ZVHX(K)*ZVHX(N)+
     1                          AKY*ZVHY(K)*ZVHY(N))*WDT
C
  164 CONTINUE
  165 CONTINUE
c
  170 CONTINUE
  180 CONTINUE
C
      DO I=1,NDIM
      DO J=1,NDIM
        AK(I,J)=AFIFI(I,J)
      ENDDO
      ENDDO
C
      DO IK=1,NDIM
         QPR(IK)=0.D0
      DO J=1,NDIM
         QPR(IK)=QPR(IK)+AK(IK,J)*TT21(J)
      ENDDO
      ENDDO
C
      DO J=1,NDIM
      IF(NEL(J,NBREL).EQ.ICV1.OR.NEL(J,NBREL).EQ.ICV2) THEN
        IF(QPR(J).LT.0.) THEN
          QUM(KK)=QUM(KK)+QPR(J)
        ELSE
          QUK(KK)=QUK(KK)+QPR(J)
        ENDIF
      ENDIF
      ENDDO
C
  134 CONTINUE
C
  723 CONTINUE
C
      RETURN
      END
C=======================================================================
C
C RACUNANJE PROTOKA DUZ KONTURE UNUTAR MREZE
C
C=========================================================================
      SUBROUTINE RKONTV(VECTJ,KONT,NEL,CORD,VG)
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /STAC/ NSTAC,NSLP,IBROJ
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /DODAT/ NDIMM
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /AXIS/ HH,RIZ,INDAX
C
      DIMENSION R12(3),EN(3),V(3)
C
      DIMENSION KONT(5,MAXLIN,*),NEL(NDIMM,*)
      DIMENSION VECTJ(2,*),CORD(3,*)
      DIMENSION NRED(4)
      DIMENSION VG(3,NET,*)
      DATA NRED/4,2,1,3/ 
C
      ONE=1.
      PI=4.*DATAN(ONE)
C
      DO 723 KK=1,NKONT
        QUK(KK)=0.D0
        QUM(KK)=0.D0
C
	DO 134 II=1,LIN(KK)
C
          NBREL=KONT(1,II,KK)
          NMAT=NEL(NDIM+1,NBREL)
          NEMA=NEL(NDIM+2,NBREL)

          IF (NMAT.LE.0) GOTO 134
          IF (NEMA.GT.0) GOTO 134
C
	  ICV1=KONT(2,II,KK)
	  ICV2=KONT(3,II,KK)

          DO KLM=1,NDIM
            IF(ICV1.EQ.NEL(KLM,NBREL)) IC1=KLM
            IF(ICV2.EQ.NEL(KLM,NBREL)) IC2=KLM
          ENDDO

          DEBLJ=1.
          IF(INDAX.EQ.1)
     +          DEBLJ=(CORD(1,NEL(IC1,NBREL))+CORD(1,NEL(IC2,NBREL)))/2.
C
          R12(1)=CORD(1,ICV1)-CORD(1,ICV2)
          R12(2)=CORD(2,ICV1)-CORD(2,ICV2)
          R12(3)=0.
          RI=DSQRT(R12(1)*R12(1)+R12(2)*R12(2)+R12(3)*R12(3))
          EN(1)=-R12(2)/RI
          EN(2)=R12(1)/RI
          EN(3)=0.
          V(1)=(VG(1,NBREL,NRED(IC1))+VG(1,NBREL,NRED(IC2)))/2.
          V(2)=(VG(2,NBREL,NRED(IC1))+VG(2,NBREL,NRED(IC2)))/2.
C         V(1)=(VECTJ(1,ICV1)+VECTJ(1,ICV2))/2.
C         V(2)=(VECTJ(2,ICV1)+VECTJ(2,ICV2))/2.
          V(3)=0.
          VN=EN(1)*V(1)+EN(2)*V(2)+EN(3)*V(3)
	  VINT=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)
	  IF(VINT.GT.1.D-10) THEN
            VI=DSQRT(VINT)
            AL=DACOS(VN/VI)*180./PI
	  ENDIF
	  QN=-VN*RI*DEBLJ
	  IF(QN.LT.0.) THEN
            QUM(KK)=QUM(KK)+QN
          ELSE
            QUK(KK)=QUK(KK)+QN
          ENDIF
C
  134   CONTINUE
C
  723 CONTINUE
C  
      IF(NKONT.GT.0) THEN
       QK=0.
       QM=0.
       DO L=1,NKONT
          QK=QK+QUK(L)
          QM=QM+QUM(L) 
       ENDDO
       IF(INDSC.EQ.0) THEN
        WRITE(IIZLAZ,7023) VVREME
        DO L=1,NKONT
         WRITE(IIZLAZ,5007) L,QUK(L),QUM(L),QUK(L)+QUM(L)
        ENDDO
        WRITE(IIZLAZ,5008) QK,QM,QK+QM
        WRITE(21,7023) VVREME
        DO L=1,NKONT
         WRITE(21,5007) L,QUK(L),QUM(L),QUK(L)+QUM(L)
        ENDDO 
        WRITE(21,5008) QK,QM,QK+QM
        ENDIF
      ENDIF
 5007 FORMAT(I6,3(1PE13.5))
 5008 FORMAT('UKUPNO',3(1PE13.5))
 7023 FORMAT(//'  F L U X   W I T H I N  C O N T O U R, TIME =',1PD13.5/
     1' KONTURA     ULAZNI       IZLAZNI      UKUPNI PROTOK-3')
C
      RETURN
      END
