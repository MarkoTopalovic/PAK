C$DEBUG      
C=========================================================================
      SUBROUTINE VRACN3D(TT1,SILE,
     1NZAD,ZADVRE,NGPSIL,MAXA,SKEF,SKEFN,KONT,
     2FZAPR,VREME,TABF,TT10,UBRZ,UBRZ0,AK,
     3VECTJ,IVECT,POVSIL,GRADJN,
     4ITFMAX,AKONST,NASLOV,ICUR,VG,GG,KOJK,ISNUMER)
      USE PPR
      USE STIFFNESS
      USE NODES
      USE ELEMENTS
      USE KONTURE
      USE PREDISCRIBED
      USE zapresil
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'paka.inc'
      INCLUDE 'mpif.h'
c      COMMON A(17000)
c      REAL A
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DETJS,DETJ,NBREL
      COMMON /KRITER/ IDOKL,MAXCUR
      COMMON /TREN13/ FPS1,FPS2,FPS3,FPS4
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /LSK1/ LSK
      COMMON /JEDANP/ INDJED,NBRF,NGL,INDTLO
      COMMON /FLUX/ KISA
      COMMON /TACNOS/ EPSTR,MAXIT,NJRAP
      COMMON /ZADPO/ LNZADJ
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /SCRATC/ ISCRC
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /NDESUK/ NDES,IDPRIT,IFORM
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT,GAMA
      COMMON /VREPER/ NPER,NTABFT
      COMMON /PROMEN/ NJUTN,INDOPT,INDPRO
      COMMON /DODAT/ NDIMM
      COMMON /DYNAM/ NDIN
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /STAC/ NSTAC,NSLP,IBROJ
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /VDP/ DT,NKORP,NN,NZAV
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
      COMMON /ICITANJE/ INPT
C
      COMMON /SILZP3/ FZAP(3,21)
C
      COMMON /PRIKAZ/ INDSC,IZPOT
C
      COMMON /PIJEZO/ CPOR(3,1000),NPIJEZ(20,100),NODP(100),NPIJ,NPOR,
     1                NPORCV(1000),NPOREL(1000),NEPOR,NPORE,NEPIJ
      COMMON /PORCEL/ RMIN(34),SMIN(34),TMIN(34)
      COMMON /NXNAST/ NXNASTRAN
      COMMON /SOLVER/ ISOLVER
      COMMON /NELDIM/ NDIMEL
      COMMON /DJERDAP/ IDJERDAP,ISPRESEK
      COMMON /STAMPAZT/ NPRINT
C
      CHARACTER*80 NASLOV
      DIMENSION TT1(*),SILE(*),ZADVRE(*)
      DIMENSION VG(3,NET,*),GG(3,NET,*)
      DIMENSION TT10(*),UBRZ(*),UBRZ0(*)
!       DIMENSION ID(1,*),NZAD(3,*),NGPSIL(12,*),MAXA(*),CORD(3,*)
      DIMENSION NZAD(3,*),NGPSIL(12,*),MAXA(*)
      DIMENSION SKEF(NDES,*),AKONST(3,5,*),SKEFN(*)
      DIMENSION POVSIL(4,*),PR1(21),ITFMAX(*)
C
      DIMENSION KONT(9,MAXLIN,*),KOJK(*)

      DIMENSION TT21(44),TT210(44)
!       ,QUK(1000),QUM(1000)
      DIMENSION R(3,3),S(3,3),T(3,3),W(3,3)
      DIMENSION RK(3,3),SK(3,3),TK(3,3),WK(3,3)
      DIMENSION FZAPR(3,*),VECTJ(3,*),GRADJN(3,*)
      DIMENSION IVECT(*),LM2(44),ICUR(*)
      DIMENSION PERM(3,4)
      DIMENSION VREME(NPER,950)
      DIMENSION TABF(2,NTABFT,*)
      DIMENSION AM1(21,21),AK(NDES,*)

      DIMENSION RS2(21)
      DIMENSION RS3(21)
      DIMENSION F36(44)
C      DIMENSION PORNIEL(100000,8)
C
      DIMENSION AFIFI1(21,21),AFIFI(21,21)
      DIMENSION PUS(21,21),AJX(20),AJY(20),AJZ(20)
      DIMENSION AJX1(20),AJY1(20),AJZ1(20),NRED(8),NREDK(8)
      DIMENSION NRED1(3),NREDT4(4),NREDT10(10),NRED26(6),NRED28(8)
      DIMENSION HT(7,3),CKL(3,8),TTE(2,3),WTET(3)
!       promenjena duzina NRED sa 8 na 10 i dodata 2 clana
      DATA NRED/8,4,2,6,7,3,1,5/  
      DATA NREDK/27,9,3,21,25,7,1,19/
      DATA NRED1/2,1,12/
      DATA NRED26/1,2,3,12,23,31/
      DATA NRED28/1,2,3,4,12,23,34,41/
      DATA NREDT4/1,2,3,4/
      DATA NREDT10/1,2,3,4,12,23,13,14,24,34/
      integer ierr, myid
      integer*8 LM2
!       ,NGPSIL
      integer Dtime(8)
C
      CHARACTER*3 STATTT
      DIMENSION NDJEL(21),NP3D1(18)
      logical OLDNEWW

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      IF (myid.ne.0) goto 1235
  
      if (.not.allocated(PORNIEL)) allocate(PORNIEL(NET,10),STAT=istat) 
      if (.not.allocated(SILEL)) allocate(SILEL(NET),STAT=istat) 
C
C      IDJERDAP=1
!        CITANJE SE RADI PRE STAMPE 
!       IF(IDJERDAP.GE.1)THEN
!       CALL VICLEAR(NP3D1,5)
! C   podprogram za citanje cvorova u kojima se stampaju filtracione sile
!               CALL VDJERDAPREADT(NP3D1,ISNUMER,II)
!        ELSEIF(IDJERDAP.EQ.-1)THEN
! C   stampaju se filtracione sile za sve lamele/sekcije
!       CALL VICLEAR(NP3D1,18)
!           DO II=1,18
!              CALL VDJERDAPREADT(NP3D1,ISNUMER,II)
!           ENDDO
!        ENDIF
C INDIKATOR ZA OSU TEZINE
 9998   IF (IOSA.EQ.0) IOSA=3

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
       CALL VINIDOK(NZAD,ICUR,CORD)

C UTVRDJIVANJE KOD KOLIKO ELEMENATA SE POJAVLJUJE JEDAN CVOR
      CALL VMVECT(IVECT)

      DO I=1,JEDN
       TT1(I)=0.D0
       UBRZ(I)=0.D0
      ENDDO      
       
C DEFINISANJE CVOROVA ZA PROTOKE
      DO KK=1,NKONT
	DO II=1,LIN(KK)
	  KOJK(KONT(2,II,KK))=KOJK(KONT(2,II,KK))+1
	  KOJK(KONT(3,II,KK))=KOJK(KONT(3,II,KK))+1
	  KOJK(KONT(4,II,KK))=KOJK(KONT(4,II,KK))+1
	  KOJK(KONT(5,II,KK))=KOJK(KONT(5,II,KK))+1
	  KOJK(KONT(6,II,KK))=KOJK(KONT(6,II,KK))+1
	  KOJK(KONT(7,II,KK))=KOJK(KONT(7,II,KK))+1
	  KOJK(KONT(8,II,KK))=KOJK(KONT(8,II,KK))+1
	  KOJK(KONT(9,II,KK))=KOJK(KONT(9,II,KK))+1
        ENDDO
      ENDDO	  
C
      R(3,1)=-0.7745966692415
      R(3,2)=0.0
      R(3,3)=0.77459666924148
      S(3,1)=-0.7745966692415
      S(3,2)=0.0
      S(3,3)=0.77459666924148
      T(3,1)=-0.7745966692415
      T(3,2)=0.0
      T(3,3)=0.77459666924148
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
      TK(3,1)=-1.0
      TK(3,2)=0.0
      TK(3,3)=1.0
      WK(3,1)=0.33333333333333
      WK(3,2)=1.33333333333333
      WK(3,3)=0.33333333333333
C 
      R(2,1)=-0.5773502691896
      R(2,2)=0.57735026918963
      S(2,1)=-0.5773502691896
      S(2,2)=0.57735026918963
      T(2,1)=-0.5773502691896
      T(2,2)=0.57735026918963
      W(2,1)=1.000000000000000
      W(2,2)=1.000000000000000

      R(1,1)=0.0
      S(1,1)=0.0
      T(1,1)=0.0
      W(1,1)=2.0
      
      PI=4.D0*DATAN(1.D0)
            
C===========================================================================
      IAXIS=6
C
      INULTO=0
       KKORAK=1
      VVREME=0.D0
C indikator LINTE=1 - linearan proracun
!       LINTE=1
C=========================================================================
          CALL DATE_AND_TIME(VALUES=Dtime)
          WRITE(*,*) 'vreme pre petlje po koracima', (Dtime(i),i=5,7)
          WRITE(3,*) 'vreme pre petlje po koracima', (Dtime(i),i=5,7)
c petlja po periodima
1235  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NPER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      DO 600 NNPER=1,NPER
      KORAK=0
c petlja po koracima !topalovic ovde pocinje 
   35 KORAK=KORAK+1
      if(myid.ne.0) goto 3456
C
C=========================================================================   
C  RAD SA STACIONARNIM STANJEM U PRVOM KORAKU I NESTACIONARNIM U OSTALIM
C  ( GEOLOSKI PROFILI ZA CERNE)
C
      IF (KKORAK.EQ.1) then
        NNSTAC=NSTAC
C       WRITE(19) KKORAK
      endif
cz uslov da prvi korak radi stacionarno
      IF (KKORAK.EQ.1.AND.NNSTAC.EQ.2) NSTAC=1
      IF (KKORAK.GT.1.AND.NNSTAC.EQ.2) NSTAC=0
C
cz      IF (KKORAK.GT.0) THEN
      TIME=VREME(NNPER,KORAK)
      VVREME=VVREME+TIME
cz      ELSE
cz       TIME=1.D7
cz       KORAK=KORAK-1
cz      ENDIF
C     STAMPANJE PRITISAKA U NEU
      ISIFRA=9000001
      INEU=49
      WRITE(INEU,9000) KKORAK,ISIFRA
 9000 FORMAT(I5,',',I7,',1,'/
     1'Pressure Face 1 Set 1'/
     2'0.0,0.0,0.0,'/
     3'0,0,0,0,0,0,0,0,0,0,'/
     4'0,0,0,0,0,0,0,0,0,0,'/
     5'0,'/
     6'0,0,3,8,'/
     7'0,0,1,')
C     KRAJ ZAGLAVLJA

      DO I=1,NPT 
       DO J=1,3
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
       
      if(myid.ne.0) goto 3456
      
	IF (NBLOCK.GT.1) THEN
       OPEN (ISCRC,FILE='ZSKLIN',FORM='UNFORMATTED',STATUS='UNKNOWN')
       REWIND ISCRC
       ELSE
        KC=NWK
      ENDIF
C===========================================================================
C ZA CVOROVE NA LINIJI PROCURIVANJA INICIJALIZACIJA ZA NESTACIONARNO
C STRUJANJE - START JE UVEK OD PRVOG NAJNIZEG CVORA PA NA VISE
C      IF (ITER.EQ.0.AND.MAXCUR.NE.0) IDOKL=2
      IF (ITER.EQ.0.AND.MAXCUR.NE.0) IDOKL=MAXCUR
C===========================================================================
      IF(NXNASTRAN.EQ.0) THEN
	IF (MAXCUR.GT.0) THEN
          if(ISOLVER.NE.-11) CALL VCLEAR(A(LSK),KC)
        ELSE
          IF (ITER.EQ.0.and.ISOLVER.NE.-11) CALL VCLEAR(A(LSK),KC)
        ENDIF
      ENDIF
       
	DO I=1,JEDN
	 SILE(I)=0.D0
	 UBRZ(I)=0.D0
	ENDDO 

C
      IF((KKORAK.EQ.1).AND.(ITER.EQ.0))THEN
!       
        IF ((ISOLVER.EQ.-11).OR.(ISOLVER.EQ.1)) THEN
C         MUMPS
          CALL sparseassembler_init(1)
        ENDIF
      ENDIF     
C GLAVNA PETLJA PO ELEMENTIMA
      IF(NXNASTRAN.EQ.1) GO TO 999
      DO 400 NBREL=1,NET
      IF(elemtip(NBREL).gt.100) THEN
        NTYPE=elemtip(NBREL)/100
        NTIMES=100
      ELSE
        NTYPE=elemtip(NBREL)/10
        NTIMES=10
      ENDIF
      NDIMEL=elemtip(NBREL)-NTIMES*NTYPE
C ODEDIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
      IBRGT=3
      IF (NDIMEL.LE.8) IBRGT=2

      NMAT=NEL(NDIM+1,NBREL)
      NEMA=NEL(NDIM+2,NBREL)

      IF (NMAT.LE.0.OR.NEMA.GT.0) THEN
         WRITE(3,*) 'NBREL,NMAT,NEMA',NBREL,NMAT,NEMA
         WRITE(*,*) 'NBREL,NMAT,NEMA',NBREL,NMAT,NEMA
!          STOP
      ENDIF
      IF (NMAT.LE.0) GOTO 400
      IF (NEMA.GT.0) GOTO 400
      
       PERM(1,1)=AKONST(3,1,NMAT)
       PERM(1,2)=0.
       PERM(1,3)=0.
       PERM(2,2)=AKONST(3,2,NMAT)
       PERM(2,1)=0.
       PERM(2,3)=0.
       PERM(3,3)=AKONST(3,3,NMAT)
       PERM(3,1)=0.
       PERM(3,2)=0.

      IF (PERM(1,1).LT.1.D-15.OR.PERM(2,2).LT.1.D-15.OR.
     1    PERM(3,3).LT.1.D-15) THEN
         WRITE(3,*) 'NBREL,NMAT,NEMA',NBREL,NMAT,NEMA
         WRITE(*,*) 'NBREL,NMAT,NEMA',NBREL,NMAT,NEMA
         WRITE(*,*) 'PERM',PERM(1,1),PERM(2,2),PERM(3,3)
         STOP
      ENDIF
C=========================================================================
      DO 130 KLM=1,NDIMEL
      CK(1,KLM)=CORD(1,NEL(KLM,NBREL))
      CK(2,KLM)=CORD(2,NEL(KLM,NBREL))
      CK(3,KLM)=CORD(3,NEL(KLM,NBREL))
C
      LM2(KLM)=ID(1,NEL(KLM,NBREL))
 130  CONTINUE
    
      if(NTYPE.eq.2) CALL VLCK2D(CK,CKL,TTE,NDIMEL)
C=======================================================================
      CALL VCLEARD(TT21,NDIM)
      CALL VPREB(TT21,TT1,NBREL)
      CALL VCLEARD(TT210,NDIM)
      CALL VPREB(TT210,TT10,NBREL)
 
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
      DO 199 I=1,NDIM     1
      RS2(I)=0.
      RS3(I)=0.
 199  CONTINUE
C
C INTEGRACIJA U GAUSOVIM TACKAMA

C		ZASIC=0.D0

C sneza 11.02.2012
c       DO II=1,NDIM
c         PR1(II)=TT21(II)-CORD(IOSA,NEL(II,NBREL))
c        ENDDO

      NGAUS=0
      ZAPREL=0.
      IBRGT1=IBRGT
      IBRGT2=IBRGT
      if(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) THEN
        IBRGT2=1
        IBRGT1=4
        IF(INDJOT.GT.0) IBRGT1=INDJOT
      ENDIF
      IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) IBRGT2=1
      IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) IBRGT1=1
      IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
        IBRGT2=1
        IBRGT1=3
      ENDIF     
      DO 180 I=1,IBRGT1
      DO 170 J=1,IBRGT2
      DO 160 L=1,IBRGT2
      NGAUS=NGAUS+1
      IF(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) THEN
       CALL VGETRST(RT,ST,TT,WTET,NGAUS,IBRGT1)
       CALL VINTERPTETRA(RT,ST,TT,0,0,NDIMEL)
       WDT=WTET(1)*WTET(2)*WTET(3)*DETJ
      ELSE IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) THEN
        CALL VINTERP1D(R(IBRGT,I),0,0,NDIMEL,AJS)  
        WDT=W(IBRGT,I)*thick(NBREL)*AJS
      ELSE IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) THEN
        CALL VINTERP2D48(R(IBRGT,J),S(IBRGT,L),0,0,NDIMEL,CKL,TTE)
        WDT=W(IBRGT,J)*W(IBRGT,L)*DETJ*thick(NBREL)
      ELSE IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
        RT=1./6
        ST=1./6
        WT=1./3
        IF(NGAUS.EQ.2) RT=2./3
        IF(NGAUS.EQ.3) ST=2./3
        CALL VINTERP2D36(RT,ST,0,0,NDIMEL,CKL,TTE)
        WDT = WT**2*thick(NBREL)*DETJ
      ELSE
        CALL VINTER3(R(IBRGT,I),S(IBRGT,J),T(IBRGT,L),0,0)
        WDT=W(IBRGT,I)*W(IBRGT,J)*W(IBRGT,L)*DETJ
      ENDIF
      
      ZAPREL=ZAPREL+WDT
      
Cs ovaj uslov nije isti za stacionarno i nestacionarno
Cs za nestacionarno proverava u svakoj iteraciji 
      IF ( (((ITER.GT.0).AND.(kkorak.EQ.1))
     1.OR.KKORAK.GT.1) .AND.(INDFS.EQ.1)) THEN
C      IF ( (((ITER.GT.0).AND.(NSTAC.EQ.1))
C     1.OR.((KKORAK.GT.0).AND.(NSTAC.EQ.0))) .AND.(INDFS.EQ.1)) THEN
C
       PR=0.D0
C okvasenost u gausovoj tacki
        DO II=1,NDIMEL
         PR=PR+H(II)*(TT21(II)-CORD(IOSA,NEL(II,NBREL)))
        ENDDO

           IF( PR.LT.0.D0 .AND. INDFS.EQ.1 )THEN
            ANPUS=AKONST(3,4,NMAT)*ZASIC
            PERM(1,1)=ZASIC*AKONST(3,1,NMAT)
            PERM(2,2)=ZASIC*AKONST(3,2,NMAT)
            PERM(3,3)=ZASIC*AKONST(3,3,NMAT)
           ELSE
            ANPUS=0.D0
            PERM(1,1)=0.0D0
            PERM(2,2)=0.0D0
            PERM(3,3)=0.0D0
           ENDIF
       ELSE
        PERM(1,1)=0.0D0
        PERM(2,2)=0.0D0
        PERM(3,3)=0.0D0
        ANPUS=0.D0
       ENDIF

      DO 165 K=1,NDIMEL
      DO 164 N=1,NDIMEL
cz S nije podeljeno sa dt; podeljeno sa dt prilikom pakovanja u matricu SKEF
         AM1(K,N)=AM1(K,N)+AKONST(3,4,NMAT)*H(K)*H(N)*WDT
         AFIFI(K,N)=AFIFI(K,N)+(AKONST(3,1,NMAT)*ZVHX(K)*ZVHX(N)+
     1AKONST(3,2,NMAT)*ZVHY(K)*ZVHY(N)+AKONST(3,3,NMAT)*ZVHZ(K)*ZVHZ(N))
     2*WDT
C
         IF (INDFS.EQ.1) THEN
         PUS(K,N)=PUS(K,N)+ANPUS*H(K)*H(N)*WDT
         AFIFI1(K,N)=AFIFI1(K,N)+(PERM(1,1)*ZVHX(K)*ZVHX(N)+
     1PERM(2,2)*ZVHY(K)*ZVHY(N)+PERM(3,3)*ZVHZ(K)*ZVHZ(N))*WDT
         ENDIF
!       if(NBREL.eq.1.AND.I.eq.IBRGT) then
!          WRITE(*,*) 'NBREL',NBREL
!          WRITE(*,*) 'H(K)',H(K)
!          WRITE(*,*) 'H(N)',H(N)
!          WRITE(*,*) 'ZVHX(K)',ZVHX(K)
!          WRITE(*,*) 'ZVHY(K)',ZVHY(K)
!          WRITE(*,*) 'ZVHZ(K)',ZVHZ(K)
!          WRITE(*,*) 'WDT',WDT
!          WRITE(*,*) 'AKONST(3,1,NMAT)',AKONST(3,1,NMAT)
!          WRITE(*,*) 'AKONST(3,2,NMAT)',AKONST(3,2,NMAT)
!          WRITE(*,*) 'AKONST(3,3,NMAT)',AKONST(3,3,NMAT)
!          write(*,*) 'AM1(',K,',',N,')',AM1(K,N)
!          write(*,*) 'AFIFI(',K,',',N,')',AFIFI(K,N)
!       ENDIF
C
  164 CONTINUE
  165 CONTINUE
  160 CONTINUE
  170 CONTINUE
  180 CONTINUE
!       WRITE(*,*) 'ZAPREL',ZAPREL
C INCIJALIZACIJA MATRICE SKEF I F36
      DO 260 I=1,NDES
      DO 258 J=1,NDES
      SKEF(I,J)=0.
 258  CONTINUE
      F36(I)=0.
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
CZ        IVRF=NGPSIL(6,1)
CZ        CALL VTIMFUN (TABF,FK1,VVREME,ITFMAX(IVRF),IVRF)
CZ        CALL VPROCEN(RS2,PR1,FK1,NEL,NBREL,CORD) 
CZ        INDPUN=1
CZ       ENDIF
C============================================================================
C POVRSINSKE SILE
C=================
      IF (INDPUN.EQ.1) GOTO 265

      DO 250 JBRPS=1,MAXSIL
      IF (NBREL.EQ.NGPSIL(1,JBRPS)) THEN
      NODE1=NGPSIL(2,JBRPS)
      NODE2=NGPSIL(3,JBRPS)
      NODE3=NGPSIL(4,JBRPS)
      NODE4=NGPSIL(5,JBRPS)
      NODE5=NGPSIL(6,JBRPS)
      NODE6=NGPSIL(7,JBRPS)
      NODE7=NGPSIL(8,JBRPS)
      NODE8=NGPSIL(9,JBRPS)
      N1=NEL(1,NBREL)
      N2=NEL(2,NBREL)
      N3=NEL(3,NBREL)
      N4=NEL(4,NBREL)
      N5=NEL(5,NBREL)
      N6=NEL(6,NBREL)
      N7=NEL(7,NBREL)
      N8=NEL(8,NBREL)
!       N9=NEL(9,NBREL)
!       N10=NEL(10,NBREL)
      IVRF=NGPSIL(10,JBRPS)
      CALL VTIMFUN (TABF,FK1,VVREME,ITFMAX(IVRF),IVRF)
      FPS1=POVSIL(1,JBRPS)*FK1
      FPS2=POVSIL(2,JBRPS)*FK1
      FPS3=POVSIL(3,JBRPS)*FK1
      FPS4=POVSIL(4,JBRPS)*FK1
C
      NUMGAU=IBRGT
      
      IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) THEN
        DO I=1,NUMGAU
          CALL VINTERP1D(R(IBRGT,I),0,0,NDIMEL,AJS)  
          WDTV=W(IBRGT,I)*thick(NBREL)*AJS
          DO K=1,NDIMEL
          DUM=TT21(K)-CORD(IOSA,NEL(K,NBREL))
          IF ( (((ITER.GT.0).AND.(NSTAC.EQ.1))
     1          .OR.((KKORAK.GT.0).AND.(NSTAC.EQ.0)))
     1          .AND.INDFS.EQ.1.AND.DUM.LT.0.) THEN
          ELSE
            RS2(K)=RS2(K)+H(K)*FPS1*WDTV
          ENDIF
        ENDDO
        ENDDO
        GOTO 250
      ENDIF
      
      IF(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) THEN
         DO NGXT = 1,3
             RT=1./6
             ST=1./6
             WT=1./3
             IF(NGXT.EQ.2) RT=2./3
             IF(NGXT.EQ.3) ST=2./3

            CALL VJACTP33(RT,ST,HT,1,NGPSIL,JBRPS,NDIMEL)
            WDTS = WT*DETJS/2
            DO K=2,8
              DO j=1,NDIMEL
                IF(NGPSIL(K,JBRPS).EQ.NEL(J,NBREL)) THEN
                  DUM=TT21(J)-CORD(IOSA,NEL(J,NBREL))
                  IF ( (((ITER.GT.0).AND.(NSTAC.EQ.1))
     1                .OR.((KKORAK.GT.0).AND.(NSTAC.EQ.0)))
     1                .AND.INDFS.EQ.1.AND.DUM.LT.0.) THEN
                  ELSE
                    RS2(J)=RS2(J)+HT(K-1,1)*FS2*WDTS
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
         ENDDO 
      
      ELSE
      
      DO 225 I=1,NUMGAU
      DO 225 J=1,NUMGAU
C
C KONSTANTNO T (ZETA)
C
      IF((NODE1.EQ.N1.AND.NODE2.EQ.N2.AND.NODE3.EQ.N3.AND.NODE4.EQ.N4).
     1OR.(NODE1.EQ.N2.AND.NODE2.EQ.N3.AND.NODE3.EQ.N4.AND.NODE4.EQ.N1).
     2OR.(NODE1.EQ.N3.AND.NODE2.EQ.N4.AND.NODE3.EQ.N1.AND.NODE4.EQ.N2).
     3OR.(NODE1.EQ.N4.AND.NODE2.EQ.N1.AND.NODE3.EQ.N2.AND.NODE4.EQ.N3).
C
     4OR.(NODE1.EQ.N1.AND.NODE2.EQ.N4.AND.NODE3.EQ.N3.AND.NODE4.EQ.N2).
     5OR.(NODE1.EQ.N4.AND.NODE2.EQ.N3.AND.NODE3.EQ.N2.AND.NODE4.EQ.N1).
     6OR.(NODE1.EQ.N3.AND.NODE2.EQ.N2.AND.NODE3.EQ.N1.AND.NODE4.EQ.N4).
     7OR.(NODE1.EQ.N2.AND.NODE2.EQ.N1.AND.NODE3.EQ.N4.AND.NODE4.EQ.N3))
     8THEN
      CALL VINTER3(R(NUMGAU,I),S(NUMGAU,J),1.D0,3,1)
       GOTO 220
      ENDIF 
C
C KONSTANTNO T (ZETA)
C 
      IF((NODE1.EQ.N5.AND.NODE2.EQ.N6.AND.NODE3.EQ.N7.AND.NODE4.EQ.N8).
     1OR.(NODE1.EQ.N6.AND.NODE2.EQ.N7.AND.NODE3.EQ.N8.AND.NODE4.EQ.N5).
     2OR.(NODE1.EQ.N7.AND.NODE2.EQ.N8.AND.NODE3.EQ.N5.AND.NODE4.EQ.N6).
     3OR.(NODE1.EQ.N8.AND.NODE2.EQ.N5.AND.NODE3.EQ.N6.AND.NODE4.EQ.N7).
C
     5OR.(NODE1.EQ.N5.AND.NODE2.EQ.N8.AND.NODE3.EQ.N7.AND.NODE4.EQ.N6).
     6OR.(NODE1.EQ.N8.AND.NODE2.EQ.N7.AND.NODE3.EQ.N6.AND.NODE4.EQ.N5).
     7OR.(NODE1.EQ.N7.AND.NODE2.EQ.N6.AND.NODE3.EQ.N5.AND.NODE4.EQ.N8).
     8OR.(NODE1.EQ.N6.AND.NODE2.EQ.N5.AND.NODE3.EQ.N8.AND.NODE4.EQ.N7))
     8THEN
      CALL VINTER3(R(NUMGAU,I),S(NUMGAU,J),-1.D0,3,1)
       GOTO 220
      ENDIF 
C 
C KONSTANTNO R (KSI)
C
      IF((NODE1.EQ.N1.AND.NODE2.EQ.N4.AND.NODE3.EQ.N8.AND.NODE4.EQ.N5).
     1OR.(NODE1.EQ.N4.AND.NODE2.EQ.N8.AND.NODE3.EQ.N5.AND.NODE4.EQ.N1).
     2OR.(NODE1.EQ.N8.AND.NODE2.EQ.N5.AND.NODE3.EQ.N1.AND.NODE4.EQ.N4).
     3OR.(NODE1.EQ.N5.AND.NODE2.EQ.N1.AND.NODE3.EQ.N4.AND.NODE4.EQ.N8).
C
     4OR.(NODE1.EQ.N1.AND.NODE2.EQ.N5.AND.NODE3.EQ.N8.AND.NODE4.EQ.N4).
     5OR.(NODE1.EQ.N5.AND.NODE2.EQ.N8.AND.NODE3.EQ.N4.AND.NODE4.EQ.N1).
     6OR.(NODE1.EQ.N8.AND.NODE2.EQ.N4.AND.NODE3.EQ.N1.AND.NODE4.EQ.N5).
     7OR.(NODE1.EQ.N4.AND.NODE2.EQ.N1.AND.NODE3.EQ.N5.AND.NODE4.EQ.N8))
     8THEN
      CALL VINTER3(1.D0,S(NUMGAU,I),T(NUMGAU,J),1,1)
       GOTO 220
      ENDIF 
C 
C KONSTANTNO R (KSI)
C
      IF((NODE1.EQ.N2.AND.NODE2.EQ.N3.AND.NODE3.EQ.N7.AND.NODE4.EQ.N6).
     1OR.(NODE1.EQ.N3.AND.NODE2.EQ.N7.AND.NODE3.EQ.N6.AND.NODE4.EQ.N2).
     2OR.(NODE1.EQ.N7.AND.NODE2.EQ.N6.AND.NODE3.EQ.N2.AND.NODE4.EQ.N3).
     3OR.(NODE1.EQ.N6.AND.NODE2.EQ.N2.AND.NODE3.EQ.N3.AND.NODE4.EQ.N7).
C
     4OR.(NODE1.EQ.N2.AND.NODE2.EQ.N6.AND.NODE3.EQ.N7.AND.NODE4.EQ.N3).
     5OR.(NODE1.EQ.N6.AND.NODE2.EQ.N7.AND.NODE3.EQ.N3.AND.NODE4.EQ.N2).
     6OR.(NODE1.EQ.N7.AND.NODE2.EQ.N3.AND.NODE3.EQ.N2.AND.NODE4.EQ.N6).
     7OR.(NODE1.EQ.N3.AND.NODE2.EQ.N2.AND.NODE3.EQ.N6.AND.NODE4.EQ.N7))
     8THEN
      CALL VINTER3(-1.D0,S(NUMGAU,I),T(NUMGAU,J),1,1)
       GOTO 220
      ENDIF 
C 
C KONSTANTNO S (ETA)
C
      IF((NODE1.EQ.N4.AND.NODE2.EQ.N3.AND.NODE3.EQ.N7.AND.NODE4.EQ.N8).
     1OR.(NODE1.EQ.N3.AND.NODE2.EQ.N7.AND.NODE3.EQ.N8.AND.NODE4.EQ.N4).
     2OR.(NODE1.EQ.N7.AND.NODE2.EQ.N8.AND.NODE3.EQ.N4.AND.NODE4.EQ.N3).
     3OR.(NODE1.EQ.N8.AND.NODE2.EQ.N4.AND.NODE3.EQ.N3.AND.NODE4.EQ.N7).
C
     4OR.(NODE1.EQ.N3.AND.NODE2.EQ.N4.AND.NODE3.EQ.N8.AND.NODE4.EQ.N7).
     5OR.(NODE1.EQ.N4.AND.NODE2.EQ.N8.AND.NODE3.EQ.N7.AND.NODE4.EQ.N3).
     6OR.(NODE1.EQ.N8.AND.NODE2.EQ.N7.AND.NODE3.EQ.N3.AND.NODE4.EQ.N4).
     7OR.(NODE1.EQ.N7.AND.NODE2.EQ.N3.AND.NODE3.EQ.N4.AND.NODE4.EQ.N8))
     8THEN
      CALL VINTER3(R(NUMGAU,I),-1.D0,T(NUMGAU,J),2,1)
       GOTO 220
      ENDIF 
C 
C KONSTANTNO S (ETA)
C
      IF((NODE1.EQ.N1.AND.NODE2.EQ.N2.AND.NODE3.EQ.N6.AND.NODE4.EQ.N5).
     1OR.(NODE1.EQ.N2.AND.NODE2.EQ.N6.AND.NODE3.EQ.N5.AND.NODE4.EQ.N1).
     2OR.(NODE1.EQ.N6.AND.NODE2.EQ.N5.AND.NODE3.EQ.N1.AND.NODE4.EQ.N2).
     3OR.(NODE1.EQ.N5.AND.NODE2.EQ.N1.AND.NODE3.EQ.N2.AND.NODE4.EQ.N6).
C
     4OR.(NODE1.EQ.N2.AND.NODE2.EQ.N1.AND.NODE3.EQ.N5.AND.NODE4.EQ.N6).
     5OR.(NODE1.EQ.N1.AND.NODE2.EQ.N5.AND.NODE3.EQ.N6.AND.NODE4.EQ.N2).
     6OR.(NODE1.EQ.N5.AND.NODE2.EQ.N6.AND.NODE3.EQ.N2.AND.NODE4.EQ.N1).
     7OR.(NODE1.EQ.N6.AND.NODE2.EQ.N2.AND.NODE3.EQ.N1.AND.NODE4.EQ.N5))
     8THEN
      CALL VINTER3(R(NUMGAU,I),1.D0,T(NUMGAU,J),2,1)
       GOTO 220
      ENDIF
 
 220    IF (KISA.EQ.1) THEN
C           DOGRADITI KASNIJE 
C KADA PADA KISA ONDA FLUX DELUJE U PRAVCU -IOSA OSE 
C TAKO DA MORA DA SE PROJEKTUJE NA NORMALU LINIJE PO KOJOJ DELUJE
C          FAK=DABS(CORD(1,NODE1)-CORD(1,NODE2))
C          WDTS=FAK
          WDTS=W(NUMGAU,I)*W(NUMGAU,J)*DETJS
        ELSE
          WDTS=W(NUMGAU,I)*W(NUMGAU,J)*DETJS
        ENDIF

       DO K=1,NDIMEL
        DUM=TT21(K)-CORD(IOSA,NEL(K,NBREL))
        IF ( (((ITER.GT.0).AND.(NSTAC.EQ.1))
     1.OR.((KKORAK.GT.0).AND.(NSTAC.EQ.0)))
     1.AND.INDFS.EQ.1.AND.DUM.LT.0.) THEN
        ELSE
         RS2(K)=RS2(K)+H(K)*FS2*WDTS
        ENDIF
       ENDDO

C     write(3,*) nbrel,i,j,R(IBRGT,I),S(IBRGT,J)
C     write(3,*) 'fs2,wdts,detjs',fs2,wdts,detjs
C     CALL Vwrr(h,ndim,'h   ')
C     CALL Vwrr(rs2,ndim,'rs2 ')
 225  CONTINUE
      ENDIF
 
      ENDIF
 250  CONTINUE

C      ENDIF
C KRAJ PETLJE AKO JE ITER>0
C=======================================================================

C PAKOVANJE POVRSINSKIH SILA SA DESNE STRANE
 265   DO I=1,NDIMEL
         F36(I)=F36(I)+RS2(I)
       ENDDO

C=========================================================================
      DO I=1,NDIMEL
      DO J=1,NDIMEL
       AK(I,J)=0.D0
      ENDDO
      ENDDO
C PAKOVANJE MATRICA AMUU,AMUW,AMWW,AKUU,AKUW,AKWW U MATRICE AMM,AC,AK
      DO I=1,NDIMEL
      DO J=1,NDIMEL
c sneza 10.02.2012
c         AK(I,J)=AFIFI(I,J)-AFIFI1(I,J)
         AK(I,J)=AFIFI(I,J)
        IF ((NSTAC.EQ.0).AND.(KKORAK.GT.0)) THEN
c sneza 10.02.2012
c         AK(I,J)=AK(I,J)+(AM1(I,J)-PUS(I,J))/TIME
         AK(I,J)=AK(I,J)+AM1(I,J)/TIME
        ENDIF
      ENDDO
      ENDDO
C============================================================================
C PAKOVANJE MATRICE AK U MATRICU SKEF
      DO 263 I=1,NDIMEL
      DO 262 J=1,NDIMEL
      SKEF(I,J)=AK(I,J)
  262 CONTINUE
  263 CONTINUE
C============================================================================
C PAKOVANJE Kfifi*fi(i-1) SA DESNE STRANE
      DO I=1,NDIMEL
       DO J=1,NDIMEL
        F36(I)=F36(I)-(AFIFI(I,J)-AFIFI1(I,J))*TT21(J)
       IF ( (NSTAC.EQ.0).AND.(KKORAK.GT.0) ) THEN
         F36(I)=F36(I)-((AM1(I,J)-PUS(I,J))*(TT21(J)-TT210(J)))/TIME
       ENDIF
       ENDDO
      ENDDO
C============================================================================
C VEKTOR PROTOKA PO CVOROVIMA
      DO I=1,NDIMEL
       N=NEL(I,NBREL)
       JJ=ID(1,N)
       IF(JJ.GT.0) THEN

         POT=0.D0
        DO J=1,NDIMEL

         POT=POT+(AFIFI(I,J)-AFIFI1(I,J))*TT21(J)
cz uporediti resenja sa ovim uslovom za nestacionarni protok
cz proveriti resenja za zasic=0.
      	IF (NSTAC.NE.1) THEN
          POT=POT+((AM1(I,J)-PUS(I,J))*(TT21(J)-TT210(J)))/TIME
        ENDIF

        ENDDO
         UBRZ(JJ)=UBRZ(JJ)+POT
       ENDIF
      ENDDO
C PAKOVANJE MATRICE ELEMENTA I VEKTRA DESNE STRANE U SISTEM
       CALL VPSKEFN(SKEF,SKEFN,NDIMEL)
       CALL VSPAKDE(SILE,F36,LM2,NDIMEL)
      IF (ISOLVER.EQ.0) THEN
        IF (MAXCUR.GT.0) THEN
          CALL VSPAKUJ(A(LSK),MAXA,SKEFN,LM2,NDIMEL)
        ELSE
cz ovaj uslov ne valja jer ne radi levu stranu za full newton      
C        IF (ITER.EQ.0.OR.NJRAP.EQ.2) 
          IF (ITER.EQ.0) 
     +                  CALL VSPAKUJ(A(LSK),MAXA,SKEFN,LM2,NDIMEL)
        ENDIF
      ELSEIF ((ISOLVER.EQ.-11).OR.(ISOLVER.EQ.1)) THEN
C      IMSL Sparse OR MUMPS
        IF(ITER.EQ.0) THEN
C         CALL VSPARSE_ELEM_LEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,1)
          CALL sparseassembler_addelemmatrix(NDES,LM2,SKEF)
        ENDIF
      ENDIF
C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE
C=======================================================================
         CALL DATE_AND_TIME(VALUES=Dtime)
         WRITE(*,*) 'vreme kraj petlje po elementima', (Dtime(i),i=5,7)
         WRITE(3,*) 'vreme kraj petlje po elementima', (Dtime(i),i=5,7)
C CVOROVI NA LINIJI PROCURIVANJA       
C=======================================================================
  	 IF (MAXCUR.GT.0) THEN
        CALL VCURI1(CORD,A(LSK),SILE,MAXA,TT1,ID,ICUR)
C       ELSE
C        IF(ITER.EQ.0) CALL VCURI1(CORD,A(LSK),SILE,MAXA,TT1,ID,ICUR)
       ENDIF
C==========================================================================
C MUMPS alociranje matrice krutosti
3456  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
           CALL VSPAKUA(A(LSK),MAXA,A(LSKE),LM2,NDES,0,
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
         CALL DATE_AND_TIME(VALUES=Dtime)
         WRITE(*,*) 'vreme posle sparse kill', (Dtime(i),i=5,7)
         WRITE(3,*) 'vreme posle sparse kill', (Dtime(i),i=5,7)

         IF(KKORAK.EQ.1) THEN
              stiff_n = JEDN
          ENDIF
         ENDIF
        ENDIF
C==========================================================================
	IF (MAXCUR.GT.0) THEN
	 IF(ISOLVER.EQ.0) THEN
          if(myid.eq.0) CALL VZADLEV(A(LSK),MAXA,A(LNZADJ),NUMZAD,NZAD)
         ELSE
          if(myid.eq.0) CALL VSPARS_ZADLEV(A(LNZADJ),NUMZAD,NZAD,
     1                                CORD,VVREME,TABF,NTABFT,ITFMAX)
         ENDIF
         CALL VRESEN(A(LSK),TT10,MAXA,JEDN,1)
        ELSE
cz proveriti i ovde uslov za full newton        
C        IF (ITER.EQ.0.OR.NJRAP.EQ.2) THEN
        IF (ITER.EQ.0) THEN
cz proveriti rad sa blokovima i zadlev
          if(myid.eq.0) then
          IF (NBLOCK.GT.1) THEN
            LLM =LRAD
            LSKE=LLM+100
            CALL VSPAKUA(A(LSK),MAXA,A(LSKE),LM2,NDES,0,
     &          A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,H,NZAD)
            CLOSE (ISCRC,STATUS='KEEP')

          ENDIF
          endif
          
          IF(ISOLVER.EQ.0) THEN
          if(myid.eq.0) CALL VZADLEV(A(LSK),MAXA,A(LNZADJ),NUMZAD,NZAD)
          ELSE
           
           if(myid.eq.0) CALL VSPARS_ZADLEV(A(LNZADJ),NUMZAD,NZAD,
     1                                  CORD,VVREME,TABF,NTABFT,ITFMAX)
          ENDIF
         CALL DATE_AND_TIME(VALUES=Dtime)
          WRITE(3,*) 'pre resen 1', (Dtime(i),i=5,7)
          CALL VRESEN(A(LSK),TT10,MAXA,JEDN,1)
         CALL DATE_AND_TIME(VALUES=Dtime)
         WRITE(3,*) 'posle resen 1', (Dtime(i),i=5,7)
        ENDIF
        ENDIF
       if(myid.ne.0) goto 4567
       TOL=1.0D-7
C=======================================================================
         CALL DATE_AND_TIME(VALUES=Dtime)
         WRITE(3,*) 'vreme pre pisanja desne strane', (Dtime(i),i=5,7)
!         WRITE(3,*) 'NUMZAD'
!         WRITE(3,*)  'Cvor, zad pot', NZAD(1,I),FK1
       DO 410 I=1,NUMZAD
c       write(iizlaz,*)'fk1=',fk1
c       do ii=1,2
c        write(iizlaz,*)TABF(1,1,Ii),TABF(2,1,Ii)
c       enddo
        CALL VTIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,I)),NZAD(3,I))
        POTZADC(I)=FK1
       IF(NZAD(2,I).EQ.6.AND.CORD(IOSA,NZAD(1,I)).GT.FK1) THEN
         POTZADC(I)=0.0
         GOTO 410
       ENDIF
!  za Grancarevo zavisnost na granicama modela
!  desna obala
       IF(NZAD(2,I).EQ.7) Then
!         WRITE(3,*)  '7, pre', NZAD(1,I),CORD(IOSA,NZAD(1,I)),FK1
        FK1=301.0+0.332*FK1
!        WRITE(3,*)  '7, posle',NZAD(1,I),CORD(IOSA,NZAD(1,I)),FK1
        POTZADC(I)=FK1
        if(CORD(IOSA,NZAD(1,I)).GT.FK1) THEN
         POTZADC(I)=0.0
         GOTO 410
        ENDIF
       ENDIF
       IF(NZAD(2,I).EQ.-1) Then
!         WRITE(3,*)  '-1, pre', NZAD(1,I),CORD(IOSA,NZAD(1,I)),FK1
        CALL VPREDF_INTERP(I,NZAD,CORD,FK1,1)
!         WRITE(3,*)  '-1,posle',NZAD(1,I),CORD(IOSA,NZAD(1,I)),FK1
        POTZADC(I)=FK1
        if(CORD(IOSA,NZAD(1,I)).GT.FK1) THEN
         POTZADC(I)=0.0
         GOTO 410
        ENDIF
       ENDIF
!  leva obala - linearno
       IF(NZAD(2,I).EQ.0) Then
!        WRITE(3,*)  '0, pre', NZAD(1,I),CORD(IOSA,NZAD(1,I)),FK1
        CALL VPREDF_INTERP(I,NZAD,CORD,FK1,2)
!         WRITE(3,*)  '0,posle',NZAD(1,I),CORD(IOSA,NZAD(1,I)),FK1
        POTZADC(I)=FK1
        IF(CORD(IOSA,NZAD(1,I)).GT.FK1) THEN
         POTZADC(I)=0.0
         GOTO 410
        ENDIF
       ENDIF
c       write(iizlaz,*)'fk1=',fk1, vvreme,ITFMAX(NZAD(3,I)),NZAD(3,I)
c       do ii=1,2
c        write(iizlaz,*)TABF(1,1,Ii),TABF(2,1,Ii)
c       enddo
       IF (NSTAC.EQ.0) THEN
C STRUJANJE NESTACIONARNO
C===========================================================================
C KONSTANTNI POTENCIJALI
        IF( (NZAD(2,I).EQ.1.OR.NZAD(2,I).EQ.6.OR.
     &        NZAD(2,I).EQ.7.OR.NZAD(2,I).EQ.0.OR.
     &        NZAD(2,I).EQ.-1.OR.
     &        (NZAD(2,I).EQ.3.AND.KKORAK.EQ.1))
     &        .AND.(ITER.EQ.0).AND.(ID(1,NZAD(1,I)).GT.0)) THEN
cz oduzimanje zadate vrednosti od prethodno izracunate             
            SILE(ID(1,NZAD(1,I)))=1.D035*
     1                            ((ZADVRE(I)*FK1)-TT1(ID(1,NZAD(1,I))))
       ENDIF
C ZPREMINSKI IZVOR, PONOR       
       IF ( (ID(1,NZAD(1,I)).GT.0) .AND. (NZAD(2,I).EQ.5)) THEN
          SILE(ID(1,NZAD(1,I)))=ZADVRE(I)*FK1
       ENDIF
cz ovo ima smisla kod izcurivanja na slobodnoj povrsini (PROVERI USLOVE?)
      IF ((KKORAK.EQ.1).AND.(NZAD(2,I).EQ.-1).AND.(INDFS.EQ.1)) 
     1      SILE(ID(1,NZAD(1,I)))=1.0D35*CORD(IOSA,NZAD(1,I))
c     
cz zadavanje promenljivog potencijala
        IF( (NZAD(2,I).EQ.2).AND.(ID(1,NZAD(1,I)).GT.0).AND.
     1     (ITER.EQ.0)) then
cc          if((CORD(IOSA,NZAD(1,I))-FK1).LE.TOL ) THEN
            SILE(ID(1,NZAD(1,I)))=1.D035*(FK1-TT1(ID(1,NZAD(1,I))))
cc         else
c proba sta raditi iznad
c            SILE(ID(1,NZAD(1,I)))=1.D035*(-TT1(ID(1,NZAD(1,I))))
cc         endif
        ENDIF
C===========================================================================
       ELSE
C===========================================================================
C KADA JE STRUJANJE STACIONARNO
C
C  ZA VDP TREBA OVAKO
C
        IF (NZAD(2,I).eq.1.OR.NZAD(2,I).EQ.6.OR.
     &        NZAD(2,I).EQ.7.OR.NZAD(2,I).EQ.0.OR.
     &        NZAD(2,I).EQ.-1) THEN
         IF ( (ID(1,NZAD(1,I)).GT.0) .AND. (ITER.EQ.0) ) THEN
          SILE(ID(1,NZAD(1,I)))=1.0D35*ZADVRE(I)*FK1
c         ELSE 
c          SILE(ID(1,NZAD(1,I)))=0.D0
         ENDIF
        ENDIF
C ZPREMINSKI IZVOR, PONOR       
          IF ( (ID(1,NZAD(1,I)).GT.0) .AND. (NZAD(2,I).EQ.5)) THEN
!               WRITE(*,*) 'ZADVRE(',I,')',ZADVRE(I)
!               WRITE(*,*) 'FK1',FK1
              SILE(ID(1,NZAD(1,I)))=SILE(ID(1,NZAD(1,I)))+ZADVRE(I)*FK1
          ENDIF

C===========================================================================
       ENDIF
  410  CONTINUE
         CALL DATE_AND_TIME(VALUES=Dtime)
         WRITE(3,*) 'vreme pre presen,2', (Dtime(i),i=5,7)

c         CALL Vwrr(SILE,JEDN,'s-p')
 4567   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL VRESEN(A(LSK),SILE,MAXA,JEDN,2)
        IF (myid.ne.0) goto 2346
c         CALL Vwrr(SILE,JEDN,'s-k')
c     OVDE UCITATI RESENJA SA DISKA U VEKTOR SILE
c
  999 IF(NXNASTRAN.EQ.1) THEN
C
            NXNAS=99    
            OPEN (NXNAS,FILE='POTNXNAS',STATUS='OLD',FORM='FORMATTED',
     1            ACCESS='SEQUENTIAL')
            DO I=1,JEDN
               READ(NXNAS,*) IDUM,SILE(I)
               write(3,4000) IDUM,SILE(I)
 4000 format(i5,f15.5)               
               IF(IDUM.NE.I) STOP 'IDUM.NE.I  PAKV3D.FOR'
            ENDDO
            CLOSE(NXNAS,STATUS='KEEP')
      ENDIF

         CALL DATE_AND_TIME(VALUES=Dtime)
         WRITE(3,*) 'vreme posle presen,2', (Dtime(i),i=5,7)
C
           DO 440 I=1,JEDN
                TT1(I)=TT1(I)+SILE(I)
 440       CONTINUE
 
!       WRITE(*,*)'ITER= ',ITER
!       WRITE(*,*)'PERIOD= ',NNPER
!       WRITE(*,*)'STEP= ',KKORAK
      IF(NXNASTRAN.EQ.1) go to 491
       IF(INDFS.EQ.0) GO TO 491
C
      CALL VKONVTF(TT1,SILE,KONVV2,1,ID,ITER)
C
2346  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KONVV2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(MAXIT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        IF (KKORAK.GT.0) THEN
C
C         CALL VCKONV1(ICUR,TT1,CORD,ID,KONV1,UBRZ)
c         proba da nasilno ne radi iteracije
cz          nasilu=0
cz          if(nstac.eq.1.and.iter.eq.0) nasilu=0
cz          write(*,*) ' nasilu',nasilu
cz          if(nasilu.eq.0) then
          IF (KONVV2*KONV1.EQ.0) THEN
            ITER=ITER+1
            WRITE(IIZLAZ,*)'ITER=',ITER
            GO TO 100
          ENDIF
         ENDIF
        
        IF (myid.ne.0) goto 2345

        IF ((NJUTN.NE.0).AND.(INDOPT.EQ.0)) RETURN

 491     CALL DATE_AND_TIME(VALUES=Dtime)
         WRITE(3,*) 'vreme pre racunanja gradijenta', (Dtime(i),i=5,7)
C======================================================================= 
C RACUNANJE GRADIJENTA POTENCIJALA, BRZINA I ZAPREMINSKIH SILA
C=======================================================================
       DO JJ=1,NPT
          FZAPR(1,JJ)=0.D0
          FZAPR(2,JJ)=0.D0
          FZAPR(3,JJ)=0.D0
	   ENDDO	
C
C RACUNANJE PROTOKA DUZ KONTURE UNUTAR MREZE /Sneza april 2013. racuna se iza izracunavanja protoka po konturi/
C
C      IF (NKONT.GT.0) THEN
C      CALL VRKONT3(TT1,ID,KONT,NEL,AKONST,QUK,QUM,CORD,KKORAK)
C      ENDIF
C
C  PETLJA PO ELEMENTIMA
C
      DO 500 NBREL=1,NET
      IF(elemtip(NBREL).gt.100) THEN
        NTYPE=elemtip(NBREL)/100
        NTIMES=100
      ELSE
        NTYPE=elemtip(NBREL)/10
        NTIMES=10
      ENDIF
      NDIMEL=elemtip(NBREL)-NTIMES*NTYPE
      
      VZAPR=0.
      NMAT=NEL(NDIM+1,NBREL)
      NEMA=NEL(NDIM+2,NBREL)

      IF (NMAT.LE.0) GOTO 500
      IF (NEMA.GT.0) GOTO 500
 
      CALL VCLEARD(TT21,NDIM)
      CALL VPREB(TT21,TT1,NBREL)
      CALL VCLEARD(TT210,NDIM)
      CALL VPREB(TT210,TT10,NBREL)
C
      DO KLM=1,NDIMEL
       CK(1,KLM)=CORD(1,NEL(KLM,NBREL))
       CK(2,KLM)=CORD(2,NEL(KLM,NBREL))
       CK(3,KLM)=CORD(3,NEL(KLM,NBREL))
      ENDDO
      
      if(NTYPE.eq.2) CALL VLCK2D(CK,CKL,TTE,NDIMEL)
!  zapremiska sila po elementu
      SILEL(NBREL)=0.D0
! 
      DO KL=1,NDIMEL
       FZAP(1,KL)=0.D0
       FZAP(2,KL)=0.D0
       FZAP(3,KL)=0.D0
      ENDDO
C      write(iizlaz,*)"pakv3d"
C
C INDIKATOR NJUTN-KOTESOVE INTEGRACIJE (0-NE,1-DA)
      KOTES=0
      NGAUS=0
      IBRGT=2
      IBRKT=2
      IF(KOTES.EQ.1) IBRKT=3
      
        IBRKT1=IBRKT
        IBRKT2=IBRKT
        if(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) THEN
          IBRKT2=1
          IBRKT1=4
        ENDIF
        IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) IBRKT2=1
        IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) IBRKT1=1
        IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
          IBRKT2=1
          IBRKT1=3
        ENDIF
!         write(*,*) 'pre petlje po gausovim tackama, NBREL',NBREL
        DO  320 I=1,IBRKT1
        DO  320 J=1,IBRKT2
        DO  320 K=1,IBRKT2
        NGAUS=NGAUS+1
        PORNIEL(NBREL,NGAUS)=0.D0
        
	IF(KOTES.EQ.1) THEN
          IF(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) then
            CALL VGETRST(RT,ST,TT,WTET,NGAUS,IBRKT1)
            CALL VINTERPTETRA(RT,ST,TT,0,0,NDIMEL)
            WDT=WTET(1)*WTET(2)*WTET(3)*DETJ
          ELSE IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) THEN
            CALL VINTERP1D(RK(IBRKT1,I),0,0,NDIMEL,AJS)
            WDT=WK(IBRKT1,I)*thick(NBREL)*AJS
          ELSE IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) THEN
            CALL VINTERP2D48(R(IBRGT,J),S(IBRGT,K),0,0,NDIMEL,CKL,TTE)
            WDT=W(IBRGT,J)*W(IBRGT,K)*DETJ*thick(NBREL)
          ELSE IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
            RT=1./6
            ST=1./6
            WT=1./3
            IF(NGAUS.EQ.2) RT=2./3
            IF(NGAUS.EQ.3) ST=2./3
            CALL VINTERP2D36(RT,ST,0,0,NDIMEL,CKL,TTE)
            WWDT = WT*thick(NBREL)*DETJ/2
          ELSE
            CALL VINTER3(RK(IBRGT,I),SK(IBRGT,J),TK(IBRGT,K),0,0)
            WDT=WK(IBRKT,I)*WK(IBRKT,J)*WK(IBRKT,K)*DETJ
          ENDIF
	ELSE
          IF(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) then
            CALL VGETRST(RT,ST,TT,WTET,NGAUS,IBRKT1)
            CALL VINTERPTETRA(RT,ST,TT,0,0,NDIMEL)
            WDT=WTET(1)*WTET(2)*WTET(3)*DETJ
          ELSE IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) THEN
            CALL VINTERP1D(R(IBRKT1,I),0,0,NDIMEL,AJS)
            WDT=W(IBRKT1,I)*thick(NBREL)*AJS
          ELSE IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) THEN
!             CALL VINTERP2D48(R(IBRGT,I),S(IBRGT,J),0,0,NDIMEL,CKL,TTE)
!             WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ*thick(NBREL)
              CALL VINTERP2D48(R(IBRGT,J),S(IBRGT,K),0,0,NDIMEL,CKL,TTE)
              WDT=W(IBRGT,J)*W(IBRGT,K)*DETJ*thick(NBREL)
          ELSE IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
            RT=1./6
            ST=1./6
            WT=1./3
            IF(NGAUS.EQ.2) RT=2./3
            IF(NGAUS.EQ.3) ST=2./3
            CALL VINTERP2D36(RT,ST,0,0,NDIMEL,CKL,TTE)
            WDT = WT*thick(NBREL)*DETJ/2
          ELSE
            CALL VINTER3(R(IBRGT,I),S(IBRGT,J),T(IBRGT,K),0,0)
            WDT=W(IBRKT,I)*W(IBRKT,J)*W(IBRKT,K)*DETJ
          ENDIF
	ENDIF
        VZAPR=VZAPR+WDT
C
C OKVASENOST U GAUS TACKI
        PR=0.D0
        DO II=1,NDIMEL
           PR=PR+H(II)*(TT21(II)-CORD(IOSA,NEL(II,NBREL)))
        ENDDO
        OKVAS=0.
        IF(PR.GT.0.D0) OKVAS=1.
c porni pritisak u gaus tacki
        PORNI=PR*GAMA
        IF(PR.GT.0.D0) PORNIEL(NBREL,NGAUS)=PORNI
C      write(iizlaz,*)NBREL,NGAUS,PORNIEL(NBREL,NGAUS)
C=======================================================================
C RACUNANJE GRADIJENTA POTENCIJALA:
C
      GRADFX=0.D0
      GRADFY=0.D0
      GRADFZ=0.D0
C ISKLJUCENO DA NE RACUNA FILTRACIONE SILE I SILE OD PORNIH PRITISAKA      
      IF(PR.GT.0.D0) THEN
C      IF(PR.GT.0.D0.OR.ISIL.EQ.0) THEN
C RACUNANJE GRADIJENTA POTENCIJALA
          DO II=1,NDIMEL
            GRADFX=GRADFX+ZVHX(II)*TT21(II)
            GRADFY=GRADFY+ZVHY(II)*TT21(II)
            GRADFZ=GRADFZ+ZVHZ(II)*TT21(II)
          ENDDO
C
C  GRADIJENT PRITISKA (FILTRACIONE SILE PO JEDINICI ZAPREMINE U GAUS TACKI)
C
        DO KL=1,NDIMEL
         IF (IFIL.EQ.1)THEN
           GRADPX=GAMA*GRADFX
           GRADPY=GAMA*GRADFY
           GRADPZ=GAMA*GRADFZ
C           WRITE(*,*) 'NBREL',NBREL
C           WRITE(*,*) 'GAMA',GAMA
C           WRITE(*,*) 'GRADPX',GRADPX
CC           WRITE(*,*) 'GRADPY',GRADPY
C           WRITE(*,*) 'GRADPZ',GRADPZ
C           WRITE(*,*) 'WDT',WDT
C           WRITE(*,*) 'H(',KL,')',H(KL)
C           IF(DABS(AKONST(3,5,NMAT)).LT.1.E-10) THEN
C                FZAP(1,KL)=FZAP(1,KL)-H(KL)*GRADPX*WDT
C                FZAP(2,KL)=FZAP(2,KL)-H(KL)*GRADPY*WDT
C                FZAP(3,KL)=FZAP(3,KL)-H(KL)*GRADPZ*WDT
C           ELSE
            IF(IPOR.EQ.0) THEN
              FZAP(1,KL)=FZAP(1,KL)-H(KL)*GRADPX*WDT
              FZAP(2,KL)=FZAP(2,KL)-H(KL)*GRADPY*WDT
              FZAP(3,KL)=FZAP(3,KL)-H(KL)*GRADPZ*WDT
            ELSE
              FZAP(1,KL)=FZAP(1,KL)-H(KL)*GRADPX*WDT*AKONST(3,5,NMAT)
              FZAP(2,KL)=FZAP(2,KL)-H(KL)*GRADPY*WDT*AKONST(3,5,NMAT)
              FZAP(3,KL)=FZAP(3,KL)-H(KL)*GRADPZ*WDT*AKONST(3,5,NMAT)
            ENDIF
C           ENDIF
         ELSEIF(IFIL.EQ.2)THEN
            IF(IPOR.EQ.0) THEN
              FZAP(1,KL)=FZAP(1,KL)+ZVHX(KL)*PORNI*WDT
              FZAP(2,KL)=FZAP(2,KL)+ZVHY(KL)*PORNI*WDT
              FZAP(3,KL)=FZAP(3,KL)+ZVHZ(KL)*PORNI*WDT
            ELSE
        FZAP(1,KL)=FZAP(1,KL)+ZVHX(KL)*PORNI*WDT*(1.D0-AKONST(3,5,NMAT))
        FZAP(2,KL)=FZAP(2,KL)+ZVHY(KL)*PORNI*WDT*(1.D0-AKONST(3,5,NMAT))
        FZAP(3,KL)=FZAP(3,KL)+ZVHZ(KL)*PORNI*WDT*(1.D0-AKONST(3,5,NMAT))
            ENDIF
C PRITISAK OD SOPSTVENE TEZINE
C           GRADPX=GAMA*DELTA(IOSA,1)*OKVAS
C           GRADPY=GAMA*DELTA(IOSA,2)*OKVAS
C           GRADPZ=GAMA*DELTA(IOSA,3)*OKVAS
         ELSEIF(IFIL.EQ.0)THEN
           GRADPX=GAMA*GRADFX
           GRADPY=GAMA*GRADFY
           GRADPZ=GAMA*GRADFZ
C           IF(DABS(AKONST(3,5,NMAT)).LT.1.E-10) THEN
C              FZAP(1,KL)=FZAP(1,KL)-H(KL)*GRADPX*WDT
C              FZAP(2,KL)=FZAP(2,KL)-H(KL)*GRADPY*WDT
C              FZAP(3,KL)=FZAP(3,KL)-H(KL)*GRADPZ*WDT
C           ELSE
            IF(IPOR.EQ.0) THEN
              FZAP(1,KL)=FZAP(1,KL)-H(KL)*GRADPX*WDT
              FZAP(2,KL)=FZAP(2,KL)-H(KL)*GRADPY*WDT
              FZAP(3,KL)=FZAP(3,KL)-H(KL)*GRADPZ*WDT
              FZAP(1,KL)=FZAP(1,KL)+ZVHX(KL)*PORNI*WDT
              FZAP(2,KL)=FZAP(2,KL)+ZVHY(KL)*PORNI*WDT
              FZAP(3,KL)=FZAP(3,KL)+ZVHZ(KL)*PORNI*WDT
            ELSE
              FZAP(1,KL)=FZAP(1,KL)-H(KL)*GRADPX*WDT*AKONST(3,5,NMAT)
              FZAP(2,KL)=FZAP(2,KL)-H(KL)*GRADPY*WDT*AKONST(3,5,NMAT)
              FZAP(3,KL)=FZAP(3,KL)-H(KL)*GRADPZ*WDT*AKONST(3,5,NMAT)
        FZAP(1,KL)=FZAP(1,KL)+ZVHX(KL)*PORNI*WDT*(1.D0-AKONST(3,5,NMAT))
        FZAP(2,KL)=FZAP(2,KL)+ZVHY(KL)*PORNI*WDT*(1.D0-AKONST(3,5,NMAT))
        FZAP(3,KL)=FZAP(3,KL)+ZVHZ(KL)*PORNI*WDT*(1.D0-AKONST(3,5,NMAT))
            ENDIF
C           ENDIF
        ENDIF
C        
        ENDDO
C         
      ENDIF

      IF ( (PR.LT.0.D0) .AND. (INDFS.EQ.1) )THEN
         AKX=(1.D0-ZASIC)*AKONST(3,1,NMAT)
         AKY=(1.D0-ZASIC)*AKONST(3,2,NMAT)
         AKZ=(1.D0-ZASIC)*AKONST(3,3,NMAT)
      ELSE
         AKX=AKONST(3,1,NMAT)
         AKY=AKONST(3,2,NMAT)
         AKZ=AKONST(3,3,NMAT)
      ENDIF

C GRADIJENT U GAUS TACKI
       AJX1(NGAUS)=GRADFX
       AJY1(NGAUS)=GRADFY
       AJZ1(NGAUS)=GRADFZ
       GG(1,NBREL,NGAUS)=GRADFX
       GG(2,NBREL,NGAUS)=GRADFY
       GG(3,NBREL,NGAUS)=GRADFZ
C BRZINA U GAUS TACKI
       AJX(NGAUS)=-AKX*GRADFX
       AJY(NGAUS)=-AKY*GRADFY
       AJZ(NGAUS)=-AKZ*GRADFZ
       VG(1,NBREL,NGAUS)=-AKX*GRADFX
       VG(2,NBREL,NGAUS)=-AKY*GRADFY
       VG(3,NBREL,NGAUS)=-AKZ*GRADFZ
C
  320   CONTINUE
!       CASE DEFAULT
!         WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NBREL)
!       END SELECT
      
      IMCVOR=0
      IVAL=0
      SELECT CASE(elemtip(NBREL))
      CASE(12)
        IVAL=2
      CASE(13)
        IVAL=2
        IMCVOR=1
      CASE(23)
        IVAL=3
      CASE(26)
        IVAL=3
        IMCVOR=1
      CASE(24)
        IVAL=4
      CASE(28)
        IVAL=4
        IMCVOR=1
      CASE(34)
        IVAL=4
      CASE(38)
        IVAL=8
      CASE(310)
        IVAL=4
        IMCVOR=1
      CASE(320)
        IVAL=8
        IMCVOR=1
      CASE DEFAULT
!         WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NBREL)
      END SELECT
C       
 	DO KK=1,NDIMEL
	  NCVOR=NEL(KK,NBREL)
	  KK1=NRED(KK)
	  SELECT CASE(elemtip(NBREL))
          CASE(12:13)
            KK1=NRED1(KK)
          CASE(23,26)
            KK1=NRED26(KK)
          CASE(24,28)
            KK1=NRED28(KK)
          CASE(34,310)
            KK1=NREDT10(KK)
          CASE DEFAULT
!          WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NBREL)
          END SELECT
	  IF(KOTES.EQ.1) KK1=NREDK(KK)
	  IF(KK.LE.IVAL.OR.IMCVOR.EQ.0) THEN
C BRZINE U CVOROVIMA
          VECTJ(1,NCVOR)=VECTJ(1,NCVOR)+AJX(KK1)/IVECT(NCVOR)
          VECTJ(2,NCVOR)=VECTJ(2,NCVOR)+AJY(KK1)/IVECT(NCVOR)
          VECTJ(3,NCVOR)=VECTJ(3,NCVOR)+AJZ(KK1)/IVECT(NCVOR)
C GRADIJENT POTENCIJALA U CVOROVIMA
	  GRADJN(1,NCVOR)=GRADJN(1,NCVOR)+AJX1(KK1)/IVECT(NCVOR)
	  GRADJN(2,NCVOR)=GRADJN(2,NCVOR)+AJY1(KK1)/IVECT(NCVOR)
	  GRADJN(3,NCVOR)=GRADJN(3,NCVOR)+AJZ1(KK1)/IVECT(NCVOR)
C FILTRACIONE SILE U CVOROVIMA	
          FZAPR(1,NCVOR)=FZAPR(1,NCVOR)+FZAP(1,KK)
          FZAPR(2,NCVOR)=FZAPR(2,NCVOR)+FZAP(2,KK)
          FZAPR(3,NCVOR)=FZAPR(3,NCVOR)+FZAP(3,KK)
c ovaj uslov nije jasan          
c          FZAPR(1,NCVOR)=FZAPR(1,NCVOR)+FZAP(1,KK)*8./VZAPR
c          FZAPR(2,NCVOR)=FZAPR(2,NCVOR)+FZAP(2,KK)*8./VZAPR
c          FZAPR(3,NCVOR)=FZAPR(3,NCVOR)+FZAP(3,KK)*8./VZAPR
          ELSE IF(KK.GT.IVAL.AND.IMCVOR.EQ.1) THEN
          KK2=KK1/10
          KK3=KK1-KK2*10
C BRZINE U CVOROVIMA
          VECTJ(1,NCVOR)=VECTJ(1,NCVOR)+(AJX(KK2)/IVECT(NCVOR))/2
     1                                 +(AJX(KK3)/IVECT(NCVOR))/2
          VECTJ(2,NCVOR)=VECTJ(2,NCVOR)+(AJY(KK2)/IVECT(NCVOR))/2
     1                                 +(AJY(KK3)/IVECT(NCVOR))/2
          VECTJ(3,NCVOR)=VECTJ(3,NCVOR)+(AJZ(KK2)/IVECT(NCVOR))/2
     1                                 +(AJZ(KK3)/IVECT(NCVOR))/2
C GRADIJENT POTENCIJALA U CVOROVIMA
          GRADJN(1,NCVOR)=GRADJN(1,NCVOR)+(AJX1(KK2)/IVECT(NCVOR))/2
     1                                   +(AJX1(KK3)/IVECT(NCVOR))/2
          GRADJN(2,NCVOR)=GRADJN(2,NCVOR)+(AJY1(KK2)/IVECT(NCVOR))/2
     1                                   +(AJY1(KK3)/IVECT(NCVOR))/2
          GRADJN(3,NCVOR)=GRADJN(3,NCVOR)+(AJZ1(KK2)/IVECT(NCVOR))/2
     1                                   +(AJZ1(KK3)/IVECT(NCVOR))/2
C FILTRACIONE SILE U CVOROVIMA  
          FZAPR(1,NCVOR)=FZAPR(1,NCVOR)+FZAP(1,KK)
          FZAPR(2,NCVOR)=FZAPR(2,NCVOR)+FZAP(2,KK)
          FZAPR(3,NCVOR)=FZAPR(3,NCVOR)+FZAP(3,KK)
        ENDIF
      ENDDO
      FXEL=0.0
      FYEL=0.0
      FZEL=0.0
      DO KK=1,NDIMEL
       FXEL=FXEL+FZAP(1,KK)
       FYEL=FYEL+FZAP(2,KK)
       FZEL=FZEL+FZAP(3,KK)
      ENDDO
      SILEL(NBREL)=dsqrt(FXEL*FXEL+FYEL*FYEL+FZEL*FZEL)/VZAPR
      IF(DABS(SILEL(NBREL)).GT.1.E-15)
     1 WRITE(INEU,9002) NBREL,SILEL(NBREL)
 9002  FORMAT(I10,',',1PE15.7,',')
C======================================================================= 
 500  CONTINUE
C
       WRITE(INEU,9001)
 9001 FORMAT('-1,0.,')
C======================================================================= 
 480   CALL VPRKONT1(TT1,ID,UBRZ,CORD,NZAD,IOSA,IVECT,ZADVRE,KOJK)
C========================================================================
C RACUNANJE PROTOKA DUZ KONTURE UNUTAR MREZE
C
!          write(*,*) 'pre RKON3V'
!        WRITE(*,*) 'pre racunanja protoka kroz kont'
      IF (NKONT.GT.0) THEN
        CALL VRKON3V(VECTJ,KONT,VG,KKORAK,ISNUMER)
      ENDIF
C
C RACUNANJE PROTOKA DUZ KONTURE UNUTAR MREZE
C
!       IF (NKONT.GT.0) THEN
! C      CALL VRKONT3(TT1,ID,KONT,NEL,AKONST,QUK,QUM,CORD,KKORAK)
!         CALL VRKONT3V(VG,KONT,QUK,QUM,CORD,KKORAK)
!       ENDIF
C
C
      CALL DATE_AND_TIME(VALUES=Dtime)
      WRITE(*,*) 'vreme pre stampanja', (Dtime(i),i=5,7)
      WRITE(3,*) 'vreme pre stampanja', (Dtime(i),i=5,7)
      WRITE(IIZLAZ,*)'PERIOD NUMBER,KKORAK ',NNPER, KKORAK         
      WRITE(*,*)'PERIOD NUMBER,KKORAK ',NNPER, KKORAK         
      IF (KKORAK.GT.0) THEN
! C Sneza 24.06.2011 dodato da stampa rezultate samo za odredjene cvorove
! C stampanje potencijala
!         CALL VSTAU09c(TT1,CORD,ID,NPT,47,1,1,KKORAK,NZAD,NUMZAD,KONT,TT1,
!      1               NEL)
! C        CALL VPRKONT1c(TT1,ID,UBRZ,CORD,NZAD,IOSA,IVECT,ZADVRE,KOJK)
! C stampanje dubine vode
! C       CALL VSTAU09c(TT1,CORD,ID,NPT,49,11,1,KKORAK,NZAD,NUMZAD,KONT,TT1)
! C stampa rezultata u neu za sve cvorove
!        CALL VSTAS09(TT1,CORD,ID,NPT,49,10,1,KKORAK,NZAD,NUMZAD,KONT)
!        CALL VSTAU09(TT1,CORD,ID,NPT,49,9,1,KKORAK,NZAD,NUMZAD,KONT,TT1)
!        CALL VSTAU09(TT1,CORD,ID,NPT,49,1,1,KKORAK,NZAD,NUMZAD,KONT,TT1)
!        CALL VSTAU09(TT1,CORD,ID,NPT,49,11,1,KKORAK,NZAD,NUMZAD,KONT,TT1)
!       CALL VSTAU09(VECTJ,CORD,ID,NPT,49,41,3,KKORAK,NZAD,NUMZAD,KONT,TT1)
!        CALL VSTAU09(GRADJN,CORD,ID,NPT,49,51,3,KKORAK,NZAD,NUMZAD,KONT,
!      1            TT1)
!       CALL VSTAU09(FZAPR,CORD,ID,NPT,49,61,3,KKORAK,NZAD,NUMZAD,KONT,TT1)
! C stampanje pornih pritisaka 
! !!! PREPRAVITI STAU35 da radi sa svim tipovima elemenata i onda odkomentarisati
!        CALL VSTAU35(PORNIEL,KKORAK)
! c stampanje zapreminskih sila
! C stampanje rezultata u svim cvorovima
!       CALL VIZLL3D(ID,TT1,KKORAK,VECTJ,QUK,QUM,VVREME,FZAPR,VG,GG,
!      1            KOTES,IDJSTAMP1,NP3D1)
!         IF (INDSC.EQ.0) THEN
!       CALL VSTAG3D(TT1,ID,NASLOV,VVREME,KKORAK,1,NPT,18,0,NET,NEL,
!      1            KORAK,VECTJ,GRADJN,NZAD,FZAPR)
!         ENDIF
!       ENDIF
! C
! 
! 
C stampanje potencijala
      IF (INDSC.EQ.0.OR.INDSC.EQ.3.OR.INDSC.EQ.4) THEN
!        WRITE(*,*)'pre stau09c '         
       CALL VSTAU09c(TT1,NPT,47,1,1,KKORAK,NZAD,NUMZAD,KONT,TT1,
     1               NEL)
C        CALL VPRKONT1c(TT1,ID,UBRZ,CORD,NZAD,IOSA,IVECT,ZADVRE,KOJK)
C stampanje dubine vode
C       CALL VSTAU09c(TT1,CORD,ID,NPT,49,11,1,KKORAK,NZAD,NUMZAD,KONT,TT1)
C stampa rezultata u neu za sve cvorove
!         WRITE(*,*)'pre STAS09 10 '         
        CALL VSTAS09(TT1,NPT,49,10,1,KKORAK,NZAD,NUMZAD,KONT,ISNUMER)
!         WRITE(*,*)'pre STAU09 9 '         
       CALL VSTAU09(TT1,NPT,49,9,1,KKORAK,NZAD,NUMZAD,KONT,TT1,ISNUMER)
!         WRITE(*,*)'pre STAU09 19 '         
       CALL VSTAU09(TT1,NPT,49,19,1,KKORAK,NZAD,NUMZAD,KONT,TT1,ISNUMER)
!         WRITE(*,*)'pre STAU09 1 '         
       CALL VSTAU09(TT1,NPT,49,1,1,KKORAK,NZAD,NUMZAD,KONT,TT1,ISNUMER)
!         WRITE(*,*)'pre STAU09 11 '         
       CALL VSTAU09(TT1,NPT,49,11,1,KKORAK,NZAD,NUMZAD,KONT,TT1,ISNUMER)
!         WRITE(*,*)'pre STAU09 41 '         
       CALL VSTAU09(VECTJ,NPT,49,41,3,KKORAK,NZAD,NUMZAD,KONT,TT1,
     1             ISNUMER)
!         WRITE(*,*)'pre STAU09 51 '         
       CALL VSTAU09(GRADJN,NPT,49,51,3,KKORAK,NZAD,NUMZAD,KONT,
     1            TT1,ISNUMER)
!         WRITE(*,*)'pre STAU09 61 '         
!       CALL VSTAU09(FZAPR,NPT,49,61,3,KKORAK,NZAD,NUMZAD,KONT,TT1)
!  stampanje sile po jedinici zapremine
!           CALL VSTAS09T(NEL,49,79,1,KKORAK,NGPSIL,MAXTQE,ISNUMER)
C stampanje pornih pritisaka 
!!! PREPRAVITI STAU35 da radi sa svim tipovima elemenata i onda odkomentarisati
!         WRITE(*,*)'pre STAU35 '         
        CALL VSTAU35(PORNIEL,KKORAK,ISNUMER)
c stampanje zapreminskih sila
C stampanje rezultata u svim cvorovima
!        WRITE(*,*)'pre IZLL3D '         
!       CALL VIZLL3D(TT1,KKORAK,VECTJ,QUK,QUM,VVREME,FZAPR,VG,GG,
      ENDIF
      CALL VIZLL3D(TT1,KKORAK,VECTJ,VVREME,FZAPR,VG,GG,
     1            KOTES,ISNUMER)
        IF (INDSC.EQ.0) THEN
!         WRITE(*,*)'pre STAG3D '         
      CALL VSTAG3D(TT1,NASLOV,VVREME,KKORAK,1,NPT,18,0,NET,
     1            KORAK,VECTJ,GRADJN,NZAD,FZAPR,ISNUMER)
        ENDIF
!       IF (INDSC.EQ.1.OR.INDSC.EQ.3) THEN
!          CALL VSTAU09c(TT1,NPT,47,1,1,KKORAK,NZAD,NUMZAD,KONT,TT1,
!      1               ISNUMER)
!       ENDIF
! C STAMPANJE REZULTATA ZA pijezometre     
      IF(INDSC.EQ.2.OR.INDSC.EQ.4) THEN
!         WRITE(*,*)'pre STAU09CT '         
          CALL VSTAU09CT(TT1,47,VVREME,KKORAK,ISNUMER)
! !        WRITE(*,*)'pre STAU09CDPV '   
! !  uradi slobodnu numeraciju ???
!          CALL VSTAU09CDPV(TT1,VECTJ,50,VVREME,KKORAK,ISNUMER)
! !          CALL VSTAU09CDP(TT1,CORD,50,VVREME,KKORAK,NEL)
! !          CALL VSTAU09CDV(VECTJ,CORD,50,VVREME,KKORAK,NEL)
! !         CALL VSTAU09CTPP(PORNIEL,CORD,47,VVREME,KKORAK,NEL)
        ENDIF
      ENDIF
! Stampanje rezultata u presecima
!  prvo se ucitavaju preseci pa stampaju rezultati za svaki presek
      IF(ISPRESEK.GT.0) THEN
         write(*,*) "stampanje rezultata preseka"
         CALL VSTAMPPRESEKNEU(TT1,VECTJ,
     +       VREME,VVREME,KKORAK,ISNUMER)
!         CLOSE(250)
!        ENDDO
      ELSEIF(ISPRESEK.LE.-1) THEN
         CALL VSTAMPPRESEKNEU(TT1,VECTJ,
     +       VREME,VVREME,KKORAK,ISNUMER)
      ENDIF
      
! 
      IF(NXNASTRAN.EQ.1) RETURN
C
      KKORAK=KKORAK+1
C uslova da sledeci korak u tom periodu ne postoji i da ide na sledeci period
 2345 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      VRNPKOR=DABS(VREME(NNPER,KORAK+1))
      CALL MPI_BCAST(VRNPKOR,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      IF (VRNPKOR.LT.1.D-10) THEN
        GOTO 600
      ELSE
cz ide na sledeci korak      
        GOTO 35
      ENDIF
  600 CONTINUE
      if(myid.eq.0) then
        IF(NSTAC.EQ.1.AND.KKORAK.EQ.2) 
     1    CALL VIZBACI(TT1,CORD,ID,NPT,NEL,NET,IOSA,NDIM,IIZLAZ)
        CALL DATE_AND_TIME(VALUES=Dtime)
        WRITE(*,*) 'vreme posle stampanja', (Dtime(i),i=5,7)
        WRITE(3,*) 'vreme posle stampanja', (Dtime(i),i=5,7)
      endif 
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      RETURN      
      END
C=======================================================================
C=======================================================================
!       SUBROUTINE VIZLL3D(TT1,KORAK,VECTJ,QUK,QUM,VVREME,FZAPR,VG,GG,
      SUBROUTINE VIZLL3D(TT1,KORAK,VECTJ,VVREME,FZAPR,VG,GG,
     1                  KOTES,ISNUMBER)
      USE NODES
      USE ELEMENTS
      USE KONTURE
      USE ppr
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
      COMMON /ICITANJE/INPT
      COMMON /PRIKAZ/ INDSC,IZPOT
      COMMON /DJERDAP/ IDJERDAP,ISPRESEK

      DIMENSION D4(4),TT1(*),VECTJ(3,*),N4(4)
      DIMENSION FZAPR(3,*),NRED(8),NREDK(8),NRED1(3)
!       DIMENSION VG(3,NET,*),GG(3,NET,*),QUK(1000),QUM(1000),NREDT(4)
      DIMENSION VG(3,NET,*),GG(3,NET,*),NREDT(4)
      DIMENSION NREDT10(10),VGT(3),NRED26(6),NRED28(8)
      DATA NRED/8,4,2,6,7,3,1,5/  
      DATA NREDK/27,9,3,21,25,7,1,19/
      DATA NRED1/2,1,3/
      DATA NRED26/1,2,3,12,23,31/
      DATA NRED28/1,2,3,4,12,23,34,41/
      DATA NREDT/1,2,3,4/
      DATA NREDT10/1,2,3,4,12,23,13,14,24,34/
      DIMENSION NP3D1(18)
C
      if (indsc.eq.0) then
!       write(*,*) "stampanje potencijala"
C      
CS  POTENCIJAL U CVOROVIMA 
CE  NODE POTENTIAL 
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2001)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6001)
      WRITE (IIZLAZ,*)'STEP NUMBER  ',KORAK

      KG = 0
   20 DO 50 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 100
      IF(ISNUMBER.EQ.0) THEN
        N4(I) = KG 
      ELSE
        N4(I) = NCVEL(KG)
      ENDIF
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
      IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5012) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      ENDIF
      GO TO 77
C
   88 IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5013) (N4(J1),D4(J1),J1=1,I)
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
CS  BRZINA U PRAVCU X1     
CE  CURRENT DESNITIES IN DIRECTION X1
C

      KG = 0
  620 DO 650 I=1,4
       KG = KG +1
       IF(KG.GT.NPT) GO TO 1600
       IF(ISNUMBER.EQ.0) THEN
         N4(I) = KG 
       ELSE
         N4(I) = NCVEL(KG)
       ENDIF
       D4(I) = VECTJ(1,KG)
  650 CONTINUE
C
 1600 I = I - 1
      IF(IFORM.EQ.1) GO TO 688
      IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5012) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      ENDIF
      GO TO 677
C
  688 IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5013) (N4(J1),D4(J1),J1=1,I)
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
CS  Brzina U PRAVCU X2     
CE  CURRENT DESNITIES IN DIRECTION X2
C

      KG = 0
  720 DO 750 I=1,4
       KG = KG +1
       IF(KG.GT.NPT) GO TO 1700
       IF(ISNUMBER.EQ.0) THEN
         N4(I) = KG 
       ELSE
         N4(I) = NCVEL(KG)
       ENDIF
       D4(I) = VECTJ(2,KG)
  750 CONTINUE
C
 1700 I = I - 1
      IF(IFORM.EQ.1) GO TO 788
      IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5012) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      ENDIF
      GO TO 777
C
  788 IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5013) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
      ENDIF
  777 IF(KG.LT.NPT) THEN
      GO TO 720
      ENDIF
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2049)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6049)
C
CS  BRZINA U PRAVCU X3     
CE  CURRENT DESNITIES IN DIRECTION X3
C

      KG = 0
  820 DO 850 I=1,4
       KG = KG +1
       IF(KG.GT.NPT) GO TO 1800
       IF(ISNUMBER.EQ.0) THEN
         N4(I) = KG 
       ELSE
         N4(I) = NCVEL(KG)
       ENDIF
       D4(I) = VECTJ(3,KG)
  850 CONTINUE
C
 1800 I = I - 1
      IF(IFORM.EQ.1) GO TO 888
      IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5012) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      ENDIF
      GO TO 877
C
  888 IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5013) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
      ENDIF
  877 IF(KG.LT.NPT) THEN
      GO TO 820
      ENDIF
!       write(*,*) "stampanje brzina"

C
CS  BRZINE
C
      NGAUS=IBRGT*IBRGT*IBRGT
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2059)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6059)
      IF (INDSC.EQ.0)THEN
       DO I=1,NET
         IF(elemtip(I).gt.100) THEN
          NTYPE=elemtip(I)/100
          NTIMES=100
        ELSE
          NTYPE=elemtip(I)/10
          NTIMES=10
        ENDIF
        NDIMEL=elemtip(I)-NTIMES*NTYPE
        IF(ISNUMBER.EQ.0) THEN
          WRITE(IIZLAZ,5000) I
        ELSE
          WRITE(IIZLAZ,5000) MCVEL(I)
        ENDIF
       DO J=1,NDIMEL
         IF(ISNUMBER.EQ.0) THEN
           NC=NEL(J,I)
         ELSE
           NC=NCVEL(NEL(J,I))
         ENDIF
         NG=NRED(J)
         SELECT CASE(elemtip(I))
          CASE(12:13)
            NG=NRED1(J)
          CASE(23:26)
            NG=NRED26(J)
          CASE(34,310)
            NG=NREDT10(J)
          CASE DEFAULT
!          WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NBREL)
         END SELECT
!          IF(NDIMEL.eq.2) NG=NRED1(J)
!          IF(NDIMEL.eq.10) NG=NREDT10(J)
	 IF(KOTES.EQ.1) NG=NREDK(J)
	     IF(INPT.EQ.1) THEN
              IF(NG.GT.10) THEN
                NG1=NG/10
                NG2=NG-NG1*10
                DO K=1,3
                  VGT(K)=(VG(K,I,NG1)+(VG(K,I,NG2)))/2
                ENDDO
                WRITE(IIZLAZ,5017) NC,(VGT(K),K=1,3)
              ELSE
                WRITE(IIZLAZ,5017) NC,(VG(K,I,NG),K=1,3)
              ENDIF
	     ELSE
              IF(NG.GT.10) THEN
                NG1=NG/10
                NG2=NG-NG1*10
                DO K=1,3
                  VGT(K)=(VG(K,I,NG1)+(VG(K,I,NG2)))/2
                ENDDO
                WRITE(IIZLAZ,5007) NC,(VGT(K),K=1,3)
              ELSE
                WRITE(IIZLAZ,5007) NC,(VG(K,I,NG),K=1,3)
              ENDIF
	     ENDIF
       ENDDO
       ENDDO
      ENDIF
C
CS  GRADIJENTI
!       write(*,*) "stampanje gradijenata"
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2069)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6069)
      IF (INDSC.EQ.0)THEN
      DO I=1,NET
         IF(elemtip(I).gt.100) THEN
          NTYPE=elemtip(I)/100
          NTIMES=100
        ELSE
          NTYPE=elemtip(I)/10
          NTIMES=10
        ENDIF
        NDIMEL=elemtip(I)-NTIMES*NTYPE
        IF(ISNUMBER.EQ.0) THEN
          WRITE(IIZLAZ,5000) I
        ELSE
          WRITE(IIZLAZ,5000) MCVEL(I)
        ENDIF
        DO J=1,NDIMEL
         IF(ISNUMBER.EQ.0) THEN
           NC=NEL(J,I)
         ELSE
           NC=NCVEL(NEL(J,I))
         ENDIF
         NG=NRED(J)
         SELECT CASE(elemtip(I))
          CASE(12:13)
            NG=NRED1(J)
          CASE(23:26)
            NG=NRED26(J)
          CASE(34,310)
            NG=NREDT10(J)
          CASE DEFAULT
!          WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NBREL)
         END SELECT
         IF(NDIMEL.eq.2) NG=NRED1(J)
         IF(NDIMEL.eq.10) NG=NREDT10(J)
	   IF(KOTES.EQ.1) NG=NREDK(J)
	   IF(INPT.EQ.1) THEN
                IF(NG.GT.10) THEN
                  NG1=NG/10
                  NG2=NG-NG1*10
                  DO K=1,3
                    VGT(K)=(GG(K,I,NG1)+(GG(K,I,NG2)))/2
                  ENDDO
                  WRITE(IIZLAZ,5017) NC,(VGT(K),K=1,3)
                 ELSE
                  WRITE(IIZLAZ,5017) NC,(GG(K,I,NG),K=1,3)
	         ENDIF
	   ELSE
	     IF(NG.GT.10) THEN
                NG1=NG/10
                NG2=NG-NG1*10
                DO K=1,3
                  VGT(K)=(GG(K,I,NG1)+(GG(K,I,NG2)))/2
                ENDDO
                WRITE(IIZLAZ,5007) NC,(VGT(K),K=1,3)
             ELSE
                WRITE(IIZLAZ,5007) NC,(GG(K,I,NG),K=1,3)
             ENDIF
	   ENDIF
        ENDDO
      ENDDO
      ENDIF
c
      endif      
!       write(*,*) "stampanje zap sila"
C
CS  ZAPREMINSKE SILE
CE  BODY FORCES
C
C	WRITE(25,*)'C  N,    FX,         FY,      FZ'
!         write(*,*) 'IZLST, IDJERDAP, NP3D1',IDJERDAP,NP3D1(1)
      IF (GAMA.GT.0.D0) THEN
!         IF(IDJERDAP.EQ.1) WRITE(25,5017) NP3D1(1)
! 	    DO NN=1,NPT
! 	     IF(IDJERDAP.EQ.1) THEN
! 	        IF(IDJSTAMP1(NN).EQ.1) THEN
! 	           IF(INPT.EQ.1) THEN
!                 WRITE(25,5017) NN,(FZAPR(KK,NN),KK=1,3)
! 	           ELSE
!                 WRITE(25,5007) NN,(FZAPR(KK,NN),KK=1,3)
! 	           ENDIF
! 	        ENDIF 
! 	     ELSE
! 	        IF(INPT.EQ.1) THEN
!              WRITE(25,5017) NN,(FZAPR(KK,NN),KK=1,3)
! 	        ELSE
!              WRITE(25,5007) NN,(FZAPR(KK,NN),KK=1,3)
! 	        ENDIF
! 	     ENDIF
! 	    ENDDO
            CALL VICLEAR(NP3D1,18)
            write(*,*)"izl, IDJERDAP", IDJERDAP
           IF(IDJERDAP.GE.1) THEN
C   podprogram za citanje cvorova u kojima se stampaju filtracione sile
             II=1
             CALL VDJERDAPREADT(NP3D1,ISNUMBER,II)
             IFILE=105
             DO IK=1,1
               IFILE=IFILE+IDJERDAP
!                do ii=1,2
!                	faktor=ii/2.0
                  faktor=1
                WRITE(IFILE,5017) NP3D1(IK)
                DO NN=1,NP3D1(IK)
                IF(ISNUMBER.EQ.0) THEN
                WRITE(IFILE,5017) IDJSTAMP1(IK,NN),
!      1          (FZAPR(KK,NN)*faktor,KK=1,3)
     1          (FZAPR(KK,IDJSTAMP1(IK,NN))*faktor,KK=1,3)
                ELSE
                WRITE(IFILE,5017) NCVEL(IDJSTAMP1(IK,NN)),
     1          (FZAPR(KK,IDJSTAMP1(IK,NN))*faktor,KK=1,3)
                ENDIF
               ENDDO
             ENDDO
           ELSEIF(IDJERDAP.EQ.-1)THEN
             DO II=1,18
               CALL VDJERDAPREADT(NP3D1,ISNUMBER,II)
             IIDJERDAP=II 
            write(*,*)"IIDJERDAP",IIDJERDAP,NP3D1(1)
             IFILE=105
             DO IK=1,1
               IFILE=IFILE+IIDJERDAP
                  faktor=1
                WRITE(IFILE,5017) NP3D1(IK)
               DO NN=1,NP3D1(IK)
                IF(ISNUMBER.EQ.0) THEN
                WRITE(IFILE,5017) IDJSTAMP1(IK,NN),
!      1          (FZAPR(KK,NN)*faktor,KK=1,3)
     1          (FZAPR(KK,IDJSTAMP1(IK,NN))*faktor,KK=1,3)
                ELSE
                WRITE(IFILE,5017) NCVEL(IDJSTAMP1(IK,NN)),
     1          (FZAPR(KK,IDJSTAMP1(IK,NN))*faktor,KK=1,3)
                ENDIF
	       ENDDO
             ENDDO
            ENDDO
!  stampa za sve cvorove u jedan fajl           
           ELSE
!             do ii=1,2
              faktor=1
             WRITE(25,5017) NPT
             DO NN=1,NPT
!                IF(INPT.EQ.1) THEN
               IF(ISNUMBER.EQ.0) THEN
                WRITE(25,5017) NN,(FZAPR(KK,NN)*faktor,KK=1,3)
               ELSE
            WRITE(25,5017) NCVEL(NN),(FZAPR(KK,NN)*faktor,KK=1,3)
               ENDIF
! 	        ELSE
!              WRITE(25,5007) NN,(FZAPR(KK,NN),KK=1,3)
! 	        ENDIF
             ENDDO
!             enddo
           ENDIF
!            DO NN=1,NPT
! 	     IF(IDJERDAP.EQ.1) THEN
! 	        IF(IDJSTAMP1(NN).EQ.5) THEN
! 	           IF(INPT.EQ.1) THEN
!                 WRITE(25,5017) NN,(FZAPR(KK,NN),KK=1,3)
! 	           ELSE
!                 WRITE(25,5007) NN,(FZAPR(KK,NN),KK=1,3)
! 	           ENDIF
! 	        ENDIF 
! 	     ELSE
! 	        IF(INPT.EQ.1) THEN
!              WRITE(25,5017) NN,(FZAPR(KK,NN),KK=1,3)
! 	        ELSE
!              WRITE(25,5007) NN,(FZAPR(KK,NN),KK=1,3)
! 	        ENDIF
! 	     ENDIF
! 	  ENDDO
      ENDIF
C
CS  PROTOK KROZ KONTURU     
CE  CURRENT FLUX WITHIN CONTOUR
C
!       write(*,*) "stampanje protoka kroz kont"
!       write(*,*) 'INDSC,NKONT', INDSC,NKONT
      IF(INDSC.EQ.0) THEN
       IF(NKONT.GT.0) THEN
         QK=0.
         QM=0.
         DO L=1,NKONT
!        write(*,*)"l,QUK(L),QUM(L)", L,QUK(L),QUM(L)
          QK=QK+QUK(L)
          QM=QM+QUM(L) 
        ENDDO
        WRITE(IIZLAZ,7023) VVREME
!         write(*,*)"NKONT", NKONT
        DO L=1,NKONT
         IF(ISNUMBER.EQ.0) THEN
         WRITE(IIZLAZ,5007) L,QUK(L),QUM(L),QUK(L)+QUM(L)
         ELSEIF(ISNUMBER.EQ.1) THEN
         LI=LINSN(L)
!          write(*,*)"L,LI", L,LI
         WRITE(IIZLAZ,5007) LI,QUK(L),QUM(L),QUK(L)+QUM(L)
         ENDIF
        ENDDO
        WRITE(IIZLAZ,5008) QK,QM,QK+QM
        WRITE(21,7023) VVREME
        DO L=1,NKONT
         IF(ISNUMBER.EQ.0) THEN
         WRITE(21,5007) L,QUK(L),QUM(L),QUK(L)+QUM(L)
         ELSEIF(ISNUMBER.EQ.1) THEN
         LI=LINSN(L)
         WRITE(21,5007) LI,QUK(L),QUM(L),QUK(L)+QUM(L)
         ENDIF
        ENDDO 
        WRITE(21,5008) QK,QM,QK+QM
       ENDIF
      ENDIF
C
      RETURN
 5000 FORMAT(/' ELEMENT',I10)     
 5002 FORMAT(4(I6,1PD13.5))
 5003 FORMAT(4(I6,3X,F10.3))
 5005 FORMAT(F10.3,10(1PE13.5))
 5007 FORMAT(I6,3(1PE13.5))
 5008 FORMAT('UKUPNO',3(1PE13.5))
 5012 FORMAT(4(I10,1PD13.5))
 5013 FORMAT(4(I10,3X,F10.3))
 5017 FORMAT(I10,3(1PE13.5))
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
 2049 FORMAT(/'     C V O R N E   B R Z I N E   U   X 3  P R A V C U '//
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
 6049 FORMAT(/'    N O D A L    V E L O C I T I E S   I N   D I R E C T'
     1,' I O N   X 3 '//
     1' NODE  VELOCITY    NODE  VELOCITY    NODE  VELOCITY    NODE VELOC
     1ITY   '/'  No.',3(15X,'No.'))
 6059 FORMAT(//'    N O D A L    V E L O C I T I E S'/
     1' NODE  VELOCITY VX  VELOCITY VY  VELOCITY VZ')
 6069 FORMAT(//'    N O D A L    G R A D I E N T S  '/
     1' NODE  GRADIENT GX  GRADIENT GY  GRADIENT GZ')
 7023 FORMAT(//'  F L U X   W I T H I N  C O N T O U R, TIME =',1PD13.5/
     1' KONTURA     ULAZNI       IZLAZNI      UKUPNI PROTOK-1')
      END
C=======================================================================
C
C RACUNANJE PROTOKA DUZ KONTURE ZA 3D ELEMENTE
C
C=========================================================================
      SUBROUTINE VRKONT3(TT1,ID,KONT,AKONST,CORD,KOR)
      USE ELEMENTS
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DETJS,DETJ,NBREL
      COMMON /STAC/ NSTAC,NSLP,IBROJ
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /DODAT/ NDIMM
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
C
      DIMENSION R(3,3),S(3,3),T(3,3),W(3,3)
      DIMENSION TT21(44)
C
      DIMENSION KONT(9,MAXLIN,*),ID(1,*)
      DIMENSION AFIFI(21,21),AM1(21,21),AK(21,21),TT1(*),CORD(3,*)
      DIMENSION AKONST(3,5,*),QPR(21)
      DIMENSION A12(3),A14(3),DN(3)
C
C*    BLOK ZA STAMPANJE PROTOKA PO ELEMENTIMA ZA GRAFIKU - PRIVREMENO
      JEDAN=1
      NULA=0
      ZERO=0.
      MJ=-1
      M8=8
C
      IN=15
      IGRAF=48
      WRITE(IGRAF,1000) KOR,IN,JEDAN
      WRITE(IGRAF,3003) 
      WRITE(IGRAF,1010) ZERO,ZERO,ZERO
      WRITE(IGRAF,1000)NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(IGRAF,1000)NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
C
C 0=ANY,1=DISPL,2=ACCEL,3=FORCE,4=STRESS,5=STRAIN,6=TEMP,OTHER=USER
C     NODAL(7) ,ELEMENTAL(8)
      WRITE(IGRAF,1000) NULA,NULA,M8,M8
C
C =1, CAN NOT LINEARY COMBINE THIS OUTPUT
C =1, COMP(0-2) ARE THE X,Y,Z COMPONENT
C =1, THIS VECTOR HAS CONTROIDAL OR NODAL OUTPUT
      WRITE(IGRAF,1000) JEDAN,NULA,JEDAN
 1000 FORMAT(10I5)
 1010 FORMAT(3(1PD12.4))
 3003 FORMAT(' PROTOK PO JEDINICI POVRSINE')
C*
      R(3,1)=-0.7745966692415
      R(3,2)=0.0
      R(3,3)=0.77459666924148
      S(3,1)=-0.7745966692415
      S(3,2)=0.0
      S(3,3)=0.77459666924148
      T(3,1)=-0.7745966692415
      T(3,2)=0.0
      T(3,3)=0.77459666924148
      W(3,1)=0.55555555555556
      W(3,2)=0.88888888888889
      W(3,3)=0.55555555555556


      R(2,1)=-0.57735026918963
      R(2,2)=0.5773502691896
      S(2,1)=-0.57735026918963
      S(2,2)=0.5773502691896
      T(2,1)=-0.57735026918963
      T(2,2)=0.5773502691896
      W(2,1)=1.000000000000000
      W(2,2)=1.000000000000000

      R(1,1)=0.0
      S(1,1)=0.0
      T(1,1)=0.0
      W(1,1)=2.0
C
! C ODEDIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
!       IBRGT=3
!       IF (NDIM.EQ.8) IBRGT=2
C
      NLM=0
      DO 723 KK=1,NKONT
        QUK(KK)=0.D0
        QUM(KK)=0.D0
C
	DO 134 II=1,LIN(KK)
	NLM=NLM+1
C	
      NBREL=KONT(1,II,KK)
      NMAT=NEL(NDIM+1,NBREL)
      NEMA=NEL(NDIM+2,NBREL)
      
      IF(elemtip(NBREL).gt.100) THEN
        NTYPE=elemtip(NBREL)/100
        NTIMES=100
      ELSE
        NTYPE=elemtip(NBREL)/10
        NTIMES=10
      ENDIF
      NDIMEL=elemtip(NBREL)-NTIMES*NTYPE
C ODEDIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
      IBRGT=3
      IF (NDIMEL.LE.8) IBRGT=2
      IBRGT1=IBRGT
      IF(elemtip(NBREL).eq.310.or.elemtip(NBREL).eq.34) then
        IBRGT=4
        IBRGT1=1
      endif

      IF (NMAT.LE.0) GOTO 134
      IF (NEMA.GT.0) GOTO 134
C
      CALL VCLEARD(TT21,NDIM)
      CALL VPREB(TT21,TT1,NBREL)
C
      DO KLM=1,NDIM
       CK(1,KLM)=CORD(1,NEL(KLM,NBREL))
       CK(2,KLM)=CORD(2,NEL(KLM,NBREL))
       CK(3,KLM)=CORD(3,NEL(KLM,NBREL))
      ENDDO
C
        ICV1=KONT(2,II,KK)
        ICV2=KONT(3,II,KK)
        ICV3=KONT(4,II,KK)
        ICV4=KONT(5,II,KK)
        ICV5=KONT(6,II,KK)
        ICV6=KONT(7,II,KK)
        ICV7=KONT(8,II,KK)
        ICV8=KONT(9,II,KK)
C
      DO 163 K=1,NDIMEL
      DO 162 N=1,NDIMEL
         AM1(K,N)=0.D0    
         AFIFI(K,N)=0.D0    
  162 CONTINUE
  163 CONTINUE
C
      DO 180 I=1,IBRGT
      DO 170 J=1,IBRGT
      DO 160 L=1,IBRGT
C     
      CALL VINTER3(R(IBRGT,I),S(IBRGT,J),T(IBRGT,L),0,0)
C  
      WDT=W(IBRGT,I)*W(IBRGT,J)*W(IBRGT,L)*DETJ
C
C OKVASENOST U GAUS TACKI
        PR=0.D0
        DO IK=1,NDIMEL
           PR=PR+H(IK)*(TT21(IK)-CORD(IOSA,NEL(IK,NBREL)))
        ENDDO
      IF ( PR.LT.0.D0 .AND. INDFS.EQ.1 )THEN
         AKX=(1.D0-ZASIC)*AKONST(3,1,NMAT)
         AKY=(1.D0-ZASIC)*AKONST(3,2,NMAT)
         AKZ=(1.D0-ZASIC)*AKONST(3,3,NMAT)
         AKS=(1.D0-ZASIC)*AKONST(3,4,NMAT)
      ELSE
         AKX=AKONST(3,1,NMAT)
         AKY=AKONST(3,2,NMAT)
         AKZ=AKONST(3,3,NMAT)
         AKS=AKONST(3,4,NMAT)
      ENDIF
C
      DO 165 K=1,NDIMEL
      DO 164 N=1,NDIMEL
C
         AM1(K,N)=AM1(K,N)+AKS*H(K)*H(N)*WDT
         AFIFI(K,N)=AFIFI(K,N)+(AKX*ZVHX(K)*ZVHX(N)+
     1                      AKY*ZVHY(K)*ZVHY(N)+AKZ*ZVHZ(K)*ZVHZ(N))*WDT
C
  164 CONTINUE
  165 CONTINUE
C
  160 CONTINUE
  170 CONTINUE
  180 CONTINUE
C
      DO I=1,NDIMEL
      DO J=1,NDIMEL
        AK(I,J)=AFIFI(I,J)
      ENDDO
      ENDDO
C
      DO IK=1,NDIMEL
        QPR(IK)=0.D0
        DO J=1,NDIMEL
          if((TT21(IK)-CORD(IOSA,NEL(IK,NBREL))).gt.0)then
              QPR(IK)=QPR(IK)+AK(IK,J)*TT21(J)
          endif
        ENDDO
      ENDDO
C
      QPRI=0.
      QPRU=0.
      DO J=1,NDIMEL
!         IF(NEL(J,NBREL).EQ.ICV1.OR.NEL(J,NBREL).EQ.ICV2.
!      &     OR.NEL(J,NBREL).EQ.ICV3.OR.NEL(J,NBREL).EQ.ICV4) THEN
           IF(QPR(J).LT.0.) THEN
              QPRI=QPRI+QPR(J)
           ELSE
              QPRU=QPRU+QPR(J)
           ENDIF
!         ENDIF
      ENDDO
      
      QUM(KK)=QUM(KK)+QPRI
      QUK(KK)=QUK(KK)+QPRU
      
C*    STAMPANJE PROTOKA PO ELEMENTU
      IF(ICV4.EQ.0) THEN
        DO I=1,3
          A12(I)=CORD(I,ICV2)-CORD(I,ICV1)
          A14(I)=CORD(I,ICV3)-CORD(I,ICV1)
        ENDDO
      ELSE
        DO I=1,3
          A12(I)=CORD(I,ICV2)-CORD(I,ICV1)
          A14(I)=CORD(I,ICV4)-CORD(I,ICV1)
        ENDDO
      ENDIF
C     DODATAK ZA PROTOK PO JEDINICI POVRSINE      
      CALL VAXBV(A12,A14,DN)
      DUM=DSQRT(DN(1)*DN(1)+DN(2)*DN(2)+DN(3)*DN(3))
      QPRUK= -QPRI/DUM
C      QPRUK= (QPRI+QPRU)/DUM
C      
      WRITE(IGRAF,5000) NLM,QPRUK
C*
  134 CONTINUE
C
  723 CONTINUE
C
C*    KRAJ BLOKA ZA STAMPANJE PROTOKA PO ELEMENTU
      WRITE(IGRAF,5000) MJ,ZERO
 5000 FORMAT(I5,3(1PD12.4))
C*
      RETURN
      END
C=======================================================================
      SUBROUTINE VREDOSLED(NRED,NREDEL,NDIMEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NRED(*),NREDEL(*)
      
      CALL VICLEAR(NRED,20)
      DO I=1,NDIMEL
        NRED(I)=NREDEL(I)
      ENDDO
      
      RETURN
      END
C=========================================================================
C
C RACUNANJE PROTOKA DUZ KONTURE UNUTAR MREZE
C
C=========================================================================
      SUBROUTINE VRKON3V(VECTJ,KONT,VG,KOR,ISNUMER)
      USE NODES
      USE ELEMENTS
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /STAC/ NSTAC,NSLP,IBROJ
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /DODAT/ NDIMM
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PRIKAZ/ INDSC,IZPOT
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
C
      DIMENSION R42(3),R13(3),R12(3),EN(3),EN1(3),V(3)
C
      DIMENSION KONT(9,MAXLIN,*)
      DIMENSION VG(3,NET,*)
      DIMENSION VECTJ(3,*)
!       DIMENSION QUK(1000),QUM(1000),NRED(10),NRED26(6),NRED28(8)
      DIMENSION NRED(10),NRED26(6),NRED28(8)
      DIMENSION NRED8(8),NRED1(3),NREDT4(4),NREDT10(10)
      DATA NRED8/8,4,2,6,7,3,1,5/
      DATA NRED1/2,1,12/
      DATA NRED26/1,2,3,12,23,31/
      DATA NRED28/1,2,3,4,12,23,34,41/
      DATA NREDT4/1,2,3,4/
      DATA NREDT10/1,2,3,4,12,23,13,14,24,34/
C
C*    BLOK ZA STAMPANJE PROTOKA PO ELEMENTIMA ZA GRAFIKU - PRIVREMENO
      JEDAN=1
      NULA=0
      ZERO=0.
      MJ=-1
      M8=8
C
      IN=15
      IGRAF=48
      WRITE(IGRAF,1000) KOR,IN,JEDAN
      WRITE(IGRAF,3003) 
      WRITE(IGRAF,1010) ZERO,ZERO,ZERO
      WRITE(IGRAF,1000)NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(IGRAF,1000)NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
C ako se stampa vrednost po elementu mora da postoji sledeci red, od verzije Femap 10.0
      WRITE(IGRAF,1000)NULA
C 0=ANY,1=DISPL,2=ACCEL,3=FORCE,4=STRESS,5=STRAIN,6=TEMP,OTHER=USER
C     NODAL(7) ,ELEMENTAL(8)
      WRITE(IGRAF,1000) NULA,NULA,M8,M8
C
C =1, CAN NOT LINEARY COMBINE THIS OUTPUT
C =1, COMP(0-2) ARE THE X,Y,Z COMPONENT
C =1, THIS VECTOR HAS CONTROIDAL OR NODAL OUTPUT
      WRITE(IGRAF,1000) JEDAN,NULA,JEDAN
 1000 FORMAT(10I5)
 1010 FORMAT(3(1PD12.4))
 3003 FORMAT(' PROTOK PO JEDINICI POVRSINE')
C*
C
      ONE=1.
      PI=4.*DATAN(ONE)
      NLM=0
C
      DO 723 KK=1,NKONT
        QUK(KK)=0.D0
        QUM(KK)=0.D0
!         write(*,*)'KK,LIN(KK)',KK,LIN(KK)
C
	DO 134 II=1,LIN(KK)
C
          NBREL=KONT(1,II,KK)
!            write(*,*)'KK,NBREL',KK,NBREL
          NMAT=NEL(NDIM+1,NBREL)
          NEMA=NEL(NDIM+2,NBREL)
          NLM=NLM+1
          
          IF(elemtip(NBREL).gt.100) THEN
            NTYPE=elemtip(NBREL)/100
            NTIMES=100
          ELSE
            NTYPE=elemtip(NBREL)/10
            NTIMES=10
          ENDIF
          NDIMEL=elemtip(NBREL)-NTIMES*NTYPE
          
          SELECT CASE(elemtip(NBREL))   
          CASE(12)
            NSCV=2
            CALL VREDOSLED(NRED,NRED1,NDIMEL)
          CASE(13)
            NSCV=3
            CALL VREDOSLED(NRED,NRED1,NDIMEL)
          CASE(23)
            NSCV=2
            CALL VREDOSLED(NRED,NRED26,NDIMEL)
          CASE(26)
            NSCV=3
            CALL VREDOSLED(NRED,NRED26,NDIMEL)
          CASE(24)
            NSCV=2
            CALL VREDOSLED(NRED,NRED28,NDIMEL)
          CASE(28)
            NSCV=3
            CALL VREDOSLED(NRED,NRED28,NDIMEL)
          CASE(34)
            NSCV=3
            CALL VREDOSLED(NRED,NREDT10,NDIMEL)
          CASE(310)
            NSCV=6
            CALL VREDOSLED(NRED,NREDT10,NDIMEL)
          CASE(38)
            NSCV=4
            CALL VREDOSLED(NRED,NRED8,NDIMEL)
          CASE DEFAULT
            WRITE(*,*) 'GRESKA, NEMA PRIKAZA ZA DATI ELEMENT'
          END SELECT
          
          IF (NMAT.LE.0) GOTO 134
          IF (NEMA.GT.0) GOTO 134
C
	  ICV1=KONT(2,II,KK)
	  ICV2=KONT(3,II,KK)
	  ICV3=KONT(4,II,KK)
	  ICV4=KONT(5,II,KK)   
	  ICV5=KONT(6,II,KK)
          ICV6=KONT(7,II,KK)
          ICV7=KONT(8,II,KK)
          ICV8=KONT(9,II,KK)
          IC1=0
          IC2=0
          IC3=0
          IC31=1
          IC32=1
          IC33=1
          IC34=1
          IC4=1
          IC41=0
          IC5=0
          IC51=1
          IC52=1
          IC6=0
          IC61=1
          IC62=1
          IC7=0
          IC71=1
          IC72=1
          IC8=0
!  za 1d elemente ucitava elemente i cvor za koji se racuna protok          
!            write(*,*) "NBREL",NBREL
!              write(*,*) "ICV1",ICV1
!             write(*,*) "ICV2",ICV2
!             write(*,*) "ICV3",ICV3
!           write(*,*) "ICV4",ICV4
!           write(*,*) "ICV5",ICV5
!           write(*,*) "ICV6",ICV6
!           write(*,*) "ICV7",ICV7
!           write(*,*) "ICV8",ICV8IC5
          DO KLM=1,NDIMEL
!             write(*,*) "NEL(",KLM,",",NBREL,")",NEL(KLM,NBREL)
            IF(ICV1.EQ.NEL(KLM,NBREL)) then
              IC1=KLM
              IF(NTYPE.EQ.1) go to 99 
            ENDIF
            IF(ICV2.EQ.NEL(KLM,NBREL)) IC2=KLM
            IF(ICV3.EQ.NEL(KLM,NBREL)) THEN
!               WRITE(*,*) 'elemtip(',NBREL,')',elemtip(NBREL)
!               IF(elemtip(NBREL).eq.13) THEN
!                 IC31=NRED(KLM)/10
!                 IC32=NRED(KLM)-IC31*10
!                 IC33=1
!                 IC34=2
!                 IC3=1
!              ELSEIF(elemtip(NBREL).eq.26.and.ICV6.eq.0) THEN
             IF(elemtip(NBREL).eq.26.and.ICV6.eq.0) THEN
                IC31=NRED(KLM)/10
                IC32=NRED(KLM)-IC31*10
                IC33=1
                IC34=2
                IC3=1
!             write(*,*) 'el',II,NBREL
             ELSE
                IC31=KLM
                IC32=1
                IC33=0
                IC34=1
                IC3=1
              ENDIF
            ENDIF
            IF(ICV4.EQ.NEL(KLM,NBREL)) THEN
              IC41=1
              IC4=KLM
            ENDIF
            IF(ICV5.EQ.NEL(KLM,NBREL)) THEN
              IC51=NRED(KLM)/10
              IC52=NRED(KLM)-IC51*10
              IC5=1
            ENDIF
            IF(ICV6.EQ.NEL(KLM,NBREL)) THEN
              IC61=NRED(KLM)/10
              IC62=NRED(KLM)-IC61*10
              IC6=1
            ENDIF
            IF(ICV7.EQ.NEL(KLM,NBREL)) THEN
              IC71=NRED(KLM)/10
              IC72=NRED(KLM)-IC71*10
              IC7=1
            ENDIF
            IF(ICV8.EQ.NEL(KLM,NBREL)) IC8=KLM
          ENDDO
 99       IF(NTYPE.EQ.1) THEN
            ICV1=NEL(1,NBREL)
            ICV2=NEL(2,NBREL)
            ICV3=NEL(3,NBREL)
          ELSEIF(NTYPE.EQ.2) THEN
            IF(IC1.EQ.1.AND.IC2.EQ.2) ICV3=NEL(3,NBREL)
            IF(IC1.EQ.2.AND.IC2.EQ.3) ICV3=NEL(1,NBREL)
            IF(IC1.EQ.3.AND.IC2.EQ.1) ICV3=NEL(2,NBREL)
          ENDIF
!              write(*,*) "ICV1,ICV2,ICV3",ICV1,ICV2,ICV3

 100      IF(NTYPE.LE.2.AND.IC5.EQ.0) THEN
            R12(1)=CORD(1,ICV1)-CORD(1,ICV2)
            R12(2)=CORD(2,ICV1)-CORD(2,ICV2)
            R12(3)=CORD(3,ICV1)-CORD(3,ICV2)
            RI=DSQRT(R12(1)*R12(1)+R12(2)*R12(2)+R12(3)*R12(3))
            EN(1)=R12(1)/RI
            EN(2)=R12(2)/RI
            EN(3)=R12(3)/RI
            IF(NTYPE.EQ.1) go to 200
!             GOTO 131
          ENDIF
C         
          IF(NSCV.eq.4) THEN
            DO I=1,3
              R13(I)=CORD(I,ICV3)-CORD(I,ICV1)
              R42(I)=CORD(I,ICV2)-CORD(I,ICV4)
            ENDDO
          ELSE
            DO I=1,3
              R13(I)=CORD(I,ICV3)-CORD(I,ICV1)
              R42(I)=CORD(I,ICV2)-CORD(I,ICV1)
            ENDDO
          ENDIF

          IF(NTYPE.EQ.2) THEN
           CALL VAXBV(R42,R13,EN1)  
           CALL VAXBV(R42,EN1,EN)  
          ELSE
           CALL VAXBV(R42,R13,EN)  
          ENDIF
          AI=DSQRT(EN(1)**2+EN(2)**2+EN(3)**2)
          EN(1)=EN(1)/AI
          EN(2)=EN(2)/AI
          EN(3)=EN(3)/AI
          AI=AI/2.

          
!   131      DO i=1,NDIMEL
!             WRITE(*,*) 'NRED(',I,')',NRED(I)
!            ENDDO
!            do i = 1, 3
!              WRITE(*,*) 'VG(1,NBREL,',I,')',VG(1,NBREL,I)
!            enddo
!             write(*,*) "IC1",IC1
!            write(*,*) "IC2",IC2
!            write(*,*) "IC31",IC31
!            write(*,*) "IC32",IC32
!            write(*,*) "IC33",IC33
!            write(*,*) "IC34",IC34
!            write(*,*) "IC4",IC4
!            write(*,*) "IC41",IC41
!            write(*,*) "IC5",IC5
!            write(*,*) "IC51",IC51
!            write(*,*) "IC52",IC52
!            write(*,*) "IC6",IC6
!            write(*,*) "IC61",IC61
!            write(*,*) "IC62",IC62
!            write(*,*) "IC7",IC7
!            write(*,*) "IC71",IC71
!            write(*,*) "IC72",IC72
!    131     write(*,*) "pre racunanja brzina"
!             V(1)=0.0
!             V(2)=0.0
!             V(3)=0.0
!  131      IF(NTYPE.EQ.3) THEN
 131           DO I=1,3
               V(I)=(VG(I,NBREL,NRED(IC1))+VG(I,NBREL,NRED(IC2))+
     1          ((VG(I,NBREL,IC31)+VG(I,NBREL,IC32)*IC33)/IC34)*IC3+
     2          VG(I,NBREL,NRED(IC4))*IC41+
     3          (VG(I,NBREL,IC51)/2+VG(I,NBREL,IC52)/2)*IC5+
     4          (VG(I,NBREL,IC61)/2+VG(I,NBREL,IC62)/2)*IC6+
     5          (VG(I,NBREL,IC71)/2+VG(I,NBREL,IC72)/2)*IC7)/NSCV
            ENDDO
!           ELSEIF (NTYPE.EQ.2) THEN
!            DO I=1,3
!               V(I)=(VG(I,NBREL,NRED(IC1))+VG(I,NBREL,NRED(IC2)))/2
! !                V(I)=(VG(I,NBREL,NRED(IC1))+VG(I,NBREL,NRED(IC2))+
! !      1         ((VG(I,NBREL,IC31)+VG(I,NBREL,IC32)*IC33)/IC34)*IC3
! !      5          )/NSCV
!            ENDDO
!           ENDIF
!           write(*,*) "posle racunanja brzina"
C      V(1)=(VECTJ(1,ICV1)+VECTJ(1,ICV2)+VECTJ(1,ICV3)+VECTJ(1,ICV4))/4.
C      V(2)=(VECTJ(2,ICV1)+VECTJ(2,ICV2)+VECTJ(2,ICV3)+VECTJ(2,ICV4))/4.
C      V(3)=(VECTJ(3,ICV1)+VECTJ(3,ICV2)+VECTJ(3,ICV3)+VECTJ(3,ICV4))/4.
!       WRITE(*,*) 'NSCV',NSCV
!         DO I=1,3
!            WRITE(*,*) 'V(',I,')=',V(I)
!            WRITE(*,*) 'EN(',I,')=',EN(I)
!         ENDDO
!            write(3,*) "EN",EN(1),EN(2),EN(3)
!            write(3,*) "V",V(1),V(2),V(3)
          VN=EN(1)*V(1)+EN(2)*V(2)+EN(3)*V(3)
          VINT=V(1)**2+V(2)**2+V(3)**2
!            write(3,*) "VI",V(1),V(2),V(3)
!            write(3,*) "VINT",VINT
           VI=0.0
          IF(VINT.GT.1.D-14) THEN
            VI=DSQRT(VINT)
            AL=DACOS(VN/VI)*180./PI
          ENDIF
 200      IF(NTYPE.EQ.1) THEN
             VN=EN(1)*VG(1,NBREL,NRED(IC1))+EN(2)*
     1        VG(2,NBREL,NRED(IC1))+EN(3)*VG(3,NBREL,NRED(IC1))
             QN=-VN*thick(NBREL)
!                 WRITE(*,*) 'VI,VG',VN,VG(1,NBREL,NRED(IC1)),EN(1)
          ELSEIF(NTYPE.EQ.2) THEN
!             QN=-VI*RI*thick(NBREL)
            QN=-VN*RI*thick(NBREL)
!              write(3,*) "QN,VN,RI",QN,VN,RI
!             write(3,*) "QN,",QN
          ELSE
            QN=-VN*AI
	  ENDIF
	  IF(QN.LT.0.) THEN
            QUM(KK)=QUM(KK)+QN
          ELSE
            QUK(KK)=QUK(KK)+QN
          ENDIF
!            WRITE(3,*) 'NBREL,thick(NBREL),VN,',NBREL,thick(NBREL),VN
C*    STAMPANJE PROTOKA PO ELEMENTU
C      
      WRITE(IGRAF,5000) NLM,-VN
C*
  134   CONTINUE
!            WRITE(*,*) 'QUM,QUK,AI',QUM(1),QUK(1),AI
C
  723 CONTINUE
C*    KRAJ BLOKA ZA STAMPANJE PROTOKA PO ELEMENTU
      WRITE(IGRAF,5000) MJ,ZERO
 5000 FORMAT(I10,3(1PD12.4))
C*
C
      QK=0.
      QM=0.
      DO L=1,NKONT
         QK=QK+QUK(L)
         QM=QM+QUM(L)
      ENDDO
       WRITE(IIZLAZ,*)' KONTURA, ULAZNI,      IZLAZNI,    UKUPNI PROTOK'
!        WRITE(*,*) "NKONT",NKONT
      DO L=1,NKONT
       IF(ISNUMER.EQ.0)THEN
       WRITE(IIZLAZ,5055) L,QUK(L),QUM(L),(QUK(L)+QUM(L))
       ELSEIF(ISNUMER.EQ.1)THEN
       LI=LINSN(L)
!        WRITE(*,*) "L,LI",L,LI
       WRITE(IIZLAZ,5055) LI,QUK(L),QUM(L),(QUK(L)+QUM(L))
       ENDIF
      ENDDO
       WRITE(IIZLAZ,*)' SUMARNI  ULAZNI,      IZLAZNI,    UKUPNI PROTOK'
       WRITE(IIZLAZ,5056) QK,QM,QK+QM
       WRITE(21,*)' KONTURA, ULAZNI,      IZLAZNI,    UKUPNI PROTOK'
      DO L=1,NKONT
       IF(ISNUMER.EQ.0)THEN
       WRITE(21,5055) L,QUK(L),QUM(L),(QUK(L)+QUM(L))
       ELSEIF(ISNUMER.EQ.1)THEN
       LI=LINSN(L)
       WRITE(21,5055) LI,QUK(L),QUM(L),(QUK(L)+QUM(L))
       ENDIF
      ENDDO
       WRITE(21,*)' SUMARNI  ULAZNI,      IZLAZNI,    UKUPNI PROTOK'
       WRITE(21,5056) QK,QM,QK+QM
C stampanje protoka po 1 konturi u POT.TXT
        If(INDSC.EQ.1)then
            WRITE(47,*) QUK(1)
            WRITE(47,*) -1
        endif
 5005 FORMAT(10X,10(1PE13.5))  	
 5055 FORMAT(I6,10(1PE13.5))  	
 5056 FORMAT(6X,10(1PE13.5))  	
C
      RETURN
      END
C=======================================================================
C
C RACUNANJE PROTOKA DUZ KONTURE ZA 3D ELEMENTE - preko brzina / Sneza april 2013./
C
C=========================================================================
!       SUBROUTINE VRKONT3V(VG,KONT,NEL,QUK,QUM,CORD,KOR)
      SUBROUTINE VRKONT3V(VG,KONT,CORD,KOR)
      USE ELEMENTS
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DETJS,DETJ,NBREL
      COMMON /STAC/ NSTAC,NSLP,IBROJ
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /DODAT/ NDIMM
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
C
      COMMON /PRIKAZ/ INDSC,IZPOT
C
      DIMENSION R42(3),R13(3),EN(3),V(3)
C
C      DIMENSION KONT(9,MAXLIN,*),NEL(NDIMM,*)
      DIMENSION KONT(9,MAXLIN,*)
      DIMENSION VG(3,NET,*),CORD(3,*)
      DIMENSION NRED(8)
      DATA NRED/8,4,2,6,7,3,1,5/  
C
C*    BLOK ZA STAMPANJE PROTOKA PO ELEMENTIMA ZA GRAFIKU - PRIVREMENO
      JEDAN=1
      NULA=0
      ZERO=0.
      MJ=-1
      M8=8
C
      IN=15
      IGRAF=48
      WRITE(IGRAF,1000) KOR,IN,JEDAN
      WRITE(IGRAF,3003) 
      WRITE(IGRAF,1010) ZERO,ZERO,ZERO
      WRITE(IGRAF,1000)NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(IGRAF,1000)NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
C
C 0=ANY,1=DISPL,2=ACCEL,3=FORCE,4=STRESS,5=STRAIN,6=TEMP,OTHER=USER
C     NODAL(7) ,ELEMENTAL(8)
      WRITE(IGRAF,1000) NULA,NULA,M8,M8
C
C =1, CAN NOT LINEARY COMBINE THIS OUTPUT
C =1, COMP(0-2) ARE THE X,Y,Z COMPONENT
C =1, THIS VECTOR HAS CONTROIDAL OR NODAL OUTPUT
      WRITE(IGRAF,1000) JEDAN,NULA,JEDAN
 1000 FORMAT(10I5)
 1010 FORMAT(3(1PD12.4))
 3003 FORMAT(' PROTOK PO JEDINICI POVRSINE')
C*
      ONE=1.
      PI=4.*DATAN(ONE)
C
      NLM=0
C
      DO 723 KK=1,NKONT
C
	DO 134 II=1,LIN(KK)
C
          NBREL=KONT(1,II,KK)
          NMAT=NEL(NDIM+1,NBREL)
          NEMA=NEL(NDIM+2,NBREL)
          NLM=NLM+1

          IF (NMAT.LE.0) GOTO 134
          IF (NEMA.GT.0) GOTO 134
CIGRAF
	  ICV1=KONT(2,II,KK)
	  ICV2=KONT(3,II,KK)
	  ICV3=KONT(4,II,KK)
	  ICV4=KONT(5,II,KK)

          DO KLM=1,NDIM
            IF(ICV1.EQ.NEL(KLM,NBREL)) IC1=KLM
            IF(ICV2.EQ.NEL(KLM,NBREL)) IC2=KLM
            IF(ICV3.EQ.NEL(KLM,NBREL)) IC3=KLM
            IF(ICV4.EQ.NEL(KLM,NBREL)) IC4=KLM
          ENDDO
C
          DO I=1,3
            R13(I)=CORD(I,ICV3)-CORD(I,ICV1)
            R42(I)=CORD(I,ICV2)-CORD(I,ICV4)
          ENDDO

          CALL VAXBV(R42,R13,EN)   
          AI=DSQRT(EN(1)*EN(1)+EN(2)*EN(2)+EN(3)*EN(3))
          EN(1)=EN(1)/AI
          EN(2)=EN(2)/AI
          EN(3)=EN(3)/AI
          AI=AI/2.
          V(1)=(VG(1,NBREL,NRED(IC1))+VG(1,NBREL,NRED(IC2))+
     1          VG(1,NBREL,NRED(IC3))+VG(1,NBREL,NRED(IC4)))/4. 
          V(2)=(VG(2,NBREL,NRED(IC1))+VG(2,NBREL,NRED(IC2))+
     1          VG(2,NBREL,NRED(IC3))+VG(2,NBREL,NRED(IC4)))/4. 
          V(3)=(VG(3,NBREL,NRED(IC1))+VG(3,NBREL,NRED(IC2))+
     1          VG(3,NBREL,NRED(IC3))+VG(3,NBREL,NRED(IC4)))/4. 
C      V(1)=(VECTJ(1,ICV1)+VECTJ(1,ICV2)+VECTJ(1,ICV3)+VECTJ(1,ICV4))/4.
C      V(2)=(VECTJ(2,ICV1)+VECTJ(2,ICV2)+VECTJ(2,ICV3)+VECTJ(2,ICV4))/4.
C      V(3)=(VECTJ(3,ICV1)+VECTJ(3,ICV2)+VECTJ(3,ICV3)+VECTJ(3,ICV4))/4.
C brzina u pravcu normale
          VN=EN(1)*V(1)+EN(2)*V(2)+EN(3)*V(3)
          VINT=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)
          IF(VINT.GT.1.D-15) THEN
C intezitet vektora brzine
            VI=DSQRT(VINT)
            AL=DACOS(VN/VI)*180./PI
          ENDIF
C
C*    STAMPANJE PROTOKA PO ELEMENTU
C      
      WRITE(IGRAF,5000) NLM,-VN
C*
  134 CONTINUE
C
  723 CONTINUE
C
C*    KRAJ BLOKA ZA STAMPANJE PROTOKA PO ELEMENTU
      WRITE(IGRAF,5000) MJ,ZERO
 5000 FORMAT(I5,3(1PD12.4))
C*
      RETURN
      END
C=======================================================================
      SUBROUTINE VSPARS_ZADLEV(NZADJ,NZADP,NZAD,CORD,VVREME,TABF,
     1                                              NTABFT,ITFMAX)
      USE STIFFNESS
      USE PREDISCRIBED
      IMPLICIT NONE
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO PRESCRIBED POTENTIAL. - UPDATING OF LEFT SIDE WITH BLOCKS
C .
C ......................................................................
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
      COMMON /CDEBUG/ IDEBUG

      DOUBLE PRECISION PENALTY_FACTOR

      INTEGER NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      INTEGER IULAZ,IIZLAZ,IPAKT,ISIL,IFIL,IPOR
      INTEGER IDEBUG

      INTEGER NZADJ,NZADP,NZAD,IOSA,ITFMAX,NTABFT
      INTEGER I,J,K, NJ, NDIJ, istat,INTERPAXIS
      DOUBLE PRECISION A,CORD,FK1,ZASIC,TABF,VVREME

      DIMENSION NZADJ(*),NZAD(3,*),CORD(3,*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*)
        
      PENALTY_FACTOR=1.0D35
      
      IF(NBLOCK.GT.1) THEN
       WRITE(IIZLAZ,*) "IMSL SPARSE WITH BLOCKS IS NOT SUPPORTED"
       STOP
      ENDIF
      
      if (.not.allocated(NZADC)) then 
        allocate(NZADC(NZADP),STAT=istat)
        CALL VICLEAR8(NZADC,NZADP)       
!         CALL VICLEAR(NZADC,NZADP)       
        DO J = 1, nonzeros
          IF(rows(J).EQ.columns(J)) THEN
            DO I=1,NZADP
              IF (NZAD(2,I).EQ.4.OR.NZAD(2,I).EQ.5) EXIT
              NJ=NZADJ(I)
              IF(NJ.eq.rows(J)) THEN
                NZADC(I)=J
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      endif

       DO 10 I=1,NZADP
         IF(IPAKT.NE.1) THEN
! ! ! ! ! ! ! ! ! ! 
!   PAKV
! ! ! ! ! ! ! ! 
          CALL VTIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,I)),NZAD(3,I))
          IF(NZAD(2,I).EQ.6.AND.CORD(IOSA,NZAD(1,I)).GT.FK1) goto 10
!  za Grancarevo zavisnost na granicama modela
!  desna obala
          IF(NZAD(2,I).EQ.7) Then
           FK1=301.0+0.332*FK1
           if(CORD(IOSA,NZAD(1,I)).GT.FK1) GOTO 10
          ENDIF
          IF(NZAD(2,I).EQ.-1) Then
! C suborutine PREDF_INTERP pakv4.for           
!           CALL VPREDF_INTERP(I,INTERPAXIS,NUMAXISPTSX,INTAXISPOINTX,
!      1                  ITFMAX,NZAD,XAXISPTCORD,CORD,TABF,VVREME,FK1,1)
!         CASE (2)
!            INTERPAXIS=2
!           CALL VPREDF_INTERP(I,INTERPAXIS,NUMAXISPTSY,INTAXISPOINTY,
!      1                  ITFMAX,NZAD,YAXISPTCORD,CORD,TABF,VVREME,FK1,1)
!         CASE (3)
!            INTERPAXIS=3
!           CALL VPREDF_INTERP(I,INTERPAXIS,NUMAXISPTSZ,INTAXISPOINTZ,
!      1                  ITFMAX,NZAD,ZAXISPTCORD,CORD,TABF,VVREME,FK1,1)
!         CASE DEFAULT
!         END SELECT
           CALL VPREDF_INTERP(I,NZAD,CORD,FK1,1)
           if(CORD(IOSA,NZAD(1,I)).GT.FK1) GOTO 10
          ENDIF
!  leva obala - linearno Grancarevo
          IF(NZAD(2,I).EQ.0) Then
           CALL VPREDF_INTERP(I,NZAD,CORD,FK1,2)
           IF(CORD(IOSA,NZAD(1,I)).GT.FK1) GOTO 10
          ENDIF
! ! ! ! ! ! ! ! ! ! 
!   PAKT
! ! ! ! ! ! ! ! 
!          ELSEIF (IPAKT.EQ.1) THEN
!           IF(NZAD(2,I).NE.0) THEN
!            CALL VTIMFUN(TABF,FK1,VVREME,ITFMAX(NZAD(3,I)),NZAD(3,I))
!           ELSE
! !  Djerdap interpolacija zadatih temperatura na spoju beton/stena    
!            IF(NZAD(3,I).EQ.1) THEN
!              INTERPAXIS=1
!              CALL VPRED_INTERP(I,INTERPAXIS,NUMAXISPTSX,INTAXISPOINTX,
!      1              ITFMAX,NZAD,XAXISPTCORD,CORD,TABF,VVREME,FK1)
!            ENDIF
!            IF(NZAD(3,I).EQ.2) THEN
!              INTERPAXIS=2
!              CALL VPRED_INTERP(I,INTERPAXIS,NUMAXISPTSY,INTAXISPOINTY,
!      1               ITFMAX,NZAD,YAXISPTCORD,CORD,TABF,VVREME,FK1)
!            ENDIF
!            IF(NZAD(3,I).EQ.3) THEN
!              INTERPAXIS=3
!              CALL VPRED_INTERP(I,INTERPAXIS,NUMAXISPTSZ,INTAXISPOINTZ,
!      1               ITFMAX,NZAD,ZAXISPTCORD,CORD,TABF,VVREME,FK1)
!            ENDIF
!           ENDIF
         ENDIF

         IF (NZAD(2,I).EQ.4.OR.NZAD(2,I).EQ.5) GOTO 10
         
         IF(NZADC(I).GT.0) stiff(NZADC(I)) = PENALTY_FACTOR
!         WRITE(3,*)"zadlev",I,NZAD(1,I),NZADC(I),stiff(NZADC(I))
   10  CONTINUE
   
      RETURN
      END

C$DEBUG      
C=========================================================================
      SUBROUTINE VRACN3DT(TT1,SILE,
     1NZAD,ZADVRE,NGPSIL,MAXA,SKEF,SKEFN,
     2FZAPR,VREME,TABF,TT10,UBRZ,UBRZ0,AK,
     3VECTJ,IVECT,POVSIL,GRADJN,
     4ITFMAX,AKONST,NASLOV,ICUR,VG,GG,INDPT,ISNUMER)
      USE mcm_database
      USE PPR
      USE STIFFNESS
      USE NODES
      USE ELEMENTS
      USE PREDISCRIBED
      USE KONTURE
      USE pflux
      USE Crtanje

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'paka.inc'
      INCLUDE 'mpif.h'
c      COMMON A(17000)
c      REAL A
      COMMON /IME/ IME
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DETJS,DETJ,NBREL
      COMMON /KRITER/ IDOKL,MAXCUR
      COMMON /TREN13/ FPS1,FPS2,FPS3,FPS4
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /LSK1/ LSK
      COMMON /JEDANP/ INDJED,NBRF,NGL,INDTLO
      COMMON /FLUX/ KISA
      COMMON /TACNOS/ EPSTR,MAXIT,NJRAP
      COMMON /ZADPO/ LNZADJ
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /SCRATC/ ISCRC
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /NDESUK/ NDES,IDPRIT,IFORM
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT,GAMA
      COMMON /VREPER/ NPER,NTABFT
      COMMON /PROMEN/ NJUTN,INDOPT,INDPRO
      COMMON /DODAT/ NDIMM
      COMMON /DYNAM/ NDIN
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /STAC/ NSTAC,NSLP,IBROJ
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /VDP/ DT,NKORP,NN,NZAV
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
      COMMON /SREDIN/ IMASS      
C
       COMMON /SILZP3/ FZAP(3,21)
C
      COMMON /PRIKAZ/ INDSC,IZPOT
C
      COMMON /PIJEZO/ CPOR(3,1000),NPIJEZ(20,100),NODP(100),NPIJ,NPOR,
     1                NPORCV(1000),NPOREL(1000),NEPOR,NPORE,NEPIJ
      COMMON /PORCEL/ RMIN(34),SMIN(34),TMIN(34)
      COMMON /NXNAST/ NXNASTRAN
      COMMON /SOLVER/ ISOLVER
      COMMON /NELDIM/ NDIMEL
      COMMON /DJERDAP/ IDJERDAP,ISPRESEK
      COMMON /ICITANJE/ INPT
      COMMON /STAMPAZT/ NPRINT
      COMMON /INDBRANA/ IBRANA
C
      CHARACTER*80 NASLOV
      DIMENSION TT1(*),SILE(*),ZADVRE(*)
      DIMENSION VG(3,NET,*),GG(3,NET,*)
      DIMENSION TT10(*),UBRZ(*),UBRZ0(*)
      DIMENSION NZAD(3,*),NGPSIL(12,*),MAXA(*)
!       DIMENSION ID(1,*),NZAD(3,*),NGPSIL(12,*),MAXA(*),CORD(3,*)
      DIMENSION SKEF(NDES,*),AKONST(3,5,*),SKEFN(*)
      DIMENSION POVSIL(4,*),PR1(21),ITFMAX(*)
C
!       DIMENSION KONT(9,MAXLIN,*),KOJK(*)

      DIMENSION TT21(44),TT210(44)
      DIMENSION R(3,3),S(3,3),T(3,3),W(3,3)
      DIMENSION RK(3,3),SK(3,3),TK(3,3),WK(3,3)
      DIMENSION FZAPR(3,*),VECTJ(3,*),GRADJN(3,*)
      DIMENSION IVECT(*),LM2(44),ICUR(*)
      DIMENSION PERM(3,4)
      DIMENSION VREME(NPER,950)
      DIMENSION TABF(2,NTABFT,*)
      DIMENSION AM1(21,21),AK(NDES,*)

      DIMENSION RS2(21)
      DIMENSION RS3(21)
      DIMENSION F36(44)
C
      DIMENSION AFIFI1(21,21),AFIFI(21,21)
      DIMENSION PUS(21,21),AJX(20),AJY(20),AJZ(20)
      DIMENSION AJX1(20),AJY1(20),AJZ1(20),NRED(8),NREDK(8)
      DIMENSION NRED1(2),NREDT4(4),NREDT10(10),NRED26(6),NRED28(8)
      DIMENSION HT(7,3),CKL(3,8),TTE(2,3),WTET(3),INDEKS(10)
      DIMENSION A12(3),A13(3),DN(3)
!       promenjena duzina NRED sa 8 na 10 i dodata 2 clana
      DATA NRED/8,4,2,6,7,3,1,5/  
      DATA NREDK/27,9,3,21,25,7,1,19/
      DATA NRED1/2,1/
      DATA NRED26/1,2,3,12,23,31/
      DATA NRED28/1,2,3,4,12,23,34,41/
      DATA NREDT4/1,2,3,4/
      DATA NREDT10/1,2,3,4,12,23,13,14,24,34/
      integer ierr, myid
      integer*8 LM2,IPROM
!       ,NGPSIL

      CHARACTER*3 STATTT
      DIMENSION NDJEL(21),NP3D1(5)
      logical OLDNEWW
      integer Dtime(8)
      CHARACTER*50 IME
      CHARACTER*53 IZIPOCT

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      IF (myid.ne.0) goto 1235
  
      if (.not.allocated(PORNIEL)) allocate(PORNIEL(NET,10),STAT=istat) 
!       if(.not.allocated(IDJSTAMP1))allocate(IDJSTAMP1(1,NPT),STAT=istat) 
      IF(MAXSIL.GT.0) THEN
       if(.not.allocated(PFLUXEL))allocate(PFLUXEL(MAXSIL),STAT=istat) 
      ENDIF
      IF(MAXTQE.GT.0) THEN
       if(.not.allocated(FCONEL))allocate(FCONEL(MAXTQE),STAT=istat) 
       if(.not.allocated(TOKOLINEL))allocate(TOKOLINEL(MAXTQE),
     1       STAT=istat) 
      ENDIF
!       DO I=1,1
!         DO JJ=1,NPT
!           IDJSTAMP1(I,JJ)=0
!         ENDDO
!       ENDDO
C
!        CALL VICLEAR(NP3D1,5)
!       CALL VICLEAR(IDJSTAMP1,NPT)
C      IDJERDAP=1
!          write(*,*) 'IDJERDAP',IDJERDAP
!       IF(IDJERDAP.GT.0)THEN
! C   program za citanje cvorova u kojima se stampaju filtracione sile
!  za Dejerdap se radi po lamelama pa uvek stampa temperature za cvrstocu
!               CALL VDJERDAPREADT(IDJSTAMP1,NP3D1,ISNUMER)
! C
!        ENDIF
C INDIKATOR ZA OSU TEZINE - ne treba za pakt
!  9998   IF (IOSA.EQ.0) IOSA=3

C INDIKATOR ZA TRAZENJE IZOHIPSI PO FAJLOVIMA:
       INDJED=0

C ZBOG PAK-a NEQ-BROJ JEDNACINA .EQ. JEDN
       NEQ=JEDN
!        ICONT=1

C IDOKL JE BROJ CVOROVA KOJI SE UZIMAJU KAO AKTIVNI NA LINIJI CURENJA 
C ICUR JE NIZ U KOME SE NALAZE CVOROVI NA LINIJI CURENJA
C
C  INICIJALIZACIJA KONVERGENCIJE ZA LINIJU CURENJA
       KONV1=1  

         MAXCUR=0
!        IDOKL=0
!        CALL VINIDOK(NZAD,ICUR,CORD)

C UTVRDJIVANJE KOD KOLIKO ELEMENATA SE POJAVLJUJE JEDAN CVOR
      CALL VMVECT(IVECT)

      DO I=1,JEDN
       TT1(I)=0.D0
!        UBRZ(I)=0.D0
      ENDDO      
       
! C DEFINISANJE CVOROVA ZA PROTOKE   ne treba za pakt
!       DO KK=1,NKONT
! 	DO II=1,LIN(KK)
! 	  KOJK(KONT(2,II,KK))=KOJK(KONT(2,II,KK))+1
! 	  KOJK(KONT(3,II,KK))=KOJK(KONT(3,II,KK))+1
! 	  KOJK(KONT(4,II,KK))=KOJK(KONT(4,II,KK))+1
! 	  KOJK(KONT(5,II,KK))=KOJK(KONT(5,II,KK))+1
! 	  KOJK(KONT(6,II,KK))=KOJK(KONT(6,II,KK))+1
! 	  KOJK(KONT(7,II,KK))=KOJK(KONT(7,II,KK))+1
! 	  KOJK(KONT(8,II,KK))=KOJK(KONT(8,II,KK))+1
! 	  KOJK(KONT(9,II,KK))=KOJK(KONT(9,II,KK))+1
!         ENDDO
!       ENDDO	  
C
      R(3,1)=-0.7745966692415
      R(3,2)=0.0
      R(3,3)=0.77459666924148
      S(3,1)=-0.7745966692415
      S(3,2)=0.0
      S(3,3)=0.77459666924148
      T(3,1)=-0.7745966692415
      T(3,2)=0.0
      T(3,3)=0.77459666924148
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
      TK(3,1)=-1.0
      TK(3,2)=0.0
      TK(3,3)=1.0
      WK(3,1)=0.33333333333333
      WK(3,2)=1.33333333333333
      WK(3,3)=0.33333333333333
C 
      R(2,1)=-0.5773502691896
      R(2,2)=0.57735026918963
      S(2,1)=-0.5773502691896
      S(2,2)=0.57735026918963
      T(2,1)=-0.5773502691896
      T(2,2)=0.57735026918963
      W(2,1)=1.000000000000000
      W(2,2)=1.000000000000000

      R(1,1)=0.0
      S(1,1)=0.0
      T(1,1)=0.0
      W(1,1)=2.0
      
      PI=4.D0*DATAN(1.D0)
            
C===========================================================================
      IAXIS=6
C
      INULTO=0
      KKORAK=1
      VVREME=0.D0
C indikator LINTE=0 - linearan proracun
      LINTE=0
C=========================================================================
c petlja po periodima
1235  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NPER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      DO 600 NNPER=1,NPER
      
      KORAK=0
      IF(mcm_kojpak.eq.7)then
      KORAK=PAKSkorBR 
      ENDIF
c petlja po koracima
35    KORAK=KORAK+1
      
      
      
      if(myid.ne.0) goto 3456
C
C=========================================================================   
C  RAD SA STACIONARNIM STANJEM U PRVOM KORAKU I NESTACIONARNIM U OSTALIM
C  ( GEOLOSKI PROFILI ZA CERNE)
C
      IF (KKORAK.EQ.1) then
        NNSTAC=NSTAC
C       WRITE(19) KKORAK
      endif
cz uslov da prvi korak radi stacionarno
      IF (KKORAK.EQ.1.AND.NNSTAC.EQ.2) NSTAC=1
      IF (KKORAK.GT.1.AND.NNSTAC.EQ.2) NSTAC=0
C
cz      IF (KKORAK.GT.0) THEN
      TIME=VREME(NNPER,KORAK)
      VVREME=VVREME+TIME
cz      ELSE
cz       TIME=1.D7
cz       KORAK=KORAK-1
cz      ENDIF

      DO I=1,NPT 
       DO J=1,3
!         VECTJ(J,I)=0.D0
        GRADJN(J,I)=0.D0
       ENDDO
      ENDDO

C     ZADAVANJE POCETNE TEMPERATURE
C      write(*,*)'INDPT',INDPT
C      write(3,*)'INDPT',INDPT
      IF (KKORAK.EQ.1.AND.NNSTAC.EQ.0) THEN
      write(*,*)'INDPT',INDPT
      write(3,*)'INDPT',INDPT
          IF(INDPT.EQ.0) THEN
C          IF(INDPT.EQ.O) THEN
             CALL VCLEAR(TT10,JEDN)
             DO I=1,JEDN
                TT10(I)=POCT
             ENDDO
          ELSE
!         UCITAVANJE TEMPERATURA KOJE CE BITI POCETNE
            IB=INDEX(IME,'.')
            IF (IB.EQ.0) THEN
             IZIPOCT=trim(IME) // '.UL'
            ELSE
             IZIPOCT=IME(1:IB-1)//'.UL'
            ENDIF
               OPEN(99,FILE=IZIPOCT,STATUS='OLD',FORM='UNFORMATTED',
     1              ACCESS='SEQUENTIAL')
             CALL VPOCETNET(TT10,NPT,99,1)
             CLOSE (99)		  
          ENDIF
      ELSE
         DO I=1,JEDN
            TT10(I)=TT1(I)
         ENDDO
      ENDIF
C TT1 - poslednja poznata temperatura   
       
C
      ITER=0
      ITER1=0
      INDSK=1
      INDPR=0
!       write(*,*)'INDJOT-pakv3dt',INDJOT
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
C      IF (ITER.EQ.0.AND.MAXCUR.NE.0) IDOKL=2
!        IF (ITER.EQ.0.AND.MAXCUR.NE.0) IDOKL=MAXCUR
C===========================================================================
      IF(NXNASTRAN.EQ.0) THEN
! 	IF (MAXCUR.GT.0) THEN
!           if(ISOLVER.NE.-11) CALL VCLEAR(A(LSK),KC)
!         ELSE
!   samo za modifikovan Njutn
          IF (ITER.EQ.0.and.ISOLVER.NE.-11) CALL VCLEAR(A(LSK),KC)
!         ENDIF
      ENDIF
       
      DO I=1,JEDN
          SILE(I)=0.D0
!           UBRZ(I)=0.D0
      ENDDO 

C
      IF((KKORAK.EQ.1).AND.(ITER.EQ.0))THEN
        IF ((ISOLVER.EQ.-11).OR.(ISOLVER.EQ.1)) THEN
C         MUMPS - full Njutn
          CALL sparseassembler_init(1)
        ENDIF
      ENDIF     
C -------------------------------------
C GLAVNA PETLJA PO ELEMENTIMA
C -------------------------------------
!       IF(NXNASTRAN.EQ.1) GO TO 999
      DO 400 NBREL=1,NET
      IF(elemtip(NBREL).gt.100) THEN
        NTYPE=elemtip(NBREL)/100
        NTIMES=100
      ELSE
        NTYPE=elemtip(NBREL)/10
        NTIMES=10
      ENDIF
      NDIMEL=elemtip(NBREL)-NTIMES*NTYPE
C ODEDIVANJE BROJA GAUSOVIH TACAKA PRILIKOM INTEGRACIJE
      IBRGT=3
      IF (NDIMEL.LE.8) IBRGT=2

      NMAT=NEL(NDIM+1,NBREL)
      NEMA=NEL(NDIM+2,NBREL)

      IF (NMAT.LE.0.OR.NEMA.GT.0) THEN
         WRITE(3,*) 'NBREL,NMAT,NEMA',NBREL,NMAT,NEMA
         WRITE(*,*) 'NBREL,NMAT,NEMA',NBREL,NMAT,NEMA
!          STOP
      ENDIF
      IF (NMAT.LE.0) GOTO 400
      IF (NEMA.GT.0) GOTO 400
      
!        PERM(1,1)=AKONST(3,1,NMAT)
!        PERM(1,2)=0.
!        PERM(1,3)=0.
!        PERM(2,2)=AKONST(3,2,NMAT)
!        PERM(2,1)=0.
!        PERM(2,3)=0.
!        PERM(3,3)=AKONST(3,3,NMAT)
!        PERM(3,1)=0.
!        PERM(3,2)=0.
! 
!       IF (PERM(1,1).LT.1.D-15.OR.PERM(2,2).LT.1.D-15.OR.
!      1    PERM(3,3).LT.1.D-15) THEN
!          WRITE(3,*) 'NBREL,NMAT,NEMA',NBREL,NMAT,NEMA
!          WRITE(*,*) 'NBREL,NMAT,NEMA',NBREL,NMAT,NEMA
!          WRITE(*,*) 'PERM',PERM(1,1),PERM(2,2),PERM(3,3)
!          STOP
!       ENDIF
C=========================================================================
      DO 130 KLM=1,NDIMEL
      CK(1,KLM)=CORD(1,NEL(KLM,NBREL))
      CK(2,KLM)=CORD(2,NEL(KLM,NBREL))
      CK(3,KLM)=CORD(3,NEL(KLM,NBREL))
C
      LM2(KLM)=ID(1,NEL(KLM,NBREL))
 130  CONTINUE
    
      if(NTYPE.eq.2) CALL VLCK2D(CK,CKL,TTE,NDIMEL)
C=======================================================================
      CALL VCLEARD(TT21,NDIMEL)
      CALL VPREB(TT21,TT1,NBREL)
      CALL VCLEARD(TT210,NDIMEL)
      CALL VPREB(TT210,TT10,NBREL)
 
C=======================================================================
      DO 163 K=1,NDIMEL
      DO 162 N=1,NDIMEL
        IF(NSTAC.EQ.0) AM1(K,N)=0.D0    
        AFIFI(K,N)=0.D0    
        AFIFI1(K,N)=0.D0    
!         PUS(K,N)=0.D0    
  162 CONTINUE
  163 CONTINUE
C Fluks I Toplotni izvor:
      DO 199 I=1,NDIM
      RS2(I)=0.
      RS3(I)=0.
 199  CONTINUE
C INCIJALIZACIJA MATRICE SKEF I F36
      DO 260 I=1,NDES
      DO 258 J=1,NDES
      SKEF(I,J)=0.
 258  CONTINUE
      F36(I)=0.
 260  CONTINUE
C
C INTEGRACIJA U GAUSOVIM TACKAMA
C
      NGAUS=0
      ZAPREL=0.
      IBRGT1=IBRGT
      IBRGT2=IBRGT
      if(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) THEN
        IBRGT2=1
        IBRGT1=4
        IF(INDJOT.GT.0) IBRGT1=INDJOT
      ENDIF
      IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) IBRGT2=1
      IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) IBRGT1=1
      IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
        IBRGT2=1
        IBRGT1=3
      ENDIF
      DO 180 I=1,IBRGT1
      DO 170 J=1,IBRGT2
      DO 160 L=1,IBRGT2
      NGAUS=NGAUS+1
      IF(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) THEN
       CALL VGETRST(RT,ST,TT,WTET,NGAUS,IBRGT1)
       CALL VINTERPTETRA(RT,ST,TT,0,0,NDIMEL)
       WDT=WTET(1)*WTET(2)*WTET(3)*DETJ
      ELSE IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) THEN
        CALL VINTERP1D(R(IBRGT,I),0,0,NDIMEL,AJS)  
        WDT=W(IBRGT,I)*thick(NBREL)*AJS
      ELSE IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) THEN
        CALL VINTERP2D48(R(IBRGT,J),S(IBRGT,L),0,0,NDIMEL,CKL,TTE)
        WDT=W(IBRGT,J)*W(IBRGT,L)*DETJ*thick(NBREL)
      ELSE IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
        RT=1./6
        ST=1./6
        WT=1./3
        IF(NGAUS.EQ.2) RT=2./3
        IF(NGAUS.EQ.3) ST=2./3
        CALL VINTERP2D36(RT,ST,0,0,NDIMEL,CKL,TTE)
        WDT = WT**2*thick(NBREL)*DETJ
      ELSE
        CALL VINTER3(R(IBRGT,I),S(IBRGT,J),T(IBRGT,L),0,0)
        WDT=W(IBRGT,I)*W(IBRGT,J)*W(IBRGT,L)*DETJ
      ENDIF
      
      ZAPREL=ZAPREL+WDT
      
      DO 165 K=1,NDIMEL
      DO 164 N=1,NDIMEL
cz S - AM1 nije podeljeno sa dt; podeljeno sa dt prilikom pakovanja u matricu SKEF
!          IF(LINTE.EQ.0) THEN
!  matrica toplotnog kapaciteta  (1.27), (1.38)
!  konstantne u j-ni (1.27) zavise od temperature
         IF(NSTAC.EQ.0.AND.IMASS.EQ.0) 
     1 AM1(K,N)=AM1(K,N)+AKONST(3,4,NMAT)*AKONST(3,5,NMAT)*H(K)*H(N)*WDT
!     (1.24), (1.36) matrica kondukcije 
!  konstantne u j-ni (1.24) zavise od temperature
         AFIFI(K,N)=AFIFI(K,N)+(AKONST(3,1,NMAT)*ZVHX(K)*ZVHX(N)+
     1AKONST(3,2,NMAT)*ZVHY(K)*ZVHY(N)+AKONST(3,3,NMAT)*ZVHZ(K)*ZVHZ(N))
     2*WDT 
!          ELSE
!          ENDIF
!       if(NBREL.eq.1.AND.I.eq.IBRGT) then
!          WRITE(*,*) 'NBREL',NBREL
!          WRITE(*,*) 'H(K)',H(K)
!          WRITE(*,*) 'H(N)',H(N)
!          WRITE(*,*) 'ZVHX(K)',ZVHX(K)
!          WRITE(*,*) 'ZVHY(K)',ZVHY(K)
!          WRITE(*,*) 'ZVHZ(K)',ZVHZ(K)
!          WRITE(*,*) 'WDT',WDT
!          WRITE(*,*) 'AKONST(3,1,NMAT)',AKONST(3,1,NMAT)
!          WRITE(*,*) 'AKONST(3,2,NMAT)',AKONST(3,2,NMAT)
!          WRITE(*,*) 'AKONST(3,3,NMAT)',AKONST(3,3,NMAT)
!           write(3,*) 'AM1(',K,',',N,')',AM1(K,N)
!           write(3,*) 'AFIFI(',K,',',N,')',AFIFI(K,N)
!       ENDIF
C
  164 CONTINUE
  165 CONTINUE
  160 CONTINUE
  170 CONTINUE
  180 CONTINUE
!       WRITE(*,*) 'ZAPREL',ZAPREL
        IF(NSTAC.EQ.0.AND.IMASS.EQ.0) THEN
       write(3,*)'NBREL,ZAPREL',NBREL,ZAPREL
            DO I=1,NDIMEL
            write(3,5000) (AM1(I,K), K=1,NDIMEL)
            ENDDO
        ENDIF

C ako je dijagonalna matrica C
       IF(NSTAC.EQ.0.AND.IMASS.EQ.1) THEN
        IF(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) THEN
       DO  K=1,NDIMEL
       bmnozi=1./4
       IF(NDIMEL.GT.4)bmnozi=1./32
       IF(K.GT.4) bmnozi=7./48
       AM1(K,K)=AM1(K,K)+AKONST(3,4,NMAT)*AKONST(3,5,NMAT)*bmnozi*ZAPREL
       ENDDO
C nije sredjeno za druge elemente       
      ELSE IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) THEN
        CALL VINTERP1D(R(IBRGT,I),0,0,NDIMEL,AJS)  
        WDT=W(IBRGT,I)*thick(NBREL)*AJS
      ELSE IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) THEN
        CALL VINTERP2D48(R(IBRGT,J),S(IBRGT,L),0,0,NDIMEL,CKL,TTE)
        WDT=W(IBRGT,J)*W(IBRGT,L)*DETJ*thick(NBREL)
      ELSE IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
        RT=1./6
        ST=1./6
        WT=1./3
        IF(NGAUS.EQ.2) RT=2./3
        IF(NGAUS.EQ.3) ST=2./3
        CALL VINTERP2D36(RT,ST,0,0,NDIMEL,CKL,TTE)
        WDT = WT**2*thick(NBREL)*DETJ
      ELSE
        CALL VINTER3(R(IBRGT,I),S(IBRGT,J),T(IBRGT,L),0,0)
        WDT=W(IBRGT,I)*W(IBRGT,J)*W(IBRGT,L)*DETJ
      ENDIF
       ENDIF

C=========================================================================
C POVRSINSKE SILE - FLUKS
C=================
!       IF (INDPUN.EQ.1) GOTO 265
!  MAXSIL - broj elemenata na kojima je zadat fluks
      DO 250 JBRPS=1,MAXSIL
      IF (NBREL.EQ.NGPSIL(1,JBRPS)) THEN
      NODE1=NGPSIL(2,JBRPS)
      NODE2=NGPSIL(3,JBRPS)
      NODE3=NGPSIL(4,JBRPS)
      NODE4=NGPSIL(5,JBRPS)
      NODE5=NGPSIL(6,JBRPS)
      NODE6=NGPSIL(7,JBRPS)
      NODE7=NGPSIL(8,JBRPS)
      NODE8=NGPSIL(9,JBRPS)
      N1=NEL(1,NBREL)
      N2=NEL(2,NBREL)
      N3=NEL(3,NBREL)
      N4=NEL(4,NBREL)
      N5=NEL(5,NBREL)
      N6=NEL(6,NBREL)
      N7=NEL(7,NBREL)
      N8=NEL(8,NBREL)
      IVRF=NGPSIL(10,JBRPS)
      CALL VTIMFUN (TABF,FK1,VVREME,ITFMAX(IVRF),IVRF)
      FPS1=POVSIL(1,JBRPS)*FK1
      FPS2=POVSIL(2,JBRPS)*FK1
      FPS3=POVSIL(3,JBRPS)*FK1
      FPS4=POVSIL(4,JBRPS)*FK1
C
      NUMGAU=IBRGT
!  nije dobro uradjen fuks za 1D elemente !!!!      
      IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) THEN
        DO I=1,NUMGAU
          CALL VINTERP1D(R(IBRGT,I),0,0,NDIMEL,AJS)  
          WDTV=W(IBRGT,I)*thick(NBREL)*AJS
          DO K=1,NDIMEL
               RS2(K)=RS2(K)+H(K)*FPS1*WDTV
        ENDDO
        ENDDO
        GOTO 250
      ENDIF
! 
!  dodati zadati fluks za 2D elemente
! 
!  fluks za tetra elemente
      IF(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) THEN
        POVRS=0.
        IJK=8
        IF(NDIMEL.EQ.4) IJK=4
!         write(IIZLAZ,*)'NBREL',NBREL,IFZRAC,IFLUXR(JBRPS)
C ZRACENJE zadato preko FLUKSA        
       IF(IFZRAC.eq.1.and.IFLUXR(JBRPS).eq.1) THEN
C odredjivanje jedinicnog vektrora normale na stranicu elementa
         DO I=1,3
           A12(I)=CORD(I,NODE2)-CORD(I,NODE1)
           A13(I)=CORD(I,NODE3)-CORD(I,NODE1)
!            A12(I)=CORD(I,NODE1)-CORD(I,NODE2)
!            A13(I)=CORD(I,NODE1)-CORD(I,NODE3)
        ENDDO
        CALL VAXBV(A12,A13,DN)
        CALL VJEDV(DN(1),DN(2),DN(3))          
c komponente jedinicnog vektora normale za zracenje
        CALL VTIMFUN (TABF,VNX,VVREME,ITFMAX(INFX),INFX)
        CALL VTIMFUN (TABF,VNY,VVREME,ITFMAX(INFY),INFY)
        CALL VTIMFUN (TABF,VNZ,VVREME,ITFMAX(INFZ),INFZ)
!         write(IIZLAZ,*)'NBREL,DN',NBREL,DN(1),DN(2),DN(3)
!         write(IIZLAZ,*)'INFX,INFY',INFX,INFY,INFZ
C proizvod I*n
        VNI=DN(1)*VNX+DN(2)*VNY+DN(3)*VNZ
!         write(IIZLAZ,*)'VNI',VNX,VNY,VNZ,VNI
        IF(VNI.LT.0) THEN
          DO NGXT = 1,3
             RT=1./6
             ST=1./6
             WT=1./3
             IF(NGXT.EQ.2) RT=2./3
             IF(NGXT.EQ.3) ST=2./3
             CALL VJACTP33(RT,ST,HT,1,NGPSIL,JBRPS,NDIMEL)
             WDTS = WT*DETJS/2
             POVRS=POVRS+WDTS
             IF(IMASS.EQ.0)THEN
             DO K=2,IJK
               DO J=1,NDIMEL
                 IF(NGPSIL(K,JBRPS).EQ.NEL(J,NBREL)) THEN
!   povrsinski fulks  j-na (1.30)                 
!            RS2(J)=RS2(J)+HT(K-1,1)*FS2*WDTS*(-VNI)*RKOREKCIJA*TIME
                  RS2(J)=RS2(J)+HT(K-1,1)*FS2*WDTS*(-VNI)*RKOREKCIJA
                 ENDIF
               ENDDO
             ENDDO
             ENDIF
          ENDDO 
          PFLUXEL(JBRPS)=-FS2*VNI*RKOREKCIJA
!           PFLUXEL(JBRPS)=-FS2*VNI*RKOREKCIJA*TIME
!            write(IIZLAZ,*)'NBREL,FS2',NBREL,FS2,RKOREKCIJA,TIME
!           write(IIZLAZ,*)'PFLUXEL',NGPSIL(1,JBRPS),PFLUXEL(JBRPS)
        ELSE
         GOTO 250
        ENDIF
       ELSE
         DO NGXT = 1,3
             RT=1./6
             ST=1./6
             WT=1./3
             IF(NGXT.EQ.2) RT=2./3
             IF(NGXT.EQ.3) ST=2./3
             CALL VJACTP33(RT,ST,HT,1,NGPSIL,JBRPS,NDIMEL)
             WDTS = WT*DETJS/2
             POVRS=POVRS+WDTS
             IF(IMASS.EQ.0)THEN
             DO K=2,IJK
              DO J=1,NDIMEL
                 IF(NGPSIL(K,JBRPS).EQ.NEL(J,NBREL)) THEN
!   povrsinski fulks  j-na (1.30)                 
                    RS2(J)=RS2(J)+HT(K-1,1)*FS2*WDTS
!           write(3,5000)  HT(K-1,1),FS2,WDTS, RS2(J)
                ENDIF
              ENDDO
            ENDDO
            ENDIF
         ENDDO 
       ENDIF
        IF(IMASS.EQ.1)THEN
         DO K=2,IJK
           IF(K.NE.5) THEN
             DO J=1,NDIMEL
               IF(NGPSIL(K,JBRPS).EQ.NEL(J,NBREL)) THEN
               bmnozi=1./3
C       IF(IJK.GT.4)bmnozi=0.
C       IF(K.GT.5) bmnozi=1./3
              IF(IJK.GT.4)bmnozi=1./32
              IF(K.GT.5) bmnozi=7./48
!  vektor konvekcije (1.31), (1.40)                    
              IF(IFZRAC.eq.1.and.IFLUXR(JBRPS).eq.1) THEN
!   povrsinski fulks  j-na (1.30)                 
                RS2(J)=RS2(J)+FS2*(-VNI)*RKOREKCIJA*bmnozi*POVRS
              ELSE
!   povrsinski fulks  j-na (1.30)                 
                RS2(J)=RS2(J)+FS2*bmnozi*POVRS
              ENDIF
            ENDIF
           ENDDO
          ENDIF
         ENDDO
        ENDIF
       
!  fluks za heksa elemente      
      ELSEIF (elemtip(NBREL).eq.38.or.elemtip(NBREL).eq.320) THEN
      
      DO 225 I=1,NUMGAU
      DO 225 J=1,NUMGAU
C
C KONSTANTNO T (ZETA)
C
      IF((NODE1.EQ.N1.AND.NODE2.EQ.N2.AND.NODE3.EQ.N3.AND.NODE4.EQ.N4).
     1OR.(NODE1.EQ.N2.AND.NODE2.EQ.N3.AND.NODE3.EQ.N4.AND.NODE4.EQ.N1).
     2OR.(NODE1.EQ.N3.AND.NODE2.EQ.N4.AND.NODE3.EQ.N1.AND.NODE4.EQ.N2).
     3OR.(NODE1.EQ.N4.AND.NODE2.EQ.N1.AND.NODE3.EQ.N2.AND.NODE4.EQ.N3).
C
     4OR.(NODE1.EQ.N1.AND.NODE2.EQ.N4.AND.NODE3.EQ.N3.AND.NODE4.EQ.N2).
     5OR.(NODE1.EQ.N4.AND.NODE2.EQ.N3.AND.NODE3.EQ.N2.AND.NODE4.EQ.N1).
     6OR.(NODE1.EQ.N3.AND.NODE2.EQ.N2.AND.NODE3.EQ.N1.AND.NODE4.EQ.N4).
     7OR.(NODE1.EQ.N2.AND.NODE2.EQ.N1.AND.NODE3.EQ.N4.AND.NODE4.EQ.N3))
     8THEN
      CALL VINTER3(R(NUMGAU,I),S(NUMGAU,J),1.D0,3,1)
       GOTO 220
      ENDIF 
C
C KONSTANTNO T (ZETA)
C 
      IF((NODE1.EQ.N5.AND.NODE2.EQ.N6.AND.NODE3.EQ.N7.AND.NODE4.EQ.N8).
     1OR.(NODE1.EQ.N6.AND.NODE2.EQ.N7.AND.NODE3.EQ.N8.AND.NODE4.EQ.N5).
     2OR.(NODE1.EQ.N7.AND.NODE2.EQ.N8.AND.NODE3.EQ.N5.AND.NODE4.EQ.N6).
     3OR.(NODE1.EQ.N8.AND.NODE2.EQ.N5.AND.NODE3.EQ.N6.AND.NODE4.EQ.N7).
C
     5OR.(NODE1.EQ.N5.AND.NODE2.EQ.N8.AND.NODE3.EQ.N7.AND.NODE4.EQ.N6).
     6OR.(NODE1.EQ.N8.AND.NODE2.EQ.N7.AND.NODE3.EQ.N6.AND.NODE4.EQ.N5).
     7OR.(NODE1.EQ.N7.AND.NODE2.EQ.N6.AND.NODE3.EQ.N5.AND.NODE4.EQ.N8).
     8OR.(NODE1.EQ.N6.AND.NODE2.EQ.N5.AND.NODE3.EQ.N8.AND.NODE4.EQ.N7))
     8THEN
      CALL VINTER3(R(NUMGAU,I),S(NUMGAU,J),-1.D0,3,1)
       GOTO 220
      ENDIF 
C 
C KONSTANTNO R (KSI)
C
      IF((NODE1.EQ.N1.AND.NODE2.EQ.N4.AND.NODE3.EQ.N8.AND.NODE4.EQ.N5).
     1OR.(NODE1.EQ.N4.AND.NODE2.EQ.N8.AND.NODE3.EQ.N5.AND.NODE4.EQ.N1).
     2OR.(NODE1.EQ.N8.AND.NODE2.EQ.N5.AND.NODE3.EQ.N1.AND.NODE4.EQ.N4).
     3OR.(NODE1.EQ.N5.AND.NODE2.EQ.N1.AND.NODE3.EQ.N4.AND.NODE4.EQ.N8).
C
     4OR.(NODE1.EQ.N1.AND.NODE2.EQ.N5.AND.NODE3.EQ.N8.AND.NODE4.EQ.N4).
     5OR.(NODE1.EQ.N5.AND.NODE2.EQ.N8.AND.NODE3.EQ.N4.AND.NODE4.EQ.N1).
     6OR.(NODE1.EQ.N8.AND.NODE2.EQ.N4.AND.NODE3.EQ.N1.AND.NODE4.EQ.N5).
     7OR.(NODE1.EQ.N4.AND.NODE2.EQ.N1.AND.NODE3.EQ.N5.AND.NODE4.EQ.N8))
     8THEN
      CALL VINTER3(1.D0,S(NUMGAU,I),T(NUMGAU,J),1,1)
       GOTO 220
      ENDIF 
C 
C KONSTANTNO R (KSI)AFIFI1
C
      IF((NODE1.EQ.N2.AND.NODE2.EQ.N3.AND.NODE3.EQ.N7.AND.NODE4.EQ.N6).
     1OR.(NODE1.EQ.N3.AND.NODE2.EQ.N7.AND.NODE3.EQ.N6.AND.NODE4.EQ.N2).
     2OR.(NODE1.EQ.N7.AND.NODE2.EQ.N6.AND.NODE3.EQ.N2.AND.NODE4.EQ.N3).
     3OR.(NODE1.EQ.N6.AND.NODE2.EQ.N2.AND.NODE3.EQ.N3.AND.NODE4.EQ.N7).
C
     4OR.(NODE1.EQ.N2.AND.NODE2.EQ.N6.AND.NODE3.EQ.N7.AND.NODE4.EQ.N3).
     5OR.(NODE1.EQ.N6.AND.NODE2.EQ.N7.AND.NODE3.EQ.N3.AND.NODE4.EQ.N2).
     6OR.(NODE1.EQ.N7.AND.NODE2.EQ.N3.AND.NODE3.EQ.N2.AND.NODE4.EQ.N6).
     7OR.(NODE1.EQ.N3.AND.NODE2.EQ.N2.AND.NODE3.EQ.N6.AND.NODE4.EQ.N7))
     8THEN
      CALL VINTER3(-1.D0,S(NUMGAU,I),T(NUMGAU,J),1,1)
       GOTO 220
      ENDIF 
C 
C KONSTANTNO S (ETA)
C
      IF((NODE1.EQ.N4.AND.NODE2.EQ.N3.AND.NODE3.EQ.N7.AND.NODE4.EQ.N8).
     1OR.(NODE1.EQ.N3.AND.NODE2.EQ.N7.AND.NODE3.EQ.N8.AND.NODE4.EQ.N4).
     2OR.(NODE1.EQ.N7.AND.NODE2.EQ.N8.AND.NODE3.EQ.N4.AND.NODE4.EQ.N3).
     3OR.(NODE1.EQ.N8.AND.NODE2.EQ.N4.AND.NODE3.EQ.N3.AND.NODE4.EQ.N7).
C
     4OR.(NODE1.EQ.N3.AND.NODE2.EQ.N4.AND.NODE3.EQ.N8.AND.NODE4.EQ.N7).
     5OR.(NODE1.EQ.N4.AND.NODE2.EQ.N8.AND.NODE3.EQ.N7.AND.NODE4.EQ.N3).
     6OR.(NODE1.EQ.N8.AND.NODE2.EQ.N7.AND.NODE3.EQ.N3.AND.NODE4.EQ.N4).
     7OR.(NODE1.EQ.N7.AND.NODE2.EQ.N3.AND.NODE3.EQ.N4.AND.NODE4.EQ.N8))
     8THEN
      CALL VINTER3(R(NUMGAU,I),-1.D0,T(NUMGAU,J),2,1)
       GOTO 220
      ENDIF 
C 
C KONSTANTNO S (ETA)
C
      IF((NODE1.EQ.N1.AND.NODE2.EQ.N2.AND.NODE3.EQ.N6.AND.NODE4.EQ.N5).
     1OR.(NODE1.EQ.N2.AND.NODE2.EQ.N6.AND.NODE3.EQ.N5.AND.NODE4.EQ.N1).
     2OR.(NODE1.EQ.N6.AND.NODE2.EQ.N5.AND.NODE3.EQ.N1.AND.NODE4.EQ.N2).
     3OR.(NODE1.EQ.N5.AND.NODE2.EQ.N1.AND.NODE3.EQ.N2.AND.NODE4.EQ.N6).
C
     4OR.(NODE1.EQ.N2.AND.NODE2.EQ.N1.AND.NODE3.EQ.N5.AND.NODE4.EQ.N6).
     5OR.(NODE1.EQ.N1.AND.NODE2.EQ.N5.AND.NODE3.EQ.N6.AND.NODE4.EQ.N2).
     6OR.(NODE1.EQ.N5.AND.NODE2.EQ.N6.AND.NODE3.EQ.N2.AND.NODE4.EQ.N1).
     7OR.(NODE1.EQ.N6.AND.NODE2.EQ.N2.AND.NODE3.EQ.N1.AND.NODE4.EQ.N5))
     8THEN
      CALL VINTER3(R(NUMGAU,I),1.D0,T(NUMGAU,J),2,1)
       GOTO 220
      ENDIF
 
 220         WDTS=W(NUMGAU,I)*W(NUMGAU,J)*DETJS
!   povrsinski fulks  j-na (1.30)                 
       DO K=1,NDIMEL
         RS2(K)=RS2(K)+H(K)*FS2*WDTS
       ENDDO

C     write(3,*) nbrel,i,j,R(IBRGT,I),S(IBRGT,J)
C     write(3,*) 'fs2,wdts,detjs',fs2,wdts,detjs
C     CALL Vwrr(h,ndim,'h   ')
C     CALL Vwrr(rs2,ndim,'rs2 ')
 225  CONTINUE
      ENDIF
 
      ENDIF
 250  CONTINUE

C      ENDIF
C KRAJ PETLJE AKO JE ITER>0
C=======================================================================
C zadata PRELAZNOST na povrsini elementa
C ----------------------------------------
!      MAXTQE - broj elemenata na kojima je zadata prelaznost
!       WRITE(3,*) "pre JBRPS,MAXTQE",JBRPS,MAXTQE
      IVODE=0
      DO 255 JBRPS=1,MAXTQE
      IF (NBREL.EQ.NELTOK(1,JBRPS)) THEN
      NODE1=NELTOK(2,JBRPS)
      NODE2=NELTOK(3,JBRPS)
      NODE3=NELTOK(4,JBRPS)
      NODE4=NELTOK(5,JBRPS)
      NODE5=NELTOK(6,JBRPS)
      NODE6=NELTOK(7,JBRPS)
      NODE7=NELTOK(8,JBRPS)
      NODE8=NELTOK(9,JBRPS)
      N1=NEL(1,NBREL)
      N2=NEL(2,NBREL)
      N3=NEL(3,NBREL)
      N4=NEL(4,NBREL)
      N5=NEL(5,NBREL)
      N6=NEL(6,NBREL)
      N7=NEL(7,NBREL)
      N8=NEL(8,NBREL)
!  Inidkator da li se radi provera da li se radi Bofang
      IPROM=NELTOK(12,JBRPS)
        IVRFH=NELTOK(10,JBRPS)
        IVRFTOK1=NELTOK(11,JBRPS)
        CALL VTIMFUN (TABF,FK1TOK,VVREME,ITFMAX(IVRFTOK1),IVRFTOK1)
!  za grafiku id funkije koeficijenta prelaza 
      FCONEL(JBRPS)=NELTOK(10,JBRPS)
!  IPROM = 3 proverava se da li je teziste elementa okvaseno      
      IF(IPROM.EQ.1.OR.IPROM.EQ.3) THEN
        IVODE=NELTOK(10,JBRPS)/1000000
        IVRFH1=NELTOK(10,JBRPS)/1000-IVODE*1000
        IVRFH2=NELTOK(10,JBRPS)-IVODE*1000000-IVRFH1*1000
! funkcija temperature okoline         
!         IF(NELTOK(11,JBRPS)/1000.GE.1) THEN
!          IVRFTOK1=NELTOK(11,JBRPS)/1000
!          IVRFTOK2=NELTOK(11,JBRPS)-IVRFTOK1*1000
!         ELSE
!          IVRFTOK1=NELTOK(11,JBRPS)
!          IVRFTOK2=0
!         ENDIF
        IVRFVODE=WATER(3,IVODE)
!  odredjivanje nivoa vode u akumulaciji         
        CALL VTIMFUN (TABF,FWATER,VVREME,ITFMAX(IVRFVODE),IVRFVODE)
!  provera da li je teziste elementa okvaseno        
!         IF(HFACE(JBRPS).GT.FWATER.OR.IVRFTOK2.EQ.0) THEN
!          IF(HFACE(JBRPS).GT.FWATER) THEN
!             CALL VTIMFUN (TABF,FK1TOK,VVREME,ITFMAX(IVRFTOK1),IVRFTOK1)
!         ELSE
!             CALL VTIMFUN (TABF,FK1TOK,VVREME,ITFMAX(IVRFTOK2),IVRFTOK2)
!         ENDIF
!  za grafiku stampanje zadate temperature okoline        
!         TOKOLINEL(JBRPS)=FK1TOK
C zavisna voda - donja voda za Djerdap        
        IF(WATER(2,IVODE).NE.0) THEN
           IPOM111=WATER(2,IVODE)
           IVRFVODE_DEP=WATER(3,IPOM111)
           CALL VTIMFUN (TABF,FWATER_DEP,VVREME,ITFMAX(IVRFVODE_DEP),
     1       IVRFVODE_DEP)
        ENDIF
!  element je okvasen; odredjivanje funkcije koeficijenta prelaza        
      ELSEIF(IPROM.EQ.2) THEN
        IVODE=NELTOK(10,JBRPS)/1000000
        IF(IVODE.LT.1) THEN
          IVRFH=NELTOK(10,JBRPS)
        ELSE
          IVRFH=NELTOK(10,JBRPS)-IVODE*1000000
          IF(WATER(2,IVODE).EQ.0) THEN
            IVRFVODE=WATER(3,IVODE)
            CALL VTIMFUN (TABF,FWATER,VVREME,ITFMAX(IVRFVODE),IVRFVODE)
          ELSE
            IVRFVODE=WATER(3,IVODE)
!               WRITE(3,*) "IVRFVODE",IVRFVODE
            CALL VTIMFUN (TABF,FWATER,VVREME,ITFMAX(IVRFVODE),IVRFVODE)
            IPOM111=WATER(2,IVODE)
            IVRFVODE_DEP=WATER(3,IPOM111)
!              WRITE(3,*) "IPOM111",IPOM111
!              WRITE(3,*) "IVRFVODE_DEP",IVRFVODE_DEP
            IF(IBRANA.EQ.0) IVRFVODE_DEP=WATER(3,WATER(2,IVODE))
            CALL VTIMFUN (TABF,FWATER_DEP,VVREME,ITFMAX(IVRFVODE_DEP),
     1       IVRFVODE_DEP)
          ENDIF
        ENDIF
!         IVRFTOK1=NELTOK(11,JBRPS)
!         CALL VTIMFUN (TABF,FK1TOK,VVREME,ITFMAX(IVRFTOK1),IVRFTOK1)
!  nema Bofanga        
!       ELSE
      ENDIF
!  cvor je okvasen, radi se Bofang. Ne proverava se da li element u vodi/vazduhu     
      IF(IPROM.EQ.2.OR.IPROM.EQ.3.AND.HFACE(JBRPS).LE.FWATER) THEN
!  DJERDAP   
        IVRFSENSOR=SENSOR(3,1)
        CALL VTIMFUN (TABF,FSENSOR(1),VVREME,ITFMAX(IVRFSENSOR)
     1                   ,IVRFSENSOR)
!  Grancarevo (ima dva senzora na vrhu i dnu akumulacije)   
!           IBRANA.EQ.1 - Grancarevo
         IF(IBRANA.EQ.1) THEN
          IF (FWATER.GE.PREKIDNAFR) THEN
             IVRFSENSOR=SENSOR(3,2)
             HSENSOR2=HSENSOR(2)
          ELSE
             IVRFSENSOR=SENSOR(3,3)      
             HSENSOR2=HSENSOR(3)
             IF(FWATER.LT.350.0) THEN
               write(*,*) 'nivo u akumulaciji ispod 350', FWATER
               write(3,*) 'nivo u akumulaciji ispod 350', FWATER
             ENDIF
          ENDIF
         ENDIF
!  kraj Grancarevo
!  odredjivanje temperature u senzorima (termometri u akumulaciji)
        CALL VTIMFUN (TABF,FSENSOR(2),VVREME,ITFMAX(IVRFSENSOR)
     1                   ,IVRFSENSOR)
        CALL VTIMFUN (TABF,FK1TOK,VVREME,ITFMAX(IVRFTOK1),IVRFTOK1)
        IF(WATER(2,IVODE).EQ.0) THEN
!  DJERDAP
!      SUBROUTINE VLININTGV(Y1,Y2,TW1,TW2,TM,IND,Y)  - pakv4.for
         IF(IBRANA.EQ.0) THEN
          CALL VLININTGV(FWATER-HSENSOR(1),FWATER-HSENSOR(2),
     1     FSENSOR(1),FSENSOR(2),FK1TOK,
     2                          0,FWATER-HFACE(JBRPS))
!  GRANCAREVO
         ELSEIF(IBRANA.EQ.1) THEN
!  interpolacija na osnovu dva senzora u akumulaciji (februar - 2018)
          IF(IBOFANG.EQ.1) THEN
            CALL VLININTGV(FWATER-HSENSOR(1),FWATER-HSENSOR2,
     1        FSENSOR(1),FSENSOR(2),FK1TOK,
     2                          0,FWATER-HFACE(JBRPS))
          ELSEIF(IBOFANG.EQ.2) THEN
!  interpolacija na osnovu senzora u akumulaciji (decembar - 2018)
            CALL VLININTGV(0,0,FSENSOR(1),FSENSOR(2),FK1TOK,
     1                               0,FWATER-HFACE(JBRPS))
          ENDIF
         ENDIF
C nema zavisne vode za Grancarevo (kod Djerdapa temp. donje vode zavisi od gornje)
        ELSE
!  donja voda Djerdap        
          CALL VLININTGV(FWATER_DEP-HSENSOR(1),
     1     FWATER_DEP-HSENSOR(2),FSENSOR(1),FSENSOR(2),FK1TOK,
     2                      1,0)
        ENDIF
      ENDIF
C
      TOKOLINEL(JBRPS)=FK1TOK
      NUMGAU=IBRGT
!  proveriti prelaznost za 1D element      
      IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) THEN
        DO I=1,NUMGAU
          CALL VINTERP1D(R(IBRGT,I),0,0,NDIMEL,AJS)  
          WDTV=W(IBRGT,I)*thick(NBREL)*AJS
C odredjivanje temperature u Gausovoj tacki i Tok-T
          TEMP=0.0
          TTT=0.0
          DO K=1,NDIMEL
            TEMP=TEMP+H(K)*TT210(K)
            TTT=TTT+H(K)*FK1TOK
!             TTT=TTT+H(K)*(FK1TOK-TT21(K))
          ENDDO
C odredjivanje koeficijenta prelaznosti HHP u zavisnosti od temperature
          IF(IPROM.EQ.1.AND.IPROM.EQ.3) THEN
            IF(HFACE(JBRPS).GT.FWATER) THEN
                CALL VTIMFUN (TABF,HHP,TEMP,ITFMAX(IVRFH1),IVRFH1)
                FCONEL(JBRPS)=IVRFH1
            ELSE
                CALL VTIMFUN (TABF,HHP,TEMP,ITFMAX(IVRFH2),IVRFH2)
                FCONEL(JBRPS)=IVRFH2
!             write(3,*) "NBREL,IVRFH2",NBREL,IVRFH2
            ENDIF
          ELSE
            CALL VTIMFUN (TABF,HHP,TEMP,ITFMAX(IVRFH),IVRFH)
            FCONEL(JBRPS)=IVRFH
!          write(3,*) "NBREL,IVRFH",NBREL,IVRFH
          ENDIF
!  promeniti ovo za 1D element        
!           DO K=1,NDIMEL
!             RS2(K)=RS2(K)+H(K)*HHP*WDTV*TTT
!           ENDDO
        ENDDO
C  (1.25)  
        DO K=1,NDIMEL
         DO N=1,NDIMEL
          AFIFI1(K,N)=AFIFI1(K,N)+HHP*H(K)*H(N)*WDTV
         ENDDO 
        ENDDO 
        GOTO 255
      ENDIF
C  
!  ugraditi prelaznost za 2D elemente
! 
!  prelaznost za tetra elemente
      IF(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) THEN
        IJK=8
        IF(NDIMEL.EQ.4) IJK=4
        CALL VICLEAR(INDEKS,NDIMEL)
            DO J=1,NDIMEL
             DO K=2,IJK
              IF(K.NE.5)THEN
              IF(NELTOK(K,JBRPS).EQ.NEL(J,NBREL)) then
               INDEKS(K)=J
!                write(3,*) 'K,INDEKS(K),J',K,INDEKS(K),J
               ENDIF
              ENDIF
             ENDDO
            ENDDO
        IBRGTP1=3
        IF(INDJOP.GT.0) IBRGTP1=INDJOP
        POVRS=0.0
        DO NGXT = 1,IBRGTP1
           IF(IBRGTP1.EQ.3)THEN
             RT=1./6
             ST=1./6
             WT=1./3
             IF(NGXT.EQ.2) RT=2./3
             IF(NGXT.EQ.3) ST=2./3
          ELSEIF(IBRGTP1.EQ.4)THEN
             RT=1./3
             ST=1./3
             WT=25./48
             SELECT CASE (NGXT)
             CASE (1) 
               WT=-27./48
             CASE (2) 
               RT=1./5
               ST=3./5
             CASE (3) 
               RT=1./5
               ST=1./5
             CASE (4) 
               RT=3./5
               ST=1./5
             CASE DEFAULT
             END SELECT
          ELSEIF(IBRGTP1.EQ.1)THEN
             RT=1./3
             ST=1./3
             WT=1.
          ENDIF
C
!               CALL VINTERPTETRAF(RT,ST,TT,KFIX,NDIMEL)
              CALL VJACTP33(RT,ST,HT,1,NELTOK,JBRPS,NDIMEL)
              WDTS = WT*DETJS/2
              POVRS=POVRS+WDTS
! C odredjivanje temperature u Gausovoj tacki i Tok-T
             TEMP=0.0
             TTT=0.0
             DO K=2,IJK
              IF(K.NE.5)THEN
!                 DO J=1,NDIMEL
!                   IF(NELTOK(K,JBRPS).EQ.NEL(J,NBREL)) THEN
                    IF(LINTE.EQ.0) THEN
                     TEMP=TEMP+HT(K-1,1)*TT210(INDEKS(K))
!  za (1.40)                    
                     TTT=TTT+HT(K-1,1)*FK1TOK
                    ELSE
                     TEMP=TEMP+HT(K-1,1)*TT21(INDEKS(K))
!  za  (1.31)                    
                     TTT=TTT+HT(K-1,1)*(FK1TOK-TT21(INDEKS(K)))
                    ENDIF
!                   ENDIF
!                 ENDDO
               ENDIF
             ENDDO
! C odredjivanje koeficijenta prelWRITE(*,*) "pre racunanja bofanga"aznosti HHP u zavisnosti od temperature
C odredjivanje koeficijenta prelaznosti HHP u zavisnosti od temperature
!              write(3,*) "NBREL,IVRFH1",NBREL,IVRFH,IVRFH1
!          write(3,*) 'TEMP',TEMP
          IF(IPROM.EQ.1.OR.IPROM.EQ.3) THEN
            IF(HFACE(JBRPS).GT.FWATER) THEN
                CALL VTIMFUN (TABF,HHP,TEMP,ITFMAX(IVRFH1),IVRFH1)
C   stampa ID funkcije za koeficijent prelaznosti
                FCONEL(JBRPS)=IVRFH1
!             write(3,*) "NBREL,IVRFH1",NBREL,IVRFH1
            ELSE
                CALL VTIMFUN (TABF,HHP,TEMP,ITFMAX(IVRFH2),IVRFH2)
                FCONEL(JBRPS)=IVRFH2
!             write(3,*) "NBREL,IVRFH2",NBREL,IVRFH2
            ENDIF
          ELSE
            CALL VTIMFUN (TABF,HHP,TEMP,ITFMAX(IVRFH),IVRFH)
            FCONEL(JBRPS)=IVRFH
!             write(3,*) "NBREL,IVRFH",NBREL,IVRFH
          ENDIF
C
!              DO K=2,8
!               IF(K.NE.5)THEN
! !                  DO J=1,NDIMEL
! !                      IF(NELTOK(K,JBRPS).EQ.NEL(J,NBREL)) THEN
! !  vektor konvekcije (1.31), (1.40)                    
!                   RS2(INDEKS(K))=RS2(INDEKS(K))+HT(K-1,1)*TTT*HHP*WDTS
! !                      ENDIF
! !                  ENDDO
!               ENDIF
!             ENDDO
! (1.25) i (1.37) - matrica konvekcije 
          IF(IMASS.EQ.0)THEN
            DO K=2,IJK
             IF(K.NE.5) THEN
!  vektor konvekcije (1.31), (1.40)                    
                  RS2(INDEKS(K))=RS2(INDEKS(K))+HT(K-1,1)*TTT*HHP*WDTS
!               DO J=1,NDIMEL
               DO N=2,IJK
                 IF(N.NE.5)THEN
!                 IF(NELTOK(N,JBRPS).EQ.NEL(J,NBREL)) THEN
               AFIFI1(INDEKS(K),INDEKS(N))=AFIFI1(INDEKS(K),INDEKS(N))+
     1                  HHP*HT(K-1,1)*HT(N-1,1)*WDTS
!          write(3,*) 'K,N,INDEKS(K),INDEKS(N)',K,N,INDEKS(K),INDEKS(N)
!                  ENDIF
                 ENDIF
               ENDDO
!               ENDDO
             ENDIF
             ENDDO 
           ENDIF
!          
         ENDDO 
        IF(IMASS.EQ.1)THEN
        DO K=2,IJK
          IF(K.NE.5) THEN
         bmnozi=1./3
C       IF(IJK.GT.4)bmnozi=0.
C       IF(K.GT.5) bmnozi=1./3
       IF(IJK.GT.4)bmnozi=1./32
       IF(K.GT.5) bmnozi=7./48
!  vektor konvekcije (1.31), (1.40)                    
          RS2(INDEKS(K))=RS2(INDEKS(K))+HHP*TTT*bmnozi*POVRS
          AFIFI1(INDEKS(K),INDEKS(K))=AFIFI1(INDEKS(K),INDEKS(K))+
     1                  HHP*bmnozi*POVRS
          ENDIF
         ENDDO
        ENDIF
!          write(3,*) 'NBREL,HHP,TTT,TEMP', NBREL,HHP,TTT,TEMP
!          DO I=1,NDIMEL
!           WRITE(3,5000) (AFIFI1(I,J),J=1,NDIMEL)
!          ENDDO
 5000 FORMAT(10(1PE12.4))
!         CALL Vwrr10(AFIFI1,NDIMEL*NDIMEL,'AFIFI1')
!          WRITE(*,*) "posle tetri"
!  prelaznost za heksa elemente         
      ELSEIF (elemtip(NBREL).eq.38.or.elemtip(NBREL).eq.320) THEN
      
      DO 2250 I=1,NUMGAU
      DO 2250 J=1,NUMGAU
C
C KONSTANTNO T (ZETA)
C
      IF((NODE1.EQ.N1.AND.NODE2.EQ.N2.AND.NODE3.EQ.N3.AND.NODE4.EQ.N4).
     1OR.(NODE1.EQ.N2.AND.NODE2.EQ.N3.AND.NODE3.EQ.N4.AND.NODE4.EQ.N1).
     2OR.(NODE1.EQ.N3.AND.NODE2.EQ.N4.AND.NODE3.EQ.N1.AND.NODE4.EQ.N2).
     3OR.(NODE1.EQ.N4.AND.NODE2.EQ.N1.AND.NODE3.EQ.N2.AND.NODE4.EQ.N3).
C
     4OR.(NODE1.EQ.N1.AND.NODE2.EQ.N4.AND.NODE3.EQ.N3.AND.NODE4.EQ.N2).
     5OR.(NODE1.EQ.N4.AND.NODE2.EQ.N3.AND.NODE3.EQ.N2.AND.NODE4.EQ.N1).
     6OR.(NODE1.EQ.N3.AND.NODE2.EQ.N2.AND.NODE3.EQ.N1.AND.NODE4.EQ.N4).
     7OR.(NODE1.EQ.N2.AND.NODE2.EQ.N1.AND.NODE3.EQ.N4.AND.NODE4.EQ.N3))
     8THEN
      CALL VINTER3(R(NUMGAU,I),S(NUMGAU,J),1.D0,3,1)
       GOTO 2200
      ENDIF 
C
C KONSTANTNO T (ZETA)
C 
      IF((NODE1.EQ.N5.AND.NODE2.EQ.N6.AND.NODE3.EQ.N7.AND.NODE4.EQ.N8).
     1OR.(NODE1.EQ.N6.AND.NODE2.EQ.N7.AND.NODE3.EQ.N8.AND.NODE4.EQ.N5).
     2OR.(NODE1.EQ.N7.AND.NODE2.EQ.N8.AND.NODE3.EQ.N5.AND.NODE4.EQ.N6).
     3OR.(NODE1.EQ.N8.AND.NODE2.EQ.N5.AND.NODE3.EQ.N6.AND.NODE4.EQ.N7).
C
     5OR.(NODE1.EQ.N5.AND.NODE2.EQ.N8.AND.NODE3.EQ.N7.AND.NODE4.EQ.N6).
     6OR.(NODE1.EQ.N8.AND.NODE2.EQ.N7.AND.NODE3.EQ.N6.AND.NODE4.EQ.N5).
     7OR.(NODE1.EQ.N7.AND.NODE2.EQ.N6.AND.NODE3.EQ.N5.AND.NODE4.EQ.N8).
     8OR.(NODE1.EQ.N6.AND.NODE2.EQ.N5.AND.NODE3.EQ.N8.AND.NODE4.EQ.N7))
     8THEN
      CALL VINTER3(R(NUMGAU,I),S(NUMGAU,J),-1.D0,3,1)
       GOTO 2200
      ENDIF 
C 
C KONSTANTNO R (KSI)
C
      IF((NODE1.EQ.N1.AND.NODE2.EQ.N4.AND.NODE3.EQ.N8.AND.NODE4.EQ.N5).
     1OR.(NODE1.EQ.N4.AND.NODE2.EQ.N8.AND.NODE3.EQ.N5.AND.NODE4.EQ.N1).
     2OR.(NODE1.EQ.N8.AND.NODE2.EQ.N5.AND.NODE3.EQ.N1.AND.NODE4.EQ.N4).
     3OR.(NODE1.EQ.N5.AND.NODE2.EQ.N1.AND.NODE3.EQ.N4.AND.NODE4.EQ.N8).
C
     4OR.(NODE1.EQ.N1.AND.NODE2.EQ.N5.AND.NODE3.EQ.N8.AND.NODE4.EQ.N4).
     5OR.(NODE1.EQ.N5.AND.NODE2.EQ.N8.AND.NODE3.EQ.N4.AND.NODE4.EQ.N1).
     6OR.(NODE1.EQ.N8.AND.NODE2.EQ.N4.AND.NODE3.EQ.N1.AND.NODE4.EQ.N5).
     7OR.(NODE1.EQ.N4.AND.NODE2.EQ.N1.AND.NODE3.EQ.N5.AND.NODE4.EQ.N8))
     8THEN
      CALL VINTER3(1.D0,S(NUMGAU,I),T(NUMGAU,J),1,1)
       GOTO 2200
      ENDIF 
C 
C KONSTANTNO R (KSI)
C
      IF((NODE1.EQ.N2.AND.NODE2.EQ.N3.AND.NODE3.EQ.N7.AND.NODE4.EQ.N6).
     1OR.(NODE1.EQ.N3.AND.NODE2.EQ.N7.AND.NODE3.EQ.N6.AND.NODE4.EQ.N2).
     2OR.(NODE1.EQ.N7.AND.NODE2.EQ.N6.AND.NODE3.EQ.N2.AND.NODE4.EQ.N3).
     3OR.(NODE1.EQ.N6.AND.NODE2.EQ.N2.AND.NODE3.EQ.N3.AND.NODE4.EQ.N7).
C
     4OR.(NODE1.EQ.N2.AND.NODE2.EQ.N6.AND.NODE3.EQ.N7.AND.NODE4.EQ.N3).
     5OR.(NODE1.EQ.N6.AND.NODE2.EQ.N7.AND.NODE3.EQ.N3.AND.NODE4.EQ.N2).
     6OR.(NODE1.EQ.N7.AND.NODE2.EQ.N3.AND.NODE3.EQ.N2.AND.NODE4.EQ.N6).
     7OR.(NODE1.EQ.N3.AND.NODE2.EQ.N2.AND.NODE3.EQ.N6.AND.NODE4.EQ.N7))
     8THEN
      CALL VINTER3(-1.D0,S(NUMGAU,I),T(NUMGAU,J),1,1)
       GOTO 2200
      ENDIF 
C 
C KONSTANTNO S (ETA)
C
      IF((NODE1.EQ.N4.AND.NODE2.EQ.N3.AND.NODE3.EQ.N7.AND.NODE4.EQ.N8).
     1OR.(NODE1.EQ.N3.AND.NODE2.EQ.N7.AND.NODE3.EQ.N8.AND.NODE4.EQ.N4).
     2OR.(NODE1.EQ.N7.AND.NODE2.EQ.N8.AND.NODE3.EQ.N4.AND.NODE4.EQ.N3).
     3OR.(NODE1.EQ.N8.AND.NODE2.EQ.N4.AND.NODE3.EQ.N3.AND.NODE4.EQ.N7).
C
     4OR.(NODE1.EQ.N3.AND.NODE2.EQ.N4.AND.NODE3.EQ.N8.AND.NODE4.EQ.N7).
     5OR.(NODE1.EQ.N4.AND.NODE2.EQ.N8.AND.NODE3.EQ.N7.AND.NODE4.EQ.N3).
     6OR.(NODE1.EQ.N8.AND.NODE2.EQ.N7.AND.NODE3.EQ.N3.AND.NODE4.EQ.N4).
     7OR.(NODE1.EQ.N7.AND.NODE2.EQ.N3.AND.NODE3.EQ.N4.AND.NODE4.EQ.N8))
     8THEN
      CALL VINTER3(R(NUMGAU,I),-1.D0,T(NUMGAU,J),2,1)
       GOTO 2200
      ENDIF 
C 
C KONSTANTNO S (ETA)
C
      IF((NODE1.EQ.N1.AND.NODE2.EQ.N2.AND.NODE3.EQ.N6.AND.NODE4.EQ.N5).
     1OR.(NODE1.EQ.N2.AND.NODE2.EQ.N6.AND.NODE3.EQ.N5.AND.NODE4.EQ.N1).
     2OR.(NODE1.EQ.N6.AND.NODE2.EQ.N5.AND.NODE3.EQ.N1.AND.NODE4.EQ.N2).
     3OR.(NODE1.EQ.N5.AND.NODE2.EQ.N1.AND.NODE3.EQ.N2.AND.NODE4.EQ.N6).
C
     4OR.(NODE1.EQ.N2.AND.NODE2.EQ.N1.AND.NODE3.EQ.N5.AND.NODE4.EQ.N6).
     5OR.(NODE1.EQ.N1.AND.NODE2.EQ.N5.AND.NODE3.EQ.N6.AND.NODE4.EQ.N2).
     6OR.(NODE1.EQ.N5.AND.NODE2.EQ.N6.AND.NODE3.EQ.N2.AND.NODE4.EQ.N1).
     7OR.(NODE1.EQ.N6.AND.NODE2.EQ.N2.AND.NODE3.EQ.N1.AND.NODE4.EQ.N5))
     8THEN
      CALL VINTER3(R(NUMGAU,I),1.D0,T(NUMGAU,J),2,1)
       GOTO 2200
      ENDIF
 2200  WDTS=W(NUMGAU,I)*W(NUMGAU,J)*DETJS
       TTT=0.0
       TEMP=0.0
       DO K=1,NDIMEL
         IF(LINTE.EQ.0) THEN
            TEMP=TEMP+H(K)*TT210(K)
!  za (1.40)                    
             TTT=TTT+H(K)*FK1TOK
           ELSE
            TEMP=TEMP+H(K)*TT21(K)
!  za  (1.31)                    
             TTT=TTT+H(K)*(FK1TOK-TT21(K))
           ENDIF
       ENDDO
C odredjivanje koeficijenta prelaznosti HHP u zavisnosti od temperature
          CALL VTIMFUN (TABF,HHP,TEMP,ITFMAX(IVRFH),IVRFH)
C vektor konvekcije (1.31), (1.40)
       DO K=1,NDIMEL
         RS2(K)=RS2(K)+H(K)*HHP*TTT*WDTS
       ENDDO
! (1.25) i (1.37) matrica konvekcije 
         DO K=1,NDIMEL
          DO N=1,NDIMEL
           AFIFI1(K,N)=AFIFI1(K,N)+HHP*H(K)*H(N)*WDTS
          ENDDO 
         ENDDO 
C          
 2250  CONTINUE
      ENDIF

      ENDIF
 255  CONTINUE

C kraj dela za prelaznost

C PAKOVANJE POVRSINSKIH SILA SA DESNE STRANE
 265   DO I=1,NDIMEL
         F36(I)=F36(I)+RS2(I)
       ENDDO
C=========================================================================
      DO I=1,NDIMEL
      DO J=1,NDIMEL
       AK(I,J)=0.D0
      ENDDO
      ENDDO
C PAKOVANJE MATRICA AMUU,AMUW,AMWW,AKUU,AKUW,AKWW U MATRICE AMM,AC,AK
      DO I=1,NDIMEL
      DO J=1,NDIMEL
        AK(I,J)=AFIFI(I,J)+AFIFI1(I,J)
! (1.28)      
        IF (NSTAC.EQ.0) AK(I,J)=AK(I,J)+AM1(I,J)/TIME
      ENDDO
      ENDDO
C============================================================================
C PAKOVANJE MATRICE AK U MATRICU SKEF
!       WRITE(*,*) "pre pakovanja AK matrice u skef"
      DO 263 I=1,NDIMEL
      DO 262 J=1,NDIMEL
        SKEF(I,J)=AK(I,J)
  262 CONTINUE
  263 CONTINUE
C============================================================================
C PAKOVANJE Kfifi*fi(i-1) SA DESNE STRANE
      IF (NSTAC.EQ.0) THEN
!  VEKTOR TOPLOTNOG KAPACITETA (1.28) ili (1.39)      
       DO I=1,NDIMEL
        DO J=1,NDIMEL
         IF(LINTE.EQ.0) THEN
! (1.39)      
           F36(I)=F36(I)+AM1(I,J)*TT210(J)/TIME
         ELSE
! (1.28)      
          F36(I)=F36(I)-(AM1(I,J)*(TT21(J)-TT210(J)))/TIME
         ENDIF
        ENDDO
       ENDDO
      ENDIF
! 
!  (1.33) vektor kondukcije      
! 
      IF (LINTE.EQ.1) THEN
       DO I=1,NDIMEL
        DO J=1,NDIMEL
          F36(I)=F36(I)-AFIFI(I,J)*TT21(J)
        ENDDO
       ENDDO
      ENDIF
C============================================================================
C VEKTOR PROTOKA PO CVOROVIMA
!       DO I=1,NDIMEL
!        N=NEL(I,NBREL)
!        JJ=ID(1,N)
! !        IF(JJ.GT.0) THEN
! ! 
! !          POT=0.D0
! ! !         DO J=1,NDIMEL
! ! ! 
! ! ! !          POT=POT+(AFIFI(I,J)-AFIFI1(I,J))*TT21(J)
! ! ! cz uporediti resenja sa ovim uslovom za nestacionarni protok
! ! ! cz proveriti resenja za zasic=0.
! ! ! !       	IF (NSTAC.NE.1) THEN
! ! ! !           POT=POT+((AM1(I,J)-PUS(I,J))*(TT21(J)-TT210(J)))/TIME
! ! ! !         ENDIF
! ! ! 
! ! !         ENDDO
! !          UBRZ(JJ)=UBRZ(JJ)+POT
! !        ENDIF
!       ENDDO
C PAKOVANJE MATRICE ELEMENTA I VEKTRA DESNE STRANE U SISTEM
!        WRITE(*,*) "pre pakovanja desne strane"
       CALL VPSKEFN(SKEF,SKEFN,NDIMEL)
!        WRITE(*,*) "posle PSKEFN"
       CALL VSPAKDE(SILE,F36,LM2,NDIMEL)
!        WRITE(*,*) "posle SPAKDE"
      IF (ISOLVER.EQ.0) THEN
        IF (MAXCUR.GT.0) THEN
          CALL VSPAKUJ(A(LSK),MAXA,SKEFN,LM2,NDIMEL)
        ELSE
cz ovaj uslov ne valja jer ne radi levu stranu za full newton      
C        IF (ITER.EQ.0.OR.NJRAP.EQ.2) 
          IF (ITER.EQ.0) 
     +                  CALL VSPAKUJ(A(LSK),MAXA,SKEFN,LM2,NDIMEL)
        ENDIF
      ELSEIF ((ISOLVER.EQ.-11).OR.(ISOLVER.EQ.1)) THEN
C      IMSL Sparse OR MUMPS
        IF(ITER.EQ.0) THEN
C         CALL VSPARSE_ELEM_LEFT(SKEF,NEL,ID,NDES,NBREL,NDIM,1)
          CALL sparseassembler_addelemmatrix(NDES,LM2,SKEF)
        ENDIF
      ENDIF
!       WRITE(*,*) "pre 400 continue"
C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE
!       WRITE(*,*) "posle petlje po elementima"
        CALL DATE_AND_TIME(VALUES=Dtime)
        WRITE(*,*) 'posle petlje po elementima', (Dtime(i),i=5,7)
        WRITE(3,*) 'posle petlje po elementima', (Dtime(i),i=5,7)
C=======================================================================
C CVOROVI NA LINIJI PROCURIVANJA       
C=======================================================================
  	 IF (MAXCUR.GT.0) THEN
        CALL VCURI1(CORD,A(LSK),SILE,MAXA,TT1,ID,ICUR)
C       ELSE
C        IF(ITER.EQ.0) CALL VCURI1(CORD,A(LSK),SILE,MAXA,TT1,ID,ICUR)
       ENDIF
C==========================================================================
C MUMPS alociranje matrice krutosti
3456  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
           CALL VSPAKUA(A(LSK),MAXA,A(LSKE),LM2,NDES,0,
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
c            CALL VSPARSE_BUILD(ISPARSE_N,ISPARSE_NZ,
c     1                DSPARSE_A,ISPARSE_ROWS,ISPARSE_COLS)
          if(myid.eq.0) CALL sparseassembler_getsparse(nonzeros,
     1                 rows,columns,stiff)

          if(myid.eq.0) CALL sparseassembler_kill()

          IF(KKORAK.EQ.1) THEN
              stiff_n = JEDN
C Komentirao Milos 20080121
C            CALL VMUMPS_INIT(ISPARSE_N,ISPARSE_NZ)
          ENDIF
         ENDIF
        ENDIF
C==========================================================================
! 	IF (MAXCUR.GT.0) THEN
! 	 IF(ISOLVER.EQ.0) THEN
!           if(myid.eq.0) CALL VZADLEV(A(LSK),MAXA,A(LNZADJ),NUMZAD,NZAD)
!          ELSE
!           if(myid.eq.0) CALL VSPARS_ZADLEV(A(LNZADJ),NUMZAD,NZAD)
!          ENDIF
!          CALL VRESEN(A(LSK),TT10,MAXA,JEDN,1)
!         ELSE
cz proveriti i ovde uslov za full newton        
C        IF (ITER.EQ.0.OR.NJRAP.EQ.2) THEN
        IF (ITER.EQ.0) THEN
cz proveriti rad sa blokovima i zadlev
          if(myid.eq.0) then
          IF (NBLOCK.GT.1) THEN
            LLM =LRAD
            LSKE=LLM+100
            CALL VSPAKUA(A(LSK),MAXA,A(LSKE),LM2,NDES,0,
     &          A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,H,NZAD)
            CLOSE (ISCRC,STATUS='KEEP')

          ENDIF
          endif
        CALL DATE_AND_TIME(VALUES=Dtime)
        WRITE(3,*) 'pre ZADLEV', (Dtime(i),i=5,7)
          
          IF(ISOLVER.EQ.0) THEN
          if(myid.eq.0) CALL VZADLEV(A(LSK),MAXA,A(LNZADJ),NUMZAD,NZAD)
          ELSE
            if(myid.eq.0) CALL VSPARS_ZADLEV(A(LNZADJ),NUMZAD,NZAD,
     1                                  CORD,VVREME,TABF,NTABFT,ITFMAX)
          ENDIF
        CALL DATE_AND_TIME(VALUES=Dtime)
        WRITE(3,*) 'posle ZADLEV', (Dtime(i),i=5,7)
          CALL VRESEN(A(LSK),TT10,MAXA,JEDN,1)
        ENDIF
        CALL DATE_AND_TIME(VALUES=Dtime)
        WRITE(3,*) 'posle RESEN', (Dtime(i),i=5,7)
!         ENDIF
       IF (myid.ne.0) goto 4567
       TOL=1.0D-7
!        WRITE(*,*) "pre petlje po zadatim vrednostima"
C=======================================================================
       DO 410 I=1,NUMZAD
c       write(iizlaz,*)'fk1=',fk1
c       do ii=1,2
c        write(iizlaz,*)TABF(1,1,Ii),TABF(2,1,Ii)
c       enddo
       IF(NZAD(2,I).NE.0) THEN
        CALL VTIMFUN(TABF,FK1,VVREME,ITFMAX(NZAD(3,I)),NZAD(3,I))
       ELSE
        IF(NZAD(3,I).EQ.1) THEN
          INTERPAXIS=1
          CALL VPRED_INTERP(I,INTERPAXIS,NUMAXISPTSX,INTAXISPOINTX,
     1                    ITFMAX,NZAD,XAXISPTCORD,CORD,TABF,VVREME,FK1)
        ENDIF
        IF(NZAD(3,I).EQ.2) THEN
          INTERPAXIS=2
          CALL VPRED_INTERP(I,INTERPAXIS,NUMAXISPTSY,INTAXISPOINTY,
     1                    ITFMAX,NZAD,YAXISPTCORD,CORD,TABF,VVREME,FK1)
        ENDIF
        IF(NZAD(3,I).EQ.3) THEN
          INTERPAXIS=3
          CALL VPRED_INTERP(I,INTERPAXIS,NUMAXISPTSZ,INTAXISPOINTZ,
     1                    ITFMAX,NZAD,ZAXISPTCORD,CORD,TABF,VVREME,FK1)
        ENDIF
!        IF(NZAD(1,I).EQ.834) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!        IF(NZAD(1,I).EQ.837) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!        IF(NZAD(1,I).EQ.932) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!        IF(NZAD(1,I).EQ.1184) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!        IF(NZAD(1,I).EQ.1468) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!        IF(NZAD(1,I).EQ.1678) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!        IF(NZAD(1,I).EQ.6045) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!        IF(NZAD(1,I).EQ.6045) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!        IF(NZAD(1,I).EQ.6683) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!        IF(NZAD(1,I).EQ.8826) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!        IF(NZAD(1,I).EQ.8818) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!  Lamela 10
!         IF(NZAD(1,I).EQ.33564) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.38505) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.33594) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.30373) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.33821) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.33853) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.33855) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.33556) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.29005) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.33862) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.33537) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.28997) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.15242) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.50789) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.50877) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.15557) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.16378) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
!         IF(NZAD(1,I).EQ.13994) write(3,*)'NZAD(1,I),FK1',NZAD(1,I),FK1
        ENDIF
c       write(iizlaz,*)'fk1=',fk1, vvreme,ITFMAX(NZAD(3,I)),NZAD(3,I)
c       do ii=1,2
c        write(iizlaz,*)TABF(1,1,Ii),TABF(2,1,Ii)
c       enddo
C===========================================================================
C KONSTANTNI POTENCIJALI
!         IF( (NZAD(2,I).EQ.1.OR.(NZAD(2,I).EQ.3.AND.KKORAK.EQ.1))
!      &                 .AND.(ITER.EQ.0).AND.(ID(1,NZAD(1,I)).GT.0)) THEN
cz oduzimanje zadate vrednosti od prethodno izracunate             
         IF(LINTE.EQ.0) THEN
            SILE(ID(1,NZAD(1,I)))=1.D035*ZADVRE(I)*FK1
         ELSE
             SILE(ID(1,NZAD(1,I)))=1.D035*
     1                   (ZADVRE(I)*FK1-TT1(ID(1,NZAD(1,I))))
         ENDIF
!         write(iizlaz,*)'desna',I,NZAD(1,I),SILE(ID(1,NZAD(1,I)))
!        ENDIF
C===========================================================================
  410  CONTINUE
 4567  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       CALL VRESEN(A(LSK),SILE,MAXA,JEDN,2)
! !  za Milosa
!          IF (KKORAK.eq.1) then
!             DO I = 1, JEDN
!                WRITE(1005,*) SILE(i)
!             END DO
!          ENDIF
!  kraj stampe za Milosa
      IF (myid.ne.0) goto 2346
      CALL DATE_AND_TIME(VALUES=Dtime)
      WRITE(3,*) 'posle RESEN-a 2', (Dtime(i),i=5,7)
c         CALL Vwrr(SILE,JEDN,'s-k')
c     OVDE UCITATI RESENJA SA DISKA U VEKTOR SILE
c
!   999 IF(NXNASTRAN.EQ.1) THEN
! C
!             NXNAS=99    
!             OPEN (NXNAS,FILE='POTNXNAS',STATUS='OLD',FORM='FORMATTED',
!      1            ACCESS='SEQUENTIAL')
!             DO I=1,JEDN
!                READ(NXNAS,*) IDUM,SILE(I)
!                write(3,4000) IDUM,SILE(I)
!  4000 format(i5,f15.5)               
!                IF(IDUM.NE.I) STOP 'IDUM.NE.I  PAKV3D.FOR'
!             ENDDO
!             CLOSE(NXNAS,STATUS='KEEP')
!       ENDIF
! 
C
           DO 440 I=1,JEDN
             IF(LINTE.EQ.0) THEN
                TT1(I)=SILE(I)
             ELSE
                TT1(I)=TT1(I)+SILE(I)
             ENDIF
 440       CONTINUE
!  Sneza 31.03.2017.
!  Ako radi stacionarnu linearnu analizu ima samo jedan prolaz, ne proverava se konvergenicja
!         WRITE(*,*)'period,korak=',NNPER,KORAK
 2346  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!         WRITE(*,*)'MPI_BCAST, linte',LINTE
      CALL MPI_BCAST(LINTE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!         WRITE(*,*)'MPI_BCAST, MAXIT',MAXIT
      CALL MPI_BCAST(MAXIT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       IF(LINTE.EQ.0) goto 491
! 
!       WRITE(*,*)'ITER= ',ITER
!       WRITE(*,*)'PERIOD= ',NNPER
!       WRITE(*,*)'STEP= ',KKORAK
!       IF(NXNASTRAN.EQ.1) go to 491
C
        WRITE(*,*)'KONVTF'
      CALL VKONVTF(TT1,SILE,KONVV2,1,ID,ITER)
C
!  2346  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(KONVV2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
C
        IF (KKORAK.GT.1) THEN
C
C        CALL VCKONV1(ICUR,TT1,CORD,ID,KONV1,UBRZ)
c         proba da nasilno ne radi iteracije
cz          nasilu=0
cz          if(nstac.eq.1.and.iter.eq.0) nasilu=0
cz          write(*,*) ' nasilu',nasilu
cz          if(nasilu.eq.0) then
          IF (KONVV2*KONV1.EQ.0) THEN
            ITER=ITER+1
            WRITE(IIZLAZ,*)'ITER=',ITER
            GO TO 100
          ENDIF
        endif 
        CALL DATE_AND_TIME(VALUES=Dtime)
        WRITE(3,*) 'pre racunanja gradijenta', (Dtime(i),i=5,7)

        IF ((NJUTN.NE.0).AND.(INDOPT.EQ.0)) RETURN

C======================================================================= 
C RACUNANJE GRADIJENTA Temperature
C=======================================================================
C  PETLJA PO ELEMENTIMA
C
 491  IF (myid.ne.0) goto 2345
!        WRITE(*,*)'pre gradijenta'
      DO 500 NBREL=1,NET
      IF(elemtip(NBREL).gt.100) THEN
        NTYPE=elemtip(NBREL)/100
        NTIMES=100
      ELSE
        NTYPE=elemtip(NBREL)/10
        NTIMES=10
      ENDIF
      NDIMEL=elemtip(NBREL)-NTIMES*NTYPE
      
!       WRITE(*,*)'NBREL=',NBREL
      
      VZAPR=0.
      NMAT=NEL(NDIM+1,NBREL)
      NEMA=NEL(NDIM+2,NBREL)

      IF (NMAT.LE.0) GOTO 500
      IF (NEMA.GT.0) GOTO 500
 
      CALL VCLEARD(TT21,NDIMEL)
      CALL VPREB(TT21,TT1,NBREL)
      CALL VCLEARD(TT210,NDIMEL)
      CALL VPREB(TT210,TT10,NBREL)
C
      DO KLM=1,NDIMEL
       CK(1,KLM)=CORD(1,NEL(KLM,NBREL))
       CK(2,KLM)=CORD(2,NEL(KLM,NBREL))
       CK(3,KLM)=CORD(3,NEL(KLM,NBREL))
      ENDDO
      
      if(NTYPE.eq.2) CALL VLCK2D(CK,CKL,TTE,NDIMEL)

       DO KL=1,NDIMEL
        FZAP(1,KL)=0.D0
        FZAP(2,KL)=0.D0
        FZAP(3,KL)=0.D0
       ENDDO
C
C INDIKATOR NJUTN-KOTESOVE INTEGRACIJE (0-NE,1-DA)
      KOTES=0
      NGAUS=0
      IBRGT=2
      IBRKT=2
      IF(KOTES.EQ.1) IBRKT=3
      
        IBRKT1=IBRKT
        IBRKT2=IBRKT
        if(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) THEN
          IBRKT2=1
          IBRKT1=4
        ENDIF
        IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) IBRKT2=1
        IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) IBRKT1=1
        IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
          IBRKT2=1
          IBRKT1=3
        ENDIF
        DO  320 I=1,IBRKT1
        DO  320 J=1,IBRKT2
        DO  320 K=1,IBRKT2
        NGAUS=NGAUS+1
        PORNIEL(NBREL,NGAUS)=0.D0
        
	IF(KOTES.EQ.1) THEN
          IF(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) then
            CALL VGETRST(RT,ST,TT,WTET,NGAUS,IBRKT1)
            CALL VINTERPTETRA(RT,ST,TT,0,0,NDIMEL)
            WDT=WTET(1)*WTET(2)*WTET(3)*DETJ
          ELSE IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) THEN
            CALL VINTERP1D(RK(IBRKT1,I),0,0,NDIMEL,AJS)
            WDT=WK(IBRKT1,I)*thick(NBREL)*AJS
          ELSE IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) THEN
            CALL VINTERP2D48(R(IBRGT,J),S(IBRGT,K),0,0,NDIMEL,CKL,TTE)
            WDT=W(IBRGT,J)*W(IBRGT,K)*DETJ*thick(NBREL)
          ELSE IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
            RT=1./6
            ST=1./6
            WT=1./3
            IF(NGAUS.EQ.2) RT=2./3
            IF(NGAUS.EQ.3) ST=2./3
            CALL VINTERP2D36(RT,ST,0,0,NDIMEL,CKL,TTE)
            WWDT = WT*thick(NBREL)*DETJ/2
          ELSE
            CALL VINTER3(RK(IBRGT,I),SK(IBRGT,J),TK(IBRGT,K),0,0)
            WDT=WK(IBRKT,I)*WK(IBRKT,J)*WK(IBRKT,K)*DETJ
          ENDIF
	ELSE
          IF(elemtip(NBREL).eq.34.or.elemtip(NBREL).eq.310) then
            CALL VGETRST(RT,ST,TT,WTET,NGAUS,IBRKT1)
            CALL VINTERPTETRA(RT,ST,TT,0,0,NDIMEL)
            WDT=WTET(1)*WTET(2)*WTET(3)*DETJ
          ELSE IF(elemtip(NBREL).eq.12.or.elemtip(NBREL).eq.13) THEN
            CALL VINTERP1D(R(IBRKT1,I),0,0,NDIMEL,AJS)
            WDT=W(IBRKT1,I)*thick(NBREL)*AJS
          ELSE IF(elemtip(NBREL).eq.24.or.elemtip(NBREL).eq.28) THEN
!             CALL VINTERP2D48(R(IBRGT,I),S(IBRGT,J),0,0,NDIMEL,CKL,TTE)
!             WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ*thick(NBREL)
              CALL VINTERP2D48(R(IBRGT,J),S(IBRGT,K),0,0,NDIMEL,CKL,TTE)
              WDT=W(IBRGT,J)*W(IBRGT,K)*DETJ*thick(NBREL)
          ELSE IF(elemtip(NBREL).eq.23.or.elemtip(NBREL).eq.26) THEN
            RT=1./6
            ST=1./6
            WT=1./3
            IF(NGAUS.EQ.2) RT=2./3
            IF(NGAUS.EQ.3) ST=2./3
            CALL VINTERP2D36(RT,ST,0,0,NDIMEL,CKL,TTE)
            WDT = WT*thick(NBREL)*DETJ/2
          ELSE
            CALL VINTER3(R(IBRGT,I),S(IBRGT,J),T(IBRGT,K),0,0)
            WDT=W(IBRKT,I)*W(IBRKT,J)*W(IBRKT,K)*DETJ
          ENDIF
	ENDIF
        VZAPR=VZAPR+WDT
C
C=======================================================================
C RACUNANJE GRADIJENTA TEMPERATURE:
C
      GRADFX=0.D0
      GRADFY=0.D0
      GRADFZ=0.D0
C RACUNANJE GRADIJENTA TEMPERATURE
          DO II=1,NDIMEL
            GRADFX=GRADFX+ZVHX(II)*TT21(II)
            GRADFY=GRADFY+ZVHY(II)*TT21(II)
            GRADFZ=GRADFZ+ZVHZ(II)*TT21(II)
          ENDDO
C
C         
!       ENDIF
         AKX=AKONST(3,1,NMAT)
         AKY=AKONST(3,2,NMAT)
         AKZ=AKONST(3,3,NMAT)
!       ENDIF
C GRADIJENT U GAUS TACKI
       AJX1(NGAUS)=GRADFX
       AJY1(NGAUS)=GRADFY
       AJZ1(NGAUS)=GRADFZ
       GG(1,NBREL,NGAUS)=GRADFX
       GG(2,NBREL,NGAUS)=GRADFY
       GG(3,NBREL,NGAUS)=GRADFZ
C BRZINA U GAUS TACKI
       AJX(NGAUS)=-AKX*GRADFX
       AJY(NGAUS)=-AKY*GRADFY
       AJZ(NGAUS)=-AKZ*GRADFZ
       VG(1,NBREL,NGAUS)=-AKX*GRADFX
       VG(2,NBREL,NGAUS)=-AKY*GRADFY
       VG(3,NBREL,NGAUS)=-AKZ*GRADFZ
C
  320   CONTINUE
!       CASE DEFAULT
!         WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NBREL)
!       END SELECT
  
      IMCVOR=0
      IVAL=0
      SELECT CASE(elemtip(NBREL))
      CASE(12)
        IVAL=2
      CASE(13)
        IVAL=2
        IMCVOR=1
      CASE(23)
        IVAL=3
      CASE(26)
        IVAL=3
        IMCVOR=1
      CASE(24)
        IVAL=4
      CASE(28)
        IVAL=4
        IMCVOR=1
      CASE(34)
        IVAL=4
      CASE(38)
        IVAL=8
      CASE(310)
        IVAL=4
        IMCVOR=1
      CASE(320)
        IVAL=8
        IMCVOR=1
      CASE DEFAULT
!         WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NBREL)
      END SELECT
C       
 	DO KK=1,NDIMEL
	  NCVOR=NEL(KK,NBREL)
	  KK1=NRED(KK)
	  SELECT CASE(elemtip(NBREL))
          CASE(12:13)
            KK1=NRED1(KK)
          CASE(23,26)
            KK1=NRED26(KK)
          CASE(24,28)
            KK1=NRED28(KK)
          CASE(34,310)
            KK1=NREDT10(KK)
          CASE DEFAULT
!          WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NBREL)
          END SELECT
	  IF(KOTES.EQ.1) KK1=NREDK(KK)
	  IF(KK.LE.IVAL.OR.IMCVOR.EQ.0) THEN
C GRADIJENT Temperature U CVOROVIMA
	  GRADJN(1,NCVOR)=GRADJN(1,NCVOR)+AJX1(KK1)/IVECT(NCVOR)
	  GRADJN(2,NCVOR)=GRADJN(2,NCVOR)+AJY1(KK1)/IVECT(NCVOR)
	  GRADJN(3,NCVOR)=GRADJN(3,NCVOR)+AJZ1(KK1)/IVECT(NCVOR)
          ELSE IF(KK.GT.IVAL.AND.IMCVOR.EQ.1) THEN
          KK2=KK1/10
          KK3=KK1-KK2*10
C GRADIJENT Temperature U CVOROVIMA
          GRADJN(1,NCVOR)=GRADJN(1,NCVOR)+(AJX1(KK2)/IVECT(NCVOR))/2
     1                                   +(AJX1(KK3)/IVECT(NCVOR))/2
          GRADJN(2,NCVOR)=GRADJN(2,NCVOR)+(AJY1(KK2)/IVECT(NCVOR))/2
     1                                   +(AJY1(KK3)/IVECT(NCVOR))/2
          GRADJN(3,NCVOR)=GRADJN(3,NCVOR)+(AJZ1(KK2)/IVECT(NCVOR))/2
     1                                   +(AJZ1(KK3)/IVECT(NCVOR))/2
          
          ENDIF
	ENDDO
C======================================================================= 
 500  CONTINUE
C
C
        CALL DATE_AND_TIME(VALUES=Dtime)
      WRITE(*,*) 'vreme pre stampanja', (Dtime(i),i=5,7)
      WRITE(3,*) 'vreme pre stampanja', (Dtime(i),i=5,7)
      WRITE(IIZLAZ,*)'PERIOD NUMBER,KKORAK ',NNPER, KKORAK        
      WRITE(*,*)'PERIOD NUMBER,KKORAK ',NNPER, KKORAK        
      IF (KKORAK.GT.0) THEN
! otvaranje fajla za stampanje zapreminskih sila
CE WRITING TEMPERATURES FOR TERMOMECHANIC
C      CALL VRESTEL(VVREME,TT1,NPT,19) 
C NPRINT definise korak u kome se stampaju temperature za PAKS
!       write(*,*)'NPRINT,KORAK',NPRINT,KORAK
!       IF(NPRINT.EQ.1)THEN
!         CALL VRESTEL(VVREME,TT1,NPT,IDJSTAMP1,NP3D1,ISNUMER)
!       ELSEIF(KORAK.EQ.NPRINT) THEN
C Sneza 24.06.2011 dodato da stampa rezultate samo za odredjene cvorove
C stampanje potencijala
!         write(*,*)'stampanje potencijala'
! 
!  privremeno za stampu
! 
        IF (INDSC.EQ.0.OR.INDSC.EQ.3.OR.INDSC.EQ.4) THEN
           IF(IDJERDAP.GT.0) THEN
            CALL VOTVTEMPDJ(KKORAK)
            write(*,*)'RESTEL,NPRINT,KORAK',NPRINT,KKORAK
            CALL VRESTEL(VVREME,TT1,NPT,ISNUMER,KKORAK)
           ELSE
            CALL VOTVTEMP(KKORAK)
            write(*,*)'RESTEL,NPRINT,KORAK',NPRINT,KKORAK
            CALL VRESTEL(VVREME,TT1,NPT,ISNUMER,KKORAK)
          ENDIF
         CALL VSTAU09c(TT1,NPT,47,1,1,KKORAK,NZAD,NUMZAD,KONT,TT1,
     1               ISNUMER)
C        CALL VPRKONT1c(TT1,ID,UBRZ,CORD,NZAD,IOSA,IVECT,ZADVRE,KOJK)
C stampanje dubine vode
C       CALL VSTAU09c(TT1,CORD,ID,NPT,49,11,1,KKORAK,NZAD,NUMZAD,KONT,TT1)
C stampa rezultata u neu za sve cvorove STAS09T
C stampanje elemenata na kojima je zadat flux
!         write(*,*)'stampanje fluksa/zracenja'
        IF(MAXSIL.GT.0) 
     1  CALL VSTAS09T(NPT,49,10,1,KKORAK,NGPSIL,MAXSIL,ISNUMER)
!  stampanje elemenata na kojima je zadata prelaznost - stampa se ID f-je koef. prelaznosti
        IF(MAXTQE.GT.0)THEN 
          CALL VSTAS09T(NPT,49,12,1,KKORAK,NGPSIL,MAXTQE,ISNUMER)
!  stampanje temperature prema Bofangu
!        WRITE(*,*)'stampanje temperature prema Bofangu'        
          CALL VSTAS09T(NPT,49,19,1,KKORAK,NGPSIL,MAXTQE,ISNUMER)
        ENDIF
! 
!        WRITE(*,*)'stampanje zadatog potencijala'        
       CALL VSTAU09T(TT1,NPT,49,9,1,KKORAK,NZAD,NUMZAD,KONT,
     1       TT1,ISNUMER)
!        WRITE(*,*)'stampanje temperatura'        
         IF(INDJOS.EQ.1) THEN
        CALL VSTAU09T(TTemp,NPT,49,1,1,KKORAK,NZAD,NUMZAD,KONT,
     1       TT1,ISNUMER)
        ELSE
        CALL VSTAU09T(TT1,NPT,49,1,1,KKORAK,NZAD,NUMZAD,KONT,
     1       TT1,ISNUMER)
        ENDIF
!  DUBINA VODE
!        CALL VSTAU09T(TT1,CORD,ID,NPT,49,11,1,KKORAK,NZAD,NUMZAD,KONT,TT1)
!  BRZINE
!       CALL VSTAU09T(VECTJ,CORD,ID,NPT,49,41,3,KKORAK,NZAD,NUMZAD,KONT,
!      1     TT1)
!        WRITE(*,*)'stampanje gradijenta'        
       CALL VSTAU09T(GRADJN,NPT,49,51,3,KKORAK,NZAD,NUMZAD,KONT,
     1            TT1,ISNUMER)
C stampanje rezultata u svim cvorovima
!       WRITE(*,*)'pre lst'        
!       WRITE(*,*)'stampanje IZLL3DT'        
      CALL VIZLL3DT(TT1,KKORAK,VECTJ,VVREME,VG,GG,
     1            KOTES,ISNUMER)
!       WRITE(*,*)'pre STAG3D'        
!       WRITE(*,*)'stampanje UNV'        
      CALL VSTAG3D(TT1,NASLOV,VVREME,KKORAK,1,NPT,18,0,NET,
     1            KORAK,VECTJ,GRADJN,NZAD,FZAPR,ISNUMER)
      ENDIF
      IF (INDSC.EQ.1.OR.INDSC.EQ.3) THEN
        CALL VSTAU09c(TT1,NPT,47,1,1,KKORAK,NZAD,NUMZAD,KONT,TT1,
     1               NEL,ISNUMER)
      ENDIF
      IF(INDSC.EQ.2.OR.INDSC.EQ.4) THEN
         WRITE(*,*)'stampanje STAU09CT'        
        CALL VSTAU09CT(TT1,47,VVREME,KKORAK,ISNUMER)
         WRITE(*,*)'posle stampanje STAU09CT'        
!        IF(NPRINT.EQ.1)THEN
!         CALL VRESTEL(VVREME,TT1,NPT,IDJSTAMP1,NP3D1)
!       WRITE(*,*)'stampanje STAU09T'        
!         CALL VSTAU09T(TT1,NPT,49,1,1,KKORAK,NZAD,NUMZAD,KONT,
!      1       T1,ISNUMER)
!        ELSEIF(KORAK.EQ.NPRINT) THEN
!       write(*,*)'TEMP,NPRINT,KORAK',NPRINT,KORAK
!        CALL VSTAU09T(TT1,NPT,49,1,1,KKORAK,NZAD,NUMZAD,KONT,
!      1       TT1,ISNUMER)
!          CALL VRESTEL(VVREME,TT1,NPT,IDJSTAMP1,NP3D1)
!        ENDIF
      ENDIF
      ENDIF
C
! Stampanje rezultata u presecima
      IF(ISPRESEK.GT.0) THEN
         CALL VSTAMPPRESEKNEU(TT1,VECTJ,
     +       VREME,VVREME,KKORAK,ISNUMER)
      ENDIF
! 
      IF(NXNASTRAN.EQ.1) RETURN
C
      KKORAK=KKORAK+1
c uslova da sledeci korak u tom periodu ne postoji i da ide na sledeci period
 2345 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      VRNPKOR=DABS(VREME(NNPER,KORAK+1))
      CALL MPI_BCAST(VRNPKOR,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      IF (VRNPKOR.LT.1.D-10) THEN
        GOTO 600
      ELSE
cz ide na sledeci korak   
          IF(mcm_kojpak.eq.7)then
            GOTO 600
             ENDIF
        GOTO 35
      ENDIF
! topalovic ovde ide kraj prtlje po koracima 
  600 CONTINUE
         CALL DATE_AND_TIME(VALUES=Dtime)
      WRITE(*,*) 'vreme posle stampanja', (Dtime(i),i=5,7)
      WRITE(3,*) 'vreme posle stampanja', (Dtime(i),i=5,7)
      if (myid.eq.0) then
        IF((NSTAC.EQ.1.AND.KKORAK.EQ.2).AND.(mcm_kojpak.NE.7)) 
     1    CALL VIZBACIT(TT1,CORD,ID,NPT,NEL,NET,NDIM,IIZLAZ)
      endif
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      RETURN      
      END
C=======================================================================
C=======================================================================
      SUBROUTINE VIZLL3DT(TT1,KORAK,VECTJ,VVREME,VG,GG,
     1                  KOTES,ISNUMBER)
      USE NODES
      USE ELEMENTS
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
      COMMON /ICITANJE/INPT
      COMMON /PRIKAZ/ INDSC,IZPOT

      DIMENSION D4(4),TT1(*),VECTJ(3,*),N4(4)
      DIMENSION NRED(8),NREDK(8),NRED1(3)
      DIMENSION VG(3,NET,*),GG(3,NET,*),NREDT(4)
      DIMENSION NREDT10(10),VGT(3),NRED26(6),NRED28(8)
      DATA NRED/8,4,2,6,7,3,1,5/  
      DATA NREDK/27,9,3,21,25,7,1,19/
      DATA NRED1/2,1,3/
      DATA NRED26/1,2,3,12,23,31/
      DATA NRED28/1,2,3,4,12,23,34,41/
      DATA NREDT/1,2,3,4/
      DATA NREDT10/1,2,3,4,12,23,13,14,24,34/
      
C
      if (indsc.eq.0) then
C      
CS  TEMPERATURA U CVOROVIMA 
CE  NODE TEMPERATURE 
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2001)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6001)
      WRITE (IIZLAZ,*)'STEP NUMBER  ',KORAK

      KG = 0
   20 DO 50 I=1,4
      KG = KG +1
      IF(KG.GT.NPT) GO TO 100
      IF(ISNUMBER.EQ.0) THEN
        N4(I) = KG 
      ELSE
        N4(I) = NCVEL(KG)
      ENDIF
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
      IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5012) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5002) (N4(J1),D4(J1),J1=1,I)
      ENDIF
      GO TO 77
C
   88 IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,5013) (N4(J1),D4(J1),J1=1,I)
      ELSE
       WRITE(IIZLAZ,5003) (N4(J1),D4(J1),J1=1,I)
      ENDIF
   77 IF(KG.LT.NPT) THEN
      GO TO 20
      ENDIF
C
CS  GRADIJENTI
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2069)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6069)
      IF (INDSC.EQ.0)THEN
      DO I=1,NET
         IF(elemtip(I).gt.100) THEN
          NTYPE=elemtip(I)/100
          NTIMES=100
        ELSE
          NTYPE=elemtip(I)/10
          NTIMES=10
        ENDIF
        NDIMEL=elemtip(I)-NTIMES*NTYPE
        IF(ISNUMBER.EQ.0) THEN
          WRITE(IIZLAZ,5000) I
        ELSE
          WRITE(IIZLAZ,5000) MCVEL(I)
        ENDIF
        DO J=1,NDIMEL
         IF(ISNUMBER.EQ.0) THEN
           NC=NEL(J,I)
         ELSE
           NC=NCVEL(NEL(J,I))
         ENDIF
         NG=NRED(J)
         SELECT CASE(elemtip(I))
          CASE(12:13)
            NG=NRED1(J)
          CASE(23:26)
            NG=NRED26(J)
          CASE(34,310)
            NG=NREDT10(J)
          CASE DEFAULT
!          WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NBREL)
         END SELECT
         IF(NDIMEL.eq.2) NG=NRED1(J)
         IF(NDIMEL.eq.10) NG=NREDT10(J)
	   IF(KOTES.EQ.1) NG=NREDK(J)
	   IF(INPT.EQ.1) THEN
                IF(NG.GT.10) THEN
                  NG1=NG/10
                  NG2=NG-NG1*10
                  DO K=1,3
                    VGT(K)=(GG(K,I,NG1)+(GG(K,I,NG2)))/2
                  ENDDO
                  WRITE(IIZLAZ,5017) NC,(VGT(K),K=1,3)
                 ELSE
                  WRITE(IIZLAZ,5017) NC,(GG(K,I,NG),K=1,3)
	         ENDIF
	   ELSE
	     IF(NG.GT.10) THEN
                NG1=NG/10
                NG2=NG-NG1*10
                DO K=1,3
                  VGT(K)=(GG(K,I,NG1)+(GG(K,I,NG2)))/2
                ENDDO
                WRITE(IIZLAZ,5007) NC,(VGT(K),K=1,3)
             ELSE
                WRITE(IIZLAZ,5007) NC,(GG(K,I,NG),K=1,3)
             ENDIF
	   ENDIF
        ENDDO
      ENDDO
      ENDIF
c
      endif      
C
      RETURN
 5000 FORMAT(/' ELEMENT',I10)     
 5002 FORMAT(4(I5,1PD13.5))
 5003 FORMAT(4(I5,3X,F10.3))
 5005 FORMAT(F10.3,10(1PE13.5))
 5007 FORMAT(I5,3(1PE13.5))
 5012 FORMAT(4(I10,1PD13.5))
 5013 FORMAT(4(I10,3X,F10.3))
 5017 FORMAT(I10,3(1PE13.5))
 2001 FORMAT(/'    T E M P E R A T U R E               '
     1,'             '//
     1' CVOR  POMERANJE   CVOR  POMERANJE   CVOR  POMERANJE   CVOR POMER 
     1NJE   '/' BROJ',3(14X,'BROJ'))
 2069 FORMAT(//'     G R A D I J E N T I  U   C V O R U   '/
     1' CVOR GRADIJENT GX GRADIJENT GY GRADIJENT GZ')
 6001 FORMAT(/'    N O D A L    T O T A L     T E M P E R A T U R E    '
     1,'                 '//
     1' NODE  DISPLAC.    NODE  DISPLAC.    NODE  DISPLAC.    NODE DISPL
     1AC.   '/'  No.',3(15X,'No.'))
 6069 FORMAT(//'    N O D A L    G R A D I E N T S  '/
     1' NODE  GRADIENT GX  GRADIENT GY  GRADIENT GZ')
      END
C=======================================================================
      SUBROUTINE VDJERDAPREAD(IDJSTAMP1,NP3D1)
       USE NODES
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       COMMON /DJERDAP/ IDJERDAP,ISPRESEK
       COMMON /ICITANJE/ INPT
       COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
       DIMENSION NP3D1(5),NDJEL(21),IDJSTAMP1(5,*)
       CHARACTER*6 IME(5)
       CHARACTER*2 IMECSV1
       CHARACTER*1 IMECSV
       
!        CALL VICLEAR(NP3D1,5)
!        KK=1
        IF (IDJERDAP.LE.10) THEN
          write(IMECSV1,20) IDJERDAP
          IME='DJ3D' // IMECSV1    
!         if(IPAKT.EQ.1)IME='PDJT' // IMECSV1 // '.NEU'
        ELSE
          write(IMECSV,21) IDJERDAP
          IME='DJ3D' // IMECSV     
!         if(IPAKT.EQ.1)IME='PDJT' // IMECSV // '.NEU'
        ENDIF
       IIFILE=90
! !        IMECSV=IME(1:IB-1)//'.CSV'
!        Write( xpresek, '(i10)' ) IDJERDAP
!        IME(1)='DJ3D'//xpresek
!        IME(2)='DJ3D2'
!        IME(3)='DJ3D3'
!        IME(4)='DJ3D4'
!        IME(5)='DJ3D5'
       icitanjedj=1
       DO KK=1,icitanjedj
          NP3D1(KK)=0
          IIFILE=IIFILE+IDJERDAP
c          OPEN (IIFILE,FILE='DJ3D1',STATUS='UNKNOWN',FORM='FORMATTED',
          OPEN (IIFILE,FILE=IME(KK),STATUS='UNKNOWN',FORM='FORMATTED',
     1      ACCESS='SEQUENTIAL')
C
!            OPEN (92,FILE='DJ3D2',STATUS='UNKNOWN',FORM='FORMATTED',
!      1      ACCESS='SEQUENTIAL')
!            OPEN (93,FILE='DJ3D3',STATUS='UNKNOWN',FORM='FORMATTED',
!      1      ACCESS='SEQUENTIAL')
!            OPEN (94,FILE='DJ3D4',STATUS='UNKNOWN',FORM='FORMATTED',
!      1      ACCESS='SEQUENTIAL')
!            OPEN (95,FILE='DJ3D5',STATUS='UNKNOWN',FORM='FORMATTED',
!      1      ACCESS='SEQUENTIAL')

           READ(IIFILE,4111,err=9999) NDJ3D1,NCVE
!           IF(STATTT.EQ.'OLD') THEN
!             OPEN (91,FILE='DJ3D1',STATUS='OLD',FORM='FORMATTED',
!      1      ACCESS='SEQUENTIAL')
!             READ(91,4111) NDJ3D1,IT1,NCVE
!             write(*,*) 'NDJ3D1,NCVE',NDJ3D1,NCVE
           DO I=1,NDJ3D1
!                 READ(IIFILE,*)
                IF(INPT.EQ.1) THEN
                   READ(IIFILE,4111) (NDJEL(J),J=1,NCVE+1)
!                    READ(IIFILE,4111) (NDJEL(J),J=IT1+1,NCVE)
                ELSE
                   READ(IIFILE,1111) (NDJEL(J),J=1,NCVE+1)
!                    READ(IIFILE,1111) (NDJEL(J),J=IT1+1,NCVE)
                ENDIF
                DO K=2,NCVE+1
                 IF(ISNUMER.EQ.0)THEN
                   IDJSTAMP1(KK,NDJEL(K))=KK
                 ELSE
                   NN=NDJEL(K)
                   N=NN-NPI+1
                   NI=NELCV(N)
                   IDJSTAMP1(KK,NI)=KK
                 ENDIF
                ENDDO
             ENDDO
             DO I=1,NPT
! odredjivanje broja cvorova za koje se pisu zapreminske sile             
                IF(IDJSTAMP1(KK,I).EQ.KK) NP3D1(KK)=NP3D1(KK)+1
             ENDDO
!              write(*,*) 'kraj citanja NP3D1',NP3D1(KK)
!           ELSE
!              STOP ' NEMA FAJLA DJ3D1'
!           ENDIF
        ENDDO
        go to 9998  
 9999     write(*,*) ' NE POSTOJI FAJL DJ3D '         
          write(3,*) ' NE POSTOJI FAJL DJ3D '
   20  FORMAT (I1)
   21  FORMAT (I2)
 1111 FORMAT(21I5)
 4111 FORMAT(21I10)
 9998 RETURN
      END
C=======================================================================
      SUBROUTINE VDJERDAPREADT(NP3D1,ISNUMER,II)
       USE PPR
       USE NODES
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       COMMON /DJERDAP/ IDJERDAP,ISPRESEK
       COMMON /ICITANJE/ INPT
       COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
       DIMENSION NDJEL(21),NP3D1(18)
!        CHARACTER*6 IME(5)
        CHARACTER*20 IME
        CHARACTER*2 IMECSV
        CHARACTER*1 IMECSV1
!         write(*,*)"IP",IP
        IIDJERDAP=IDJERDAP
        IF(IDJERDAP.EQ.-1) THEN
         IIDJERDAP=II
        if (allocated(IDJSTAMP1)) deallocate(IDJSTAMP1)
        ENDIF
        IF (IIDJERDAP.LT.10) THEN
          write(IMECSV1,20) IIDJERDAP
          IME='DJ3D' // IMECSV1    
!         if(IPAKT.EQ.1)IME='PDJT' // IMECSV1 // '.NEU'
        ELSE
          write(IMECSV,21) IIDJERDAP
          IME='DJ3D' // IMECSV     
!         if(IPAKT.EQ.1)IME='PDJT' // IMECSV // '.NEU'
        ENDIF
        write(*,*) "ime", IME
!        KK=1
       IIFILE=90
        KK=1
!        DO KK=1,1
          NP3D1(KK)=0
          IIFILE=IIFILE+IIDJERDAP
c          OPEN (IIFILE,FILE='DJ3D1',STATUS='UNKNOWN',FORM='FORMATTED',
          OPEN (IIFILE,FILE=IME,STATUS='UNKNOWN',FORM='FORMATTED',
     1      ACCESS='SEQUENTIAL')
C
!  6.8.2018. promenjeno da cita cvorove u kojima se stampaju sile, ne elemente sa cvorovima za koje se stampaju
!           READ(IIFILE,4111,err=9999) NDJ3D1,IT1,NCVE
          READ(IIFILE,4001,err=9999)NP3D1(KK)
      if(.not.allocated(IDJSTAMP1))allocate
     1   (IDJSTAMP1(1,NP3D1(KK)),STAT=istat) 
      DO I=1,1
        DO JJ=1,NP3D1(KK)
          IDJSTAMP1(I,JJ)=0
        ENDDO
      ENDDO
!             OPEN (91,FILE='DJ3D1',STATUS='OLD',FORM='FORMATTED',
!      1      ACCESS='SEQUENTIAL')
!             READ(91,4111) NDJ3D1,IT1,NCVE
             write(*,*) 'NP3D1(KK)',NP3D1(KK)
           DO I=1,NP3D1(KK)
!                 READ(IIFILE,*)
                IF(INPT.EQ.1) THEN
                   READ(IIFILE,4001) IDJSTAMP1(KK,I)
!                    READ(IIFILE,4111) (NDJEL(J),J=1,IT1)
!                    READ(IIFILE,4111) (NDJEL(J),J=IT1+1,NCVE)
                ELSE
                   READ(IIFILE,1001) IDJSTAMP1(KK,I)
!                    READ(IIFILE,1111) (NDJEL(J),J=1,IT1)
!                    READ(IIFILE,1111) (NDJEL(J),J=IT1+1,NCVE)
                ENDIF
!                 DO K=1,NCVE
                  IF(ISNUMER.EQ.0)THEN
!                    IDJSTAMP1(KK,NDJEL(K))=KK
                  ELSE
!                     NN=NDJEL(K)
                    NN=IDJSTAMP1(KK,I)
                    N=NN-NPI+1
                    NI=NELCV(N)
                    IDJSTAMP1(KK,I)=NI
! !              write(*,*) "I,k,NN",I,K,NN,NI              
                  ENDIF
!              write(*,*) 'IDJSTAMP1',I, IDJSTAMP1(KK,I)
!                 ENDDO
            ENDDO
!              DO I=1,NPT
!                 IF(IDJSTAMP1(KK,I).EQ.KK) NP3D1(KK)=NP3D1(KK)+1
!              ENDDO
             write(*,*) 'kraj citanja NP3D1',NP3D1(KK)
!           ELSE
!              STOP ' NEMA FAJLA DJ3D1'
!           ENDIF
!         ENDDO
        go to 9998  
 9999     write(*,*) ' NE POSTOJI FAJL DJ3DT '         
          write(3,*) ' NE POSTOJI FAJL DJ3D '
   20  FORMAT (I1)
   21  FORMAT (I2)
 1001 FORMAT(I5)
 1111 FORMAT(14I5)
 4001 FORMAT(I10)
 4111 FORMAT(13I10)
 9998 RETURN
      END
C=======================================================================
       SUBROUTINE VREADPRESEK (ISNUMER)
        USE PRESEK
        USE ELEMENTS
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        COMMON /DJERDAP/ IDJERDAP,ISPRESEK
        COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
!         DIMENSION NP3D1(5),NDJEL(21)
        CHARACTER*250 ACOZ
        CHARACTER*20 IME
        CHARACTER*2 IMECSV
        CHARACTER*1 IMECSV1,PREF
!         CHARACTER*4 PRESECI(18,3)
!        DATA PRESECI/2,1,12/
       IPIJEZ=250
!         select case (IP)
!          case(1)
!           PREF='O'
!          case(2)
!           PREF='S'
!          case(3)
!           PREF='R'
!          case(4)
!           PREF='B'
!         end select       
        write(*,*) 'IDJERDAP',IDJERDAP
        IF (IDJERDAP.GE.10) THEN
          write(IMECSV,21) IDJERDAP
          IME='L' // IMECSV // '.PRE'    
!           IME='L' // IMECSV // PREF // '.DAT'    
        ELSEIF (IDJERDAP.GE.1) THEN
          write(IMECSV1,20) IDJERDAP
          IME='L' // IMECSV1 // '.PRE'    
!           IME='L0' // IMECSV1 // PREF // '.DAT'    
        ENDIF
         write(*,*) "IME ",IME
!         IPIJEZ=IPIJEZ
        OPEN (IPIJEZ,FILE=IME,STATUS='UNKNOWN',FORM='FORMATTED',
     1      ACCESS='SEQUENTIAL')
!         IF(IPAKT.EQ.1)THEN
!        OPEN (IPIJEZ,FILE='PIJEZT.DAT',STATUS='UNKNOWN',
!      1      FORM='FORMATTED',ACCESS='SEQUENTIAL')
!         ENDIF
!          IF(IP.GT.1)Then
!            deallocate (NP_ID)
!            deallocate (NP_ELEMENT)
!            deallocate (NP_COORDS)
!            deallocate (NELP)
!          ENDIF
!      
         READ(IPIJEZ,355) IPRES
         if (.not.allocated(NPRESEK)) 
     1              allocate(NPRESEK(IPRES),STAT=istat)
         if (.not.allocated(NPELEM)) 
     1              allocate(NPELEM(IPRES),STAT=istat)
         READ(IPIJEZ,356) (NPRESEK(I),I=1,IPRES)
         READ(IPIJEZ,356) (NPELEM(I),I=1,IPRES)
         MAXNP=NPRESEK(1)
         MAXELP=NPELEM(1)
         IF(IPRES.GT.1) THEN
           do INP=2,IPRES
             IF(MAXNP.LT.NPRESEK(INP)) MAXNP=NPRESEK(INP)
             IF(MAXELP.LT.NPELEM(INP)) MAXELP=NPELEM(INP)
           enddo
         ENDIF
!           write(*,*)"IPRES,MAXNP,MAXELP",IPRES,MAXNP,MAXELP
!            if (NPRESEK.eq.0) STOP 'nisu definisani cvorovi u preseku'
              if (.not.allocated(NP_ID)) 
     1              allocate(NP_ID(MAXNP,IPRES),STAT=istat)
              if (.not.allocated(NP_ELEMENT)) allocate(NP_ELEMENT
     1                                 (MAXNP,IPRES),STAT=istat)
              if (.not.allocated(NP_COORDS)) allocate(NP_COORDS
     1                               (6,MAXNP,IPRES),STAT=istat)
              if (.not.allocated(NPROP)) allocate(NPROP
     1                                 (MAXNP,IPRES),STAT=istat)
              if (.not.allocated(NELP)) allocate(NELP
     1                             (5,MAXELP,IPRES),STAT=istat)
         do INP=1,IPRES
!            write(*,*)"NPRESEK(INP),NPELEM",NPRESEK(INP),NPELEM(INP)
           DO I=1,NPRESEK(INP)
            IF(IPAKT.eq.0) THEN
             READ(IPIJEZ,360) NP_ID(I,INP),NP_ELEMENT(I,INP),
     1             (NP_COORDS(jj,I,INP),jj=1,6),NPROP(I,INP)
            ELSEIF(IPAKT.eq.1) THEN
             READ(IPIJEZ,358) NP_ID(I,INP),NP_ELEMENT(I,INP),
     1             (NP_COORDS(jj,I,INP),jj=1,6),NPROP(I,INP)
            ENDIF
!        write(3,*)"NP_ID(I)",NP_ID(I,INP),NP_ELEMENT(I,INP),NPROP(I,INP)
            IF(ISNUMER.EQ.1) THEN
               NNMP=NP_ELEMENT(I,INP)
               JJ=NNMP-NMI+1
            IF(NNMP.LT.NMI.OR.NNMP.GT.NMA.OR.MELCV(JJ).EQ.0) THEN
            write(*,*) 'pre', INP,NP_ID(I,INP),I,NNMP,'van modela'
                STOP 
                ENDIF
                  NP_ELEMENT(I,INP)=MELCV(JJ)
              ENDIF
            ENDDO
! C elementi u kojima se stampaju vrednosti za prikaz rezultata
!             if (NPELEM.eq.0) STOP 'nisu definisani elementi u preseku'
            DO II=1,NPELEM(INP)
               READ(IPIJEZ,1000) ACOZ
               READ(ACOZ,359) (NELP(jj,II,INP),jj=1,4)
               READ(ACOZ,361) NELP(5,II,INP)
!                READ(IPIJEZ,359) (NELP(jj,II,INP),jj=1,4)
!                READ(IPIJEZ,361) NELP(5,II,INP)
!           write(3,*)"NELP read",II,(NELP(ji,II,INP),ji=1,5)
! !             write(3,*)"dpoint",DPOINT_ID(I),DP_ELEMENT(I)
            ENDDO
          enddo
!            endif
   20  format (I1)
   21  format (I2)
  355  FORMAT(I10)
  356  FORMAT(23I10)
  358  FORMAT(2I10,6F10.2,I10)
  360  FORMAT(2I10,6F10.2,10X,I10)
  359  FORMAT(4I10)
  361  FORMAT(210X,I5)
 1000  FORMAT(A250)
 1117  FORMAT(21I10,2I5,E10.3,I5)
       RETURN
      END
C=======================================================================
       SUBROUTINE VREADPRESEKALL (ISNUMER)
        USE PRESEK
        USE ELEMENTS
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        COMMON /DJERDAP/ IDJERDAP,ISPRESEK
        COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
!         DIMENSION NP3D1(5),NDJEL(21)
        CHARACTER*250 ACOZ
        CHARACTER*20 IME
        CHARACTER*2 IMECSV
        CHARACTER*1 IMECSV1,PREF
!         CHARACTER*4 PRESECI(18,3)
!        DATA PRESECI/2,1,12/
       IPIJEZ=250
        IF(ISPRESEK.EQ.-1) IME='L-ALL.PRE'
        IF(ISPRESEK.EQ.-2) IME='V-ALL.PRE'
         write(*,*) "IME ",IME
        OPEN (IPIJEZ,FILE=IME,STATUS='UNKNOWN',FORM='FORMATTED',
     1      ACCESS='SEQUENTIAL')
         READ(IPIJEZ,355) IPRES
         IIPRES=IPRES
!          IPRES=IIPRES
         if (.not.allocated(NPRESEK)) 
     1              allocate(NPRESEK(IIPRES),STAT=istat)
         if (.not.allocated(NPELEM)) 
     1              allocate(NPELEM(IIPRES),STAT=istat)
!          if(IPRES.GT.24)THEN
!           READ(IPIJEZ,356) (NPRESEK(I),I=1,24)
!           READ(IPIJEZ,356) (NPELEM(I),I=1,24)
!           READ(IPIJEZ,356) (NPRESEK(I),I=25,IIPRES)
!           READ(IPIJEZ,356) (NPELEM(I),I=25,IIPRES)
!          ELSE
          READ(IPIJEZ,1010) (NPRESEK(I),I=1,IIPRES)
          READ(IPIJEZ,1010) (NPELEM(I),I=1,IIPRES)
!          ENDIF
         MAXNP=NPRESEK(1)
         MAXELP=NPELEM(1)
         IF(IIPRES.GT.1) THEN
           do INP=2,IIPRES
             IF(MAXNP.LT.NPRESEK(INP)) MAXNP=NPRESEK(INP)
             IF(MAXELP.LT.NPELEM(INP)) MAXELP=NPELEM(INP)
           enddo
         ENDIF
!           write(*,*)"IPRES,MAXNP,MAXELP",IIPRES,MAXNP,MAXELP
!            if (NPRESEK.eq.0) STOP 'nisu definisani cvorovi u preseku'
              if (.not.allocated(NP_ID)) 
     1              allocate(NP_ID(MAXNP,IIPRES),STAT=istat)
              if (.not.allocated(NP_ELEMENT)) allocate(NP_ELEMENT
     1                                 (MAXNP,IIPRES),STAT=istat)
              if (.not.allocated(NP_COORDS)) allocate(NP_COORDS
     1                               (6,MAXNP,IIPRES),STAT=istat)
              if (.not.allocated(NPROP)) allocate(NPROP
     1                                 (MAXNP,IIPRES),STAT=istat)
              if (.not.allocated(NELP)) allocate(NELP
     1                             (5,MAXELP,IIPRES),STAT=istat)
         do INP=1,IIPRES
!             write(*,*)"IPRES",INP
!             write(*,*)"NPRESEK(INP),NPELEM",NPRESEK(INP),NPELEM(INP)
           DO I=1,NPRESEK(INP)
            IF(IPAKT.eq.0) THEN
             READ(IPIJEZ,360) NP_ID(I,INP),NP_ELEMENT(I,INP),
     1             (NP_COORDS(jj,I,INP),jj=1,6),NPROP(I,INP)
            ELSEIF(IPAKT.eq.1) THEN
             READ(IPIJEZ,358) NP_ID(I,INP),NP_ELEMENT(I,INP),
     1             (NP_COORDS(jj,I,INP),jj=1,6),NPROP(I,INP)
            ENDIF
!        write(3,*)"NP_ID(I)",NP_ID(I,INP),NP_ELEMENT(I,INP),NPROP(I,INP)
            IF(ISNUMER.EQ.1) THEN
               NNMP=NP_ELEMENT(I,INP)
               JJ=NNMP-NMI+1
            IF(NNMP.LT.NMI.OR.NNMP.GT.NMA.OR.MELCV(JJ).EQ.0) THEN
            write(*,*) 'pre', INP,NP_ID(I,INP),I,NNMP,'van modela'
                STOP 
                ENDIF
                  NP_ELEMENT(I,INP)=MELCV(JJ)
              ENDIF
            ENDDO
! C elementi u kojima se stampaju vrednosti za prikaz rezultata
!             if (NPELEM.eq.0) STOP 'nisu definisani elementi u preseku'
            DO II=1,NPELEM(INP)
               READ(IPIJEZ,1000) ACOZ
               READ(ACOZ,359) (NELP(jj,II,INP),jj=1,4)
               READ(ACOZ,361) NELP(5,II,INP)
!               write(*,*)"NPELEM",NPELEM(INP),II
!                READ(IPIJEZ,359) (NELP(jj,II,INP),jj=1,4)
!                READ(IPIJEZ,361) NELP(5,II,INP)
            ENDDO
          enddo
!            endif
   20  format (I1)
   21  format (I2)
  355  FORMAT(I10)
  356  FORMAT(24I10)
  358  FORMAT(2I10,6F10.2,I10)
  360  FORMAT(2I10,6F10.2,10X,I10)
  359  FORMAT(4I10)
  361  FORMAT(210X,I5)
 1000  FORMAT(A250)
 1010  FORMAT(10I10)
 1117  FORMAT(21I10,2I5,E10.3,I5)
       RETURN
      END
C=======================================================================
