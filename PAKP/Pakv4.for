C==========================================================================
      SUBROUTINE PERIOD(VREME)
      USE MESURMENTPOINTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*250 ACOZ
      COMMON /VREPER/ NPER,NTABFT
      COMMON /VDP/ DT,NKORP,NN,NZAV
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /SRPSKI/ ISRPS
      COMMON /PRIKAZ/ INDSC,IZPOT
      DIMENSION VREME(NPER,950)

      DO I=1,NPER
      DO J=1,950
       VREME(I,J)=0.
      ENDDO    
      ENDDO
C
      NKOR=0
C
      DO I=1,NPER
      CALL ISPITA(ACOZ)
      READ(ACOZ,1004) NKORP,DT
 1004 FORMAT (I5,D10.3)
       IF(INDSC.EQ.2.OR.INDSC.EQ.1) go to 20
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,4002) NKORP,DT
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,5002) NKORP,DT
 4002 FORMAT(//
     111X,'BROJ KORAKA PO PERIODU....................... NKORP =',I5/
     111X,'VREDNOST KORAKA.............................. DT  =',1PD10.3)
 5002 FORMAT(//
     111X,'NUMBER OF STEPS IN PERIOD ................... NKORP =',I5/
     111X,'TIME STEP IN PERIOD ......................... DT  =',1PD10.3)
 20      CONTINUE
      DO J=1,NKORP
        VREME(I,J)=DT
      ENDDO
C
      NKOR=NKOR+NKORP
C
      ENDDO
      BRKORAKA=NKOR
      RETURN
      END
C==========================================================================
C==========================================================================
      SUBROUTINE ULTABF (TABF,ITFMAX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /VREPER/ NPER,NTABFT
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /SRPSKI/ ISRPS
      COMMON /PRIKAZ/ INDSC,IZPOT

      CHARACTER*250 ACOZ

      DIMENSION TABF(2,NTABFT,*),ITFMAX(*)

      DO I=1,NTABFT
       CALL ISPITA(ACOZ)
       READ(ACOZ,1000) IBR,IMAX
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,4000) IBR,IMAX
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,5000) IBR,IMAX
       ITFMAX(IBR)=IMAX
       DO J=1,IMAX
       CALL ISPITA(ACOZ)
       READ(ACOZ,1004) TABF(1,IBR,J),TABF(2,IBR,J)
      IF(INDSC.EQ.2.OR.INDSC.EQ.1) go to 20
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,4002) TABF(1,IBR,J),TABF(2,IBR,J)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,5002) TABF(1,IBR,J),TABF(2,IBR,J)
 20      continue
       ENDDO
      ENDDO
      RETURN
 1004 FORMAT (2D10.3)
 1000 FORMAT (2I5)
 4002 FORMAT(//
     111X,'VREME (ARGUMENT)................................t =',1PD10.3/
     111X,'VREDNOST VREMENSKE FUNKCIJE..................f(t) =',1PD10.3)
 4000 FORMAT(//
     111X,'REDNI BROJ VREMENSKE FUNKCIJE (IBR) .................. =',I5/
     111X,'MAKSIM. BR. TACAKA ZA DATU VREM. FUNK. (IMAX) ........ =',I5)
 5002 FORMAT(//
     111X,'TIME .......................................... t =',1PD10.3/
     111X,'TIME FUNCTION VALUE ........................ f(t) =',1PD10.3)
 5000 FORMAT(//
     111X,'TIME FUNCTION NUMBER ........................... (IBR) =',I5/
     111X,'MAKSIMAL NUMBER OF POINTS FOR TIME FUNCTION ... (IMAX) =',I5)
      END
C==========================================================================
C==========================================================================
      SUBROUTINE TIMFUN(TABF,SILAF,UVREME,IBRT,IFUN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /VREPER/ NPER,NTABFT
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT

      DIMENSION TABF(2,NTABFT,*)

      IF(IFUN.EQ.0) THEN
c         stop 'stop - (IFUN.EQ.0)'        
        RETURN
      ENDIF
!       WRITE(3,*) 'IFUN',IFUN
      
      TOL=1.D-20
      IF (DABS(UVREME).LT.TOL) UVREME=1.D-20
       VREML=UVREME+UVREME/1.D10
       VREMD=UVREME-UVREME/1.D10
C
      IF(IBRT.EQ.1) THEN
        SILAF=TABF(2,IFUN,1)
        RETURN
      ENDIF

      DO I=1,IBRT-1
      IF((VREML.GT.TABF(1,IFUN,I)).AND.(TABF(1,IFUN,I+1).GT.VREMD)) THEN
        SRT=TABF(1,IFUN,I+1)-TABF(1,IFUN,I) 
        SRF=TABF(2,IFUN,I+1)-TABF(2,IFUN,I) 
        DELTAT=UVREME-TABF(1,IFUN,I)
        SILAF=TABF(2,IFUN,I)+(DELTAT/SRT)*SRF
        RETURN
       ENDIF
      ENDDO
      WRITE(IIZLAZ,*)'GRESKA U ULAZNIM PODACIMA ZA VREM. FUNKCIJU',IFUN
      WRITE(IIZLAZ,*)'VREME JE IZVAN OPSEGA VREMENSKE FUNKCIJE', UVREME
      WRITE(*,*)'VREME JE IZVAN OPSEGA VREMENSKE FUNKCIJE',IFUN, UVREME
      STOP
      END
C==========================================================================
C==========================================================================
      SUBROUTINE MVECT(IVECT)
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /DODAT/ NDIMM

      DIMENSION IVECT(*)   
!       DIMENSION NEL(NDIMM,*)    

      DO I=1,NPT
       IVECT(I)=0
      ENDDO      


      DO I=1,NET
      IF(elemtip(I).gt.100) THEN
        NTYPE=elemtip(I)/100
        NTIMES=100
      ELSE
        NTYPE=elemtip(I)/10
        NTIMES=10
      ENDIF
      NDIMEL=elemtip(I)-NTIMES*NTYPE
       DO J=1,NDIMEL
        IF (NEL(J,I).NE.0) THEN
         IVECT(NEL(J,I))=IVECT(NEL(J,I))+1
        ENDIF
       ENDDO
      ENDDO      

      END
C======================================================================

C=======================================================================
      SUBROUTINE JEDNA1(A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO EQUALIZING 2 REAL VECTORS IN ACCORDANCE WITH TERM :
CS.    P R O G R A M
CS        ZA IZJEDNACAVANJE 2 REALNA VEKTORA U SAGLASNOSTI SA IZRAZOM :
C .
C .         A(I)=B(I)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),B(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' JEDNA1'
      DO 10 I=1,N
   10 A(I)=B(I)
      RETURN
      END
C=======================================================================

C=======================================================================
      SUBROUTINE DETER3(GRAD,DETG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CALCULATE DETERMINANT OF MATRIX 3 X 3
CS.   P R O G R A M
CS.      ZA RACUNANJE DETERMINANTE MATRICE 3 X 3
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION GRAD(3,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' DETER3'
      DETG=GRAD(1,1)*(GRAD(2,2)*GRAD(3,3)-GRAD(2,3)*GRAD(3,2))+
     1     GRAD(1,2)*(GRAD(2,3)*GRAD(3,1)-GRAD(2,1)*GRAD(3,3))+
     2     GRAD(1,3)*(GRAD(2,1)*GRAD(3,2)-GRAD(2,2)*GRAD(3,1))
      RETURN
      END
C=======================================================================

C=======================================================================
      SUBROUTINE MINV3(XJ,DET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INVERSE MATRIX 3 X 3
CS.   P R O G R A M
CS.      ZA INVERTOVANJE MATRICE 3 X 3
C .
C ......................................................................
C
      COMMON /ELEMEN/ ELAST(6,6),XJJ(3,3),ALFA(3),TEMP0,DETT,NLM,KK
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /SRPSKI/ ISRPS
      DIMENSION XJ1(3,3),XJ(3,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' MINV3 '
C
      CALL DETER3(XJ,DET)
      IF(DET.GT.1.D-10) GOTO 10
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2000) DET,NLM,NGE
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6000) DET,NLM,NGE
      PRINT *, ' DET < 0'
      IF(DABS(DET).LT.1.D-10) STOP 'MINV3'
   10 XJ1(1,1)= (XJ(2,2)*XJ(3,3)-XJ(3,2)*XJ(2,3))/DET
      XJ1(2,1)=-(XJ(2,1)*XJ(3,3)-XJ(3,1)*XJ(2,3))/DET
      XJ1(3,1)= (XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))/DET
      XJ1(1,2)=-(XJ(1,2)*XJ(3,3)-XJ(3,2)*XJ(1,3))/DET
      XJ1(2,2)= (XJ(1,1)*XJ(3,3)-XJ(3,1)*XJ(1,3))/DET
      XJ1(3,2)=-(XJ(1,1)*XJ(3,2)-XJ(3,1)*XJ(1,2))/DET
      XJ1(1,3)= (XJ(1,2)*XJ(2,3)-XJ(2,2)*XJ(1,3))/DET
      XJ1(2,3)=-(XJ(1,1)*XJ(2,3)-XJ(2,1)*XJ(1,3))/DET
      XJ1(3,3)= (XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2))/DET
      CALL JEDNA1(XJ,XJ1,9)
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(/
     1' DETERMINANTA MANJA ILI JEDNAKA NULI, DET = ',1PD10.3/
     1' ELEMENT =',I8,'     GRUPA =',I5)
C-----------------------------------------------------------------------
 6000 FORMAT(/
     1' ZERO OR NEGATIVE DETERMINANTE, DET = ',1PD10.3/
     1' ELEMENT =',I8,'     GROUP =',I5)
C-----------------------------------------------------------------------
      END
C===================================================================

C==========================================================================
      SUBROUTINE MINV(A,N,D,L,M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    INVERTOVANJE MATRICE
CE    INVERT OF MATRIX
CS    UPOTREBA
CS    CALL MINV(A,N,D,L,M)
CS    A=ULAZNA MATRICA KOJA SE ZAMENJUJE INVERTOVANOM
CS    N=RED MATRICE A
CS    D=DETERMINANTA MATRICE A
CS    L=RADNI VEKTOR DUZINE N
CS    M=RADNI VEKTOR DUZINE N
CS    UPOTREBLJEN JE STANDARDNI GAUSS-
CS    JORDAN-OV  ALGORITAM.RACUNA SE
CS    DETERMINANTA.DETERMINANTA NULA
CS    OZNACAVA DA JE MATRICA SINGULARNA.
C
CE    INVERZION OF MATRIX
C
      DIMENSION L(*),M(*)
C
      DIMENSION A(*)
C
CS    NACI NAJVECI ELEMENT
CE    SEARCHING MAX. MEMBER
C
      D=1.0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
C   10 IF(DABS(BIGA)-DABS(A(IJ))) 15,20,20
      IF(DABS(BIGA)-DABS(A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
C
CS    IZMENITI VRSTE
CE    CHANGE ROW
C
      J=L(K)
      IF(J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI)=HOLD
C
CS    IZMENITI KOLONE
CE    CHANGE COLUMN
C
   35 I=M(K)
      IF(I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI)=HOLD
C
CS    PODELITI KOLONU SA NEGATIVNIM STOZEROM-
CS    VREDNOST STOZERNOG ELEMENTA JE U BIGA.
C
   45 IF(BIGA) 48,46,48
   46 D=0.0
      RETURN
   48 DO 55 I=1,N
      IF(I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
C
CS    REDUKOVATI MATRICU
CE    REDUCTION OF MATRIX
C
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
   65 CONTINUE
C
C     PODELITI VRSTU STOZEROM
C
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
C
C     PROIZVOD STOZERA
C
      D=D*BIGA
C
C     ZAMENITI STOZER RECIPROCNOM VREDNOSCU
C
      A(KK)=1.0/BIGA
   80 CONTINUE
C
CS    POSLEDNJA IZMENA VRSTA I KOLONA
CE    LAST       
C
      K=N
  100 K=(K-1)
      IF(K) 150,150,105
  105 I=L(K)
      IF(I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI)=HOLD
  120 J=M(K)
      IF(J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI)=HOLD
      GO TO 100
  150 RETURN
      END
C==========================================================================
C==========================================================================
      SUBROUTINE PREB(TT21,TT1,NBREL)
      USE NODES
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /DODAT/ NDIMM
      COMMON /NDESUK/ NDES,IDPRIT,IFORM
      COMMON /NELDIM/ NDIMEL

!       DIMENSION ID(1,*),NEL(NDIMM,*),TT21(*),TT1(*)
      DIMENSION TT21(*),TT1(*)

      DO KLM=1,NDIMEL
         JED=ID(1,NEL(KLM,NBREL))
         IF(JED.NE.0) TT21(KLM)=TT1(JED)
      ENDDO
      END
C==========================================================================
C==========================================================================
      SUBROUTINE SSTAM(TT1,CORD,ID,NPROM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      
      DIMENSION TT1(*),CORD(3,*),ID(1,*)

      TOL=1.D-8


      DO I=1,NPT
       IF (DABS(CORD(1,I)-2.25D-3).LT.TOL) THEN
        IF(ID(NPROM,I).EQ.0) THEN
          S=0.D0
        ELSE
          S=TT1(ID(NPROM,I))
         ENDIF
       WRITE(IIZLAZ,1000)I,CORD(2,I),S
       ENDIF
      ENDDO      

1000  FORMAT(I5,2D13.5)

      END
C==========================================================================
C==========================================================================
      SUBROUTINE OPTIM1(TT1,SILE,NEL,ID,NZAD,ZADVRE,NGPSIL,
     &MAXA,CORD,SKEF,SKEFN,KONT,FZAPR,VREME,TABF,TT10,UBRZ,UBRZ0,BRZ,
     &BRZ0,AK,VECTJ,IVECT,POVSIL,GRADJN,ITFMAX,AKONST,NASLOV,
     &ICUR,PRIV,PRIV1,VG,GG,KOJK)
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /NDESUK/ NDES,IDPRIT,IFORM
      COMMON /VREPER/ NPER,NTABFT
      COMMON /PROMEN/ NJUTN,INDOPT,INDPRO
      COMMON /DODAT/ NDIMM
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /VDP/ DT,NKORP,NN,NZAV
      COMMON /CERNE/ IGRESK,IPROTK,IZSILE
C
      CHARACTER*80 NASLOV
      DIMENSION TT1(*),SILE(*),ZADVRE(*),CORD(3,*)
      DIMENSION TT10(*),UBRZ(*),UBRZ0(*),BRZ(*),BRZ0(*)
      DIMENSION NEL(NDIMM,*),ID(1,*),NZAD(3,*),NGPSIL(12,*),MAXA(*)
      DIMENSION SKEF(NDES,*),AKONST(3,5,*),SKEFN(*)
      DIMENSION POVSIL(2,*),ITFMAX(*)
      DIMENSION KONT(9,MAXLIN,*),KOJK(*)   
      DIMENSION VG(3,NET,*),GG(3,NET,*)
C
      DIMENSION FZAPR(2,*),VECTJ(2,*),GRADJN(2,*)
      DIMENSION IVECT(*),ICUR(*)
      DIMENSION VREME(NPER,950),AK(NDES,*)
      DIMENSION TABF(2,NTABFT,*),PRIV1(*),PRIV(*)
      DIMENSION IPOM(200)
      INTEGER*8 NGPSIL
C==========================================================================
C
	OPEN(IGRESK,FILE='GRESKA')
C
	DO 1001 L=1,NN+1
C	
      DO M=1,JEDN
       UBRZ0(M)=0.D0
       PRIV(M)=0.D0
      ENDDO      
C
        M=NJUTN
        XMIN=AKONST(1,INDPRO,M)
        XMAX=AKONST(2,INDPRO,M)
        ODNOS=AKONST(3,4,M)
C
C LINEARNA VEZA KMIN I KMAX
C
 	IF (NZAV.EQ.1) THEN
	DK=(XMAX-XMIN)/NN
	DDK=DK*L
	XSR=XMIN+DDK-DK
C
	ELSE
C
C EKSPONENCIJALNA VEZA KMIN I KMAX
C
	DK=(DLOG10(XMAX)-DLOG10(XMIN))/NN
	DDK=DK*L
	XSR=10.D0**(DLOG10(XMIN)+DDK-DK)
 	ENDIF
C
        TOL=1.D-2
       INDOPT=0
       ITER1=0
C
      DO 300 NNPER=1,NKORP
C 
       VVREME=DT*NNPER
       FIZAD=0.D0
       KBR=0      
       DO I=1,NUMZAD
        IF(NZAD(2,I).EQ.-2) THEN
        CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,I)),NZAD(3,I))
         FIZAD=FIZAD+ZADVRE(I)*FK1
         KBR=KBR+1
         IPOM(KBR)=ID(1,NZAD(1,I))
        ENDIF
       ENDDO
C
       IF (KBR.EQ.0) THEN
        WRITE(*,*)'NIJE ZADATO U KOM CVORU SE TRAZI POREDJENJE!!!'
        STOP
       ENDIF
C
       FIZAD=FIZAD/KBR
C
         AKONST(3,INDPRO,M)=XSR
C
 	IF (ODNOS.GT.1.D-10.AND.INDPRO.EQ.1) 
     1  AKONST(3,2,M)=AKONST(3,INDPRO,M)/ODNOS
 	IF (ODNOS.GT.1.D-10.AND.INDPRO.EQ.2) 
     1  AKONST(3,1,M)=AKONST(3,INDPRO,M)*ODNOS
C
      CALL RACUN(TT1,SILE,NEL,ID,NZAD,ZADVRE,NGPSIL,
     &MAXA,CORD,SKEF,SKEFN,KONT,FZAPR,VREME,TABF,TT10,UBRZ,UBRZ0,
     &AK,VECTJ,IVECT,POVSIL,GRADJN,ITFMAX,AKONST,NASLOV,
     &ICUR,VG,GG,KOJK)
         SRE=0.D0
         DO I=1,KBR
         SRE=SRE+TT1(IPOM(I))
         ENDDO
         FSR=SRE/KBR
C
       BRZ(NNPER)=FIZAD
       UBRZ0(NNPER)=FSR
C
	IF (NNPER.GT.1) THEN
	AA=UBRZ0(NNPER-1)
	BB=UBRZ0(NNPER)
	CC=BRZ(NNPER-1)
	DD=BRZ(NNPER)
	R=DABS(((AA+BB)/2.D0*DT)-((CC+DD)/2.D0*DT))
C
        PRIV(NNPER)=PRIV(NNPER-1)+R
	ENDIF
C
 300   CONTINUE
C
	GRESKA=PRIV(NKORP)
C
       IF(L.EQ.1.AND.NZAV.EQ.0) WRITE(IGRESK,*) '1   LOG'
       IF(L.EQ.1.AND.NZAV.EQ.1) WRITE(IGRESK,*) '1   LIN'
       IF(L.EQ.1) WRITE(IGRESK,*) 'PERMEABILITY     ERROR'
       WRITE(IGRESK,2003) XSR,GRESKA
 2003  FORMAT(1E12.5,1E12.5)
C
	BRZ0(L)=XSR
	PRIV1(L)=GRESKA
 1001 CONTINUE
C
	DO J=1,NN+1
	IF(PRIV1(J).LT.GRESKA) THEN
	GRESKA=PRIV1(J)
	XSR=BRZ0(J)
	ENDIF
	ENDDO
C
       WRITE(IIZLAZ,5000) XSR
       INDOPT=1
       AKONST(3,INDPRO,M)=XSR
C
 5000 FORMAT(//78('*')/,
     1'...............OPTIMAL PERMEABILITY K=',1PD12.5,'...............
     1...............'
     1/,78('*')/)
C
	CLOSE(IGRESK)
C
        END
C==========================================================================
C==========================================================================
      SUBROUTINE PRKONT1(TT1,ID,UBRZ,CORD,NZAD,IOSA,IVECT,ZADVRE,KOJK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /PRIKAZ/ INDSC,IZPOT
      DIMENSION TT1(*),ID(1,*),UBRZ(*),CORD(3,*),NZAD(3,*),IVECT(*)
      DIMENSION ZADVRE(*),KOJK(*)


         TULAZ=0.D0
         TIZLAZ=0.D0
         PROT=0.D0
         FSIZ1=0.D0
         FSIZ2=0.D0
         FSIS1=0.D0
         FSIS2=0.D0
         SPOLJ1=0.D0
         SPOLJ2=0.D0
        ULPROT=0.D0
        ZLPROT=0.D0
        AKORIT=0.D0
        PUL=0.
        PIZ=0.
          WRITE(IIZLAZ,*)'PROTOCI PO CVOROVIMA RAZLICITI OD NULE:'

C=====================================================
          DO II=1,NUMZAD          
           JJEDN=ID(1,NZAD(1,II))
           IF (JJEDN.NE.0) THEN
             Q=UBRZ(JJEDN)
             IF (DABS(ZADVRE(II)-(320.)).LT.1.E-5) AKORIT=AKORIT+Q
             IF (Q.GT.0.D0) THEN
               ULPROT=ULPROT+Q
             ELSE
               ZLPROT=ZLPROT+Q
             ENDIF
             IF (DABS(Q).GT.1.D-10) THEN
               if(indsc.eq.0) WRITE(IIZLAZ,1000)NZAD(1,II),Q
             ENDIF  
           ENDIF
          ENDDO
        WRITE(IIZLAZ,*)'PROTOK PO KONTURI KORITA= ',AKORIT
        WRITE(IIZLAZ,*)'ULAZNI PROTOK ZA ZADATE CVOROVE= ',ULPROT
        WRITE(IIZLAZ,*)'IZLAZNI PROTOK ZA ZADATE CVOROVE= ',ZLPROT
        WRITE(IIZLAZ,*)
        WRITE(IIZLAZ,*)

	  DO 100 I=1,NPT
           JJEDN=ID(1,I)
          IF (JJEDN.NE.0) THEN
           Q=UBRZ(JJEDN)
           
           IF (KOJK(I).GT.0) THEN
            IF (Q.GT.0.D0) THEN
              PUL=PUL+Q
            ELSE
              PIZ=PIZ+Q
            ENDIF
             IF (DABS(Q).GT.1.D-10) THEN
               if(indsc.eq.0) WRITE(IIZLAZ,1000)I,Q
             ENDIF  
           ENDIF
C OVO PROVERITI           
c           IF (IVECT(I).LE.2) THEN
           IF (IVECT(I).LE.NDIM-1) THEN
             IF (DABS(Q).GT.1.D-10) THEN
               if(indsc.eq.0) WRITE(IIZLAZ,*)I,IVECT(I),Q
             ENDIF  
            IF (Q.GT.0.D0) THEN
              SPOLJ1=SPOLJ1+Q
            ELSE
              SPOLJ2=SPOLJ2+Q
            ENDIF
           ENDIF
           IF (Q.GT.0.D0) THEN
             TULAZ=TULAZ+Q
           ELSE
             TIZLAZ=TIZLAZ+Q
           ENDIF
           PROT=PROT+Q
           
C=====================================================
C SPECIJALNI DODATNI USLOV ZA CERNI:
           IF((TT1(JJEDN)-CORD(IOSA,I)).LT.1.D-10) THEN
            IF (Q.GT.0.D0) THEN
              FSIZ1=FSIZ1+Q
            ELSE
              FSIZ2=FSIZ2+Q
            ENDIF

           ELSE

            IF (Q.GT.0.D0) THEN
              FSIS1=FSIS1+Q
            ELSE
              FSIS2=FSIS2+Q
            ENDIF

           ENDIF
          ENDIF
 100      CONTINUE
C

        WRITE(IIZLAZ,*)'UKUPNI PROTOK ZA SVE CVOROVE= ',PROT
        WRITE(IIZLAZ,*)'UKUPNI ULAZNI PROTOK ZA SVE CVOROVE= ',TULAZ
        WRITE(IIZLAZ,*)'UKUPNI IZLAZNI PROTOK ZA SVE CVOROVE= ',TIZLAZ
        WRITE(IIZLAZ,*)'UKUPNI ULAZNI PROTOK IZNAD SL. POVRS.= ',FSIZ1
        WRITE(IIZLAZ,*)'UKUPNI IZLAZNI PROTOK IZNAD SL. POVRS.= ',FSIZ2
        WRITE(IIZLAZ,*)'UKUPNI ULAZNI PROTOK  ISPOD SL. POVRS.= ',FSIS1
        WRITE(IIZLAZ,*)'UKUPNI IZLAZNI PROTOK ISPOD SL. POVRS.= ',FSIS2
        WRITE(IIZLAZ,*)
        WRITE(IIZLAZ,*)'ULAZNI PROTOK ZA SPOLJASNJU KONTURU= ',SPOLJ1
        WRITE(IIZLAZ,*)'IZLAZNI PROTOK ZA SPOLJASNJU KONTURU= ',SPOLJ2
        WRITE(IIZLAZ,*)
        WRITE(IIZLAZ,*)'ULAZNI PROTOK ZA KONTURNE CVOROVE= ',PUL
        WRITE(IIZLAZ,*)'IZLAZNI PROTOK ZA KONTURNE CVOROVE= ',PIZ
         
1000    FORMAT('CVOR=',I5,5X,'FLUX=',D13.5)          

       END
C==========================================================================
C==========================================================================
      SUBROUTINE PRKONT1c(TT1,ID,UBRZ,CORD,NZAD,IOSA,IVECT,ZADVRE,KOJK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      DIMENSION TT1(*),ID(1,*),UBRZ(*),CORD(3,*),NZAD(3,*),IVECT(*)
      DIMENSION ZADVRE(*),KOJK(*)


         TULAZ=0.D0
         TIZLAZ=0.D0
         PROT=0.D0
         FSIZ1=0.D0
         FSIZ2=0.D0
         FSIS1=0.D0
         FSIS2=0.D0
         SPOLJ1=0.D0
         SPOLJ2=0.D0
        ULPROT=0.D0
        ZLPROT=0.D0
        AKORIT=0.D0
        PUL=0.
        PIZ=0.
C         WRITE(IIZLAZ,*)'PROTOCI PO CVOROVIMA RAZLICITI OD NULE:'

C=====================================================
          DO II=1,NUMZAD          
           JJEDN=ID(1,NZAD(1,II))
           IF (JJEDN.NE.0) THEN
             Q=UBRZ(JJEDN)
             IF (DABS(ZADVRE(II)-(320.)).LT.1.E-5) AKORIT=AKORIT+Q
             IF (Q.GT.0.D0) THEN
               ULPROT=ULPROT+Q
             ELSE
               ZLPROT=ZLPROT+Q
             ENDIF
C             IF (DABS(Q).GT.1.D-10) THEN
C               WRITE(IIZLAZ,1000)NZAD(1,II),Q
C             ENDIF  
           ENDIF
          ENDDO
C        WRITE(IIZLAZ,*)'PROTOK PO KONTURI KORITA= ',AKORIT
        WRITE(47,*) ULPROT
        WRITE(47,*) -1
C        WRITE(49,*)'ULAZNI PROTOK ZA ZADATE CVOROVE= ',ULPROT
C        WRITE(49,*)'IZLAZNI PROTOK ZA ZADATE CVOROVE= ',ZLPROT
C        WRITE(IIZLAZ,*)
C        WRITE(IIZLAZ,*)

C	  DO 100 I=1,NPT
C           JJEDN=ID(1,I)
C          IF (JJEDN.NE.0) THEN
C           Q=UBRZ(JJEDN)
C           
C           IF (KOJK(I).GT.0) THEN
C            IF (Q.GT.0.D0) THEN
C              PUL=PUL+Q
C            ELSE
C              PIZ=PIZ+Q
C            ENDIF
C             IF (DABS(Q).GT.1.D-10) THEN
C               WRITE(IIZLAZ,1000)I,Q
C             ENDIF  
C           ENDIF
C OVO PROVERITI           
Cc           IF (IVECT(I).LE.2) THEN
C           IF (IVECT(I).LE.NDIM-1) THEN
C             IF (DABS(Q).GT.1.D-10) THEN
C               WRITE(IIZLAZ,*)I,IVECT(I),Q
C             ENDIF  
C            IF (Q.GT.0.D0) THEN
C              SPOLJ1=SPOLJ1+Q
C            ELSE
C              SPOLJ2=SPOLJ2+Q
C            ENDIF
C           ENDIF
C           IF (Q.GT.0.D0) THEN
C             TULAZ=TULAZ+Q
C           ELSE
C             TIZLAZ=TIZLAZ+Q
C           ENDIF
C           PROT=PROT+Q
C           
CC=====================================================
CC SPECIJALNI DODATNI USLOV ZA CERNI:
C           IF((TT1(JJEDN)-CORD(IOSA,I)).LT.1.D-10) THEN
C            IF (Q.GT.0.D0) THEN
C              FSIZ1=FSIZ1+Q
C            ELSE
C              FSIZ2=FSIZ2+Q
C            ENDIF
C
C           ELSE
C
C            IF (Q.GT.0.D0) THEN
C              FSIS1=FSIS1+Q
C            ELSE
C              FSIS2=FSIS2+Q
C            ENDIF
C
C           ENDIF
C          ENDIF
C 100      CONTINUE
C

C        WRITE(IIZLAZ,*)'UKUPNI PROTOK ZA SVE CVOROVE= ',PROT
C        WRITE(IIZLAZ,*)'UKUPNI ULAZNI PROTOK ZA SVE CVOROVE= ',TULAZ
C        WRITE(IIZLAZ,*)'UKUPNI IZLAZNI PROTOK ZA SVE CVOROVE= ',TIZLAZ
C        WRITE(IIZLAZ,*)'UKUPNI ULAZNI PROTOK IZNAD SL. POVRS.= ',FSIZ1
C        WRITE(IIZLAZ,*)'UKUPNI IZLAZNI PROTOK IZNAD SL. POVRS.= ',FSIZ2
C        WRITE(IIZLAZ,*)'UKUPNI ULAZNI PROTOK  ISPOD SL. POVRS.= ',FSIS1
C        WRITE(IIZLAZ,*)'UKUPNI IZLAZNI PROTOK ISPOD SL. POVRS.= ',FSIS2
C        WRITE(IIZLAZ,*)
C        WRITE(IIZLAZ,*)'ULAZNI PROTOK ZA SPOLJASNJU KONTURU= ',SPOLJ1
C        WRITE(IIZLAZ,*)'IZLAZNI PROTOK ZA SPOLJASNJU KONTURU= ',SPOLJ2
C        WRITE(IIZLAZ,*)
C        WRITE(IIZLAZ,*)'ULAZNI PROTOK ZA KONTURNE CVOROVE= ',PUL
C        WRITE(IIZLAZ,*)'IZLAZNI PROTOK ZA KONTURNE CVOROVE= ',PIZ
C         
1000    FORMAT('CVOR=',I5,5X,'FLUX=',D13.5)          

       END
C==========================================================================
      SUBROUTINE CURI1(CORD,ALEVO,SILE,MAXA,TT1,ID,ICUR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /DODAT/ NDIMM
      COMMON /KRITER/ IDOKL,MAXCUR
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR

      DIMENSION CORD(3,*),ALEVO(*),SILE(*),TT1(*)
      DIMENSION MAXA(*),ID(1,*),ICUR(*)

	IF (MAXCUR.EQ.0) RETURN

         DO 10 I=1,IDOKL
            NODE=ICUR(I)
            JJ=ID(1,NODE)
            IF(JJ.LE.0) GO TO 10
            H=CORD(IOSA,NODE)
            DH=H-TT1(JJ)
         WRITE(3,*) 'CVOR, COORD, FI', NODE,H,TT1(JJ)
            IF(DH.LT.1.D-6) THEN
              SILE(JJ)=1.0D35*DH
              ALEVO(MAXA(JJ))=1.0D35
            ENDIF
   10    CONTINUE
      RETURN
      END
C==========================================================================
C==========================================================================
      SUBROUTINE CKONV1(ICUR,TT1,CORD,ID,KONVER,UBRZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /KRITER/ IDOKL,MAXCUR
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR

      DIMENSION TT1(*),CORD(3,*),UBRZ(*)
      DIMENSION ICUR(*),ID(1,*)

	IF (MAXCUR.EQ.0) RETURN

        MAX=IDOKL
        KONVER=1
        IF (IDOKL.EQ.MAXCUR) RETURN

        PRPLUS=0.D0
        DO I=1,JEDN
         IF (UBRZ(I).GT.0.D0) PRPLUS=PRPLUS+UBRZ(I)
        ENDDO
       
        TOL=PRPLUS*0.01D0

         NODE=ICUR(IDOKL+1)
         JJ=ID(1,NODE)
          IF ( TT1(JJ).GT.CORD(IOSA,NODE) ) THEN
            WRITE(IIZLAZ,*)'NODE=',NODE
            MAX=IDOKL+1
            KONVER=0
         ENDIF

       IDOKL=MAX
      RETURN
      END
C==========================================================================
C==========================================================================
      SUBROUTINE INIDOK(NZAD,ICUR,CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /KRITER/ IDOKL,MAXCUR
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
      COMMON /SRPSKI/ ISRPS

      DIMENSION NZAD(3,*),ICUR(*)
      DIMENSION CORD(3,*)

C POTPROGRAM ZA DEFINISANJE CVOROVA NA LINIJI CURENJA U NIZU ICUR(*)
C ODREDJIVANJE MAKSIMALNOG BROJA CVOROVA NA LINIJI CURENJA -MAXCUR
C I INICIJALIZACIJA ZA POCETAK PRORACUNA SA UZETIM PRVIM CVOROM 
C KAO POCETNIM ZA LINIJU CURENJA

        IDOKL=0
C
        K=0     
        DO I=1,NUMZAD
         IF (NZAD(2,I).EQ.4) THEN
            NODE=NZAD(1,I)
              K=K+1
              ICUR(K)=NODE
         ENDIF
        ENDDO   

        MAXCUR=K

        DO I=1,K-1
         DO J=I+1,K
          NODE1=ICUR(I)
          NODE2=ICUR(J)
          IF (CORD(IOSA,NODE1).GT.CORD(IOSA,NODE2)) THEN
            IP=ICUR(I)
            ICUR(I)=ICUR(J)
            ICUR(J)=IP
          ENDIF
         ENDDO
        ENDDO

C          IDOKL=0
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,*) 'CVOROVI NA LINIJI CURENJA POREDJANI PO VISINI'
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,*) 'NODES ON SEEPAGE LINE '
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,*) 'CVOR         VERTIKALNA(Y/Z)-KOORDINATA'
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,*) 'NODE         VERTICAL(Y/Z)-COORDINATE'
         DO I=1,MAXCUR
          WRITE(IIZLAZ,1000)ICUR(I),4,1,CORD(IOSA,ICUR(I))
         ENDDO
C       CALL IWRR(ICUR,MAXCUR,'ICUR=')
 1000   FORMAT (3I6,F10.3)
      
      END
C==========================================================================
C=========================================================================
      SUBROUTINE PROCEN(RS2,PR,FAK,NEL,NBREL,CORD) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /DODAT/ NDIMM

      DIMENSION RS2(*),PR(*),L(4)
      DIMENSION NEL(NDIMM,*)
      DIMENSION CORD(3,*)

        KBRPOT=0
        DO I=1,4
           L(I)=1
          IF (PR(I).LT.0.D0) L(I)=-1
          RS2(I)=0.D0
           IF (L(I).EQ.1) KBRPOT=KBRPOT+1
        ENDDO
       
         PROC=0.25D0*KBRPOT
         XL=1.D10
         XD=-1.D10
         DO I=1,4
          NODE=NEL(I,NBREL)
          IF (CORD(1,NODE).LT.XL) XL=CORD(1,NODE)
          IF (CORD(1,NODE).GT.XD) XD=CORD(1,NODE)
         ENDDO

C	  PROC=8.D0
          POVPAD=(XD-XL)
         
          
      IF ( (L(1).EQ.-1).AND.(L(2).EQ.-1)
     &.AND.(L(3).EQ.-1).AND.(L(4).EQ.-1)) THEN
       RETURN
      ELSE IF ( (L(1).EQ.1).AND.(L(2).EQ.1)
     &.AND.(L(3).EQ.1).AND.(L(4).EQ.1)) THEN 
       RETURN
      ELSE
        DO I=1,4
         IF (L(I).EQ.1) RS2(I)=FAK*POVPAD/KBRPOT
        ENDDO
      ENDIF
       
      END
C==========================================================================
C=========================================================================
      SUBROUTINE PSKEFN(SKEF,SKEFN,NDIMEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /NDESUK/ NDES,IDPRIT,IFORM
       DIMENSION SKEF(NDES,*),SKEFN(*)

	 K=0
       DO I=1,NDIMEL
        DO J=I,NDIMEL
          K=K+1
         	SKEFN(K)=SKEF(I,J)
        ENDDO
       ENDDO

      END
C==========================================================================

C======================================================================
      SUBROUTINE DELJIV(IBROJ,IDELI,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO TEST DIVIDE OF NUMBER WITHOUT RESIDUE
CS.    P R O G R A M
CS.        ZA ISPITIVANJE DELJIVOSTI BROJA BEZ OSTATKA
C .
CE.            IBROJ  - BROJ KOJI DELIMO
CS.            IBROJ  - BROJ KOJI DELIMO
C .
CE.            IDELI  - BROJ KOJIM DELIMO
CS.            IDELI  - BROJ KOJIM DELIMO
C .
CE.            IND    - INDIKATOR OF DIVIDE
CE.                     =0; WITHOUT RESIDUE
CE.                     =0; WITH RESIDUE
CS.            IND    - INDIKATOR DELJIVOSTI
CS.                     =0; BEZ OSTATKA
CS.                     =0; SA OSTATKOM
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' DELJIV'
      BROJ=IBROJ
      DELI=IDELI
      REZ=BROJ/DELI
      IREZ=REZ
      OST=REZ-IREZ
      TOL=0.0D00
      IND=0
      IF(DABS(OST).GT.TOL) IND=1
      RETURN
      END
C=======================================================================

C======================================================================
      SUBROUTINE CHECKJ(NELG,NBREL,CORDG,CORD,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /SRPSKI/ ISRPS
      COMMON /DODAT/ NDIMM

      DIMENSION P(3,21),XJ(3,3),CORDG(3,*),CORD(3,*)
      DIMENSION NELG(NDIMM,*)
      DIMENSION CK(3,21),H(21)
C
      R=0.D0
      S=0.D0
      T=0.D0

  10  DO KLM=1,NDIM
       NODE=NELG(KLM,NBREL)
       IF(NODE.GT.NPT) THEN
          CK(1,KLM)=CORDG(1,NODE-NPT)
          CK(2,KLM)=CORDG(2,NODE-NPT)
          CK(3,KLM)=CORDG(3,NODE-NPT)
       ELSE
          CK(1,KLM)=CORD(1,NODE)
          CK(2,KLM)=CORD(2,NODE)
          CK(3,KLM)=CORD(3,NODE)
       ENDIF
      ENDDO
C
      RP=1.0+R
      SP=1.0+S
      TP=1.0+T
      RM=1.0-R
      TM=1.0-T
      SM=1.0-S
      RR=1.0-R*R
      SS=1.0-S*S
      TT=1.0-T*T
      DO 82 I=1,21
      H(I)=0.
      DO 82 J=1,3
   82 P(J,I)=0.
C
C     INTERPOLACIJSKE FUNKCIJE I NJIHOVI IZVODI
CE    INTERPOLATION FUNCTIONS AND DERIVATIVES
C
CS    PRVIH 8 CVOROVA
CE    FIRST 8 NODES
C
      H(1)=0.125*RP*SP*TP
      H(2)=0.125*RM*SP*TP
      H(3)=0.125*RM*SM*TP
      H(4)=0.125*RP*SM*TP
      H(5)=0.125*RP*SP*TM
      H(6)=0.125*RM*SP*TM
      H(7)=0.125*RM*SM*TM
      H(8)=0.125*RP*SM*TM
C
      P(1,1)=0.125*SP*TP
      P(1,2)=-P(1,1)
      P(1,3)=-0.125*SM*TP
      P(1,4)=-P(1,3)
      P(1,5)=0.125*SP*TM
      P(1,6)=-P(1,5)
      P(1,7)=-0.125*SM*TM
      P(1,8)=-P(1,7)
C
      P(2,1)=0.125*RP*TP
      P(2,2)=0.125*RM*TP
      P(2,3)=-P(2,2)
      P(2,4)=-P(2,1)
      P(2,5)=0.125*RP*TM
      P(2,6)=0.125*RM*TM
      P(2,7)=-P(2,6)
      P(2,8)=-P(2,5)
C
      P(3,1)=0.125*RP*SP
      P(3,2)=0.125*RM*SP
      P(3,3)=0.125*RM*SM
      P(3,4)=0.125*RP*SM
      P(3,5)=-P(3,1)
      P(3,6)=-P(3,2)
      P(3,7)=-P(3,3)
      P(3,8)=-P(3,4)
C
C      IF (NDIM.EQ.8) GO TO 50
C
CS    JAKOBIJAN U TACKI R,S,T
CE    JACOBIAN AT POINT R,S,T
C
      DO 61 I=1,3
      DO 61 J=1,3
      XJ(I,J)=0.
      DO 61 KK=1,NDIM
      XJ(I,J)=XJ(I,J)+P(I,KK)*CK(J,KK)
   61 CONTINUE
C     WRITE(3,8888)((XJ(I,J),J=1,3),I=1,3)
C8888 FORMAT(' X',3D11.3)
C
CS     DETERMINANTA JAKOBIJANA U TACKI R,S,T
CE     JACOBIAN DETERMINANT AT POINT R,S,T
C
      CALL MINV3(XJ,DET1)
C
C     WRITE(3,7888) ((XJ(I,J),J=1,3),I=1,3)
C7888 FORMAT(' X-',3D11.3)
      IF (DET1.GT.1.D-15) GO TO 77
      IF (IND.EQ.0) THEN
C      WRITE(IIZLAZ,2000) NBREL,0,R,S,T,DET1
      IF (DET1.LT.0.D0) THEN
      DO I=1,4
       IPRIV=NELG(I,NBREL)
       NELG(I,NBREL)=NELG(I+4,NBREL)
       NELG(I+4,NBREL)=IPRIV
      ENDDO
      ELSE
       WRITE(IIZLAZ,200) (NELG(J,NBREL),J=1,NDIM)
       STOP
      ENDIF
      ELSE
      WRITE(IIZLAZ,2000) NBREL+NET,0,R,S,T,DET1
       STOP
      ENDIF
      GOTO 10
C

 200  FORMAT(8(I10))
C
C
 77   IF( IND.EQ.1) WRITE(IIZLAZ,201)NBREL+NET
      RETURN
C

  201  FORMAT('DET ZA ELEMENT ',I10,' JE OK')
C
 2000 FORMAT(' ** GRESKA **: JAKOBIJAN JEDNAK ILI MANJI OD NULE',
     1       ' ZA ELEMENT No.',I10/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',1PD12.5)
C 6000 FORMAT(' ** ERROR **: JACOBIAN EQUAL OR LESS THEN ZERO',
C     1       ' FOR ELEMENT No.',I10/
C     1       9X,'KFIX=',I5/
C     2       12X,'R=',F10.5/
C     3       12X,'S=',F10.5/
C     4       12X,'T=',F10.5/
C     5       10X,'DET=',F10.5)
C

      END
C=======================================================================
      SUBROUTINE ULTABFPT(NTAV,FN,NBR,MAXT,IBR,IMAX,ITMAX,TFMX,TMAX,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*80 ACOZ
CE READ INPUT DATA FOR FUNCTIONS GIVEN AS TABLE
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /SRPSKI/ ISRPS
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      DIMENSION NTAV(*),FN(2,NBR,*)
      DO 100 N=1,NBR
         CALL ISPITA(ACOZ)
         IF(INDFOR.EQ.1) READ(1,*) IBR,IMAX
         IF(INDFOR.EQ.2) READ(ACOZ,1000) IBR,IMAX
         IF(IMAX.LE.MAXT) GO TO 20
         RETURN
   20    NTAV(IBR)=IMAX
         IF(ISRPS.EQ.0) WRITE(IIZLAZ,2000) IBR,NTAV(IBR)
         IF(ISRPS.EQ.1) WRITE(IIZLAZ,6000) IBR,NTAV(IBR)
 2000 FORMAT(//11X,'FUNKCIJA BR.=',I3,5X,'BROJ TACAKA NA KRIVOJ=',I3/)
 6000 FORMAT(//11X,'FUNCTION No.=',I3,5X,'NUMBER POINTS ON CURVE=',I3/)
         DO 25 J=1,IMAX
            CALL ISPITA(ACOZ)
            IF(INDFOR.EQ.1) READ(1,*)(FN(I,IBR,J),I=1,2)
            IF(INDFOR.EQ.2) READ(ACOZ,1001)(FN(I,IBR,J),I=1,2)
  25     CONTINUE
         IF(IND.EQ.0) GO TO 30
         IF(FN(1,IBR,IMAX).GE.TMAX) GO TO 30
         ITMAX=1
         TFMX=FN(1,IBR,IMAX)
         RETURN
   30    IF(ISRPS.EQ.0) WRITE(IIZLAZ,2001)
         IF(ISRPS.EQ.1) WRITE(IIZLAZ,6001)
 2001 FORMAT(11X,'ARGUMENT',6X,'FUNKCIJA')
 6001 FORMAT(11X,'ARGUMENT',6X,'FUNCTION')
      WRITE(IIZLAZ,5002)((FN(I,IBR,J),I=1,2),J=1,IMAX)
 5002 FORMAT(9X,1PD12.5,2X,1PD12.5)
  100 CONTINUE
      IBR=0
      RETURN
 1000 FORMAT(10I5)
 1001 FORMAT(2F10.3)
      END
C=======================================================================
      SUBROUTINE BOFANG(Y1,Y2,TW1,TW2,TM,TS,C,Y)
      USE PREDISCRIBED
!       INTEGER*8 
      DOUBLE PRECISION FTOK,Y1,Y2,TW1,TW2,TM,TS,C,Y
      
      C=(TW1*exp(-Alpha*Y2)-TW2*exp(-Alpha*Y1))/((1-exp(-Alpha*Y1))*
     1               exp(-Alpha*Y2)-(1-exp(-Alpha*Y2))*exp(-Alpha*Y1))
      TS=(TW2*(1-exp(-Alpha*Y1))-TW1*(1-exp(-Alpha*Y2)))/((1-exp
     1(-Alpha*Y1))*exp(-Alpha*Y2)-(1-exp(-Alpha*Y2))*exp(-Alpha*Y1))
      TM=C+(TS-C)*exp(-Alpha*Y)
      TM=Rd_BOFANG*TM

      RETURN
      END
C=======================================================================
      SUBROUTINE LININTGV(Y1,Y2,TW1,TW2,TM,IND,Y)
      USE PREDISCRIBED
      COMMON /INDBRANA/ IBRANA
      INTEGER*4 IND,IBRANA
      DOUBLE PRECISION FTOK,Y1,Y2,TW1,TW2,TM,TS,Y,CKOEF,DELILAC
      
!         WRITE(3,*) "pakv4, Ibrana",IBRANA
!  Djerdap
      IF(IBRANA.EQ.0) THEN
        IF(Y.LT.Y1) THEN
          W=Y/Y1
          TS=0.5*(TM+TW1)
          IF(TS.LE.0) TS=0
          TM=TS*(1-W)+TW1*W
        ELSEIF(Y.GE.Y1.AND.Y.LE.Y2) THEN
          W=(Y-Y1)/(Y2-Y1)
          TM=TW1*(1-W)+TW2*W
        ELSE
          TM=TW2
        ENDIF
        IF(IND.EQ.1) THEN
          TS=0.5*(TM+TW1)
          W=D_BOFANG/Y1
          TD=TS*(1-W)+TW1*W
          TM=0.5*(TS+TD)
!           WRITE(3,*) "D,Y1,TS,TW1,TM",D_BOFANG,Y1,TS,TW1,TM
        ENDIF  
!  Grancarevo
      ELSEIF(IBRANA.EQ.1) THEN
         IF(IBOFANG.EQ.1) THEN
           DELILAC=0.0
           DELILAC=(1-exp(-Alpha*Y1))*exp(-Alpha*Y2)
           DELILAC=DELILAC-(1-exp(-Alpha*Y2))*exp(-Alpha*Y1)
           CKOEF=(TW1*exp(-Alpha*Y2)-TW2*exp(-Alpha*Y1))/DELILAC
           TS=TW2*(1-exp(-Alpha*Y1))-TW1*(1-exp(-Alpha*Y2))
           TS=TS/DELILAC
           TM=CKOEF+(TS-CKOEF)*exp(-Alpha*Y)
        ELSEIF(IBOFANG.EQ.2) THEN
!  TW2 - temperatura na povrsini akumulacije        
!  TW1 - temperatura na dnu akumulacije        
           TM=TW1+(TW2-TW1)*exp(-Alpha*Y)
        ENDIF
!  Granc arevo kraj 
!         WRITE(3,*) "Alpha,DELILAC,CKOEF,TS",Alpha,DELILAC,CKOEF,TS
!         WRITE(3,*) "Y1,Y2,TW1,TW2,TM,IND,Y",Y1,Y2,TW1,TW2,TM,IND,Y
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE TEMP_DEMPENDANCE(TM,TS,C)
      USE PREDISCRIBED
!       INTEGER*8 
      DOUBLE PRECISION TM,TS,C
      
      TM=C+((TS-C)/(Alpha*D_BOFANG))*(1-exp(-Alpha*D_BOFANG))
            
      RETURN
      END
C=======================================================================
      SUBROUTINE PRED_INTERP(I,INTERPAXIS,NUMAXISPOINTS,INTAXISPOINT,
     1                     ITFMAX,NZAD,AXISPTCORD,CORD,TABF,VVREME,FK1)
      integer*4 NZAD(3,*),ITFMAX(*),I,INTERPAXIS
      integer*8 ILOC,INDLOC,NSEG,NUMAXISPOINTS,INTAXISPOINT(2,*)
                         
      double precision AXISPTCORD(*),CORD(3,*),VVREME,FK1,TABF,Fx1,Fx2
      
      DO ILOC=2,NUMAXISPOINTS
         IF(AXISPTCORD(1).GE.CORD(INTERPAXIS,NZAD(1,I))) THEN
            INDLOC=1
            EXIT
         ENDIF
         IF(AXISPTCORD(NUMAXISPOINTS).LE.CORD(INTERPAXIS,
     1                                      NZAD(1,I))) THEN
            INDLOC=3
            EXIT
         ENDIF
         IF(AXISPTCORD(ILOC-1).LE.CORD(INTERPAXIS,NZAD(1,I)).AND.
     1         AXISPTCORD(ILOC).GE.CORD(INTERPAXIS,NZAD(1,I))) THEN
            INDLOC=2
            NSEG=ILOC-1
            EXIT
         ENDIF
        ENDDO
!         WRITE(3,*) "INDLOC",INDLOC
        SELECT CASE (INDLOC)
        CASE (1)
!          WRITE(3,*) "INTAXISPOINT(2,1)",INTAXISPOINT(2,1)
         CALL TIMFUN(TABF,FK1,VVREME,ITFMAX(INTAXISPOINT(2,1))
     1                                         ,INTAXISPOINT(2,1))
!          WRITE(3,*) "FK1",FK1
        CASE (2)
!          WRITE(3,*) "INTAXISPOINT(2,",NSEG,")",INTAXISPOINT(2,NSEG)
!          WRITE(3,*) "INTAXISPOINT(2,",NSEG+1,")",INTAXISPOINT(2,NSEG+1)
         CALL TIMFUN(TABF,Fx1,VVREME,ITFMAX(INTAXISPOINT(2,NSEG))
     1                                      ,INTAXISPOINT(2,NSEG))
         CALL TIMFUN(TABF,Fx2,VVREME,ITFMAX(INTAXISPOINT(2,NSEG+1))
     1                                      ,INTAXISPOINT(2,NSEG+1))
         FK1=Fx1+(CORD(INTERPAXIS,NZAD(1,I))-AXISPTCORD(NSEG))*
     1         ((Fx2-Fx1)/(AXISPTCORD(NSEG+1)-AXISPTCORD(NSEG)))
!          WRITE(3,*) "FK1,Fx1,Fx2",FK1,Fx1,Fx2
        CASE (3)
!          WRITE(3,*) "INTAXISPOINT(2,1)",INTAXISPOINT(2,1)
         CALL TIMFUN(TABF,FK1,VVREME,ITFMAX(INTAXISPOINT(2,
     1                NUMAXISPOINTS)),INTAXISPOINT(2,NUMAXISPOINTS))
!          WRITE(3,*) "FK1",FK1
        CASE DEFAULT
         PRINT*, "OUT OF INTERPOLATION SEGMENTS"
         STOP
        END SELECT
!           WRITE(3,*) "FK1,Fx1,Fx2",FK1,Fx1,Fx2

      RETURN
      END
C=======================================================================
      SUBROUTINE PREDF_INTERP(I,NZAD,CORD,FK1,IND)
      integer*4 NZAD(3,*),I,IAXIS,IND
!       integer*8 ILOC,INDLOC,NSEG,NUMAXISPOINTS
                         
      double precision CORD(3,*),FK1
      double precision Fkf1,Fx2,Fn1,Fakt,Akonst,Cin,y1,y2
      
!       DO ILOC=2,NUMAXISPOINTS
!          IF(AXISPTCORD(1).GE.CORD(INTERPAXIS,NZAD(1,I))) THEN
!             INDLOC=1
!             EXIT
!          ENDIF
!          IF(AXISPTCORD(NUMAXISPOINTS).LE.CORD(INTERPAXIS,
!      1                                      NZAD(1,I))) THEN
!             INDLOC=3
!             EXIT
!          ENDIF
        IAXIS=2
        IF(IND.EQ.1) THEN
          Fakt=80.0
          Akonst=301.0
          Cin=0.332
          y1=125.0
          y2=500.0
        ELSEIF(IND.EQ.2) THEN
          Fakt=30.0
          Akonst=0.0
          Cin=1.0
          y1=150.0
          y2=500.0
        ENDIF
!         IF(AXISPTCORD(1).LE.CORD(INTERPAXIS,NZAD(1,I)).AND.
!      1         AXISPTCORD(2).GE.CORD(INTERPAXIS,NZAD(1,I))) THEN
!             INDLOC=2
!             NSEG=ILOC-1
!          Fkf1=Fakt/(AXISPTCORD(1)-AXISPTCORD(2))
!          CALL TIMFUN(TABF,Fx2,VVREME,ITFMAX(INTAXISPOINT(2,2))
!      1                                      ,INTAXISPOINT(2,2))
         Fkf1=Fakt/(y1-y2)
         Fn1=Akonst+Cin*FK1-Fkf1*y1
!          Write(3,*) 'Fkf1,Fn1,FK1', Fkf1,Fn1,FK1
         FK1=Fkf1*CORD(IAXIS,NZAD(1,I))+Fn1
!         ELSE
!          PRINT*, "OUT OF INTERPOLATION SEGMENTS"
!         STOP
! !             EXIT
!         ENDIF
!        ENDDO
!         WRITE(3,*) "INDLOC",INDLOC

      RETURN
      END
