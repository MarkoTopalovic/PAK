C$DEBUG      
C==========================================================================
C==========================================================================
      SUBROUTINE VMAXATE(MAXA,MHT,NE,NTE,JEDN,LM,NWK)
      USE NODES
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS    PODPROGRAM ZA FORMIRANJE VEKTORA VISINA STUBOVA I MAXA
CS    KONACNO SE SMESTAJU U ISTI PROSTOR
CE    PROGRAM TO DETERMINE COLUMN HEIGHTS VECTOR AND MAXA
C
      COMMON /DODAT/ NDIMM
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
C
      DIMENSION MAXA(*),MHT(*),LM(*)
      COMMON /CDEBUG/ IDEBUG
      integer*8 LM

      IF(IDEBUG.GT.0) WRITE(*,*) 'MAXATE'
     
C
CS    PETLJA PO ELEMENTIMA
CE    ELEMENT LOOP
C

C      DO I=1,36
C       LM(I)=0
C      ENDDO

      DO I=1,JEDN+1
       MHT(I)=0
       MAXA(I)=0
      ENDDO
!        write(*,*) "u MAXATE, pre petlje po el"

      DO 100 NLM=1,NE
      IF(elemtip(NLM).gt.100) THEN
        NTYPE=elemtip(NLM)/100
        NTIMES=100
      ELSE
        NTYPE=elemtip(NLM)/10
        NTIMES=10
      ENDIF
      NDIMEL=elemtip(NLM)-NTIMES*NTYPE
         KK=0
         DO 2 I=1,NDIMEL
            IF(NEL(I,NLM).EQ.0) GO TO 2
            N=NEL(I,NLM)
!         write(*,*) "NEL(I,NLM)",I,NLM,NEL(I,NLM),ID(J,N)
               DO 1 J=1,1
                  IF(ID(J,N).LE.0) GO TO 1
                  KK=KK+1
                  LM(KK)=ID(J,N)
    1          CONTINUE
C            ENDIF
    2    CONTINUE
C         CALL VIWRR(LM,KK,'  LM')
C
         LS=JEDN+1
         DO 10 I=1,KK
            IF (LM(I).LT.LS) LS=LM(I)
   10    CONTINUE

C         WRITE(IIZLAZ,*)'LS=',LS
C
         DO 20 I=1,KK
            II=LM(I)
            ME=II-LS
            IF(ME.GT.MHT(II)) MHT(II)=ME
C         WRITE(IIZLAZ,*)'MHT(II),II',MHT(II),II
   20    CONTINUE
C
  100 CONTINUE
!        write(*,*) "u MAXATE, posle petlje po el"
C
CS    VEKTOR MAXA
CE    VECTOR MAXA
C
      MAXA(1)=1
      MAXA(2)=2
      DO 200 I=2,JEDN
       MAXA(I+1)=MAXA(I)+MHT(I)+1
 200  CONTINUE
      NWK=MAXA(JEDN+1)-1
      LS = JEDN+1
      DO 210 I=1,LS
 210  MHT(I)=MAXA(I)



      RETURN
      END
C==========================================================================
C==========================================================================
      SUBROUTINE VOTVORI
CS OTVARANJE DATOTEKA
CE OPEN FILES
      COMMON /IME/ IME
      COMMON /SRPSKI/ ISRPS
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /CERNE/ IGRESK,IPROTK,IZSILE
      CHARACTER*50 IME
C      LOGICAL OLDNEW
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' OTVORI'
C
CS ULAZNA DATOTEKA
CE INPUT FILE
C
      IF(ISRPS.EQ.0)
     *WRITE(*,2000)
      IF(ISRPS.EQ.1)
     *WRITE(*,6000)
 2000 FORMAT('   UNETI IME ULAZNE DATOTEKE / "pakv.dat" :')
 6000 FORMAT('   ENTER NAME OF INPUT FILE / "pakv.dat" :')
      READ (*,910) IME
      IF(IME.EQ.'                    ') IME = 'pakv.dat            '
  910 FORMAT (A)
      CALL VIMENA(IME)
      OPEN (IULAZ,FILE=IME,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='SEQUENTIAL')
C
CS  SCRATCH DATOTEKE
CE  SCRATCH FILES
C
      OPEN (15,STATUS='SCRATCH',FORM='UNFORMATTED',
     1 ACCESS='SEQUENTIAL')
C otvara ZTEMP datoteku samo za pakt     
      IF(IPAKT.EQ.1)THEN
      OPEN (19,FILE='ZITEMP',STATUS='UNKNOWN',FORM='UNFORMATTED',
     1 ACCESS='SEQUENTIAL')
      ENDIF
      OPEN (20,STATUS='SCRATCH',FORM='UNFORMATTED',
     1 ACCESS='SEQUENTIAL')
C
      RETURN
      END
C======================================================================
      SUBROUTINE VISPITA(ACOZ)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO TEST INPUT CARD
CS.    P R O G R A M
CS.        ZA ISPITIVANJE ULAZNE KARTICE
CS.        ZA BROJANJE UCITANIH KARTICA - KARTIC=KARTIC
C .
CS.        AKO U PRVE DVE KOLONE KARTICE STOJI: 'C ' ILI 'C-' ILI 'C*'
CS.        KARTICA SE IGNORISE
C .
C ......................................................................
C
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
       COMMON /BROJUK/ INDFOR,NULAZ
C      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
C     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG

      CHARACTER*250 ACOZ


      IF(IDEBUG.GT.0) PRINT *, ' ISPITA'

C   10 KARTIC = KARTIC + 1
  10  READ(IULAZ,1000) ACOZ
      IF(ACOZ(1:2).EQ.'C '.OR.ACOZ(1:2).EQ.'C*'.OR.ACOZ(1:2).EQ.'C-'.OR.
     1ACOZ(1:2).EQ.'c ') GO TO 10
      IF(INDFOR.EQ.2) RETURN
      BACKSPACE 1
      RETURN
 1000 FORMAT(A250)
      END

C.......................................................................
C.                                                                     .
C=======================================================================
      SUBROUTINE VZATVOR
CS ZATVARANJE DATOTEKA
CE CLOSE FILES
      COMMON /PRIMER/ IPRBR,INDIZL,INDGRA
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' ZATVOR'
C
      CLOSE (1,STATUS='KEEP')
      IF(INDIZL.EQ.1) CLOSE (IIZLAZ,STATUS='KEEP')
      IF(INDGRA.EQ.1) THEN
        CLOSE (18,STATUS='KEEP')
      ENDIF
      RETURN
      END
C======================================================================
      SUBROUTINE VADDSTF(A,B,C,S,P,JDIAG,LD,NEL,INDSK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C.... ASSEMBLE GLOBAL ARRAYS
C
      DIMENSION A(*),B(*),JDIAG(*),P(*),S(NEL,*),LD(*),C(*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' ADDSTF'
      DO 200 J = 1,NEL
      K = LD(J)
      IF(K.EQ.0) GO TO 200
      B(K) = B(K) + P(J)
      IF(INDSK.EQ.0) GO TO 200
      L = JDIAG(K) - K
      DO 100 I = 1,NEL
      M = LD(I)
      IF(M.GT.K.OR.M.EQ.0) GO TO 100
      M = L + M
      A(M) = A(M) + S(I,J)
      C(M) = C(M) + S(J,I)
  100 CONTINUE
  200 CONTINUE
      RETURN
      END    
C======================================================================
      SUBROUTINE VUACTCL(A,C,B,JDIAG,NEQ,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION A(*),B(*),JDIAG(*),C(*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' UACTCL'
C
C.... UNSYMMETRIC,ACTIVE COLUMN PROFILE EQUATION SOLVER
C
      IF(KKK.EQ.1)THEN
        DO 800 J = 1,NEQ
        JD = JDIAG(J)
  800   C(JD)=1.D0
      ENDIF
C
C.... FACTOR A TO UT*D*U REDUCE B TO Y
C
      JR = 0
      DO 300 J = 1,NEQ
      JD = JDIAG(J)
      JH = JD - JR
      IF(JH.LE.1) GO TO 300
      IS = J + 1 - JH
      IE = J - 1
      IF(KKK.EQ.2) GO TO 250
      K = JR + 1
      ID = 0
C
C.... REDUCE ALL EQUATIONS EXCEPT DIAGONAL
C
      DO 200 I = IS,IE
      IR = ID
      ID = JDIAG(I)
      IH = MIN0(ID - IR - 1,I - IS)
      IF(IH.EQ.0) GO TO 150
      A(K) = A(K) - DOT(A(K - IH),C(ID - IH),IH)
      C(K) = C(K) - DOT(C(K - IH),A(ID - IH),IH)
  150 IF(DABS(A(ID)).GT.1.D-20) C(K) = C(K)/A(ID)
C 150 IF(DABS(A(ID)).GT.1.D-40) C(K) = C(K)/A(ID)
  200 K = K + 1
C
C.... REDUCE DIAGONAL TERM
C
      A(JD) = A(JD) - DOT(A(JR + 1),C(JR + 1),JH - 1)
C
C.... FORWARD REDUCE THE R.H.S.
C
  250 IF(KKK.EQ.2) B(J) = B(J) - DOT(C(JR + 1),B(IS),JH - 1)
  300 JR = JD
      IF(KKK.EQ.1) RETURN
C
C.... BACKSUBSTITUTION
C
      J = NEQ
      JD = JDIAG(J)
  500 IF(DABS(A(JD)).GT.1.D-20) B(J) = B(J)/A(JD)
C 500 IF(DABS(A(JD)).GT.1.D-40) B(J) = B(J)/A(JD)
      D = B(J)
      J = J - 1
      IF(J.LE.0) RETURN
      JR = JDIAG(J)
      IF(JD - JR.LE.1) GO TO 700
      IS = J - JD + JR + 2
      K = JR - IS + 1
      DO 600 I = IS,J
      B(I) = B(I) - A(I + K)*D
  600 CONTINUE
  700 JD = JR
      GO TO 500
      END                 
C======================================================================
      FUNCTION VDOT(A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C...  SCALAR PRODUCT      
      DIMENSION A(*),B(*)
      DOT=0.D0
      DO 100 I=1,N
      DOT = DOT + A(I)*B(I)
  100 CONTINUE
      RETURN
      END
C======================================================================
      SUBROUTINE VKONVTF(TT1,F,KONV,N,ID,ITER)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /TACNOS/ EPSTR,MAXIT,NJRAP
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /SRPSKI/ ISRPS

      DIMENSION TT1(*),F(*),ID(1,*)
      KONV=1

      DTMOD=0.
      TMAX=0.
   
      DO 10 J=1,NPT
      I=ID(N,J)
      IF (I.EQ.0) GO TO 10
      TMAX=TMAX+TT1(I)**2
      DTMOD=DTMOD+F(I)**2
 10   CONTINUE       

      IF (DTMOD.GT.1.D-10) DTMOD=DSQRT(DTMOD)
      IF (TMAX.GT.1.D-10) TMAX=DSQRT(TMAX)

      TOLT=TMAX*EPSTR

      IF (DTMOD.GT.TOLT) KONV=0
      IF (TMAX.LT.1.D-10) THEN
      KONV=1
      ELSE 
      ENE=DTMOD/TMAX
       IF(ISRPS.EQ.0)THEN
         WRITE(IIZLAZ,1000) ITER,DTMOD,ITER,ENE
         WRITE(*,1000) ITER,DTMOD,ITER,ENE
       END IF
        IF(ISRPS.EQ.1)THEN
         WRITE(IIZLAZ,2000) ITER,DTMOD,ITER,ENE
         WRITE(*,2000) ITER,DTMOD,ITER,ENE
        END IF
      ENDIF
 1000 FORMAT(/'      ITERACIJA = ',I5,5X,'NORMA  = ',1PD12.4/
     1        '        NORMA (',I3,') / NORMA ( 0 )  = ',1PD12.4)
 2000 FORMAT(/'      ITERATION = ',I5,5X,'NORM   = ',1PD12.4/
     1        '          NORM (',I3,') / NORM ( 0 )  = ',1PD12.4)
      END
C======================================================================

C======================================================================
      SUBROUTINE VINTERP1D(R,KFIX,NPARAM,NDIM,AJS)
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DET,DET1,NBREL
      COMMON /TREN1/ FPS1,FPS2,FP1,FP2,VNN,PERM,CORDY
      COMMON /VISKOZ/ AMI
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP

      DIMENSION ZVHR1D(3)
      DIMENSION HP(4)
      DIMENSION TT21(44),TT210(44),H2N(2,18),H2NP(2,4)
      DIMENSION B(4,18),BT(18,4),GRADH(2,4),GRADHT(4,2),ZHP(4,2)
      DIMENSION GH2NT(18,2)
      DIMENSION H2NN(2,18)
      DIMENSION PERM(2,2),CORDY(4)

      DIMENSION AJ(1,3)

      RP=1.D0+R
      RM=1.D0-R
      RR=1.D0-R*R

      IF (NDIM.GT.2) THEN
         H(3)=RR
      ELSE
         H(3)=0.D0
      ENDIF
      H(2)=0.5*RM-0.5*H(3)
      H(1)=0.5*RP-0.5*H(3)

      IF (NDIM.GT.2) THEN
         ZVHR1D(3)=-2.D0*R
      ELSE
         ZVHR1D(3)=0.D0
      ENDIF
      ZVHR1D(2)=-0.5-0.5*ZVHR1D(3)
      ZVHR1D(1)=0.5-0.5*ZVHR1D(3)
      
!       AJ(1,1)=DOT(ZVHR1D,CK(1,1),NDIM)
!       AJ(1,2)=DOT(ZVHR1D,CK(2,1),NDIM)
!       AJ(1,3)=DOT(ZVHR1D,CK(3,1),NDIM)

      DO I=1,3
        AJ(1,I)=0.D0
      ENDDO
      DO I=1,NDIM
        DO J=1,3
          AJ(1,J)=AJ(1,J)+ZVHR1D(I)*CK(J,I)
        enddo
      ENDDO

!       DO 61 I=1,3
!       AJ(1,I)=0.
!       DO 61 KK=1,NDIM
!       AJ(1,I)=AJ(1,I)+ZVHR1D(KK)*CK(I,KK)
!    61 CONTINUE
      
      AJS = DSQRT(AJ(1,1)**2+AJ(1,2)**2+AJ(1,3)**2)

      IF(KFIX.GT.0) GO TO 70
      
      X=CK(1,2)-CK(1,1)
      Y=CK(2,2)-CK(2,1)
      Z=CK(3,2)-CK(3,1)
      DL=DSQRT(X**2+Y**2+Z**2)
      ABH=DSQRT(X**2+Y**2)
      BCH=DSQRT(Y**2+Z**2)
      CAH=DSQRT(Z**2+X**2)
      if(ABH.eq.0) then
        if(X.gt.0) then
            CSA=1
        else
            CSA=-1
        endif
        if(Y.gt.0) then
            SA=1
        else
            SA=-1
        endif
      ELSE
        CSA=X/ABH
        SA=Y/ABH
        CSDA=ABH/DL
      ENDIF
      if(BCH.eq.0) then
        if(Y.gt.0) then
            CSB=1
        else
            CSB=-1
        endif
        if(Z.gt.0) then
            SB=1
        else
            SB=-1
        endif
      ELSE
        CSB=Y/BCH
        SB=Z/BCH
        CSDB=BCH/DL
      ENDIF
      if(CAH.eq.0) then
        if(X.gt.0) then
            CSC=1
        else
            CSC=-1
        endif
        if(Z.gt.0) then
            SC=1
        else
            SC=-1
        endif
      ELSE
        CSC=X/CAH
        SC=Z/CAH
        CSDG=CAH/DL
      ENDIF

      DO 20 I=1,NDIM
         if(ABH.eq.0) then
            ZVHX(I)=(1/AJS)*ZVHR1D(I)*CSA*CSC
         else
            ZVHX(I)=(1/AJS)*ZVHR1D(I)*CSA*CSDA
         endif
         if(BCH.eq.0) then
            ZVHY(I)=(1/AJS)*ZVHR1D(I)*SA*CSB
         else
            ZVHY(I)=(1/AJS)*ZVHR1D(I)*CSB*CSDB
         endif
         if(CAH.eq.0) then
            ZVHZ(I)=(1/AJS)*ZVHR1D(I)*SC*SB
         else
            ZVHZ(I)=(1/AJS)*ZVHR1D(I)*SC*CSDG
         endif
   20 CONTINUE  

   70 IF (NPARAM.EQ.1) THEN
        IF (NDIM.EQ.4) PJ=0.5D0*(FPS1+FPS2)
        FS2=PJ
      ENDIF
C
      END
C======================================================================
C======================================================================
      SUBROUTINE VGETRST(R,S,T,W,NGR,NGIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION W(3)
      
      PET=5.
      IF(NGIP.EQ.1) THEN
          R = 0.25
          S = 0.25
          T = 0.25 
          W(1)= 1./6.
          W(2)= 1.
          W(3)= 1.
      ENDIF
      IF(NGIP.EQ.4) THEN
          AA=(5.+3.*DSQRT(PET))/20.
          BB=(5.-   DSQRT(PET))/20.
          IF(NGR.EQ.1) THEN
              R=BB
              S=BB
              T=BB
          ENDIF
          IF(NGR.EQ.2) THEN
              R=AA
              S=BB
              T=BB
          ENDIF
          IF(NGR.EQ.3) THEN
              R=BB
              S=AA
              T=BB
          ENDIF
          IF(NGR.EQ.4) THEN
              R=BB
              S=BB
              T=AA
          ENDIF
          W(1)=1./24.
          W(2)=1.
          W(3)=1.
      ENDIF
      IF(NGIP.EQ.5) THEN
        IF(NGR.EQ.1) THEN
            R= 1./4
            S= 1./4
            T= 1./4
            W(1)=-4./30
        ENDIF
        IF(NGR.EQ.2) THEN
            R= 1./6
            S= 1./6
            T= 1./6
            W(1)= 9./120
        ENDIF
        IF(NGR.EQ.3) THEN
            R= 1./2
            S= 1./6
            T= 1./6
            W(1)= 9./120
        ENDIF
        IF(NGR.EQ.4) THEN
            R= 1./6
            S= 1./2
            T= 1./6
            W(1)= 9./120
        ENDIF
        IF(NGR.EQ.5) THEN
            R= 1./6
            S= 1./6
            T= 1./2
            W(1)= 9./120
        ENDIF
        W(2)= 1.
        W(3)= 1.
      ENDIF
               
      RETURN
      END
C======================================================================
C======================================================================
      SUBROUTINE VJACTP33(R,S,HT,IND,NGPSIL,JBRPS,NDIM)
      USE NODES
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     POVRSINSKI JAKOBIJAN NA GRANICI 3/D TETRA ELEMENTA
C
C
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DET,DET1,NBREL
      COMMON /TREN13/ FPS1,FPS2,FPS3,FPS4
      COMMON /VISKOZ/ AMI
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /CDEBUG/ IDEBUG
      DIMENSION XJ(3,3),HT(7,3),NGPSIL(12,*)
!       INTEGER*8 NGPSIL
C
      IF(IDEBUG.GT.0) PRINT *, ' JACTP33'
      
      HT(1,1)=1.-R-S
      HT(2,1)=R
      HT(3,1)=S
      HT(4,1)=0.
      HT(1,2)=-1.
      HT(2,2)=1.
      HT(3,2)=0.
      HT(1,3)=-1.
      HT(2,3)=0.
      HT(3,3)=1.
C
C     KOREKCIJA FUNKCIJA KADA JE BROJ CVOROVA VECI OD 3
C
   26 IF(NDIM.eq.10) THEN
C       CETVRTI CVOR
        HT(5,1)=4.*R*(1.-R-S)
        HT(5,2)=4.*(1.-R-S)-4.*R
        HT(5,3)=-4.*R
C       PETI CVOR
        HT(6,1)=4.*R*S
        HT(6,2)=4.*S
        HT(6,3)=4.*R
C       SESTI CVOR
        HT(7,1)=4.*S*(1.-R-S)
        HT(7,2)=-4.*S
        HT(7,3)=4.*(1.-R-S)-4.*S
        
C
C
C       KOREKCIJA PRVA TRI CVORA
C
        HT(1,1)=HT(1,1)-.5*HT(5,1)
        HT(2,1)=HT(2,1)-.5*HT(5,1)
        HT(1,2)=HT(1,2)-.5*HT(5,2)
        HT(2,2)=HT(2,2)-.5*HT(5,2)
        HT(1,3)=HT(1,3)-.5*HT(5,3)
        HT(2,3)=HT(2,3)-.5*HT(5,3)
        
        HT(2,1)=HT(2,1)-.5*HT(6,1)
        HT(3,1)=HT(3,1)-.5*HT(6,1)
        HT(2,2)=HT(2,2)-.5*HT(6,2)
        HT(3,2)=HT(3,2)-.5*HT(6,2)
        HT(2,3)=HT(2,3)-.5*HT(6,3)
        HT(3,3)=HT(3,3)-.5*HT(6,3)
    
        HT(3,1)=HT(3,1)-.5*HT(7,1)
        HT(1,1)=HT(1,1)-.5*HT(7,1)
        HT(3,2)=HT(3,2)-.5*HT(7,2)
        HT(1,2)=HT(1,2)-.5*HT(7,2)
        HT(3,3)=HT(3,3)-.5*HT(7,3)
        HT(1,3)=HT(1,3)-.5*HT(7,3)
      
      ENDIF
C
C     JAKOBIJAN
C
  400 DO 60 I=1,3
         XJ(1,I)=0.D0
         XJ(2,I)=0.D0
      DO 60 J=1,7
         NJ=NGPSIL(J+1,JBRPS)
         IF(NJ.EQ.0) GO TO 60
         XJ(1,I)=XJ(1,I)+HT(J,2)*CORD(I,NJ)
         XJ(2,I)=XJ(2,I)+HT(J,3)*CORD(I,NJ)
   60 CONTINUE
C
C     DETERMINANTA JAKOBIJANA
C
      DETX=XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2)
      DETY=XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1)
      DETZ=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
      DET=DSQRT(DETX**2+DETY**2+DETZ**2)
      IF(DET.GT.1.0D-10) GO TO 70
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2000) JBRPS,DET
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6000) JBRPS,DET
      STOP ' PROGRAM STOP - Pakv2.for - JACTP33'
      
!    70 IF (NDIM.EQ.10) PRJ=(FPS1+FPS2+FPS3)/3
   70 PRJ=(FPS1+FPS2+FPS3)/3
      FS2=PRJ
      
      IF(IND.EQ.1) RETURN
      CINA=DETX/DET
      CINB=DETY/DET
      CING=DETZ/DET
      
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(/' ','DETERMINANTA JAKOBIJANA MANJA OD NULE'/
     1' ','3-D PRITISAK BROJ =',I10/' ','DETERMINANTA =',1PD12.5)
 6000 FORMAT(/' ','ZERO OR NEGATIVE DETERMINANTE'/
     1' ','3-D PRESSURE NUMBER =',I10/' ','DETERMINANT =',1PD12.5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C======================================================================
      SUBROUTINE VINTERPTETRA(R,S,T,KFIX,NPARAM,NDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       FOR INTERPOLATION FUNCTIONS AND JACOBIAN MATRIX IN CURRENT
CE.       INTEGRATION POINT (R,S,T - ARE NATURAL COORDINATES)
C .
C ......................................................................
C
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DET,DET1,NBREL
      COMMON /TREN13/ FPS1,FPS2,FPS3,FPS4
      COMMON /VISKOZ/ AMI
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /CDEBUG/ IDEBUG
      DIMENSION P(NDIM,4),XJ(3,3),IPERM(8),PJ(3,NDIM)
      DATA IPERM/2,3,4,1,6,7,8,5/
C
      IF(IDEBUG.GT.0) PRINT *, ' INTERPTETRA'
C
C     INTERPOLACIJSKE FINKCIJE I NJIHOVI IZVODI
C
C      WRITE(3,1001) NCVE
C      DO I=1,47
C        WRITE(3,1001) (NOP(I,J), J=1,10)
C      ENDDO
C 1001 FORMAT(10I5)
C      STOP
C
C     PRVIH 4 CVOROVA
C
      IF(KFIX.eq.4.or.KFIX.eq.1.or.KFIX.eq.2.or.KFIX.eq.3) 
     1   WRITE(*,*) 'iz interp tetra'
      H(1)= 1-R-S-T
      H(2)= R
      H(3)= S
      H(4)= T
C
      P(1,2)=-1.
      P(2,2)= 1.
      P(3,2)= 0.
      P(4,2)= 0.
C
      P(1,3)=-1.
      P(2,3)= 0.
      P(3,3)= 1.
      P(4,3)= 0.
C
      P(1,4)=-1.
      P(2,4)= 0.
      P(3,4)= 0.
      P(4,4)= 1.
C
      IF(NDIM.EQ.4) GO TO 50
      IF(NDIM.NE.10) STOP 'NDIM.NE.10'
!       I=0
!    55 I=I+1
!       NN=I+4
!       IF(NOP(NLM,NN).EQ.0) GO TO 55
!       GO TO (5,6,7,8,9,10),I
C
C     STEPENE SLOBODE ZA CVOROVE PREKO 4
C
C     PETI CVOR
    5 H(5)= 4.*R*(1.-R-S-T)
      P(5,2)= 4.*(1.-R-S-T)-4.*R
      P(5,3)=-4.*R
      P(5,4)=-4.*R
!       GO TO 55
C     SESTI CVOR
    6 H(6)= 4.*R*S
      P(6,2)= 4.*S
      P(6,3)= 4.*R
      P(6,4)= 0.
!       GO TO 55
C     SEDMI CVOR
    7 H(7)= 4.*S*(1.-R-S-T)
      P(7,2)=-4.*S
      P(7,3)= 4.*(1.-R-S-T)-4.*S
      P(7,4)=-4.*S
!       GO TO 55
C     OSMI CVOR
    8 H(8)= 4.*T*(1.-R-S-T)
      P(8,2)=-4.*T
      P(8,3)=-4.*T
      P(8,4)= 4.*(1.-R-S-T)-4.*T
!       GO TO 55
C     DEVETI CVOR
    9 H(9)= 4.*R*T
      P(9,2)= 4.*T
      P(9,3)= 0.
      P(9,4)= 4.*R
!       GO TO 55
C     DESETI CVOR
   10 H(10)= 4.*S*T
      P(10,2)= 0.
      P(10,3)= 4.*T
      P(10,4)= 4.*S
C
C     KOREKCIJE PRVIH 4 FUNKCIJA AKO SU UPOTREBLJENI CVOROVI PREKO 4
C
C  MEDJUCVOROVI OD 5 DO 10
      H(1)=H(1)-0.5*(H(5)+H(7)+H(8))
      H(2)=H(2)-0.5*(H(5)+H(6)+H(9))
      H(3)=H(3)-0.5*(H(6)+H(7)+H(10))
      H(4)=H(4)-0.5*(H(8)+H(9)+H(10))
      DO 210 J=2,4
        P(1,J)=P(1,J)-0.5*(P(5,J)+P(7,J)+P(8,J))
        P(2,J)=P(2,J)-0.5*(P(5,J)+P(6,J)+P(9,J))
        P(3,J)=P(3,J)-0.5*(P(6,J)+P(7,J)+P(10,J))
        P(4,J)=P(4,J)-0.5*(P(8,J)+P(9,J)+P(10,J))
  210 CONTINUE
C
C     JAKOBIJAN U TACKI R,S,T
C
   50 DO 60 I=1,3
      DO 60 J=1,3
      XJ(I,J)=0.
      DO 60 KM=1,NDIM
      XJ(I,J)=XJ(I,J)+P(KM,I+1)*CK(J,KM)
   60 CONTINUE

      CALL VMINV3(XJ,DET1)
      
   77 DO 85 I=1,3
      DO 85 JJ=1,NDIM
      PJ(I,JJ)=0.
      DO 85 K=1,3
      PJ(I,JJ)=PJ(I,JJ) + XJ(I,K)*P(JJ,K+1)
   85 CONTINUE
C
      DO 123 J=1,NDIM
      ZVHX(J)=PJ(1,J) 
      ZVHY(J)=PJ(2,J) 
 123  ZVHZ(J)=PJ(3,J) 
C
C     DERERMINANTA JAKOBIJANA U TACKI R,S,T
C
      IF(DET1.GT.1.D-13) RETURN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) NLM,KFIX,R,S,T,DET1
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) NLM,KFIX,R,S,T,DET1
        write(izlaz,*) 'NDIM',NDIM

      STOP 'PROGRAM STOP - PAK32 - JACTE33A'
C
 1000 FORMAT(4I10,3(1PE13.5))
 1100 FORMAT(8I10)
 1200 FORMAT(' h',4(1PE13.5))
C-----------------------------------------------------------------------
 2000 FORMAT(' ** GRESKA **: NEGATIVNA ILI NULA DETERMINATA JAKOBIJANA',
     1       ' ZA ELEMENT BR.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(/' ','ZERO OR NEGATIVE JACOBIAN DETERMINANTE'/
     1' ','ELEMENT NUM. =',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',D12.5)
C-----------------------------------------------------------------------
      END
C======================================================================
C======================================================================
      SUBROUTINE VINTERPTETRAF(R,S,T,KFIX,NDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       FOR INTERPOLATION FUNCTIONS AND JACOBIAN MATRIX IN CURRENT
CE.       INTEGRATION POINT (R,S,T - ARE NATURAL COORDINATES)
C         Prescribed flux or covnection
C .
C ......................................................................
C
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DET,DET1,NBREL
      COMMON /TREN13/ FPS1,FPS2,FPS3,FPS4
      COMMON /VISKOZ/ AMI
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      DIMENSION P(NDIM,4),XJ(3,3),IPERM(8),PJ(3,NDIM)
      DATA IPERM/2,3,4,1,6,7,8,5/
C
      IF(IDEBUG.GT.0) PRINT *, ' INTERPTETRAF'
C
C     INTERPOLACIJSKE FINKCIJE I NJIHOVI IZVODI
C
C      WRITE(3,1001) NCVE
C      DO I=1,47
C        WRITE(3,1001) (NOP(I,J), J=1,10)
C      ENDDO
C 1001 FORMAT(10I5)
C      STOP
C
C     PRVIH 4 CVOROVA
C
!       IF(KFIX.eq.4.or.KFIX.eq.1.or.KFIX.eq.2.or.KFIX.eq.3) 
!      1   WRITE(*,*) 'iz interp tetra, KFIX', KFIX
      H(1)= 1-R-S-T
      H(2)= R
      H(3)= S
      H(4)= T
C
C
      IF(NDIM.EQ.4) GO TO 50
      IF(NDIM.NE.10) STOP 'NDIM.NE.10'
!       I=0
!    55 I=I+1
!       NN=I+4
!       IF(NOP(NLM,NN).EQ.0) GO TO 55
!       GO TO (5,6,7,8,9,10),I
C
C     STEPENE SLOBODE ZA CVOROVE PREKO 4
C
C     PETI CVOR
    5 H(5)= 4.*R*(1.-R-S-T)
C     SESTI CVOR
    6 H(6)= 4.*R*S
C     SEDMI CVOR
    7 H(7)= 4.*S*(1.-R-S-T)
C     OSMI CVOR
    8 H(8)= 4.*T*(1.-R-S-T)
C     DEVETI CVOR
    9 H(9)= 4.*R*T
C     DESETI CVOR
   10 H(10)= 4.*S*T
C
C     KOREKCIJE PRVIH 4 FUNKCIJA AKO SU UPOTREBLJENI CVOROVI PREKO 4
C     MEDJUCVOROVI OD 5 DO 10
      H(1)=H(1)-0.5*(H(5)+H(7)+H(8))
      H(2)=H(2)-0.5*(H(5)+H(6)+H(9))
      H(3)=H(3)-0.5*(H(6)+H(7)+H(10))
      H(4)=H(4)-0.5*(H(8)+H(9)+H(10))
C Izvodi interpolacionih funkcija
      P(1,2)=-1.
      P(2,2)= 1.
      P(3,2)= 0.
      P(4,2)= 0.
C
      P(1,3)=-1.
      P(2,3)= 0.
      P(3,3)= 1.
      P(4,3)= 0.
C
      P(1,4)=-1.
      P(2,4)= 0.
      P(3,4)= 0.
      P(4,4)= 1.
C     PETI CVOR
      P(5,2)= 4.*(1.-R-S-T)-4.*R
      P(5,3)=-4.*R
      P(5,4)=-4.*R
C     SESTI CVOR
      P(6,2)= 4.*S
      P(6,3)= 4.*R
      P(6,4)= 0.
C     SEDMI CVOR
      P(7,2)=-4.*S
      P(7,3)= 4.*(1.-R-S-T)-4.*S
      P(7,4)=-4.*S
C     OSMI CVOR
      P(8,2)=-4.*T
      P(8,3)=-4.*T
      P(8,4)= 4.*(1.-R-S-T)-4.*T
C     DEVETI CVOR
      P(9,2)= 4.*T
      P(9,3)= 0.
      P(9,4)= 4.*R
C     DESETI CVOR
      P(10,2)= 0.
      P(10,3)= 4.*T
      P(10,4)= 4.*S
C
      DO 210 J=2,4
        P(1,J)=P(1,J)-0.5*(P(5,J)+P(7,J)+P(8,J))
        P(2,J)=P(2,J)-0.5*(P(5,J)+P(6,J)+P(9,J))
        P(3,J)=P(3,J)-0.5*(P(6,J)+P(7,J)+P(10,J))
        P(4,J)=P(4,J)-0.5*(P(8,J)+P(9,J)+P(10,J))
  210 CONTINUE
C
C     JAKOBIJAN U TACKI R,S,T
C
   50 DO 60 I=1,3
      DO 60 J=1,3
      XJ(I,J)=0.
      DO 60 KM=1,NDIM
      XJ(I,J)=XJ(I,J)+P(KM,I+1)*CK(J,KM)
   60 CONTINUE
! 
!       write(3,*) "NBREL,NDIM",NBREL,NDIM
       SELECT CASE (KFIX)
!  s=0
       CASE (1)
         DETX=XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1)
         DETY=XJ(1,1)*XJ(3,3)-XJ(1,3)*XJ(3,1)
         DETZ=XJ(2,1)*XJ(3,3)-XJ(2,3)*XJ(3,1)
         DET1=DSQRT(DETX**2+DETY**2+DETZ**2)
         IF(DET1.GT.1.0D-10) GO TO 77
         IF(ISRPS.EQ.0)
     1   WRITE(IZLAZ,2000) JBRPS,DET1
         IF(ISRPS.EQ.1)
     1   WRITE(IZLAZ,6000) JBRPS,DET1
         STOP ' PROGRAM STOP - Pakv2.for - intertetraf'
!  t=0
       CASE (2)
         DETX=XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2)
         DETY=XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1)
         DETZ=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
         DET1=DSQRT(DETX**2+DETY**2+DETZ**2)
         IF(DET1.GT.1.0D-10) GO TO 77
         IF(ISRPS.EQ.0)
     1   WRITE(IZLAZ,2000) JBRPS,DET1
         IF(ISRPS.EQ.1)
     1   WRITE(IZLAZ,6000) JBRPS,DET1
         STOP ' PROGRAM STOP - Pakv2.for - intertetraf'
!  r=0
       CASE (3)
         DETX=XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2)
         DETY=XJ(2,1)*XJ(3,3)-XJ(2,3)*XJ(3,1)
         DETZ=XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1)
         DET1=DSQRT(DETX**2+DETY**2+DETZ**2)
         IF(DET1.GT.1.0D-10) GO TO 77
         IF(ISRPS.EQ.0)
     1   WRITE(IZLAZ,2000) JBRPS,DET1
         IF(ISRPS.EQ.1)
     1   WRITE(IZLAZ,6000) JBRPS,DET1
         STOP ' PROGRAM STOP - Pakv2.for - intertetraf'
!  r+s+t=1
       CASE (4)
         CALL VMINV3(XJ,DET1)
       CASE DEFAULT
!         WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NBREL)
       END SELECT
!       CALL VMINV3(XJ,DET1)
      
   77 DO 85 I=1,3
      DO 85 JJ=1,NDIM
      PJ(I,JJ)=0.
      DO 85 K=1,3
      PJ(I,JJ)=PJ(I,JJ) + XJ(I,K)*P(JJ,K+1)
   85 CONTINUE
C
      DO 123 J=1,NDIM
      ZVHX(J)=PJ(1,J) 
      ZVHY(J)=PJ(2,J) 
 123  ZVHZ(J)=PJ(3,J) 
C
C     DERERMINANTA JAKOBIJANA U TACKI R,S,T
C
      IF(DET1.GT.1.D-13) RETURN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) NBREL,KFIX,R,S,T,DET1
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) NBREL,KFIX,R,S,T,DET1

      STOP 'PROGRAM STOP - PAKV2 - INTERTETRAF'
C
 1000 FORMAT(4I10,3(1PE13.5))
 1100 FORMAT(8I10)
 1200 FORMAT(' h',4(1PE13.5))
C-----------------------------------------------------------------------
 2000 FORMAT(' ** GRESKA **: NEGATIVNA ILI NULA DETERMINATA JAKOBIJANA',
     1       ' ZA ELEMENT BR.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(/' ','ZERO OR NEGATIVE JACOBIAN DETERMINANTE'/
     1' ','ELEMENT NUM. =',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',D12.5)
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE VMINV2(XJJ,DETT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO INVERSE MATRIX 2 X 2
CS.   P R O G R A M
CS.      ZA INVERTOVANJE MATRICE 2 X 2
C .
C ......................................................................
C
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DETJS,DETJ,NBREL
      DIMENSION XJJ(2,2)
C
      IF(IDEBUG.GT.0) PRINT *, ' MINV2 '
C
      DETT=XJJ(1,1)*XJJ(2,2)-XJJ(1,2)*XJJ(2,1)
      IF(DETT.GT.1.D-10) GOTO 10
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) DETT,NBREL
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) DETT,NBREL
      PRINT *, ' DET < 0'
      IF(DABS(DETT).LT.1.D-10) STOP 'MINV2'
   10 XJP=XJJ(1,1)
      XJJ(1,1)=XJJ(2,2)/DETT
      XJJ(2,2)=XJP/DETT
      XJJ(1,2)=-XJJ(1,2)/DETT
      XJJ(2,1)=-XJJ(2,1)/DETT
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(/
     1' DETERMINANTA MANJA ILI JEDNAKA NULI, DET = ',1PD10.3/
     1' ELEMENT =',I8)
C-----------------------------------------------------------------------
 6000 FORMAT(/
     1' ZERO OR NEGATIVE DETERMINANTE, DET = ',1PD10.3/
     1' ELEMENT =',I8)
C-----------------------------------------------------------------------
      END
C======================================================================
      
C======================================================================      
      SUBROUTINE VLCK2D(CK,CKL,TTE,NDIMEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      DIMENSION CK(3,21),CKL(3,8),A12(3),A13(3)
      DIMENSION EN(3),Y(3),TTE(2,3)
      
      DO I=1,8
       DO j=1,3
        CKL(J,I)=0.D0
       ENDDO
      ENDDO
      
      DO 620 I=1,3
         A12(I) = CK(I,2) - CK(I,1)
  620    A13(I) = CK(I,3) - CK(I,1)
C        VEKTOR NORMALE
         EN(1) = A12(2)*A13(3) - A12(3)*A13(2)
         EN(2) = A12(3)*A13(1) - A12(1)*A13(3)
         EN(3) = A12(1)*A13(2) - A12(2)*A13(1)
C        VEKTOR Y
         Y(1) = EN(2)*A12(3) - EN(3)*A12(2)
         Y(2) = EN(3)*A12(1) - EN(1)*A12(3)
         Y(3) = EN(1)*A12(2) - EN(2)*A12(1)
         YI = DSQRT(Y(1)*Y(1)+Y(2)*Y(2)+Y(3)*Y(3))
         A12I = DSQRT(A12(1)*A12(1)+A12(2)*A12(2)+A12(3)*A12(3))
         DO 640 I=1,3
            TTE(1,I) = A12(I)/A12I
  640       TTE(2,I) = Y(I)/YI
C           FORMIRANJE LOKALNIH KOORDINATA
            CKL(1,1)=0.0D0
            CKL(2,1)=0.0D0
            DO 660 I=2,NDIMEL
                DO 670 J=1,2
                    DO 670 K =1,3
  670       CKL(J,I) = CKL(J,I) + TTE(J,K)*(CK(K,I) - CK(K,1))
  660    CONTINUE
  
        RETURN
        END
C======================================================================
      SUBROUTINE VINTERP2D36(R,S,KFIX,NPARAM,NDIMEL,CKL,TTE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DETJS,DETJ,NBREL
      COMMON /TREN1/ FPS1,FPS2,FP1,FP2,VNN,PERM,CORDY
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /VISKOZ/ AMI
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP

      DIMENSION ZVHR(9),ZVHS(9)
      DIMENSION HP(4),X(9),Y(9),Z(9),CKL(3,8),TTE(2,3)
      DIMENSION ZVHXL(9),ZVHYL(9)
      DIMENSION TT21(44),TT210(44),H2N(2,18),H2NP(2,4)
      DIMENSION B(4,18),BT(18,4),GRADH(2,4),GRADHT(4,2),ZHP(4,2)
      DIMENSION GH2NT(18,2)
      DIMENSION H2NN(2,18)
      DIMENSION PERM(2,2),CORDY(4)

      DIMENSION AJ(2,2)
      
      DO 130 KLM=1,NDIMEL
       X(KLM)=CKL(1,KLM)
       Y(KLM)=CKL(2,KLM)
  130 CONTINUE

      RP=1.D0-R-S

      IF (NDIMEL.GT.3) THEN
         H(6)=4*S*RP
         H(5)=4*R*S
         H(4)=4*R*RP
         H(3)=S-0.5*(H(5)+H(6))
         H(2)=R-0.5*(H(4)+H(5))
         H(1)=RP-0.5*(H(4)+H(6))
      ELSE
         H(3)=S
         H(2)=R
         H(1)=RP
      ENDIF

      DEBLJ=1.D0
      IF (INDAX.EQ.1) THEN
      DEBLJ=0.D0
      DO I=1,NDIMEL
        DEBLJ=DEBLJ+H(I)*X(I)
      ENDDO
      ENDIF

      IF (NDIMEL.GT.3) THEN
         ZVHR(6)=-4*S
         ZVHS(6)=4*RP-4*S
         ZVHR(5)=4*S
         ZVHS(5)=4*R
         ZVHR(4)=4*RP-4*R
         ZVHS(4)=-4*R
         ZVHR(3)=-0.5*(ZVHR(5)+ZVHR(6))
         ZVHS(3)=1.0-0.5*(ZVHS(5)+ZVHS(6))
         ZVHR(2)=1.0-0.5*(ZVHR(4)+ZVHR(5))
         ZVHS(2)=-0.5*(ZVHS(4)+ZVHS(5))
         ZVHR(1)=-1.0-0.5*(ZVHR(4)+ZVHR(6))
         ZVHS(1)=-1.0-0.5*(ZVHS(4)+ZVHS(6))
      ELSE
         ZVHR(3)=0.0
         ZVHS(3)=1.0
         ZVHR(2)=1.0
         ZVHS(2)=0.0
         ZVHR(1)=-1.0
         ZVHS(1)=-1.0
      ENDIF

      AJ(1,1)=DOT(ZVHR,X,NDIMEL)
      AJ(1,2)=DOT(ZVHR,Y,NDIMEL)
      AJ(2,1)=DOT(ZVHS,X,NDIMEL)
      AJ(2,2)=DOT(ZVHS,Y,NDIMEL)

      IF(KFIX.GT.0) GO TO 70

      DETJ=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)

      IF (DETJ.LT.0.0D0) THEN
      WRITE(*,*)'DETERMINANTA MANJA OD NULE!!!'
      STOP
      ENDIF

      DO 20 I=1,NDIMEL
      ZVHXL(I)=(ZVHR(I)*AJ(2,2)-ZVHS(I)*AJ(1,2))/DETJ
      ZVHYL(I)=(ZVHS(I)*AJ(1,1)-ZVHR(I)*AJ(2,1))/DETJ
   20 CONTINUE
    
      DO I=1,NDIMEL
        ZVHX(I)=ZVHXL(I)*TTE(1,1)+ZVHYL(I)*TTE(2,1)
        ZVHY(I)=ZVHXL(I)*TTE(1,2)+ZVHYL(I)*TTE(2,2)
        ZVHZ(I)=ZVHXL(I)*TTE(1,3)+ZVHYL(I)*TTE(2,3)
      ENDDO
   
      RETURN

   70 IF (NPARAM.EQ.1) THEN

        IF (NDIMEL.EQ.4) PJ=0.5D0*(FPS1+FPS2)

        FS2=PJ

        IF (INDAX.EQ.1.AND.HH.GT.0.D0) THEN
           FS2=PJ/(2.D0*RIZ*3.14159D0*HH)
        ENDIF

      ENDIF
CE     LINE JACOBIAN DETERMINANT
      GO TO (71,72),KFIX
CE     CONSTANT KSI
   71 DET2= AJ(2,1)*AJ(2,1)+AJ(2,2)*AJ(2,2)
      GO TO 74
CE     CONSTANT ETA
   72 DET2= AJ(1,1)*AJ(1,1)+AJ(1,2)*AJ(1,2)
   74 DET2=DSQRT(DET2)
      DETJS=DET2
      IF ( DET2.GT.1.D-15) RETURN
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2000) NBREL,KFIX,R,S,DET2
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6000) NBREL,KFIX,R,S,DET2
      STOP
C
 2000 FORMAT(' ** GRESKA **: JAKOBIJAN JEDNAK ILI MANJI OD NULE',
     1       ' ZA ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     5       10X,'DET=',F10.5)
 6000 FORMAT(' ** ERROR **: JACOBIAN EQUAL OR LESS THEN ZERO',
     1       ' FOR ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     5       10X,'DET=',F10.5)
C

      END
C======================================================================
C======================================================================
      SUBROUTINE VINTERP2D48(R,S,KFIX,NPARAM,NDIMEL,CKL,TTE)
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DETJS,DETJ,NBREL
      COMMON /TREN1/ FPS1,FPS2,FP1,FP2,VNN,PERM,CORDY
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /VISKOZ/ AMI
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP

      DIMENSION ZVHR(9),ZVHS(9)
      DIMENSION ZVHXL(9),ZVHYL(9),TTE(2,3)
      DIMENSION HP(4),X(9),Y(9)
      DIMENSION TT21(44),TT210(44),H2N(2,18),H2NP(2,4)
      DIMENSION B(4,18),BT(18,4),GRADH(2,4),GRADHT(4,2),ZHP(4,2)
      DIMENSION GH2NT(18,2)
      DIMENSION H2NN(2,18)
      DIMENSION PERM(2,2),CORDY(4)

      DIMENSION AJ(2,2),XJ(2,2),XJJ(2,2),XJ3D(3,3),CKL(3,8)

      RP=1.D0+R
      SP=1.D0+S
      RM=1.D0-R
      SM=1.D0-S
      RR=1.D0-R*R
      SS=1.D0-S*S

      IF (NDIMEL.GT.4) THEN
         IF (NDIMEL.EQ.8) THEN
           H(9)=0.D0
         ELSE
           H(9)=RR*SS
         ENDIF 
         H(8)=0.5*RP*SS-0.5*H(9)
         H(7)=0.5*RR*SM-0.5*H(9)
         H(6)=0.5*RM*SS-0.5*H(9)
         H(5)=0.5*RR*SP-0.5*H(9)
         H(4)=0.25*RP*SM-0.5*(H(7)+H(8))-0.25*H(9)
         H(3)=0.25*RM*SM-0.5*(H(6)+H(7))-0.25*H(9)
         H(2)=0.25*RM*SP-0.5*(H(5)+H(6))-0.25*H(9)
         H(1)=0.25*RP*SP-0.5*(H(5)+H(8))-0.25*H(9)

      ELSE

         H(4)=0.25*RP*SM
         H(3)=0.25*RM*SM
         H(2)=0.25*RM*SP
         H(1)=0.25*RP*SP

      ENDIF

!       DEBLJ=1.D0
!       IF (INDAX.EQ.1) THEN
!       DEBLJ=0.D0
!       DO I=1,NDIMEL
!         DEBLJ=DEBLJ+H(I)*CK(1,I)
!       ENDDO
!       ENDIF

      IF (NDIMEL.GT.4) THEN
         IF (NDIMEL.EQ.8) THEN
           ZVHR(9)=0.D0
           ZVHS(9)=0.D0
         ELSE
           ZVHR(9)=-2.D0*R*SS
           ZVHS(9)=-2.D0*S*RR
         ENDIF
         ZVHR(8)=0.5*SS-0.5*ZVHR(9)
         ZVHS(8)=-RP*S-0.5*ZVHS(9)
         ZVHR(7)=-R*SM-0.5*ZVHR(9)
         ZVHS(7)=-0.5*RR-0.5*ZVHS(9)
         ZVHR(6)=-0.5*SS-0.5*ZVHR(9)
         ZVHS(6)=-RM*S-0.5*ZVHS(9)
         ZVHR(5)=-R*SP-0.5*ZVHR(9)
         ZVHS(5)=0.5*RR-0.5*ZVHS(9)
         ZVHR(4)=0.25*SM-0.5*(ZVHR(7)+ZVHR(8))-0.25*ZVHR(9)
         ZVHS(4)=-0.25*RP-0.5*(ZVHS(7)+ZVHS(8))-0.25*ZVHS(9)
         ZVHR(3)=-0.25*SM-0.5*(ZVHR(6)+ZVHR(7))-0.25*ZVHR(9)
         ZVHS(3)=-0.25*RM-0.5*(ZVHS(6)+ZVHS(7))-0.25*ZVHS(9)
         ZVHR(2)=-0.25*SP-0.5*(ZVHR(5)+ZVHR(6))-0.25*ZVHR(9)
         ZVHS(2)=0.25*RM-0.5*(ZVHS(5)+ZVHS(6))-0.25*ZVHS(9)
         ZVHR(1)=0.25*SP-0.5*(ZVHR(5)+ZVHR(8))-0.25*ZVHR(9)
         ZVHS(1)=0.25*RP-0.5*(ZVHS(5)+ZVHS(8))-0.25*ZVHS(9)
         

      ELSE

         ZVHR(4)=0.25*SM
         ZVHS(4)=-0.25*RP
         ZVHR(3)=-0.25*SM
         ZVHS(3)=-0.25*RM
         ZVHR(2)=-0.25*SP
         ZVHS(2)=0.25*RM
         ZVHR(1)=0.25*SP
         ZVHS(1)=0.25*RP

      ENDIF

   50 DO 62 I=1,2
      XJ(1,I)=0.0D0
      XJ(2,I)=0.0D0
      DO 60 J=1,NDIMEL
      XJ(1,I)=XJ(1,I)+ZVHR(J)*CKL(I,J)
      XJ(2,I)=XJ(2,I)+ZVHS(J)*CKL(I,J)
   60 CONTINUE
      XJJ(1,I)=XJ(1,I)
      XJJ(2,I)=XJ(2,I)
   62 CONTINUE
      
      IF(KFIX.GT.0) GO TO 70
CS    DETERMINANTA JAKOBIJANA I INVERZNI JAKOBIJAN
      CALL VMINV2(XJ,DETJ)

      IF (DETJ.LT.0.0D0) THEN
      WRITE(*,*)'DETERMINANTA MANJA OD NULE!!!'
      STOP
      ENDIF

      DO 20 I=1,NDIMEL
        ZVHXL(I)=(ZVHR(I)*XJ(1,1)+ZVHS(I)*XJ(2,1))/DETJ
        ZVHYL(I)=(ZVHS(I)*XJ(2,2)+ZVHR(I)*XJ(1,2))/DETJ
   20 CONTINUE
      
      DO I=1,NDIMEL
        ZVHX(I)=ZVHXL(I)*TTE(1,1)+ZVHYL(I)*TTE(2,1)
        ZVHY(I)=ZVHXL(I)*TTE(1,2)+ZVHYL(I)*TTE(2,2)
        ZVHZ(I)=ZVHXL(I)*TTE(1,3)+ZVHYL(I)*TTE(2,3)
      ENDDO

      RETURN

   70 IF (NPARAM.EQ.1) THEN

        IF (NDIMEL.EQ.4) PJ=0.5D0*(FPS1+FPS2)

        FS2=PJ

        IF (INDAX.EQ.1.AND.HH.GT.0.D0) THEN
           FS2=PJ/(2.D0*RIZ*3.14159D0*HH)
        ENDIF

      ENDIF
CE     LINE JACOBIAN DETERMINANT
      GO TO (71,72),KFIX
CE     CONSTANT KSI
   71 DET2= AJ(2,1)*AJ(2,1)+AJ(2,2)*AJ(2,2)
      GO TO 74
CE     CONSTANT ETA
   72 DET2= AJ(1,1)*AJ(1,1)+AJ(1,2)*AJ(1,2)
   74 DET2=DSQRT(DET2)
      DETJS=DET2
      IF ( DET2.GT.1.D-15) RETURN
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2000) NBREL,KFIX,R,S,DET2
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6000) NBREL,KFIX,R,S,DET2
      STOP
C
 2000 FORMAT(' ** GRESKA **: JAKOBIJAN JEDNAK ILI MANJI OD NULE',
     1       ' ZA ELEMENT No.',I10/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     5       10X,'DET=',F10.5)
 6000 FORMAT(' ** ERROR **: JACOBIAN EQUAL OR LESS THEN ZERO',
     1       ' FOR ELEMENT No.',I10/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     5       10X,'DET=',F10.5)
C

      END
C======================================================================

C======================================================================
      SUBROUTINE VINTER3(R,S,T,KFIX,NPARAM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DET,DET1,NBREL
      COMMON /CVORV/ NOD9(13),NND9
      COMMON /TREN13/ FPS1,FPS2,FPS3,FPS4
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /VISKOZ/ AMI
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /SRPSKI/ ISRPS
      COMMON /NELDIM/ NDIMEL

      DIMENSION P(3,21),PJ(3,21),XJ(3,3),IPERM(8)
C
      DATA IPERM/2,3,4,1,6,7,8,5/
      
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
      IF (NDIMEL.EQ.8) GO TO 50
C
CS    STEPENE SLOBODE ZA CVOROVE PREKO 8
CE     DEGREES OF FREADOM FOR NODES OVER 8
C
      I=0
    2 I=I+1
      IF (I.GT.NND9) GO TO 30
      NN=NOD9(I)-8
      GO TO (9,10,11,12,13,14,15,16,17,18,19,20,21),NN
C
C
    9 H(9)=0.25*RR*SP*TP
      P(1,9)=-0.50*R*SP*TP
      P(2,9)=0.25*RR*TP
      P(3,9)=0.25*RR*SP
      GO TO 2
   10 H(10)=0.25*RM*SS*TP
      P(1,10)=-0.25*SS*TP
      P(2,10)=-0.50*RM*S*TP
      P(3,10)=0.25*RM*SS
      GO TO 2
   11 H(11)=0.25*RR*SM*TP
      P(1,11)=-0.50*R*SM*TP
      P(2,11)=-0.25*RR*TP
      P(3,11)=0.25*RR*SM
      GO TO 2
   12 H(12)=0.25*RP*SS*TP
      P(1,12)=0.25*SS*TP
      P(2,12)=-0.50*RP*S*TP
      P(3,12)=0.25*RP*SS
      GO TO 2
   13 H(13)=0.25*RR*SP*TM
      P(1,13)=-0.50*R*SP*TM
      P(2,13)=0.25*RR*TM
      P(3,13)=-0.25*RR*SP
      GO TO 2
   14 H(14)=0.25*RM*SS*TM
      P(1,14)=-0.25*SS*TM
      P(2,14)=-0.50*RM*S*TM
      P(3,14)=-0.25*RM*SS
      GO TO 2
   15 H(15)=0.25*RR*SM*TM
      P(1,15)=-0.50*R*SM*TM
      P(2,15)=-0.25*RR*TM
      P(3,15)=-0.25*RR*SM
      GO TO 2
   16 H(16)=0.25*RP*SS*TM
      P(1,16)=0.25*SS*TM
      P(2,16)=-0.50*RP*S*TM
      P(3,16)=-0.25*RP*SS
      GO TO 2
   17 H(17)=0.25*RP*SP*TT
      P(1,17)=0.25*SP*TT
      P(2,17)=0.25*RP*TT
      P(3,17)=-0.50*RP*SP*T
      GO TO 2
   18 H(18)=0.25*RM*SP*TT
      P(1,18)=-0.25*SP*TT
      P(2,18)=0.25*RM*TT
      P(3,18)=-0.50*RM*SP*T
      GO TO 2
   19 H(19)=0.25*RM*SM*TT
      P(1,19)=-0.25*SM*TT
      P(2,19)=-0.25*RM*TT
      P(3,19)=-0.50*RM*SM*T
      GO TO 2
   20 H(20)=0.25*RP*SM*TT
      P(1,20)=0.25*SM*TT
      P(2,20)=-0.25*RP*TT
      P(3,20)=-0.50*RP*SM*T
      GO TO 2
   21 H(21)=RR*SS*TT
      P(1,21)=-2.0*R*SS*TT
      P(2,21)=-2.0*S*RR*TT
      P(3,21)=-2.0*T*RR*SS
      GO TO 2
C
CS    KOREKCIJE PRVIH 20 FINKCIJA AKO JE UPOTREBLJEN CVOR 21
CE    CORECTION OF FIRST 20 FUNCTIONS IF NODE 21 EXISTS
C
   30 IN=NOD9(NND9)
      IF(IN.NE.21) GO TO 40
      DO 36 I=1,8
      H(I)=H(I)-0.125*H(21)
      DO 36 J=1,3
   36 P(J,I)=P(J,I)-0.125*P(J,21)
      IF(NND9.EQ.1) GO TO 51
      DO 37 I=1,NND9-1
      IN=NOD9(I)
      H(IN)=H(IN)-0.25*H(21)
      DO 37 J=1,3
   37 P(J,IN)=P(J,IN)-0.25*P(J,21)
C
CS    KOREKCIJE PRVIH 8 FUNKCIJA AKO SU UPOTREBLJENI CVOROVI PREKO 8
CE    CORECTION OF FIRST 8 FUNCTIONS IF NODES OVER 8 EXISTS
C
   40 IH=0
   41 IH=IH+1
      IF(IH.GT.NND9) GO TO 50
      IN=NOD9(IH)
      IF(IN.GT.16) GO TO 46
      I1=IN-8
      I2=IPERM(I1)
C
      H(I1)=H(I1)-0.5*H(IN)
      H(I2)=H(I2)-0.5*H(IN)
      H(IH+8)=H(IN)
      DO 45 J=1,3
      P(J,I1)=P(J,I1)-0.5*P(J,IN)
      P(J,I2)=P(J,I2)-0.5*P(J,IN)
   45 P(J,IH+8)=P(J,IN)
      GO TO 41
C
   46 IF(IN.EQ.21) GO TO 51
      I1=IN-16
      I2=I1+4
      H(I1)=H(I1)-0.5*H(IN)
      H(I2)=H(I2)-0.5*H(IN)
      H(IH+8)=H(IN)
      DO 47 J=1,3
      P(J,I1)=P(J,I1)-0.5*P(J,IN)
      P(J,I2)=P(J,I2)-0.5*P(J,IN)
   47 P(J,IH+8)=P(J,IN)
      GO TO 41
C
   51 H(NND9+8)=H(21)
      DO 39 J=1,3
   39 P(J,NND9+8)=P(J,21)
C
CS    JAKOBIJAN U TACKI R,S,T
CE    JACOBIAN AT POINT R,S,T
C
   50 DO 61 I=1,3
      DO 61 J=1,3
      XJ(I,J)=0.
      DO 61 KK=1,NDIMEL
      XJ(I,J)=XJ(I,J)+P(I,KK)*CK(J,KK)
   61 CONTINUE
C     WRITE(3,8888)((XJ(I,J),J=1,3),I=1,3)
C8888 FORMAT(' X',3D11.3)
      IF(KFIX.GT.0) GO TO 70
C
CS     DETERMINANTA JAKOBIJANA U TACKI R,S,T
CE     JACOBIAN DETERMINANT AT POINT R,S,T
C
C      CALL VMINV(XJ,3,DET1,LE,ME)
      CALL VMINV3(XJ,DET1)
C
C     WRITE(3,7888) ((XJ(I,J),J=1,3),I=1,3)
C7888 FORMAT(' X-',3D11.3)
      IF (DET1.GT.1.D-15) GO TO 77
      WRITE(IIZLAZ,2000) NBREL,KFIX,R,S,T,DET1
      STOP
C
   77 DO 85 I=1,3
      DO 85 JJ=1,NDIMEL
      PJ(I,JJ)=0.
      DO 85 K=1,3
      PJ(I,JJ)=PJ(I,JJ) + XJ(I,K)*P(K,JJ)
   85 CONTINUE
C
      DO 123 J=1,NDIMEL
      ZVHX(J)=PJ(1,J) 
      ZVHY(J)=PJ(2,J) 
 123  ZVHZ(J)=PJ(3,J) 
C
C
      RETURN
C
C
CS     DETERMINATA POVRSINSKOG JAKOBIJANA
CE     SURFACE JACOBIAN DETERMINANT
C
C OVO VAZI SAMO ZA OSMOCVORNE ELEMENTE:
   70  IF (NPARAM.EQ.1) THEN
       IF (NDIMEL.EQ.8) PRJ=0.25D0*(FPS1+FPS2+FPS3+FPS4)
       FS2=PRJ
C
      IF (HH.GT.0.D0.AND.RIZ.GT.0.D0) THEN
        FS2=PRJ/(2.D0*RIZ*3.14159D0*HH)
      ENDIF
C
       ENDIF			       
C
       GO TO (71,72,73),KFIX
CS     KONSTANTNO KSI
CE     CONSTANT KSI
   71 DET=(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))**2+(XJ(3,1)*XJ(2,3)-
     1XJ(3,3)*XJ(2,1))**2+(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))**2
      GO TO 74
CS     KONSTANTNO ETA
CE     CONSTANT ETA
   72 DET=(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))**2+(XJ(1,1)*XJ(3,3)-
     1XJ(1,3)*XJ(3,1))**2+(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))**2
      GO TO 74
CS     KONSTANTNO ZETA
CE     CONSTANT ZETA
   73 DET=(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))**2+(XJ(1,1)*
     1XJ(2,3)-XJ(1,3)*XJ(2,1))**2+(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))**2
   74 DET=DSQRT(DET)
      IF ( DET.GT.1.D-15) RETURN
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2000) NBREL,KFIX,R,S,T,DET
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6000) NBREL,KFIX,R,S,T,DET
      STOP
C
 2000 FORMAT(' ** GRESKA **: JAKOBIJAN JEDNAK ILI MANJI OD NULE',
     1       ' ZA ELEMENT No.',I10/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
 6000 FORMAT(' ** ERROR **: JACOBIAN EQUAL OR LESS THEN ZERO',
     1       ' FOR ELEMENT No.',I10/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
C

      END

C======================================================================

C======================================================================
      SUBROUTINE VOTVIZL
CS OTVARANJE IZLAZNE DATOTEKE
CE OPEN OUTPUT FILE
      COMMON /SRPSKI/ ISRPS
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /IMEULZ/ PAKLST,PAKUNV,ZSILE,PAKNEU,IDUZIN,ZTEMP
      CHARACTER *54 PAKLST,PAKUNV,ZSILE,PAKNEU,ZTEMP
      CHARACTER*1 IME*50,STAT*3
      LOGICAL OLDNEW
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' OTVIZL'
C
C
CS IZLAZNA DATOTEKA
CE OUTPUT FILE
C
    5 CONTINUE
      IF(ISRPS.EQ.1)
     *WRITE(*,*)'   ENTER NAME OF OUTPUT FILE /"',PAKLST(1:IDUZIN),'"'
      IF(ISRPS.EQ.0)
     *WRITE(*,*)'   UNETI IME IZLAZNE DATOTEKE /"',PAKLST(1:IDUZIN),'"'
      READ (*,910) IME
      IF(IME.EQ.'                    ') IME = PAKLST
  910 FORMAT (A)
C
   10 STAT='NEW'
      INQUIRE(FILE=IME,EXIST=OLDNEW)
      IF(OLDNEW) STAT='OLD'
      IF(STAT.EQ.'NEW') THEN
      OPEN (IIZLAZ,FILE=IME,STATUS='NEW',FORM='FORMATTED',
     1 ACCESS='SEQUENTIAL')
                        ELSE
      OPEN (IIZLAZ,FILE=IME,STATUS='OLD',FORM='FORMATTED',
     1 ACCESS='SEQUENTIAL')
                        ENDIF
C
      IND=0
      IF(STAT.EQ.'OLD') CALL VBRIS (IME,IIZLAZ,IND)
      IF(IND.EQ.1)GO TO 10
      IF(IND.EQ.2)GO TO 5
C
      RETURN
      END
C======================================================================
      SUBROUTINE VBRIS (IME,IUN,IND)
CS BRISANJE FILE-A
CE DELETE FILES
      CHARACTER*1 IME(20),CH*1
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' BRIS  '
      IND=2
      IF(ISRPS.EQ.0)
     *WRITE(*,2000) (IME(I),I=1,20)
      IF(ISRPS.EQ.1)
     *WRITE(*,6000) (IME(I),I=1,20)
 2000 FORMAT ('   FILE    ',20A1,' VEC POSTOJI'//
     1 '                 <ENTER>     PREBRISATI'/
     1 '                   "N"       OBICI ')
 6000 FORMAT ('   FILE    ',20A1,' ALREADY EXISTS'//
     1 '                 PRESS "ENTER" TO DELETE OR'/
     1 '                 KEY   "N"     TO BYPASS')
      READ(*,910) CH
  910 FORMAT(A)
      IF(CH.EQ.' ') CH = 'D'
      IF(CH.EQ.'D'.OR.CH.EQ.'d') THEN
                    CLOSE (IUN,STATUS='DELETE')
                    IND=1
                        ELSE
                        CLOSE (IUN,STATUS='KEEP')
                    ENDIF
      RETURN
      END
C======================================================================
      SUBROUTINE VBRISZST (IME,IUN,IND)
CS BRISANJE FILE-A
CE DELETE FILES
      CHARACTER*1 IME(20),CH*1
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' BRIS  '
      IND=2
!       IF(ISRPS.EQ.0)
!      *WRITE(*,2000) (IME(I),I=1,20)
!       IF(ISRPS.EQ.1)
!      *WRITE(*,6000) (IME(I),I=1,20)
!  2000 FORMAT ('   FILE    ',20A1,' VEC POSTOJI'//
!      1 '                 <ENTER>     PREBRISATI'/
!      1 '                   "N"       OBICI ')
!  6000 FORMAT ('   FILE    ',20A1,' ALREADY EXISTS'//
!      1 '                 PRESS "ENTER" TO DELETE OR'/
!      1 '                 KEY   "N"     TO BYPASS')
!       READ(*,910) CH
      CH=' '
  910 FORMAT(A)
      IF(CH.EQ.' ') CH = 'D'
      IF(CH.EQ.'D'.OR.CH.EQ.'d') THEN
                    CLOSE (IUN,STATUS='DELETE')
                    IND=1
                        ELSE
                        CLOSE (IUN,STATUS='KEEP')
                    ENDIF
      RETURN
      END
C======================================================================
C======================================================================
      SUBROUTINE VWRR10(A,N,CHAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /CDEBUG/ IDEBUG
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE REAL VECTOR IN OUTPUT FILE
CS.    P R O G R A M
CS        ZA ZAPISIVANJE REALNOG VEKTORA U IZLAZNI FILE
C .
C ......................................................................
C
      CHARACTER*10 CHAR
      DIMENSION A(*)
      IF(IDEBUG.GT.0) PRINT *, ' WRR   '
C
      WRITE(IIZLAZ,5010) CHAR
      WRITE(IIZLAZ,5000) (A(I),I=1,N)
      RETURN
C
 5010 FORMAT(A10)
 5000 FORMAT(10(1PD12.4))
      END
C=======================================================================
C======================================================================
      SUBROUTINE VWRR(A,N,CHAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /CDEBUG/ IDEBUG
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE REAL VECTOR IN OUTPUT FILE
CS.    P R O G R A M
CS        ZA ZAPISIVANJE REALNOG VEKTORA U IZLAZNI FILE
C .
C ......................................................................
C
      CHARACTER*4 CHAR
      DIMENSION A(*)
      IF(IDEBUG.GT.0) PRINT *, ' WRR   '
C
      WRITE(IIZLAZ,5010) CHAR
      WRITE(IIZLAZ,5000) (A(I),I=1,N)
      RETURN
C
 5010 FORMAT(A5)
 5000 FORMAT(4(1PD18.9))
      END
C=======================================================================
      SUBROUTINE VIWRR(M,N,CHAR)
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE INTEGER VECTOR IN OUTPUT FILE
CS.    P R O G R A M
CS        ZA ZAPISIVANJE CELOBROJNOG VEKTORA U IZLAZNI FILE
C .
C ......................................................................
C
C
      CHARACTER*4 CHAR
      DIMENSION M(*)
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' IWRR  '
C
      WRITE(IIZLAZ,5010) CHAR
 
      WRITE(IIZLAZ,5000) (M(I),I=1,N)
      RETURN
C
 5010 FORMAT(A5)
 5000 FORMAT(14I5)
      END
C=======================================================================
      SUBROUTINE VULAZT1(ISNUMER)
      USE NODES
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PRIKAZ/ INDSC,IZPOT
      COMMON /SRPSKI/ ISRPS
      COMMON /ICITANJE/INPT
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      CHARACTER*250 ACOZ
!       DIMENSION ID(1,*),CORD(3,*)
      JEDN=0
      ICVOR=0
      NPA=0
      NPI=0
      write(*,*) 'ulazt1,ISNUMER', ISNUMER
      DO 10 I=1,NPT
      CALL VISPITA(ACOZ)
!       IF(INPT.EQ.1)THEN
!        READ(ACOZ,1017) N,NJEDN,(CORD(J,N),J=1,3),KORC
!       ELSE
!        READ(ACOZ,1007) N,NJEDN,(CORD(J,N),J=1,3),KORC
!       ENDIF
      IF(INPT.EQ.1)THEN
       READ(ACOZ,1017) N,NJEDN,(CORD(J,I),J=1,3),KORC
      ELSE
       READ(ACOZ,1007) N,NJEDN,(CORD(J,I),J=1,3),KORC
      ENDIF
!       IF(N.GT.NPT) STOP 'NE RADI SLOBODNA NUMERACIJA (N.GT.NPT)'
      IF((N.GT.NPT).AND.(ISNUMER.EQ.0)) 
     * STOP 'ISNUMER=0 (N.GT.NPT)'
C==========================================================================
C SAMO PRIVREMENO ZBOG GLEDANJA SLOBODNE POVRSINE U PRAVCU OSE Y
C MORA SE PROMENITI OSA Z SA OSOM Y
C      READ(ACOZ,1007) N,NJEDN,X,Y,Z,KORC
C      CORD(1,N)=X
C      CORD(2,N)=Z
C      CORD(3,N)=-Y
C==========================================================================
        IF (NJEDN.EQ.0) THEN
          JEDN=JEDN+1
          ID(1,I)=JEDN
        ELSE
          ID(1,I)=0
        ENDIF
      IF(ISNUMER.EQ.1) THEN 
!       ELSEIF(ISNUMER.EQ.1) THEN 
!         ICVOR=ICVOR+1
        IF(I.EQ.1) THEN
          NPA=N
          NPI=N
        ELSE
          IF(NPA.LT.N) NPA=N
          IF(NPI.GT.N) NPI=N
        ENDIF
        IF (NJEDN.EQ.0) NCVEL(I)=N
!           JEDN=JEDN+1
!           ID(1,I)=JEDN
!           NCVEL(I)=N
!         ELSE
!           ID(1,I)=0
!         ENDIF
      ENDIF
   10 CONTINUE
!       do ii=1,NPT
!            write(IIZLAZ,*) ii,ID(1,ii),NCVEL(ii)
!       enddo
      IF(INDSC.EQ.2.OR.INDSC.EQ.1) RETURN
!       Write(*,*) "indsc ulazt1", INDSC
 1007 FORMAT(I5,3X,I2,3F10.5,I5)
 1017 FORMAT(I10,3X,I2,3F10.5,I5)
C 1008 FORMAT(I5,3X,I2,3E12.5,I5)
 1008 FORMAT(I5,2X,I6,3X,3(F10.3),I5)
 1018 FORMAT(I10,2X,I6,3X,3(F10.3),I5)
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2000)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6000)
 2000 FORMAT(//
     *'CVOR, D.O.F., K O O R D I N A T E   C V O R O V A'/)
 6000 FORMAT(//
     *'NODE, D.O.F., C O O R D I N A T E S   O F   N O D E'/)
      DO 12 I=1,NPT
      IF (ID(1,I).EQ.0) JJ=1
!       IF (ID(1,I).NE.0) JJ=0
      IF (ID(1,I).NE.0) JJ=ID(1,I)
      IF(INPT.EQ.1) THEN
       IF(ISNUMER.EQ.0) THEN 
        WRITE(IIZLAZ,1018) I,JJ,(CORD(J,I),J=1,3),KORC
       ELSE
      WRITE(IIZLAZ,1018) NCVEL(I),ID(1,I),(CORD(J,I),J=1,3),KORC
       ENDIF
      ELSE
       IF(ISNUMER.EQ.0) THEN 
        WRITE(IIZLAZ,1008) I,JJ,(CORD(J,I),J=1,3),KORC
       ELSE
       WRITE(IIZLAZ,1008) NCVEL(I),JJ,(CORD(J,I),J=1,3),KORC
       ENDIF
      ENDIF
   12 CONTINUE
      RETURN
      END
C==========================================================================
      SUBROUTINE VULAZT1D()
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /TIPEL/ NP2DMX
      COMMON /DODAT/ NDIMM
      COMMON /PRIKAZ/ INDSC,IZPOT
      COMMON /SRPSKI/ ISRPS
      CHARACTER*250 ACOZ
C
      if(.not.allocated(thick)) then
        allocate(thick(NET),STAT=istat)
        if(istat.ne.0) stop 'error allocating thick'
      endif
      
      DO 20 J=1,NET
      CALL VISPITA(ACOZ)
      IF(NP2DMX.EQ.2) THEN
       READ(ACOZ,1007) NN,(NEL(I,NN),I=1,2),NEL(NDIMM-1,NN),
     1                          NEL(NDIMM,NN),thick(j)
      ELSE
       READ(ACOZ,1008) NN,(NEL(I,NN),I=1,3),NEL(NDIMM-1,NN),
     1                          NEL(NDIMM,NN),thick(j)
       ENDIF
      IF(NN.GT.NET) STOP 'NE RADI SLOBODNA NUMERACIJA (NN.GT.NET)'
   20 CONTINUE
      IF(INDSC.EQ.2.OR.INDSC.EQ.1) RETURN
!       IF(INDSC.EQ.2) RETURN
 1007 FORMAT(5I5,F10.5)
 1008 FORMAT(6I5,F10.5)
 1009 FORMAT(4I6,F10.5)
 1010 FORMAT(5I6,F10.5)
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,1104)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6104)
 1104 FORMAT(/
     1'ELEMENT,  C  V  O  R  O  V  I   , MATERIJAL, POPRECNI PRESEK'/)
 6104 FORMAT(/
     1'ELEMENT,    N  O  D  E  S    ,MATERIAL NUMBER, CROSS SECTION'/)
      DO 22 J=1,NET
       IF(NP2DMX.EQ.2) THEN
        WRITE(IIZLAZ,1009) J,(NEL(I,J),I=1,2),NEL(NDIMM-1,J),
     1                        thick(j)
       ELSE
        WRITE(IIZLAZ,1009) J,(NEL(I,J),I=1,3),NEL(NDIMM-1,J),
     1                        thick(j)
       ENDIF
   22 CONTINUE
      RETURN
      END
C==========================================================================
      SUBROUTINE VULAZT2()
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /TIPEL/ NP2DMX
      COMMON /DODAT/ NDIMM
      COMMON /PRIKAZ/ INDSC,IZPOT
      COMMON /SRPSKI/ ISRPS
      CHARACTER*250 ACOZ
!       DIMENSION NEL(NDIMM,*)
C
      DO 20 J=1,NET
      CALL VISPITA(ACOZ)
      READ(ACOZ,1008) NN,(NEL(I,NN),I=1,4),NEL(NDIMM-1,NN),NEL(NDIMM,NN)
       IF(NP2DMX.GT.4) THEN
         CALL VISPITA(ACOZ)
         READ(ACOZ,1008) (NEL(I,NN),I=5,NDIMM-2)
       ENDIF
      IF(NN.GT.NET) STOP 'NE RADI SLOBODNA NUMERACIJA (NN.GT.NET)'
   20 CONTINUE
      IF(INDSC.EQ.2.OR.INDSC.EQ.1) RETURN
!       IF(INDSC.EQ.2) RETURN
 1008 FORMAT(12I5)
 1009 FORMAT(24I6)
 1010 FORMAT(6X,24I6)
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,1104)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6104)
 1104 FORMAT(/
     1'ELEMENT,   C   V   O   R   O   V   I     , MATERIJAL'/)
 6104 FORMAT(/
     1'ELEMENT,    N   O   D   E   S     ,MATERIAL NUMBER'/)
      DO 22 J=1,NET
       WRITE(IIZLAZ,1009) J,(NEL(I,J),I=1,4),NEL(NDIMM-1,J),NEL(NDIMM,J)
        IF(NP2DMX.GT.4) WRITE(IIZLAZ,1010) (NEL(I,J),I=5,NDIMM-2)
   22 CONTINUE
      RETURN
      END
C==========================================================================

C==========================================================================
      SUBROUTINE VULAZR3(ISNUMER)
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /TIPEL/ NP2DMX
      COMMON /DODAT/ NDIMM
      COMMON /SRPSKI/ ISRPS
      COMMON /ICITANJE/INPT
      COMMON /PRIKAZ/ INDSC,IZPOT
      CHARACTER*250 ACOZ
!       DIMENSION NEL(NDIMM,*)
C
      if(.not.allocated(thick)) then
        allocate(thick(NET),STAT=istat)
        if(istat.ne.0) stop 'error allocating thick'
      endif
      if(.not.allocated(elemtip)) then
        allocate(elemtip(NET),STAT=istat)
        if(istat.ne.0) stop 'error allocating elemtip'
      endif
      DO i=1,6
       eltypes(i)=0
      ENDDO
      IELEM=0
      NMA=0
      NMI=0
      DO 20 J=1,NET
      IELEM=IELEM+1
      CALL VISPITA(ACOZ)
      IF(INPT.EQ.1) THEN
        READ(ACOZ,1006) NETIP
      ELSE
        READ(ACOZ,1005) NETIP
      ENDIF
!       WRITE(*,*) 'NETIP',NETIP
      if(NETIP.EQ.0) READ(ACOZ,1006) NETIP
      if(NETIP.gt.100) then
        NTMP=NETIP/100
      ELSE
        NTMP=NETIP/10
      ENDIF
      eltypes(NTMP)=NETIP
      NETIP=NTMP
!  nema slobodne numeracije elemenata
      IF(ISNUMER.EQ.0) THEN
      IF(INPT.EQ.1) THEN
c      READ(ACOZ,1118) NN,(NEL(I,NN),I=1,8),NEL(NDIMM-1,NN)
!       READ(ACOZ,1118) NN,(NEL(I,NN),I=1,8),NEL(NDIMM-1,NN),NEL(NDIMM,NN)
      SELECT CASE (NETIP)
       CASE(1)
        READ(ACOZ,1117) NN,(NEL(I,NN),I=1,NDIMM-2),NEL(NDIMM-1,NN),
     1                  NEL(NDIMM,NN),thick(j),elemtip(j)
       CASE(2)
        READ(ACOZ,1117) NN,(NEL(I,NN),I=1,NDIMM-2),NEL(NDIMM-1,NN),
     1                  NEL(NDIMM,NN),thick(j),elemtip(j)
       CASE(3)
        READ(ACOZ,1117) NN,(NEL(I,NN),I=1,NDIMM-2),NEL(NDIMM-1,NN),
     1                  NEL(NDIMM,NN),tmp,elemtip(j)
       CASE DEFAULT
        WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',NETIP
        STOP
      END SELECT
      ELSE
      SELECT CASE (NETIP)
       CASE(1)
!       READ(ACOZ,1008) NN,(NEL(I,NN),I=1,8),NEL(NDIMM-1,NN),NEL(NDIMM,NN)
        READ(ACOZ,1007) NN,(NEL(I,NN),I=1,NDIMM-2),NEL(NDIMM-1,NN),
     1                  NEL(NDIMM,NN),thick(j),elemtip(j)
       CASE(2)
        READ(ACOZ,1007) NN,(NEL(I,NN),I=1,NDIMM-2),NEL(NDIMM-1,NN),
     1                  NEL(NDIMM,NN),thick(j),elemtip(j)
       CASE(3)
        READ(ACOZ,1007) NN,(NEL(I,NN),I=1,NDIMM-2),NEL(NDIMM-1,NN),
     1                  NEL(NDIMM,NN),tmp,elemtip(j)
       CASE DEFAULT
        WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',NETIP
        STOP
       END SELECT
      ENDIF
      IF(NN.GT.NET) STOP 'ISNUMER=0 (NN.GT.NET)'
      ELSEIF(ISNUMER.EQ.1) THEN
       IF(INPT.EQ.1) THEN
c      READ(ACOZ,1118) NN,(NEL(I,NN),I=1,8),NEL(NDIMM-1,NN)
!       READ(ACOZ,1118) NN,(NEL(I,NN),I=1,8),NEL(NDIMM-1,NN),NEL(NDIMM,NN)
      SELECT CASE (NETIP)
       CASE(1)
        READ(ACOZ,1117) NN,(NEL(I,j),I=1,NDIMM-2),NEL(NDIMM-1,j),
     1                  NEL(NDIMM,j),thick(j),elemtip(j)
       CASE(2)
        READ(ACOZ,1117) NN,(NEL(I,j),I=1,NDIMM-2),NEL(NDIMM-1,j),
     1                  NEL(NDIMM,j),thick(j),elemtip(j)
       CASE(3)
        READ(ACOZ,1117) NN,(NEL(I,j),I=1,NDIMM-2),NEL(NDIMM-1,j),
     1                  NEL(NDIMM,j),tmp,elemtip(j)
       CASE DEFAULT
        WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',NETIP
        STOP
      END SELECT
      ELSE
      SELECT CASE (NETIP)
       CASE(1)
!       READ(ACOZ,1008) NN,(NEL(I,NN),I=1,8),NEL(NDIMM-1,NN),NEL(NDIMM,NN)
        READ(ACOZ,1007) NN,(NEL(I,j),I=1,NDIMM-2),NEL(NDIMM-1,j),
     1                  NEL(NDIMM,j),thick(j),elemtip(j)
       CASE(2)
        READ(ACOZ,1007) NN,(NEL(I,j),I=1,NDIMM-2),NEL(NDIMM-1,j),
     1                  NEL(NDIMM,j),thick(j),elemtip(j)
       CASE(3)
        READ(ACOZ,1007) NN,(NEL(I,j),I=1,NDIMM-2),NEL(NDIMM-1,j),
     1                  NEL(NDIMM,j),tmp,elemtip(j)
       CASE DEFAULT
        WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',NETIP
        STOP
       END SELECT
      ENDIF
         MCVEL(J)=NN
!          NI=NN
         IF(J.EQ.1) THEN
            NMA=NN
            NMI=NN
         ELSE
            IF(NMA.LT.NN) NMA=NN
            IF(NMI.GT.NN) NMI=NN
         ENDIF
         NN=J
      ENDIF
!        IF(NP2DMX.GT.8) THEN
!          CALL VISPITA(ACOZ)
!          IF(INPT.EQ.1) THEN
!            READ(ACOZ,1018) (NEL(I,NN),I=9,NDIMM-2)
!          ELSE
!            READ(ACOZ,1008) (NEL(I,NN),I=9,NDIMM-2)
!          ENDIF
!        ENDIF
!       IF(NN.GT.NET) STOP 'NE RADI SLOBODNA NUMERACIJA (NN.GT.NET)'
   20 CONTINUE
      numeltip=0
      do i=1,6
       if(eltypes(i).ne.0) then
       numeltip=numeltip+1
        do j=1,i
        if(i.gt.1.and.eltypes(j).eq.0) then
          eltypes(j)=eltypes(i)
          eltypes(i)=0
          exit
        endif
        enddo
       endif
      enddo
      IF(INDSC.EQ.2.OR.INDSC.EQ.1) RETURN
!       do i=1,6
!         write(*,*) 'eltypes(',i,')',eltypes(i)
!       enddo

 1005 FORMAT(125X,I5)
 1006 FORMAT(230X,I5)
 1007 FORMAT(23I5,E10.3,I5)
 1008 FORMAT(23I5,I5,I10)
 1009 FORMAT(24I6)
 1010 FORMAT(6X,24I6)
 1017 FORMAT(9I10,I5)
 1117 FORMAT(21I10,2I5,E10.3,I5)
 1018 FORMAT(9I10,I5)
 1118 FORMAT(9I10,2I5)
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,1104)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6104)
 1104 FORMAT(/
     111X,'ELEMENT,   C   V   O   R   O   V   I     , MATERIJAL,
     2POPRECNI PRESEK,ELTIP'/)
 6104 FORMAT(/
     111X,'ELEMENT,    N   O   D   E   S             ,MATERIAL NUMBER,
     2CROSS SECTION,ELTIP'/)
      DO 22 J=1,NET
       IF(INPT.EQ.1) THEN
       WRITE(IIZLAZ,1117) J,(NEL(I,J),I=1,NDIMM-2),NEL(NDIMM-1,J),
     1                  NEL(NDIMM,J),thick(j),elemtip(j)
       ELSE
       WRITE(IIZLAZ,1007) J,(NEL(I,J),I=1,NDIMM-2),NEL(NDIMM-1,J),
     1                  NEL(NDIMM,J),thick(j),elemtip(j)
       ENDIF
   22 CONTINUE
      RETURN
      END
C==========================================================================

C==========================================================================
      SUBROUTINE VULAZT3(NZAD,ZADVRE,NZADJ,ISNUMER)
      USE NODES
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /SRPSKI/ ISRPS
      COMMON /ICITANJE/INPT
      COMMON /PRIKAZ/ INDSC,IZPOT
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      CHARACTER*250 ACOZ
      DIMENSION NZAD(3,*),ZADVRE(*),NZADJ(*)
      INTEGER NN,NI,N
C
      DO 28 I=1,NUMZAD
      CALL VISPITA(ACOZ)
      IF(INPT.EQ.1) THEN
        READ(ACOZ,1022) NN
      ELSE
        READ(ACOZ,1012) NN
      ENDIF
      IF(ISNUMER.EQ.0) THEN
         NI=NN
         IF(NN.LE.0.OR.NN.GT.NPT) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) NN
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) NN
            STOP ' PROGRAM STOP - PAKV2 - ULAZT3'
         ENDIF
      ELSE
         N=NN-NPI+1
         NI=NELCV(N)
         IF(NI.EQ.0.OR.NN.LT.NPI.OR.NN.GT.NPA) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) NN
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) NN
            STOP ' PROGRAM STOP - PAKV2 - ULAZT3'
         ENDIF
      ENDIF
      NZAD(1,I)=NI
      IF(INPT.EQ.1) THEN
        READ(ACOZ,1122) NZAD(2,I),NZAD(3,I),ZADVRE(I)
      ELSE
        READ(ACOZ,1112) NZAD(2,I),NZAD(3,I),ZADVRE(I)
      ENDIF
      NZADJ(I)=ID(1,NZAD(1,I))
  28  CONTINUE
      IF(INDSC.EQ.2.OR.INDSC.EQ.1) RETURN
 1012 FORMAT(I5)
 1022 FORMAT(I10)
 1112 FORMAT(5X,2I5,F10.3)
 1122 FORMAT(10X,2I5,F10.3)
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,1101)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6101)
 1101 FORMAT(10X,
     1'BROJ CVORA, IND. PROMENLJIVE, VREMEN. FUNK., ZAD. VREDNOST'/)
 6101 FORMAT(10X,
     1'NODE NUMB., IND. QUANTITIES , TIME FUNCTION,  V A L U E'/)
      DO 32 I=1,NUMZAD
      WRITE(IIZLAZ,1013) NZAD(1,I),NZAD(2,I),NZAD(3,I),ZADVRE(I)
!       NZADJ(I)=ID(1,NZAD(1,I))
!       write(*,*) 'iz ulazt3 ID(1,NZAD(1,',I,'))',ID(1,NZAD(1,I))
  32  CONTINUE
 1013 FORMAT (I20,I10,I14,10X,F10.3)
 2100 FORMAT(//' GRESKA U UCITAVANJU ZADATIH POT/TEMP'/
     1' CVOR',I10,' NE POSTOJI'//' PROGRAM STOP - PAK051 - FSPOLJ'//)
 6100 FORMAT(//' ERROR IN INPUT DATA ABOUT PRESCRIBE POT/TEM'/
     1' NODE',I10,' NOT EXIST'//' PROGRAM STOP - PAK051 - FSPOLJ'//)
      RETURN
      END
C==========================================================================

C==========================================================================
      SUBROUTINE VULAZAXIS()
      USE PREDISCRIBED
      COMMON /PRIKAZ/ INDSC,IZPOT
      COMMON /ICITANJE/ INPT
      CHARACTER*250 ACOZ
      INTEGER*8 I
C
      if(NUMAXISPTSX.gt.0) then
      DO I=1,NUMAXISPTSX
        CALL VISPITA(ACOZ)
!         IF(INPT.EQ.1) THEN
            READ(ACOZ,1022) INTAXISPOINTX(1,I),
     1              INTAXISPOINTX(2,I),XAXISPTCORD(I)
!         ELSE
!             READ(ACOZ,1012) INTAXISPOINTX(1,I),
!      1              INTAXISPOINTX(2,I),XAXISPTCORD(I)
!         ENDIF
      ENDDO
      endif
      if(NUMAXISPTSY.gt.0) then
      DO I=1,NUMAXISPTSY
        CALL VISPITA(ACOZ)
!         IF(INPT.EQ.1) THEN
            READ(ACOZ,1022) INTAXISPOINTY(1,I),
     1              INTAXISPOINTY(2,I),YAXISPTCORD(I)
!         ELSE
!             READ(ACOZ,1012) INTAXISPOINTY(1,I),
!      1              INTAXISPOINTY(2,I),YAXISPTCORD(I)
!         ENDIF
      ENDDO
      endif
      if(NUMAXISPTSZ.gt.0) then
      DO I=1,NUMAXISPTSZ
        CALL VISPITA(ACOZ)
!         IF(INPT.EQ.1) THEN
            READ(ACOZ,1022) INTAXISPOINTZ(1,I),
     1              INTAXISPOINTZ(2,I),ZAXISPTCORD(I)
!         ELSE
!             READ(ACOZ,1012) INTAXISPOINTZ(1,I),
!      1              INTAXISPOINTZ(2,I),ZAXISPTCORD(I)
!         ENDIF
      ENDDO
      endif
      IF(INDSC.EQ.2.OR.INDSC.EQ.1) RETURN
!       IF(INDSC.EQ.2) RETURN
 1012 FORMAT(2I5,F10.3)
 1022 FORMAT(2I10,F10.3)
      IF(ISRPS.EQ.0)
     *WRITE(3,1101)
      IF(ISRPS.EQ.1)
     *WRITE(3,6101)
 1101 FORMAT(10X,
     1'BROJ TACKE, BROJ VREMEN. FUNK., KOORDINATA'/)
 6101 FORMAT(10X,
     1'POINT NUMB., ID OF TIME FUNCTION,  COORDINATE'/)
      DO I=1,NUMAXISPTSX
      WRITE(3,1013) INTAXISPOINTX(1,I),INTAXISPOINTX(2,I),
     1 XAXISPTCORD(I)
      ENDDO
 1013 FORMAT (2I10,F10.3)
      RETURN
      END
C==========================================================================

C==========================================================================
      SUBROUTINE VULAZKO(AK,NUMMAT,INDPK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PROMEN/ NJUTN,INDOPT,INDPRO
      COMMON /VDP/ DT,NKORP,NN,NZAV
      COMMON /SRPSKI/ ISRPS
      CHARACTER*250 ACOZ
      DIMENSION AK(3,5,*)
      DIMENSION INDPK(*)

      DO I=1,NUMMAT
       CALL VISPITA(ACOZ)
       READ(ACOZ,1014) AK(1,1,I),AK(1,2,I),AK(1,3,I),
C     1AK(2,1,I),AK(2,2,I),AK(2,3,I),AK(3,4,I),N,NAV,INDPK(I)
     1AK(2,1,I),AK(2,2,I),AK(2,3,I),AK(3,4,I),INDPK(I),N,NAV
       IF (AK(2,1,I).GT.1.D-15) THEN
        NJUTN=I
        INDPRO=1
	NN=N
	NZAV=NAV
       ENDIF
       IF (AK(2,2,I).GT.1.D-15) THEN
        NJUTN=I
        INDPRO=2
	NN=N
	NZAV=NAV
       ENDIF
       IF (AK(2,3,I).GT.1.D-15) THEN
        NJUTN=I
        INDPRO=3
	NN=N
	NZAV=NAV
       ENDIF
       AK(3,1,I)=AK(1,1,I)
       AK(3,2,I)=AK(1,2,I)
       AK(3,3,I)=AK(1,3,I)
	IF (N.GE.1) THEN
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,6016)I,AK(1,1,I),I,AK(1,2,I),I,AK(1,3,I),
     1I,AK(2,1,I),I,AK(2,2,I),I,AK(2,3,I),I,AK(3,4,I),NZAV,NN
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6016)I,AK(1,1,I),I,AK(1,2,I),I,AK(1,3,I),
     1I,AK(2,1,I),I,AK(2,2,I),I,AK(2,3,I),I,AK(3,4,I),NZAV,NN
	ELSE
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,6015)I,AK(1,1,I),I,AK(1,2,I),I,AK(1,3,I),
     1I,AK(2,1,I),I,AK(2,2,I),I,AK(2,3,I),I,AK(3,4,I)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6015)I,AK(1,1,I),I,AK(1,2,I),I,AK(1,3,I),
     1I,AK(2,1,I),I,AK(2,2,I),I,AK(2,3,I),I,AK(3,4,I)
	ENDIF
      ENDDO
      RETURN
 1014 FORMAT(7F10.3,3I5)
C
 6015 FORMAT(////
     111X,'PERMEABILITY Kxmin FOR MATERIAL ',I3,' ...Kxmin= ',1PD10.3//
     111X,'PERMEABILITY Kymin FOR MATERIAL ',I3,' ...Kymin= ',1PD10.3//
     111X,'STORAGEmin FOR MATERIAL ',I3,' ...........Smin = ',1PD10.3///
     111X,'PERMEABILITY Kxmax FOR MATERIAL ',I3,' ...Kxmax= ',1PD10.3//
     111X,'PERMEABILITY Kymax FOR MATERIAL ',I3,' ...Kymax= ',1PD10.3//
     111X,'STORAGEmax FOR MATERIAL ',I3,' ...........Smax = ',1PD10.3//
     111X,'KX/KY FOR MATERIAL ',I3,' ................KX/KY= ',1PD10.3//)
C
 6016 FORMAT(////
     111X,'PERMEABILITY Kxmin FOR MATERIAL ',I3,' ...Kxmin= ',1PD10.3//
     111X,'PERMEABILITY Kymin FOR MATERIAL ',I3,' ...Kymin= ',1PD10.3//
     111X,'STORAGEmin FOR MATERIAL ',I3,' ...........Smin = ',1PD10.3///
     111X,'PERMEABILITY Kxmax FOR MATERIAL ',I3,' ...Kxmax= ',1PD10.3//
     111X,'PERMEABILITY Kymax FOR MATERIAL ',I3,' ...Kymax= ',1PD10.3//
     111X,'STORAGEmax FOR MATERIAL ',I3,' ...........Smax = ',1PD10.3//
     111X,'KX/KY FOR MATERIAL ',I3,' ................KX/KY= ',1PD10.3//
     111X,'LINEAR OR LOGARITHM RELATION Kmin-Kmax.....NZAV= ',I5//
     111X,'NUMBER OF POINTS FOR RELATION Kmin-Kmax......NN= ',I5//)
      
      END
C=======================================================================
      FUNCTION VDELTA(I,J)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CRONECKER DELTA
CS.   P R O G R A M
CS.      ZA KRONEKER DELTA
C .
C ......................................................................
C
      DELTA=0.D0
      IF(I.EQ.J) DELTA=1.D0
      RETURN
      END
C=======================================================================
      SUBROUTINE VZAGLAV
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT


      IF(IDEBUG.GT.0) PRINT *, ' ZAGLAV'
C
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2050)
      WRITE(IIZLAZ,6100)
      WRITE(IIZLAZ,2103)
      WRITE(IIZLAZ,2115)
      RETURN
C----------------------------------------------------------------------
 2050 FORMAT('1'/////////33X,'P R O G R A M'/
     1 10X,'ZA STACIONARNU I NELINEARNU ANALIZU STRUJANJA'/
     1 19X,'FLUIDA KROZ POROZNE SREDINE'/
     1 26X,'METODOM KONACNIH ELEMENATA')
 6100 FORMAT(/////////////
     1 11X,' PPPPPPPPP         AAAA        KK      KK     PPPPPPPPP'/
     2 11X,' PPPPPPPPPP       AAAAAA       KK     KK      PPPPPPPPPP'/
     3 11X,' PP      PP      AA    AA      KK    KK       PP      PP'/
     4 11X,' PP      PP     AA      AA     KK   KK        PP      PP'/
     5 11X,' PPPPPPPPPP     AA      AA     KK  KK         PPPPPPPPPP'/
     6 11X,' PPPPPPPPP      AA      AA     KKKKK          PPPPPPPPP'/
     7 11X,' PP             AAAAAAAAAA     KKKKKKK        PP'/
     8 11X,' PP             AAAAAAAAAA     KK    KK       PP'/
     1 11X,' PP             AA      AA     KK     KK      PP'/
     2 11X,' PP             AA      AA     KK      KK     PP')
 2103 FORMAT(//////3X,'---------------------------  VERZIJA  1.00  ',
     1                '---------------------------')
 2115 FORMAT(//
     1 30X,'MASINSKI FAKULTET'/
     1 30X,'UL. S. JANJICA 6'/
     1 30X,'34000 KRAGUJEVAC'/
     1 30X,'SRBIJA-JUGOSLAVIJA'///)
       END
C----------------------------------------------------------------------
C=======================================================================
      SUBROUTINE VMAXAPR(A,V,R,MAXA,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.  P R O G R A M
CE.     TO MULTIPLY MATRIX STORED BY COLUMNS AND VECTOR
CS.   P R O G R A M
CS.      ZA MNOZENJE MATRICE SLOZENE PO STUPCIMA I VEKTORA
C .
CE.       A - MATRIX STORE BY COLUMNS
CE.       V - VECTOR
CE.       R - R=R+A*V
CE.    MAXA - ADDRESSES OF DIAGONAL ELEMENTS IN MATRIX A
CE.       N - NUMBER OF ROW IN A, ALSO NUMBER OF MEMBERS IN VECTORS V
CE.           AND R
CS.       A - MATRICA UREDJENA PO STUPCIMA
CS.       V - VEKTOR KOJI SE MNOZI MATRICOM
CS.       R - PROIZVOD R=R+A*V
CS.    MAXA - VEKTOR ADRESA DIJAGONALNIH CLANOVA MATRICE A (DIMENZ. N+1)
CS.       N - BROJ VRSTA U MATRICI A, TJ. CLANOVA U VEKTORIMA V I R
C .
C ......................................................................
C
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),MAXA(*),V(*),R(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' MAXAPR'
      DO 10 IV =1,N
C
CE       RIGHT FROM DIAGONAL
CS       DESNO OD DIJAGONALE
C
C         WRITE(IIZLAZ,*)'IV=',IV
         DO 1 I=IV,N
         IF (I.EQ.1) THEN
            II=1
            R(IV)=R(IV)+A(1)*V(1)
            GOTO 1
         ENDIF
C            IMH = MAXA(I+1)-MAXA(I)-1
            IMH = MAXA(I)-MAXA(I-1)-1
            IF(IMH.GT.I) GO TO 1
            II=MAXA(I)-(I-IV)
            IF(II.LE.MAXA(I-1)) GO TO 1
            R(IV)=R(IV)+A(II)*V(I)
    1    CONTINUE
C
CE       LEFT FROM DIAGONAL
CS       LEVO OD DIJAGONALE
C
         IF(IV.EQ.1) GO TO 10
         IMH = MAXA(IV) - MAXA(IV-1) - 1
         IF(IMH.EQ.0) GO TO 10
         IP=IV-IMH
         IK=IV-1
         DO 2 I=IP,IK
            II=MAXA(IV-1)+1+(I-IP)
   2    R(IV)=R(IV)+A(II)*V(I)
   10 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VAXBV( A, B, C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS  VEKTORSKI PROIZVOD
C        C = A X B
C
      DIMENSION A(*),B(*),C(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' AXBV  '
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      RETURN
      END
C======================================================================
      SUBROUTINE VINTERST(R,S,T)
   
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      COMMON /TRENT3/ ZVHX(21),ZVHY(21),ZVHZ(21),CK(3,21),H(21),
     1 FS2,DET,DET1,NBREL
      COMMON /CVORV/ NOD9(13),NND9
      COMMON /TREN13/ FPS1,FPS2,FPS3,FPS4
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /VISKOZ/ AMI
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /SRPSKI/ ISRPS
      COMMON /PIJEZO/ CPOR(3,1000),NPIJEZ(20,100),NODP(100),NPIJ,NPOR,
     1                NPORCV(1000),NPOREL(1000),NEPOR,NPORE,NEPIJ

      DIMENSION P(3,21),PJ(3,21),XJ(3,3),IPERM(8)
C
      DATA IPERM/2,3,4,1,6,7,8,5/
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
C     INTERPOLACIJSKE FUNKCIJE 
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
c      XT=0.
c      YT=0.
c      ZT=0.
c      DO 86 JJ=1,NDIM
c          XT=XT + H(JJ)*CPOR(1,JJ)
c          YT=YT + H(JJ)*CPOR(2,JJ)
c          ZT=ZT + H(JJ)*CPOR(3,JJ)
c   86 CONTINUE
C
C
      END

C======================================================================

C======================================================================
      SUBROUTINE VINTERP(R,S,KFIX,NPARAM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      COMMON /TRENUT/ TT21,H,HP,ZVHX,ZVHY,HV2,HV3,ZVXT,ZVYT,DETJ,
     1DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,B,BT,GRADH,GRADHT,TT210,
     2DEBLJ,ZHP,H2N,H2NP,GH2NT,H2NN
      COMMON /TREN1/ FPS1,FPS2,FP1,FP2,VNN,PERM,CORDY
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /VISKOZ/ AMI
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP

      DIMENSION ZVHR(9),ZVHS(9)
      DIMENSION H(9),ZVHX(9),ZVHY(9)
      DIMENSION HP(4),X(9),Y(9)
      DIMENSION TT21(44),TT210(44),H2N(2,18),H2NP(2,4)
      DIMENSION B(4,18),BT(18,4),GRADH(2,4),GRADHT(4,2),ZHP(4,2)
      DIMENSION GH2NT(18,2)
      DIMENSION H2NN(2,18)
      DIMENSION PERM(2,2),CORDY(4)

      DIMENSION AJ(2,2)

      RP=1.D0+R
      SP=1.D0+S
      RM=1.D0-R
      SM=1.D0-S
      RR=1.D0-R*R
      SS=1.D0-S*S

      IF (NDIM.GT.4) THEN
         IF (NDIM.EQ.8) THEN
           H(9)=0.D0
         ELSE
           H(9)=RR*SS
         ENDIF 
         H(8)=0.5*RP*SS-0.5*H(9)
         H(7)=0.5*RR*SM-0.5*H(9)
         H(6)=0.5*RM*SS-0.5*H(9)
         H(5)=0.5*RR*SP-0.5*H(9)
         H(4)=0.25*RP*SM-0.5*(H(7)+H(8))-0.25*H(9)
         H(3)=0.25*RM*SM-0.5*(H(6)+H(7))-0.25*H(9)
         H(2)=0.25*RM*SP-0.5*(H(5)+H(6))-0.25*H(9)
         H(1)=0.25*RP*SP-0.5*(H(5)+H(8))-0.25*H(9)

      ELSE

         H(4)=0.25*RP*SM
         H(3)=0.25*RM*SM
         H(2)=0.25*RM*SP
         H(1)=0.25*RP*SP

      ENDIF

      DEBLJ=1.D0
      IF (INDAX.EQ.1) THEN
      DEBLJ=0.D0
      DO I=1,NDIM
        DEBLJ=DEBLJ+H(I)*X(I)
      ENDDO
      ENDIF

      IF (NDIM.GT.4) THEN
         IF (NDIM.EQ.8) THEN
           ZVHR(9)=0.D0
           ZVHS(9)=0.D0
         ELSE
           ZVHR(9)=-2.D0*R*SS
           ZVHS(9)=-2.D0*S*RR
         ENDIF
         ZVHR(8)=0.5*SS-0.5*ZVHR(9)
         ZVHS(8)=-RP*S-0.5*ZVHS(9)
         ZVHR(7)=-R*SM-0.5*ZVHR(9)
         ZVHS(7)=-0.5*RR-0.5*ZVHS(9)
         ZVHR(6)=-0.5*SS-0.5*ZVHR(9)
         ZVHS(6)=-RM*S-0.5*ZVHS(9)
         ZVHR(5)=-R*SP-0.5*ZVHR(9)
         ZVHS(5)=0.5*RR-0.5*ZVHS(9)
         ZVHR(4)=0.25*SM-0.5*(ZVHR(7)+ZVHR(8))-0.25*ZVHR(9)
         ZVHS(4)=-0.25*RP-0.5*(ZVHS(7)+ZVHS(8))-0.25*ZVHS(9)
         ZVHR(3)=-0.25*SM-0.5*(ZVHR(6)+ZVHR(7))-0.25*ZVHR(9)
         ZVHS(3)=-0.25*RM-0.5*(ZVHS(6)+ZVHS(7))-0.25*ZVHS(9)
         ZVHR(2)=-0.25*SP-0.5*(ZVHR(5)+ZVHR(6))-0.25*ZVHR(9)
         ZVHS(2)=0.25*RM-0.5*(ZVHS(5)+ZVHS(6))-0.25*ZVHS(9)
         ZVHR(1)=0.25*SP-0.5*(ZVHR(5)+ZVHR(8))-0.25*ZVHR(9)
         ZVHS(1)=0.25*RP-0.5*(ZVHS(5)+ZVHS(8))-0.25*ZVHS(9)
         

      ELSE

         ZVHR(4)=0.25*SM
         ZVHS(4)=-0.25*RP
         ZVHR(3)=-0.25*SM
         ZVHS(3)=-0.25*RM
         ZVHR(2)=-0.25*SP
         ZVHS(2)=0.25*RM
         ZVHR(1)=0.25*SP
         ZVHS(1)=0.25*RP

      ENDIF

      AJ(1,1)=DOT(ZVHR,X,NDIM)
      AJ(1,2)=DOT(ZVHR,Y,NDIM)
      AJ(2,1)=DOT(ZVHS,X,NDIM)
      AJ(2,2)=DOT(ZVHS,Y,NDIM)
      IF(KFIX.GT.0) GO TO 70

      DETJ=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)

      IF (DETJ.LT.0.0D0) THEN
      WRITE(*,*)'DETERMINANTA MANJA OD NULE!!!'
      STOP
      ENDIF

      DO 20 I=1,NDIM
      ZVHX(I)=(ZVHR(I)*AJ(2,2)-ZVHS(I)*AJ(1,2))/DETJ
      ZVHY(I)=(ZVHS(I)*AJ(1,1)-ZVHR(I)*AJ(2,1))/DETJ
   20 CONTINUE  

   70 IF (NPARAM.EQ.1) THEN

        IF (NDIM.EQ.4) PJ=0.5D0*(FPS1+FPS2)

        FS2=PJ

        IF (INDAX.EQ.1.AND.HH.GT.0.D0) THEN
           FS2=PJ/(2.D0*RIZ*3.14159D0*HH)
        ENDIF

      ENDIF
CE     LINE JACOBIAN DETERMINANT
      GO TO (71,72),KFIX
CE     CONSTANT KSI
   71 DET2= AJ(2,1)*AJ(2,1)+AJ(2,2)*AJ(2,2)
      GO TO 74
CE     CONSTANT ETA
   72 DET2= AJ(1,1)*AJ(1,1)+AJ(1,2)*AJ(1,2)
   74 DET2=DSQRT(DET2)
      DETJS=DET2
      IF ( DET2.GT.1.D-15) RETURN
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2000) NBREL,KFIX,R,S,DET2
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6000) NBREL,KFIX,R,S,DET2
      STOP
C
 2000 FORMAT(' ** GRESKA **: JAKOBIJAN JEDNAK ILI MANJI OD NULE',
     1       ' ZA ELEMENT No.',I10/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     5       10X,'DET=',F10.5)
 6000 FORMAT(' ** ERROR **: JACOBIAN EQUAL OR LESS THEN ZERO',
     1       ' FOR ELEMENT No.',I10/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     5       10X,'DET=',F10.5)
C

      END
C======================================================================
C=======================================================================

      SUBROUTINE VJEDV(A,B,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' JEDV  '
C
CS  IZRACUNAVANJE JEDINICNOG VEKTORA
C        D^2 = A^2 + B^2 + C^2
C        A=A/D, B=B/D, C=C/D
C
      D=DSQRT(A*A+B*B+C*C)
      A=A/D
      B=B/D
      C=C/D
      RETURN
      END
C=======================================================================

! C=======================================================================
!       SUBROUTINE VTRAZI (N,INCVOR,isrednji)
!       USE NODE
!       
!       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!       COMMON /CDEBUG/ IDEBUG
!       integer*8 idonji, igornji,isrednji  
!       
!       idonji=1
!       igornji=N
!       DO WHILE (idonji.LE.igornji)
!          isrednji=(igornji+idonji)/2
!          IF(ID_NODE(isrednji).EQ.INCVORE) RETURN
!          IF(INCVORE.LT.ID_NODE(isrednji)) THEN
!             igornji=isrednji-1
!          ELSE
!             idonji=srednji+1
!          ENDIF
!       END DO
!       STOP 'nije nadjen id cvora' 
!       END
! C=======================================================================
C======================================================================
      SUBROUTINE VVEZACV(NP)
      USE NODES
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CONECT NUMBERS NODES IN FREE NUMERATION
CS.   P R O G R A M
CS.        ZA VEZU BROJEVA CVOROVA I NJIHOVOG REDNOG BROJA U VEKTORU
CS.        BROJEVA CVOROVA U SLOBODNOJ NUMERACIJI
C .
CE.    V E C T O R S
CE.               NCVEL(NP) - VECTOR NUMBER OF NODES
CE.        NELCV(NPA-NPI+1) - VECTOR CONECT NUMBER OF NODES
CS.    V E K T O R I
CS.               NCVEL(NP) - VEKTOR BROJEVA CVOROVA
CS.        NELCV(NPA-NPI+1) - VEKTOR VEZE BROJEVA CVOROVA I NJIHOVOG
CS.                           REDNOG BROJA U VEKTORU BROJEVA CVOROVA
C .
C ......................................................................
C
!        DIMENSION NCVEL(*),NELCV(*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, 'PAKV2 - VEZACV'
      DO 10 I=1,NP
         J=NCVEL(I)
         J=J-NPI+1
         NELCV(J)=I
   10 CONTINUE
      RETURN
      END
C======================================================================
C=======================================================================
      SUBROUTINE VVEZAEL(NE)
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CONECT NUMBERS ELEMENTS IN FREE NUMERATION
CS.   P R O G R A M
CS.      ZA POVEZIVANJE BROJEVA ELEMENATA U SLOBODNOJ NUMERACIJI
C .
CE.    V E C T O R S                                                   
CE.               MCVEL(NE) - VECTOR NUMBER OF ELEMENTS       
CE.        MELCV(NMA-NMI+1) - VEKCOR CONECT NUMBER OF ELEMENT
CS.    V E K T O R I                                                   
CS.               MCVEL(NE) - VEKTOR BROJEVA ELEMENATA                
CS.        MELCV(NMA-NMI+1) - VNDIMM,IELEMEKTOR VEZE BROJEVA ELEMENATA I NJIHOVOG  
CS.                           REDNOG BROJA U VEKTORU BROJEVA ELEMENATA  
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
!       DIMENSION MCVEL(*),MELCV(*)
C
      IF(IDEBUG.GT.0) PRINT *, 'PAKV2 - VEZAEL'
      DO 10 I=1,NE
         J=MCVEL(I)
         J=J-NMI+1
         MELCV(J)=I
   10 CONTINUE
      RETURN
      END
C======================================================================
C=======================================================================
      SUBROUTINE VVEZACE(NE,NDIMM)
      USE NODES
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO CONECT NUMBERS NODES ELEMENTS IN FREE NUMERATION
CS.   P R O G R A M
CS.      ZA POVEZIVANJE BROJEVA CVOROVA ELEMENATA U SLOBODNOJ NUMERACIJI
C .
C ......................................................................
C
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
!       DIMENSION NELCV(*),NEL(NE,*)
C
      IF(IDEBUG.GT.0) PRINT *, 'PAKV2 - VEZACE'
!        write(*,*) "NDIMM",NDIMM
!        write(*,*) " EL(1,2)",NEL(1,2)
      DO 20 I=1,NE
         DO 10 J=1,NDIMM
            N=NEL(J,I)
            IF(N.LE.0)GO TO 10
            NN=N
            IF(N.EQ.0) GO TO 10
            IF(N.GE.NPI.AND.N.LE.NPA) THEN
               N=NN-NPI+1
            ENDIF
            K=NELCV(N)
!       write(*,*) "I,J,NEL(J,I),NELCV(N)", I,J, NEL(J,I),NEL(1,2)
            IF(K.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) I,NN
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) I,NN
      STOP ' PROGRAM STOP - PAKV2 - VEZACE'
            ENDIF
            NEL(J,I)=K
   10    CONTINUE
   20 CONTINUE
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//' GRESKA U UCITAVANJU ULAZNIH PODATAKA O ELEMENTIMA'/
     1,I5,' CVOR',I5,' NE POSTOJI'//' PROGRAM STOP - PAKV2 - VEZACE'//)
C-----------------------------------------------------------------------
 6000 FORMAT(//' ERROR IN INPUT DATA ABOUT ELEMENTS'/,I5,
     1' NODE',I5,' NOT EXIST'//' PROGRAM STOP - PAKV2 - VEZACE'//)
C-----------------------------------------------------------------------
      END
C======================================================================
