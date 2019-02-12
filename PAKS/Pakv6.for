C=======================================================================
C
CE    BLOCKS
CS    BLOKOVI
C
C   SUBROUTINE   SPAKUJ
C                SPAKUA
C                RESEN
C                RESENA
C                ZADLEV
C                BLKDEF
C                RESENB
C                CLEARB
C                RBLOCK
C                WBLOCK
C                READDB
C                WRITEB
C
C=======================================================================
      SUBROUTINE VSPAKUJ(SK,MAXA,SKE,LM,ND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO PLACE MATRICES OF ELEMENTS INTO SYSTEM - MAIN
CS.   P R O G R A M
CS.      ZA RAZMESTANJE MATRICA ELEMENATA U SISTEM - GLAVNI
C .
C ......................................................................
C
      include 'paka.inc'
c      COMMON A(17000)
c      REAL A
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION SKE(*),LM(*),MAXA(*),SK(*)
      integer*8 LM
C
      LCMPC=1
      LMPC=1
      ISKDSK=0

      IF(IDEBUG.GT.0) PRINT *, ' SPAKUJ'
      IF(NBLOCK.EQ.1) THEN
CE       WITHOUT BLOCKS
CS       BEZ BLOKOVA
         CALL VSPAKUA(SK,MAXA,SKE,LM,ND,0,
     &               A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,A(LCMPC),A(LMPC))
         IF(ISKDSK.NE.0)
     &   WRITE(ISCRC)ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*(ND+1)/2)
      ELSE
CE       WITH BLOCKS
CS       SA BLOKOVIMA
         WRITE(ISCRC)ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*(ND+1)/2)
      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VSPAKUA(SK,MAXA,SKE,LM,ND,INDD,
     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO PLACE MATRICES OF ELEMENTS INTO SYSTEM - WITHOUT BLOCKS
CS.   P R O G R A M
CS.      ZA RAZMESTANJE MATRICA ELEMENATA U SISTEM - BEZ BLOKOVA
C .
C ......................................................................
C
C      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /CDEBUG/ IDEBUG
      DIMENSION SKE(*),LM(*),MAXA(*),SK(*),MNQ(*),LREC(*),
     &          CMPC(1,*),MPC(1,*)
      integer*8 LM
      
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUA'
C
C
CS  PETLJA PO BLOKOVIMA
      DO 20 KB0=1,NBLOCK
      IF(NBLOCK.EQ.1)THEN
        MNQ0=1
        MNQ1=JEDN
        MXMN=0
        IF(INDD.EQ.1)GO TO 9
        GO TO 11
      ENDIF
C
C
      MNQ0=MNQ(KB0)
      MNQ1=MNQ(KB0+1)-1
      MXMN=MAXA(MNQ(KB0))-1
      LDB=MAXA(MNQ(KB0+1))-MAXA(MNQ(KB0))
      CALL VRBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
CS  PETLJA PO ELEMENTIMA
    9   REWIND ISCRC
   10   READ(ISCRC,END=15,ERR=999)
     &  ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*(ND+1)/2)
C-----------------------------------------------
   11 NDI=0
      DO 200 I=1,ND
      II=LM(I)
C
C
      IF(II.LT.0)THEN
        IIP=-II
        ICM=MPC(1,IIP)
        DO 320 L=1,NEZAV
          II=MPC(L+1,IIP)
          IF(II.LT.MNQ0.OR.(II.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 320
          CMI=CMPC(ICM,L)
          MI=MAXA(II)-MXMN
          KS=I
          DO 310 J=1,ND
            JJ=LM(J)
            IF(JJ)303,310,307
  303       JJP=-JJ
            JCM=MPC(1,JJP)
            KSS=KS
            IF(J.GE.I)KSS=J+NDI
              DO 318 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 318
                IJ=II-JJ
                IF(IJ)318,314,314
  314           KK=MI+IJ
                CMJ=CMPC(JCM,K)
                SK(KK)=SK(KK)+CMI*CMJ*SKE(KSS)
  318         CONTINUE
              GO TO 310
C
C
  307         IJ=II-JJ
              IF(IJ)310,311,311
  311         KK=MI+IJ
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
              SK(KK)=SK(KK)+CMI*SKE(KSS)
  310         KS=KS+ND-J
  320   CONTINUE
        GO TO 200
      ENDIF
      IF(II.LT.MNQ0.OR.(II.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 200
      MI=MAXA(II)-MXMN
      KS=I
      DO 220 J=1,ND
      JJ=LM(J)
      IF(JJ)420,220,110
C
C
  420       JJP=-JJ
            JCM=MPC(1,JJP)
            KSS=KS
            IF(J.GE.I)KSS=J+NDI
              DO 418 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 418
                CMJ=CMPC(JCM,K)
C
                IJ=II-JJ
                IF(IJ)418,415,415
  415           KK=MI+IJ
                SK(KK)=SK(KK)+CMJ*SKE(KSS)
  418           CONTINUE
                GO TO 220
C
C
  110 IJ=II-JJ
      IF(IJ)220,210,210
  210 KK=MI+IJ
      KSS=KS
      IF(J.GE.I)KSS=J+NDI
      SK(KK)=SK(KK)+SKE(KSS)
  220 KS=KS+ND-J
  200 NDI=NDI+ND-I
C
      IF(INDD.EQ.1.OR.NBLOCK.GT.1) GO TO 10
   15 IF(NBLOCK.GT.1) CALL VWBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
C
   20 CONTINUE
      RETURN
999   PRINT *,'ERROR: reading element stifness matrix from disk'
      STOP
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VRESEN(B,V,MAXA,NN,KKK)
      USE STIFFNESS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO SOLUTION OF SYSTEM OF EQUATIONS
CS.   P R O G R A M
CS.      ZA RESAVANJE SISTEMA JEDNACINA
C .
C ......................................................................
C
      include 'paka.inc'
      INCLUDE 'mpif.h'
c      COMMON A(17000)
c      REAL A
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
C      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILIISK,ILISE,
C     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /SOLVER/ ISOLVER


      DIMENSION B(*),V(*),MAXA(*)
      integer ierr, myid
      integer*8 i

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      IF (myid.ne.0) goto 10
C

!         CALL Vwrr(B,nwk,'LEV')
!         CALL Vwrr(V,nn,'DES')

      IF(IDEBUG.GT.0) PRINT *, ' RESEN'
      IF(ISRPS.EQ.0)
     1WRITE(*,2000) KKK
      IF(ISRPS.EQ.1)
     1WRITE(*,6000) KKK
C
  10  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ISOLVER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      IF(IABS(ISOLVER).GT.O) THEN
         IF(ISOLVER.EQ.1) THEN
C        CALL VICCGMA(B,V,MAXA,NWK,NN,KKK,IIZLAZ,TOLG,ALFAG)
         ELSE
c            if(imumps.eq.1) then
c                 if(kkk.eq.2) then
               CALL dmumps1(rows,columns,stiff,V,nonzeros,stiff_n,kkk)
c                    IF (myid.ne.0) goto 20
c                endif
c            endif
         ENDIF
      ELSE
        if(myid.eq.0) then
        IF(NBLOCK.EQ.1) THEN
CE       WITHOUT BLOCKS
CS       BEZ BLOKOVA
C        IF(JPS.GT.1.AND.JPBR.LT.JPS1) THEN
C           CALL VRESENP(B,V,MAXA,NN,JEDNP,IIZLAZ,KKK)
C        ELSE
            CALL VRESENA(B,V,MAXA,NN,IIZLAZ,KKK)
C            CALL VRESEN_(B,V,MAXA,NN,NWK,KKK)
C        ENDIF
        ELSE
CE       WITH BLOCKS
CS       SA BLOKOVIMA
         CALL VRESENB(B,V,MAXA,A(LMNQ),A(LICPL),A(LLREC),NN,KC,
     1   NBLOCK,LR,IBLK,IIZLAZ,KKK)
        ENDIF
        endif
C........................
      ENDIF
!       CALL Vwrr(V,nn,'rsh posle resavanja')
C        CALL Vwrr(B,nn,'sk-kraj')
C........................
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(' *** RESAVANJE SISTEMA JEDNACINA (',I1,') ***')
C----------------------------------------------------------------------
 6000 FORMAT(' *** SOLUTION EQUATIONS SYSTEM (',I1,') ***')
C----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VRESENP(A,V,MAXA,NN,JEDNP,IIZLAZ,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE------- RESEN FOR ONE PART SYSTEM -------------
CS------- RESEN ZA JEDNODELNI SITEM -------------
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /RSN   / DETER,IPIVOT,IDETER
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /CDEBUG/ IDEBUG
      COMMON /SRPSKI/ ISRPS
      DIMENSION A(*),V(*),MAXA(*)
C
C OVDE DODATO:
      IPIVOT=1
      IDETER=0

      IF(IDEBUG.GT.0) PRINT *, ' RESENP'
      IF(JEDNP.EQ.0) RETURN
      IF(KKK-2)40,150,950
  40  DETER=1.D0
      IF(IPIVOT.EQ.-1)IPIVOT=1
      AMAXK=A(1)
      AMINK=A(1)
      AMAXF=A(1)
      AMINF=A(1)
C
C     FAKTORIZACIJA BLOKA 1 PREKO BLOKA 1
C
      DO 140 N=1,JEDNP
        KN=MAXA(N)
        KL=KN+1
        KU=MAXA(N+1)-1
        KH=KU-KL
CZILE
        IF(DABS(A(KN)).LT.1.D-15) THEN
      IF(ISRPS.EQ.0)
     1PRINT *, ' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0'
      IF(ISRPS.EQ.1)
     1PRINT *, ' DIAGONAL STIFFNESS MATRIX = 0'
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2001)N
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6001)N
          A(KN)=1.D0
          DO 10 IO=KL,KU
   10     A(IO)=0.D0
        ELSE
          IF(A(KN).GT.AMAXK) AMAXK=A(KN)
          IF(A(KN).LT.AMINK) AMINK=A(KN)
        ENDIF
CZILE
      IF(KH)110,90,50
  50    K=N-KH
        IC=0
        KLT=KU
        DO 80 J=1,KH
        IC=IC+1
        KLT=KLT-1
        KI=MAXA(K)
        ND=MAXA(K+1)-KI-1
        IF(ND)80,80,60
  60    KK=MIN0(IC,ND)
        C=0.D0
        DO 70 L=1,KK
  70    C=C+A(KI+L)*A(KLT+L)
        A(KLT)=A(KLT)-C
  80    K=K+1
  90    K=N
        B=0.D0
        DO 100 KK=KL,KU
        K=K-1
        KI=MAXA(K)
        C=A(KK)/A(KI)
        B=B+C*A(KK)
 100    A(KK)=C
        A(KN)=A(KN)-B
C
      IF(A(KN).GT.AMAXF) AMAXF=A(KN)
      IF(A(KN).LT.AMINF) AMINF=A(KN)
C
 110  IF(A(KN))115,120,125
 115  IF(IPIVOT.EQ.1)IPIVOT=-1
      IF(IPIVOT)125,120,125
 120    CONTINUE
      PRINT *, ' PIVOT < 0'
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2000)N,A(KN)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6000)N,A(KN)
      IF(KOR.GT.0.OR.ITER.GT.0) GO TO 125
        STOP
 125  IF(IDETER.EQ.0) GO TO 140
C     DETER=DETER*A(KN)
 140   CONTINUE
C
C     FAKTORIZACIJA BLOKA 2A PREKO BLOKA 1  I BLOKA 2B
C
      DO 540 N=JEDNP+1,NN
        KN=MAXA(N)
        KLA=KN+N-JEDNP
        KLB=KN+1
        KU=MAXA(N+1)-1
        KUB=KLB+N-JEDNP-1
        KHA=KU -KLA
        KHB=KUB-KLB
        KH =KU -KLB
C....  BLOK 2A
      IF(KHA)540,545,450
 450    K=JEDNP+1-KHA
        IC=0
        KLT=KU
        DO 480 J=1,KHA
        IC=IC+1
        KLT=KLT-1
        KI=MAXA(K)
        ND=MAXA(K+1)-KI-1
        IF(ND)480,480,460
 460    KK=MIN0(IC,ND)
        C=0.D0
        DO 470 L=1,KK
 470    C=C+A(KI+L)*A(KLT+L)
        A(KLT)=A(KLT)-C
 480    K=K+1
C....  BLOK 2B
 545  IF(KHB)490,490,550
 550    K=N-KHB
        KLT=KUB
        KJJ=K-JEDNP-1
        DO 580 J=1,KHB
        KJJ=KJJ+1
        KLT=KLT-1
        KI=MAXA(K)
        ND=MAXA(K+1)-KI-KJJ
        IF(ND)580,580,560
 560    KK=MIN0(ND,KHA+1)
        C=0.D0
        DO 570 L=KJJ,KJJ+KK-1
 570    C=C+A(KI+L)*A(KLT+L)
        A(KLT)=A(KLT)-C
 580    K=K+1
C
 490    K=JEDNP+1
        B=0.D0
        DO 500 KK=KLA,KU
        K=K-1
        KI=MAXA(K)
        C=A(KK)/A(KI)
        B=B+C*A(KK)
 500    A(KK)=C
        A(KN)=A(KN)-B
 540   CONTINUE
       RETURN
CE     REDUCE  FREE  VECTOR
CS     REDUKOVANJE SLOBODNOG VEKTORA
C
C      STANDARDNA REDUKCIJA
C
 150    DO 180 N=1,JEDNP
        KL=MAXA(N)+1
        KU=MAXA(N+1)-1
        IF(KU-KL)180,160,160
 160    K=N
        C=0.
        DO 170 KK=KL,KU
        K=K-1
 170    C=C+A(KK)*V(K)
        V(N)=V(N)-C
 180    CONTINUE
C
C      MODIFIKOVANA REDUKCIJA UPOTREBOM BLOKA 2A
C
        DO 780 N=JEDNP+1,NN
        KL=MAXA(N)+N-JEDNP
        KU=MAXA(N+1)-1
        IF(KU-KL)780,760,760
 760    K=JEDNP+1
        C=0.
        DO 770 KK=KL,KU
        K=K-1
 770    C=C+A(KK)*V(K)
        V(N)=V(N)-C
 780    CONTINUE
CP
        IF(KKK.EQ.2)RETURN
CP
C       ZAMENA UNAZAD
 950    DO 200 N=1,JEDNP
        K=MAXA(N)
 200    V(N)=V(N)/A(K)
CCC***        IF(JEDNP.EQ.1)RETURN
C.....   KOREKCIJA OD POMERANJA GLOBALNIH CVOROVA 
        N=NN
        DO 230 L=JEDNP+1,NN
        KL=MAXA(N)+N-JEDNP
        KU=MAXA(N+1)-1
        IF(KU-KL)230,210,210
 210    K=JEDNP+1
        DO 220 KK=KL,KU
        K=K-1
 220    V(K)=V(K)-A(KK)*V(N)
 230    N=N-1
C
        N=JEDNP
        DO 930 L=2,JEDNP
        KL=MAXA(N)+1
        KU=MAXA(N+1)-1
        IF(KU-KL)930,910,910
 910    K=N
        DO 920 KK=KL,KU
        K=K-1
 920    V(K)=V(K)-A(KK)*V(N)
 930    N=N-1
CZILE
CS.     CISCENJE NUMERICKIH GRESAKA ZA POMERANJA
CE.     CLEANING NUMERICAL ERRORS FOR DISPLACEMENTS
C        CALL VCISTIP(V,JEDNP)
CZILE
        RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'//
     1  ' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I4,//' ','PIVOT=',
     2  D20.12)
 2001 FORMAT(' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0 ZA JEDNACINU',I6)
C-----------------------------------------------------------------------
 6000 FORMAT(//' ','MATRIX OF SYSTEM IS NOT POSITIVE DEFINITE'//
     1  ' ','PIVOT IS NOT POSITIVE FOR EQUATION NO.',I4,//' ','PIVOT=',
     2  D20.12)
 6001 FORMAT(' DIAGONAL STIFFNESS MATRIX = 0 FOR EQUATION',I6)
C-----------------------------------------------------------------------
        END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VRESENA(A,V,MAXA,NN,IIZLAZ,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE       TO L*D*LT FAKTORIZE WITHOUT BLOCKS
CS.   P R O G R A M
CS.      ZA L*D*LT FAKTORIZACIJU BEZ BLOKOVA
C .
C ......................................................................
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /RSN   / DETER,IPIVOT,IDETER
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /CDEBUG/ IDEBUG
      COMMON /SRPSKI/ ISRPS
      DIMENSION A(*),V(*),MAXA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' RESENA'
      IPIVOT=0
      IDETER=0
      NEQ=JEDN
      IF(NN.EQ.0) RETURN
      IF(KKK-2)40,150,150
  40  DETER=1.D0
      IF(IPIVOT.EQ.-1)IPIVOT=1
      AMAXK=A(1)
      AMINK=A(1)
      AMAXF=A(1)
      AMINF=A(1)
      DO 140 N=1,NN
        KN=MAXA(N)
        KL=KN+1
        KU=MAXA(N+1)-1
        KH=KU-KL
CZILE
        IF(DABS(A(KN)).LT.1.D-15.AND.N.LE.NEQ) THEN
      IF(ISRPS.EQ.0)
     1PRINT *, ' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0',N
      IF(ISRPS.EQ.1)
     1PRINT *, ' DIAGONAL STIFFNESS MATRIX = 0',N
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2001)N
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6001)N
          A(KN)=1.D0
          DO 10 IO=KL,KU
   10     A(IO)=0.D0
        ELSE
          IF(A(KN).GT.AMAXK) AMAXK=A(KN)
          IF(A(KN).LT.AMINK) AMINK=A(KN)
        ENDIF
CZILE
      IF(KH)110,90,50
  50    K=N-KH
        IC=0
        KLT=KU
        DO 80 J=1,KH
        IC=IC+1
        KLT=KLT-1
        KI=MAXA(K)
        ND=MAXA(K+1)-KI-1
        IF(ND)80,80,60
  60    KK=MIN0(IC,ND)
        C=0.D0
        DO 70 L=1,KK
  70    C=C+A(KI+L)*A(KLT+L)
        A(KLT)=A(KLT)-C
  80    K=K+1
  90    K=N
        B=0.D0
        DO 100 KK=KL,KU
        K=K-1
        KI=MAXA(K)
        C=A(KK)/A(KI)
        B=B+C*A(KK)
 100    A(KK)=C
        A(KN)=A(KN)-B
C
      IF(A(KN).GT.AMAXF) AMAXF=A(KN)
      IF(A(KN).LT.AMINF) AMINF=A(KN)
C
 110  IF(A(KN))115,120,125
CC 115  IF(IPIVOT.EQ.1.AND.(ICONT.EQ.0.OR.(ICONT.NE.0.AND.KN.LE.NEQ)))
 115  IF(IPIVOT.EQ.1.AND.(ICONT.EQ.0.OR.(ICONT.NE.0.AND.N.LE.NEQ)))
     $   IPIVOT=-1
      IF(IPIVOT)125,116,125
CC 116  IF(ICONT.NE.0.AND.KN.GT.NEQ)GO TO 125
 116  IF(ICONT.NE.0.AND.N.GT.NEQ)GO TO 125
 120    CONTINUE
      PRINT *, ' PIVOT < 0'
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2000)N,A(KN)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6000)N,A(KN)
CFICA      IF(KOR.GT.0.OR.ITER.GT.0) GO TO 125
CZILE        STOP
 125  IF(IDETER.EQ.0) GO TO 140
C     DETER=DETER*A(KN)
 140   CONTINUE
C
C       WRITE(3,*) 'AMAXK,AMINK',AMAXK,AMINK
C       WRITE(3,*) 'AMAXF,AMINF',AMAXF,AMINF
       RETURN
CE     REDUCE  FREE  VECTOR
CS     REDUKOVANJE SLOBODNOG VEKTORA
 150    CONTINUE
CZILE
CS.     CISCENJE NUMERICKIH GRESAKA ZA POMERANJA
CE.     CLEANING NUMERICAL ERRORS FOR DISPLACEMENTS
C        CALL VCISTIN(V,NN)
CZILE
        DO 180 N=1,NN
        KL=MAXA(N)+1
        KU=MAXA(N+1)-1
        IF(KU-KL)180,160,160
 160    K=N
        C=0.D0
        DO 170 KK=KL,KU
        K=K-1
 170    C=C+A(KK)*V(K)
        V(N)=V(N)-C
 180    CONTINUE
C       ZAMENA UNAZAD
        DO 200 N=1,NN
        K=MAXA(N)
CZILE
C        IF(DABS(A(K)-1.D0).LT.1.D-15.AND.N.LE.NEQ) V(N)=0.D0
CZILE
 200    V(N)=V(N)/A(K)
        IF(NN.EQ.1)RETURN
        N=NN
        DO 230 L=2,NN
        KL=MAXA(N)+1
        KU=MAXA(N+1)-1
        IF(KU-KL)230,210,210
 210    K=N
        DO 220 KK=KL,KU
        K=K-1
 220    V(K)=V(K)-A(KK)*V(N)
 230    N=N-1
CZILE
CS.     CISCENJE NUMERICKIH GRESAKA ZA POMERANJA
CE.     CLEANING NUMERICAL ERRORS FOR DISPLACEMENTS
C        CALL VCISTIP(V,NN)
CZILE
        RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'//
     1  ' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I4,//' ','PIVOT=',
     2  D20.12)
 2001 FORMAT(' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0 ZA JEDNACINU',I6)
C-----------------------------------------------------------------------
 6000 FORMAT(//' ','MATRIX OF SYSTEM IS NOT POSITIVE DEFINITE'//
     1  ' ','PIVOT IS NOT POSITIVE FOR EQUATION NO.',I4,//' ','PIVOT=',
     2  D20.12)
 6001 FORMAT(' DIAGONAL STIFFNESS MATRIX = 0 FOR EQUATION',I6)
C-----------------------------------------------------------------------
        END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VRESENB(A,V,MAXA,MNQ,ICPL,LREC,NN,KC,NBLOCK,LR,IBLK,
     1  IIZLAZ,KKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE       TO L*D*LT FAKTORIZE WITH BLOCKS
CS.   P R O G R A M
CS.      ZA L*D*LT FAKTORIZACIJU SA BLOKOVIMA
C .
C ......................................................................
C
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ITERBR/ ITER
      COMMON /RSN   / DETER,IPIVOT,IDETER
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /SRPSKI/ ISRPS
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*),V(*),MAXA(*),MNQ(*),ICPL(*),LREC(*)
C
C================
C OVDE DODATO:
      IPIVOT=1
      IDETER=0
      ICONT=1
      NEQ=JEDN
C================

      IF(IDEBUG.GT.0) PRINT *, ' RESENB'
      KC1=KC+1
      IF(KKK-2)40,150,150
   40 DETER=1.D0
      IF(IPIVOT.EQ.-1)IPIVOT=1
      DO 145 NB=1,NBLOCK
      IF(ISRPS.EQ.0)
     1WRITE(*,2002) NB,NBLOCK
      IF(ISRPS.EQ.1)
     1WRITE(*,6002) NB,NBLOCK
      CALL VRBLOCK(A(1),LREC,NB,LR,MAXA(MNQ(NB+1))-MAXA(MNQ(NB)),IBLK)
      IF(NB.EQ.1)THEN
        AMAXK=A(1)
        AMINK=A(1)
        AMAXF=A(1)
        AMINF=A(1)
      ENDIF
      KB=ICPL(NB)
      KHB=NB-KB+1
      KI0=1
      IF(KHB.GT.1) KI0=KC1
CE    LOOP OVER COUPLED BLOCKS
CS    PETLJA PO SPREGNUTIM BLOKOVIMA
        DO 142 JB=1,KHB
        IF(KB.NE.NB)
     1CALL VRBLOCK(A(KC1),LREC,KB,LR,MAXA(MNQ(KB+1))-MAXA(MNQ(KB)),IBLK)
CE    LOOP OVER EQUATIONS OF BLOCK NB
CS    PETLJA PO JEDNACINAMA BLOKA NB
      DO 140 N=MNQ(NB),MNQ(NB+1)-1
        KN=MAXA(N)-MAXA(MNQ(NB))+1
        KL=KN+1
        KU=MAXA(N+1)-MAXA(MNQ(NB))
        KH=KU-KL
CZILE
        IF(DABS(A(KN)).LT.1.D-15.AND.N.LE.NEQ) THEN
      IF(ISRPS.EQ.0)
     1PRINT *, ' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0 ZA JEDNACINU',N
      IF(ISRPS.EQ.1)
     1PRINT *, ' DIAGONAL STIFFNESS MATRIX = 0 FOR EQUATION',N
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2001)N
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6001)N
          A(KN)=1.D0
          DO 10 IO=KL,KU
   10     A(IO)=0.D0
        ELSE
          IF(A(KN).GT.AMAXK) AMAXK=A(KN)
          IF(A(KN).LT.AMINK) AMINK=A(KN)
        ENDIF
CZILE
      IF(KB.LT.NB) THEN
                   IF(KH) 140,140,50
                   ELSE
                   IF(KH) 110,90,50
                   ENDIF
C
   50 K1=N-KH
CE    AKO JE K1.GE.MNQ(KB+1)  JEDNACINA N NIJE POVEZANA SA BLOKOM KB
CS    AKO JE K1.GE.MNQ(KB+1)  JEDNACINA N NIJE POVEZANA SA BLOKOM KB
         IF(K1-MNQ(KB+1)) 52,140,140
CE    K  JE PRVA JEDNACINA IZ BLOKA KB SA KOJOM JE SPREGNUTA J-NA N
CS    K  JE PRVA JEDNACINA IZ BLOKA KB SA KOJOM JE SPREGNUTA J-NA N
   52    K=MAX0(K1,MNQ(KB))
         IC=K-K1
CE     KLT  ADRESA PRVOG CLANA G (KASNIJE L)  KOJI JE POVEZAN SA  KB
CS     KLT  ADRESA PRVOG CLANA G (KASNIJE L)  KOJI JE POVEZAN SA  KB
       KLT=KU-IC
CE     KH   BROJ CLANOVA POVEZANIH SA  KB
CS     KH   BROJ CLANOVA POVEZANIH SA  KB
          IF(KB.LT.NB) THEN
          KH=MNQ(KB+1)-K
                ELSE
                KI0=1
                KH=N-K
                ENDIF
      DO 80 J=1,KH
      IC=IC+1
      KLT=KLT-1
      KI=MAXA(K)
      ND=MAXA(K+1)-KI-1
      IF(ND)80,80,60
   60 KK=MIN0(IC,ND)
CE    KIB  ADRESSES OF DIAGONAL ELEMENT OF BLOCK KB
CS    KIB  ADRESA DIJAGONALNOG CLANA BLOKA KB
       KIB=KI-MAXA(MNQ(KB))+KI0
      C=0.0D0
      DO 70 L=1,KK
   70 C=C+A(KIB+L)*A(KLT+L)
      A(KLT)=A(KLT)-C
   80 K=K+1
C
   90 IF(KB.LT.NB)GO TO 140
      K=N
      B=0.0D0
      DO 100 KK=KL,KU
      K=K-1
      C=A(KK)/V(K)
      B=B+C*A(KK)
  100 A(KK)=C
      A(KN)=A(KN)-B
C
      IF(A(KN).GT.AMAXF) AMAXF=A(KN)
      IF(A(KN).LT.AMINF) AMINF=A(KN)
C
 110  IF(A(KN))115,120,125
 115  IF(IPIVOT.EQ.1.AND.(ICONT.EQ.0.OR.(ICONT.NE.0.AND.N.LE.NEQ)))
     $   IPIVOT=-1
      IF(IPIVOT)125,116,125
 116  IF(ICONT.NE.0.AND.N.GT.NEQ)GO TO 125
 120  CONTINUE
      PRINT *, ' PIVOT < 0'
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2000)N,A(KN)
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6000)N,A(KN)
      IF(KOR.GT.0.OR.ITER.GT.0) GO TO 125
C      STOP
  125 IF(IDETER.EQ.0) GO TO 139
C     DETER=DETER*A(KN)
  139 V(N)=A(KN)
  140 CONTINUE
  142 KB=KB+1
      CALL VWBLOCK(A(1),LREC,NB,LR,MAXA(MNQ(NB+1))-MAXA(MNQ(NB)),IBLK)
  145 CONTINUE
CE    PAMTI DIJAGONALNE CLANOVE U DRUGOJ POLOVINI VEKTORA  A
CS    PAMTI DIJAGONALNE CLANOVE U DRUGOJ POLOVINI VEKTORA  A
      DO 147 N=1,NN
  147 A(KC1+N)=V(N)
      RETURN
C
CE    REDUKOVANJE SLOBODNOG VEKTORA
CS    REDUKOVANJE SLOBODNOG VEKTORA
  150 DO 185 NB=1,NBLOCK
      CALL VRBLOCK(A(1),LREC,NB,LR,MAXA(MNQ(NB+1))-MAXA(MNQ(NB)),IBLK)
      DO 180 N=MNQ(NB),MNQ(NB+1)-1
      KL=MAXA(N)-MAXA(MNQ(NB))+2
      KU=MAXA(N+1)-MAXA(MNQ(NB))
      IF(KU-KL)180,160,160
  160 K=N
      C=0.0D0
      DO 170 KK=KL,KU
      K=K-1
  170 C=C+A(KK)*V(K)
      V(N)=V(N)-C
  180 CONTINUE
  185 CONTINUE
C
CE    BACKSUBSTITUTION
CS    ZAMENA UNAZAD
C
      DO 200 N=1,NN
  200 V(N)=V(N)/A(KC1+N)
      IF(NN.EQ.1)RETURN
C
      N=NN
      NB=NBLOCK
      DO 235 LB=1,NBLOCK
      CALL VRBLOCK(A(1),LREC,NB,LR,MAXA(MNQ(NB+1))-MAXA(MNQ(NB)),IBLK)
      L0=MNQ(NB)
      IF(NB.EQ.1) L0=2
      DO 230 L=L0,MNQ(NB+1)-1
      KL=MAXA(N)-MAXA(MNQ(NB))+2
      KU=MAXA(N+1)-MAXA(MNQ(NB))
      IF (KU-KL)230,210,210
  210 K=N
      DO 220 KK=KL,KU
      K=K-1
  220 V(K)=V(K)-A(KK)*V(N)
  230 N=N-1
  235 NB=NB-1
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(//' ','MATRICA SISTEMA NIJE POZITIVNO DEFINITNA'
     1//' ','PIVOT NIJE POZITIVAN ZA JEDNACINU BR.',I4,//' ','PIVOT=',
     2D20.12)
 2001 FORMAT(' DIJAGONALNI CLAN MATRICE KRUTOSTI = 0 ZA JEDNACINU',I6)
 2002 FORMAT(' *** *** BLOK MATRICE SISTEMA',I5,' /',I4,' ***')
C-----------------------------------------------------------------------
 6000 FORMAT(//' ',' MATRIX OF SYSTEM IS NOT POSITIVE DEFINITE'//
     1' ','PIVOT IS NOT POSITIVE FOR EQUATION NO.',I4,//' ','PIVOT=',
     2D20.12)
 6001 FORMAT(' DIAGONAL STIFFNESS MATRIX = 0 FOR EQUATION',I6)
 6002 FORMAT(' *** *** BLOCK OF SYSTEM MATRIX',I5,' /',I4,' ***')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VZADLEV(B,MAXA,NZADJ,NZADP,NZAD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO PRESCRIBED DISPLACEM. - UPDATING OF LEFT SIDE WITH BLOCKS
CS.   P R O G R A M
CS.      ZA ZADATA POMERANJA - KOREKCIJA LEVE STRANE SA BLOKOVIMA
C .
C ......................................................................
C
      include 'paka.inc'
c      COMMON A(17000)
c      REAL A
C      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),MAXA(*),NZADJ(*),IA(1),NZAD(3,*)
      EQUIVALENCE (A(1),IA(1))
      MNQ(NN)=IA(LMNQ+NN-1)
C
      IF(IDEBUG.GT.0) PRINT *, ' ZADLEB'
      VBR=1.0D35
      DO 20 NB=1,NBLOCK
      IF(NBLOCK.EQ.1)THEN
        MNQ0=1
        MNQ1=JEDN
        MXMN=0
        GO TO 11
      ENDIF
C
C
      MNQ0=MNQ(NB)
      MNQ1=MNQ(NB+1)-1
      MXMN=MAXA(MNQ(NB))-1
      LDB=MAXA(MNQ(NB+1))-MAXA(MNQ(NB))
      CALL VRBLOCK(B,A(LLREC),NB,LR,LDB,IBLK)
   11   DO 10 I=1,NZADP
          IF (NZAD(2,I).EQ.4.OR.NZAD(2,I).EQ.5) GOTO 10
          NJ=NZADJ(I)
          IF(NJ.LT.MNQ0.OR.(NJ.GT.MNQ1.AND.NBLOCK.GT.1)) GO TO 10
          NDIJ=MAXA(NJ)-MXMN
          B(NDIJ) = VBR
   10   CONTINUE
      IF(NBLOCK.EQ.1)RETURN
      CALL VWBLOCK(B,A(LLREC),NB,LR,LDB,IBLK)
   20 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VCLEARB(B,MAXA,MNQ,LREC,NBLOCK,LR,IBLK,NUL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CLEAR A FLOATING-POINT ARRAY WITH BLOCKS
CS.    P R O G R A M
CS.        ZA BRISANJE REALNIH VEKTORA SA BLOKOVIMA
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),MAXA(*),MNQ(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' CLEARB'
      CALL VCLEAR(B,NUL)
      IF(NBLOCK.EQ.1)RETURN
      DO 10 NB=1,NBLOCK
  10  CALL VWBLOCK(B,LREC,NB,LR,MAXA(MNQ(NB+1))-MAXA(MNQ(NB)),IBLK)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VRBLOCK(B,LREC,NB,LR,KC,IBLK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO READING OF BLOCK  NB
CS.   P R O G R A M
CS.      ZA CITANJE  BLOKA  NB
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' RBLOCK'
      I=0
      DO 10 N=LREC(NB),LREC(NB+1)-1
      LL=LR
      IF(N.EQ.LREC(NB+1)-1) LL=KC-I
      READ(IBLK,REC=N) (B(I+J),J=1,LL)
   10 I=I+LR
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VWBLOCK(B,LREC,NB,LR,KC,IBLK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO WRITE OF BLOCK  NB
CS.   P R O G R A M
CS.      ZA ZAPISIVANJE BLOKA  NB
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WBLOCK'
      I=0
      DO 10 N=LREC(NB),LREC(NB+1)-1
      LL=LR
      IF(N.EQ.LREC(NB+1)-1) LL=KC-I
      WRITE(IBLK,REC=N) (B(I+J),J=1,LL)
   10 I=I+LR
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VREADDB(B,MAXA,MNQ,LREC,NBLOCK,LR,IBLK,IU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO TRANSFER SEQUENTIAL FILE INTO BLOCKS
CS.   P R O G R A M
CS.      ZA PREBACIVANJE SEKVENCIJALNE DATOTEKE U BLOKOVE
C .
C ......................................................................
C
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),MAXA(*),MNQ(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' READDB'
      DO 10 NB=1,NBLOCK
      LL=MAXA(MNQ(NB+1))-MAXA(MNQ(NB))
      CALL VREADDD(B,LL,IPODS,IU,LDUZI)
   10 CALL VWBLOCK(B,LREC,NB,LR,LL,IBLK)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VWRITEB(B,MAXA,MNQ,LREC,NBLOCK,LR,IBLK,IU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO TRANSFER BLOCKS INTO SEQUENTIAL FILE
CS.   P R O G R A M
CS.      ZA PREBACIVANJE  BLOKOVA U SEKVENCIJALNU DATOTEKU
C .
C ......................................................................
C
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /CDEBUG/ IDEBUG
      DIMENSION B(*),MAXA(*),MNQ(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WRITEB'
      DO 10 NB=1,NBLOCK
      LL=MAXA(MNQ(NB+1))-MAXA(MNQ(NB))
      CALL VRBLOCK(B,LREC,NB,LR,LL,IBLK)
   10 CALL VWRITDD(B,LL,IPODS,IU,LDUZI)
      RETURN
      END
C=======================================================================
      SUBROUTINE VSILPAK(F,FE,LM,ND,CMPC,MPC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE      ASSEMBLE INTERNAL FORCES
CS      RAZMESTANJE UNUTRASNJIH SILA
C
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      DIMENSION F(*),FE(*),LM(*),CMPC(MMP,*),MPC(NEZA1,*)
C
      DO 50 L=1,ND
      LL=LM(L)
      IF(LL.GT.0)THEN
        F(LL)=F(LL)+FE(L)
C
      ELSEIF(LL.LT.0)THEN
        LLP=-LL
        LCM=MPC(1,LLP)
        DO 20 K=1,NEZAV
          LL=MPC(K+1,LLP)
          IF(LL.EQ.0)GO TO 20
          CML=CMPC(LCM,K)
          F(LL)=F(LL)+CML*FE(L)
   20   CONTINUE
      ENDIF
C
   50 CONTINUE
      RETURN
      END
C=======================================================================
      FUNCTION CONDOF(U,CMPC,MPC,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE      TO CALCULATE CONSTRAIND DOF
CS      IZRACUNAVANJE ZAVISNOG STEPENA SLOBODE (MULTIPOINT CONSTRAINT)
C
      COMMON /CDEBUG/ IDEBUG
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      DIMENSION U(*),CMPC(MMP,*),MPC(NEZA1,*)
      IF(IDEBUG.GT.0) PRINT *, ' CONDOF'
C
      IIP=-II
      ICM=MPC(1,IIP)
      CONDOF=0.D0
      DO 20 L=1,NEZAV
      II=MPC(L+1,IIP)
   20 IF(II.GT.0) CONDOF=CONDOF+U(II)*CMPC(ICM,L)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VSWRR(A,MAXA,N,CHAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.       PRINT PROFILED MATRIX
CS        STAMPANJE PROFILISANE MATRICE
C .
C ......................................................................
C
C      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
C     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      CHARACTER*4 CHAR
      COMMON /CDEBUG/ IDEBUG
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      DIMENSION A(*),MAXA(*),V(30)
C
      IF(IDEBUG.GT.0) PRINT *, ' SWRR'
      WRITE(IIZLAZ,5010) CHAR
      DO 20 I=1,N
      DO 10 J=I,N
        V(J)=0.D0
        IJ=J-I
        KK=MAXA(J)+IJ
        IF(KK.LT.MAXA(J+1)) V(J)=A(KK)
10    CONTINUE
        WRITE(IIZLAZ,5000) (V(J),J=1,N)
        V(I)=0.D0
20    CONTINUE
      RETURN
C
 5010 FORMAT(A4)
 5000 FORMAT(30G10.3)
      END
C=======================================================================
C
C=======================================================================
      
C==========================================================================
      SUBROUTINE VBLOCK1
      USE MATRIXINIT
      include 'paka.inc'
c      COMMON A(30000)
c      REAL A
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SRPSKI/ ISRPS
      COMMON /DUPLAP/ IDVA
      COMMON /SCRATC/ ISCRC
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /LSK1/ LSK

      ISCRC=10

      IBLK =46
      NBMAX=500
      IREST=0
CS    LDUZI - DUZINA STAZE KOD JEDINICE SA DIREKTNIM PRISTUPOM
      LDUZI=5120
CZILESK
      IF((LMAX+(NWK)*IDVA).LT.(MTOT-1)) THEN
        NBLOCK=1
        LMNQ=LMAX
        LICPL=LMAX
        LLREC=LMAX
        LR=LDUZI/8
      ELSE
        LAMAX=LMAX+(NWK)*IDVA
        LMNQ=LMAX
        LICPL=LMNQ+NBMAX+1
        LLREC=LICPL+NBMAX
        LMAX=LLREC+NBMAX+1
        CALL VDELJIV(LMAX,2,INDL)
        IF(INDL.EQ.0) LMAX=LMAX+1

CE      LR  JE DUZINA ZAPISA NE U BAJTIMA VEC KAO BROJ CLANOVA
CE      EKVIVALENTNOG VEKTORA U DUPLOJ PRECIZNOSTI
CS      LR  JE DUZINA ZAPISA NE U BAJTIMA VEC KAO BROJ CLANOVA
CS      EKVIVALENTNOG VEKTORA U DUPLOJ PRECIZNOSTI
        LR=LDUZI/8
        MTMM=MTOT-1
        CALL VBLKDEF(MAXA,A(LMNQ),A(LICPL),A(LLREC),
     1              NBLOCK,KC,NBMAX,LR,JEDN,LMAX,MTMM,IDVA,IIZLAZ)

C        WRITE(IIZLAZ,*)'BROJ BLOKOVA=',NBLOCK
C        STOP
CE      OPEN FILES FOR BLOCKS
CS      OTVARANJE FILE-A ZA BLOKOVE
        CALL VBLOKOP(IBLK,LDUZI)
      ENDIF
      write(*,*) 'NBLOCK=',NBLOCK
C      NPODS(JPBR,49)=NBLOCK
C
      LSK=LMAX
      IF(NBLOCK.EQ.1) THEN
         LRAD=LSK+NWK*IDVA
         CALL VCLEAR(A(LSK),NWK)
C         NPODS(JPBR,35)=LMAX13+1
CZILESK         IF(IREST.NE.2)CALL VWRITDD(A(LSK),NWK,IPODS,LMAX13,LDUZI)
      ELSE
         LRAD=MTOT-1
         CALL VCLEAR(A(LSK),KC)
         NUL=KC
       CALL VCLEARB(A(LSK),MAXA,A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,NUL)
         
C         NPODS(JPBR,35)=LMAX13+1
C         IF(IREST.NE.2.AND.(IREST.NE.1.OR.METOD.NE.-1))
C     +     CALL VBLKZAP(A(LSK),MAXA,A(LMNQ))
C         NPODS(JPBR,60)=LMAX13+1
C         IF(IREST.NE.2.AND.(IREST.NE.1.OR.METOD.NE.-1))
C     +     CALL VBLKZAP(A(LSK),MAXA,A(LMNQ))
      ENDIF
      CALL VDELJIV(LRAD,2,INDL)
      IF(INDL.EQ.0) LRAD=LRAD+1
C            LMAX=LRAD
      LMM=LRAD+1
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2000) JEDN,NWK,LMAXM,LMM,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6000) JEDN,NWK,LMAXM,LMM,MTOT
      IF(NBLOCK.GT.1) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2010) LAMAX,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6010) LAMAX,MTOT
      ENDIF
      IF(IREST.NE.2.AND.LMM.GT.MTOT) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2010) LMM,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6010) LMM,MTOT
      STOP ' BLOCKS'
      ENDIF

C-----------------------------------------------------------------------
 2050 FORMAT(///'1',///6X,
     1'PODACI O DIMENZIJAMA MATRICA SISTEMA I RADNOG VEKTORA - A'/
     16X,57('-')///)
 2060 FORMAT(//6X,
     1'Z A   P O D S T R U K T U R U   B R O J ..... JPBR =',I5///)
 2070 FORMAT(//6X,'Z A   K O N T U R N E   C V O R O V E'///)
 2000 FORMAT(
     111X,'BROJ JEDNACINA .............................. JEDN =',I10// 
     211X,'BROJ ELEMENATA U MATRICI SISTEMA ............. NWK =',I10// 
     311X,'MAKSIMALNA VISINA STUBCA U MATRICI SISTEMA ... MHC =',I10// 
     311X,'ZAUZET PROSTOR U RADNOM VEKTORU ............. LMAX =',I10// 
     411X,'RASPOLOZIV PROSTOR U RADNOM VEKTORU ......... MTOT =',I10///)
 2001 FORMAT(
     111X,'BROJ JEDNACINA .............................. JEDN =',I10// 
     211X,'BROJ ELEMENATA U MATRICI SISTEMA ............. NWK =',I10// 
     411X,'RASPOLOZIV PROSTOR U RADNOM VEKTORU ......... MTOT =',I10///)
 2010 FORMAT(///' NEDOVOLJNA DIMENZIJA ZA SMESTANJE MATRICA I VEKTORA ',
     1'SISTEMA :'/' POTREBNA DIMENZIJA,    LMAXM=',I10/' RASPOLOZIVA ',
     2'DIMENZIJA, MTOT =',I10//
     3' TREBA POVECATI RADNI VEKTOR A()'/' PROGRAM STOP')
 2100 FORMAT(//' NEDOVOLJNO MEMORIJE ZA RAD PO BLOKOVIMA'/
     1' POVECATI DUZINU RADNOG VEKTORA NA MTOT =',I10)
 2200 FORMAT(//' PROGRAM NE MOZE RADITI SA PODSTRUKTURAMA PO BLOKOVIMA'/
     1' POVECATI DUZINU RADNOG VEKTORA ILI BROJ PODSTRUKTURA')
 2300 FORMAT(//' PROGRAM NE MOZE RADITI DINAMICKU ANALIZU PO BLOKOVIMA')
 2400 FORMAT(//' PROGRAM NE MOZE RACUNATI KRITICNE SILE PO BLOKOVIMA')
C-----------------------------------------------------------------------
 6050 FORMAT(///'1',///6X,
     1'DATA ABOUT DIMENSIONS OF MATRIX OF SYSTEM AND WORKING VECTOR A'/
     16X,62('-')///)
 6060 FORMAT(//6X,
     1'F O R   S U B S T R U C T U R E   N U M B E R .... JPBR =',I5///)
 6070 FORMAT(//6X,'F O R   C O N T O U R   N O D E S'///)
 6000 FORMAT(
     111X,'NUMBER OF EQUATIONS ......................... JEDN =',I10//
     211X,'NUMBER OF ELEMENTS INTO MATRIX OF SYSTEM ..... NWK =',I20//
     211X,'MAXIMUM HEIGHT COLUMN INTO MATRIX OF SYSTEM .. MHC =',I20//
     311X,'REQUIRED DIMENSION INTO WORKING VECTOR ...... LMAX =',I20//
     411X,'MAXIMUM TOTAL STORAGE AVAILABLE ............. MTOT =',I10///)
 6001 FORMAT(
     111X,'NUMBER OF EQUATIONS ......................... JEDN =',I10//
     211X,'NUMBER OF ELEMENTS INTO MATRIX OF SYSTEM ..... NWK =',I10//
     411X,'MAXIMUM TOTAL STORAGE AVAILABLE ............. MTOT =',I10///)
 6010 FORMAT(///' UNSUFFICIENT DIMENSION OF WORKING VECTOR:'/
     1' REQUIRED DIMENSION LMAXM =',I10/
     2' AVAILABLE DIMENSION MTOT =',I10)
 6100 FORMAT(//' UNSUFFICIENT MEMORY FOR WORK PER BLOCKS'/
     1' REQUIRED DIMENSION OF WORKING VECTOR, MTOT =',I10)
 6200 FORMAT(//' PROGRAM CAN NOT BE WORKED WITH SUBSTRUCTURES AND BLOCKS
     1')
 6300 FORMAT(//' PROGRAM CAN NOT BE WORKED DYNAMIC ANALISYS WITH BLOCKS 
     1')
 6400 FORMAT(//' PROGRAM CAN NOT BE SOLVED CRITICAL LOADS WITH BLOCKS')
C-----------------------------------------------------------------------
      END
C
C======================================================================
      SUBROUTINE VBLOKOP(IBLK,LDUZI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO OPEN FILE FOR BLOCKS
CS.   P R O G R A M
CS.      ZA FILE ZA BLOKOVE - OTVARA SE SAMO PO POTREBI
C .
C ......................................................................
C
C      CHARACTER*20    FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
C     1                FITEMP,FIBLOK,FIBACK
C      COMMON /PAKJED/ FIULAZ,FIZLAZ,FIGRAF,FIELEM,FIPODS,FIFILE,FMODES,
C     1                FITEMP,FIBLOK,FIBACK
C      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
C     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C      COMMON /MASINA/ INDPC,ICRTA
      COMMON /CDEBUG/ IDEBUG

      CHARACTER*20    FIBLOK
      FIBLOK='ZIBLOK              '
C
      IF(IDEBUG.GT.0) PRINT *, ' BLOKOP'
C      IF(INDPC.EQ.2) THEN
C         OPEN (IBLK,ACCESS='DIRECT',
C     1         STATUS='UNKNOWN',FORM='UNFORMATTED',RECL=LDUZI)
C      ELSE
         OPEN (IBLK,FILE=FIBLOK,ACCESS='DIRECT',
     1         STATUS='UNKNOWN',FORM='UNFORMATTED',RECL=LDUZI)
C      ENDIF
      RETURN
      END
C======================================================================
C=======================================================================
      SUBROUTINE VBLKDEF( MAXA,MNQ,ICPL,LREC,
     1NBLOCK,KC,NBMAX,LR,NEQ,LMAX,MTOT,IDVA,IIZLAZ)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO DEFINE OF WORKS WITH BLOCKS
CS.   P R O G R A M
CS.      ZA DEFINISANJE RADA SA BLOKOVIMA
C .
C ......................................................................
C
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
C      COMMON /ULAZNI/ IULAZ,IIZLAZ
      DIMENSION MAXA(*),MNQ(*),ICPL(*),LREC(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' BLKDEF'
      KC=(MTOT-LMAX)/2/IDVA
      IF(KC.LT.NEQ) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2000)NEQ,KC
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6000)NEQ,KC
         STOP
      ENDIF
CE    CONNECTION OF BLOCKS WITH EQUATIONS
CS    VEZA BLOKOVA SA JEDNACINAMA
      NB=0
      N=1
      M1=0
    5 NB=NB+1
      IF(NB.GT.NBMAX) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2010)NBMAX
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6010)NBMAX
         STOP
      ENDIF
      MNQ(NB)=N
      M=0
   10 M=M+MAXA(N+1)-MAXA(N)
      IF(M.GT.KC)GO TO 5
      IF(N.EQ.NEQ)GO TO 30
         IF(M.GT.M1) M1=M
         N=N+1
         IF(M.EQ.KC)GO TO 5
         GO TO 10
   30 NBLOCK=NB
      MNQ(NBLOCK+1)=NEQ+1
CE    VEKTOR OF COUPLED BLOCKS    ICPL
CS    VEKTOR SPREGNUTIH BLOKOVA   ICPL
      DO 50 NB=1,NBLOCK
        K=MNQ(NB)
        DO 40 N=MNQ(NB),MNQ(NB+1)-1
        KH=MAXA(N+1)-MAXA(N)-1
   40   IF(N-KH.LT.K) K=N-KH
      KB=0
   42 KB=KB+1
      IF(K.LT.MNQ(KB+1))GO TO 45
      GO TO 42
   45 ICPL(NB)=KB
   50 CONTINUE
CE    POINTER OF WRITE FOR BLOCK    LREC
CS    POLOZAJ ZAPISA BLOKA    LREC
      LREC(1)=1
      DO 70 NB=1,NBLOCK
      LTH=MAXA(MNQ(NB+1))-MAXA(MNQ(NB))
      L1=LTH/LR
      E=FLOAT(LTH)/FLOAT(LR)
      IF(E.GT.FLOAT(L1)) L1=L1+1
      LREC(NB+1)=LREC(NB)+L1
   70 CONTINUE
C
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2030) NBLOCK,M1
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6030) NBLOCK,M1
      RETURN
C-----------------------------------------------------------------------
 2030 FORMAT('1',///6X,
     1'PODACI O BLOKOVIMA MATRICE KRUTOSTI'/
     16X,35('-')///
     111X,'BROJ BLOKOVA .............................. NBLOCK =',I9///
     211X,'MAKSIMALNA DUZINA BLOKA ....................... M1 =',I9)
 2000 FORMAT(/1X,'NEDOVOLJNA DIMENZIJA BLOKA :'
     1/1X,'MINIMALNA POTREBNA DIMENZIJA, JEDN =',I10
     2/1X,'RASPOLOZIVA DIMENZIJA,          KC =',I10)
 2010 FORMAT(/1X,'NEDOVOLJAN BROJ BLOKOVA :'
     1/1X,'RASPOLOZIVI BROJ BLOKOVA,    NBMAX =',I10)
C-----------------------------------------------------------------------
 6030 FORMAT('1',///6X,
     1'DATA ABOUT BLOCKS OF STIFFNESS MATRIX'/
     16X,37('-')///
     111X,'NUMBER OF BLOCKS .......................... NBLOCK =',I9///
     211X,'MAXIMUM LENGTH OF BLOCKS ...................... M1 =',I9)
 6000 FORMAT(/1X,'UNSUFFICIENT DIMENSION OF BLOCK:'
     1/1X,'MINIMUM REQUIRED DIMENSION  JEDN =',I10
     2/1X,'AVAILABLE DIMENSION ........  KC =',I10)
 6010 FORMAT(/1X,'UNSUFFICINET NUMBER OF BLOCKS :'
     1/1X,'AVAILABLE NUMBER OF BLOCKS  NBMAX =',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C=======================================================================
      SUBROUTINE VCLEARD(A,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
      IF(IDEBUG.GT.0) PRINT *, ' CLEAR'
      DO 10 I=1,N
   10 A(I)=0.0D0
      RETURN
      END
C=======================================================================

C=======================================================================
      SUBROUTINE VCLEAR(A,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CLEAR A FLOATING-POINT ARRAY
CS.    P R O G R A M
CS.        ZA BRISANJE REALNIH VEKTORA
C .
CE.    I=1,N  (N - LENGTH OF VECTOR -  A)
CE.         A(I) - CLEAR VECTOR
CS.    I=1,N  (N - DUZINA VEKTORA -  A)
CS.         A(I) - VEKTOR KOJI SE BRISE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' CLEAR'
      DO 10 I=1,N
   10 A(I)=0.0D0
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE VICLEAR(IA,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CLEAR INTEGER ARRAY
CS.    P R O G R A M
CS.        ZA BRISANJE CELOBROJNIH VEKTORA
C .
CE.    I=1,N  (N - LENGTH OF VECTOR - IA)
CE.        IA(I) - CLEAR VECTOR
CS.    I=1,N  (N - DUZINA VEKTORA - IA)
CS.        IA(I) - VEKTOR KOJI SE BRISE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ICLEAR'
      DO 10 I=1,N
   10 IA(I)=0
      RETURN
      END
C=======================================================================
      SUBROUTINE VICLEAR8(IA,N)

      integer*8 IA,N
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CLEAR INTEGER ARRAY
CS.    P R O G R A M
CS.        ZA BRISANJE CELOBROJNIH VEKTORA
C .
CE.    I=1,N  (N - LENGTH OF VECTOR - IA)
CE.        IA(I) - CLEAR VECTOR
CS.    I=1,N  (N - DUZINA VEKTORA - IA)
CS.        IA(I) - VEKTOR KOJI SE BRISE
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' ICLEAR8'
      DO 10 I=1,N
   10 IA(I)=0
      RETURN
      END
C=======================================================================
      SUBROUTINE VERROR(I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   PROGRAM
CE.      TO PRINT ERRORS
CS.   PROGRAM
CS.      ZA STAMPANJE GRESAKA I UPOZORENJA
C .
C ......................................................................
C
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
C      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
C     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' ERROR '
C
      GO TO (1,2,3),I
C
    1 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2000) LMAX,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6000) LMAX,MTOT
      STOP
C
    2 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2000) LMAX,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6000) LMAX,MTOT
      STOP
C
    3 CONTINUE
      IF(ISRPS.EQ.0)
     1WRITE(IIZLAZ,2000) LMAX,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IIZLAZ,6000) LMAX,MTOT
      STOP
C-----------------------------------------------------------------------
 2000 FORMAT(///' NEDOVOLJNA DIMENZIJA U OSNOVNOM RADNOM VEKTORU A ZA UC
     1ITAVANJE ULAZNIH PODATAKA'/' POTREBNA DIMENZIJA, LMAX=',I10/
     2' RASPOLOZIVA DIMENZIJA, MTOT=',I10)
C-----------------------------------------------------------------------
 6000 FORMAT(///'  ISUFFICIENT DIMENSION OF WORKING VECTOR FOR READING I
     1NPUT DATA'/
     1' REQUIRED STORAGE  .......  LMAX =',I10/
     2' MAXIMUM STORAGE AVAILABLE  MTOT =',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C=======================================================================
      SUBROUTINE VIREADD(IA,N,II,LS,LD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO READ INTEGER VECTOR ON DIRECT ACCESS FILE
CS.    P R O G R A M
CS        ZA CITANJE CELOBROJNOG VEKTORA SA FILE SA DIREKTNIM PRISTUPOM
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IA(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' IREADD'
      IF(N.EQ.0) RETURN
      NK=N*4/LD
      NN=NK*LD/4
      IF(N.GT.NN) NK=NK+1
      DO 10 K=1,NK
         IK=(K-1)*LD/4+1
         NN=K*LD/4
         IF(K.EQ.NK) NN=N
         LS=LS+1
         READ(II,REC=LS) (IA(I),I=IK,NN)
   10 CONTINUE
C      IF(II.EQ.8) THEN
C         WRITE(3,*) 'IRN,LS,NK,LD,II',N,LS,NK,LD,II
C         CALL VIWRR(IA,N,'IRN ')
C      ENDIF
      RETURN
      END   
C======================================================================
      SUBROUTINE VWRITDD(A,N,II,LS,LD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO WRITE REAL VECTOR AT DIRECT ACCESS DISK
CS.    P R O G R A M
CS        ZA ZAPISIVANJE REALNOG VEKTORA NA DISK SA DIREKTNIM PRISTUPOM
C .
C ......................................................................
C
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' WRITDD'
C      IF(ISKDSK.EQ.0.AND.N.EQ.NWK)RETURN
      IF(N.EQ.0) RETURN
      NK=N*8/LD
      NN=NK*LD/8
      IF(N.GT.NN) NK=NK+1
C      PRINT *, 'II,N,LD,NK',II,N,LD,NK
      DO 10 K=1,NK
         IK=(K-1)*LD/8+1
         NN=K*LD/8
         IF(K.EQ.NK) NN=N
         LS=LS+1
C         PRINT *,'K,LS,IK,NN',K,LS,IK,NN
         WRITE(II,REC=LS) (A(I),I=IK,NN)
   10 CONTINUE
C      IF(II.EQ.8) THEN
C         WRITE(3,*) 'WN,LS,NK,LD,II',N,LS,NK,LD,II
C         CALL VWRR(A,N,'WN  ')
C      ENDIF
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE VREADDD(A,N,II,LS,LD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO READ REAL VECTOR ON DIRECT ACCESS FILE
CS.    P R O G R A M
CS        ZA CITANJE REALNOG VEKTORA SA FILE SA DIREKTNIM PRISTUPOM
C .
C ......................................................................
C
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' READDD'
C      IF(ISKDSK.EQ.0.AND.N.EQ.NWK)RETURN
      IF(N.EQ.0) RETURN
      NK=N*8/LD
      NN=NK*LD/8
      IF(N.GT.NN) NK=NK+1
      DO 10 K=1,NK
         IK=(K-1)*LD/8+1
         NN=K*LD/8
         IF(K.EQ.NK) NN=N
         LS=LS+1
         READ(II,REC=LS) (A(I),I=IK,NN)
   10 CONTINUE
C      IF(II.EQ.8) THEN
C         WRITE(3,*) 'RN,LS,NK,LD,II',N,LS,NK,LD,II
C         CALL VWRR(A,N,'RN  ')
C      ENDIF
      RETURN
      END
C=======================================================================      