C======================================================================
CE        ALL SUBROUTINES IN THIS FILE ARE PROGRAMD BY DRAKCE
C======================================================================
      SUBROUTINE ISPAKUJ2(IROW,NAXA,MAXA,NWK,JEDN)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      include 'paka.inc'
      
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /SKDISK/ ISKDSK
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /CDEBUG/ IDEBUG
      COMMON /PREDOT/ IPREDOT, NNWK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
C
      INTEGER*8 MAXA,NWK,I8,NZERO8,J8
      DIMENSION MAXA(JEDN+1),LM(100),IROW(*),NAXA(JEDN+1),IND(JEDN)
     +,LMS(100)
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUJ'
      NBLOCK=1
      IF(NBLOCK.EQ.1) THEN
CE       WITHOUT BLOCKS
CS       BEZ BLOKOVA
CS     UNULJIVANJE MATRICE
CE     CLEARING MATRIX ISK 
      DO 30 I=1, JEDN+1
         NAXA(I) = 1
C         NAXA(I) = I
C         IROW(I) = I
  30  CONTINUE
      NNWK = 2*JEDN
C      NAXA(JEDN+1) = JEDN+1
      REWIND(IDRAKCE)
 1000    FORMAT(101I5)
 
C     PREDTACKANJE - POMOCU KOJEG UTVRDJUJEMO TEORETSKI MAKSIMUM BROJA
C                    NE NULA CLANOVA (NEKI CLANOVI SE EVIDENTIRAJU VISE PUTA)
C                    FORMIRAJU SE VISINE SA NE NULA CLANOVIMA (CUVAJU SE U NAXA)
      IPREDOT = 1      
      DO 10 J=1,NELUK
         READ(IDRAKCE) ND,(LM(I),I=1,ND)
         DO I = 1,ND
           DO K = I,ND
             IF(LM(I).LT.LM(K))THEN
               ILM = LM(I)
               LM(I) = LM(K)
               LM(K) = ILM
             END IF
           END DO
         END DO
         LMNQ=1
         LLREC=1
         LR=1
C
C         CALL IIISWR(LM,ND,' LM ')  
C
         CALL ISPAKUA2(IROW,NAXA,MAXA,LM,ND,0,
     &               A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,A(LCMPC),A(LMPC),
     &               NWK,JEDN)
  10  CONTINUE
        CALL VREM(6)
C     ZAVRSENO PRED TACKANJE I FORMIRA SE NAXA KAO DA SU SVI ELEMENTI
C     NENULA CLANOVI POMOCU VISINA
      IPOM = 1
      NNZERO = 0
      DO I=1,JEDN+1
         IPOM = NNZERO + 1
         NNZERO = NAXA(I)+NNZERO
         NAXA(I) = IPOM
      ENDDO
      IF(IPOM+IROWS.GT.MTOT) CALL ERROR(1)
C     UNULJAVA SE IROW NIZ JER NEKI CLANOVI MOGU BITI RAZLICITI OD NULE
      DO I8=1,IPOM
        IROW(I8)=0
      ENDDO
      REWIND(IDRAKCE)     
        CALL VREM(6)
C     POCINJE TACKANJE - FORMIRANJE IROW NIZA
      IPREDOT = 2
      DO 20 J=1,NELUK
         READ(IDRAKCE) ND,(LM(I),I=1,ND)
         LMNQ=1
         LLREC=1
         LR=1
C
C         CALL IIISWR(LM,ND,' LM ')  
C
         CALL ISPAKUA2(IROW,NAXA,MAXA,LM,ND,0,
     &               A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,A(LCMPC),A(LMPC),
     &               NWK,JEDN)
  20  CONTINUE
        CALL VREM(6)
C     FORMIRAN NIZ IROW - IZBACUJU SE NULA CLANOVI (ZBOG DUPLIRANJA)
      NWK = 1
      J8 = 1
      I8 = NAXA(JEDN+1)
      NN = 1
    9 NN = NN+1
      IF(NN.GT.JEDN+1)GOTO 99
      IF(IROW(NAXA(NN)-1).GT.0)GOTO 9
      DO I8=NAXA(NN-1),NAXA(NN)-1
        IF(IROW(I8).EQ.0)GOTO 99
      ENDDO
   99 J8=I8
      NWK = J8
      DO I=NN,JEDN
         DO I8=NAXA(I),NAXA(I+1)-1
           IF(IROW(I8).GT.0)THEN
             IROW(J8)=IROW(I8)
             J8=J8+1
             IROW(I8) = 0
           ENDIF
         ENDDO
         NAXA(I) = NWK
         NWK = J8
      ENDDO
      NWKMAX = NAXA(JEDN+1)-1
      NAXA(JEDN+1)=NWK        
C     SORTIRANJE IROWS NIZA (NIJE POTREBNO)
C                KORISTI VISE OD 25% VREMENA
      ISORT = 0
      IF(ISORT.EQ.1)THEN
        CALL VREM(6)
        J8 = 1
        IF(ICCGG.EQ.1)THEN
          DO 100 I = 1, JEDN
            DO 101 I8 = NAXA(I), NAXA(I+1) - 1
              IND(IROW(I8)) = 1         
  101       CONTINUE
            DO 102 I8 = I, JEDN
              IF(IND(I8).EQ.1)THEN
                IROW(J8) = I8
                IND(I8) = 0
                J8 = J8 + 1
              ENDIF
  102       CONTINUE
  100     CONTINUE
        ELSE
          DO 200 I = 1, JEDN
            DO 201 I8 = NAXA(I), NAXA(I+1) - 1
              IND(IROW(I8)) = 1         
  201       CONTINUE
            DO 202 I8 = I, 1, -1
              IF(IND(I8).EQ.1)THEN
                IROW(J8) = I8
                IND(I8) = 0
                J8 = J8 + 1
              ENDIF
  202       CONTINUE
  200     CONTINUE
        ENDIF
      ENDIF
      NZERO = MAXA(JEDN+1)-NAXA(JEDN+1)
      NNZERO = NAXA(JEDN+1)-1
      NWK = MAXA(JEDN+1)-1

        CALL VREM(6)
      CALL ICOUNT2(NWK,JEDN,NZERO8,NEED1,NEED2,NEED3,NNZERO,NAXA,
     1             NWKMAX,NELUK)
C PROVERI TREBA DA SE IZBACI      NNZERO=NWK-NZERO8
      ELSE
CE       WITH BLOCKS
CS       SA BLOKOVIMA
         STOP 'NERADI SA BLOKOVIMA'
C         WRITE(ISCRC)ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*(ND+1)/2)
      ENDIF
      RETURN
      END
C=========================================================================================
      SUBROUTINE ICOUNT2(NWK,NN,NZERO,NEED1,NEED2,NEED3,NNZERO,NAXA,
     1                   NWKMAX,NELUK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      INTEGER*8 NWK,I8,NZERO
      DIMENSION NAXA(NN+1)
      NZERO=0
      NNZERO = NAXA(NN + 1) - 1
      NZERO = NWK - NNZERO
      A1=FLOAT(NZERO)/FLOAT(NWK)*100.
      A2=FLOAT(NNZERO)/FLOAT(NWK)*100.
      NEED1=NNZERO*2+NN*5
      NEED2=NNZERO
      NEED3=NEED1+NEED2
      IF(ISRPS.EQ.0)
     1WRITE(*,20)NN,NWK,NZERO,A1,NNZERO,A2,NEED1,NEED2,NEED3,
     1           NWKMAX,NELUK
      IF(ISRPS.EQ.1)
     1WRITE(*,30)NN,NWK,NZERO,A1,NNZERO,A2,NEED1,NEED2,NEED3,
     1           NWKMAX,NELUK
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,20)NN,NWK,NZERO,A1,NNZERO,A2,NEED1,NEED2,NEED3,
     1           NWKMAX,NELUK
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,30)NN,NWK,NZERO,A1,NNZERO,A2,NEED1,NEED2,NEED3,
     1           NWKMAX,NELUK
20    FORMAT(//' DIMENZIJA MATRICE    : ',I12/
     *         ' UKUPAN BROJ CLANOVA  : ',I12/
     *         ' BROJ CLANOVA = 0     : ',I12,'  ILI ',F6.3,' %'/
     *         ' BROJ NE NULA CLANOVA : ',I12,'  ILI ',F6.3,' %'/
     *         ' POTREBNA VELICINA GLAVNOG VEKTORA  : ',I12/
     *         ' POTREBNA VELICINA POMOCNOG VEKTORA : ',I12/
     *' UKUPNA POTREBNA VELICINA GLAVNOG I POMOCNOG VEKTORA : ',I9/
     *' VELICINA VEKTORA ZA ODREDJIVANJE NENULA CLANOVA : ',I9/
     *         ' BROJ ELEMENATA       : ',I12//)
30    FORMAT(//' DIMENSION OF SYSTEM MATRIX: ',I12/
     *         ' TOTAL NUMBER OF MEMBERS   : ',I12/
     *         ' NUMBER OF MEMBERS = 0     : ',I12,'  OR ',F6.3,' %'/
     *         ' NUMBER OF NONZERO MEMBERS : ',I12,'  OR ',F6.3,' %'/
     *         ' REQUIRED DIMENSION IN MAIN ARRAY     : ',I12/
     *         ' REQUIRED DIMENSION IN AUXILIARY ARRAY: ',I12/
     *' TOTAL REQUIRED DIMENSION IN MAIN AND AUXILIARY ARRAY: ',I12/
     *' SIZE OF TEMPORARY ARRAY NEEDED FOR NONZERO ELEMENTS : ',I9/
     *         ' NUMBER OF ELEMENTS        : ',I12//)
      RETURN
      END
C=======================================================================
      SUBROUTINE ISPAKUA2(IROW,NAXA,MAXA,LM,ND,INDD,
     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC,NWK,JEDN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.      TO ADD ELEMENT STIFFNESS MATRICES TO GLOBAL STIFFNESS MATRIX
CS.   P R O G R A M
CS.      ZA RAZMESTANJE MATRICA ELEMENATA U SISTEM - BEZ BLOKOVA
C .
C ......................................................................
C
      include 'paka.inc'
      
!      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /CDEBUG/ IDEBUG
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /PROBAS/ IILS
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      INTEGER*8 MAXA,NWK,MI,KK,KKNWK
      DIMENSION LM(*),MAXA(JEDN+1),MNQ(*),LREC(*),
     &          CMPC(MMP,*),MPC(NEZA1,*),IROW(*), NAXA(JEDN+1)
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUA2'
C
C
CS  PETLJA PO BLOKOVIMA
      DO 20 KB0=1,NBLOCK
      IF(NBLOCK.EQ.1)THEN
C      IF(IILS.NE.-1) CALL KORIGC(A(LCORUL),A(LCMPC),MMP,NP,NEZAV)
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
      CALL RBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
CS  PETLJA PO ELEMENTIMA
    9   REWIND ISCRC
   10   IF(ICCGG.EQ.2) THEN
           READ(ISCRC,END=15,ERR=999)
     &     ND,(LM(I),I=1,ND)
        ELSE
           READ(ISCRC,END=15,ERR=999)
     &     ND,(LM(I),I=1,ND)
        ENDIF
C-----------------------------------------------
   11 NDI=0
99999 FORMAT(5I5)
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
            IF(ICCGG.EQ.2) THEN
              KSSU=(I-1)*ND+J
              KSSL=(J-1)*ND+I
            ELSE
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
            ENDIF
              DO 318 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 318
                IJ=II-JJ
                IF(IJ)318,314,314
  314           CMJ=CMPC(JCM,K)
                IF(ICCGG.EQ.2) THEN
                  KK=MI-IJ
                  KKNWK=KK+NWK
                ELSE          
                  KK=MI+IJ
                ENDIF
                CALL ADDRC(II, JJ, MAXA, NAXA, JEDN, IROW, ICCGG)
  318         CONTINUE
              GO TO 310
C
C
  307         IJ=II-JJ
              IF(IJ)310,311,311
  311         IF(ICCGG.EQ.2) THEN
                KK=MI-IJ
                KKNWK=KK+NWK
                KSSU=(I-1)*ND+J
                KSSL=(J-1)*ND+I
              ELSE
                KK=MI+IJ
                KSS=KS
                IF(J.GE.I)KSS=J+NDI
              ENDIF
              CALL ADDRC(II, JJ, MAXA, NAXA, JEDN, IROW, ICCGG)
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
            IF(ICCGG.EQ.2) THEN
              KSSU=(I-1)*ND+J
              KSSL=(J-1)*ND+I
            ELSE
              KSS=KS
              IF(J.GE.I)KSS=J+NDI
            ENDIF
              DO 418 K=1,NEZAV
                JJ=MPC(K+1,JJP)
                IF(JJ.EQ.0)GO TO 418
                CMJ=CMPC(JCM,K)
C
                IJ=II-JJ
                IF(IJ)418,415,415
  415           IF(ICCGG.EQ.2) THEN
                  KK=MI-IJ
                  KKNWK=KK+NWK
                ELSE
                  KK=MI+IJ
                ENDIF
                CALL ADDRC(II, JJ, MAXA, NAXA, JEDN, IROW, ICCGG)
  418           CONTINUE
                GO TO 220
C
C
  110 IJ=II-JJ
      IF(IJ)220,210,210
  210 IF(ICCGG.EQ.2) THEN
         KK=MI-IJ
         KKNWK=KK+NWK
         KSSU=(I-1)*ND+J
         KSSL=(J-1)*ND+I
      ELSE
         KK=MI+IJ
         KSS=KS
         IF(J.GE.I)KSS=J+NDI
      ENDIF
      CALL ADDRC(II, JJ, MAXA, NAXA, JEDN, IROW, ICCGG)
  220 KS=KS+ND-J
  200 NDI=NDI+ND-I
C
      IF(INDD.EQ.1.OR.NBLOCK.GT.1) GO TO 10
   15 IF(NBLOCK.GT.1) CALL WBLOCK(SK,LREC,KB0,LR,LDB,IBLK)
C
   20 CONTINUE
      RETURN
999   PRINT *,'ERROR: reading element stifness matrix from disk'
      STOP
      END
C=========================================================================================
       
      SUBROUTINE ADDRC(IC, IR, MAXA, NAXA, JEDN, IROW, ICCGG)
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /PREDOT/ IPREDOT, NNWK
      INTEGER*8 MAXA,I8
	DIMENSION MAXA(JEDN+1), NAXA(JEDN+1), IROW(NNWK+2) !(NAXA(JEDN+1))
!	IC = 1
!   10	IC = IC + 1
!      IF (MAXA(IC).LE.KK) GO TO 10
!      IC = IC - 1
!      IR = IC + MAXA(IC) - KK
      IF(IPREDOT.EQ.1)THEN
         IF ((ICCGG.EQ.1).AND.(IR.LT.IC))
     +   NAXA(IR) = NAXA(IR) + 1
         IF ((ICCGG.EQ.-1).AND.(IR.LT.IC))
     +   NAXA(IC) = NAXA(IC) + 1
      ELSE
         IF (ICCGG.EQ.1) THEN
            DO 11 I8 = NAXA(IR), NAXA(IR+1)-1
               IF(IROW(I8).EQ.0) IROW(I8) = IC
               IF(IROW(I8).EQ.IC)RETURN
   11       CONTINUE
         ELSE
            DO 12 I8 = NAXA(IC), NAXA(IC+1)-1
               IF(IROW(I8).EQ.0) IROW(I8) = IR
               IF(IROW(I8).EQ.IR)RETURN
   12       CONTINUE
         ENDIF
      ENDIF
         RETURN
      IF (ICCGG.EQ.1) THEN
        DO 20 I8 = NAXA(IR), NAXA(IR+1)-1
            IF(IROW(I8).EQ.IC)RETURN
20      CONTINUE     
        LMAX = LMAX+1
        IF(LMAX.GT.MTOT) CALL ERROR(1)
        DO 40 I = JEDN+1, IR+1, -1
            IROW(NAXA(I)) = IROW(NAXA(I-1))
            NAXA(I) = NAXA(I) + 1
40      CONTINUE
        I8 = NAXA(IR)
        IROW(I8) = IC
      ELSE
        DO 21 I8 = NAXA(IC), NAXA(IC+1)-1
            IF(IROW(I8).EQ.IR)RETURN
21      CONTINUE     
        LMAX = LMAX+1
        IF(LMAX.GT.MTOT) CALL ERROR(1)
        DO 41 I = JEDN+1, IC+1, -1
            IROW(NAXA(I)) = IROW(NAXA(I-1))
            NAXA(I) = NAXA(I) + 1
41      CONTINUE
        I8 = NAXA(IC)
        IROW(I8) = IR
      ENDIF
      RETURN
      END
    
C==========================================================================  
      SUBROUTINE STIROW(IROW, MAXA, JEDN, IDDOT, ICCGG, IMUMPS)
      COMMON /IMEULZ/ PAKLST,PAKUNV,PAKNEU
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      DIMENSION MAXA(JEDN+1), IROW((MAXA(JEDN+1)-1)*(IMUMPS+1))
      CHARACTER*40 FNAME,FNAME1
      CHARACTER *24 PAKLST,PAKUNV,PAKNEU
      IB=INDEX(PAKLST,'.')
      IUNIT=9999
      WRITE(FNAME1,11)JEDN,ICCGG,IMUMPS,IDDOT
11    FORMAT('_',I10,'_',I2,'_',I1,'_',I1,'.ROW')
      FNAME=PAKLST(1:IB-1)//FNAME1
      OPEN(UNIT=IUNIT,NAME=FNAME)
      WRITE(IUNIT,*)'JEDN =', JEDN, ' ICCGG =', ICCGG, 
     +              ' IMUMPS =', IMUMPS, ' NWK =', MAXA(JEDN+1)-1
      WRITE(IUNIT,*)'MAXA = '
      WRITE(IUNIT,51)(MAXA(I),I=1,JEDN+1)
      DO 10 I=1, JEDN
        WRITE(IUNIT,61)I
        WRITE(IUNIT,51)(IROW(J),J=MAXA(I),MAXA(I+1)-1)
10    CONTINUE
      IF(IMUMPS.EQ.1)THEN
        NWK = MAXA(JEDN+1)-1
        DO 20 I=1, JEDN
          WRITE(IUNIT,61)I
          WRITE(IUNIT,51)(IROW(J),J=MAXA(I)+NWK,MAXA(I+1)-1+NWK)
20      CONTINUE
      END IF
      CLOSE(IUNIT)
51    FORMAT(5000I11)
61    FORMAT('COL',I10)
62    FORMAT('ROW',I10)
      RETURN
      END
C==========================================================================  

      SUBROUTINE FORMCOL(ICOL,MAXA,JEDN)
      DIMENSION MAXA(JEDN+1), ICOL(MAXA(JEDN+1)-1)
      DO 10 I=1,JEDN
        DO 20 J=MAXA(I),MAXA(I+1)-1
            ICOL(J) = I
20      CONTINUE     
10    CONTINUE     
      RETURN
      END
C==========================================================================  
      
