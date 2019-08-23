C======================================================================
CE        U ovom fajlu su potprogrami koji vrse novo tackanje u RB drvo
C======================================================================
C=======================================================================
      SUBROUTINE SPAKUJMT(SK,MAXA,SKE,LM,ND)
      USE DRAKCE8
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'paka.inc'
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG

      DIMENSION SKE(*),LM(*),MAXA(*),SK(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUJMT'
      IF(NBLOCK.EQ.1) THEN
!           CALL ISPAKGMT(SK,AIROWS,MAXA,SKE,LM,ND,0,
!     &               A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,A(LCMPC),A(LMPC))
      IF(ICCGG.EQ.2) THEN
      CALL sparseassembler_addnonsymmatrix
     1 (ND,LM,SKE,NEZAV,NMPC,A(LCMPC),A(LMPC))        
      ELSE !IF(ICCGG.EQ.2) THEN         
      CALL sparseassembler_addelemmatrix
     1 (ND,LM,SKE,NEZAV,NMPC,A(LCMPC),A(LMPC))
      ENDIF !IF(ICCGG.EQ.2) THEN
      
      ELSE
          STOP 'NE RADI S BLOKOVIMA'
      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ISPAKGMT(SK,IROW,MAXA,SKE,LM,ND,INDD,
     &                  MNQ,LREC,NBLOCK,LR,IBLK,CMPC,MPC)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'paka.inc'
      
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
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
      DIMENSION LM(*),MAXA(*),MNQ(*),LREC(*),
     &          CMPC(MMP,*),MPC(NEZA1,*),IROW(*),SK(*),SKE(*)
      IF(IDEBUG.GT.0) PRINT *, ' SPAKUA'
C

      
!      
!        MNQ0=1
!        MNQ1=JEDN
!        MXMN=0
!        IF(INDD.EQ.1)THEN
!        REWIND ISCRC
!   10   CONTINUE
!           READ(ISCRC,END=15,ERR=999)
!     &     ND,(LM(I),I=1,ND),(SKE(I),I=1,ND*(ND+1)/2)       
!        ENDIF
!C-----------------------------------------------
!   11 NDI=0
!
!      DO 200 I=1,ND
!      II=LM(I)
!
!      IF(II.LT.0)THEN !pocetak vezanih pomeranja
!        IIP=-II
!        ICM=MPC(1,IIP)
!        DO 320 L=1,NEZAV
!          II=MPC(L+1,IIP)
!          IF(II.LT.MNQ0) GO TO 320
!          CMI=CMPC(ICM,L)
!          !MI=MAXA(II)-MXMN GDE LI SE OVO KORISTI
!          KS=I
!          DO 310 J=1,ND
!            JJ=LM(J)
!            IF(JJ)303,310,307
!  303       JJP=-JJ
!            JCM=MPC(1,JJP)
!              KSS=KS
!              IF(J.GE.I)KSS=J+NDI
!              
!              DO K=1,NEZAV
!                JJ=MPC(K+1,JJP)
!                IF(JJ.GE.0) THEN
!                IJ=II-JJ
!                IF(IJ.GE.0) THEN
!                  CMJ=CMPC(JCM,K)
!                ! NADJOH MI
!                  KK=MI+IJ+1
!                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
!                  
! 1666               CONTINUE
!                    KK=KK-1
!                    IF(IROW(KK).NE.II-IJ)GOTO 1666
!                    
!                  SK(KK)=SK(KK)+CMI*CMJ*SKE(KSS)
!                ENDIF
!                ENDIF
!               ENDDO 
!  !318         CONTINUE
!              
!              
!              GO TO 310
!C
!C
!  307         IJ=II-JJ
!              IF(IJ)310,311,311
!  311         CONTINUE
!                ! KK=MI+IJ+1 EVO GA MI
!                !  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
! 6616               CONTINUE
!                    KK=KK-1
!                    IF(IROW(KK).NE.II-IJ)GOTO 6616
!                KSS=KS
!                IF(J.GE.I)KSS=J+NDI
!                SK(KK)=SK(KK)+CMI*SKE(KSS)
!  310         KS=KS+ND-J
!  320   CONTINUE
!        GO TO 200
!      ENDIF      !IF(II.LT.0)THEN kraj vezanih pomeranja
!      
!      
!      
!      
!      
!      IF(II.LT.MNQ0) GO TO 200
!      MI=MAXA(II)-MXMN
!      KS=I
!      IVRS=0
!      
!      
!      DO 220 J=1,ND
!      JJ=LM(J)
!      IF(JJ.GT.0)IVRS=IVRS+1
!      IF(JJ)420,220,110
!C
!C
!  420       JJP=-JJ
!            JCM=MPC(1,JJP)
!              KSS=KS
!              IF(J.GE.I)KSS=J+NDI
!              DO K=1,NEZAV
!                JJ=MPC(K+1,JJP)
!                IF(JJ.GE.0)THEN
!                CMJ=CMPC(JCM,K)
!                IJ=II-JJ
!                IF(IJ.GE.0) THEN
!                  KK=MI+IJ+1
!                  IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
! 2666               CONTINUE
!                    KK=KK-1
!                    IF(IROW(KK).NE.II-IJ)GOTO 2666
!                  SK(KK)=SK(KK)+CMJ*SKE(KSS)
!                ENDIF
!                ENDIF
!              ENDDO
!              
!                GO TO 220
!C
!  110 IJ=II-JJ
!      IF(IJ.GE.0) THEN
!
!         KK=MI+IJ+1
!      IF(KK.GT.MAXA(II+1))KK=MAXA(II+1)
! 6626     CONTINUE
!          KK=KK-1
!          IF(IROW(KK).NE.II-IJ)GOTO 6626
!         KSS=KS
!         IF(J.GE.I)KSS=J+NDI
!         SK(KK)=SK(KK)+SKE(KSS)
!      ENDIF
!      
!  220 KS=KS+ND-J
!      
!      
!      
!  200 NDI=NDI+ND-I
!C
!      IF(INDD.EQ.1) GO TO 10
!   15 IF(NBLOCK.GT.1) CALL WBLOCK(SK,LREC,KB0,LR,LDB,IBLK)  
      RETURN
999   PRINT *,'ERROR: reading element stifness matrix from disk'
      STOP
      END
C===================================================================

C=======================================================================
      SUBROUTINE BUSYMATRICA()
      USE MATRICA
      USE STIFFNESS
      USE DRAKCE8
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   OVO JE SUBRUTINA KOJA SE KORISTI ZA FORMIRANJE MATRICE KRUTOSTI    
CS    UMESTO DRAKCETOVOG TACKANJA. PROBLEM JE STO SE POZIVA NA VISE MESTA    
C .   NAKON PETLJI PO GRUPAMA ELEMENATA, TAKO DA JE TESKA ZA TESTIRANJE
C ......................................................................
      include 'paka.inc'
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /SOPSVR/ ISOPS,ISTYP,NSOPV,ISTSV,IPROV,IPROL
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' BUSYMATRICA'  
      if(.not.allocated(rows)) then
            call sparseassembler_getnz(nonzeros)
            if (nonzeros.ge.2147483647) stop 'izmeni imaxa u 64bit'
! 2147483647 =(2**32)-1)/2 Maximum value for a variable of type int.
! pored imaxa svuda gde se koristi neki clan niza mora da bude int64  
            
            allocate(rows(nonzeros),STAT=istat)
            if(istat.ne.0) stop 'error allocating rows'
            allocate(iirows(nonzeros),STAT=istat)
            allocate(columns(nonzeros),STAT=istat)
            if(istat.ne.0) stop 'error allocating columns'
            allocate(iicolumns(nonzeros),STAT=istat)
            allocate(IMAXA(JEDN+1),STAT=istat)
            if(istat.ne.0) stop 'error allocating IMAXA'
            allocate(ALSK(nonzeros),STAT=istat)
            if(istat.ne.0) stop 'error allocating ALSK'
             ALLOCATE (AIROWS(nonzeros*2), STAT = istat)
         if(istat.ne.0) stop 'error allocating AIROWS'
              IF(NDIN.GT.0.OR.ISOPS.GT.0) THEN 
                IF(IMASS.GE.1) THEN
              ALLOCATE (ALSM(nonzeros), STAT = iAllocateStatus)
          IF (iAllocateStatus /= 0) write(3,*)'ALSM Not enough memory *'
          IF (iAllocateStatus /= 0) STOP '*** ALSM Not enough memory *'
                ENDIF
             ENDIF
             IF(NDIN.GT.0) THEN
                IF(IDAMP.GT.0) THEN
                   ALLOCATE (ALSC(nonzeros), STAT = iAllocateStatus)
          IF (iAllocateStatus /= 0) write(3,*)'ALSC Not enough memory *'
          IF (iAllocateStatus /= 0) STOP '*** ALSC Not enough memory *'
                ENDIF
             ENDIF
          endif

      CALL sparseassembler_getsparse(nonzeros,rows,columns,ALSK,IMAXA)
      IMAXA(JEDN+1) = nonzeros+1
      do i=1,nonzeros
            iirows(i)= rows(i)
            iicolumns(i)=columns(i)
            AIROWS(i)= rows(i)
            AIROWS(nonzeros+i)=columns(i)
      enddo
   
          CALL sparseassembler_kill()
          
      RETURN
      END