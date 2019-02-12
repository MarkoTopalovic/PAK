C======================================================================
C
C======================================================================
      SUBROUTINE VTGRMAT(AK,NUMMAT,NETIP,II,NDIM,NASLOV,GAMA)
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO PRINT MATERIALS AND PROPERTIES IN NEU FILE
C .
C ......................................................................
C
      CHARACTER*80 NASLOV
      DIMENSION AK(3,5,*)
C
      IND=-1
      ITYP=100
	VERSION=10.0
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
      WRITE(II,2002) NASLOV
      WRITE(II,5200) VERSION
      WRITE(II,5100) IND
C
      NULA=0
      ZERO=0.
C
      EX=0.
      EY=0.
      EZ=0.
      GX=0.
      GY=0.
      GZ=0.
      VX=0.
      VY=0.
      VZ=0.

      IND=-1
      ITYP=601
      ITYPE1=-601
      IJEDAN=1
      IDESET=10
      IDPET=25
      IDVESTA=200
      IPEDESET=50
      ISEDAM=70
      ISHAPEBEAM=5
      ISPET=75
      APRESEK=1.0
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
c      if (nummat.gt.120) stop 'nummat.gt.120 - color - TRGRMAT'
C stampanje materijala
      DO I=1,NUMMAT
        if (nummat.gt.120) then
          ICOL=I-120
        else
          ICOL=I
         endif
c         ICOL=124-(I-1)*2
c         ICOL=120-(I-1)*3
         WRITE(II,5600) I,ITYPE1,ICOL,NULA,NULA,I,NULA
         WRITE(II,2000) I
         WRITE(II,5600) IDESET
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) IDPET
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA
         WRITE(II,5600) IDVESTA
         WRITE(II,5500) EX,EY,EZ,GX,GY,GZ,VX,VY,VZ,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,AK(3,1,I),AK(3,2,I),AK(3,3,I),ZERO,
     +                  ZERO,ZERO,AK(3,4,I),GAMA
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5600) IPEDESET
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) ISEDAM
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
      ENDDO
      WRITE(II,5100) IND
C
        IND=-1
        ITYP=402
        WRITE(II,5100) IND
        WRITE(II,5100) ITYP
      DO J=1,numeltip
        SELECT CASE (eltypes(J))
        CASE (12)
          IPRO=5
        CASE (13)
          IPRO=37
        CASE (23)
          IPRO=19
        CASE (24)
          IPRO=19
        CASE (26)
          IPRO=20
        CASE (28)
          IPRO=20
        CASE (34)
          IPRO=25
        CASE (310)
          IPRO=26
        CASE (38)
          IPRO=25
        CASE (320)
          IPRO=26
        CASE DEFAULT
          WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',eltypes(J)
        END SELECT
C stampanje propertija - stampa za svaki properti materijal, tj. bro properija=NUMMAT*numeltip
        DO I=1,NUMMAT
         ICOL=I
c         ICOL=124-(I-1)*2
c         ICOL=120-(I-1)*3
C broj propertija = ltypes(J)+I
         WRITE(II,5600) eltypes(J)*100+I,ICOL,I,IPRO,I,NULA
         WRITE(II,2001) I
C od 3 linije  se razlikuje properti u zavisnosti od tipa konacnog elementa
              SELECT CASE (IPRO)
C      grede bez i sa medjucvorovima
                 CASE (5, 37)
			    WRITE(II,5600) NULA,ISHAPEBEAM,IJEDAN,NULA
			    WRITE(II,5600) IDESET
			    WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
			    WRITE(II,5600) NULA,NULA
			    WRITE(II,5600) ISPET
			    WRITE(II,5500) APRESEK,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) APRESEK,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5600) NULA
			    WRITE(II,5600) NULA
C 
C 2D elementi
			CASE (19, 20)
			    WRITE(II,5600) NULA,NULA,NULA,NULA
			    WRITE(II,5600) IDESET
			    WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
			    WRITE(II,5600) NULA,NULA
			    WRITE(II,5600) ISPET
			    WRITE(II,5500) APRESEK,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5600) NULA
			    WRITE(II,5600) NULA
C
C 3D elementi			
			CASE (25, 26)
			    WRITE(II,5600) NULA,NULA,NULA,NULA
			    WRITE(II,5600) IDESET
			    WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
			    WRITE(II,5600) NULA,NULA
			    WRITE(II,5600) ISPET
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5600) NULA
			    WRITE(II,5600) NULA
			CASE DEFAULT
			END SELECT
        ENDDO
C        WRITE(II,5100) IND
      enddo
        WRITE(II,5100) IND
      RETURN
C
 5100 FORMAT(I5)
 5500 FORMAT(10(1PE13.5))
 5600 FORMAT(10I5)
 5700 FORMAT(I4,1(1PE13.5),3I4)
 2000 FORMAT(' MAT',I4)
 2001 FORMAT(' PRO',I4)
 2002 FORMAT(A80)
 5200 FORMAT(F4.1,',')
      END
C======================================================================
C
C======================================================================
      SUBROUTINE VTGRMATP(AK,NUMMAT,II,NASLOV,GAMA,IVER)
      USE PRESEK
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO PRINT MATERIALS AND PROPERTIES IN NEU FILE
C .
C ......................................................................
C
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /DJERDAP/ IDJERDAP,ISPRESEK
      CHARACTER*80 NASLOV
      CHARACTER*30 IMEPR
      DIMENSION AK(3,5,*)
C
      IND=-1
      ITYP=100  
      VERSION=10.0
      IF(IVER.eq.44) VERSION=4.4
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
      WRITE(II,2002) NASLOV
      WRITE(II,5200) VERSION
      WRITE(II,5100) IND
C
      NULA=0
      ZERO=0.
C
      EX=0.
      EY=0.
      EZ=0.
      GX=0.
      GY=0.
      GZ=0.
      VX=0.
      VY=0.
      VZ=0.
C
      IND=-1
      IF(IVER.eq.10) THEN
!  verzija 10.0
        ITYP=601
        ITYPE1=-601
      ELSEIF(IVER.eq.44) THEN
        ITYP=401
        ITYPE1=1
      ENDIF
      IJEDAN=1
      IDESET=10
      IDPET=25
      IDVESTA=200
      IPEDESET=50
      ISEDAM=70
      ICPET=45
      IDTRI=93
      ISHAPEBEAM=5
      ISPET=75
      APRESEK=1.0
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
c      if (nummat.gt.120) stop 'nummat.gt.120 - color - TRGRMAT'
C stampanje materijala
      DO I=1,NUMMAT
        if (nummat.gt.120) then
          ICOL=I-120
        else
          ICOL=I
         endif
         IF(IVER.EQ.44) THEN
         WRITE(II,5600) I,ICOL,NULA,I,NULA
         WRITE(II,2000) I
         WRITE(II,5500) EX,EY,EZ
         WRITE(II,5500) GX,GY,GZ
         WRITE(II,5500) VX,VY,VZ
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO
         WRITE(II,5500) AK(3,1,I),AK(3,2,I),AK(3,3,I),ZERO,ZERO
         WRITE(II,5500) ZERO
         WRITE(II,5500) AK(3,4,I),GAMA,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA
         WRITE(II,5500) ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO
         WRITE(II,5600) NULA,NULA,NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA
         ELSEIF(IVER.EQ.10) THEN
c         ICOL=124-(I-1)*2
c         ICOL=120-(I-1)*3
         WRITE(II,5600) I,ITYPE1,ICOL,NULA,NULA,I,NULA
         WRITE(II,2000) I
         WRITE(II,5600) IDESET
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) IDPET
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA
         WRITE(II,5600) IDVESTA
         WRITE(II,5500) EX,EY,EZ,GX,GY,GZ,VX,VY,VZ,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,AK(3,1,I),AK(3,2,I),AK(3,3,I),ZERO,
     +                  ZERO,ZERO,AK(3,4,I),GAMA
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     +                  ZERO,ZERO
         WRITE(II,5600) IPEDESET
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) ISEDAM
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
         WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                  NULA,NULA
       ENDIF
      ENDDO
      WRITE(II,5100) IND
C
        IND=-1
        ITYP=402
        WRITE(II,5100) IND
        WRITE(II,5100) ITYP
!       DO J=1,numeltip
!         SELECT CASE (eltypes(J))
!         CASE (12)
!           IPRO=5
!         CASE (13)
!           IPRO=37
!         CASE (23)
          IPRO=19
!         CASE (24)
!           IPRO=19
!         CASE (26)
!           IPRO=20
!         CASE (28)IMEPR
!           IPRO=20
!         CASE (34)
!           IPRO=25
!         CASE (310)
!           IPRO=26
!         CASE (38)
!           IPRO=25
!         CASE (320)
!           IPRO=26
!         CASE DEFAULT
!           WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',eltypes(J)
!         END SELECT
C stampanje propertija - stampa za svaki properti materijal, tj. bro properija=NUMMAT*numeltip
        NUMMATT=NUMMAT
        IF(IPAKT.EQ.0.AND.ISPRESEK.EQ.-2) NUMMATT=13
        DO I=1,NUMMATT
         ICOL=I
c         ICOL=124-(I-1)*2
c         ICOL=120-(I-1)*3
C broj propertija = eltypes(J)+I
C         ieltypesp=2300+I
         ieltypesp=I
         WRITE(II,5600) ieltypesp,ICOL,I,IPRO,I,NULA
! stampanje naziva propertija         
!  nazivi propertija za filtraciju
        IF (IPAKT.EQ.0.AND.IDJERDAP.NE.0)THEN
         SELECT CASE (ICOL)
          CASE  (1)
             IMEPR="BE-L8"
          CASE  (2)
             IMEPR="BE-L9"
          CASE  (3)
             IMEPR='BE-L10'
          CASE  (4)
             IMEPR='BE-L11'
          CASE  (5)
             IMEPR='BE-L12'
          CASE  (6)     
             IMEPR='BE-L13'
          CASE  (7)
             IMEPR='BE-L14'
          CASE  (8) 
             IMEPR='BE-TPH'
          CASE  (9)   
             IMEPR='IZ-L8P'
          CASE (10)  
             IMEPR='IZ-L9P    '
          CASE (11)   
             IMEPR='IZ-L10P   '
          CASE (12)  
             IMEPR='IZ-L11P   '
          CASE (13)  
             IMEPR='IZ-L12P   '
          CASE (14)  
             IMEPR='IZ-L13P   '
          CASE (15)  
             IMEPR='IZ-L14P   '
          CASE (16)  
             IMEPR='IZ-L8D    '
          CASE (17)  
             IMEPR='IZ-L9D    '
          CASE (18)  
             IMEPR='IZ-L10D   '
          CASE (19)  
             IMEPR='IZ-L11D   '
          CASE (20)  
             IMEPR='IZ-L12D   '
          CASE (21)  
             IMEPR='IZ-L13D   '
          CASE (22)  
             IMEPR='IZ-L14D   '
          CASE (23)  
             IMEPR='SM-1      '
          CASE (24)  
             IMEPR='SM-1P     '
          CASE (25)  
             IMEPR='SM-2      '
          CASE (26)   
             IMEPR='SM-2P     '
          CASE (27)  
             IMEPR='SM-3      '
          CASE (28)   
             IMEPR='SM-4      '
          CASE (29)   
             IMEPR='SM-5      '
          CASE (30)   
             IMEPR='GT        '
          CASE (31) 
             IMEPR='DTG-Du-8/9'
          CASE (32) 
             IMEPR='DTG-Su-8/9'
          CASE (33)  
             IMEPR='DTG-Du-9    '
          CASE (34)  
             IMEPR='DTG-Su-9    '
          CASE (35)  
             IMEPR='DTG-Du-9/10 '
          CASE (36)  
             IMEPR='DTG-Su-9/10 '
          CASE (37)  
             IMEPR='DTG-Du-10   '
          CASE (38)  
             IMEPR='DTG-Su-10   '
          CASE (39)  
             IMEPR='DTG-Du-11   '
          CASE (40)  
             IMEPR='DTG-Su-11   '
          CASE (41)  
             IMEPR='DTG-Du-11/12'
          CASE (42)  
             IMEPR='DTG-Su-11/12'
          CASE (43) 
             IMEPR='DTG-Du-12   '
          CASE (44)  
             IMEPR='DTG-Su-12   '
          CASE (45)  
             IMEPR='DTG-Su-12/13'
          CASE (46)  
             IMEPR='BE-L1       '
          CASE (47)  
             IMEPR='BE-L2            '
          CASE (48) 
             IMEPR='BE-L3   '
          CASE (49) 
             IMEPR='BE-L4   '
          CASE (50) 
             IMEPR='BE-L5   '
          CASE (51)  
             IMEPR='BE-L6   '
          CASE (52)  
             IMEPR='BE-L7   '
          CASE (53)  
             IMEPR='BE-RZ   '
          CASE (54)  
             IMEPR='BE-RZB  '
          CASE (55)   
             IMEPR='IZ-L1P  '
          CASE (56)  
             IMEPR='IZ-L2P  '
          CASE (57)   
             IMEPR='IZ-L3P  '
          CASE (58)   
             IMEPR='IZ-L4P  '
          CASE (59)  
             IMEPR='IZ-L5P  '
          CASE (60)  
             IMEPR='IZ-L6P  '
          CASE (61)  
             IMEPR='IZ-L7P  '
          CASE (62)  
             IMEPR='IZ-L1D  '
          CASE (63)  
             IMEPR='IZ-L2D  '
          CASE (64)  
             IMEPR='IZ-L3D  '
          CASE (65)  
             IMEPR='IZ-L4D  '
          CASE (66)   
             IMEPR='IZ-L5D  '
          CASE (67)   
             IMEPR='IZ-L6D  '
          CASE (68)  
             IMEPR='IZ-L7D  '
          CASE (69)  
             IMEPR='DTG-S0-1'
          CASE (70)  
             IMEPR='DTG-Du-1'
          CASE (71)  
             IMEPR='DTG-Dn-1'
          CASE (72)  
             IMEPR='DTG-Su-1'
          CASE (73) 
             IMEPR='DTG-Du-1/2  '
          CASE (74)  
             IMEPR='DTG-Dn-1/2  '
          CASE (75)  
             IMEPR='DTG-Su-1/2  '
          CASE (76)  
             IMEPR='DTG-Du-2    '
          CASE (77)  
             IMEPR='DTG-Dn-2    '
          CASE (78)  
             IMEPR='DTG-Su-2    '
          CASE (79)   
             IMEPR='DTG-Du-2/3  '
          CASE (80)   
             IMEPR='DTG-Dn-2/3  '
          CASE (81)  
             IMEPR='DTG-Su-2/3  '
          CASE (82)   
             IMEPR='DTG-Du-3    '
          CASE (83)  
             IMEPR='DTG-Dn-3    '
          CASE (84)  
             IMEPR='DTG-Su-3    '
          CASE (85)  
             IMEPR='DTG-Du-3/4  '
          CASE (86)   
             IMEPR='DTG-Dn-3/4  '
          CASE (87)   
             IMEPR='DTG-Su-3/4  '
          CASE (88)  
             IMEPR='DTG-Du-4    '
          CASE (89)  
             IMEPR='DTG-Su-4    '
          CASE (90)  
             IMEPR='DTG-Du-4/5  '
          CASE (91)  
             IMEPR='DTG-Dn-4/5  '
          CASE (92)  
             IMEPR='DTG-Su-4/5  '
          CASE (93)  
             IMEPR='DTG-Du-5    '
          CASE (94)  
             IMEPR='DTG-Dn-5    '
          CASE (95) 
             IMEPR='DTG-Su-5    '
          CASE (96)   
             IMEPR='DTG-Du-5/6  '
          CASE (97)   
             IMEPR='DTG-Dn-5/6  '
          CASE (98)  
             IMEPR='DTG-Su-5/6  '
          CASE (99)  
             IMEPR='DTG-Du-6    '
          CASE(100)   
             IMEPR='DTG-Su-6    '
          CASE(101)   
             IMEPR='DTG-Du-6/7  '
          CASE(102)  
             IMEPR='DTG-Dn-6/7  '
          CASE(103)  
             IMEPR='DTG-Su-6/7  '
          CASE(104)  
             IMEPR='DTG-Du-7    '
          CASE(105)  
             IMEPR='DTG-Dn-7    '
          CASE(106)   
             IMEPR='DTG-Su-7    '
          CASE(107)  
             IMEPR='DTG-Du-7/8  '
          CASE(108)  
             IMEPR='DTG-Dn-7/8  '
          CASE(109)  
             IMEPR='DTG-Su-7/8  '
          CASE(110)  
             IMEPR='DTG-DG-L1-6 '
          CASE(111)  
             IMEPR='DTG-DG-L7   '
          CASE(112)  
             IMEPR='BE-S3'
          CASE(113)  
             IMEPR='BE-S2'
          CASE(114)   
             IMEPR='BE-S1'
          CASE(115)   
             IMEPR='BE-MB'
          CASE(116)  
             IMEPR='BE-BP    '
          CASE(117)  
             IMEPR='BE-RZ    '
          CASE(118)  
             IMEPR='IZ-S3-1P '
          CASE(119)  
             IMEPR='IZ-S3-2P '
          CASE(120)   
             IMEPR='IZ-S3-3P '
          CASE(121)  
             IMEPR='IZ-S3-4P '
          CASE(122)  
             IMEPR='IZ-S2-1P '
          CASE(123)  
             IMEPR='IZ-S2-2P '
          CASE(124)  
             IMEPR='IZ-S2-3P '
          CASE(125)  
             IMEPR='IZ-S2-4P '
          CASE(126) 
             IMEPR='IZ-S1-1P '
          CASE(127) 
             IMEPR='IZ-S1-2P '
          CASE(128) 
             IMEPR='IZ-S1-3P '
          CASE(129)  
             IMEPR='IZ-S1-4P '
          CASE(130)  
             IMEPR='IZ-MB-1P '
          CASE(131)   
             IMEPR='IZ-MB-2P '
          CASE(132)  
             IMEPR='IZ-MB-3P '
          CASE(133)  
             IMEPR='IZ-S3-1D '
          CASE(134)  
             IMEPR='IZ-S3-2D '
          CASE(135)  
             IMEPR='IZ-S3-3D '
          CASE(136)   
             IMEPR='IZ-S3-4D '
          CASE(137)  
             IMEPR='IZ-S2-1D '
          CASE(138)   
             IMEPR='IZ-S2-2D '
          CASE(139)   
             IMEPR='IZ-S2-3D '
          CASE(140)  
             IMEPR='IZ-S2-4D '
          CASE(141)  
             IMEPR='IZ-S1-1D '
          CASE(142)  
             IMEPR='IZ-S1-2D '
          CASE(143)   
             IMEPR='IZ-S1-3D '
          CASE(144)  
             IMEPR='IZ-S1-4D '
          CASE(145) 
             IMEPR='IZ-MB-1D '
          CASE(146)  
             IMEPR='IZ-MB-2D '
          CASE(147) 
             IMEPR='IZ-MB-3D '
          CASE(148)  
             IMEPR='IB-B1    '
          CASE(149)  
             IMEPR='IB-B2    '
          CASE(150)  
             IMEPR='IB-B3    '
          CASE(151)  
             IMEPR='IB-B4    '
          CASE(152)  
             IMEPR='IB-B5    '
          CASE(153)  
             IMEPR='IB-B6    '
          CASE(154)  
             IMEPR='IB-B7    '
          CASE(155)  
             IMEPR='IB-B8    '
          CASE(156)   
             IMEPR='IB-B9   '
          CASE(157)  
             IMEPR='IB-B10   '
          CASE(158)  
             IMEPR='IB-B11   '
          CASE(159)  
             IMEPR='IB-B12   '
          CASE(160)  
             IMEPR='IB-B13   '
          CASE(161) 
             IMEPR='IB-B14   '
          CASE(162)  
             IMEPR='IB-B15   '
          CASE(163)  
             IMEPR='IB-B16   '
          CASE(164)  
             IMEPR='IB-B17   '
          CASE(165)  
             IMEPR='IB-B18   '
          CASE(166)  
             IMEPR='IB-B19   '
          CASE(167)   
             IMEPR='IB-B20   '
          CASE(168)  
             IMEPR='IB-B21   '
          CASE(169)  
             IMEPR='IB-B22   '
          CASE(170)  
             IMEPR='IB-B23   '
          CASE(171) 
             IMEPR='IB-B24   '
          CASE(172)  
             IMEPR='IB-B25   '
          CASE(173)  
             IMEPR='IB-B26   '
          CASE(174)  
             IMEPR='IB-B27   '
          CASE(175)  
             IMEPR='IB-B28   '
          CASE(176)   
             IMEPR='IB-B29   '
          CASE(177)  
             IMEPR='IB-B30   '
          CASE(178)  
             IMEPR='IB-B31   '
          CASE(179)  
             IMEPR='IB-B32   '
          CASE(180)  
             IMEPR='IB-B33   '
          CASE(181)  
             IMEPR='IB-B34   '
          CASE(182)  
             IMEPR='IB-B35   '
          CASE(183)  
             IMEPR='IB-B36   '
          CASE(184)  
             IMEPR='FP-7/8-IG/24   '
          CASE(185)  
             IMEPR='FP-7/8-Uz24/30 '
          CASE(186)  
             IMEPR='FP-7/8-Uz30/40 '
          CASE(187)   
             IMEPR='FP-7/8-Uz40/50 '
          CASE(188)   
             IMEPR='FP-7/8-Uz50/OB '
          CASE(189)  
             IMEPR='FP-7/8-OB/S    '
          CASE(190)  
             IMEPR='FP-7/8-P/S-G   '
          CASE(191)  
             IMEPR='FP-7/8-S       '
          CASE(192)  
             IMEPR='FP-7/8-P/S-D   '
          CASE(193)  
             IMEPR='FS-7/8-S/P-G   '
          CASE(194)  
             IMEPR='FS-7/8-Gore    '
          CASE(195)   
             IMEPR='FS-7/8-S/P-D   '
          CASE(196)   
             IMEPR='FS-7/8-Dole    '
          CASE(197)  
             IMEPR='FS-7/8-S       '
          CASE(198)   
             IMEPR='FP-8/9-IG/24   '
          CASE(199)   
             IMEPR='FP-8/9-Uz24/30 '
          CASE(200)   
             IMEPR='FP-8/9-Uz30/40 '
          CASE(201)  
             IMEPR='FP-8/9-Uz40/50 '
          CASE(202)  
             IMEPR='FP-8/9-Uz50/OB '
          CASE(203)  
             IMEPR='FP-8/9-OB/S    '
          CASE(204)  
             IMEPR='FP-8/9-P/S     '
          CASE(205)  
             IMEPR='FP-8/9-S       '
          CASE(206)  
             IMEPR='FS-8/9-S/P     '
          CASE(207)  
             IMEPR='FS-8/9-Gore    '
          CASE(208)  
             IMEPR='FS-8/9-Dole    '
          CASE(209)  
             IMEPR='FS-8/9-S       '
          CASE(210)  
             IMEPR='FP-9/10-IG/24  '
          CASE(211)  
             IMEPR='FP-9/10-Uz24/30'
          CASE(212)  
             IMEPR='FP-9/10-Uz30/40'
          CASE(213)   
             IMEPR='FP-9/10-Uz40/50'
          CASE(214)  
             IMEPR='FP-9/10-Uz50/OB'
          CASE(215)   
             IMEPR='FP-9/10-OB/S   '
          CASE(216)  
             IMEPR='FP-9/10-P/S    '
          CASE(217)   
             IMEPR='FP-9/10-S      '
          CASE(218)  
             IMEPR='FS-9/10-S/P    '
          CASE(219)  
             IMEPR='FS-9/10-Gore   '
          CASE(220)  
             IMEPR='FS-9/10-Dole   '
          CASE(221)  
             IMEPR='FS-9/10-S      '
          CASE(222)   
             IMEPR='FP-10/11-IG/24 '
          CASE(223)   
             IMEPR='FP-10/11-Uz24/30 '
          CASE(224) 
             IMEPR='FP-10/11-Uz30/40 '
          CASE(225) 
             IMEPR='FP-10/11-Uz40/50 '
          CASE(226) 
             IMEPR='FP-10/11-Uz50/OB '
          CASE(227)  
             IMEPR='FP-10/11-OB/S    '
          CASE(228)   
             IMEPR='FP-10/11-P/S     '
          CASE(229)  
             IMEPR='FP-10/11-S       '
          CASE(230) 
             IMEPR='FS-10/11-S/P     '
          CASE(231)  
             IMEPR='FS-10/11-Gore    '
          CASE(232)  
             IMEPR='FS-10/11-Dole    '
          CASE(233) 
             IMEPR='FS-10/11-S       '
          CASE(234) 
             IMEPR='FP-11/12-IG/24   '
          CASE(235) 
             IMEPR='FP-11/12-Uz24/30 '
          CASE(236)  
             IMEPR='FP-11/12-Uz30/40 '
          CASE(237)  
             IMEPR='FP-11/12-Uz40/50 '
          CASE(238)  
             IMEPR='FP-11/12-Uz50/OB '
          CASE(239) 
             IMEPR='FP-11/12-OB/S    '
          CASE(240)  
             IMEPR='FP-11/12-P/S     '
          CASE(241)  
             IMEPR='FP-11/12-S       '
          CASE(242) 
             IMEPR='FS-11/12-S/P     '
          CASE(243) 
             IMEPR='FS-11/12-Gore    '
          CASE(244)  
             IMEPR='FS-11/12-Dole    '
          CASE(245)  
             IMEPR='FS-11/12-S       '
          CASE(246)  
             IMEPR='FP-12/13-IG/24   '
          CASE(247)  
             IMEPR='FP-12/13-Uz24/30 '
          CASE(248)  
             IMEPR='FP-12/13-Uz30/40 '
          CASE(249)   
             IMEPR='FP-12/13-Uz40/50 '
          CASE(250)  
             IMEPR='FP-12/13-Uz50/OB '
          CASE(251)   
             IMEPR='FP-12/13-OB/S    '
          CASE(252)   
             IMEPR='FP-12/13-P/S-G   '
          CASE(253)  
             IMEPR='FP-12/13-P/S-D   '
          CASE(254)   
             IMEPR='FP-12/13-S       '
          CASE(255)   
             IMEPR='FS-12/13-S/P-G   '
          CASE(256)  
             IMEPR='FS-12/13-Gore    '
          CASE(257)   
             IMEPR='FS-12/13-S/P-D   '
          CASE(258)   
             IMEPR='FS-12/13-Dole    '
          CASE(259)  
             IMEPR='FS-12/13-S       '
          CASE(260)   
             IMEPR='FP-13/14-IG/24   '
          CASE(261)  
             IMEPR='FP-13/14-Uz24/30 '
          CASE(262)  
             IMEPR='FP-13/14-Uz30/40 '
          CASE(263)  
             IMEPR='FP-13/14-Uz40/50 '
          CASE(264)  
             IMEPR='FP-13/14-Uz50/OB '
          CASE(265)  
             IMEPR='FP-13/14-OB/S    '
          CASE(266)   
             IMEPR='FP-13/14-P/S     '
          CASE(267)  
             IMEPR='FS-13/14-S/P     '
          CASE(268)  
             IMEPR='FS-13/14-Gore    '
          CASE(269)  
             IMEPR='FS-13/14-Dole    '
          CASE(270)  
             IMEPR='L8-S-P/S         '
          CASE(271)  
             IMEPR='L9-S-P/S         '
          CASE(272)  
             IMEPR='L10-S-P/S        '
          CASE(273)  
             IMEPR='L11-S-P/S        '
          CASE(274)   
             IMEPR='L12-S-P/S        '
          CASE(275)   
             IMEPR='L13-S-P/S        '
          CASE(276)  
             IMEPR='L14-S-P/S        '
          CASE(277)  
             IMEPR='FP-14/15-IG/24   '
          CASE(278)  
             IMEPR='FP-14/15-Uz24/30 '
          CASE(279)  
             IMEPR='FP-14/15-Uz30/40 '
          CASE(280)  
             IMEPR='FP-14/15-Uz40/50 '
          CASE(281)  
             IMEPR='FP-14/15-Uz50/OB '
          CASE(282)  
             IMEPR='FP-14/15-OB/S    '
          CASE(283)  
             IMEPR='FP-14/15-P/S     '
          CASE(284)   
             IMEPR='FS-14/15-S/P     '
          CASE(285)   
             IMEPR='FS-14/15-Gore    '
          CASE(286)  
             IMEPR='FS-14/15-Dole    '
          CASE(287)  
             IMEPR='P4-U-RE '
          CASE(288)   
             IMEPR='P4-U-RA '
          CASE(289)  
             IMEPR='P5-U-RE '
          CASE(290)  
             IMEPR='P5-U-RA '
          CASE(291)  
             IMEPR='P6-U-RE '
          CASE(292)  
             IMEPR='P6-U-RA '
          CASE(293)  
             IMEPR='P7-U-RE '
          CASE(294)  
             IMEPR='P7-U-RA '
          CASE(295)   
             IMEPR='P8-U-RE '
          CASE(296)   
             IMEPR='P8-U-RA '
          CASE(297)  
             IMEPR='ZA-2D   '
          CASE(298)  
             IMEPR='PR-2D   '
          CASE(299) 
             IMEPR='FP-0/1-Uz24/30   '
          CASE(300) 
             IMEPR='FP-0/1-Uz30/40   '
          CASE(301)  
             IMEPR='FP-0/1-Uz40/50   '
          CASE(302)  
             IMEPR='FP-0/1-Uz50/OB   '
          CASE(303)  
             IMEPR='FP-0/1-OB/S      '
          CASE(304)  
             IMEPR='FS-0/1-Gore      '
          CASE(305)  
             IMEPR='FS-0/1-Niz1.8/31 '
          CASE(306)  
             IMEPR='FS-0/1-Dole      '
          CASE(307)  
             IMEPR='FP-0/1-Dole      '
          CASE(308)  
             IMEPR='FP-0/1-Uz1.8/24  '
          CASE(309)  
             IMEPR='FP-0/1-38        '
          CASE(310)  
             IMEPR='FP-0/1-34        '
          CASE(311)  
             IMEPR='FP-0/1-P/S       '
          CASE(312)  
             IMEPR='FS-0/1-S/P       '
          CASE(313)   
             IMEPR='FP-1/2-IG/24     '
          CASE(314)  
             IMEPR='FP-1/2-Uz24/30'
          CASE(315)  
             IMEPR='FP-1/2-Uz30/40'
          CASE(316)  
             IMEPR='FP-1/2-Uz40/50'
          CASE(317)   
             IMEPR='FP-1/2-Uz50/OB'
          CASE(318)  
             IMEPR='FP-1/2-OB/S   '
          CASE(319)  
             IMEPR='FP-1/2-P/S    '
          CASE(320)  
             IMEPR='FS-1/2-S/P    '
          CASE(321)   
             IMEPR='FS-1/2-Gore   '
          CASE(322)   
             IMEPR='FS-1/2-Dole   '
          CASE(323)  
             IMEPR='FP-1/2-S      '
          CASE(324) 
             IMEPR='FS-1/2-S      '
          CASE(325)  
             IMEPR='FP-2/3-IG/24  '
          CASE(326)  
             IMEPR='FP-2/3-Uz24/30'
          CASE(327)  
             IMEPR='FP-2/3-Uz30/40'
          CASE(328) 
             IMEPR='FP-2/3-Uz40/50'
          CASE(329)  
             IMEPR='FP-2/3-Uz50/OB'
          CASE(330)  
             IMEPR='FP-2/3-OB/S   '
          CASE(331) 
             IMEPR='FP-2/3-P/S    '
          CASE(332) 
             IMEPR='FP-2/3-S      '
          CASE(333)  
             IMEPR='FS-2/3-S/P    '
          CASE(334)  
             IMEPR='FS-2/3-Gore   '
          CASE(335)  
             IMEPR='FS-2/3-Dole   '
          CASE(336)  
             IMEPR='FS-2/3-S      '
          CASE(337)   
             IMEPR='FP-3/4-IG/24  '
          CASE(338)  
             IMEPR='FP-3/4-Uz24/30'
          CASE(339)   
             IMEPR='FP-3/4-Uz30/40'
          CASE(340)   
             IMEPR='FP-3/4-Uz40/50'
          CASE(341)  
             IMEPR='FP-3/4-Uz50/OB'
          CASE(342)  
             IMEPR='FP-3/4-OB/S   '
          CASE(343)  
             IMEPR='FP-3/4-P/S    '
          CASE(344) 
             IMEPR='FP-3/4-S      '
          CASE(345)  
             IMEPR='FS-3/4-S/P    '
          CASE(346)  
             IMEPR='FS-3/4-Gore   '
          CASE(347)  
             IMEPR='FS-3/4-Dole   '
          CASE(348)  
             IMEPR='FS-3/4-S      '
          CASE(349) 
             IMEPR='FP-4/5-IG/24  '
          CASE(350)  
             IMEPR='FP-4/5-Uz24/30'
          CASE(351)  
             IMEPR='FP-4/5-Uz30/40'
          CASE(352)  
             IMEPR='FP-4/5-Uz40/50'
          CASE(353)   
             IMEPR='FP-4/5-Uz50/OB'
          CASE(354)  
             IMEPR='FP-4/5-OB/S   '
          CASE(355)  
             IMEPR='FP-4/5-P/S    '
          CASE(356)  
             IMEPR='FP-4/5-S   '
          CASE(357)  
             IMEPR='FS-4/5-S/P '
          CASE(358)  
             IMEPR='FS-4/5-Gore'
          CASE(359) 
             IMEPR='FS-4/5-Dole'
          CASE(360)  
             IMEPR='FS-4/5-S   '
          CASE(361)  
             IMEPR='FP-5/6-IG/24   '
          CASE(362)  
             IMEPR='FP-5/6-Uz24/30 '
          CASE(363)   
             IMEPR='FP-5/6-Uz30/40 '
          CASE(364)  
             IMEPR='FP-5/6-Uz40/50 '
          CASE(365)  
             IMEPR='FP-5/6-Uz50/OB '
          CASE(366)  
             IMEPR='FP-5/6-OB/S    '
          CASE(367)  
             IMEPR='FP-5/6-P/S     '
          CASE(368)  
             IMEPR='FP-5/6-S       '
          CASE(369)  
             IMEPR='FS-5/6-S/P     '
          CASE(370)    
             IMEPR='FS-5/6-Gore    '
          CASE(371)  
             IMEPR='FS-5/6-Dole    '
          CASE(372)  
             IMEPR='FS-5/6-S       '
          CASE(373)  
             IMEPR='FP-6/7-IG/24   '
          CASE(374)  
             IMEPR='FP-6/7-Uz24/30 '
          CASE(375)  
             IMEPR='FP-6/7-Uz30/40 '
          CASE(376) 
             IMEPR='FP-6/7-Uz40/50 '
          CASE(377) 
             IMEPR='FP-6/7-Uz50/OB '
          CASE(378)  
             IMEPR='FP-6/7-OB/S    '
          CASE(379)  
             IMEPR='FP-6/7-P/S     '
          CASE(380)  
             IMEPR='FP-6/7-S       '
          CASE(381)  
             IMEPR='FS-6/7-S/P     '
          CASE(382)  
             IMEPR='FS-6/7-Gore    '
          CASE(383)  
             IMEPR='FS-6/7-Dole    '
          CASE(384)  
             IMEPR='FS-6/7-S       '
          CASE(385)  
             IMEPR='L1-S-P/S       '
          CASE(386) 
             IMEPR='L2-S-P/S       '
          CASE(387)  
             IMEPR='L3-S-P/S       '
          CASE(388)  
             IMEPR='L4-S-P/S       '
          CASE(389)  
             IMEPR='L5-S-P/S       '
          CASE(390)  
             IMEPR='L6-S-P/S       '
          CASE(391) 
             IMEPR='L7-S-P/S       '
          CASE(392)  
             IMEPR='P1-U-RE        '
          CASE(393)  
             IMEPR='P1-U-RA        '
          CASE(394)  
             IMEPR='P2-U-RE        '
          CASE(395)  
             IMEPR='P2-U-RA        '
          CASE(396)  
             IMEPR='P3-U-RE        '
          CASE(397) 
             IMEPR='P3-U-RA        '
          CASE(398)   
             IMEPR='F-II/III-U42/70.5'
          CASE(399)   
             IMEPR='F-II/III-U25.2/42'
          CASE(400)  
             IMEPR='F-II/III-U18.5/25.2'
          CASE(401) 
             IMEPR='F-II/III-IG      '
          CASE(402) 
             IMEPR='F-II/III-U10.9/19.4'
          CASE(403) 
             IMEPR='F-II/III-SG      '
          CASE(404)   
             IMEPR='F-II/III-MG      '
          CASE(405) 
             IMEPR='F-II/III-Dole    '
          CASE(406) 
             IMEPR='F-II/III-DG      '
          CASE(407)  
             IMEPR='F-II/III-N       '
          CASE(408)  
             IMEPR='F-II/III-U22/48  '
          CASE(409)  
             IMEPR='F-II/III-Gore    '
          CASE(410)  
             IMEPR='F-I/II-U42/70.5  '
          CASE(411)  
             IMEPR='F-I/II-U25.2/42  '
          CASE(412)  
             IMEPR='F-I/II-U18.5/25.2'
          CASE(413)  
             IMEPR='F-I/II-IG        '
          CASE(414)  
             IMEPR='F-I/II-10.9/19.4 '
          CASE(415)  
             IMEPR='F-I/II-SG        '
          CASE(416)  
             IMEPR='F-I/II-MG        '
          CASE(417)  
             IMEPR='F-I/II-Dole      '
          CASE(418)  
             IMEPR='F-I/II-DG        '
          CASE(419) 
             IMEPR='F-I/II-N         '
          CASE(420)
             IMEPR='F-I/II-U22/48    '
          CASE(421) 
             IMEPR='F-I/II-Gore      '
          CASE(422) 
             IMEPR='F-I/MB-U42/70.5  '
          CASE(423)  
             IMEPR='F-I/MB-U27.9/42  '
          CASE(424) 
             IMEPR='F-I/MB-U18.5/27.9'
          CASE(425)  
             IMEPR='F-I/MB-IG        '
          CASE(426) 
             IMEPR='F-I/MB-10.9/19.4 '
          CASE(427) 
             IMEPR='F-I/MB-SG        '
          CASE(428)  
             IMEPR='F-I/MB-MG        '
          CASE(429)  
             IMEPR='F-I/MB-Dole      '
          CASE(430) 
             IMEPR='F-I/MB-DG        '
          CASE(431)  
             IMEPR='F-I/MB-N         '
          CASE(432)  
             IMEPR='F-I/MB-N5.6/23.6 '
          CASE(433)  
             IMEPR='F-I/MB-N23.6/48.1'
          CASE(434)  
             IMEPR='IZ-I/MB-1        '
          CASE(435)  
             IMEPR='IZ-I/MB-2        '
          CASE(436)    
             IMEPR='IZ-I/MB-3        '
          CASE(437)  
             IMEPR='IZ-I/MB-4        '
          CASE(438)  
             IMEPR='IZ-I/MB-5        '
          CASE(439)   
             IMEPR='IZ-I/MB-6        '
          CASE(440)  
             IMEPR='F-MB/BP-IG       '
          CASE(441)  
             IMEPR='F-MB/BP-U42/72.5 '
          CASE(442)   
             IMEPR='F-MB/BP-U15.6/42 '
          CASE(443)   
             IMEPR='F-MB/BP-G        '
          CASE(444)  
             IMEPR='F-MB/BP-Dole     '
          CASE(445)  
             IMEPR='F-MB/BP-N23.7/43.4'
          CASE(446)   
             IMEPR='F-MB/BP-Gore'
          CASE(447)   
             IMEPR='E-Z-MB-S    '
          CASE(448)  
             IMEPR='E-Z-S-S     '
          CASE(449)  
             IMEPR='L8-6B       '
          CASE(450)  
             IMEPR='L8-5B       '
          CASE(451) 
             IMEPR='L8-4B       '
          CASE(452) 
             IMEPR='L8-3B       '
          CASE(453)  
             IMEPR='L8-2B       '
          CASE(454)  
             IMEPR='L8-1B       '
          CASE(455)  
             IMEPR='D-8-2B      '
          CASE(456)  
             IMEPR='D-8-1B      '
          CASE(457)   
             IMEPR='V-8-1B      '
          CASE(458)    
             IMEPR='V-8-2B      '
          CASE(459)   
             IMEPR='V-8-3B      '
          CASE(460)  
             IMEPR='Su-8/9-2B   '
          CASE(461) 
             IMEPR='L9-4B       '
          CASE(462) 
             IMEPR='L9-3B       '
          CASE(463)  
             IMEPR='L9-2B       '
          CASE(464)  
             IMEPR='L9-1B       '
          CASE(465)  
             IMEPR='D-9-3B      '
          CASE(466) 
             IMEPR='D-9-2B      '
          CASE(467)  
             IMEPR='D-9-1B      '
          CASE(468)  
             IMEPR='Su-8/9-1B   '
          CASE(469) 
             IMEPR='Su-9-2B     '
          CASE(470) 
             IMEPR='Su-9-1B     '
          CASE(471) 
             IMEPR='Su-9/10-2B  '
          CASE(472)  
             IMEPR='L10-4B      '
          CASE(473)  
             IMEPR='L10-3B      '
          CASE(474)  
             IMEPR='L10-2B      '
          CASE(475) 
             IMEPR='L10-1B      '
          CASE(476) 
             IMEPR='D-10-3B     '
          CASE(477)  
             IMEPR='D-10-2B     '
          CASE(478)  
             IMEPR='D-10-1B     '
          CASE(479) 
             IMEPR='Su-9/10-1B  '
          CASE(480) 
             IMEPR='Su-10-2B    '
          CASE(481) 
             IMEPR='Su-10-1B    '
          CASE(482) 
             IMEPR='L11-3B      '
          CASE(483)  
             IMEPR='L11-2B      '
          CASE(484) 
             IMEPR='L11-1B      '
          CASE(485)  
             IMEPR='Du-11-1B    '
          CASE(486) 
             IMEPR='Du-11-2B    '
          CASE(487)  
             IMEPR='Du-11-3B    '
          CASE(488)  
             IMEPR='Du-11-4B    '
          CASE(489) 
             IMEPR='Du-11-5B    '
          CASE(490)  
             IMEPR='D-11-4B     '
          CASE(491)  
             IMEPR='D-11-3B     '
          CASE(492) 
             IMEPR='D-11-2B     '
          CASE(493) 
             IMEPR='D-11-1B     '
          CASE(494) 
             IMEPR='Su-11-1B    '
          CASE(495)  
             IMEPR='Su-11/12-3B '
          CASE(496)  
             IMEPR='Su-11/12-2B '
          CASE(497)  
             IMEPR='L12-4B      '
          CASE(498)  
             IMEPR='L12-3B      '
          CASE(499)  
             IMEPR='L12-2B      '
          CASE(500) 
             IMEPR='L12-1B      '
          CASE(501) 
             IMEPR='D-12-3B     '
          CASE(502)   
             IMEPR='D-12-2B     '
          CASE(503)  
             IMEPR='D-12-1B     '
          CASE(504)  
             IMEPR='Su-11/12-1B '
          CASE(505)  
             IMEPR='Su-12-2B    '
          CASE(506)  
             IMEPR='Su-12-1B    '
          CASE(507) 
             IMEPR='Su-12/13-1B '
          CASE(508) 
             IMEPR='L7-4B'
          CASE(509)   
             IMEPR='L7-3B'
          CASE(510)  
             IMEPR='L7-2B'
          CASE(511)  
             IMEPR='L7-1B'
          CASE(512) 
             IMEPR='B1   '
          CASE(513) 
             IMEPR='B2   '
          CASE(514)  
             IMEPR='B3   '
          CASE(515)  
             IMEPR='B4   '
          CASE(516)   
             IMEPR='B5   '
          CASE(517) 
             IMEPR='B6   '
          CASE(518)  
             IMEPR='B7   '
          CASE(519)   
             IMEPR='B8   '
          CASE(520)  
             IMEPR='B9   '
          CASE(521)  
             IMEPR='B10  '
          CASE(522)  
             IMEPR='B11  '
          CASE(523)  
             IMEPR='B12  '
          CASE(524)  
             IMEPR='B13  '
          CASE(525)  
             IMEPR='B14  '
          CASE(526)  
             IMEPR='B15  '
          CASE(527)  
             IMEPR='B16  '
          CASE(528)  
             IMEPR='B17  '
          CASE(529)  
             IMEPR='B18  '
          CASE(530)  
             IMEPR='B19  '
          CASE(531)  
             IMEPR='B20  '
          CASE(532)  
             IMEPR='B21  '
          CASE(533)  
             IMEPR='B22  '
          CASE(534)  
             IMEPR='B23  '
          CASE(535) 
             IMEPR='B24  '
          CASE(536)  
             IMEPR='B25  '
          CASE(537) 
             IMEPR='B26  '
          CASE(538)  
             IMEPR='B27  '
          CASE(539)  
             IMEPR='B28  '
          CASE(540)  
             IMEPR='B29  '
          CASE(541)   
             IMEPR='B30  '
          CASE(542)  
             IMEPR='B31  '
          CASE(543)   
             IMEPR='B32  '
          CASE(544)  
             IMEPR='B33  '
          CASE(545) 
             IMEPR='B34  '
          CASE(546) 
             IMEPR='B35  '
          CASE(547)  
             IMEPR='B36  '
         CASE DEFAULT
         END SELECT
        ELSEIF (IPAKT.EQ.1.AND.IDJERDAP.NE.0)THEN
         SELECT CASE (I)
           CASE  (1)
              IMEPR='beton'
           CASE  (2)
              IMEPR='SM-1 '
           CASE  (3)
              IMEPR="SM-1'"
           CASE  (4) 
              IMEPR='SM-2 '
           CASE  (5) 
              IMEPR="SM-2'"
           CASE  (6)
              IMEPR='SM-3 '
           CASE  (7) 
              IMEPR='SM-4 '
           CASE  (8) 
              IMEPR='SM-5 '
           CASE  (9) 
              IMEPR='injekciona zavesa '         
         CASE DEFAULT
         END SELECT
        ELSEIF (IPAKT.EQ.0.AND.ISPRESEK.EQ.-2)THEN
!          write(*,*) 'Vlasina properti'
         SELECT CASE (I)
          CASE  (1)
             IMEPR="Glineno jezgro"
          CASE  (2)
             IMEPR="Glineni sut 1"
          CASE  (3)
             IMEPR='Glineni sut 2'
          CASE  (4)
             IMEPR='Kameni sut 1'
          CASE  (5)
             IMEPR='Beton'
          CASE  (6)     
             IMEPR='Injekciona zavesa'
          CASE  (7)
             IMEPR='Stenska masa'
          CASE  (8) 
             IMEPR='Malo sispucali skriljci'
          CASE  (9)   
             IMEPR='Drenaza'
          CASE (10)  
             IMEPR='Drenaza - kameni nabacaj'
          CASE (11)   
             IMEPR='Deponija'
          CASE (12)  
             IMEPR='Nasip'
          CASE (13)  
             IMEPR='Kameni sut 2'
!           CASE  (1)
!              IMEPR="K-J-1"
!           CASE  (2)
!              IMEPR="K-J-2"
!           CASE  (3)
!              IMEPR='K-J-3'
!           CASE  (4)
!              IMEPR='K-J-4'
!           CASE  (5)
!              IMEPR='K-J-5'
!           CASE  (6)     
!              IMEPR='K-J-6'
!           CASE  (7)
!              IMEPR='K-J-7'
!           CASE  (8) 
!              IMEPR='K-J-8'
!           CASE  (9)   
!              IMEPR='K-J-9'
!           CASE (10)  
!              IMEPR='K-J-10'
!           CASE (11)   
!              IMEPR='K-J-11'
!           CASE (12)  
!              IMEPR='K-J-12'
!           CASE (13)  
!              IMEPR='K-G1'
!           CASE (14)  
!              IMEPR='K-G2'
!           CASE (15)  
!              IMEPR='K-K1'
!           CASE (16)  
!              IMEPR='K-K2'
!           CASE (17)  
!              IMEPR='K-D'
!           CASE (18)  
!              IMEPR='K-N'
!           CASE (19)  
!              IMEPR='K-B-IG'
!           CASE (20)  
!              IMEPR='K-B-OT'
!           CASE (21)  
!              IMEPR='K-IZ-P1'
!           CASE (22)  
!              IMEPR='K-IZ-P2'
!           CASE (23)  
!              IMEPR='K-IZ-P3'
!           CASE (24)  
!              IMEPR='K-IZ-P4'
!           CASE (25)  
!              IMEPR='K-IZ-P5'
!           CASE (26)   
!              IMEPR='K-IZ-P6'
!           CASE (27)  
!              IMEPR='K-IZ-P7'
!           CASE (28)   
!              IMEPR='K-IZ-P8'
!           CASE (29)   
!              IMEPR='K-IZ-P9'
!           CASE (30)   
!              IMEPR='K-IZ-P10'
!           CASE (31) 
!              IMEPR='K-IZ-P11'
!           CASE (32) 
!              IMEPR='K-IZ-P12'
!           CASE (33)  
!              IMEPR='K-IZ-D1'
!           CASE (34)  
!              IMEPR='K-IZ-D2'
!           CASE (35)  
!              IMEPR='K-IZ-D3'
!           CASE (36)  
!              IMEPR='K-IZ-D4'
!           CASE (37)  
!              IMEPR='K-IZ-D5'
!           CASE (38)  
!              IMEPR='K-IZ-D6'
!           CASE (39)  
!              IMEPR='K-IZ-D7'
!           CASE (40)  
!              IMEPR='K-IZ-D8'
!           CASE (41)  
!              IMEPR='K-IZ-D9'
!           CASE (42)  
!              IMEPR='K-IZ-D10'
!           CASE (43) 
!              IMEPR='K-IZ-D11'
!           CASE (44)  
!              IMEPR='K-IZ-D12'
!           CASE (45)  
!              IMEPR='K-S-1'
!           CASE (46)  
!              IMEPR='K-S-2'
!           CASE (47)  
!              IMEPR='K-S-3'
!           CASE (48) 
!              IMEPR='K-S-4'
!           CASE (49) 
!              IMEPR='K-D-1'
!           CASE (50) 
!              IMEPR='K-D-2'
!           CASE (51)  
!              IMEPR='K-Z'
         CASE DEFAULT
         END SELECT
        ENDIF
         WRITE(II,2005) IMEPR
!          WRITE(II,2001) I
C od 3 linije  se razlikuje properti u zavisnosti od tipa konacnog elementa
              SELECT CASE (IPRO)
! C      grede bez i sa medjucvorovima
!                  CASE (5, 37)
! 			    WRITE(II,5600) NULA,ISHAPEBEAM,IJEDAN,NULA
! 			    WRITE(II,5600) IDESET
! 			    WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
! 			    WRITE(II,5600) NULA,NULA
! 			    WRITE(II,5600) ISPET
! 			    WRITE(II,5500) APRESEK,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) APRESEK,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5600) NULA
! 			    WRITE(II,5600) NULA
C 
C 2D elementi
			CASE (19, 20)
!  v4.4			
                   IF (IVER.EQ.44) THEN
		          WRITE(II,5600) NULA,NULA,NULA,NULA
                          WRITE(II,5600) ICPET
                WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
                WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
                WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
                WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
                WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
                WRITE(II,5600) NULA,NULA,NULA,NULA,NULA
			    WRITE(II,5600) IDTRI
			    WRITE(II,5500) APRESEK,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO
                      ELSEIF(IVER.EQ.10)THEN
C v10
			    WRITE(II,5600) NULA,NULA,NULA,NULA
			    WRITE(II,5600) IDESET
			    WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
			    WRITE(II,5600) NULA,NULA
			    WRITE(II,5600) ISPET
			    WRITE(II,5500) APRESEK,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
			    WRITE(II,5600) NULA
			    WRITE(II,5600) NULA
                      ENDIF
C 3D elementi			
! 			CASE (25, 26)
! 			    WRITE(II,5600) NULA,NULA,NULA,NULA
! 			    WRITE(II,5600) IDESET
! 			    WRITE(II,5600) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
! 			    WRITE(II,5600) NULA,NULA
! 			    WRITE(II,5600) ISPET
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,Z2002ERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5500) ZERO,ZERO,ZERO,ZERO,ZERO
! 			    WRITE(II,5600) NULA
! 			    WRITE(II,5600) NULA
			CASE DEFAULT
                     END SELECT
        ENDDO
C        WRITE(II,5100) IND
!       enddo
        WRITE(II,5100) IND
      RETURN
C
 5100 FORMAT(I5)
 5500 FORMAT(10(1PE13.5))
 5600 FORMAT(10I5)
 5700 FORMAT(I4,1(1PE13.5),3I4)
 2000 FORMAT(' MAT',I4)
 2001 FORMAT(' PRO',I4)
 2002 FORMAT(A80)
 2005 FORMAT(A30)
 5200 FORMAT(F4.1,',')
      END
C======================================================================
C
C======================================================================
      SUBROUTINE VTGRAUK(NP,II,ISNUMER)
      USE NODES
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO PRINT COORDINATES IN UNIVERSAL FILE
CS.    P R O G R A M
CS.       ZA STAMPANJE KOORDINATA CVOROVA U UNIVERZALNI FILE
C .
C ......................................................................
C
!       DIMENSION CORD(3,*)
!       DIMENSION ID(1,*)
      DIMENSION IDD(6)
C
      IDD(2)=0
      IDD(3)=0
      IDD(4)=0
      IDD(5)=0
      IDD(6)=0
C
      NULA=0
C
      IND=-1
      ITYP=403
      ICOL=46
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
      DO 10 I=1,NP
         IF(ISNUMER.EQ.0) THEN
           NI=I
         ELSE
           NI=NCVEL(I)
         ENDIF
         DO 20 IJ=1,1
            IF(ID(IJ,I).GT.0) THEN
               IDD(IJ)=0
            ELSE
               IDD(IJ)=1
            ENDIF
   20    CONTINUE
         WRITE(II,5000) NI,NULA,NULA,NULA,ICOL,(IDD(J),J=1,6),
     +                  (CORD(J,I),J=1,3),NULA
   10 CONTINUE
      WRITE(II,5100) IND
      RETURN
C
 5000 FORMAT(I10,10I3,3(1PE13.5),I2)
 5100 FORMAT(I5)
      END
C======================================================================
      SUBROUTINE VTGRAUKP(II,IPR)
      USE PRESEK
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO PRINT COORDINATES CROSS-SECTION IN UNIVERSAL FILE
CS.    P R O G R A M
CS.       ZA STAMPANJE KOORDINATA CVOROVA PRESEKA U UNIVERZALNI FILE
C .
C ......................................................................
C
!       DIMENSION CORD(3,*)
!       DIMENSION ID(1,*)
      DIMENSION IDD(6)
C
      IDD(1)=0
      IDD(2)=0
      IDD(3)=0
      IDD(4)=0
      IDD(5)=0
      IDD(6)=0
C
      NULA=0
C
      IND=-1
      ITYP=403
      ICOL=46
      WRITE(II,5100) IND
      WRITE(II,5100) ITYP
      DO 10 I=1,NPRESEK(IPR)
!          IF(ISNUMER.EQ.0) THEN
           NI=NP_ID(I,IPR)
!          ELSE
!            NI=NCVEL(I)
!          ENDIF
!          DO 20 IJ=1,1
!             IF(ID(IJ,I).GT.0) THEN
!                IDD(IJ)=0
!             ELSE
!                IDD(IJ)=1
!             ENDIF
!    20    CONTINUE
         WRITE(II,5000) NI,NULA,NULA,NULA,ICOL,(IDD(J),J=1,6),
     +                  (NP_COORDS(J,I,IPR),J=1,3),NULA
   10 CONTINUE
      WRITE(II,5100) IND
      RETURN
C
 5000 FORMAT(I10,10I3,3(1PE13.5),I2)
 5100 FORMAT(I5)
      END
C=====================================================================
C=====================================================================
      SUBROUTINE VTGRAU2(NEL,NCVE,NE,IGRAF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C .......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT 2/D ELEMENTS DATA IN NEUTRAL GRAPHICS FILE
CS.   P R O G R A M
CS.       ZA STAMPANJE 2/D ELEMENATA U NEUTRALNI FILE
C .
C .......................................................................
C 
      COMMON /MATER/ NUMMAT
      COMMON /DODAT/ NDIMM
      DIMENSION NEL(NDIMM,*)
      DIMENSION FIZ(14)         
C      COMMON /CDEBUG/ IDEBUG
C
C      IF(IDEBUG.GT.0) PRINT *, ' TGRAU2'
C
C     FIZICKE OSOBINE
C
      NULA=0
      ZERO=0.
      ONE=1.
      JEDAN=1
      INDPR=25
      INDPD=2
      I11=11
      I2=2
      I4=4
      I8=8
      I14=14
      I48=48
C
C     E L E M E N T I
C
      NNCVE=NCVE
      IF(NCVE.LT.8) NCVE=4
      IF(NCVE.EQ.9) NCVE=8
C     GRAFICKI OPIS ELEMENTA:
C     MEMBRANE - SA 4 CVOROVA = 13, SA 8 CVOROVA = 14
C     PLANE STRAIN - SA 4 CVOROVA = 19, SA 8 CVOROVA = 20
C     AXISYMMETRIC - SA 4 CVOROVA = 23, SA 8 CVOROVA = 24
C     PLATE - SA 4 CVOROVA = 17, SA 8 CVOROVA = 18
C     LAMINATE - SA 4 CVOROVA = 21, SA 8 CVOROVA = 22
      IFGD=19
      IF(NCVE.EQ.8) IFGD=20
C     VRSTA 2D ELEMENTA: SA 4 CVOROVA = 4, SA 8 CVOROVA = 5
      IFDI=4
      IF(NCVE.EQ.8) IFDI=5
C     TABELA FIZICKIH OSOBINA
C      IPTN=ISUMGR
C     BROJ CVOROVA NA ELEMENTU
      NNODS=NCVE
      IND=-1
      ITYP=404
      WRITE(IGRAF,1100) IND
      WRITE(IGRAF,1100) ITYP
      DO 10 I=1,NE
C        TABELA MATERIJALA
         MPTN=NEL(NNCVE+1,I)
         NEMA=NEL(NNCVE+2,I)
         IF(NEMA.GT.0) MPTN=NUMMAT+1
C        BOJA  
         ICOL=MPTN
c         ICOL=124-(MPTN-1)*2
c         ICOL=120-(MPTN-1)*3
         WRITE(IGRAF,1000) I,ICOL,MPTN,IFGD,IFDI,MPTN,NULA,NULA
         IF(NCVE.EQ.4) THEN
            WRITE(IGRAF,1001) (NEL(J,I),J=1,4),NULA,NULA,NULA,NULA,
     +                        NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                        NULA,NULA
         ELSE
            WRITE(IGRAF,1001) (NEL(J,I),J=1,8),
     +                        NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                        NULA,NULA
         ENDIF
         WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         WRITE(IGRAF,1003) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                     NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
   10 CONTINUE
      WRITE(IGRAF,1100) IND
C      ISUMEL=ISUMEL+NE
      NCVE=NNCVE
      RETURN
C
 1000 FORMAT(I10,7I4)
 1001 FORMAT(10I10)
 1011 FORMAT(10I2)
 1002 FORMAT(3F3.0)
 1003 FORMAT(16I2)
 1100 FORMAT(I5)
 1200 FORMAT(6(1PE13.5))
      END
C=====================================================================
      SUBROUTINE VTGRAU3(NCVE,NE,IGRAF,ISNUMER)
      USE NODES
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C .......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT (3/D) ELEMENTS DATA IN NEUTRAL GRAPHICS FILE
CS.   P R O G R A M
CS.       ZA STAMPANJE (3/D) ELEMENATA U NEUTRALNI FILE
C .
C .......................................................................
C 
      COMMON /MATER/ NUMMAT
      COMMON /DODAT/ NDIMM
!       DIMENSION NEL(NDIMM,*)
      INTEGER IELNODE
      DIMENSION FIZ(14),IELNODE(20)     
C      COMMON /CDEBUG/ IDEBUG
C
C      IF(IDEBUG.GT.0) PRINT *, ' TGRAU3'
C
C     FIZICKE OSOBINE
C
      NULA=0
      ZERO=0.
      ONE=1.
      JEDAN=1
      INDPR=25
      INDPD=2
      I11=11
      I2=2
      I4=4
      I8=8
      I14=14
      I48=48
C
C     E L E M E N T I
C
      NNCVE=NCVE
      IND=-1
      ITYP=404
      WRITE(IGRAF,1100) IND
      WRITE(IGRAF,1100) ITYP
      DO 10 I=1,NE
!         NTYPE=elemtip(NE)/10
!         NDIM=elemtip(NE)-10*NTYPE
C     GRAFICKI OPIS ELEMENTA: SA 8 CVOROVA = 25, SA 20 CVOROVA = 26
C     VRSTA 3/D ELEMENTA: SA 8 CVOROVA = 8, SA 20 CVOROVA = 12
      SELECT CASE (elemtip(I))
      CASE (12)
        IFGD=5
        IFDI=0
      CASE (13)
        IFGD=37
        IFDI=1
      CASE (23)
        IFGD=19
        IFDI=2
      CASE (24)
        IFGD=19
        IFDI=4
      CASE (26)
        IFGD=20
        IFDI=3
      CASE (28)
        IFGD=20
        IFDI=5
      CASE (34)
        IFGD=25
        IFDI=6
      CASE (310)
        IFGD=26
        IFDI=10
      CASE (38)
        IFGD=25
        IFDI=8
      CASE (320)
        IFGD=26
        IFDI=12
      CASE DEFAULT
        WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(I)
      END SELECT
C     TABELA FIZICKIH OSOBINA
C      IPTN=ISUMGR
C     BROJ CVOROVA NA ELEMENTU
      
C        TABELA MATERIJALA
         MPTN=NEL(NNCVE+1,I)
         NEMA=NEL(NNCVE+2,I)
         IF(NEMA.GT.0) MPTN=NUMMAT+1
C        BOJA  
         ICOL=MPTN
         IPROP=elemtip(I)*100+MPTN
c         ICOL=124-(MPTN-1)*2
c         ICOL=120-(MPTN-1)*3
C         WRITE(IGRAF,1000) I,ICOL,elemtip(I)+MPTN,IFGD,IFDI,MPTN,NULA,NULA
         IF(ISNUMER.EQ.0)THEN
          WRITE(IGRAF,1111) I,ICOL,IPROP,IFGD,IFDI,MPTN,NULA,NULA
         ELSEIF(ISNUMER.EQ.1)THEN
         WRITE(IGRAF,1111) MCVEL(I),ICOL,IPROP,IFGD,IFDI,MPTN,NULA,NULA
         ENDIF
!  slobodna numeracija         
         IF(ISNUMER.EQ.0)THEN
           do J=1,20
             IELNODE(J)=NEL(J,I)
           ENDDO
         ELSEIF(ISNUMER.EQ.1)THEN
           do J=1,20
             NN=NEL(J,I)
             IELNODE(J)=NN
             IF(NN.NE.0) IELNODE(J)=NCVEL(NN)
           ENDDO
         ENDIF
         SELECT CASE (elemtip(I))
         CASE (12)
!             WRITE(IGRAF,1001) (NEL(J,I),J=1,2),NULA,NULA,NULA,NULA,NULA,
            WRITE(IGRAF,1001) (IELNODE(J),J=1,2),NULA,NULA,NULA,NULA,NULA,
     +                     NULA,NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                     NULA,NULA
            WRITE(IGRAF,1002) ONE,ZERO,ZERO
         CASE (13)
!             WRITE(IGRAF,1001) (NEL(J,I),J=1,3),NULA,NULA,NULA,NULA,NULA,
          WRITE(IGRAF,1001) (IELNODE(J),J=1,3),NULA,NULA,NULA,NULA,NULA,
     +                     NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                     NULA,NULA
            WRITE(IGRAF,1002) ONE,ZERO,ZERO
         CASE (23)
!             WRITE(IGRAF,1001) (NEL(J,I),J=1,3),NULA,NULA,NULA,NULA,
           WRITE(IGRAF,1001) (IELNODE(J),J=1,3),NULA,NULA,NULA,NULA,
     +                        NULA,NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                        NULA,NULA
            WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         CASE (24)
!             WRITE(IGRAF,1001) (NEL(J,I),J=1,4),NULA,NULA,NULA,NULA,
            WRITE(IGRAF,1001) (IELNODE(J),J=1,4),NULA,NULA,NULA,NULA,
     +                        NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                        NULA,NULA
            WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         CASE (26)
!             WRITE(IGRAF,1001) (NEL(J,I),J=1,3),NULA,(NEL(J,I),J=4,6),NULA,NULA
        WRITE(IGRAF,1001) (IELNODE(J),J=1,3),NULA,(IELNODE(J),J=4,6)
     +        ,NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA
     +                        ,NULA,NULA,NULA,NULA
            WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         CASE (28)
!             WRITE(IGRAF,1001) (NEL(J,I),J=1,8),
            WRITE(IGRAF,1001) (IELNODE(J),J=1,8),
     +                        NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                        NULA,NULA
            WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         CASE (34)
!             WRITE(IGRAF,1001) NEL(1,I),NEL(2,I),NEL(3,I),NULA,NEL(4,I),
            WRITE(IGRAF,1001) (IELNODE(J),J=1,3),NULA,IELNODE(4),
     +                        NULA,NULA,NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                        NULA,NULA
            WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         CASE (310)
!              WRITE(IGRAF,1001) NEL(1,I),NEL(2,I),NEL(3,I),NULA,NEL(4,I),
!      +                        NULA,NULA,NULA,NEL(5,I),NEL(6,I)
!              WRITE(IGRAF,1001) NEL(7,I),NULA,NEL(8,I),NEL(9,I),NEL(10,I)
             WRITE(IGRAF,1001) IELNODE(1),IELNODE(2),IELNODE(3),NULA,
     +                        IELNODE(4),NULA,NULA,NULA,IELNODE(5)
     +                        ,IELNODE(6)
            WRITE(IGRAF,1001) IELNODE(7),NULA,IELNODE(8),IELNODE(9),
     +                        IELNODE(10),NULA,NULA,NULA,NULA,NULA
!      +                        ,NULA,NULA,NULA,NULA,NULA
            WRITE(IGRAF,1002) ZERO,ZERO,ZERO
          CASE (38)
!             WRITE(IGRAF,1001) (NEL(J,I),J=5,8),(NEL(J,I),J=1,4),
            WRITE(IGRAF,1001) (IELNODE(J),J=5,8),(IELNODE(J),J=1,4),
     +                        NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                        NULA,NULA
            WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         CASE (320)
!             WRITE(IGRAF,1001) (NEL(J,I),J=5,8),(NEL(J,I),J=1,4),
!      +                        NEL(17,I),NEL(18,I)
!             WRITE(IGRAF,1001) NEL(19,I),NEL(20,I),(NEL(J,I),J=9,12),
!      +                        (NEL(J,I),J=13,16)
            WRITE(IGRAF,1001) (IELNODE(J),J=5,8),(IELNODE(J),J=1,4),
     +                        IELNODE(17),IELNODE(18)
            WRITE(IGRAF,1001) IELNODE(19),IELNODE(20),
     +                        (IELNODE(J),J=9,12),
     +                        (IELNODE(J),J=13,16)
            WRITE(IGRAF,1002) ZERO,ZERO,ZERO
          CASE DEFAULT
          WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(I)
         END SELECT

         WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         WRITE(IGRAF,1003) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                     NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
   10 CONTINUE
      WRITE(IGRAF,1100) IND
C      ISUMEL=ISUMEL+NE
      NCVE=NNCVE
      RETURN
C
 1000 FORMAT(I10,7I4)
 1111 FORMAT(I10,7I6)
 1001 FORMAT(10I10)
 1011 FORMAT(10I2)
 1002 FORMAT(3F3.0)
 1003 FORMAT(16I2)
 1100 FORMAT(I5)
 1200 FORMAT(6(1PE13.5))
      END
C=======================================================================
C=====================================================================
      SUBROUTINE VTGRAU3P(IGRAF,IPR)
      USE PRESEK
!       USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C .......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT (2/D) ELEMENTS DATA IN NEUTRAL GRAPHICS FILE
CS.   P R O G R A M
CS.       ZA STAMPANJE (2/D) ELEMENATA preseka U NEUTRALNI FILE
C .
C .......................................................................
C 
      COMMON /MATER/ NUMMAT
      COMMON /DJERDAP/ IDJERDAP,ISPRESEK
!       COMMON /DODAT/ NDIMM
!       DIMENSION NEL(NDIMM,*)
      INTEGER IELNODE
      DIMENSION IELNODE(20)     
C      COMMON /CDEBUG/ IDEBUG
C
C      IF(IDEBUG.GT.0) PRINT *, ' TGRAU3P'
C
C     FIZICKE OSOBINE
C
      NULA=0
      ZERO=0.
      ONE=1.
      JEDAN=1
      INDPR=25
      INDPD=2
      I11=11
      I2=2
      I4=4
      I8=8
      I14=14
      I48=48
      NCVE=3
C
C     E L E M E N T I
C
      NNCVE=NCVE
      IND=-1
      ITYP=404
      WRITE(IGRAF,1100) IND
      WRITE(IGRAF,1100) ITYP
      DO 10 I=1,NPELEM(IPR)
!         NTYPE=elemtip(NE)/10
!         NDIM=elemtip(NE)-10*NTYPE
C     GRAFICKI OPIS ELEMENTA: SA 8 CVOROVA = 25, SA 20 CVOROVA = 26
C     VRSTA 3/D ELEMENTA: SA 8 CVOROVA = 8, SA 20 CVOROVA = 12
!       SELECT CASE (elemtip(I))
!       CASE (12)
!         IFGD=5
!         IFDI=0
!       CASE (13)
!         IFGD=37
!         IFDI=1
!       CASE (23)
         IFGD=19
         IFDI=2
!       CASE (24)
!         IFGD=19
!         IFDI=4
!       CASE (26)
!         IFGD=20
!         IFDI=3
!       CASE (28)
!         IFGD=20
!         IFDI=5
!       CASE (34)
!         IFGD=25
!         IFDI=6
!       CASE (310)
!         IFGD=26
!         IFDI=10
!       CASE (38)
!         IFGD=25
!         IFDI=8
!       CASE (320)
!         IFGD=26
!         IFDI=12
!       CASE DEFAULT
!         WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(I)
!       END SELECT
C     TABELA FIZICKIH OSOBINA
C      IPTN=ISUMGR
C     BROJ CVOROVA NA ELEMENTU
      
C        TABELA MATERIJALA
!          MPTN=NEL(NNCVE+1,I)
!          NEMA=NEL(NNCVE+2,I)
!          IF(NEMA.GT.0) MPTN=NUMMAT+1
!          elemtip(I)=23
!          MPTN=1
!  stampanje propertija na osnovu zadnjeg cvora elemenata
!          write(*,*) "ISPRESEK",ISPRESEK
         MPTN=NELP(5,I,IPR)
         IF(ISPRESEK.GE.-1) THEN
          IF(NELP(5,I,IPR).EQ.1) THEN
           MPTN=1
          elseif (NELP(5,I,IPR).EQ.9) THEN
           MPTN=9
          else
            DO J=2,4
               ncvorpr=NELP(J,I,IPR)
               MPTN=NPROP(ncvorpr,IPR)
               IF(MPTN.NE.1.AND.MPTN.NE.9) GO TO 20
            ENDDO
          ENDIF
         ELSEIF(ISPRESEK.EQ.-2) THEN
            ncvorpr=NELP(2,I,IPR)
            MPTN1=NPROP(ncvorpr,IPR)
            ncvorpr=NELP(3,I,IPR)
            MPTN2=NPROP(ncvorpr,IPR)
            ncvorpr=NELP(4,I,IPR)
            MPTN3=NPROP(ncvorpr,IPR)
            IF(MPTN.NE.7) GO TO 20
            IF((MPTN1.EQ.MPTN2).OR.(MPTN1.EQ.MPTN3)) MPTN=MPTN1
            IF(MPTN2.EQ.MPTN3) MPTN=MPTN2
            IF(MPTN.EQ.48) THEN
              MPTN=8
            ELSE
              MPTN=7
            ENDIF
         ENDIF
!          write(*,*) 'MPTN',MPTN,MPTN1,MPTN2
C        BOJA  
   20    ICOL=124-(MPTN-1)*2
!          IPROP=1
         elemtipe=23
         IPROP=MPTN
!           write(*,*) 'I,MPTN,IPROP',I,MPTN,IPROP
!          IPROP=elemtipe*100+MPTN
c         ICOL=124-(MPTN-1)*2
c         ICOL=120-(MPTN-1)*3
C         WRITE(IGRAF,1000) I,ICOL,elemtip(I)+MPTN,IFGD,IFDI,MPTN,NULA,NULA
!          IF(ISNUMER.EQ.0)THEN
          WRITE(IGRAF,1111) I,ICOL,IPROP,IFGD,IFDI,MPTN,NULA,NULA
!          ELSEIF(ISNUMER.EQ.1)THEN
!          WRITE(IGRAF,1111) MCVEL(I),ICOL,IPROP,IFGD,IFDI,MPTN,NULA,NULA
!          ENDIF
!  slobodna numeracija         
!          IF(ISNUMER.EQ.0)THEN
           do J=1,3
             IELNODE(J)=NELP(J+1,I,IPR)
           ENDDO
!          ELSEIF(ISNUMER.EQ.1)THEN
!            do J=1,20
!              NN=NEL(J,I)
!              IELNODE(J)=NN
!              IF(NN.NE.0) IELNODE(J)=NCVEL(NN)
!            ENDDO
!          ENDIF
!          SELECT CASE (elemtip(I))
!          CASE (12)
! !             WRITE(IGRAF,1001) (NEL(J,I),J=1,2),NULA,NULA,NULA,NULA,NULA,
!             WRITE(IGRAF,1001) (IELNODE(J),J=1,2),NULA,NULA,NULA,NULA,NULA,
!      +                     NULA,NULA,NULA
!             WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
!      +                     NULA,NULA
!             WRITE(IGRAF,1002) ONE,ZERO,ZERO
!          CASE (13)
! !             WRITE(IGRAF,1001) (NEL(J,I),J=1,3),NULA,NULA,NULA,NULA,NULA,
!           WRITE(IGRAF,1001) (IELNODE(J),J=1,3),NULA,NULA,NULA,NULA,NULA,
!      +                     NULA,NULA
!             WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
!      +                     NULA,NULA
!             WRITE(IGRAF,1002) ONE,ZERO,ZERO
!          CASE (23)
!             WRITE(IGRAF,1001) (NEL(J,I),J=1,3),NULA,NULA,NULA,NULA,
           WRITE(IGRAF,1001) (IELNODE(J),J=1,3),NULA,NULA,NULA,NULA,
     +                        NULA,NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                        NULA,NULA
            WRITE(IGRAF,1002) ZERO,ZERO,ZERO
!          CASE (24)
! !             WRITE(IGRAF,1001) (NEL(J,I),J=1,4),NULA,NULA,NULA,NULA,
!             WRITE(IGRAF,1001) (IELNODE(J),J=1,4),NULA,NULA,NULA,NULA,
!      +                        NULA,NULA
!             WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
!      +                        NULA,NULA
!             WRITE(IGRAF,1002) ZERO,ZERO,ZERO
!          CASE (26)
! !             WRITE(IGRAF,1001) (NEL(J,I),J=1,3),NULA,(NEL(J,I),J=4,6),NULA,NULA
!         WRITE(IGRAF,1001) (IELNODE(J),J=1,3),NULA,(IELNODE(J),J=4,6)
!      +        ,NULA,NULA
!             WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA
!      +                        ,NULA,NULA,NULA,NULA
!             WRITE(IGRAF,1002) ZERO,ZERO,ZERO
!          CASE (28)
! !             WRITE(IGRAF,1001) (NEL(J,I),J=1,8),
!             WRITE(IGRAF,1001) (IELNODE(J),J=1,8),
!      +                        NULA,NULA
!             WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
!      +                        NULA,NULA
!             WRITE(IGRAF,1002) ZERO,ZERO,ZERO
!          CASE (34)
! !             WRITE(IGRAF,1001) NEL(1,I),NEL(2,I),NEL(3,I),NULA,NEL(4,I),
!             WRITE(IGRAF,1001) (IELNODE(J),J=1,4),
!      +                        NULA,NULA,NULA,NULA,NULA
!             WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
!      +                        NULA,NULA
!             WRITE(IGRAF,1002) ZERO,ZERO,ZERO
!          CASE (310)
!              WRITE(IGRAF,1001) IELNODE(1),IELNODE(2),IELNODE(3),NULA,
!      +                        IELNODE(4),NULA,NULA,NULA,IELNODE(5)
!      +                        ,IELNODE(6)
!             WRITE(IGRAF,1001) IELNODE(7),NULA,IELNODE(8),IELNODE(9),
!      +                        IELNODE(10),NULA,NULA,NULA,NULA,NULA
! !      +                        ,NULA,NULA,NULA,NULA,NULA
!             WRITE(IGRAF,1002) ZERO,ZERO,ZERO
!           CASE (38)
! !             WRITE(IGRAF,1001) (NEL(J,I),J=5,8),(NEL(J,I),J=1,4),
!             WRITE(IGRAF,1001) (IELNODE(J),J=5,8),(IELNODE(J),J=1,4),
!      +                        NULA,NULA
!             WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
!      +                        NULA,NULA
!             WRITE(IGRAF,1002) ZERO,ZERO,ZERO
!          CASE (320)
!             WRITE(IGRAF,1001) (IELNODE(J),J=5,8),(IELNODE(J),J=1,4),
!      +                        IELNODE(17),IELNODE(18)
!             WRITE(IGRAF,1001) IELNODE(19),IELNODE(20),
!      +                        (IELNODE(J),J=9,12),
!      +                        (IELNODE(J),J=13,16)
!             WRITE(IGRAF,1002) ZERO,ZERO,ZERO
!           CASE DEFAULT
!           WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(I)
!          END SELECT

         WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         WRITE(IGRAF,1003) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                     NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
   10 CONTINUE
      WRITE(IGRAF,1100) IND
C      ISUMEL=ISUMEL+NE
      NCVE=NNCVE
      RETURN
C
 1000 FORMAT(I10,7I4)
 1111 FORMAT(I10,7I6)
 1001 FORMAT(10I10)
 1011 FORMAT(10I2)
 1002 FORMAT(3F3.0)
 1003 FORMAT(16I2)
 1100 FORMAT(I5)
 1200 FORMAT(6(1PE13.5))
      END
C=======================================================================
      SUBROUTINE VSTAU09(RTH,NPP,II,IND,IKOL,KOR,NZAD,NUMZAD,
     +                  KONT,POT,ISNUMBER)
      USE NODES
      USE ELEMENTS
      USE KONTURE
      USE PREDISCRIBED
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
      COMMON /VOPTIM/ NKONT,MAXLIN
      DIMENSION RTH(IKOL,*),FSP(6),IN(8),NZAD(3,*)
      DIMENSION KONT(9,MAXLIN,*),POT(*)
!       DIMENSION KONT1(40000)
C
      FSP(3)=0.
C
      JEDAN=1
      NULA=0
      ZERO=0.
      ONE=1.
      MJ=-1
      M451=451
      M1=1
      M2=2
      M3=3
      M4=4
      M5=5
      M6=6
      M7=7
      M8=8
      M9=9
C      II=49
C
      IN(1)=IND
      IN(2)=IND+1
      IN(3)=IND+2
      IN(4)=IND+3
      IN(5)=IND+4
      IN(6)=IND+5
      IN(7)=IND+6
      IN(8)=IND+7
C
C PRIVREMENO
!       IF(IN(1).EQ.10) THEN
!       DO I=1,40000
!          KONT1(I)=0
!       ENDDO
!       ENDIF
C      
      IP=0
C
   30 IP=IP+1
      WRITE(II,1000) KOR,IN(IP),JEDAN
      IF(IN(1).EQ.1.AND.IP.EQ.1) WRITE(II,3001) 
      IF(IN(1).EQ.9.AND.IP.EQ.1) WRITE(II,3002) 
      IF(IN(1).EQ.10.AND.IP.EQ.1) WRITE(II,3003) 
      IF(IN(1).EQ.11.AND.IP.EQ.1) WRITE(II,3004) 
      IF(IN(1).EQ.19.AND.IP.EQ.1) WRITE(II,3019) 
      IF(IN(1).EQ.41.AND.IP.EQ.1) WRITE(II,3041) 
      IF(IN(1).EQ.41.AND.IP.EQ.2) WRITE(II,3042) 
      IF(IN(1).EQ.41.AND.IP.EQ.3) WRITE(II,3043) 
      IF(IN(1).EQ.41.AND.IP.EQ.4) WRITE(II,3044) 
      IF(IN(1).EQ.51.AND.IP.EQ.1) WRITE(II,3051) 
      IF(IN(1).EQ.51.AND.IP.EQ.2) WRITE(II,3052) 
      IF(IN(1).EQ.51.AND.IP.EQ.3) WRITE(II,3053) 
      IF(IN(1).EQ.51.AND.IP.EQ.4) WRITE(II,3054) 
      IF(IN(1).EQ.61.AND.IP.EQ.1) WRITE(II,3061) 
      IF(IN(1).EQ.61.AND.IP.EQ.2) WRITE(II,3062) 
      IF(IN(1).EQ.61.AND.IP.EQ.3) WRITE(II,3063) 
      IF(IN(1).EQ.61.AND.IP.EQ.4) WRITE(II,3064) 
      WRITE(II,1010) ZERO,ZERO,ZERO
C OVO SU NULE ZA TEMPERATURE
      IF(IP.EQ.1.AND.IN(1).LT.40)
     +WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.1.AND.IN(1).GT.40)
     +WRITE(II,1000)IN(2),IN(3),IN(4),NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.2)
     +  WRITE(II,1000)IN(2),NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.3)
     +  WRITE(II,1000)NULA,IN(3),NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.4)
     +  WRITE(II,1000)NULA,NULA,IN(4),NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
C
C 0=ANY,1=DISPL,2=ACCEL,3=FORCE,4=STRESS,5=STRAIN,6=TEMP,OTHER=USER
C     NODAL(7) ,ELEMENTAL(8)
C dodat sledeci red za verziju neu fajla 10.0 /Sneza 26.12.2016.
      WRITE(II,1000) NULA
      IF(IN(1).EQ.1) WRITE(II,1000) NULA,NULA,M6,M7
      IF(IN(1).EQ.9) WRITE(II,1000) NULA,NULA,M7,M7
      IF(IN(1).EQ.10) WRITE(II,1000) NULA,NULA,M8,M8
      IF(IN(1).EQ.11) WRITE(II,1000) NULA,NULA,M6,M7
      IF(IN(1).EQ.19) WRITE(II,1000) NULA,NULA,M7,M7
      IF(IN(1).EQ.41.OR.IN(1).EQ.51.OR.IN(1).EQ.61)
     1WRITE(II,1000) NULA,NULA,M3,M7
C
C =1, CAN NOT LINEARY COMBINE THIS OUTPUT
C =1, COMP(0-2) ARE THE X,Y,Z COMPONENT
C =1, THIS VECTOR HAS CONTROIDAL OR NODAL OUTPUT
      IF(IN(1).GT.40) THEN
         IF(IP.EQ.1) THEN
            WRITE(II,1000) JEDAN,JEDAN,JEDAN
         ELSE
            WRITE(II,1000) NULA,JEDAN,JEDAN
         ENDIF
      ELSE
         WRITE(II,1000) JEDAN,NULA,JEDAN
      ENDIF
      DO 10 I=1,NPP
         DO 20 J=1,IKOL
            FSP(J) = 0.0D0
            IF(IKOL.EQ.1) THEN
               IF(IN(1).EQ.1) THEN
                  IF(ID(J,I).EQ.0) GO TO 20
                  K = ID(J,I)
!                   write(*,*) "I,K,NCVEL(I",I,K,NCVEL(I)
                  IF(K.GT.0) THEN
                     IF(INDFS.GT.0) THEN
c                       POTENCIJAL ISPOD SLOBODNE POVRSINE
                        IF((POT(K)-0.001).GT.CORD(IOSA,I)) then
                            FSP(J)=RTH(J,K)
                        else 
                            go to 10
                        endif
                     ELSE   
                        IF((POT(K)-0.001).GT.CORD(IOSA,I)) then
                            FSP(J)=RTH(J,K)
                        else 
                            go to 10
                        endif
c  Sneza 10.9.2015  'stampa vrednost 0 iznad slobodne povrsine
C                       FSP(J)=RTH(J,K) 
                     ENDIF
                  ENDIF
               ELSEIF(IN(1).EQ.11) THEN
                  IF(ID(J,I).EQ.0) GO TO 20
                  K = ID(J,I)
                  IF(K.GT.0) THEN
C    provera da li je cvor okvasen
                IF((IPAKT.EQ.0).AND.((POT(K)-0.0001)
     1             .LT.CORD(IOSA,I))) goto 10
                 FSP(J)=PPC
C                 DUBINA VODE ILI PRITISAK VODENOG STUBA
                  FSP(J)=(RTH(1,K)-CORD(IOSA,I))
!                      IF(FSP(J).GT.0.) FSP(J)=0.     
                ENDIF
               ELSEIF(IN(1).EQ.9) THEN
                 FSP(J)=1.0D0
                 DO L=1,NUMZAD 
                    IF(NZAD(1,L).EQ.I) GO TO 1201
                 ENDDO
                 FSP(J)=0.D0
 1201            CONTINUE
               ELSEIF(IN(1).EQ.19) THEN
!                  FSP(J)=0.0
                 DO L=1,NUMZAD 
                    FSP(J)=POTZADC(L)
                    IF(NZAD(1,L).EQ.I) GO TO 1210
                 ENDDO
                 FSP(J)=0.D0
 1210            CONTINUE
               ELSE
                 DO LL=1,NKONT
                    FSP(J)=LL
                 DO L=1,LIN(LL)
                 DO LI=2,9 
C                 DO LI=2,NNODE+1 
                    IF(KONT(LI,L,LL).EQ.I) THEN
!                        IEL=KONT(1,L,LL)
                      IF(ISNUMBER.EQ.0) THEN
                        IEL=KONT(1,L,LL)
                      ELSE
                        IEL=MCVEL(KONT(1,L,LL))
                      ENDIF
!                        IF(KONT1(IEL).EQ.0) THEN
!                           KONT1(IEL)=1
!                           GO TO 1202
!                        ELSE
!                           GO TO 10
!                        ENDIF 
                    ENDIF
                 ENDDO
                 ENDDO
                 ENDDO
                 FSP(J)=0.D0
 1202            CONTINUE
               ENDIF
            ELSE
!  stampa brzine i gradijente samo za okvasene cvorove 19.5.2018            
               K = ID(1,I)
!                IF(INDFS.GT.0) THEN
                  IF(K.GT.0) THEN
          IF((POT(K)-0.001).GT.CORD(IOSA,I)) FSP(J)=RTH(J,K)
c                  ELSE
c                     FSP(J)=RTH(J,I)
                  ENDIF
!                ELSE   
!                   FSP(J)=RTH(J,K)
!                ENDIF
             ENDIF
   20    CONTINUE
         VAL=FSP(1)
         IF(IP.EQ.1.AND.IN(1).GT.40)
     +      VAL=DSQRT(FSP(1)*FSP(1)+FSP(2)*FSP(2)+FSP(3)*FSP(3))
         IF(IP.EQ.2) VAL=FSP(1)
         IF(IP.EQ.3) VAL=FSP(2)
         IF(IP.EQ.4) VAL=FSP(3)
czn         IF(DABS(VAL).LT.1.E-12) GO TO 10
         IF(IN(1).EQ.10) THEN
           IF(ISNUMBER.EQ.0) THEN
            WRITE(II,5000) IEL,VAL
           ELSE
            WRITE(II,5000) MCVEL(IEL),VAL
           ENDIF
         ELSE
           IF(ISNUMBER.EQ.0) THEN
            WRITE(II,5000) I,VAL
           ELSE
!              write(*,*) "I,NCVEL(I)",I,NCVEL(I)
            WRITE(II,5000) NCVEL(I),VAL
           ENDIF
         ENDIF  
   10 CONTINUE
      WRITE(II,5000) MJ,ZERO
      IF(IP.LT.4.AND.IN(1).GT.40) GO TO 30
      RETURN
 1000 FORMAT(10I5)
 1001 FORMAT(I10,F12.4)
 3001 FORMAT('TOTAL POTENTIAL')
 3002 FORMAT('PRESCRIBED POTENTIAL')
 3003 FORMAT('PRESCRIBED SURFACES')
 3004 FORMAT('DEPTH OF WATER')
 3019 FORMAT('PRESCRIBED POTENTIAL2')
 3041 FORMAT('TOTAL VELOCITY')
 3042 FORMAT('VX VELOCITY')
 3043 FORMAT('VY VELOCITY')
 3044 FORMAT('VZ VELOCITY')
 3051 FORMAT('TOTAL GRADIENT')
 3052 FORMAT('GX GRADIENT')
 3053 FORMAT('GY GRADIENT')
 3054 FORMAT('GZ GRADIENT')
 3061 FORMAT('TOTAL FORCE')
 3062 FORMAT('FX FORCE')
 3063 FORMAT('FY FORCE')
 3064 FORMAT('FZ FORCE')
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(3(1PD12.4))
 5000 FORMAT(I10,3(1PD12.4))
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE VNEUTRA(VREME,NPER,NASLOV,IND,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        WITH LOOP OVER TIME PERIODS AND STEPS FOR NEUTRAL FILE
CS.    P R O G R A M
CS.        SA PETLJOM PO VREMENSKIM PERIODIMA I KORACIMA ZA NEUTRAL 
C .
CE.    V A R I A B L E S
CE.        NPER  - TOTAL NUMBER OF PERIODS WITH CONSTANT TIME STEP,
CE.                SEE CARD /3/
C .
C ......................................................................
C
      CHARACTER*80 NASLOV
      DIMENSION VREME(NPER,950)
      DIMENSION IS(1)
C
      MJ=-1
      IF(IND.EQ.2) GO TO 999
C
      JEDAN=1
      NULA=0
      ZERO=0.
      ONE=1.
      M2=2
      M3=3
      M4=4
      M5=5
      M8=8
      M9=9
C
      WRITE(II,1000) MJ
      M450=450
      WRITE(II,1001) M450
C
      MM=19
C      IF(NGENL.EQ.0.AND.NDIN.EQ.0.AND.ISOPS.EQ.0) MM=1
C
CE    BASIC LOOP OVER TIME PERIODS
CS    OSNOVNA PETLJA PO VREMENSKIM PERIODIMA
C
      KKORAK=1
      VVREME=0.
      DO 100 NNPER=1,NPER
C
CE    BASIC LOOP OVER TIME STEPS
CS    OSNOVNA PETLJA PO VREMENSKIM KORACIMA
C
         KORAK=0
   35    KORAK=KORAK+1
         TIME=VREME(NNPER,KORAK)
         VVREME=VVREME+TIME
C
         KORBR=KKORAK
      WRITE(II,1000) KORBR
      NDT=2
      IF (NPER.EQ.1.AND.DABS(VREME(1,2)).LT.1.D-10) NDT=1
      IF(NDT.EQ.1) THEN
         WRITE(II,1005) KORBR
         WRITE(II,1000) M9,MM
         WRITE(II,1010) ZERO
      ELSE
         WRITE(II,1006) KORBR,VVREME
         WRITE(II,1000) M9,MM
         WRITE(II,1010) VVREME
      ENDIF
      WRITE(II,1000) JEDAN
      WRITE(II,2000) NASLOV
C
      KKORAK=KKORAK+1
      IF (DABS(VREME(NNPER,KORAK+1)).LT.1.D-10) THEN
        GOTO 100
      ELSE
        GOTO 35
      ENDIF
C
  100 CONTINUE
C
      WRITE(II,1000) MJ
C
      WRITE(II,1000) MJ
      M451=451
      WRITE(II,1001) M451
      RETURN
C
  999 WRITE(II,1000) MJ
      CLOSE(II,STATUS='KEEP')
      RETURN
C
 1000 FORMAT(10I5)
 1001 FORMAT(I10,F12.4)
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(3(1PD12.4))
 2000 FORMAT(A80)
      END

C======================================================================
      SUBROUTINE VSTAS09(RTH,NPP,II,IND,IKOL,KOR,NZAD,NUMZAD,
     +                  KONT,ISNUMBER)
      USE NODES
      USE ELEMENTS
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /VOPTIM/ NKONT,MAXLIN
      DIMENSION RTH(IKOL,*),FSP(6),IN(8),NZAD(3,*)
      DIMENSION KONT(9,MAXLIN,*)
!       DIMENSION KONT1(40000)
C
      JEDAN=1
      NULA=0
      ZERO=0.
      ONE=1.
      MJ=-1
      M451=451
      M1=1
      M2=2
      M3=3
      M4=4
      M5=5
      M6=6
      M7=7
      M8=8
      M9=9
C      II=49
C
      IN(1)=IND
      IN(2)=IND+1
      IN(3)=IND+2
      IN(4)=IND+3
      IN(5)=IND+4
      IN(6)=IND+5
      IN(7)=IND+6
      IN(8)=IND+7
C
C PRIVREMENO
!       IF(IN(1).EQ.10) THEN
!       DO I=1,40000
!          KONT1(I)=0
!       ENDDO
!       ENDIF
C      
      IP=0
C
   30 IP=IP+1
      WRITE(II,1000) KOR,IN(IP),JEDAN
      IF(IN(1).EQ.1.AND.IP.EQ.1) WRITE(II,3001) 
      IF(IN(1).EQ.9.AND.IP.EQ.1) WRITE(II,3002) 
      IF(IN(1).EQ.10.AND.IP.EQ.1) WRITE(II,3003) 
      IF(IN(1).EQ.41.AND.IP.EQ.1) WRITE(II,3041) 
      IF(IN(1).EQ.41.AND.IP.EQ.2) WRITE(II,3042) 
      IF(IN(1).EQ.41.AND.IP.EQ.3) WRITE(II,3043) 
      IF(IN(1).EQ.41.AND.IP.EQ.4) WRITE(II,3044) 
      IF(IN(1).EQ.51.AND.IP.EQ.1) WRITE(II,3051) 
      IF(IN(1).EQ.51.AND.IP.EQ.2) WRITE(II,3052) 
      IF(IN(1).EQ.51.AND.IP.EQ.3) WRITE(II,3053) 
      IF(IN(1).EQ.51.AND.IP.EQ.4) WRITE(II,3054) 
      WRITE(II,1010) ZERO,ZERO,ZERO
C OVO SU NULE ZA TEMPERATURE
      IF(IP.EQ.1.AND.IN(1).LT.40)
     +WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.1.AND.IN(1).GT.40)
     +WRITE(II,1000)IN(2),IN(3),IN(4),NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.2)
     +  WRITE(II,1000)IN(2),NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.3)
     +  WRITE(II,1000)NULA,IN(3),NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.4)
     +  WRITE(II,1000)NULA,NULA,IN(4),NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
C
C 0=ANY,1=DISPL,2=ACCEL,3=FORCE,4=STRESS,5=STRAIN,6=TEMP,OTHER=USER
C     NODAL(7) ,ELEMENTAL(8)
C dodat sledeci red za verziju neu fajla 10.0 /Sneza 26.12.2016.
      WRITE(II,1000) NULA
      IF(IN(1).EQ.1) WRITE(II,1000) NULA,NULA,M6,M7
      IF(IN(1).EQ.9) WRITE(II,1000) NULA,NULA,M7,M7
      IF(IN(1).EQ.10) WRITE(II,1000) NULA,NULA,M8,M8
      IF(IN(1).EQ.41.OR.IN(1).EQ.51) WRITE(II,1000) NULA,NULA,M3,M7
C
C =1, CAN NOT LINEARY COMBINE THIS OUTPUT
C =1, COMP(0-2) ARE THE X,Y,Z COMPONENT
C =1, THIS VECTOR HAS CONTROIDAL OR NODAL OUTPUT
      IF(IN(1).GT.40) THEN
         IF(IP.EQ.1) THEN
            WRITE(II,1000) JEDAN,JEDAN,JEDAN
         ELSE
            WRITE(II,1000) NULA,JEDAN,JEDAN
         ENDIF
      ELSE
         WRITE(II,1000) JEDAN,NULA,JEDAN
      ENDIF
C      DO 10 I=1,NPP
         DO 20 J=1,IKOL
            FSP(J) = 0.0D0
            IF(IKOL.EQ.1) THEN
               IF(IN(1).EQ.1) THEN
                  IF(ID(J,I).EQ.0) GO TO 20
                  K = ID(J,I)
!                   IF(K.GT.0) THEN
C OVO JE SAMO ZA 3D SLOBODNA POVRSINA
C                     IF(RTH(1,K).GT.CORD(3,I)) FSP(J)=RTH(1,K)
                     FSP(J)=RTH(1,K)
!                   ENDIF
               ELSEIF(IN(1).EQ.9) THEN
                 FSP(J)=1.0D0
                 DO L=1,NUMZAD 
!  koliko je I u ovom uslovu?                 
                    IF(NZAD(1,L).EQ.I) GO TO 1201
                 ENDDO
                 FSP(J)=0.D0
 1201            CONTINUE
               ELSE
                 DO LL=1,NKONT
                    FSP(J)=LL
                 DO L=1,LIN(LL) 
C                 DO LI=2,5 
C                 DO LI=2,NNODE+1 
C                    IF(KONT(LI,L,LL).EQ.I) THEN
                  IF(ISNUMBER.EQ.0) THEN
                    IEL=KONT(1,L,LL)
                  ELSE
                    IEL=MCVEL(KONT(1,L,LL))
                  ENDIF
                  VAL=FSP(1)
                  WRITE(II,5000) IEL,VAL
C                       IF(KONT1(IEL).EQ.0) THEN
C                          KONT1(IEL)=1
C                          GO TO 1202
C                       ELSE
C                          GO TO 10
C                       ENDIF 
C                    ENDIF
C                 ENDDO
                 ENDDO
                 ENDDO
                 FSP(J)=0.D0
 1202            CONTINUE
               ENDIF
            ELSE
               FSP(J)=RTH(J,I)
            ENDIF
   20    CONTINUE
         VAL=FSP(1)
         IF(IP.EQ.1.AND.IN(1).GT.40)
     +      VAL=DSQRT(FSP(1)*FSP(1)+FSP(2)*FSP(2)+FSP(3)*FSP(3))
         IF(IP.EQ.2) VAL=FSP(1)
         IF(IP.EQ.3) VAL=FSP(2)
         IF(IP.EQ.4) VAL=FSP(3)
C         IF(DABS(VAL).LT.1.E-12) GO TO 10
C         IF(IN(1).EQ.10) THEN
C            WRITE(II,5000) IEL,VAL
C         ELSE
C            WRITE(II,5000) I,VAL
C         ENDIF  
C   10 CONTINUE
      WRITE(II,5000) MJ,ZERO
      IF(IP.LT.4.AND.IN(1).GT.40) GO TO 30
      RETURN
 1000 FORMAT(10I5)
 1001 FORMAT(I10,F12.4)
 3001 FORMAT('TOTAL POTENTIAL')
 3002 FORMAT('PRESCRIBED POTENTIAL')
 3003 FORMAT('PRESCRIBED SURFACES')
 3041 FORMAT('TOTAL VELOCITY')
 3042 FORMAT('VX VELOCITY')
 3043 FORMAT('VY VELOCITY')
 3044 FORMAT('VZ VELOCITY')
 3051 FORMAT('TOTAL GRADIENT')
 3052 FORMAT('GX GRADIENT')
 3053 FORMAT('GY GRADIENT')
 3054 FORMAT('GZ GRADIENT')
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(3(1PD12.4))
 5000 FORMAT(I10,3(1PD12.4))
      END
C=====================================================================
      SUBROUTINE VTGRAU3K(KONT,NDIM,IGRAF)
      USE ELEMENTS
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C .......................................................................
C .
CE.   P R O G R A M
CE.       TO PRINTOUT 2/D ELEMENTS DATA IN NEUTRAL GRAPHICS FILE
CS.   P R O G R A M
CS.       ZA STAMPANJE 2/D ELEMENATA U NEUTRALNI FILE
C .
C .......................................................................
C 
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /MATER/ NUMMAT
      COMMON /DODAT/ NDIMM
      DIMENSION KONT(9,MAXLIN,*)
      DIMENSION FIZ(14)         
C      COMMON /CDEBUG/ IDEBUG
C
C      IF(IDEBUG.GT.0) PRINT *, ' TGRAU3K'
C
C     FIZICKE OSOBINE
C
      NULA=0
      ZERO=0.
      ONE=1.
      JEDAN=1
      INDPR=25
      INDPD=2
      I11=11
      I2=2
      I4=4
      I8=8
      I14=14
      I48=48
C
C     E L E M E N T I
C
      NCVE=4
C     GRAFICKI OPIS ELEMENTA:
C     MEMBRANE - SA 4 CVOROVA = 13, SA 8 CVOROVA = 14
C     PLANE STRAIN - SA 4 CVOROVA = 19, SA 8 CVOROVA = 20
C     AXISYMMETRIC - SA 4 CVOROVA = 23, SA 8 CVOROVA = 24
C     PLATE - SA 4 CVOROVA = 17, SA 8 CVOROVA = 18
C     LAMINATE - SA 4 CVOROVA = 21, SA 8 CVOROVA = 22
      IFGD=19
      IF(NCVE.EQ.8) IFGD=20
C     VRSTA 2D ELEMENTA: SA 4 CVOROVA = 4, SA 8 CVOROVA = 5
      IFDI=4
      IF(NCVE.EQ.8) IFDI=5
C     TABELA FIZICKIH OSOBINA
C      IPTN=ISUMGR
C     BROJ CVOROVA NA ELEMENTU
      NNODS=NCVE
      IND=-1
      ITYP=404
      WRITE(IGRAF,1100) IND
      WRITE(IGRAF,1100) ITYP
      NLM=0
      DO 10 KK=1,NKONT
      DO 10 II=1,LIN(KK)
         NLM=NLM+1
C	
           NBREL=KONT(1,II,KK)
C          TABELA MATERIJALA****
           MPTN=NEL(NDIM+1,NBREL)
           NEMA=NEL(NDIM+2,NBREL)
           IF(NEMA.GT.0) MPTN=NUMMAT+1
C odredjivanje topologije elementa
C     VRSTA 2D ELEMENTA: SA 4 CVOROVA = 4, SA 8 CVOROVA = 5 SA 3cvora=2
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
            IFDI=0
            goto 10
          CASE(13)
            IFDI=0
            goto 10
          CASE(23)
            NSCV=2
            goto 10
          CASE(26)
            NSCV=3
            goto 10
          CASE(24)
            IFDI=4
            goto 10
          CASE(28)
            IFDI=4
            goto 10
          CASE(34)
            IFDI=3
          CASE(310)
            IFDI=3
          CASE(38)
            IFDI=4
          CASE DEFAULT
            WRITE(*,*) 'GRESKA, NEMA PRIKAZA ZA DATI ELEMENT'
          END SELECT

C        BOJA  
         ICOL=MPTN
C properti ID
c         ICOL=124-(MPTN-1)*2
c         ICOL=120-(MPTN-1)*3
         WRITE(IGRAF,1000) NLM,ICOL,MPTN,IFGD,IFDI,MPTN,NULA,NULA
C         IF(NCVE.EQ.4) THEN
            WRITE(IGRAF,1001) (KONT(J,II,KK),J=2,5),NULA,NULA,NULA,NULA,
     +                        NULA,NULA
            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                        NULA,NULA
C         ELSE
C            WRITE(IGRAF,1001) (NEL(J,I),J=1,8),
C     +                        NULA,NULA
C            WRITE(IGRAF,1011) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
C     +                        NULA,NULA
C         ENDIF
         WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         WRITE(IGRAF,1002) ZERO,ZERO,ZERO
         WRITE(IGRAF,1003) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,
     +                     NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
   10 CONTINUE
      WRITE(IGRAF,1100) IND
C      ISUMEL=ISUMEL+NE
C      NCVE=NNCVE
      RETURN
C
 1000 FORMAT(I10,7I4)
 1001 FORMAT(10I10)
 1011 FORMAT(10I2)
 1002 FORMAT(3F3.0)
 1003 FORMAT(16I2)
 1100 FORMAT(I6)
 1200 FORMAT(6(1PE13.5))
      END
C=======================================================================
      SUBROUTINE VTGRMATK(KONT,AKONST,NDIM)
      USE ELEMENTS
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /DODAT/ NDIMM
      DIMENSION KONT(9,MAXLIN,*)
      DIMENSION AKONST(3,5,*)
C*    BLOK ZA STAMPANJE MATERIJALA PO ELEMENTIMA ZA GRAFIKU - PRIVREMENO
      JEDAN=1
      NULA=0
      ZERO=0.
      MJ=-1
      M8=8
C
      IN=14
      IGRAF=48
      WRITE(IGRAF,1000) JEDAN,IN,JEDAN
      WRITE(IGRAF,3003) 
      WRITE(IGRAF,1010) ZERO,ZERO,ZERO
      WRITE(IGRAF,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(IGRAF,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
C od verzije femapa 10.0 za vrednosti u elementu postoji jos jedan red
      WRITE(IGRAF,1000) NULA
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
 3003 FORMAT(' KOEFICIJENT FILTRACIJE')
C*

      NLM=0
      DO 723 KK=1,NKONT
C
	DO 134 II=1,LIN(KK)
	NLM=NLM+1
C	
      NBREL=KONT(1,II,KK)
      NMAT=NEL(NDIM+1,NBREL)
      NEMA=NEL(NDIM+2,NBREL)

      AKX=AKONST(3,1,NMAT)
C*    STAMPANJE MATERIJALA PO ELEMENTU
      IF(ISNUMER.EQ.0) THEN
        WRITE(IGRAF,5000) NLM,AKX
      ELSEIF(ISNUMER.EQ.1) THEN
        WRITE(IGRAF,5000) MCVEL(NLM),AKX
      ENDIF
C*
  134 CONTINUE
C
  723 CONTINUE
C
C*    KRAJ BLOKA ZA STAMPANJE PROTOKA PO ELEMENTU
      WRITE(IGRAF,5000) MJ,ZERO
 5000 FORMAT(I10,3(1PD12.4))
C*
      RETURN
      END
C=======================================================================
      SUBROUTINE VTGRMATS(AKONST,NDIM,NE,ISNUMER)
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /DODAT/ NDIMM
!       DIMENSION NEL(NDIMM,*)
      DIMENSION AKONST(3,5,*)
C*    BLOK ZA STAMPANJE MATERIJALA PO ELEMENTIMA ZA GRAFIKU - PRIVREMENO
      JEDAN=1
      NULA=0
      ZERO=0.
      MJ=-1
      M8=8
C
      IN=13
      IGRAF=49
      WRITE(IGRAF,1000) JEDAN,IN,JEDAN
      WRITE(IGRAF,3003) 
      WRITE(IGRAF,1010) ZERO,ZERO,ZERO
      WRITE(IGRAF,1000)NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(IGRAF,1000)NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(IGRAF,1000)NULA
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
 3003 FORMAT(' KOEFICIJENT FILTRACIJE')
C*

	DO 134 NBREL=1,NE
C	
      NMAT=NEL(NDIM+1,NBREL)

      AKX=AKONST(3,1,NMAT)
C*    STAMPANJE MATERIJALA PO ELEMENTU
      IF(ISNUMER.EQ.0) THEN
        WRITE(IGRAF,5000) NBREL,AKX
      ELSEIF(ISNUMER.EQ.1) THEN
        WRITE(IGRAF,5000) MCVEL(NBREL),AKX
      ENDIF
C*
  134 CONTINUE
C
C*    KRAJ BLOKA ZA STAMPANJE PROTOKA PO ELEMENTU
      WRITE(IGRAF,5000) MJ,ZERO
 5000 FORMAT(I10,3(1PD12.4))
C*
      RETURN
      END
C=======================================================================
C stampanje rezultata u odredjenim cvorovima
C=======================================================================
      SUBROUTINE VSTAU09c(RTH,NPP,II,IND,IKOL,KOR,NZAD,NUMZAD,
     +                  KONT,POT,ISNUMBER)
      USE NODES
      USE ELEMENTS
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      CHARACTER*6 PCIME(34)
      
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
      COMMON /VOPTIM/ NKONT,MAXLIN
      DIMENSION RTH(IKOL,*),FSP(6),IN(8),NZAD(3,*)
      DIMENSION KONT(9,MAXLIN,*),POT(*)
      DIMENSION H(8)
      COMMON /PIJEZO/ CPOR(3,1000),NPIJEZ(20,100),NODP(100),NPIJ,NPOR,
     1                NPORCV(1000),NPOREL(1000),NEPOR,NPORE,NEPIJ
      COMMON /PORCEL/ RMIN(34),SMIN(34),TMIN(34)
      COMMON /DODAT/ NDIMM
C
      IF(NEPIJ.EQ.0) GO TO 100
      
C
      PCIME(1)='P07-04'
      PCIME(2)='P07-05'
      PCIME(3)='P07-07'
      PCIME(4)='P07-08'
      PCIME(5)='P07-10'
      PCIME(6)='P07-11'
      PCIME(7)='P07-12'
      PCIME(8)='P07-13'
      PCIME(9)='P11-01'
      PCIME(10)='P11-04'
      PCIME(11)='P11-05'
      PCIME(12)='P11-07'
      PCIME(13)='P11-08'
      PCIME(14)='P11-09'
      PCIME(15)='P11-10'
      PCIME(16)='P11-11'
      PCIME(17)='P11-12'
      PCIME(18)='P11-13'
      PCIME(19)='P11-14'
      PCIME(20)='P11-15'
      PCIME(21)='P11-16'
      PCIME(22)='P11-17'
      PCIME(23)='P11-18'
      PCIME(24)='P16-01'
      PCIME(25)='P16-04'
      PCIME(26)='P16-05'
      PCIME(27)='P16-06'
      PCIME(28)='P16-07'
      PCIME(29)='P16-08'
      PCIME(30)='P16-09'
      PCIME(31)='P16-10'
      PCIME(32)='P16-11'
      PCIME(33)='P16-12'
      PCIME(34)='P16-13'
C
      MJ=-1

      if(npij.gt.0) then
         WRITE(II,5000) MJ
         write(3,*) ' UCITAVANJE PODATAKA O PIJEZOMETRIMA'
         DO 12 I=1,NPIJ
            WRITE(3,*) ' PIJEZOMETAR',I,' BROJ CVOROVA',NODP(I)
            WRITE(3,*) ' CVOR,       X         Y       Z'
            DO 22 J=1,NODP(I)
               ji=NPIJEZ(J,I)
               WRITE(3,4001) ji,CORD(1,ji),CORD(2,ji),CORD(3,ji)
               if(j.gt.1) then
                  dif=CORD(3,NPIJEZ(J-1,I))-CORD(3,NPIJEZ(J,I))
                  if(dif.lt.0.) then
                     write(3,*) ' greska - pijezometar',i,' cvorovi',
     1                          NPIJEZ(J-1,I),NPIJEZ(J,I)
                     write(*,*) ' greska - pijezometar',i,' cvorovi',
     1                          NPIJEZ(J-1,I),NPIJEZ(J,I)
                     stop
                  endif
               endif
   22       CONTINUE
   12    CONTINUE
 4001    format(i10,3f10.1)
   
         DO 11 I=1,NPIJ
            DO 21 J=1,NODP(I)
               IF(RTH(1,NPIJEZ(J,I)).GT.CORD(IOSA,NPIJEZ(J,I))) THEN
                  WRITE(II,5000) I+1,RTH(1,NPIJEZ(J,I))
                  GO TO 11
               ENDIF
   21       CONTINUE
   11    CONTINUE
      endif
C
  100 IF(NEPOR.EQ.0) RETURN
C porne celije
C      if(npor.gt.0) then
C         WRITE(II,5000) MJ
C         DO 13 I=1,npor
C            WRITE(II,5000) nporcv(i),RTH(1,nporcv(i))
C   13    CONTINUE
C         WRITE(II,5000) MJ
C      endif
C porne celije ekstrapolacija
      if(NPORE.gt.0) then
         WRITE(II,5000) MJ
         DO 14 I=1,npore
C ODREDJIVANJE POTENCIJALA U PORNOJ CELIJI
              RPP=1.0+RMIN(I)
              SPP=1.0+SMIN(I)
              TPP=1.0+TMIN(I)
              RM=1.0-RMIN(I)
              TM=1.0-TMIN(I)
              SM=1.0-SMIN(I)
              RR=1.0-RMIN(I)*RMIN(I)
              SS=1.0-SMIN(I)*SMIN(I)
              TT=1.0-TMIN(I)*TMIN(I)
              DO IH=1,8
                H(IH)=0.
              ENDDO
              H(1)=0.125*RPP*SPP*TPP
              H(2)=0.125*RM*SPP*TPP
              H(3)=0.125*RM*SM*TPP
              H(4)=0.125*RPP*SM*TPP
              H(5)=0.125*RPP*SPP*TM
              H(6)=0.125*RM*SPP*TM
              H(7)=0.125*RM*SM*TM
              H(8)=0.125*RPP*SM*TM
            PPC=0.
            XPC=0.
            YPC=0.
            ZPC=0.
            DO JJH=1,8
              ncvor=NEL(JJH,NPOREL(I))          
              PPC=PPC + H(JJH)*RTH(1,ncvor)
              XPC=XPC+H(JJH)*CORD(1,ncvor)
              YPC=YPC+H(JJH)*CORD(2,ncvor)
              ZPC=ZPC+H(JJH)*CORD(3,ncvor)
              if(ppc.lt.zpc) ppc=0.D0                
C              WRITE(II,5001) I,ncvor,NPOREL(I),H(JJH),RTH(1,ncvor)
           enddo
C
            WRITE(II,5010) PCIME(I),PPC,XPC,YPC,ZPC
   14    CONTINUE
         WRITE(II,5000) MJ
      endif
C
      RETURN
 1000 FORMAT(10I6)
 1001 FORMAT(I10,F12.4)
 3001 FORMAT('TOTAL POTENTIAL')
 3004 FORMAT('DEPTH OF WATER')
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(3(1PD12.4))
 1015 FORMAT('P-',I2,I10,1PD12.4)
 1016 FORMAT('P07-',I2,I10,1PD12.4)
 1017 FORMAT('P11-',I2,I10,1PD12.4)
 1018 FORMAT('P16-',I2,I10,1PD12.4)
 5000 FORMAT(I10,3(1PE12.4))
 5010 FORMAT(A6,4(1PE12.4))
 5001 FORMAT(3(I10),2(1PE12.4))
      END
C stampanje rezultata u odredjenim cvorovima
C=======================================================================
C Stampanje Temperatura na mestu termometara
C=======================================================================
      SUBROUTINE VSTAU09CT(RTH,II,VREME,KORAK,ISNUMBER)
      USE NODES
      USE ELEMENTS
      USE MESURMENTPOINTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /DODAT/ NDIMM
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
C      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
C      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
!       COMMON /VOPTIM/ NKONT,MAXLIN,LIN(10)
      DIMENSION RTH(1,*)
!       DIMENSION NEL(NDIMM,*)
      DIMENSION H(10)
!       COMMON /PIJEZO/ CPOR(3,1000),NPIJEZ(20,100),NODP(100),NPIJ,NPOR,
!      1                NPORCV(1000),NPOREL(1000),NEPOR,NPORE,NEPIJ
!       COMMON /PORCEL/ RMIN(34),SMIN(34),TMIN(34)
C
!       IF(NEPIJ.EQ.0) GO TO 100
C
      MJ=-1
C       
      IF(MAX_MPOINTS.LE.0) RETURN
         WRITE(II,5000) MJ

 4001    format(6f10.2)
C porne celije
C      if(npor.gt.0) then
C         WRITE(II,5000) MJ
C         DO 13 I=1,npor
C            WRITE(II,5000) nporcv(i),RTH(1,nporcv(i))
C   13    CONTINUE
C         WRITE(II,5000) MJ
C      endif
C Termometri ekstrapolacija
!          WRITE(51,*) "Pijez., Potencijal"
         WRITE(II,5000) MJ
         WRITE(II,1010) VREME
!          WRITE(II,*) "ISNUMBER",ISNUMBER
         DO 14 I=1,MAX_MPOINTS
C ODREDJIVANJE POTENCIJALA na mestu termometra
C za heksa elemente
C              RPP=1.0+RMIN(I)
C              SPP=1.0+SMIN(I)
C              TPP=1.0+TMIN(I)
C              RM=1.0-RMIN(I)
C              TM=1.0-TMIN(I)
C              SM=1.0-SMIN(I)
C              RR=1.0-RMIN(I)*RMIN(I)
C              SS=1.0-SMIN(I)*SMIN(I)
C              TT=1.0-TMIN(I)*TMIN(I)
C              DO IH=1,8
C                H(IH)=0.
C              ENDDO
C              H(1)=0.125*RPP*SPP*TPP
C              H(2)=0.125*RM*SPP*TPP
C              H(3)=0.125*RM*SM*TPP
C              H(4)=0.125*RPP*SM*TPP
C              H(5)=0.125*RPP*SPP*TM
C              H(6)=0.125*RM*SPP*TM
C              H(7)=0.125*RM*SM*TM
C              H(8)=0.125*RPP*SM*TM
C TETRA ELEMENTI SA 10 CVOROVA
              RMIN=MP_COORDS(4,I)
              SMIN=MP_COORDS(5,I)
              TMIN=MP_COORDS(6,I)
              RST=1.0-RMIN-SMIN-TMIN
              DO IH=1,10
                H(IH)=0.
              ENDDO
              H(5)=4*RMIN*RST
              H(6)=4*RMIN*SMIN
              H(7)=4*SMIN*RST
              H(8)=4*TMIN*RST
              H(9)=4*RMIN*TMIN
              H(10)=4*SMIN*TMIN
              H(1)=RST-0.5*(H(5)+H(7)+H(8))
              H(2)=RMIN-0.5*(H(5)+H(6)+H(9))
              H(3)=SMIN-0.5*(H(6)+H(7)+H(10))
              H(4)=TMIN-0.5*(H(8)+H(9)+H(10))
              PPC=0.
              XPC=0.
              YPC=0.
              ZPC=0.
            DO JJH=1,10
              ncvor=NEL(JJH,MP_ELEMENT(I)) 
              
              PPC=PPC + H(JJH)*RTH(1,ncvor)
              XPC=XPC+H(JJH)*CORD(1,ncvor)
              YPC=YPC+H(JJH)*CORD(2,ncvor)
              ZPC=ZPC+H(JJH)*CORD(3,ncvor)
             IF(ISNUMBER.EQ.0) THEN
              WRITE(II,5001) I,ncvor,MP_ELEMENT(I),H(JJH),RTH(1,ncvor)
!               write(II,*) 'I,ncvor,ELEMENT',I,ncvor,MP_ELEMENT(I)
             ELSEIF(ISNUMBER.EQ.1) THEN
              ncvor1=NCVEL(ncvor)
              eleme1=MCVEL(MP_ELEMENT(I))
              WRITE(II,5001) I,ncvor1,MCVEL(MP_ELEMENT(I)),
     1                       H(JJH),RTH(1,ncvor)
!               write(II,*) 'I,ncvor,ELEM',I,ncvor,MCVEL(MP_ELEMENT(I))
!              write(*,*) 'I,ncvor1,EL',I,ncvor1,MP_ELEMENT(I),eleme1
             ENDIF
           enddo
C  NEMA ZA TEMPERATURE OVOG USLOVA       if(ppc.lt.zpc) ppc=0.D0                
C
           IF((IPAKT.EQ.0).AND.((PPC-0.0001).LT.MP_COORDS(IOSA,I)))THEN
                PPC=0.0
!             WRITE(51,5011) MPOINT_ID(I),PPC
           ENDIF
           MP_RESULTS(KORAK,I+1)=PPC
           WRITE(II,5010) MPOINT_ID(I),PPC,XPC,YPC,ZPC
           IF (IPAKT.EQ.0)  WRITE(51,5011) MPOINT_ID(I),PPC
   14    CONTINUE
         WRITE(II,5000) MJ
         MP_VREME(KORAK)=VREME
!          write(3,*) 'KORAK',KORAK
!           write(3,*) 'MP_RESULTS(KORAK,2)',MP_RESULTS(KORAK,2)
!          if(korak.eq.BRKORAKA) then
!             do ikbr=1,BRKORAKA
!               write(3,*) 'MP_RESULTS(',ikbr,',2)',MP_RESULTS(ikbr,2)
!             enddo
!          endif
C
      RETURN
 1000 FORMAT(10I6)
 1001 FORMAT(I10,F12.4)
 3001 FORMAT('TOTAL TEMPERATURE')
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(1PD12.4)
 1015 FORMAT('P-',I2,I10,1PD12.4)
 1016 FORMAT('P07-',I2,I10,1PD12.4)
 1017 FORMAT('P11-',I2,I10,1PD12.4)
 1018 FORMAT('P16-',I2,I10,1PD12.4)
 5000 FORMAT(I10,3(1PE12.4))
 5010 FORMAT(A15,4(1PE12.4))
 5011 FORMAT(A15,",",(1PE12.4,","))
 5001 FORMAT(3(I10),2(1PE12.4))
      END
C=======================================================================
C Stampanje POTENCIJALA U DODATNIM CVOROVIMA ZA CRTANJE REZULTATA
C=======================================================================
      SUBROUTINE VSTAU09CDP(RTH,II,VREME,KORAK,ISNUMBER)
      USE NODES
      USE ELEMENTS
      USE MESURMENTPOINTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /DODAT/ NDIMM
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
C      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
C      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
!       COMMON /VOPTIM/ NKONT,MAXLIN,LIN(10)
      DIMENSION RTH(1,*)
!       DIMENSION NEL(NDIMM,*)
      DIMENSION H(10)
C
      MJ=-1
C       
      IF(MAX_DPOINTS.LE.0) RETURN
         WRITE(II,5000) MJ

 4001    format(6f10.2)
C Dodatne tacke - ekstrapolacija
         WRITE(II,5000) MJ
         WRITE(II,1010) VREME
         WRITE(II,*) "POTENCIJAL"
         DO 14 I=1,MAX_DPOINTS
C ODREDJIVANJE POTENCIJALA na mestu termometra
C za heksa elemente
C              RPP=1.0+RMIN(I)
C              SPP=1.0+SMIN(I)
C              TPP=1.0+TMIN(I)
C              RM=1.0-RMIN(I)
C              TM=1.0-TMIN(I)
C              SM=1.0-SMIN(I)
C              RR=1.0-RMIN(I)*RMIN(I)
C              SS=1.0-SMIN(I)*SMIN(I)
C              TT=1.0-TMIN(I)*TMIN(I)
C              DO IH=1,8
C                H(IH)=0.
C              ENDDO
C              H(1)=0.125*RPP*SPP*TPP
C              H(2)=0.125*RM*SPP*TPP
C              H(3)=0.125*RM*SM*TPP
C              H(4)=0.125*RPP*SM*TPP
C              H(5)=0.125*RPP*SPP*TM
C              H(6)=0.125*RM*SPP*TM
C              H(7)=0.125*RM*SM*TM
C              H(8)=0.125*RPP*SM*TM
C TETRA ELEMENTI SA 10 CVOROVA
              RMIN=DP_COORDS(4,I)
              SMIN=DP_COORDS(5,I)
              TMIN=DP_COORDS(6,I)
              RST=1.0-RMIN-SMIN-TMIN
              DO IH=1,10
                H(IH)=0.
              ENDDO
              H(5)=4*RMIN*RST
              H(6)=4*RMIN*SMIN
              H(7)=4*SMIN*RST
              H(8)=4*TMIN*RST
              H(9)=4*RMIN*TMIN
              H(10)=4*SMIN*TMIN
              H(1)=RST-0.5*(H(5)+H(7)+H(8))
              H(2)=RMIN-0.5*(H(5)+H(6)+H(9))
              H(3)=SMIN-0.5*(H(6)+H(7)+H(10))
              H(4)=TMIN-0.5*(H(8)+H(9)+H(10))
            PPC=0.
            XPC=0.
            YPC=0.
            ZPC=0.
            DO JJH=1,10
              ncvor=NEL(JJH,DP_ELEMENT(I))          
              PPC=PPC + H(JJH)*RTH(1,ncvor)
              XPC=XPC+H(JJH)*CORD(1,ncvor)
              YPC=YPC+H(JJH)*CORD(2,ncvor)
              ZPC=ZPC+H(JJH)*CORD(3,ncvor)
           enddo
            if(PPC.lt.ZPC) PPC=0.D0                
!               WRITE(II,5001) I,ncvor,DP_ELEMENT(I),H(JJH),RTH(1,ncvor)
C
         IF((IPAKT.EQ.0).AND.((PPC-0.0001).LT.DP_COORDS(3,I))) THEN
              PPC=0.0
            ENDIF
            DP_RESULTS(KORAK,1,I+1)=PPC
            WRITE(II,5010) DPOINT_ID(I),PPC
!             WRITE(II,5010) DPOINT_ID(I),PPC,XPC,YPC,ZPC
   14    CONTINUE
         WRITE(II,5000) MJ
!          DP_VREME(KORAK)=VREME
!          write(3,*) 'KORAK',KORAK
!           write(3,*) 'MP_RESULTS(KORAK,2)',MP_RESULTS(KORAK,2)
!          if(korak.eq.BRKORAKA) then
!             do ikbr=1,BRKORAKA
!               write(3,*) 'MP_RESULTS(',ikbr,',2)',MP_RESULTS(ikbr,2)
!             enddo
!          endif
C
      RETURN
 1000 FORMAT(10I6)
 1001 FORMAT(I10,F12.4)
 3001 FORMAT('POTENTIAL')
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(1PD12.4)
 1015 FORMAT('P-',I2,I10,1PD12.4)
 1016 FORMAT('P07-',I2,I10,1PD12.4)
 1017 FORMAT('P11-',I2,I10,1PD12.4)
 1018 FORMAT('P16-',I2,I10,1PD12.4)
 5000 FORMAT(I10,3(1PE12.4))
 5010 FORMAT(A15,4(1PE12.4))
 5001 FORMAT(3(I10),2(1PE12.4))
      END
C=======================================================================
C=======================================================================
C Stampanje brzine U DODATNIM TACKAMA ZA REZULTATE
C=======================================================================
      SUBROUTINE VSTAU09CDV(RTH,II,VREME,KORAK,NEL,ISNUMBER)
      USE NODES
      USE MESURMENTPOINTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /DODAT/ NDIMM
C      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
C      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
!       COMMON /VOPTIM/ NKONT,MAXLIN,LIN(10)
      DIMENSION RTH(3,*)
      DIMENSION NEL(NDIMM,*)
      DIMENSION H(10)
!       COMMON /PIJEZO/ CPOR(3,1000),NPIJEZ(20,100),NODP(100),NPIJ,NPOR,
!      1                NPORCV(1000),NPOREL(1000),NEPOR,NPORE,NEPIJ
!       COMMON /PORCEL/ RMIN(34),SMIN(34),TMIN(34)
C
!       IF(NEPIJ.EQ.0) GO TO 100
C
      MJ=-1
C       
      IF(MAX_DPOINTS.LE.0) RETURN
         WRITE(II,5000) MJ

 4001    format(6f10.2)
C porne celije
C      if(npor.gt.0) then
C         WRITE(II,5000) MJ
C         DO 13 I=1,npor
C            WRITE(II,5000) nporcv(i),RTH(1,nporcv(i))
C   13    CONTINUE
C         WRITE(II,5000) MJ
C      endif
C Termometri ekstrapolacija
         WRITE(II,5000) MJ
         WRITE(II,1010) VREME
         WRITE(II,*) "Brzine"
         DO 14 I=1,MAX_DPOINTS
C ODREDJIVANJE POTENCIJALA na mestu termometra
C za heksa elemente
C              RPP=1.0+RMIN(I)
C              SPP=1.0+SMIN(I)
C              TPP=1.0+TMIN(I)
C              RM=1.0-RMIN(I)
C              TM=1.0-TMIN(I)
C              SM=1.0-SMIN(I)
C              RR=1.0-RMIN(I)*RMIN(I)
C              SS=1.0-SMIN(I)*SMIN(I)
C              TT=1.0-TMIN(I)*TMIN(I)
C              DO IH=1,8
C                H(IH)=0.
C              ENDDO
C              H(1)=0.125*RPP*SPP*TPP
C              H(2)=0.125*RM*SPP*TPP
C              H(3)=0.125*RM*SM*TPP
C              H(4)=0.125*RPP*SM*TPP
C              H(5)=0.125*RPP*SPP*TM
C              H(6)=0.125*RM*SPP*TM
C              H(7)=0.125*RM*SM*TM
C              H(8)=0.125*RPP*SM*TM
C TETRA ELEMENTI SA 10 CVOROVA
              RMIN=DP_COORDS(4,I)
              SMIN=DP_COORDS(5,I)
              TMIN=DP_COORDS(6,I)
              RST=1.0-RMIN-SMIN-TMIN
              DO IH=1,10
                H(IH)=0.
              ENDDO
              H(5)=4*RMIN*RST
              H(6)=4*RMIN*SMIN
              H(7)=4*SMIN*RST
              H(8)=4*TMIN*RST
              H(9)=4*RMIN*TMIN
              H(10)=4*SMIN*TMIN
              H(1)=RST-0.5*(H(5)+H(7)+H(8))
              H(2)=RMIN-0.5*(H(5)+H(6)+H(9))
              H(3)=SMIN-0.5*(H(6)+H(7)+H(10))
              H(4)=TMIN-0.5*(H(8)+H(9)+H(10))
            PPC=0.
            XPC=0.
            YPC=0.
            ZPC=0.
            DO JJH=1,10
              ncvor=NEL(JJH,DP_ELEMENT(I))          
              PPC=PPC + H(JJH)*RTH(2,ncvor)
              XPC=XPC+H(JJH)*CORD(1,ncvor)
              YPC=YPC+H(JJH)*CORD(2,ncvor)
              ZPC=ZPC+H(JJH)*CORD(3,ncvor)
!               if(PPC.lt.ZPC) PPC=0.D0                
C              WRITE(II,5001) I,ncvor,DP_ELEMENT(I),H(JJH),RTH(1,ncvor)
           enddo
C
            DP_RESULTS(KORAK,2,I+1)=PPC
            WRITE(II,5010) DPOINT_ID(I),PPC
!             WRITE(II,5010) DPOINT_ID(I),PPC,XPC,YPC,ZPC
   14    CONTINUE
         WRITE(II,5000) MJ
!          DP_VREME(KORAK)=VREME
!          write(3,*) 'KORAK',KORAK
!           write(3,*) 'DP_RESULTS(KORAK,2)',DP_RESULTS(KORAK,2)
!          if(korak.eq.BRKORAKA) then
!             do ikbr=1,BRKORAKA
!               write(3,*) 'DP_RESULTS(',ikbr,',2)',DP_RESULTS(ikbr,2)
!             enddo
!          endif
C
      RETURN
 1000 FORMAT(10I6)
 1001 FORMAT(I10,F12.4)
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(1PD12.4)
 1015 FORMAT('P-',I2,I10,1PD12.4)
 1016 FORMAT('P07-',I2,I10,1PD12.4)
 1017 FORMAT('P11-',I2,I10,1PD12.4)
 1018 FORMAT('P16-',I2,I10,1PD12.4)
 5000 FORMAT(I10,3(1PE12.4))
 5010 FORMAT(A10,4(1PE12.4))
 5001 FORMAT(3(I10),2(1PE12.4))
      END
C=======================================================================
C Stampanje POTENCIJALA i BRZINE U DODATNIM CVOROVIMA ZA CRTANJE REZULTATA
C=======================================================================
      SUBROUTINE VSTAU09CDPV(RTH,RTH1,II,VREME,KORAK,ISNUMBER)
      USE NODES
      USE ELEMENTS
      USE MESURMENTPOINTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /DODAT/ NDIMM
C      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
C      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
!       COMMON /VOPTIM/ NKONT,MAXLIN,LIN(10)
      DIMENSION RTH(1,*),RTH1(3,*)
!       DIMENSION NEL(NDIMM,*)
      DIMENSION H(10)
C
      MJ=-1
C       
      IF(MAX_DPOINTS.LE.0) RETURN
!          WRITE(II,5000) MJ

 4001    format(6f10.2)
C Dodatne tacke - ekstrapolacija
!          WRITE(II,5000) MJ
!          WRITE(II,1010) VREME
         WRITE(II,*) "ID tacke, POTENCIJAL, BRZINA"
!          WRITE(3,*) "ID tacke, POTENCIJAL, BRZINA"
        DO 14 I=1,MAX_DPOINTS
C ODREDJIVANJE POTENCIJALA na mestu termometra
C za heksa elemente
C              RPP=1.0+RMIN(I)
C              SPP=1.0+SMIN(I)
C              TPP=1.0+TMIN(I)
C              RM=1.0-RMIN(I)
C              TM=1.0-TMIN(I)
C              SM=1.0-SMIN(I)
C              RR=1.0-RMIN(I)*RMIN(I)
C              SS=1.0-SMIN(I)*SMIN(I)
C              TT=1.0-TMIN(I)*TMIN(I)
C              DO IH=1,8
C                H(IH)=0.
C              ENDDO
C              H(1)=0.125*RPP*SPP*TPP
C              H(2)=0.125*RM*SPP*TPP
C              H(3)=0.125*RM*SM*TPP
C              H(4)=0.125*RPP*SM*TPP
C              H(5)=0.125*RPP*SPP*TM
C              H(6)=0.125*RM*SPP*TM
C              H(7)=0.125*RM*SM*TM
C              H(8)=0.125*RPP*SM*TM
C TETRA ELEMENTI SA 10 CVOROVA
              RMIN=DP_COORDS(4,I)
              SMIN=DP_COORDS(5,I)
              TMIN=DP_COORDS(6,I)
              RST=1.0-RMIN-SMIN-TMIN
              DO IH=1,10
                H(IH)=0.
              ENDDO
              H(5)=4*RMIN*RST
              H(6)=4*RMIN*SMIN
              H(7)=4*SMIN*RST
              H(8)=4*TMIN*RST
              H(9)=4*RMIN*TMIN
              H(10)=4*SMIN*TMIN
              H(1)=RST-0.5*(H(5)+H(7)+H(8))
              H(2)=RMIN-0.5*(H(5)+H(6)+H(9))
              H(3)=SMIN-0.5*(H(6)+H(7)+H(10))
              H(4)=TMIN-0.5*(H(8)+H(9)+H(10))
            PPC=0.
            PPC1=0.
            PPC2=0.
            PPC3=0.
            XPC=0.
            YPC=0.
            ZPC=0.
            DO JJH=1,10
              ncvor=NEL(JJH,DP_ELEMENT(I))          
              PPC=PPC + H(JJH)*RTH(1,ncvor)
              PPC1=PPC1 + H(JJH)*RTH1(1,ncvor)
              PPC2=PPC2 + H(JJH)*RTH1(2,ncvor)
              PPC3=PPC3 + H(JJH)*RTH1(3,ncvor)
               XPC=XPC+H(JJH)*CORD(1,ncvor)
               YPC=YPC+H(JJH)*CORD(2,ncvor)
               ZPC=ZPC+H(JJH)*CORD(3,ncvor)
           enddo
           PPCV=SQRT(PPC1*PPC1+PPC2*PPC2+PPC3*PPC3)
            if(PPC.lt.DP_COORDS(3,I)) THEN
              PPC=DP_COORDS(3,I)                
              PPCV=0.D0                
            endif
!             WRITE(3,*) I,DP_ELEMENT(I),XPC,YPC,ZPC,PPC,PPCV
C
            DP_RESULTS(KORAK,1,I+1)=PPC
            DP_RESULTS(KORAK,2,I+1)=PPCV
            WRITE(II,5010) DPOINT_ID(I),PPC,PPCV
!             WRITE(II,5010) DPOINT_ID(I),PPC,XPC,YPC,ZPC
   14    CONTINUE
!          WRITE(II,5000) MJ
!          DP_VREME(KORAK)=VREME
!          write(3,*) 'KORAK',KORAK
!           write(3,*) 'MP_RESULTS(KORAK,2)',MP_RESULTS(KORAK,2)
!          if(korak.eq.BRKORAKA) then
!             do ikbr=1,BRKORAKA
!               write(3,*) 'MP_RESULTS(',ikbr,',2)',MP_RESULTS(ikbr,2)
!             enddo
!          endif
C
      RETURN
 1000 FORMAT(10I6)
 1001 FORMAT(I10,F12.4)
 3001 FORMAT('POTENTIAL')
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(1PD12.4)
 5010 FORMAT(A15,",",2(1PE12.4,","))
 5001 FORMAT(3(I10),2(1PE12.4))
      END
C=======================================================================
C=======================================================================
      SUBROUTINE VSTAU35(SIGMA,KOR,ISNUMBER)
      USE ELEMENTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /DODAT/ NDIMM
      COMMON /NIDEAS/ IDEAS
      COMMON /SRPSKI/ ISRPS
C
      DIMENSION SIGMA(NET,8),S(7),N(20)
C
      NNCVE=NDIM
      NCVE=NDIM
      IF(NDIM.LT.20) NDIM=8
      IF(NDIM.EQ.21) NDIM=20
C      
      JEDAN=1
      NULA=0
      ZERO=0.
      ONE=1.
      MJ=-1
      M451=451
      M2=2
      M3=3
      M4=4
      M5=5
      M8=8
      M9=9
C
      II=49
      NGAUSX=2
      NGAUSY=2
      NGAUSZ=2
      ISTNA=3
C
c / DA LI TREBA OVO      NPROS=MODPR2(NMODM)
c      NGXYZ=NGS12*NPROS*MXS
C
   88 IND=60031
      NND=70031
      CALL VNAPUN3(N,NDIM,NND)
      KK=1
      INP=5
      GO TO 10
C
C
   10 WRITE(II,1000) KOR,IND,JEDAN
      IF(KK.EQ.1.AND.INP.NE.0) THEN
         IF(IND.EQ.60031) WRITE(II,60031) 
      ENDIF   
      WRITE(II,1010) ZERO,ZERO,ZERO
      IF(ISTNA.EQ.0) THEN
       WRITE(II,1000) N(1),N(2),N(3),N(4),N(5),N(6),N(7),N(8),N(9),N(10)
       WRITE(II,1000) 
     +       N(11),N(12),N(13),N(14),N(15),N(16),N(17),N(18),N(19),N(20)
      ELSE
        WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
        WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      ENDIF
      INDM=MOD(IND,100)
C dodat sledeci red za verziju neu fajla 10.0 /Sneza 26.12.2016.
      WRITE(II,1000) NULA
      IF(INDM.GE.10.AND.INDM.LT.50) WRITE(II,1000) NULA,NULA,M4,M8
      IF(INDM.GE.50.AND.INDM.LT.99) WRITE(II,1000) NULA,NULA,M5,M8
      IF(ISTNA.EQ.0) THEN
         WRITE(II,1000) NULA,NULA,JEDAN
      ELSE
         WRITE(II,1000) JEDAN,NULA,JEDAN
      ENDIF
C      
      DO 40 NLM=1,NET
          IF(elemtip(NLM).gt.100) THEN
            NTYPE=elemtip(NLM)/100
            NTIMES=100
          ELSE
            NTYPE=elemtip(NLM)/10
            NTIMES=10
          ENDIF
          NDIMEL=elemtip(NLM)-NTIMES*NTYPE
C
          NMM=NLM
!           IF(ISNUMBER.EQ.1) NMM=MCVEL(NLM)
C         MAT=NMAT(NLM)
          MAT=NEL(NDIM+1,NLM)
          NGS12=8
        SELECT CASE (elemtip(NLM))
        CASE (12:13)
          NGS12=2
        CASE (23,26)
          NGS12=3
        CASE (24,28,34,310)
          NGS12=4
         CASE DEFAULT
C          WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NLM)
        END SELECT
C
          CALL VSREDK3(SIGMA,S,NLM,NGS12,NET)
        IF(ISNUMBER.EQ.1) NMM=MCVEL(NLM)
        IF(KK.EQ.1.AND.INP.NE.0) THEN
           IF(IND.EQ.60031.AND.S(7).GE.0.D0) then
            WRITE(II,5000) NMM,S(7)
           ELSE
              WRITE(II,5000) NMM,0.D0
           ENDIF
         ENDIF
   40 CONTINUE
      WRITE(II,5500) MJ,ZERO
      IF(ISTNA.GT.0) GO TO 31
C
      DO 30 NC=1,NDIMEL
      KTMP=1
      DO 20 NLM=1,NET
C  preskace stampanje kod 1D i 2D elemenata
!           WRITE(*,*) NC,NLM,elemtip(NLM)
       SELECT CASE (elemtip(NLM))
         CASE (12:28)
          GOTO 20
         CASE (34,310)
          IPODT=3
          NGAUSX=4
         CASE DEFAULT
C          WRITE(*,*) 'ELEMENT TYPE UNKNOWN NETIP=',elemtip(NLM)
      END SELECT
        IF(KTMP.eq.1) THEN
          IF(elemtip(NLM).gt.100) THEN
            NTYPE=elemtip(NLM)/100
            NTIMES=100
          ELSE
            NTYPE=elemtip(NLM)/10
            NTIMES=10
          ENDIF
          NDIMEL=elemtip(NLM)-NTIMES*NTYPE
          WRITE(II,1000) KOR,N(NC),JEDAN
          IF(KK.EQ.1.AND.INP.NE.0) THEN
            IF(IND.EQ.60031) WRITE(II,70031) NC 
          ENDIF   
          WRITE(II,1010) ZERO,ZERO,ZERO
          WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
          WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
          IF(INDM.GE.10.AND.INDM.LT.50) WRITE(II,1000) NULA,NULA,M4,M8
          IF(INDM.GE.50.AND.INDM.LT.99) WRITE(II,1000) NULA,NULA,M5,M8
          WRITE(II,1000) NULA,NULA,NULA
          KTMP=0
        ENDIF
        
        IF(NDIMEL.EQ.8.AND.NGAUSX.EQ.2) THEN
          IF(NC.EQ.1) IJ=NGAUSY*NGAUSZ+NGAUSZ+1
          IF(NC.EQ.2) IJ=NGAUSZ+1
          IF(NC.EQ.3) IJ=1
          IF(NC.EQ.4) IJ=NGAUSY*NGAUSZ+1
          IF(NC.EQ.5) IJ=NGAUSY*NGAUSZ+NGAUSZ+2
          IF(NC.EQ.6) IJ=NGAUSZ+2
          IF(NC.EQ.7) IJ=2
          IF(NC.EQ.8) IJ=NGAUSY*NGAUSZ+2
        ENDIF
        IF(NDIMEL.EQ.8.AND.NGAUSX.EQ.3) THEN
          IF(NC.EQ.1) IJ=2*NGAUSY*NGAUSZ+2*NGAUSZ+1
          IF(NC.EQ.2) IJ=2*NGAUSZ+1
          IF(NC.EQ.3) IJ=1
          IF(NC.EQ.4) IJ=2*NGAUSY*NGAUSZ+1
          IF(NC.EQ.5) IJ=2*NGAUSY*NGAUSZ+2*NGAUSZ+3
          IF(NC.EQ.6) IJ=2*NGAUSZ+3
          IF(NC.EQ.7) IJ=3
          IF(NC.EQ.8) IJ=2*NGAUSY*NGAUSZ+3
        ENDIF
        IF(NDIMEL.EQ.20.AND.NGAUSX.EQ.2) THEN
          IF(NC.EQ.1) IJ=NGAUSY*NGAUSZ+NGAUSZ+1
          IF(NC.EQ.2) IJ=NGAUSZ+1
          IF(NC.EQ.3) IJ=1
          IF(NC.EQ.4) IJ=NGAUSY*NGAUSZ+1
          IF(NC.EQ.5) IJ=NGAUSY*NGAUSZ+NGAUSZ+2
          IF(NC.EQ.6) IJ=NGAUSZ+2
          IF(NC.EQ.7) IJ=2
          IF(NC.EQ.8) IJ=NGAUSY*NGAUSZ+2
          IF(NC.EQ.9) THEN
              IJ=-(NGAUSY*NGAUSZ+NGAUSZ+1)
              JI=NGAUSZ+1
          ENDIF
          IF(NC.EQ.10) THEN
              IJ=-(NGAUSZ+1)
              JI=1
          ENDIF
          IF(NC.EQ.11) THEN
              IJ=-1
              JI=NGAUSY*NGAUSZ+1
          ENDIF
          IF(NC.EQ.12) THEN
              IJ=-(NGAUSY*NGAUSZ+1)
              JI=NGAUSY*NGAUSZ+NGAUSZ+1
          ENDIF
          IF(NC.EQ.13) THEN
              IJ=-(NGAUSY*NGAUSZ+NGAUSZ+1)
              JI=NGAUSY*NGAUSZ+NGAUSZ+2
          ENDIF
          IF(NC.EQ.14) THEN
              IJ=-(NGAUSZ+1)
              JI=NGAUSZ+2
          ENDIF
          IF(NC.EQ.15) THEN
              IJ=-1
              JI=2
          ENDIF
          IF(NC.EQ.16) THEN
              IJ=-(NGAUSY*NGAUSZ+1)
              JI=NGAUSY*NGAUSZ+2
          ENDIF
          IF(NC.EQ.17) THEN
              IJ=-(NGAUSY*NGAUSZ+NGAUSZ+2)
              JI=NGAUSZ+2
          ENDIF
          IF(NC.EQ.18) THEN
              IJ=-(NGAUSZ+2)
              JI=2
          ENDIF
          IF(NC.EQ.19) THEN
              IJ=-2
              JI=NGAUSY*NGAUSZ+2
          ENDIF
          IF(NC.EQ.20) THEN
              IJ=-(NGAUSY*NGAUSZ+2)
              JI=NGAUSY*NGAUSZ+NGAUSZ+2
          ENDIF
        ENDIF
        IF(NDIMEL.EQ.20.AND.NGAUSX.EQ.3) THEN
          IF(NC.EQ.1) IJ=2*NGAUSY*NGAUSZ+2*NGAUSZ+1
          IF(NC.EQ.2) IJ=2*NGAUSZ+1
          IF(NC.EQ.3) IJ=1
          IF(NC.EQ.4) IJ=2*NGAUSY*NGAUSZ+1
          IF(NC.EQ.5) IJ=2*NGAUSY*NGAUSZ+2*NGAUSZ+3
          IF(NC.EQ.6) IJ=2*NGAUSZ+3
          IF(NC.EQ.7) IJ=3
          IF(NC.EQ.8) IJ=2*NGAUSY*NGAUSZ+3
          IF(NC.EQ.9) IJ=NGAUSY*NGAUSZ+2*NGAUSZ+1
          IF(NC.EQ.10) IJ=NGAUSZ+1
          IF(NC.EQ.11) IJ=NGAUSY*NGAUSZ+1
          IF(NC.EQ.12) IJ=2*NGAUSY*NGAUSZ+NGAUSZ+1
          IF(NC.EQ.13) IJ=2*NGAUSY*NGAUSZ+2*NGAUSZ+2
          IF(NC.EQ.14) IJ=2*NGAUSZ+2
          IF(NC.EQ.15) IJ=2
          IF(NC.EQ.16) IJ=2*NGAUSY*NGAUSZ+2
          IF(NC.EQ.17) IJ=NGAUSY*NGAUSZ+2*NGAUSZ+3
          IF(NC.EQ.18) IJ=NGAUSZ+3
          IF(NC.EQ.19) IJ=NGAUSY*NGAUSZ+3
          IF(NC.EQ.20) IJ=2*NGAUSY*NGAUSZ+NGAUSZ+3
        ENDIF
c    
       IF(IPODT.EQ.3) THEN
C          WRITE(*,*) 'NDIMEL,NGAUSX,NLM',NDIMEL,NGAUSX,NLM
       IF(NDIMEL.EQ.4.AND.NGAUSX.EQ.1) THEN
         IF(NC.EQ.1) IJ=1
         IF(NC.EQ.2) IJ=1
         IF(NC.EQ.3) IJ=1
         IF(NC.EQ.4) IJ=1
       ENDIF
       IF(NDIMEL.EQ.10.AND.NGAUSX.EQ.4) THEN
         IF(NC.EQ.1) IJ=1
         IF(NC.EQ.2) IJ=2
         IF(NC.EQ.3) IJ=3
         IF(NC.EQ.4) IJ=4
         IF(NC.EQ.5) THEN
            IJ=-1
            JI= 2
         ENDIF
         IF(NC.EQ.6) THEN
            IJ=-2
            JI= 3
         ENDIF
         IF(NC.EQ.7) THEN
            IJ=-1
            JI= 3
         ENDIF
         IF(NC.EQ.8) THEN
            IJ=-1
            JI= 4
         ENDIF
         IF(NC.EQ.9) THEN
            IJ=-2
            JI= 4
         ENDIF
         IF(NC.EQ.10) THEN
            IJ=-3
            JI= 4
         ENDIF
       ENDIF
       IF(NDIMEL.EQ.10.AND.NGAUSX.EQ.5) THEN
         IF(NC.EQ.1) IJ=2
         IF(NC.EQ.2) IJ=3
         IF(NC.EQ.3) IJ=4
         IF(NC.EQ.4) IJ=5
         IF(NC.EQ.5) THEN
            IJ=-3
            JI= 2
         ENDIF
         IF(NC.EQ.6) THEN
            IJ=-4
            JI= 3
         ENDIF
         IF(NC.EQ.7) THEN
            IJ=-2
            JI= 4
         ENDIF
         IF(NC.EQ.8) THEN
            IJ=-2
            JI= 5
         ENDIF
         IF(NC.EQ.9) THEN
            IJ=-3
            JI= 5
         ENDIF
         IF(NC.EQ.10) THEN
            IJ=-5
            JI= 4
         ENDIF
       ENDIF
       ENDIF
C
          NMM=NLM
          IF(ISNUMBER.EQ.1) NMM=MCVEL(NLM)
            IF(KK.EQ.1) THEN
               IF(IND.EQ.60031.AND.IJ.GT.0) THEN
                 IF(SIGMA(NLM,IJ).GE.0)WRITE(II,5000) NMM,SIGMA(NLM,IJ)
               ELSEIF(IND.EQ.60031.AND.IJ.LT.0) THEN
                  POLA=(SIGMA(NLM,-IJ)+SIGMA(NLM,JI))/2
                  IF(SIGMA(NLM,IJ).GE.0)WRITE(II,5000) NMM,POLA
               ENDIF
            ENDIF  
C
   20 CONTINUE
      WRITE(II,5500) MJ,ZERO
   30 CONTINUE
C
   31 CONTINUE
C      IF(KK.EQ.1.AND.INP.NE.0) THEN
C         IF(IND.EQ.60031) GO TO 182
C      ENDIF
  999 RETURN
 1000 FORMAT(10(I10))
 1001 FORMAT(I10,',',F12.4,',')
 1500 FORMAT(I5)
60031 FORMAT('PORNI PRITISAK')
70031 FORMAT('SOLID NODE',I10,' PORNI PRITISAK')
C
!  1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PE12.4)
 1010 FORMAT(3(1PE12.4,','))
 5000 FORMAT(I10,3(1PE12.4))
!  5500 FORMAT(I2,',',3(1PE12.4,','))
 5500 FORMAT(I10,3(1PD12.4))
 
      END
C=======================================================================
C=======================================================================
      SUBROUTINE VNAPUN3(N,NCVE,NND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     
      DIMENSION N(*)
C
      CALL VICLEAR(N,20)
      IDEL=200
      NND=NND-IDEL
      DO 10 I=1,NCVE
         N(I)=NND+I*IDEL
   10 CONTINUE
      RETURN
      END
C=======================================================================
C=====================================================================
      SUBROUTINE VSREDK25(SIGMA,S,NPROS,NGS12,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     
      DIMENSION SIGMA(*),S(7)
C
      MINUS=0
      IF(IND.LT.0) THEN
         MINUS=1
         IND=-IND
      ENDIF
      CALL VCLEAR(S,7)
      DO 10 I=1,NGS12
         GO TO (1,2,3,4,5)IND
    1    S(1)=S(1)+SIGMA(1+(I-1)*NPROS)/NGS12
         GO TO 10
    2    S(2)=S(2)+SIGMA(2+(I-1)*NPROS)/NGS12
         GO TO 10
    3    S(3)=S(3)+SIGMA(3+(I-1)*NPROS)/NGS12
         GO TO 10
    4    S(4)=S(4)+SIGMA(4+(I-1)*NPROS)/NGS12
         GO TO 10
    5    S(1)=S(1)+SIGMA(1+(I-1)*NPROS)/NGS12
         S(2)=S(2)+SIGMA(2+(I-1)*NPROS)/NGS12
         S(3)=S(3)+SIGMA(3+(I-1)*NPROS)/NGS12
         S(4)=S(4)+SIGMA(4+(I-1)*NPROS)/NGS12
   10 CONTINUE
C
      IF(IND.EQ.5) THEN
         IF(MINUS.EQ.0) THEN
            XMY=S(1)-S(2)
            YMZ=S(2)-S(4)
            ZMX=S(4)-S(1)
            EFEK=0.5D0*(XMY*XMY+YMZ*YMZ+ZMX*ZMX)+
     1           3.D0*(S(3)*S(3))
            IF(EFEK.GT.1.D-19) S(5)=DSQRT(EFEK)
         ELSE
            SR=(S(1)+S(2)+S(4))/3.
            SD1=S(1)-SR
            SD2=S(2)-SR
            SD3=S(4)-SR
            EFEK=2.*(SD1*SD1+SD2*SD2+SD3*SD3+
     1               S(4)*S(4))/3.
            IF(EFEK.GT.1.D-19) S(5)=DSQRT(EFEK)
         ENDIF
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE VSREDK3(SIGMA,S,NLM,NGS12,NE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     
      DIMENSION SIGMA(NE,*),S(7)
C
      CALL VCLEAR(S,7)
      DO 10 I=1,NGS12
         S(1)=S(1)+SIGMA(NLM,I)
         S(2)=S(2)+SIGMA(NLM,I)
         S(3)=S(3)+SIGMA(NLM,I)
         S(4)=S(4)+SIGMA(NLM,I)
         S(5)=S(5)+SIGMA(NLM,I)
         S(6)=S(6)+SIGMA(NLM,I)
   10 CONTINUE
      S(1)=S(1)/NGS12
      S(2)=S(2)/NGS12
      S(3)=S(3)/NGS12
      S(4)=S(4)/NGS12
      S(5)=S(5)/NGS12
      S(6)=S(6)/NGS12
      S(7)=S(1)
c      XMY=S(1)-S(2)
c      YMZ=S(2)-S(3)
c      ZMX=S(3)-S(1)
c      EFEK=0.5D0*(XMY*XMY+YMZ*YMZ+ZMX*ZMX)+
c     1     3.D0*(S(4)*S(4)+S(5)*S(5)+S(6)*S(6))
c      IF(EFEK.GT.1.D-19) S(7)=DSQRT(EFEK)

      RETURN
      END
C======================================================================
C ZA TEMPERATURE
C======================================================================
      SUBROUTINE VSTAS09T(NPP,II,IND,IKOL,KOR,NGPSIL,MAXSIL,ISNUMBER)
       USE NODES
       USE ELEMENTS
       USE PREDISCRIBED
       USE pflux
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
!       COMMON /VOPTIM/ NKONT,MAXLIN
      DIMENSION FSP(6),IN(8),NGPSIL(12,*)
      integer NGPSIL
!       DIMENSION KONT(9,MAXLIN,*)
!       DIMENSION KONT1(40000)
C
      JEDAN=1
      NULA=0
      ZERO=0.
      ONE=1.
      MJ=-1
      M451=451
      M1=1
      M2=2
      M3=3
      M4=4
      M5=5
      M6=6
      M7=7
      M8=8
      M9=9
C      II=49
C
      IN(1)=IND
      IN(2)=IND+1
      IN(3)=IND+2
      IN(4)=IND+3
      IN(5)=IND+4
      IN(6)=IND+5
      IN(7)=IND+6
      IN(8)=IND+7
C
C PRIVREMENO
!       IF(IN(1).EQ.10) THEN
!       DO I=1,40000
!          KONT1(I)=0
!       ENDDO
!       ENDIF
C      
!       Write(*,*) "neu, ind, maxsil",IND,MAXSIL
      IP=0
C
   30 IP=IP+1
      WRITE(II,1000) KOR,IN(IP),JEDAN
      IF(IN(1).EQ.1.AND.IP.EQ.1) WRITE(II,3001) 
      IF(IN(1).EQ.19.AND.IP.EQ.1) WRITE(II,3002) 
      IF(IN(1).EQ.10.AND.IP.EQ.1) WRITE(II,3003) 
      IF(IN(1).EQ.12.AND.IP.EQ.1) WRITE(II,3012) 
      IF(IN(1).EQ.41.AND.IP.EQ.1) WRITE(II,3041) 
      IF(IN(1).EQ.41.AND.IP.EQ.2) WRITE(II,3042) 
      IF(IN(1).EQ.41.AND.IP.EQ.3) WRITE(II,3043) 
      IF(IN(1).EQ.41.AND.IP.EQ.4) WRITE(II,3044) 
      IF(IN(1).EQ.51.AND.IP.EQ.1) WRITE(II,3051) 
      IF(IN(1).EQ.51.AND.IP.EQ.2) WRITE(II,3052) 
      IF(IN(1).EQ.51.AND.IP.EQ.3) WRITE(II,3053) 
      IF(IN(1).EQ.51.AND.IP.EQ.4) WRITE(II,3054) 
      WRITE(II,1010) ZERO,ZERO,ZERO
C OVO SU NULE ZA TEMPERATURE
      IF(IP.EQ.1.AND.IN(1).LT.40)
     +WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.1.AND.IN(1).GT.40)
     +WRITE(II,1000)IN(2),IN(3),IN(4),NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.2)
     +  WRITE(II,1000)IN(2),NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.3)
     +  WRITE(II,1000)NULA,IN(3),NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.4)
     +  WRITE(II,1000)NULA,NULA,IN(4),NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
C
C 0=ANY,1=DISPL,2=ACCEL,3=FORCE,4=STRESS,5=STRAIN,6=TEMP,OTHER=USER
C     NODAL(7) ,ELEMENTAL(8)
C dodat sledeci red za verziju neu fajla 10.0 /Sneza 26.12.2016.
      WRITE(II,1000) NULA
      IF(IN(1).EQ.1) WRITE(II,1000) NULA,NULA,M6,M7
      IF(IN(1).EQ.19) WRITE(II,1000) NULA,NULA,M8,M8
      IF(IN(1).EQ.10) WRITE(II,1000) NULA,NULA,M8,M8
      IF(IN(1).EQ.12) WRITE(II,1000) NULA,NULA,M8,M8
      IF(IN(1).EQ.41.OR.IN(1).EQ.51) WRITE(II,1000) NULA,NULA,M3,M7
C
C =1, CAN NOT LINEARY COMBINE THIS OUTPUT
C =1, COMP(0-2) ARE THE X,Y,Z COMPONENT
C =1, THIS VECTOR HAS CONTROIDAL OR NODAL OUTPUT
      IF(IN(1).GT.40) THEN
         IF(IP.EQ.1) THEN
            WRITE(II,1000) JEDAN,JEDAN,JEDAN
         ELSE
            WRITE(II,1000) NULA,JEDAN,JEDAN
         ENDIF
      ELSE
         WRITE(II,1000) JEDAN,NULA,JEDAN
      ENDIF
C      DO 10 I=1,NPP
         DO 20 J=1,IKOL
            FSP(J) = 0.0D0
            IF(IKOL.EQ.1) THEN
               IF(IN(1).EQ.1) THEN
                  IF(ID(J,I).EQ.0) GO TO 20
                  K = ID(J,I)
                  IF(K.GT.0) THEN
C OVO JE SAMO ZA 3D SLOBODNA POVRSINA
C                     IF(RTH(1,K).GT.CORD(3,I)) FSP(J)=RTH(1,K)
!                      FSP(J)=RTH(1,K)
                  ENDIF
!                ELSEIF(IN(1).EQ.9) THEN
!                  FSP(J)=1.0D0
!                  DO L=1,NUMZAD 
!                     IF(NZAD(1,L).EQ.I) GO TO 1201
!                  ENDDO
!                  FSP(J)=0.D0
!  1201            CONTINUE
               ELSE
!                  Write(*,*) "neu, ind, maxsil",IND,MAXSIL
!                  DO LL=1,NKONT 
                 DO LL=1,MAXSIL
!                     FSP(J)=LL
!                  DO L=1,LIN(LL) 
C                 DO LI=2,5 
C                 DO LI=2,NNODE+1 
C                    IF(KONT(LI,L,LL).EQ.I) THEN
C   zadat flux-zracenje
                   IF(IND.EQ.10) THEN
                     IF(ISNUMBER.EQ.0) THEN
                       IEL=NGPSIL(1,LL)
                     ELSE
                       IEL=MCVEL(NGPSIL(1,LL))
                     ENDIF
                     VAL=PFLUXEL(LL)
                   ELSEIF(IND.EQ.12) THEN
C   zadata  prelaznost, stampa ID funkcije za koeficijent prelaznosti
                     IF(ISNUMBER.EQ.0) THEN
                       IEL=NELTOK(1,LL)
                     ELSE
                       IEL=MCVEL(NELTOK(1,LL))
                     ENDIF
!                    Write(*,*) IEL
                    VAL=FCONEL(LL)
!                    Write(*,*) IEL,VAL
                   ELSEIF(IND.EQ.19) THEN
C   zadata  prelaznost, stampa temperaturu okoline na elementu
                     IF(ISNUMBER.EQ.0) THEN
                       IEL=NELTOK(1,LL)
                     ELSE
                       IEL=MCVEL(NELTOK(1,LL))
                     ENDIF
!                    Write(*,*) IEL
                     VAL=TOKOLINEL(LL)
!                    Write(*,*) IEL,VAL
                   ENDIF
                  WRITE(II,5000) IEL,VAL
C                       IF(KONT1(IEL).EQ.0) THEN
C                          KONT1(IEL)=1
C                          GO TO 1202
C                       ELSE
C                          GO TO 10
C                       ENDIF 
C                    ENDIF
C                 ENDDO
!                  ENDDO
                 ENDDO
                 FSP(J)=0.D0
 1202            CONTINUE
               ENDIF
            ELSE
!              FSP(J)=RTH(J,I)
            ENDIF
   20    CONTINUE
!          VAL=FSP(1)
!          IF(IP.EQ.1.AND.IN(1).GT.40)
!      +      VAL=DSQRT(FSP(1)*FSP(1)+FSP(2)*FSP(2)+FSP(3)*FSP(3))
!          IF(IP.EQ.2) VAL=FSP(1)
!          IF(IP.EQ.3) VAL=FSP(2)
!          IF(IP.EQ.4) VAL=FSP(3)
! C         IF(DABS(VAL).LT.1.E-12) GO TO 10
!          IF(IN(1).EQ.10) THEN
!             WRITE(II,5000) IEL,VAL
! C         ELSE
! C            WRITE(II,5000) I,VAL
!          ENDIF  
C   10 CONTINUE
      WRITE(II,5000) MJ,ZERO
      IF(IP.LT.4.AND.IN(1).GT.40) GO TO 30
      RETURN
 1000 FORMAT(10I5)
 1001 FORMAT(I10,F12.4)
 3001 FORMAT('TEMPERATURE')
 3002 FORMAT('TEMPERATURE-BOFANG CONVECTION')
!  3003 FORMAT('PRESCRIBED SURFACES')
 3003 FORMAT('PRESCRIBED RADIATION-FLUX')
 3012 FORMAT('PRESCRIBED CONVECTION')
 3041 FORMAT('TOTAL VELOCITY')
 3042 FORMAT('VX VELOCITY')
 3043 FORMAT('VY VELOCITY')
 3044 FORMAT('VZ VELOCITY')
 3051 FORMAT('TOTAL GRADIENT')
 3052 FORMAT('GX GRADIENT')
 3053 FORMAT('GY GRADIENT')
 3054 FORMAT('GZ GRADIENT')
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(3(1PD12.4))
 5000 FORMAT(I10,3(1PD12.4))
      END
C=====================================================================
C=======================================================================
      SUBROUTINE VSTAU09T(RTH,NPP,II,IND,IKOL,KOR,NZAD,NUMZAD,
     +                  KONT,POT,ISNUMBER)
      USE NODES
      USE ELEMENTS
      USE KONTURE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
      COMMON /VOPTIM/ NKONT,MAXLIN
      DIMENSION RTH(IKOL,*),FSP(6),IN(8),NZAD(3,*)
      DIMENSION KONT(9,MAXLIN,*),POT(*)
!       DIMENSION KONT1(40000)
C
      FSP(3)=0.
C
      JEDAN=1
      NULA=0
      ZERO=0.
      ONE=1.
      MJ=-1
      M451=451
      M1=1
      M2=2
      M3=3
      M4=4
      M5=5
      M6=6
      M7=7
      M8=8
      M9=9
C      II=49
C
      IN(1)=IND
      IN(2)=IND+1
      IN(3)=IND+2
      IN(4)=IND+3
      IN(5)=IND+4
      IN(6)=IND+5
      IN(7)=IND+6
      IN(8)=IND+7
C
C PRIVREMENO
      IF(IN(1).EQ.10) THEN
!       DO I=1,40000
!          KONT1(I)=0
!       ENDDO
      ENDIF
C      
      IP=0
C
   30 IP=IP+1
      WRITE(II,1000) KOR,IN(IP),JEDAN
      IF(IN(1).EQ.1.AND.IP.EQ.1) WRITE(II,3001) 
      IF(IN(1).EQ.9.AND.IP.EQ.1) WRITE(II,3002) 
!       IF(IN(1).EQ.10.AND.IP.EQ.1) WRITE(II,3003) 
      IF(IN(1).EQ.51.AND.IP.EQ.1) WRITE(II,3051) 
      IF(IN(1).EQ.51.AND.IP.EQ.2) WRITE(II,3052) 
      IF(IN(1).EQ.51.AND.IP.EQ.3) WRITE(II,3053) 
      IF(IN(1).EQ.51.AND.IP.EQ.4) WRITE(II,3054) 
      WRITE(II,1010) ZERO,ZERO,ZERO
C OVO SU NULE ZA TEMPERATURE
      IF(IP.EQ.1.AND.IN(1).LT.40)
     +WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.1.AND.IN(1).GT.40)
     +WRITE(II,1000)IN(2),IN(3),IN(4),NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.2)
     +  WRITE(II,1000)IN(2),NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.3)
     +  WRITE(II,1000)NULA,IN(3),NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.4)
     +  WRITE(II,1000)NULA,NULA,IN(4),NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
C
C 0=ANY,1=DISPL,2=ACCEL,3=FORCE,4=STRESS,5=STRAIN,6=TEMP,OTHER=USER
C     NODAL(7) ,ELEMENTAL(8)
C dodat sledeci red za verziju neu fajla 10.0 /Sneza 26.12.2016.
      WRITE(II,1000) NULA
      IF(IN(1).EQ.1) WRITE(II,1000) NULA,NULA,M6,M7
      IF(IN(1).EQ.9) WRITE(II,1000) NULA,NULA,M7,M7
      IF(IN(1).EQ.10) WRITE(II,1000) NULA,NULA,M8,M8
      IF(IN(1).EQ.11) WRITE(II,1000) NULA,NULA,M6,M7
      IF(IN(1).EQ.41.OR.IN(1).EQ.51.OR.IN(1).EQ.61)
     1WRITE(II,1000) NULA,NULA,M3,M7
C
C =1, CAN NOT LINEARY COMBINE THIS OUTPUT
C =1, COMP(0-2) ARE THE X,Y,Z COMPONENT
C =1, THIS VECTOR HAS CONTROIDAL OR NODAL OUTPUT
      IF(IN(1).GT.40) THEN
         IF(IP.EQ.1) THEN
            WRITE(II,1000) JEDAN,JEDAN,JEDAN
         ELSE
            WRITE(II,1000) NULA,JEDAN,JEDAN
         ENDIF
      ELSE
         WRITE(II,1000) JEDAN,NULA,JEDAN
      ENDIF
      DO 10 I=1,NPP
         DO 20 J=1,IKOL
            FSP(J) = 0.0D0
            IF(IKOL.EQ.1) THEN
               IF(IN(1).EQ.1) THEN
                  IF(ID(J,I).EQ.0) GO TO 20
                  K = ID(J,I)
                  IF(K.GT.0) THEN
!                      IF(INDFS.GT.0) THEN
c                       POTENCIJAL ISPOD SLOBODNE POVRSINE
!                         IF(POT(K).GT.CORD(IOSA,I)) then
                            FSP(J)=RTH(J,K)
!                         else 
!                             go to 10
!                         endif
!                      ELSE   
!                         IF(POT(K).GT.CORD(IOSA,I)) then
!                             FSP(J)=RTH(J,K)
!                         else 
!                             go to 10
!                         endif
c  Sneza 10.9.2015  'stampa vrednost 0 iznad slobodne povrsine
C                       FSP(J)=RTH(J,K) 
!                      ENDIF
                  ENDIF
!                ELSEIF(IN(1).EQ.11) THEN
!                   IF(ID(J,I).EQ.0) GO TO 20
!                   K = ID(J,I)
!                   IF(K.GT.0) THEN
! C                    DUBINA VODE ILI PRITISAK VODENOG STUBA
!                      FSP(J)=-(RTH(1,K)-CORD(IOSA,I))
!                      IF(FSP(J).GT.0.) FSP(J)=0.     
!                   ENDIF
               ELSEIF(IN(1).EQ.9) THEN
                 FSP(J)=1.0D0
                 DO L=1,NUMZAD 
                    IF(NZAD(1,L).EQ.I) GO TO 1201
                 ENDDO
                 FSP(J)=0.D0
 1201            CONTINUE
               ELSE
                 DO LL=1,NKONT
                    FSP(J)=LL
                 DO L=1,LIN(LL)
                 DO LI=2,9 
C                 DO LI=2,NNODE+1 
                    IF(KONT(LI,L,LL).EQ.I) THEN
                     IF(ISNUMBER.EQ.0) THEN
                       IEL=KONT(1,L,LL)
                     ELSE
                       IEL=MCVEL(KONT(1,L,LL))
                     ENDIF
!                        IF(KONT1(IEL).EQ.0) THEN
!                           KONT1(IEL)=1
!                           GO TO 1202
!                        ELSE
!                           GO TO 10
!                        ENDIF 
                    ENDIF
                 ENDDO
                 ENDDO
                 ENDDO
                 FSP(J)=0.D0
 1202            CONTINUE
               ENDIF
            ELSE
               IF(INDFS.GT.0) THEN
                  K = ID(1,I)
                  IF(K.GT.0) THEN
!                      IF(POT(K).GT.CORD(IOSA,I)) FSP(J)=RTH(J,I)
c                  ELSE
!                      FSP(J)=RTH(J,I)
                  ENDIF
               ELSE   
                  FSP(J)=RTH(J,I)
               ENDIF
            ENDIF
   20    CONTINUE
         VAL=FSP(1)
         IF(IP.EQ.1.AND.IN(1).GT.40)
     +      VAL=DSQRT(FSP(1)*FSP(1)+FSP(2)*FSP(2)+FSP(3)*FSP(3))
         IF(IP.EQ.2) VAL=FSP(1)
         IF(IP.EQ.3) VAL=FSP(2)
         IF(IP.EQ.4) VAL=FSP(3)
czn         IF(DABS(VAL).LT.1.E-12) GO TO 10
         IF(IN(1).EQ.10) THEN
            WRITE(II,5000) IEL,VAL
         ELSE
            IF(ISNUMBER.EQ.0) THEN
              WRITE(II,5000) I,VAL
            ELSE
              WRITE(II,5000) NCVEL(I),VAL
            ENDIF
         ENDIF  
   10 CONTINUE
      WRITE(II,5000) MJ,ZERO
      IF(IP.LT.4.AND.IN(1).GT.40) GO TO 30
      RETURN
 1000 FORMAT(10I5)
 1001 FORMAT(I10,F12.4)
 3001 FORMAT('TEMPERATURE')
 3002 FORMAT('PRESCRIBED TEMPERATURE')
!  3003 FORMAT('PRESCRIBED SURFACES')
 3051 FORMAT('TOTAL GRADIENT')
 3052 FORMAT('GX GRADIENT')
 3053 FORMAT('GY GRADIENT')
 3054 FORMAT('GZ GRADIENT')
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(3(1PD12.4))
 5000 FORMAT(I10,3(1PE12.4))
      END
C=======================================================================
       SUBROUTINE VSTAMPPRESEKNEU(TT1,VECTJ,
     +       VREME,VVREME,KKORAK,ISNUMER)
       USE ppr
       USE PRESEK
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT,GAMA
       COMMON /VREPER/ NPER,NTABFT
       COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
       COMMON /DJERDAP/ IDJERDAP,ISPRESEK
       DIMENSION TT1(*),VECTJ(3,*)
       DIMENSION VREME(NPER,950)

       DO IPR=1,IPRES
       IFILE=300
          IFILE=IFILE+IPR
!           write(*,*) "stampanje potencijala preseka"
       CALL VSTAU09PRESEK(TT1,IFILE,1,1,KKORAK,TT1,ISNUMER,10,IPR)
        if(IPAKT.EQ.0) THEN
!           write(*,*) "stampanje dubine vode"
      CALL VSTAU09PRESEK(TT1,IFILE,11,1,KKORAK,TT1,ISNUMER,10,IPR)
!           write(*,*) "stampanje brzina preseka"
        CALL VSTAU09PRESEK(VECTJ,IFILE,41,3,KKORAK,TT1,
     1             ISNUMER,10,IPR)
!           write(*,*) "stampanje pornih pritisaka preseka"
         CALL VSTAU09CTPP(TT1,IFILE,VVREME,KKORAK,ISNUMER,10,IPR)
         ENDIF
        ENDDO
!        
       RETURN  
       END
C=======================================================================
C Stampanje POTENCIJALA i BRZINE U Presecima ZA CRTANJE REZULTATA
C=======================================================================
      SUBROUTINE VSTAU09PRESEK(RTH,II,IND,IKOL,KOR,
     +                  POT,ISNUMBER,IVER,IPR)
      USE NODES
      USE ELEMENTS
      USE PRESEK
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      DIMENSION RTH(IKOL,*),FSP(6),IN(8),H(10)
      DIMENSION POT(*)
C
      FSP(3)=0.
C
      JEDAN=1
      NULA=0
      ZERO=0.
      ONE=1.
      MJ=-1
      M451=451
      M1=1
      M2=2
      M3=3
      M4=4
      M5=5
      M6=6
      M7=7
      M8=8
      M9=9
C      II=49
!       write(*,*) "ifile", II,IND,IKOL,KOR,NPRESEK
C
      IN(1)=IND
      IN(2)=IND+1
      IN(3)=IND+2
      IN(4)=IND+3
      IN(5)=IND+4
      IN(6)=IND+5
      IN(7)=IND+6
      IN(8)=IND+7
C
      IP=0
C
   30 IP=IP+1
      WRITE(II,1000) KOR,IN(IP),JEDAN
       if(IPAKT.EQ.0) THEN
      IF(IN(1).EQ.1.AND.IP.EQ.1) WRITE(II,3001) 
      IF(IN(1).EQ.11.AND.IP.EQ.1) WRITE(II,3004) 
      ELSE
      IF(IN(1).EQ.1.AND.IP.EQ.1) WRITE(II,3002) 
      ENDIF
      IF(IN(1).EQ.41.AND.IP.EQ.1) WRITE(II,3041) 
      IF(IN(1).EQ.41.AND.IP.EQ.2) WRITE(II,3042) 
      IF(IN(1).EQ.41.AND.IP.EQ.3) WRITE(II,3043) 
      IF(IN(1).EQ.41.AND.IP.EQ.4) WRITE(II,3044) 
      WRITE(II,1010) ZERO,ZERO,ZERO
C OVO SU NULE ZA TEMPERATURE
      IF(IP.EQ.1.AND.IN(1).LT.40)
     +WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.1.AND.IN(1).GT.40)
     +WRITE(II,1000)IN(2),IN(3),IN(4),NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.2)
     +  WRITE(II,1000)IN(2),NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.3)
     +  WRITE(II,1000)NULA,IN(3),NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      IF(IP.EQ.4)
     +  WRITE(II,1000)NULA,NULA,IN(4),NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
C
C 0=ANY,1=DISPL,2=ACCEL,3=FORCE,4=STRESS,5=STRAIN,6=TEMP,OTHER=USER
C     NODAL(7) ,ELEMENTAL(8)
C dodat sledeci red za verziju neu fajla 10.0 /Sneza 26.12.2016.
      IF(IVER.EQ.10)  WRITE(II,1000) NULA
      IF(IN(1).EQ.1) WRITE(II,1000) NULA,NULA,M6,M7
      IF(IN(1).EQ.9) WRITE(II,1000) NULA,NULA,M7,M7
      IF(IN(1).EQ.10) WRITE(II,1000) NULA,NULA,M8,M8
      IF(IN(1).EQ.11) WRITE(II,1000) NULA,NULA,M6,M7
      IF(IN(1).EQ.19) WRITE(II,1000) NULA,NULA,M7,M7
      IF(IN(1).EQ.41.OR.IN(1).EQ.51.OR.IN(1).EQ.61)
     1WRITE(II,1000) NULA,NULA,M3,M7
C
C =1, CAN NOT LINEARY COMBINE THIS OUTPUT
C =1, COMP(0-2) ARE THE X,Y,Z COMPONENT
C =1, THIS VECTOR HAS CONTROIDAL OR NODAL OUTPUT
      IF(IN(1).GT.40) THEN
         IF(IP.EQ.1) THEN
            WRITE(II,1000) JEDAN,JEDAN,JEDAN
         ELSE
            WRITE(II,1000) NULA,JEDAN,JEDAN
         ENDIF
      ELSE
         WRITE(II,1000) JEDAN,NULA,JEDAN
      ENDIF
C ODREDJIVANJE POTENCIJALA/BRZINE u tackama preseka
      DO 10 I=1,NPRESEK(IPR)
!       write(*,*) "I",I
         DO 20 J=1,IKOL
            FSP(J) = 0.0D0
C TETRA ELEMENTI SA 10 CVOROVA
            RMIN=NP_COORDS(4,I,IPR)
            SMIN=NP_COORDS(5,I,IPR)
            TMIN=NP_COORDS(6,I,IPR)
            RST=1.0-RMIN-SMIN-TMIN
!        write(*,*) "RMIN",RMIN,SMIN,TMIN
            DO IH=1,10
              H(IH)=0.
            ENDDO
            H(5)=4*RMIN*RST
            H(6)=4*RMIN*SMIN
            H(7)=4*SMIN*RST
            H(8)=4*TMIN*RST
            H(9)=4*RMIN*TMIN
            H(10)=4*SMIN*TMIN
            H(1)=RST-0.5*(H(5)+H(7)+H(8))
            H(2)=RMIN-0.5*(H(5)+H(6)+H(9))
            H(3)=SMIN-0.5*(H(6)+H(7)+H(10))
            H(4)=TMIN-0.5*(H(8)+H(9)+H(10))
            PPC=0.0D0
            VPC=0.0D0
            XPC=0.0D0
            YPC=0.0D0
            ZPC=0.0D0
            IF(IKOL.EQ.1) THEN
               IF(IN(1).EQ.1) THEN
!             write(3,*)"I,NP_ELEMENT(I)", I,NP_ELEMENT(I,IPR)       
!                write(3,*) RMIN,SMIN,TMIN      
                 DO JJH=1,10
                   ncvor=NEL(JJH,NP_ELEMENT(I,IPR))  
                   PPC=PPC + H(JJH)*RTH(1,ncvor)
                   XPC=XPC+H(JJH)*CORD(1,ncvor)
                   YPC=YPC+H(JJH)*CORD(2,ncvor)
                   ZPC=ZPC+H(JJH)*CORD(3,ncvor)
!             write(3,*)"ncvor", JJH,ncvor,H(JJH),RTH(1,ncvor)  
!             write(3,*) PPC
                 enddo
C    provera da li je cvor preseka okvasen
                IF((IPAKT.EQ.0).AND.((PPC-0.0001)
     1             .LT.NP_COORDS(IOSA,I,IPR))) goto 10
                 FSP(J)=PPC
!                  write(3,*) FSP(J),XPC,YPC,ZPC
              ELSEIF(IN(1).EQ.11) THEN
                 DO JJH=1,10
                   ncvor=NEL(JJH,NP_ELEMENT(I,IPR))  
                   PPC=PPC + H(JJH)*RTH(1,ncvor)
                   XPC=XPC+H(JJH)*CORD(1,ncvor)
                   YPC=YPC+H(JJH)*CORD(2,ncvor)
                   ZPC=ZPC+H(JJH)*CORD(3,ncvor)
                 enddo
C    provera da li je cvor preseka okvasen
                IF((IPAKT.EQ.0).AND.((PPC-0.0001)
     1             .LT.NP_COORDS(IOSA,I,IPR))) goto 10
                 FSP(J)=PPC
!                      DUBINA VODE ILI PRITISAK VODENOG STUBA
                 FSP(J)=(PPC-NP_COORDS(IOSA,I,IPR))
               ENDIF
            ELSE
!  stampa brzine samo za okvasene cvorove 19.5.2018            
!               write(*,*)"J,NP_ELEMENT(I)", J,NP_ELEMENT(I)       
                 DO JJH=1,10
                   ncvor=NEL(JJH,NP_ELEMENT(I,IPR))          
!               write(*,*)"JJH,ncvor", JJH,ncvor       
                   PPC=PPC + H(JJH)*POT(ncvor)
                   VPC=VPC + H(JJH)*RTH(J,ncvor)
                   XPC=XPC+H(JJH)*CORD(1,ncvor)
                   YPC=YPC+H(JJH)*CORD(2,ncvor)
                   ZPC=ZPC+H(JJH)*CORD(3,ncvor)
                 enddo
                 FSP(J)=VPC
            IF((PPC-0.0001).LT.NP_COORDS(IOSA,I,IPR)) goto 20
!               IF((PPC-0.0001).LT.NP_COORDS(IOSA,I)) FSP(J)=0.0
!           IF((POT(K)-0.001).GT.CORD(IOSA,I)) FSP(J)=RTH(J,K)
             ENDIF
   20    CONTINUE
         VAL=FSP(1)
         IF(IP.EQ.1.AND.IN(1).GT.40)
     +      VAL=DSQRT(FSP(1)*FSP(1)+FSP(2)*FSP(2)+FSP(3)*FSP(3))
         IF(IP.EQ.2) VAL=FSP(1)
         IF(IP.EQ.3) VAL=FSP(2)
         IF(IP.EQ.4) VAL=FSP(3)
czn         IF(DABS(VAL).LT.1.E-12) GO TO 10
!          IF(IN(1).EQ.10) THEN
!            IF(ISNUMBER.EQ.0) THEN
!             WRITE(II,5000) IEL,VAL
!            ELSE
!             WRITE(II,5000) MCVEL(IEL),VAL
!            ENDIF
!          ELSE
!  u presecima je uvek numeracija od 1 do NPRESEK(IPR)
!            IF(ISNUMBER.EQ.0) THEN
            WRITE(II,5000) I,VAL
!            ELSE
! !              write(*,*) "I,NCVEL(I)",I,NCVEL(I)
!             WRITE(II,5000) NCVEL(I),VAL
!            ENDIF
!          ENDIF  
   10 CONTINUE
      WRITE(II,5000) MJ,ZERO
      IF(IP.LT.4.AND.IN(1).GT.40) GO TO 30
      RETURN
 1000 FORMAT(10I5)
 1001 FORMAT(I10,F12.4)
 3001 FORMAT('TOTAL POTENTIAL')
 3002 FORMAT('TEMPERATURE')
 3004 FORMAT('DEPTH OF WATER')
 3041 FORMAT('TOTAL VELOCITY')
 3042 FORMAT('VX VELOCITY')
 3043 FORMAT('VY VELOCITY')
 3044 FORMAT('VZ VELOCITY')
 1005 FORMAT('PAK CASE',I5)
 1006 FORMAT('CASE',I5,' TIME',1PD12.4)
 1010 FORMAT(3(1PD12.4))
 5000 FORMAT(I10,3(1PD12.4))
      END
C=======================================================================
C Stampanje pornih pritisaka u elementima preseka
C=======================================================================
      SUBROUTINE VSTAU09CTPP(RTH,II,VREME,KOR,ISNUMBER,IVER,IPR)
      USE NODES
      USE ELEMENTS
      USE PRESEK
      USE ppr
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
!       COMMON /DODAT/ NDIMM
C      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
!       COMMON /VOPTIM/ NKONT,MAXLIN,LIN(10)
      DIMENSION RTH(1,*)
      DIMENSION H(10)
!       COMMON /PIJEZO/ CPOR(3,1000),NPIJEZ(20,100),NODP(100),NPIJ,NPOR,
!      1                NPORCV(1000),NPOREL(1000),NEPOR,NPORE,NEPIJ
!       COMMON /PORCEL/ RMIN(34),SMIN(34),TMIN(34)
C
!       IF(NEPIJ.EQ.0) GO TO 100
C
      MJ=-1
C       
C porne celije
      JEDAN=1
      NULA=0
      ZERO=0.
      ONE=1.
      MJ=-1
      M451=451
      M1=1
      M2=2
      M3=3
      M4=4
      M5=5
      M6=6
      M7=7
      M8=8
      M9=9
C      II=49
!       write(*,*) "ifile", II,IND,IKOL,KOR,NPRESEK
C
      IND=60
C
      WRITE(II,1000) KOR,IND,JEDAN
      WRITE(II,3001) 
      WRITE(II,1010) ZERO,ZERO,ZERO
C OVO SU NULE ZA TEMPERATURE
      WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
      WRITE(II,1000) NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA,NULA
C
C 0=ANY,1=DISPL,2=ACCEL,3=FORCE,4=STRESS,5=STRAIN,6=TEMP,OTHER=USER
C     NODAL(7) ,ELEMENTAL(8)
C dodat sledeci red za verziju neu fajla 10.0 /Sneza 26.12.2016.
      IF(IVER.EQ.10)  WRITE(II,1000) NULA
      WRITE(II,1000) NULA,NULA,M6,M7
C
C =1, CAN NOT LINEARY COMBINE THIS OUTPUT
C =1, COMP(0-2) ARE THE X,Y,Z COMPONENTNUMMAT
C =1, THIS VECTOR HAS CONTROIDAL OR NODAL OUTPUT
         WRITE(II,1000) JEDAN,NULA,JEDAN
! 
         DO 14 I=1,NPRESEK(IPR)
C ODREDJIVANJE pornih pritisaka u cvorovima preseka
C za heksa elemente
C              RPP=1.0+RMIN(I)
C              SPP=1.0+SMIN(I)
C              TPP=1.0+TMIN(I)
C              RM=1.0-RMIN(I)
C              TM=1.0-TMIN(I)
C              SM=1.0-SMIN(I)
C              RR=1.0-RMIN(I)*RMIN(I)
C              SS=1.0-SMIN(I)*SMIN(I)
C              TT=1.0-TMIN(I)*TMIN(I)
C              DO IH=1,8
C                H(IH)=0.
C              ENDDO
C              H(1)=0.125*RPP*SPP*TPP
C              H(2)=0.125*RM*SPP*TPP
C              H(3)=0.125*RM*SM*TPP
C              H(4)=0.125*RPP*SM*TPP
C              H(5)=0.125*RPP*SPP*TM
C              H(6)=0.125*RM*SPP*TM
C              H(7)=0.125*RM*SM*TM
C              H(8)=0.125*RPP*SM*TM
C TETRA ELEMENTI SA 10 CVOROVA
              RMIN=NP_COORDS(4,I,IPR)
              SMIN=NP_COORDS(5,I,IPR)
              TMIN=NP_COORDS(6,I,IPR)
              RST=1.0-RMIN-SMIN-TMIN
              DO IH=1,10
                H(IH)=0.
              ENDDO
              H(5)=4*RMIN*RST
              H(6)=4*RMIN*SMIN
              H(7)=4*SMIN*RST
              H(8)=4*TMIN*RST
              H(9)=4*RMIN*TMIN
              H(10)=4*SMIN*TMIN
              H(1)=RST-0.5*(H(5)+H(7)+H(8))
              H(2)=RMIN-0.5*(H(5)+H(6)+H(9))
              H(3)=SMIN-0.5*(H(6)+H(7)+H(10))
              H(4)=TMIN-0.5*(H(8)+H(9)+H(10))
            PPC=0.
            XPC=0.
            YPC=0.
            ZPC=0.
        NGAUS=0    
        dgaust=0.0
        minr=50.0
        ngnaj=1
        IBRGT1=4
        DO L=1,IBRGT1
          NGAUS=NGAUS+1
          CALL VGETRST(RT,ST,TT,WTET,NGAUS,IBRGT1)
        dgaust=(RMIN-RT)*(RMIN-RT)+(SMIN-ST)*(SMIN-ST)
         dgaust=dgaust+(TMIN-TT)*(TMIN-TT)
         dgaus=sqrt(dgaus)
          IF (dgaus.LT.minr) then
            minr=dgaus
            ngnaj=NGAUS
          endif
        ENDDO
!             DO JJH=1,10
!               ncvor=NEL(JJH,MP_ELEMENT(I))          
!               PPC=PPC + H(JJH)*RTH(1,ncvor)
!               XPC=XPC+H(JJH)*CORD(1,ncvor)
!               YPC=YPC+H(JJH)*CORD(2,ncvor)
!               ZPC=ZPC+H(JJH)*CORD(3,ncvor)
!            enddo
C
!             MP_RESULTS(KORAK,I+1)=PPC
        IF(PORNIEL(NP_ELEMENT(I,IPR),ngnaj).GT.0.D0)
     1   WRITE(II,5010) NP_ID(I,IPR),PORNIEL(NP_ELEMENT(I,IPR),ngnaj)
   14    CONTINUE
         WRITE(II,5000) MJ,ZERO
!          MP_VREME(KORAK)=VREME
!          write(3,*) 'KORAK',KORAK
!           write(3,*) 'MP_RESULTS(KORAK,2)',MP_RESULTS(KORAK,2)
!          if(korak.eq.BRKORAKA) then
!             do ikbr=1,BRKORAKA
!               write(3,*) 'MP_RESULTS(',ikbr,',2)',MP_RESULTS(ikbr,2)
!             enddo
!          endif
C
      RETURN
 1000 FORMAT(10I6)
 1001 FORMAT(I10,F12.4)
 1005 FORMAT('PAK CASE',I5)
 1010 FORMAT(3(1PD12.4))
 3001 FORMAT('PORE PRESSURE')
 5000 FORMAT(I10,3(1PE12.4))
 5010 FORMAT(I10,4(1PE12.4))
 5001 FORMAT(3(I10),2(1PE12.4))
      END
C=======================================================================
