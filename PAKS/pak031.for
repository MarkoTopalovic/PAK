C=======================================================================
C
C=======================================================================
      SUBROUTINE UMOD51(FUNMAT,MAT,KARTI,MATE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO READ DATA FOR MATERIAL MODEL NUMBER 51
C .
CE.             MAT - MATERIAL NUMBER
C .
CE.   FUNMAT(MAT,1) - 
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      DIMENSION FUNMAT(MATE,*),IDUM(4)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' UMOD51'
C CARD 6B
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (IDUM(I),I=1,4)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1100) (IDUM(I),I=1,4)
      FUNMAT(MAT,30)=IDUM(1)
      FUNMAT(MAT,31)=IDUM(2)
      FUNMAT(MAT,32)=IDUM(3)
      FUNMAT(MAT,33)=IDUM(4)
C CARD 6C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(MAT,I),I=1,7)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(MAT,I),I=1,7)
C CARD 6D
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(MAT,I),I=8,13)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(MAT,I),I=8,13)
C CARD 6E
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(MAT,I),I=14,19)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(MAT,I),I=14,19)
C CARD 6F
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(MAT,I),I=20,22)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(MAT,I),I=20,22)
C CARD 6G
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(MAT,I),I=23,25)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(MAT,I),I=23,25)
C CARD 6H
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (FUNMAT(MAT,I),I=26,29)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (FUNMAT(MAT,I),I=26,29)
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) MAT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) MAT
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2012) (IDUM(I),I=1,4)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6012) (IDUM(I),I=1,4)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2013) (FUNMAT(MAT,I),I=1,7)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6013) (FUNMAT(MAT,I),I=1,7)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2014) (FUNMAT(MAT,I),I=8,13)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6014) (FUNMAT(MAT,I),I=8,13)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2015) (FUNMAT(MAT,I),I=14,19)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6015) (FUNMAT(MAT,I),I=14,19)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2016) (FUNMAT(MAT,I),I=20,22)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6016) (FUNMAT(MAT,I),I=20,22)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2017) (FUNMAT(MAT,I),I=23,25)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6017) (FUNMAT(MAT,I),I=23,25)
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2018) (FUNMAT(MAT,I),I=26,29)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6018) (FUNMAT(MAT,I),I=26,29)
      RETURN
C
 1100 FORMAT(14I5)
 1000 FORMAT(7F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(6X,
     1'MODEL MATERIJALA BROJ =    51 (SPREGNUTI PLASTICNI SA RAZARANJEM)
     2'///11X,'MATERIJAL BROJ =',I5)
 2012 FORMAT(//' CARD 6B'/4I5)
 2013 FORMAT(//' CARD 6C'/7(1PD11.4))
 2014 FORMAT(//' CARD 6D'/7(1PD11.4))
 2015 FORMAT(//' CARD 6E'/7(1PD11.4))
 2016 FORMAT(//' CARD 6F'/7(1PD11.4))
 2017 FORMAT(//' CARD 6G'/7(1PD11.4))
 2018 FORMAT(//' CARD 6H'/7(1PD11.4))
C-----------------------------------------------------------------------
 6000 FORMAT(6X,
     1'MATERIAL MODEL NUMBER =    51 (COUPLED PLASTIC-DAMAGED MODEL)'///
     211X,'MATERIAL NUMBER =',I5)
 6012 FORMAT(//' CARD 6B'/4I5)
 6013 FORMAT(//' CARD 6C'/7(1PD11.4))
 6014 FORMAT(//' CARD 6D'/7(1PD11.4))
 6015 FORMAT(//' CARD 6E'/7(1PD11.4))
 6016 FORMAT(//' CARD 6F'/7(1PD11.4))
 6017 FORMAT(//' CARD 6G'/7(1PD11.4))
 6018 FORMAT(//' CARD 6H'/7(1PD11.4))
C-----------------------------------------------------------------------
      END
