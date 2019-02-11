C======================================================================
      SUBROUTINE ISOHIP(X,Y,IPOINT,Y0,XX,YY,ZZ,AKK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /JEDANP/ INDJED,NBRF,NGL,INDTLO
      DIMENSION IPOINT(*),X(*),Y(*),Y0(*)

C      NBRF=44

	IF (INDJED.EQ.0) THEN
	  NGL=1
        CALL ISJEDN(X,Y,Y0,IPOINT,NBRF,NGL)
C        OPEN(23,FILE='ISOHYP')
C        WRITE(23,*)'        XX        YY        ZZ       AKK        Y0'
C        DO I=1,NBRF
C         WRITE(23,200) I,IPOINT(I)
C        ENDDO
        INDJED=1
      ENDIF

      MIN=0
      AMIN=1.D10
      DO I=1,NGL
       RXY=RASTXY(XX,YY,X(I),Y(I))
       IF(RXY.LT.AMIN) THEN
         MIN=I
         AMIN=RXY
       ENDIF
      ENDDO
             
	AKK=FITMAX(MIN,Y0,ZZ,NBRF,IPOINT,IPOZ)

C ANTILOGARITMOVANJE ZA PERMEABILNOST
       AKK=10.D0**(AKK)
C PERMEABILNOST U M/S:
       AKK=AKK/60.D0
C PONOVO VRACANJE NA LOGARITAM:
C       AKK=DLOG(AKK)
       
C      WRITE(23,100) XX,YY,ZZ,AKK,IPOZ,MIN
100   FORMAT(4(F10.3),2I5)
200   FORMAT(I5)
      END
C==========================================================================
C==========================================================================
      FUNCTION FITMAX(MIN,Y0,Z,NBRF,IPOINT,IP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y0(*),IPOINT(*),B(6)

C OVO SU KOEFICIJENTI POLINOMA ZA MAX DOBIJENI IZ MATLAB-A

       B(1)=-6.02781
       B(2)=44.47922
       B(3)=-45.70244
       B(4)=-116.93644
       B(5)=26.98260
       B(6)=64.56207
  

      NN=0
      DO I=1,NBRF
       IF ((MIN.GT.NN).AND.(MIN.LE.NN+IPOINT(I))) THEN
        IP=I
        GOTO 10
       ENDIF
       NN=NN+IPOINT(I)
      ENDDO

10     Y=Y0(IP)
       Y1=Y/1000.D0
       Z1=Z/1000.D0
       FITMAX=B(1)+B(2)*Y1+B(3)*Z1+B(4)*Y1**2+B(5)*Y1*Z1+B(6)*Z1**2          
        
      END
C==========================================================================
C==========================================================================
      FUNCTION RASTXY(X1,Y1,X2,Y2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       RASTXY=DSQRT((X2-X1)**2+(Y2-Y1)**2)
      END
C==========================================================================
C==========================================================================
      SUBROUTINE ISJEDN(X,Y,Y0,IPOINT,NBRF,NGL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      CHARACTER *9 IMEISO
      CHARACTER *5 IMEIS1
      CHARACTER *4 DODAT
	DIMENSION IPOINT(*),X(*),Y(*),Y0(*)
      
      IISO=23
      OPEN(22,FILE='ISOHIP.TXT')
C

C      OPEN(27,FILE='POLY.TXT')

C       DO I=1,6
C        READ(27,*)BB
C        B(I)=BB
C       ENDDO

C       CLOSE(27)
C

	CALL PRESYZ(Y0)
       DODAT='.TXT'
      DO J=1,NBRF
       READ(22,110) IMEIS1
       IMEISO=IMEIS1 // DODAT
       OPEN(IISO,FILE=IMEISO)
       READ (IISO,100) NP
        IPOINT(J)=NP
        DO I=NGL,NGL+NP-1
  	   READ(IISO,200) X(I),Y(I)
        ENDDO
        NGL=NGL+NP
       CLOSE(IISO)
       IMEISO='         '
      ENDDO
C
       NGL=NGL-1
       CLOSE(22)



 100   FORMAT(I5)
 110   FORMAT(A5)
 200   FORMAT(2F10.0)

      END
C======================================================================
C==========================================================================
      SUBROUTINE PRESYZ(Y0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

	DIMENSION Y0(*)

       Y0(1)=207.68430314
       Y0(2)=251.139
       Y0(3)=260.09571624
       Y0(4)=204.08796580
       Y0(5)=201.29538162
       Y0(6)=270.05510136
       Y0(7)=198.06351358
       Y0(8)=194.03668989
       Y0(9)=277.54294740
       Y0(10)=286.69464807
       Y0(11)=186.27250611
       Y0(12)=294.81009942
       Y0(13)=178.49513311
       Y0(14)=300.09637279
       Y0(15)=172.29609457
       Y0(16)=308.94501605
       Y0(17)=164.03280735
       Y0(18)=317.59078371
       Y0(19)=159.24523380
       Y0(20)=325.63220137
       Y0(21)=153.35418611
       Y0(22)=332.70951540
       Y0(23)=140.32480260
       Y0(24)=339.95542592
       Y0(25)=132.07422100
       Y0(26)=345.51255603
       Y0(27)=125.91739947
       Y0(28)=349.94827032
       Y0(29)=120.10795822
       Y0(30)=351.79694084
       Y0(31)=119.56305772
       Y0(32)=119.56305772
       Y0(33)=117.82974005
       Y0(34)=358.80682089
       Y0(35)=115.88649503
       Y0(36)=109.32575044
       Y0(37)=366.24660556 
       Y0(38)=379.357
       Y0(39)=92.14353962
       Y0(40)=395.02073547 
       Y0(41)=80.28820634
       Y0(42)=408.74017670 
       Y0(43)=66.217
       Y0(44)=-20.0D0


      END
C======================================================================

C======================================================================
      SUBROUTINE NGXYZ(NEL,NBREL,H,CORD,XX,YY,ZZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /DODAT/ NDIMM
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET

      DIMENSION NEL(NDIMM,*),CORD(3,*),H(*)
C
       XX=0.D0
       YY=0.D0
       ZZ=0.D0
       DO I=1,NDIM
        NODE=NEL(I,NBREL)
        XX=XX+H(I)*CORD(1,NODE)
        YY=YY+H(I)*CORD(2,NODE)
        ZZ=ZZ+H(I)*CORD(3,NODE)
       ENDDO

      END
C======================================================================
