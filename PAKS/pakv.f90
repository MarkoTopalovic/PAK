SUBROUTINE VTKPRINT(A)

COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,IOPGL(6),KOSI,NDIN,ITEST
 COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,LMXAU,LAPRS
 
DIMENSION A(*)
CALL PAKSVTK(A(LCORD),A(LCVEL),ICVEL,NP,88)


END





subroutine PAKSVTK(CORD,NCVEL,ICVEL,NPP,II)
!************************************************************************
!
!    Purpose: vrite VTK files for Paraview 
!
!  Called by: 
!
!       Date: 0
!
!     Errors: 
!
!      Notes: 
!
!************************************************************************
!
      USE mcm_database
      USE CVOROVI
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        character(len=9) extension
        character(len=65) paks_filename_vtk
        character(len=3) prefix
        integer  i
        character(len=7) charParNumb,charParNumbDouble
        integer mcm_np1
        integer ICVEL, NI, NCVEL, J, NNP, II, IDEAS, IDD, NPP

        CHARACTER*20    VTKIME
        INTEGER IVTKCOUNTER
        COMMON /VTKVALUES/ VTKIME,IVTKCOUNTER
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,IOPGL(6),KOSI,NDIN,ITEST
      COMMON /NIDEAS/ IDEAS
      DIMENSION CORD(NPP,*),NCVEL(*)
      DIMENSION IDD(6)
!      ALLOCATE (VTKCVOROVI(NP,3))
      IF (ALLOCIRANICVOROVI.EQ.(.FALSE.)) THEN
      ALLOCATE (VTKCVOROVI(NP,3))
      ALLOCIRANICVOROVI = .TRUE.
      END IF

      NULA=0
	 JEDAN=1

prefix = 'FEM'
!
if (IVTKCOUNTER.lt.10) then
	write(unit=extension,fmt=1000) IVTKCOUNTER
else if (IVTKCOUNTER.lt.100) then
	write(unit=extension,fmt=1010) IVTKCOUNTER
else if (IVTKCOUNTER.lt.1000) then
	write(unit=extension,fmt=1020) IVTKCOUNTER	

else if (IVTKCOUNTER.lt.10000) then
	write(unit=extension,fmt=1030) IVTKCOUNTER
else
	write(unit=extension,fmt=1040) IVTKCOUNTER
end if
1000 format(I1,'.vtk')
1010 format(I2,'.vtk')
1020 format(I3,'.vtk')
1030 format(I4,'.vtk')
1040 format(I5,'.vtk')

if (NP.lt.10) then
	write(unit=charParNumb,fmt=2000) NP
	write(unit=charParNumbDouble,fmt=2000) 2*NP
else if (NP.lt.100) then
		write(unit=charParNumb,fmt=2010) NP
	write(unit=charParNumbDouble,fmt=2010) 2*NP
else
		write(unit=charParNumb,fmt=2020) NP
	write(unit=charParNumbDouble,fmt=2020) 2*NP
end if
2000 format(I1)
2010 format(I2)
2020 format(I3)

 

!    paks_filename_vtk = prefix//IME(1:(ipaks_filelen-4))//extension
paks_filename_vtk = prefix//VTKIME(1:(ipaks_filelen))//extension
mcm_np1 = NP

   open(unit=88,file=paks_filename_vtk,status='unknown',form='formatted')

WRITE(88, 81)
WRITE(88, 82)
WRITE(88, 83)
WRITE(88, 84)

NPcvorova = NP

DO  I=1,NP
      DO J=1,3
      IF (ALLOCIRANAPOMERANJA.EQ.(.TRUE.)) THEN
        VTKCVOROVI(I,J) = CORD(I,J)+VTKPOMERANJA(I,J)
        ELSE
        VTKCVOROVI(I,J) = CORD(I,J)
        END IF
    END DO
END DO
WRITE(88,*)"POINTS ",charParNumb," double"
DO I=1,NP
         WRITE(88,882)(VTKCVOROVI(I,J),J=1,3)
ENDDO
IF (mcm_np.GT.NP) THEN ! ako nije znaci da nije jos ucitao mcm ulaznu datoteku
DO I=1,NP
        ! cuvanje starih koordinata pre azuriranja u paku
         par(mcm_np-NP+I)%xold(1) = par(mcm_np-NP+I)%x(1)
         par(mcm_np-NP+I)%xold(2) = par(mcm_np-NP+I)%x(2)
         par(mcm_np-NP+I)%xold(3) = par(mcm_np-NP+I)%x(3)
        ! prenos koordinata iz paka u mcm Topalovic
         par(mcm_np-NP+I)%xmove(1) = VTKCVOROVI(I,1)
         par(mcm_np-NP+I)%xmove(2) = VTKCVOROVI(I,2)
         par(mcm_np-NP+I)%xmove(3) = VTKCVOROVI(I,3)
        
ENDDO
END IF
!
!CLOSE(UNIT=88)
!!
  81 FORMAT ("# vtk DataFile Version 3.0")
  82 FORMAT ("TEST")
  83 FORMAT ("ASCII")
  84 FORMAT ("DATASET UNSTRUCTURED_GRID")

 882 FORMAT (F15.8,F15.8,F15.8)

end subroutine PAKSVTK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PAKSVTKELEM(NEL,MCVEL,ICVEL,NMAT,THID,IGRAF)
        USE ELEMENTI
        USE CVOROVI
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /SUMELE/ ISUMEL,ISUMGR
      COMMON /SRPSKI/ ISRPS
      COMMON /NIDEAS/ IDEAS
      COMMON /STAPRO/ THIDP(100),NPROPR(100)
      DIMENSION NEL(NE,*),MCVEL(*),NMAT(*),THID(*)
      DIMENSION FIZ(14)         
      COMMON /CDEBUG/ IDEBUG
      IF (ALLOCIRANIELEMENTI.EQ.(.FALSE.)) THEN
      ALLOCATE (VTKELEMENTI(NE,4))
      
      IF(IDEBUG.GT.0) PRINT *, ' TGRAU8'

        WRITE(88,*)"CELLS ",NE," ",(NE*5)
        DO  I=1,NE
            WRITE(88,1001) 4,(NEL(I,J)-1,J=1,4)
        END DO
        WRITE(88,*)"CELL_TYPES ",NE
        DO  I=1,NE
            WRITE(88,*) 9
        END DO
      CLOSE(UNIT=88)
              DO  I=1,NE
              DO J=1,4
            VTKELEMENTI(I,J) = NEL(I,J)
            END DO
              END DO
              
              ALLOCIRANIELEMENTI = .TRUE.
      END IF
      RETURN
 1001 FORMAT(5(I6))
 
      END
      
      SUBROUTINE VTKPRINTMODULEELEMENT
      USE ELEMENTI
      USE CVOROVI
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8

      IF(IDEBUG.GT.0) PRINT *, ' TGRAU8'
  
      IF (ALLOCIRANIELEMENTI.EQ.(.TRUE.)) THEN

        WRITE(88,*)"CELLS ",NE," ",(NE*5)
        DO  I=1,NE
            WRITE(88,1002) 4,(VTKELEMENTI(I,J)-1,J=1,4)
        END DO
        WRITE(88,*)"CELL_TYPES ",NE
        DO  I=1,NE
            WRITE(88,*) 9
        END DO
      CLOSE(UNIT=88)
      END IF        
      RETURN
 1002 FORMAT(5(I6))
      
      END