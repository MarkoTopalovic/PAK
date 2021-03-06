SUBROUTINE mcm_geometry_spvp
!**********************************************************************
!
!    Purpose : This module writes out the .geo file
!
!  Called by : init_state_output
!
!       Date : 22-07-99
!
!     Errors : 
!
!      Notes :  The .geo file contain the coordinates of all nodes
!               and the repartition of differents parts ( that is to
!               say the different pieces )
!
!**********************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
LOGICAL fexist
INTEGER :: mat,ind=1,i,j
REAL(kind=real_acc) :: xnode,ynode,znode,zero=0
CHARACTER(LEN=25) :: file_write
INTEGER, DIMENSION(30) :: spmaterial, vpmaterial
!
do i=1,30
 spmaterial(i) = 0
 vpmaterial(i) = 0
enddo
!
!
! if old time file exists then delete it
!
mcm_file_time=mcm_fileout(1:mcm_filelen(2))//"_time.txt"
inquire(file=mcm_file_time,exist=fexist)
if(fexist) then
	open(unit=3,file=mcm_file_time,status='old',form='formatted')
	close(unit=3,status='delete')
end if
!
file_write=mcm_fileout(1:mcm_filelen(2))//".geo"
open(unit=2,file=file_write,status='unknown',form='formatted')
	!
	! The following lines must allways be exactly the same
	!
	write(2,*)"frame geometry"
	write(2,*)"node id assign"
	write(2,*)"element id assign"
	write(2,*)"coordinates"
	write(2,fmt='(I7)') mcm_np  !write number of nodes
	!
	!
	! Next loop writes the coordinates of all nodes
	!
	select case(mcm_ndim)
		case(2)
			do j=1,mcm_np
			    mat=real(par(j)%mat)
				write(2,fmt='(E12.5,E12.5,E12.5)') par(j)%x(1),par(j)%x(2),0.0_d
				if(j.lt.mcm_svp) then
				 spmaterial(mat)=spmaterial(mat)+1
				else
				 vpmaterial(mat)=vpmaterial(mat)+1
				endif
			end do
		case(3)
			do j=1,mcm_np
				mat=real(par(j)%mat)
				write(2,fmt='(E12.5,E12.5,E12.5)') par(j)%x(1),par(j)%x(2),par(j)%x(3)
				vpmaterial(mat)=vpmaterial(mat)+1
			end do
	end select
	!
	! Write the repartition of the nodes into the different parts
	!
	do j=1,mcm_nummat
		write(2,fmt='(A5,I1)') "part ",j
		write(2,fmt='(A20,I1)')"points list of part ",j
		write(2,*)"point"
		write(2,fmt='(I8)') vpmaterial(j)
		do i=mcm_svp,mcm_evp
		 if(par(i)%mat.eq.j) then
			write(2,fmt='(I8)') i
		 endif
		end do
	end do
    
	do j=1,mcm_nummat
		write(2,fmt='(A5,I1)') "part ",mcm_nummat+j
		write(2,fmt='(A20,I1)')"points list of part ",mcm_nummat+j
		write(2,*)"point"
		write(2,fmt='(I8)') spmaterial(j)
		do i=mcm_ssp,mcm_esp
		 if(par(i)%mat.eq.j) then
			write(2,fmt='(I8)') i
		 endif
		end do
	end do

	close(unit=2)

END SUBROUTINE mcm_geometry_spvp
