subroutine mcm_state_output
!************************************************************************
!
!    Purpose: Control routine for writing state plots 
!
!  Called by: initial, solution
!
!       Date: 07-08-2002
!
!     Errors: If using ensight, stop if not 2D or 3D colocated.
!
!      Notes: Routine added to allow coice of writing MCMGUI files
!             or Ensight files.
!
!************************************************************************
!
use mcm_database
!
implicit none
!
select case (mcm_state_opt)
 case(1)
  !
  ! Ensight 5.0 output
  !
  select case (mcm_disctype)
   case(0,1)
    ! colocated
    select case (mcm_ndim)
     case(1)
      !
      write(*,1000)
      write(13,1000)
      stop
	  !  
     case(2)
	  !
      call mcm_write_ensight_2d
	  call mcm_write_res
	  !
	 case(3)
	  ! 
	  call mcm_write_ensight_3d
	  call mcm_write_res
	 end select
   case default
    !
    write(*,2000)
    write(13,2000)
    stop
  end select
  mcm_istate = mcm_istate + 1
  !

 case(2)
  !
  ! Ensight CASE file output
  !
  select case (mcm_disctype)
   case(0,1)
    ! colocated
    select case (mcm_ndim)
     case(1)
      !
      write(*,1000)
      write(13,1000)
      stop
	  !  
     case(2)
	  !
	  call mcm_write_case
      call mcm_write_case_2d
	  !
	 case(3)
	  !
	  call mcm_write_case 
	  call mcm_write_case_3d
	  ! Tijana added - write VTKs
	  call mcm_write_vtk	  
	 end select
   case default
    !
    write(*,2000)
    write(13,2000)
    stop
  end select
  !mcm_istate = mcm_istate + 1
  !
 case(3)
  !
  ! LSDYNA d3plot format
  !
  call mcm_d3plot_state
  !
 case(4)
  !
  ! Topalovic VTK
  !
  if ((mcm_kojpak.eq.5).and.(mcm_initialized.eq.(.false.)).and.((mcm_ptime.ge.mcm_endtime).or.(mcm_istate.eq.0))) then
  ! todo topalovic
  ! privremeno stavljeno false da stampa cesto vtkove, kada je true stampa samo 1 na kraju mcm proracuna u pak-mcm povezivanju
      mcm_initialized = .true.
      call WriteVtk
      call WriteVtkmat(1)
      !call WriteVtkmat(2)
      mcm_istate = mcm_istate + 1
  else if (mcm_kojpak.eq.4) then
      call WriteVtk
      call WriteVtkmat(1)
      !call WriteVtkmat(2)
      mcm_istate = mcm_istate + 1
  end if
  
 case default
  ! write MCMGUI output files
  call mcm_stateplot
  !

  
end select
!
return
!
1000 format('Ensight output not supported for 1D calculations')
2000 format('Ensight output not supported for stress-velocity points')
!
end subroutine mcm_state_output