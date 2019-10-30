subroutine mcm_move
!************************************************************************
!
!    Purpose: update particle velocity and position
!
!  Called by: solution
!
!       Date: 09-08-2002
!
!     Errors: 
!
!      Notes: Originally this routine updated both the velocity and position
!             The velocity update has been moved to a separate routine
!             so that the velocity can be modified/used before the position
!             update without requiring the call to be placed in this routine.
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j,icounter
logical :: belowPlane
real :: mcm_dyingPlane
!
do i=1,mcm_ndim
 mcm_coord_maxmin(1,i) =  1.0e+20_d
 mcm_coord_maxmin(2,i) = -1.0e+20_d
enddo
!
! Update positions and store max and min values for linked list.
!
do i=1,mcm_np
    if(par(i)%lifeStatus.eq.1) then
        par(i)%x(abs(mcm_plane_particle)) = par(i)%x(abs(mcm_plane_particle))+ par(i)%v(abs(mcm_plane_particle))*mcm_dt
        else
 do j=1,mcm_ndim
  par(i)%x(j) = par(i)%x(j) + par(i)%v(j)*mcm_dt
  !
  if(par(i)%x(j).lt.mcm_coord_maxmin(1,j)) mcm_coord_maxmin(1,j) = par(i)%x(j)
  if(par(i)%x(j).gt.mcm_coord_maxmin(2,j)) mcm_coord_maxmin(2,j) = par(i)%x(j)
  !
 enddo
 endif
enddo
!
if(mcm_axopt.eq.4) then
 ! ensure that particles do not pass through axis
 do i=1,mcm_np
  if(par(i)%x(1).le.0.2_d*par(i)%h) then
   par(i)%x(1) = par(i)%x(1) - par(i)%v(1)*mcm_dt
   par(i)%v(1) = 0.0_d
  endif
 enddo
endif
!
! If symmetry planes are active then elastically reflect any particle that has
!  passed through an active symmetry plane
!
if(mcm_boundary) call mcm_check_sym_pen

if((mcm_birth_death.eq.1).or.(mcm_birth_death.eq.3))then 
    if (mcm_np.lt.(mcm_max_np-newBornMax2))then
    icounter = 0
 do i=1,mcm_np

     if(par(i)%lifeStatus.eq.1) then         
        belowPlane = .false.
        
        
        
        if(mcm_bPO.gt.0)then !particle gets reborn if it falls below the plane
        
        
        if (mcm_plane_particle.eq.1)then
            if(par(i)%x(1).le.mcm_birthPlane(1)) then
                 belowPlane = .true. !particle is below x birth plane and needs to be copied
            endif
        endif
        if (mcm_plane_particle.eq.2)then
            if(par(i)%x(2).le.mcm_birthPlane(2)) then
                 belowPlane = .true. !particle is below y birth plane and needs to be copied
            endif
        endif
        if (mcm_plane_particle.eq.3)then
            if(par(i)%x(3).le.mcm_birthPlane(3)) then
                 belowPlane = .true. !particle is below x birth plane and needs to be copied
            endif
        endif
        
        else !particle gets reborn if it passes by the plane
            
            if (mcm_plane_particle.eq.1)then
            if(par(i)%x(1).ge.mcm_birthPlane(1)) then
                 belowPlane = .true. !particle passed x birth plane and needs to be copied
            endif
        endif
        if (mcm_plane_particle.eq.2)then
            if(par(i)%x(2).ge.mcm_birthPlane(2)) then
                 belowPlane = .true. !particle passed y birth plane and needs to be copied
            endif
        endif
        if (mcm_plane_particle.eq.3)then
            if(par(i)%x(3).ge.mcm_birthPlane(3)) then
                 belowPlane = .true. !particle passed x birth plane and needs to be copied
            endif
        endif
            
        endif
        
        
        
        
         
         if (belowPlane) then
          icounter = icounter + 1
          par(mcm_np + icounter)=par(i)
          par(mcm_np + icounter)%x(1) = par(mcm_np + icounter)%xzero(1)
          par(mcm_np + icounter)%x(2) = par(mcm_np + icounter)%xzero(2) 
          par(mcm_np + icounter)%x(3) = par(mcm_np + icounter)%xzero(3)
          
          par(mcm_np + icounter)%p = 0.0_d                   ! pressure
          par(mcm_np + icounter)%e = 0.0_d                   ! TOTAL particle internal energy
          par(mcm_np + icounter)%etry = 0.0_d                ! trial particle internal energy
          par(mcm_np + icounter)%einc = 0.0_d 
          
          par(i)%p = 0.0_d                                   ! pressure
          par(i)%e = 0.0_d                                   ! TOTAL particle internal energy
          par(i)%etry = 0.0_d                                ! trial particle internal energy
          par(i)%einc = 0.0_d 
          
       par(i)%lifeStatus = 0
       endif
      endif
 enddo
 mcm_np = mcm_np + icounter
 mcm_esp=mcm_np
 mcm_evp=mcm_np
else
    do i=1,mcm_np
        par(i)%lifeStatus = 0
        enddo
end if
end if

if(mcm_birth_death.ge.2)then 
icounter = 0
 do i=1,mcm_np
     if(mcm_dPO.gt.0)then !particle dies if it falls below the plane
    
         
         
        if((par(i)%x(mcm_dPO).le.(mcm_deathPlane(mcm_dPO)+par(i)%h)).and.(par(i)%x(mcm_dPO).ge.mcm_deathPlane(mcm_dPO))) then
            par(i)%lifeStatus = 2 !particle is dying, but its not dead yet
        endif
        
        
        if(par(i)%x(mcm_dPO).le.mcm_deathPlane(mcm_dPO)) then
            icounter = icounter + 1
            par(i)%active = .false.
            par(i)%delpointer = i+1
        endif
     else !particle dies if it goes by the plane
         mcm_dyingPlane = mcm_deathPlane(abs(mcm_dPO))-par(i)%h
         if((par(i)%x(abs(mcm_dPO)).ge.mcm_dyingPlane).and.(par(i)%x(abs(mcm_dPO)).le.mcm_deathPlane(abs(mcm_dPO))))then
            par(i)%lifeStatus = 2 !particle is dying, but its not dead yet
        endif
         if(par(i)%x(abs(mcm_dPO)).ge.mcm_deathPlane(abs(mcm_dPO))) then
            icounter = icounter + 1
            par(i)%active = .false.
            par(i)%delpointer = i+1
         endif
    endif
 enddo
 
 do i=1,mcm_np
    if(par(i)%active.eq.(.false.)) then
    do j=i+1,mcm_np
        if(par(j)%active.eq.(.true.)) then
            par(i)=par(j)
            par(j)%active = .false.
            exit
        endif
    enddo
    endif
 enddo
 
 
 mcm_np = mcm_np - icounter
 mcm_esp=mcm_np
 mcm_evp=mcm_np

end if

end subroutine mcm_move