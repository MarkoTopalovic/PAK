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
!
do i=1,mcm_ndim
 mcm_coord_maxmin(1,i) =  1.0e+20_d
 mcm_coord_maxmin(2,i) = -1.0e+20_d
enddo
!
! Update positions and store max and min values for linked list.
!
do i=1,mcm_np
 do j=1,mcm_ndim
  par(i)%x(j) = par(i)%x(j) + par(i)%v(j)*mcm_dt
  !
  if(par(i)%x(j).lt.mcm_coord_maxmin(1,j)) mcm_coord_maxmin(1,j) = par(i)%x(j)
  if(par(i)%x(j).gt.mcm_coord_maxmin(2,j)) mcm_coord_maxmin(2,j) = par(i)%x(j)
  !
 enddo
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

if(mcm_emitter.eq.(.true.))then ! hardcoded
    if (mcm_np.lt.(mcm_max_np-100))then
    icounter = 0
 do i=1,mcm_np
      !if((par(i)%x(2).le.199).and.(par(i)%mat.eq.2).and.(par(i)%newborn.eq.(.true.))) then !casa v.1.0
     if ((i.gt.825).and.(i.lt.844))then
         par(i)%newborn = .false. !casa v.2.0
         endif
      if((par(i)%x(2).le.0.01).and.(par(i)%mat.eq.2).and.(par(i)%newborn.eq.(.true.))) then !casa v.2.0
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
          
       par(i)%newborn = .false.
      endif
 enddo
 mcm_np = mcm_np + icounter
 mcm_esp=mcm_np
 mcm_evp=mcm_np
else
    do i=1,mcm_np
        par(i)%newborn = .false.
        enddo
end if
end if

end subroutine mcm_move