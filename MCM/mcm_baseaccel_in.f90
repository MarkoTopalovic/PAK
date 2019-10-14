subroutine mcm_baseaccel_in
!************************************************************************
!
!    Purpose: Read base acceleration values or set to zero
!
!  Called by: getinput
!
!       Date: 9/5/2005
!
!     Errors: 
!
!      Notes: 
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer lcount
character (len=80):: mssg,textline
!
if(mcm_baseaccel) then
 mssg='Error reading base accelerations card'
 call mcm_gttxsg(textline,lcount)
 !
 read(Unit=textline,fmt=100,err=500) mcm_base_a(1), mcm_base_a(2), mcm_base_a(3)
 !
endif
!
if(mcm_nthpx.ne.1) mcm_base_a(1) = 0.0_d
if(mcm_nthpy.ne.1) mcm_base_a(2) = 0.0_d
if(mcm_nthpz.ne.1) mcm_base_a(3) = 0.0_d
!
return
!
100 format(3e20.0)
500 call mcm_termin(textline,mssg,lcount,1)
!
    end subroutine mcm_baseaccel_in
    
    
subroutine mcm_birth_death_in
!************************************************************************
!
!    Purpose: Read birt and death particles
!
!  Called by: getinput
!
!       Date: 10/04/2019
!
!     Errors: 
!
!      Notes: 
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer lcount,i
integer startNewBorn
integer endNewBorn
character (len=80):: mssg,textline
!
if(mcm_birth_death.eq.1) then
 mssg='Error reading birth particle card'
 call mcm_gttxsg(textline,lcount)
 !
 read(Unit=textline,fmt=100,err=500) startNewBorn, endNewBorn
 newBornMax2 = 2*(endNewBorn - startNewBorn)
 !
endif
!
do i=1,mcm_np
 if ((i.ge.startNewBorn).and.(i.le.endNewBorn))then
  par(i)%newborn = .true.
 else
     par(i)%newborn = .false.
 endif
  par(i)%contactparticle = 0 ! initially not in contact with other material
enddo
!
if (mcm_plane_particle.gt.0) then
    call mcm_gttxsg(textline,lcount)
 !
    read(Unit=textline,fmt=101,err=500) mcm_birthPlane(1), mcm_birthPlane(2), mcm_birthPlane(3)
    endif
return
!
100 format(2i10)
    101 format(3e20.0)
500 call mcm_termin(textline,mssg,lcount,1)
!
end subroutine mcm_birth_death_in
