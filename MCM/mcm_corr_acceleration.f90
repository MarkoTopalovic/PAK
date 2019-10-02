      subroutine mcm_corr_acceleration
!************************************************************************
!
!    Purpose: Calculates particle acceleration
!
!  Called by: MOMENTUM
!
!       Date: 12-08-2002
!
!     Errors: 
!
!      Notes: Uses Bonet's mixed correction
!             
!             
!
!************************************************************************
!
use mcm_database
use mcm_math_func
!
implicit none
!
integer:: nn,nn1,i,j,l,k,m,n

REAL(kind=real_acc)    :: Volj, dwdx(mcm_ndim),deltasig(3,3,3),havg,w,dwdr,rhoi
real(kind=real_acc), dimension(3,3) :: grad_v
real(kind=real_acc), dimension(3,3) :: sigmai, qi
!
real(kind=real_acc) :: massj, rhoj, hj, holdj
real(kind=real_acc), dimension(3) :: xj
real(kind=real_acc), dimension(3,3) :: sigmaj, qj
!
real(kind=real_acc) :: sum_w, sum_gradw(3), corr_gradw(3), B(3,3), invB(3,3), dsig
!real(kind=real_acc), dimension(3) :: deltapsi,normalprojection, tangentprojection
!REAL(kind=real_acc)    :: deltapsilength,bndnormlength,dotproduct, normalmagnitude,tangentmagnitude
!
!	loop over all paticles or all velocity particles for the 
!	noncollocated discretization,   particle interactions
!
do i=mcm_svp,mcm_evp		!svp=start velocity point, evp=end velocity point
 par(i)%a = 0.0_d
 !________________________________________________________
 !
 ! Add acceleration due to contact force term if active
 ! Add acceleration due to dyna coupling - tom - 14-11-02
 !________________________________________________________
 !
 if(mcm_contacttype.GT.0) then
  do l = 1,mcm_ndim
   par(i)%a(l) = par(i)%a(l) + par(i)%repulsion(l)
  enddo
 endif
 !________________
 !
 ! End Contact
 !________________
 !
 !
 rhoi = par(i)%rho
 sigmai = par(i)%sigma
 qi = par(i)%q
 !
 !-------------------------------------------------------------------------------------------------
 ! Loop over all neighbours to calculate the sum_w and sum_gradw terms for this neighbourhood
 !
 sum_w = 0
 sum_gradw = 0
 do k=1,par(i)%nnbr + par(i)%g_nnbr
  !
  call mcm_get_j_moment_info(i,k,xj,massj,rhoj,hj,holdj,sigmaj,qj)
  !
  Volj = massj/rhoj
  havg = 0.5_d*(par(i)%h+hj)
  !
  ! Note: kernel centred at j particle
  call mcm_kernel(w,xj,par(i)%x,havg)
  call mcm_gradw(dwdx,dwdr,xj,par(i)%x,havg)
  !call mcm_kernel(w,par(i)%x,xj,havg)
  !call mcm_gradw(dwdx,dwdr,par(i)%x,xj,havg)
  !
  sum_w = sum_w + Volj*w
  do l=1,mcm_ndim
   sum_gradw(l) = sum_gradw(l) + Volj*dwdx(l)
  enddo
  !
 enddo
 !
 !-------------------------------------------------------------------------------------------------
 ! Now calculate B tensor
 !
 B = 0.0_d
 !
 do k=1,par(i)%nnbr + par(i)%g_nnbr
  !
  call mcm_get_j_moment_info(i,k,xj,massj,rhoj,hj,holdj,sigmaj,qj)
  !
  Volj = massj/rhoj
  havg = 0.5_d*(par(i)%h+hj)
  !
  ! Note: kernel centred at j particle
  call mcm_kernel(w,xj,par(i)%x,havg)
  call mcm_gradw(dwdx,dwdr,xj,par(i)%x,havg)
  !call mcm_kernel(w,par(i)%x,xj,havg)
  !call mcm_gradw(dwdx,dwdr,par(i)%x,xj,havg)
  !
  do m=1,mcm_ndim
   corr_gradw(m) = (dwdx(m)*sum_w - w*sum_gradw(m)) / sum_w**2
  enddo
  !
  do m=1,mcm_ndim
   do n = 1,mcm_ndim
    B(n,m) = B(n,m) + Volj * (xj(n) - par(i)%x(n)) * corr_gradw(m)
   enddo
  enddo
  !
 enddo
 !
 ! Calculate inverse of B tensor
 !
 invB = inverse(B)
 !
 !-------------------------------------------------------------------------------------------------
 !
 ! Now calculate momentum equation
 !
 deltasig = 0.0_d
 do k=1,par(i)%nnbr + par(i)%g_nnbr
  !
  call mcm_get_j_moment_info(i,k,xj,massj,rhoj,hj,holdj,sigmaj,qj)
  !
  Volj = massj/rhoj
  havg = 0.5_d*(par(i)%h+hj)
  !
  ! Note: kernel centred at j particle
  call mcm_kernel(w,xj,par(i)%x,havg)
  call mcm_gradw(dwdx,dwdr,xj,par(i)%x,havg)
  !call mcm_kernel(w,par(i)%x,xj,havg)
  !call mcm_gradw(dwdx,dwdr,par(i)%x,xj,havg)
  !
  do m=1,mcm_ndim
   corr_gradw(m) = (dwdx(m)*sum_w - w*sum_gradw(m)) / sum_w**2
  enddo
  !
  do n=1,mcm_ndim
   do m=1,mcm_ndim
    ! Sum form of momentum equation
    dsig = (sigmaj(m,n)-qj(m,n))/(rhoj**2) + (sigmai(m,n)-qi(m,n))/(rhoi**2)
	!
    if(k.eq.par(i)%nnbr) dsig = 0.0_d
	!
    do l=1,mcm_ndim
     deltasig(m,n,l)=deltasig(m,n,l) + dsig*corr_gradw(l)
    enddo 
   enddo
  enddo
 enddo
 !
 do n=1,mcm_ndim
  do m=1,mcm_ndim 
   do l=1,mcm_ndim 
    par(i)%a(n) = par(i)%a(n) + massj * deltasig(n,m,l) * invB(l,m) 
!        deltapsi(l) = 0.0 ! inicijalizacija za racunanje normale
   enddo 
  enddo
 enddo
 
!    if(mcm_contacttype.EQ.0) then
! ! calculation of kernel repulsion for friction
! ! isto kao i ovo gore samo se ne racuna za sve cestice vec samo za kontakt cestice
!  do l=1,mcm_ndim
!    par(i)%kernelrepulsion(l) = 0.0
!  enddo
!  
!  deltapsilength = 0.0
!  do k=1,par(i)%nnbr + par(i)%g_nnbr
!  !
!   

!     do j=1,par(i)%g_ncont
!      if (mcm_g_contlist(j,i).eq.k) then ! k je kontaktna ghost cestica

!      end if
!     enddo
!     do j=1,par(i)%ncont
!      if (mcm_contlist(j,i).eq.k) then ! k je kontaktna realna cestica

!      end if
!     enddo
!     
!      if (contactparticle) then
!          call mcm_get_j_moment_info(i,k,xj,massj,rhoj,hj,holdj,sigmaj,qj)
!          !
!          Volj = massj/rhoj
!          havg = 0.5_d*(par(i)%h+hj)
!          !
!          ! Note: kernel centred at j particle
!          call mcm_kernel(w,xj,par(i)%x,havg)
!          call mcm_gradw(dwdx,dwdr,xj,par(i)%x,havg)
!          !call mcm_kernel(w,par(i)%x,xj,havg)
!          !call mcm_gradw(dwdx,dwdr,par(i)%x,xj,havg)
!          !
!          do m=1,mcm_ndim
!           corr_gradw(m) = (dwdx(m)*sum_w - w*sum_gradw(m)) / sum_w**2
!          enddo
!          !
!          do n=1,mcm_ndim
!           do m=1,mcm_ndim
!            ! Sum form of momentum equation
!            dsig = (sigmaj(m,n)-qj(m,n))/(rhoj**2) + (sigmai(m,n)-qi(m,n))/(rhoi**2)
!	        !
!            if(k.eq.par(i)%nnbr) dsig = 0.0_d
!	        !
!            do l=1,mcm_ndim
!             par(i)%kernelrepulsion(l) = par(i)%kernelrepulsion(l) + dsig*corr_gradw(l)
!            enddo 
!           enddo
!          enddo
!     end if !if (contactparticle) then
!     
! enddo !do k=1,par(i)%nnbr + par(i)%g_nnbr
!  end if !if(mcm_contacttype.EQ.0)
!  
!   !friction section
!    !racunanje normale na osnovu Randles LIbersky rada iz 1996
!  ! Smoothed Particle Hydrodynamics: Some recent improvements and applications
!  ! jednacina 29 psi je materijal 
!  do l=1,mcm_ndim
!   ! sum form of basic SPH momentum equation
!    deltapsi(l) = deltapsi(l) - par(k)%mat*massj*dwdx(l)
!  enddo
!    deltapsilength = sqrt(deltapsi(1)**2+deltapsi(2)**2+deltapsi(3)**2)
!    do k=1,mcm_ndim
!     par(i)%bndnorm(k) = deltapsi(k)/deltapsilength
!    enddo
!   ! sada imamo sve komponente normale  bndnorm
!  bndnormlength = par(i)%bndnorm(1)**2+par(i)%bndnorm(2)**2+par(i)%bndnorm(3)**2 ! |v|**2
!  ! i duzinu normale
!  
!  dotproduct = 0.0 ! inicijalizacija
!  if(mcm_contacttype.EQ.0) then ! za kontakt repulziju
!    dotproduct = par(i)%kernelrepulsion(1)*par(i)%bndnorm(1)+par(i)%kernelrepulsion(2)*par(i)%bndnorm(2)+par(i)%kernelrepulsion(3)*par(i)%bndnorm(3)
!    do k=1,mcm_ndim
!     normalprojection(k) = (dotproduct/bndnormlength)*par(i)%bndnorm(k) !normalna komponenta kontaktne sile u cestici i
!     tangentprojection(k) = par(i)%kernelrepulsion(k)-normalprojection(k)! tangentna komponenta kontaktne sile u cestici i
!    enddo
!  else !if(mcm_contacttype.GT.0)  za Campbellovu repulzivnu silu
!    dotproduct = par(i)%repulsion(1)*par(i)%bndnorm(1)+par(i)%repulsion(2)*par(i)%bndnorm(2)+par(i)%repulsion(3)*par(i)%bndnorm(3)
!    do k=1,mcm_ndim
!     normalprojection(k) = (dotproduct/bndnormlength)*par(i)%bndnorm(k) !normalna komponenta kontaktne sile u cestici i
!     tangentprojection(k) = par(i)%repulsion(k)-normalprojection(k)! tangentna komponenta kontaktne sile u cestici i
!    enddo
!  end if !if(mcm_contacttype.EQ.0)
!  
!  ! sada imamo sve komponente projekcija kontaktne "sile" (ubrzanja) na normalu
!  ! treba da odredimo intentitet normalne i tangentne komponente
!  ! da bi proverili da li je zadovoljen uslov trenja
!    
!    normalmagnitude = sqrt(normalprojection(1)**2 + normalprojection(2)**2 + normalprojection(3)**2)
!    tangentmagnitude = sqrt(tangentprojection(1)**2 + tangentprojection(2)**2 + tangentprojection(3)**2)
!    
!!    if (tangentmagnitude.lt.(normalmagnitude*0.5)) then !hardcoded koeficijent trenja
!!    do l=1,mcm_ndim
!!       par(i)%a(l) = 0 ! staticko trenje ne dozvoljava kretranje
!!    enddo
!!    else
!!    do l=1,mcm_ndim
!!       par(i)%a(l) = par(i)%a(l) - normalprojection(l)*0.3 ! hardcoded dinamicki koeficijent trenja
!!    enddo
!!    end if
!    ! end of friction section
!    
! 
 !
enddo
! 
end subroutine mcm_corr_acceleration
