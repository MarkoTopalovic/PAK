      subroutine mcm_oldacceleration
!************************************************************************
!
!    Purpose: Calculates particle acceleration
!
!  Called by: MOMENTUM
!
!       Date: 15-01-99
!
!     Errors: 
!
!      Notes: 
!             
!             
!
!************************************************************************
!
use mcm_database
USE ZFLUID
!
implicit none
!
integer:: nn,nn1,i,j,l,k,m,n,kn,INODEBROJAC
logical:: contactparticle
REAL(kind=real_acc)    :: Volj, dwdx(mcm_ndim),deltasig(mcm_ndim),havg,dwdr,rhoi, deltasigcont(mcm_ndim)
real(kind=real_acc), dimension(3,3) :: grad_v
real(kind=real_acc), dimension(3,3) :: sigmai, qi
!
real(kind=real_acc) :: massj, rhoj, hj, holdj, factor
real(kind=real_acc), dimension(3) :: xj
real(kind=real_acc), dimension(3,3) :: sigmaj, qj
real(kind=real_acc) :: inv_w_deltap, epsilon, fij, ratio,npow
real(kind=real_acc), dimension(3) :: deltapsi,normalprojection, tangentprojection
REAL(kind=real_acc)    :: deltapsilength,bndnormlength,dotproduct, normalmagnitude,tangentmagnitude
!
!	loop over all paticles or all velocity particles for the 
!	noncollocated discretization,   particle interactions
!
! Contact force for writing to disk
do i = 1, mcm_nummat
 mcm_mat(i)%xcf = 0.0
 mcm_mat(i)%ycf = 0.0
 mcm_mat(i)%zcf = 0.0
enddo
!
call mcm_normals ! RACUNANJE NORMALA ZA TRENJE
INODEBROJAC =1
!
!
IF (ALLOCIRANECVORNESILE.EQ.(.FALSE.)) THEN
    ALLOCATE ( NODEBROJ(NPTI))
    ALLOCATE (CVORNESILE(NPTI,3))
    ALLOCIRANECVORNESILE = .TRUE.
END IF
!
!
do i=mcm_svp,mcm_evp		!svp=start velocity point, evp=end velocity point
 par(i)%a = 0.0_d
 !________________________________________________________
 !
 ! Add acceleration due to contact force term if active
 !________________________________________________________
 !
 !  if(mcm_contacttype.EQ.0) then
 ! calculation of kernel repulsion for friction
 ! isto kao i ovo gore samo se ne racuna za sve cestice vec samo za kontakt cestice
 ! do l=1,mcm_ndim
 !   par(i)%kernelrepulsion(l) = 0.0
 ! enddo
  ! deltapsilength = 0.0
 !    do j=1,par(i)%g_ncont
 !     if (mcm_g_contlist(j,i).eq.k) then ! k je kontaktna ghost cestica
 !        contactparticle = .true.
 !     end if
 !    enddo
 !    do j=1,par(i)%ncont
 !     if (mcm_contlist(j,i).eq.k) then ! k je kontaktna realna cestica
 !       contactparticle = .true.
 !     end if
 !    enddo
!      do l=1,mcm_ndim
!        par(i)%kernelrepulsion(l) = par(i)%kernelrepulsion(l) + massj * deltasig(l)
!      enddo
!     end if !if (contactparticle)  
!  enddo !k=1,par(i)%nnbr + par(i)%g_nnbr    
! end if !if(mcm_contacttype.EQ.0)
 !if(mcm_contacttype.EQ.0) then ! za Campbellovu repulzivnu silu
!    dotproduct = par(i)%kernelrepulsion(1)*par(i)%bndnorm(1)+par(i)%kernelrepulsion(2)*par(i)%bndnorm(2)+par(i)%kernelrepulsion(3)*par(i)%bndnorm(3)
!    do k=1,mcm_ndim
!     normalprojection(k) = (dotproduct/bndnormlength)*par(i)%bndnorm(k) !normalna komponenta kontaktne sile u cestici i
!     tangentprojection(k) = par(i)%kernelrepulsion(k)-normalprojection(k)! tangentna komponenta kontaktne sile u cestici i
!    enddo
!  !if(mcm_contacttype.GT.0) za kontakt repulziju
!    dotproduct = par(i)%repulsion(1)*par(i)%bndnorm(1)+par(i)%repulsion(2)*par(i)%bndnorm(2)+par(i)%repulsion(3)*par(i)%bndnorm(3)
!    do k=1,mcm_ndim
!     normalprojection(k) = (dotproduct/bndnormlength)*par(i)%bndnorm(k) !normalna komponenta kontaktne sile u cestici i
!     tangentprojection(k) = par(i)%repulsion(k)-normalprojection(k)! tangentna komponenta kontaktne sile u cestici i
!    enddo
! else !if(mcm_contacttype.GT.0) za kontakt repulziju
!    dotproduct = par(i)%repulsion(1)*par(i)%bndnorm(1)+par(i)%repulsion(2)*par(i)%bndnorm(2)+par(i)%repulsion(3)*par(i)%bndnorm(3)
!    do k=1,mcm_ndim
!     normalprojection(k) = (dotproduct/bndnormlength)*par(i)%bndnorm(k) !normalna komponenta kontaktne sile u cestici i
!     tangentprojection(k) = par(i)%repulsion(k)-normalprojection(k)! tangentna komponenta kontaktne sile u cestici i
!    enddo
!  end if !if(mcm_contacttype.EQ.0)
  
 if(mcm_contacttype.GT.0) then
  do l = 1,mcm_ndim
      !par(i)%a(l) = par(i)%a(l) + par(i)%repulsion(l)
    !par(i)%contactforce(l) = 0.0_d
    !par(i)%normalforce(l) = 0.0_d
    !par(i)%repulsion(l) = 0.0_d
      !     normalprojection(l) = (dotproduct/par(i)%bndnormlength)*par(i)%bndnorm(l) !normalna komponenta kontaktne sile u cestici i
!     !normalprojection(l) = (dotproduct/bndnormlength)*(-1)*par(i)%bndnorm(l) !normalna komponenta kontaktne sile u cestici i
!     tangentprojection(l) = deltasigcont(l)-normalprojection(l)! tangentna komponenta kontaktne sile u cestici i
!        !par(i)%contactforce(l) = par(i)%contactforce(l)+massj *tangentprojection(l)
!        !par(i)%normalforce(l) = par(i)%normalforce(l)+factor*par(i)%bndnorm(l)
!        !par(i)%normalforce(l) = factor*par(i)%bndnorm(l)
!        !par(i)%repulsion(l) = par(i)%repulsion(l) + massj * deltasigcont(l)
      !        par(i)%contactforce(l) = par(i)%contactforce(l)+tangentprojection(l)
!        par(i)%normalforce(l) = par(i)%normalforce(l)+normalprojection(l)
  bndnormlength = par(i)%bndnorm(1)**2+par(i)%bndnorm(2)**2+par(i)%bndnorm(3)**2
    dotproduct = par(i)%repulsion(1)*par(i)%bndnorm(1)+par(i)%repulsion(2)*par(i)%bndnorm(2)+par(i)%repulsion(3)*par(i)%bndnorm(3)
    do k=1,mcm_ndim
     normalprojection(k) = (dotproduct/bndnormlength)*par(i)%bndnorm(k) !normalna komponenta kontaktne sile u cestici i
     tangentprojection(k) = par(i)%repulsion(k)-normalprojection(k)! tangentna komponenta kontaktne sile u cestici i
    enddo
    ! sada imamo sve komponente projekcija kontaktne "sile" (ubrzanja) na normalu
  ! treba da odredimo intentitet normalne i tangentne komponente
  ! da bi proverili da li je zadovoljen uslov trenja
    normalmagnitude = sqrt(normalprojection(1)**2 + normalprojection(2)**2 + normalprojection(3)**2)
    tangentmagnitude = sqrt(tangentprojection(1)**2 + tangentprojection(2)**2 + tangentprojection(3)**2)
!    if (tangentmagnitude.lt.(normalmagnitude*0.5)) then !hardcoded koeficijent trenja
!    do l=1,mcm_ndim
!       par(i)%a(l) = 0 ! staticko trenje ne dozvoljava kretranje
!    enddo
!    else
!    do l=1,mcm_ndim
!!       par(i)%a(l) = par(i)%a(l) - normalprojection(l)*0.3 ! hardcoded dinamicki koeficijent trenja
!    enddo
!    end if
    if (tangentmagnitude.lt.(normalmagnitude*mcm_staticfriction)) then
    par(i)%a(l) = par(i)%a(l) + par(i)%repulsion(l)
    else
    par(i)%a(l) = par(i)%a(l) + normalprojection(l) + tangentprojection(l) * mcm_dynamicfriction
    !
    end if
  enddo
  mcm_mat(par(i)%mat)%xcf = mcm_mat(par(i)%mat)%xcf + par(i)%mass*par(i)%repulsion(1)
  mcm_mat(par(i)%mat)%ycf = mcm_mat(par(i)%mat)%ycf + par(i)%mass*par(i)%repulsion(2)
  mcm_mat(par(i)%mat)%zcf = mcm_mat(par(i)%mat)%zcf + par(i)%mass*par(i)%repulsion(3)
  ! za Campbellovu repulzivnu silu
! ! ovo se kopira u paks
  NODEBROJ(INODEBROJAC) = i
  CVORNESILE(INODEBROJAC,1) = par(i)%repulsion(1)
  CVORNESILE(INODEBROJAC,2) = par(i)%repulsion(2)
  CVORNESILE(INODEBROJAC,3) = par(i)%repulsion(3)
! ! ovo se kopira u paks
  ! za Campbellovu repulzivnu silu
! Campbelova sila se ne ponasa lepo, bolje ponasanje daje default kernel  
 endif
! OVO JE INICIJALIZACIJA ZFLUID MODULA KOJIME SE PRENOSE SILE U PAK
if((par(i).mat.eq.mcm_nummat).and.(mcm_contacttype.EQ.0).and.(mcm_kojpak.ne.4)) then
    NODEBROJ(INODEBROJAC) = i
    CVORNESILE(INODEBROJAC,1) = 0
    CVORNESILE(INODEBROJAC,2) = 0
    CVORNESILE(INODEBROJAC,3) = 0
end if
 !
 rhoi = par(i)%rho
 sigmai = par(i)%sigma
 qi = par(i)%q
 ! dotproduct = 0.0_d
 do k=1,par(i)%nnbr + par(i)%g_nnbr
  !contactparticle = .false.
  call mcm_get_j_moment_info(i,k,xj,massj,rhoj,hj,holdj,sigmaj,qj)
  !
  Volj = massj/rhoj
  havg = 0.5_d*(par(i)%h+hj)
  deltasig = 0.0_d
  deltasigcont = 0.0_d
  normalprojection = 0.0_d
  tangentprojection = 0.0_d
  !dotproduct = 0.0_d
  call mcm_gradw(dwdx,dwdr,par(i)%x,xj,havg)
  j = mcm_nbrlist(k,i)
  do n=1,mcm_ndim
   do m=1,mcm_ndim
    ! Sum form of basic SPH momentum equation
    if(par(i)%mat.eq.par(j)%mat) then ! ISTI MATERIJAL
    deltasig(n) = deltasig(n)+ ( (sigmaj(m,n)-qj(m,n))/(rhoj**2)  +          &
	                             (sigmai(m,n)-qi(m,n))/(rhoi**2)) * dwdx(m)
   else ! RAZLICITI KONTAKT MATERIJAL
   deltasigcont(n) = deltasigcont(n)+ ( (sigmaj(m,n)-qj(m,n))/(rhoj**2)  +          &
	                             (sigmai(m,n)-qi(m,n))/(rhoi**2)) * dwdx(m)
   !dotproduct = dotproduct + deltasigcont(n)*par(i)%bndnorm(n) 
   !contactparticle = .true.
   end if
   enddo!do m=1,mcm_ndim
  enddo!do n=1,mcm_ndim
  !
!   RASTAVLJANJE KONTAKTNE SILE OD CESTICE (RAZLICITOG MATERIJALA) NA NORMALNU I TANGENTNU KOMPONENTU
!   ZA POSMATRANU CESTICU I IMAMO K KONTAKTNIH CESTICA, AKO JE RAZLICIT MATERIJAL deltasigcont(n) <> 0
  bndnormlength = par(i)%bndnorm(1)**2+par(i)%bndnorm(2)**2+par(i)%bndnorm(3)**2
  if (bndnormlength.gt.0) then
   dotproduct = deltasigcont(1)*par(i)%bndnorm(1)+deltasigcont(2)*par(i)%bndnorm(2)+deltasigcont(3)*par(i)%bndnorm(3)
   !   factor = massj*(dotproduct/par(i)%bndnormlength)
    do l=1,mcm_ndim
     normalprojection(l) = (dotproduct/bndnormlength)*par(i)%bndnorm(l) !normalna komponenta kontaktne sile u cestici i
     tangentprojection(l) = deltasigcont(l)-normalprojection(l)! tangentna komponenta kontaktne sile u cestici i
    enddo !l=1,mcm_ndim
!   ZA PROVERU USLOVA TRENJA POTREBNI SU INTENYITETI KOMPONENTI  
    normalmagnitude = sqrt(normalprojection(1)**2 + normalprojection(2)**2 + normalprojection(3)**2)
    tangentmagnitude = sqrt(tangentprojection(1)**2 + tangentprojection(2)**2 + tangentprojection(3)**2)
!   DODAVANJE KONTAKTNE SILE U ZAVISNOSTI OD USLOVA TRENJA      
  do l=1,mcm_ndim
!   sum form of basic SPH momentum equation
!   par(i)%a(l) = par(i)%a(l) + massj * deltasig(l)
    if (tangentmagnitude.lt.(normalmagnitude*mcm_staticfriction)) then
    par(i)%a(l) = par(i)%a(l) + massj * deltasig(l) + massj * deltasigcont(l)
    else
    par(i)%a(l) = par(i)%a(l) + massj * deltasig(l) + massj * normalprojection(l) + massj * tangentprojection(l) * mcm_dynamicfriction
    !   ovde je pretpostavka da se cestica krece u suprotnom pravcu od tangentne komponente sile reakcije
    !   tom kretanju se suprotsavlja sila dinamickog trenja miK*Fn
    !   Fn = Tf/miS dakle u mcm_dynamicfriction figurise i 1/miS
    end if
  enddo !l=1,mcm_ndim
  endif !if (bndnormallength.gt.0) then
!******************************************** P A K *******************************************
! ovo se kopira u paks
! CESTICAMA KOPIRANIM IZ PAKA SU JE DODELJEN POSLEDNJI MATERIJAL  
    if((par(i).mat.eq.mcm_nummat).and.(mcm_kojpak.ne.4)) then
        NODEBROJ(INODEBROJAC) = i
        CVORNESILE(INODEBROJAC,1) = CVORNESILE(INODEBROJAC,1) + massj * deltasigcont(1)
        CVORNESILE(INODEBROJAC,2) = CVORNESILE(INODEBROJAC,2) + massj * deltasigcont(2)
        CVORNESILE(INODEBROJAC,3) = CVORNESILE(INODEBROJAC,3) + massj * deltasigcont(3)
    end if
!    
 enddo !k=1,par(i)%nnbr + par(i)%g_nnbr ! KRAJ PETLJE PO KONTAKTIMA
! ovo se kopira u paks
!******************************************** P A K *******************************************  
    if(par(i).mat.eq.mcm_nummat) then
        INODEBROJAC = INODEBROJAC+1 ! brojac za pak ZFLUID modul
    end if
enddo ! do i=mcm_svp,mcm_evp
! 
end subroutine mcm_oldacceleration
!
!
subroutine mcm_normals
!!!!!!!!!!!!!**********************************************************************************************
!
!           ZA RACUNANJE NORMALA ZA TRETMAN TRENJA U KONTAKTU
!           PRVO SE ZA SVE CESTICE NADJE DEFAULT NORMALA NA OSNOVU "BOJE" A POTOM SE NORMALE UNAPREDJUJU
!           (AKO JE TO MOGUCE)
!
!!!!!!!!!!!!!!!!!!!!***************************************************************************************
use mcm_database
!
implicit none
!
integer:: nn,nn1,i,j,l,k,m,n,kn,idebug
!
real(kind=real_acc), dimension(3) :: deltapsi,normalprojection, tangentprojection
REAL(kind=real_acc)    :: deltapsilength,bndnormlength,dotproduct, normalmagnitude,tangentmagnitude
!
real(kind=real_acc) :: distance,distance1,distance2, distance3, tempdistance
integer::   closestid1,closestid2,closestid3, tempid
integer, dimension(10) :: initialmaterialparticle
real(kind=real_acc), dimension(3) :: vectora, vectorb
real(kind=real_acc), dimension(3) :: xj
real(kind=real_acc) :: massj, rhoj, hj, holdj
real(kind=real_acc), dimension(3,3) :: sigmaj, qj
real(kind=real_acc)    :: Volj, dwdx(mcm_ndim),havg,w,dwdr,rhoi
logical:: normalfound
logical:: onecontact
logical:: twocontact
logical:: threecontact
do i = 1, mcm_nummat
!     initialmaterialparticle(i) = 0
 averagematnormal(i,1) = 0.0_d
 averagematnormal(i,2) = 0.0_d
 averagematnormal(i,3) = 0.0_d
enddo
normalfound = .false.
!
!
!racunanje normale na osnovu Randles LIbersky rada iz 1996
  ! Smoothed Particle Hydrodynamics: Some recent improvements and applications
  ! jednacina 29 psi je materijal 
  ! DEFAULT NORMALA POSLE SE U NAREDNIM PETLJAMA POBOLJSAVA AKO JE TO MOGUCE
do i=mcm_svp,mcm_evp		!svp=start velocity point, evp=end velocity point
  deltapsi(1) = 0.0_d
  deltapsi(2) = 0.0_d
  deltapsi(3) = 0.0_d
  idebug = par(i)%nnbr + par(i)%g_nnbr
  do k=1,par(i)%nnbr + par(i)%g_nnbr
  !
  call mcm_get_j_moment_info(i,k,xj,massj,rhoj,hj,holdj,sigmaj,qj)
    Volj = massj/rhoj
  havg = 0.5_d*(par(i)%h+hj)
 call mcm_kernel(w,xj,par(i)%x,havg)
  call mcm_gradw(dwdx,dwdr,xj,par(i)%x,havg)
    do l=1,mcm_ndim
   ! sum form of basic SPH momentum equation
    deltapsi(l) = deltapsi(l) - par(k)%mat*massj*dwdx(l)
    enddo
  end do !k=1,par(i)%nnbr + par(i)%g_nnbr
    deltapsilength = sqrt(deltapsi(1)**2+deltapsi(2)**2+deltapsi(3)**2)
    do l=1,mcm_ndim
    if (deltapsilength.gt.0) then
     par(i)%bndnorm(l) = deltapsi(l)/deltapsilength
     !  bndnormlength = par(i)%bndnorm(1)**2+par(i)%bndnorm(2)**2+par(i)%bndnorm(3)**2 ! |v|**2
  ! i duzinu normale
     !averagematnormal(par(i)%mat,l) = averagematnormal(par(i)%mat,l) + (par(i)%bndnorm(l))
    end if 
    enddo
!   ! sada imamo sve komponente normale  bndnorm
enddo ! do i=mcm_svp,mcm_evp 
!  averagematnormal SLUZI ZA USKLADJIVANJE SMERA, POTREBNO JE DA SVE NORMALE BUDU U ISTOM SMERU (PROSECNE NORMALE)
!  
!do i=mcm_svp,mcm_evp
!!  !provera pravca ako je skalarni proizvod veci od 0 onda su istog pravca
!  dotproduct = 0.0 ! inicijalizacija
!!!   ! svi vektori normale treba da budu istog smera kao srednji mega vektor (nisam ga skalirao/delio sa brojem cestica) jer mi samo treba pravac
!  dotproduct = averagematnormal(par(i)%mat,1)*par(i)%bndnorm(1) + averagematnormal(par(i)%mat,2)*par(i)%bndnorm(2) + averagematnormal(par(i)%mat,3)*par(i)%bndnorm(3)
!   dotproduct = par(initialmaterialparticle(par(i).mat))%bndnorm(1)*par(i)%bndnorm(1) + par(initialmaterialparticle(par(i).mat))%bndnorm(2)*par(i)%bndnorm(2) + par(initialmaterialparticle(par(i).mat))%bndnorm(3)*par(i)%bndnorm(3)
! svi vektori treba da budu istog smera kao vektor normale prve cestice od tog materijala
!      if (dotproduct.lt.0) then ! vektor suprotnog smera
!      par(i)%bndnorm(1) = -par(i)%bndnorm(1)
!      par(i)%bndnorm(2) = -par(i)%bndnorm(2)
!      par(i)%bndnorm(3) = -par(i)%bndnorm(3)
!      end if
!
!enddo ! do i=mcm_svp,mcm_evp  
!kraj 
!racunanja normale na osnovu Randles LIbersky rada iz 1996
! Smoothed Particle Hydrodynamics: Some recent improvements and applications
! 
!  
!*****************************************************************************
!   STARO RESENJE KOJE NIJE DAVALO DOVOLJNO DOBRE NORMALE
!  ! racunanje normale 2.0 na osnovu 3 tacke drugog materijala 
! 
!do i=mcm_svp,mcm_evp 
!if (initialmaterialparticle(par(i).mat).eq.0)then
!  distance = 0.0
!  distance1 = 0.0
!  distance2  = 0.0
!  distance3 = 0.0
! distance1 = 10000.0
!  distance2  = 10000.0
!  distance3 = 10000.0
!  tempdistance = 0.0
!  closestid1 = 0
!  closestid2 = 0
!  closestid3 = 0
!  tempid = 0
!    onecontact = .false.
!    twocontact = .false.
!    threecontact = .false.
!    
!  ! trazi najblizu cesticu od razlicitog materijala u skupu kontaktnih cestica  
!    do k=1,par(i)%nnbr + par(i)%g_nnbr
!    j = mcm_nbrlist(k,i)
!    if (par(i)%mat.ne.par(j)%mat) then
!        distance =sqrt((par(i)%x(1)-par(j)%x(1))**2+(par(i)%x(2)-par(j)%x(2))**2+(par(i)%x(3)-par(j)%x(3))**2)
!        if ((distance.lt.distance1).or.(onecontact.eq.(.false.)))then ! tekuca cestica k je najbliza
!        onecontact = .true.
!        closestid1 = j ! id najblize cestice od razlicitog materijala
!        distance1 = distance
!    end if
!    end if
!    end do !do k=1,par(i)%nnbr + par(i)%g_nnbr
!
!if (onecontact.eq.(.true.))then ! ako nema najblizu cesticu nema potrebe da trazim drugu, ako ima prvu, drugu trazim po istom principu kao i prvu
!    do k=1,par(i)%nnbr + par(i)%g_nnbr
!        j = mcm_nbrlist(k,i)
!        if ((par(i)%mat.ne.par(j)%mat).and.(j.ne.closestid1)) then ! skup svih susednih cestica od razlicitog materijala osim najblize cestice
!        distance =sqrt((par(i)%x(1)-par(j)%x(1))**2+(par(i)%x(2)-par(j)%x(2))**2+(par(i)%x(3)-par(j)%x(3))**2)
!***************************************
      !  closestid3 = closestid2 !druga najbliza cestica postaje treca 
      !distance3 = distance2
      !closestid2 = closestid1 ! stara prva cestica postaje druga
      !distance2 = distance1
      !closestid1 = j ! id najblize cestice od razlicitog materijala
      !distance1 = distance
      !else if (distance.lt.distance2)then ! tekuca cestica je bliza od druge ali nije bliza od trece
      !closestid3 = closestid2 !druga najbliza cestica postaje treca 
      !distance3 = distance2
      !closestid2 = j ! id druge najblize cestice od razlicitog materijala
      !distance2 = distance
      !  else if (distance.lt.distance3)then ! tekuca cestica je bliza samo od trece cestice
      !      closestid3 = j ! id druge najblize cestice od razlicitog materijala
      !    distance3 = distance
      !  end if
!par(i)%bndnorm(4) = par(i)%bndnorm(1)*par(closestid1)%x(1) + par(i)%bndnorm(2)*par(closestid1)%x(2) + par(i)%bndnorm(3)*par(closestid1)%x(3)
!************************************************************
!        if ((distance.lt.distance2).or.(twocontact.eq.(.false.))) then ! tekuca cestica k je druga najbliza
!            twocontact = .true.      
!            closestid2 = j ! id druge najblize cestice od razlicitog materijala
!            distance2 = distance
!        end if
!        end if
!    end do !do k=1,par(i)%nnbr + par(i)%g_nnbr
!end if
!!  
! if ((onecontact.eq.(.true.)).and.(twocontact.eq.(.true.)))then ! da bi trazio trecu najblizu moram da postoje prve dve
!    do k=1,par(i)%nnbr + par(i)%g_nnbr
!        j = mcm_nbrlist(k,i)
!        if ((par(i)%mat.ne.par(j)%mat).and.(j.ne.closestid1).and.(j.ne.closestid2)) then ! skup svih susednih cestica od razlicitog materijala osim najblize cestice
!        distance =sqrt((par(i)%x(1)-par(j)%x(1))**2+(par(i)%x(2)-par(j)%x(2))**2+(par(i)%x(3)-par(j)%x(3))**2)
!        if ((distance.lt.distance2).or.(threecontact.eq.(.false.))) then ! tekuca cestica k je druga najbliza
!            threecontact = .true.      
!            closestid3 = j ! id druge najblize cestice od razlicitog materijala
!            distance3 = distance
!        end if
!        end if
!    end do !do k=1,par(i)%nnbr + par(i)%g_nnbr
!end if 
!!  
!  if (threecontact.eq.(.true.)) then
!!!  ! za posmatranu cesticu i sada imamo 3 najblize cestice od razlicitog materijala 
!!!  !ako uslov nije ispunjen cestica nema 3 bliske cestice od razlicitog materijala
! ako jeste tu cesticu uzimamo za referentnu i njenu normalu dodeljujemo svim cesticama od tog materijala
!  initialmaterialparticle(par(i).mat) = i
!!!  ! od te 3 cestice napravimo 2 vektora
!   vectora(1) = par(closestid1)%x(1)-par(closestid3)%x(1)
!   vectora(2) = par(closestid1)%x(2)-par(closestid3)%x(2)
!   vectora(3) = par(closestid1)%x(3)-par(closestid3)%x(3)
!  
!   vectorb(1) = par(closestid2)%x(1)-par(closestid3)%x(1)
!   vectorb(2) = par(closestid2)%x(2)-par(closestid3)%x(2)
!   vectorb(3) = par(closestid2)%x(3)-par(closestid3)%x(3)
!!!
!!!  !normala na posmatranu cesticu je vektorski proizvod ova 2 vektkora
!  par(i)%bndnorm(1) = vectora(2)*vectorb(3)-vectora(3)*vectorb(2)
!  par(i)%bndnorm(2) = vectora(3)*vectorb(1)-vectora(1)*vectorb(3)
!  par(i)%bndnorm(3) = vectora(1)*vectorb(2)-vectora(2)*vectorb(1)
!!!  
!  bndnormlength = par(i)%bndnorm(1)**2+par(i)%bndnorm(2)**2+par(i)%bndnorm(3)**2
!  bndnormlength = sqrt(bndnormlength)
!  
!  if (bndnormlength.ne.0) then
!  
!  par(i)%bndnorm(1) = par(i)%bndnorm(1)/bndnormlength
!  par(i)%bndnorm(2) = par(i)%bndnorm(2)/bndnormlength
!  par(i)%bndnorm(3) = par(i)%bndnorm(3)/bndnormlength
!
!  end if
!!
!    end if
!!  
!enddo ! do i=mcm_svp,mcm_evp
!
!
!  ! racunanje normale 3.0 na osnovu 3 tacke koje su prenete iz paka (poslednji materijal)
!   posto su cestice nastale od ljuske ove cestice su u jednoj povrsi 
!   ovo izgleda radi samo za slucaj kada je mcm povezan sa pakom i kada su pak elementi ljuska
if (mcm_kojpak.eq.5)then
do i=mcm_svp,mcm_evp  
  distance = 0.0
  distance1 = 0.0
  distance2  = 0.0
  distance3 = 0.0
  tempdistance = 0.0
  closestid1 = 0
  closestid2 = 0
  closestid3 = 0
  tempid = 0
  onecontact = .false.
  twocontact = .false.
  threecontact = .false.
!   
  if (par(i)%mat.eq.mcm_nummat) then
!   
! trazi najblizu cesticu od istog materijala u skupu kontaktnih cestica  
      idebug = par(i)%nnbr + par(i)%g_nnbr
    do k=1,par(i)%nnbr + par(i)%g_nnbr
    j = mcm_nbrlist(k,i)
    if (par(i)%mat.eq.par(j)%mat) then
        distance =sqrt((par(i)%x(1)-par(j)%x(1))**2+(par(i)%x(2)-par(j)%x(2))**2+(par(i)%x(3)-par(j)%x(3))**2)
        if ((distance.lt.distance1).or.(onecontact.eq.(.false.)))then ! tekuca cestica k je najbliza
        onecontact = .true.
        closestid1 = j ! id najblize cestice od razlicitog materijala
        distance1 = distance
    end if
    end if
    end do !do k=1,par(i)%nnbr + par(i)%g_nnbr
!
if (onecontact.eq.(.true.))then ! ako nema najblizu cesticu nema potrebe da trazim drugu, ako ima prvu, drugu trazim po istom principu kao i prvu
    idebug = par(i)%nnbr + par(i)%g_nnbr
    do k=1,par(i)%nnbr + par(i)%g_nnbr
        j = mcm_nbrlist(k,i)
        if ((par(i)%mat.eq.par(j)%mat).and.(j.ne.closestid1)) then ! skup svih susednih cestica od razlicitog materijala osim najblize cestice
        distance =sqrt((par(i)%x(1)-par(j)%x(1))**2+(par(i)%x(2)-par(j)%x(2))**2+(par(i)%x(3)-par(j)%x(3))**2)
        if ((distance.lt.distance2).or.(twocontact.eq.(.false.))) then ! tekuca cestica k je druga najbliza
            twocontact = .true.      
            closestid2 = j ! id druge najblize cestice od razlicitog materijala
            distance2 = distance
        end if
        end if
    end do !do k=1,par(i)%nnbr + par(i)%g_nnbr
end if
!  
 if ((onecontact.eq.(.true.)).and.(twocontact.eq.(.true.)))then ! da bi trazio trecu najblizu moram da postoje prve dve
    do k=1,par(i)%nnbr + par(i)%g_nnbr
        j = mcm_nbrlist(k,i)
        if ((par(i)%mat.eq.par(j)%mat).and.(j.ne.closestid1).and.(j.ne.closestid2)) then ! skup svih susednih cestica od razlicitog materijala osim najblize cestice
        distance =sqrt((par(i)%x(1)-par(j)%x(1))**2+(par(i)%x(2)-par(j)%x(2))**2+(par(i)%x(3)-par(j)%x(3))**2)
        if ((distance.lt.distance2).or.(threecontact.eq.(.false.))) then ! tekuca cestica k je druga najbliza
            threecontact = .true.      
            closestid3 = j ! id druge najblize cestice od razlicitog materijala
            distance3 = distance
        end if
        end if
    end do !do k=1,par(i)%nnbr + par(i)%g_nnbr
end if 
!   AKO NEMA 3 NAJBLIZE CESTICE VEC SAMO 2 UZIMAMO POSMATRANU CESTICU KAO TRECU
if ((threecontact.eq.(.false.)).and.(twocontact.eq.(.true.))) then
    threecontact = .true.
    closestid3 = i
end if
!
! 
!  
  if (threecontact.eq.(.true.)) then
!!  ! za posmatranu cesticu i treba da postoje 2 najblize cestice od istog materijala 
!!  !ako uslov nije ispunjen nesto je otkinulo cesticu
!!  ! od ove 3 cestice napravimo 2 vektora
   vectora(1) = par(closestid1)%x(1)-par(closestid3)%x(1)
   vectora(2) = par(closestid1)%x(2)-par(closestid3)%x(2)
   vectora(3) = par(closestid1)%x(3)-par(closestid3)%x(3)
!  
   vectorb(1) = par(closestid2)%x(1)-par(closestid3)%x(1)
   vectorb(2) = par(closestid2)%x(2)-par(closestid3)%x(2)
   vectorb(3) = par(closestid2)%x(3)-par(closestid3)%x(3)
!!
!!  !normala na posmatranu cesticu je vektorski proizvod ova 2 vektkora
  par(i)%bndnorm(1) = vectora(2)*vectorb(3)-vectora(3)*vectorb(2)
  par(i)%bndnorm(2) = vectora(3)*vectorb(1)-vectora(1)*vectorb(3)
  par(i)%bndnorm(3) = vectora(1)*vectorb(2)-vectora(2)*vectorb(1)
!!  
  bndnormlength = par(i)%bndnorm(1)**2+par(i)%bndnorm(2)**2+par(i)%bndnorm(3)**2
  bndnormlength = sqrt(bndnormlength)
!  
  if (bndnormlength.ne.0) then 
      par(i)%bndnorm(1) = par(i)%bndnorm(1)/bndnormlength
      par(i)%bndnorm(2) = par(i)%bndnorm(2)/bndnormlength
      par(i)%bndnorm(3) = par(i)%bndnorm(3)/bndnormlength
  end if
!
 end if
!    
end if
! 
enddo ! do i=mcm_svp,mcm_evp


!
! RACUNANJE NORMALE 4.0 
!ZA CESTICU PESKA SE TRAZI NAJBLIZA CESTICA PLOCE 
! I ONDA SE CESTICI PESKA DODELJUJE NJENA NORMALA
do i=mcm_svp,mcm_evp  
  distance = 0.0
  distance1 = 0.0
  tempdistance = 0.0
  closestid1 = 0
  onecontact = .false.
!    
  if (par(i)%mat.ne.mcm_nummat) then ! nije ploca onda mora da je pesak
!    
  ! trazi najblizu cesticu od ploce 
      idebug = par(i)%nnbr + par(i)%g_nnbr
    do k=1,par(i)%nnbr + par(i)%g_nnbr
    j = mcm_nbrlist(k,i)
        if (par(i)%mat.ne.par(j)%mat) then
            distance =sqrt((par(i)%x(1)-par(j)%x(1))**2+(par(i)%x(2)-par(j)%x(2))**2+(par(i)%x(3)-par(j)%x(3))**2)
            if ((distance.lt.distance1).or.(onecontact.eq.(.false.)))then ! tekuca cestica k je najbliza
            onecontact = .true.
            closestid1 = j ! id najblize cestice od razlicitog materijala
            distance1 = distance
            end if
        end if
    end do !do k=1,par(i)%nnbr + par(i)%g_nnbr

  if (onecontact.eq.(.true.)) then
  par(i)%bndnorm(1) = par(closestid1)%bndnorm(1)
  par(i)%bndnorm(2) = par(closestid1)%bndnorm(2)
  par(i)%bndnorm(3) = par(closestid1)%bndnorm(3)
  end if
!
end if  
!  
enddo ! do i=mcm_svp,mcm_evp
!
endif
do i=mcm_svp,mcm_evp
!
! USKLADJIVANJE SMERA PREMA Z OSI    
  if (par(i)%bndnorm(3).lt.0) then
   par(i)%bndnorm(1) = -par(i)%bndnorm(1)
   par(i)%bndnorm(2) = -par(i)%bndnorm(2)
   par(i)%bndnorm(3) = -par(i)%bndnorm(3)    
  end if
!
  !  par(i)%bndnorm(1) = par(initialmaterialparticle(par(i).mat))%bndnorm(1)
  !  par(i)%bndnorm(2) = par(initialmaterialparticle(par(i).mat))%bndnorm(2)
  !  par(i)%bndnorm(3) = par(initialmaterialparticle(par(i).mat))%bndnorm(3)
  !  par(i)%bndnorm(4) = par(initialmaterialparticle(par(i).mat))%bndnorm(4)
  !  par(i)%bndnormlength = par(initialmaterialparticle(par(i).mat))%bndnormlength 
enddo ! do i=mcm_svp,mcm_evp
!
end subroutine mcm_normals