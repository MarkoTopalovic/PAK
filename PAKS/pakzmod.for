      MODULE DRAKCE8
        INTEGER*8 nwk8
        INTEGER*8, DIMENSION(:), ALLOCATABLE :: MAXA8
        INTEGER*4, DIMENSION(:), ALLOCATABLE :: IVRS
        INTEGER*4, DIMENSION(:), ALLOCATABLE :: ISK
        INTEGER*4, DIMENSION(:), ALLOCATABLE :: AIROWS
      END MODULE
      
      MODULE FSIDENT
        INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: IDENT
      END MODULE
      
      MODULE PLAST3D
        REAL*8, DIMENSION(:), ALLOCATABLE :: PLAST
        REAL*8, DIMENSION(:), ALLOCATABLE :: PLASTMM !za Micunov materijalni model
        LOGICAL ALLOCIRANMM
        REAL*8, DIMENSION(:), ALLOCATABLE :: PLAS1
        REAL*8, DIMENSION(:), ALLOCATABLE :: PLAS0
        INTEGER*4 iBrisiBrojac
      END MODULE
      
      MODULE ELEMENTI
        INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: VTKELEMENTI
      END MODULE
!     ovo sluzi za citanje / pisanje mesto fajla ali samo za 1 slucaj koji je pravio problem      
      MODULE WSTAZK
        LOGICAL ALLOCIRANAMATRICA
        REAL*8, DIMENSION(:), ALLOCATABLE :: RTWRITE
      END MODULE
!     ovo je novi modul za cuvanje matrice K      
      MODULE MATRICA
        LOGICAL ALLOCIRANAMATRICAK
        REAL*8, DIMENSION(:), ALLOCATABLE :: ALSK
        REAL*8, DIMENSION(:), ALLOCATABLE :: ALSM
        REAL*8, DIMENSION(:), ALLOCATABLE :: ALSC
      END MODULE
!
!     podaci o preseku za koji se crtaju rezultati za djerdap
      module PRESEKS
        integer*4 :: IPRES
        integer*4,dimension(:),allocatable :: NPRESEK
        integer*4,dimension(:,:),allocatable :: NP_ELEMENT
        integer*4,dimension(:,:),allocatable :: NP_ID
        double precision,dimension(:,:,:),allocatable :: NP_COORDS
        integer*4,dimension(:),allocatable :: NPELEM
        integer*4,dimension(:,:,:),allocatable :: NELP
        integer*4,dimension(:,:),allocatable :: NPROP
        integer*4,dimension(:,:),allocatable :: NPRO3
        double precision,dimension(:),allocatable :: UPRIS
      end module