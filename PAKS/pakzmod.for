      MODULE DRAKCE8
        INTEGER*8 nwk8
        INTEGER*8, DIMENSION(:), ALLOCATABLE :: MAXA8
        INTEGER*4, DIMENSION(:), ALLOCATABLE :: ISK
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
      END MODULE
      
      MODULE ELEMENTI
        INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: VTKELEMENTI
      END MODULE
!     ovo sluzi za citanje / pisanje mesto fajla ali samo za 1 slucaj koji je pravio problem      
      MODULE WSTAZK
        LOGICAL ALLOCIRANAMATRICA
        REAL*8, DIMENSION(:), ALLOCATABLE :: RTWRITE
      END MODULE      
      