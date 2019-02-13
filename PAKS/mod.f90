      MODULE MUMPSM
        INCLUDE 'dmumps_struc.h'
        TYPE (DMUMPS_STRUC) mumps_par
      END MODULE
      
      module ppr
        double precision,dimension(:,:),allocatable :: PORNIEL
        integer,dimension(:,:),allocatable :: IDJSTAMP1
      end module ppr
      
      module STIFFNESS
        integer*8 :: stiff_n
        integer*8 :: nonzeros
        integer*8,dimension(:),allocatable :: rows
        integer*8,dimension(:),allocatable :: columns
        double precision,dimension(:),allocatable :: stiff
      end module STIFFNESS
      
      module ELEMENTS
        double precision,dimension(:),allocatable :: thick
        integer,dimension(:),allocatable :: elemtip
        integer :: numeltip
        integer :: NMA
        integer :: NMI
        integer,dimension(6) :: eltypes
        integer,dimension(:,:),allocatable :: NEL
        integer,dimension(:),allocatable :: MCVEL
        integer,dimension(:),allocatable :: MELCV
      end module ELEMENTS
      
      module NODES
!         integer :: ISNUMER
        integer :: NPA
        integer :: NPI
        integer,dimension(:,:),allocatable :: ID
        integer,dimension(:),allocatable :: NCVEL
        integer,dimension(:),allocatable :: NELCV
        double precision, dimension(:,:),allocatable :: CORD
      end module NODES
      
      module MATRIXINIT
        integer*8,dimension(:),allocatable :: MAXA
        integer*8,dimension(:),allocatable :: MHT
      end module MATRIXINIT
      
      module PREDISCRIBED
        integer*8,dimension(:),allocatable :: NZADC
        double precision,dimension(:),allocatable :: POTZADC
        integer,dimension(:,:),allocatable :: NELTOK
        integer,dimension(:,:),allocatable :: NELR
        integer*8 :: MAXTQE
        integer*8 :: MAXER
        integer*8 :: NWATERS
        integer*8 :: NSENSORS
        integer*8 :: IBOFANG
        double precision :: Alpha
        double precision :: Rd_BOFANG
        double precision :: D_BOFANG
        integer*8,dimension(3) :: IDSENSOR
        double precision,dimension(3) :: DISTSENSOR
        integer*8,dimension(:,:),allocatable :: WATER
        integer*8,dimension(:,:),allocatable :: SENSOR
        double precision,dimension(:),allocatable :: HSENSOR
        double precision :: TPOC
        double precision,dimension(:),allocatable :: HFACE
        double precision,dimension(3) :: FSENSOR
        integer*8 :: NUMAXISPTSX
        integer*8 :: NUMAXISPTSY
        integer*8 :: NUMAXISPTSZ
        integer*8,dimension(:,:),allocatable :: INTAXISPOINTX
        integer*8,dimension(:,:),allocatable :: INTAXISPOINTY
        integer*8,dimension(:,:),allocatable :: INTAXISPOINTZ
        double precision,dimension(:),allocatable :: XAXISPTCORD
        double precision,dimension(:),allocatable :: YAXISPTCORD
        double precision,dimension(:),allocatable :: ZAXISPTCORD
        integer*8 :: IFZRAC
        integer*8 :: INFX
        integer*8 :: INFY
        integer*8 :: INFZ
        integer*8,dimension(:),allocatable :: IFLUXR
        double precision :: RKOREKCIJA
        double precision :: PREKIDNAFR
      end module PREDISCRIBED
      
      module MESURMENTPOINTS
        integer*8 :: BRKORAKA
        integer*8 :: MAX_MPOINTS
!         character*10,dimension(:),allocatable :: MPOINT_ID
!  promenjena duzina termometra za Grancarevo
        character*15,dimension(:),allocatable :: MPOINT_ID
        integer,dimension(:),allocatable :: MP_ELEMENT
        double precision,dimension(:,:),allocatable :: MP_COORDS
        double precision,dimension(:),allocatable :: MP_VREME
        double precision,dimension(:,:),allocatable :: MP_RESULTS
        double precision,dimension(:),allocatable :: MP_RESULTS_NIZ
        integer*8 :: MAX_DPOINTS
        character*10,dimension(:),allocatable :: DPOINT_ID
        integer*8,dimension(:),allocatable :: DP_ELEMENT
        double precision,dimension(:,:),allocatable :: DP_COORDS
        double precision,dimension(:,:,:),allocatable :: DP_RESULTS
      end module MESURMENTPOINTS
      
      module RESULTS
        double precision,dimension(:),allocatable :: TT1
        double precision,dimension(:),allocatable :: TT10
        double precision,dimension(:),allocatable :: PRIV
        double precision,dimension(:),allocatable :: PRIV1
        double precision,dimension(:),allocatable :: UBRZ0
        double precision,dimension(:),allocatable :: UBRZ
        double precision,dimension(:),allocatable :: BRZ0
        double precision,dimension(:),allocatable :: BRZ
        double precision,dimension(:,:),allocatable :: SKEF
        double precision,dimension(:),allocatable :: SKEFN
        double precision,dimension(:,:),allocatable :: AK
        double precision,dimension(:,:),allocatable :: DEFOR
        double precision,dimension(:,:,:),allocatable :: VG
        double precision,dimension(:,:,:),allocatable :: GG
        double precision,dimension(:,:),allocatable :: VECTJ
        double precision,dimension(:,:),allocatable :: POMER
        integer*8,dimension(:),allocatable :: IVECT
        double precision,dimension(:),allocatable :: SILE
      end module RESULTS
      
      module KONTURE
        integer*8,dimension(:),allocatable :: LIN
        integer,dimension(:),allocatable :: LINSN
        double precision,dimension(:),allocatable :: QUK
        double precision,dimension(:),allocatable :: QUM
      end module KONTURE
 
      module pflux
        double precision,dimension(:),allocatable :: PFLUXEL
        double precision,dimension(:),allocatable :: FCONEL
        double precision,dimension(:),allocatable :: TOKOLINEL
      end module pflux
      
      module PRESEK
        integer*4,dimension(:),allocatable :: NPRESEK
        integer*4 :: IPRES
        integer*4,dimension(:,:),allocatable :: NP_ELEMENT
        integer*4,dimension(:,:),allocatable :: NP_ID
        double precision,dimension(:,:,:),allocatable :: NP_COORDS
        integer*4,dimension(:),allocatable :: NPELEM
        integer*4,dimension(:,:,:),allocatable :: NELP
        integer*4,dimension(:,:),allocatable :: NPROP
!         integer :: NUMMAT
!         double precision,dimension(:,:),allocatable :: EPRESEK
!         integer*8,dimension(:),allocatable :: EL_TYPE
      end module PRESEK
      module zapresil
        double precision,dimension(:),allocatable :: SILEL
      end module zapresil
      module Crtanje
        double precision,dimension(:),allocatable :: TTemp
    end module Crtanje
    
    module TEMPCVOROVI
        LOGICAL ALLOCIRANEtemperature
        double precision,dimension(:),allocatable :: TEMPuCVORU
    end module TEMPCVOROVI
