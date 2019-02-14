C$DEBUG 
      ! AKO JE EXE ONDA JE PROGRAM AKO JE LIB ONDA JE SUBROUTINE
      !PROGRAM GLAVNI
      SUBROUTINE VBRISI
      WRITE (*,*) "BRISI"
      RETURN
      END
      SUBROUTINE VPAKV
      use mcm_database
      USE NODES
      USE ELEMENTS
      USE MATRIXINIT
      USE PREDISCRIBED
      USE RESULTS
      USE MESURMENTPOINTS
      USE KONTURE
      USE PRESEK
      USE TEMPCVOROVI
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'paka.inc'
      INCLUDE 'mpif.h'
!      PARAMETER (NTOT = 200000000)
      COMMON /IME/ IME
      COMMON /SRPSKI/ ISRPS
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /CDEBUG/ IDEBUG
      COMMON /RADNIV/ MAXVEC
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /JEDANP/ INDJED,NBRF,NGL,INDTLO
      COMMON /TACNOS/ EPSTR,MAXIT,NJRAP
      COMMON /PROMEN/ NJUTN,INDOPT,INDPRO
      COMMON /PENALL/ PENALT,PRESS,IBRGT
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /NDESUK/ NDES,IDPRIT,IFORM
      COMMON /TIPEL/ NP2DMX
      COMMON /BROJUK/ INDFOR,NULAZ
      COMMON /PRIMER/ IPRBR,INDIZL,INDGRA
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /POCETN/ IPOCU,IPOCV,IPOCP,IPOCT,POCU,POCV,POCP,POCT,GAMA
      COMMON /VREPER/ NPER,NTABFT
      COMMON /REPERM/ MREPER(4)
      COMMON /AXIS/ HH,RIZ,INDAX
      COMMON /DODAT/ NDIMM
      COMMON /DYNAM/ NDIN
      COMMON /ELECTR/ INDJOT,INDFS,INDJOP
      COMMON /STAC/ NSTAC,NSLP,IBROJ
      COMMON /MATER/ NUMMAT
      COMMON /VOPTIM/ NKONT,MAXLIN
      COMMON /VDP/ DT,NKORP,NN,NZAV
      COMMON /CERNE/ IGRESK,IPROTK,IZSILE
      COMMON /DUPLAP/ IDVA
      COMMON /FLUX/ KISA
      COMMON /ZADPO/ LNZADJ
      COMMON /OSATEZ/ ZASIC,IOSA,ISIL,IFIL,IPOR
      COMMON /SREDIN/ IMASS
C
      COMMON /PRIKAZ/ INDSC,IZPOT
      COMMON /PIJEZO/ CPOR(3,1000),NPIJEZ(20,100),NODP(100),NPIJ,NPOR,
     1                NPORCV(1000),NPOREL(1000),NEPOR,NPORE,NEPIJ
      COMMON /ICITANJE/ INPT
      COMMON /SOLVER/ ISOLVER
      COMMON /NXNAST/ NXNASTRAN
      COMMON /DJERDAP/ IDJERDAP,ISPRESEK
      COMMON /INDBRANA/ IBRANA
      COMMON /STAMPAZT/ NPRINT
C
C    DEFINISANJE MAKSIMALNE RADNE MEMORIJE U VEKTORU A
C
!      COMMON A(NTOT)
!      REAL A

      DIMENSION IA(1)
      EQUIVALENCE(A(1),IA(1))
C
C
C
C
      CHARACTER*80 NASLOV
      CHARACTER*80 KRAJ
      CHARACTER*250 ACOZ
      CHARACTER*50 IME
      CHARACTER*54 IMECSV
      CHARACTER*20 IMEPR
      CHARACTER*2 IMECSV2
      CHARACTER*1 IMECSV1,PREF
      CHARACTER*53 IZIPOCT
      
C
      CHARACTER*5 PRTOK(10)
      


      DIMENSION LM2(36)
C      DIMENSION IBFK(5,3),FAKP(5,3),TOPM(5),TMNM(5)
      DATA PRTOK /'FLUX1','FLUX2','FLUX3','FLUX4','FLUX5','FLUX6',
     *            'FLUX7','FLUX8','FLUX9','FLU10'/
      integer ierr, myid
      integer*8 LM2,INDWATER
      integer Dtime(8),ISNUMER
      INTEGER*4 IBRANA
      ALLOCIRANEtemperature = .FALSE.
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      !return
      IF (myid.ne.0) goto 1234
       MTOT=NTOTAL
       MAXVEC=NTOTAL
C
C IDVA - INDIKATOR ZA NUMERICKU TACNOST
C IDVA=1 JEDNOSTRUKA TACNOST
C IDVA=2 DVOSTRUKA TACNOST
C
      IDVA=2
C      IDODAV=1
      IDODAV=0
C       ISRPS=0
      ISRPS=1
C      IULAZ=5
      IULAZ=1
C KOD NDP-a JE IIZLAZ=3 A IULAZ=1
      IIZLAZ=3
C      IIZLAZ=4
      IDEBUG=0
      INDFOR=2
      NBLOCK=1
      INDIZL=0
      INDGRA=0
      IFORM=0
      IPRBR=0
      INDJOT=0
      NGE=1
C
      IGRESK=35
      IPROTK=21
      IZSILE=25
C
C=======================================================================
      CALL VOTVORI
      CALL DATE_AND_TIME(VALUES=Dtime)
      WRITE(*,*) 'posle otvaranja dat fajla', (Dtime(i),i=5,7)
      CALL VISPITA(ACOZ)                                       
      READ(ACOZ,1000) NASLOV
 1000 FORMAT(A80)
C=======================================================================
      CALL VISPITA(ACOZ)                                       
      READ(ACOZ,1001) INDFOR,INPT,NXNASTRAN,IPAKT,IDJERDAP,ISPRESEK,
     * IBRANA
 1001 FORMAT(14I5)
      write (*,*) 'IBRANA',IBRANA
C=======================================================================
      CALL VISPITA(ACOZ)                                       
      IF(INPT.EQ.1) THEN
       READ(ACOZ,1022) NPT,NGET,NMATT,NDIN,NPER,NPRINT,NSTAC,NULAZ,
     1       ISNUMER
      ELSE
       READ(ACOZ,1002) NPT,NGET,NMATT,NDIN,NPER,NPRINT,NSTAC,NULAZ,
     1       ISNUMER
      ENDIF
      if(NSTAC.EQ.2) THEN
        isave=1
      ELSE
        isave=1
      ENDIF
!       write(*,*)"NSTAC,isave",NSTAC,isave
 1002 FORMAT(16I5)
 1022 FORMAT(I10,15I5)
      IF(NPER.EQ.0) NPER = 1
      IF(NGET.EQ.0) NGET = 1
      IF(NMATT.EQ.0) NMATT=1
      IF(NPRINT.EQ.0) NPRINT=1
      IF(NULAZ.EQ.0) NULAZ = 1
      CALL VOTVIZL
      CALL VZAGLAV
      CALL VOTVGRA
      IF(IPAKT.EQ.0) CALL VOTVSIL
      WRITE(3,*) 'posle otvaranja dat fajla', (Dtime(i),i=5,7)
!  otvara se u pakv3d      
!       IF(IPAKT.EQ.1) CALL VOTVTEMP
      CALL VOTVNEU
      WRITE(IIZLAZ,3019) NASLOV
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,3020) INDFOR,INPT,NXNASTRAN,IPAKT,IDJERDAP,ISPRESEK,
     * IBRANA
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,7020) INDFOR,INPT,NXNASTRAN,IPAKT,IDJERDAP,ISPRESEK,
     * IBRANA
 3019 FORMAT(//78('*')/,A80,/78('*')///)

 3020 FORMAT(//
     111X,'NACIN UCITAVANJA ULAZNIH PODATAKA ............. INDFOR =',I5/
     116X,'EQ.0; INDFOR = 1'/
     116X,'EQ.1; U SLOBODNOM FORMATU'/
     116X,'EQ.2; U OPISANOM FORMATU'///
     111X,'BROJ CVOROVA VECI OD 99999 ....................   INPT =',I5/
     116X,'EQ.0; NE'/
     116X,'EQ.1; DA'///
     111X,'BROJ CITANJE POTENCIJALA SA DISKA .......... NXNASTRAN =',I5/
     116X,'EQ.0; NE'/
     116X,'EQ.1; DA'///
     111X,'IZVRSAVANJE PROGRAMA PAKP/PAKT ................. IPAKT =',I5/
     116X,'EQ.0; PAKP'/
     116X,'EQ.1; PAKT'///
     111X,'STAMPANJE ZAP.SILA/TEMP U CVOROVIMA LAMELE...IDJERDAP =',I5/
     116X,'EQ.0; STAMPA U SVIM CVOROVIMA'/
     116X,'GE.1; STAMPA U CVOROVIMA ZA LAMELU'/
     111X,'STAMPANJE REZULTATA U PRESECIMA .............ISPRESEK =',I5/
     116X,'EQ.1; stampa se presek za lamelu IDJERDAP'/
     116X,'EQ.0; NE STAMPAJU SE'/
     116X,'EQ.-1; STAMPAUJU SE ZA SVE PRESEKE za Djerdap'/
     111X,'Indikator brane za interpolaciju temp ..........IBRANA =',I5/
     116X,'EQ.0; Djerdap'/
     116X,'EQ.1; Grancarevo'///)
 7020 FORMAT(//
     111X,'FORMAT INDIKATOR .............................. INDFOR =',I5/
     116X,'EQ.0; INDFOR = 1'/
     116X,'EQ.1; FREE FORMAT'/
     116X,'EQ.2; FIXED FORMAT, AS DESCRIBED'///
     111X,'Num. NODE GREAT FROM 99999 ....................   INPT =',I5/
     116X,'EQ.0; NO'/
     116X,'EQ.1; YES'///
     111X,'BROJ CITANJE POTENCIJALA SA DISKA .......... NXNASTRAN =',I5/
     116X,'EQ.0; NO'/
     116X,'EQ.1; YES'///
     111X,'START PROGRAM PAKP/PAKT ........................ IPAKT =',I5/
     116X,'EQ.0; PAKP'/
     116X,'EQ.1; PAKT'///
     111X,'WRITE FIL.FORCE/TEMP IN NODES FOR LAMELE ....IDJERDAP =',I5/
     116X,'EQ.0; IN ALL NODES'/
     116X,'GE.1; ONLY FOR LAMELE'/
     111X,'WRITE RESULTS IN CROSS SECTION ..............ISPRESEK =',I5/
     116X,'EQ.0; NO'/
     116X,'EQ.-1; all cross section for dam Djerdap'/
     111X,'INDICATOR for interpolation temp     ..........IBRANA =',I5/
     116X,'EQ.0; Djerdap'/
     116X,'EQ.1; Grancarevo'///)
C
      IF(ISRPS.EQ.0)
C     *WRITE(IIZLAZ,2001) NPT,NGET,NMATT,NDIN,NPER,NPRINT,NSTAC,IANIZ
     *WRITE(IIZLAZ,2001) NPT,NPER,NPRINT,NSTAC,ISNUMER
      IF(ISRPS.EQ.1)
C     *WRITE(IIZLAZ,6001) NPT,NGET,NMATT,NDIN,NPER,NPRINT,NSTAC,IANIZ
     *WRITE(IIZLAZ,2001) NPT,NPER,NPRINT,NSTAC,ISNUMER
 2001 FORMAT(6X,'O S N O V N I    P O D A C I    O    P R O B L E M U'
     1/6X,51('-')///
     111X,'UKUPAN BROJ CVORNIH TACAKA ....................... NP =',I10/
     116X,'EQ.0; PREKIDA SE IZVRSAVANJE PROGRAMA'///
C     211X,'BROJ GRUPA ELEMENATA ............................ NGET =',I5/
C     216X,'EQ.0; POSTAJE "1"; (MAX. 10 GRUPA)'///
C     311X,'BROJ RAZLICITIH MATERIJALA ..................... NMATT =',I5/
C     316X,'EQ.0; POSTAJE "1"'///
C     411X,'INDIKATOR DINAMIKE .............................. NDIN =',I5/
C     416X,'EQ.0; STATIKA '/
C     416X,'EQ.1; DINAMIKA'///
     511X,'BROJ PERIODA SA KONSTANTNIM VREMENSKIM KORACIMA . NPER =',I5/
     516X,'EQ.0; POSTAJE "1"; (MAX. 16 PERIODA)'///
     611X,'DEFINISANJE STAMPARSKOG KORAKA ................ NPRINT =',I5/
     616X,'EQ.0; POSTAJE "1"'///
     411X,'INDIKATOR STACIONARNOSTI ....................... NSTAC =',I5/
     416X,'EQ.0; NESTACIONARNA ANALIZA'/
     416X,'EQ.1; STACIONARNA ANALIZA'/
     416X,'EQ.2; PRVI KORAK STACIONARAN A DRUGI NESTACIONARNI'///
     411X,'INDIKATOR SLOBODNE NUMERACIJE ................ ISNUMER =',I5/
     416X,'EQ.0; NE'/
     416X,'EQ.1; DA')
 6001 FORMAT(6X,'B A S I C    D A T A    F O R   T H E   P R O B L E M'
     1/6X,53('-')///
     111X,'TOTAL NUMBER OF NODAL POINTS ..................... NP =',I10/
     116X,'EQ.0; PROGRAM STOP'///
C     211X,'NUMBER OF ELEMENT GROUPS ........................ NGET =',I5/
C     216X,'EQ.0; DEFAULT SET "1"'///
C     311X,'NUMBER OF DIFFERENT MATERIALS .................. NMATT =',I5/
C     316X,'EQ.0; DEFAULT SET "1"'///
C     411X,'DYNAMICS INDICATOR .............................. NDIN =',I5/
C     416X,'EQ.0; STATIC'/
C     416X,'EQ.1; DYNAMIC'///
     511X,'NUMBER OF CONSTANT TIME STEP PERIODS ............ NPER =',I5/
     516X,'EQ.0; DEFAULT SET "1"; (MAX. 16 PERIODS)'///
     611X,'OUTPUT PRINTING INTERVAL ...................... NPRINT =',I5/
     616X,'EQ.0; DEFAULT SET "1"'///
     411X,'INDICATOR OF STATIONARITY ...................... NSTAC =',I5/
     416X,'EQ.0; TRANSIENT ANALYSIS'/
     416X,'EQ.1; STATIONARY ANALYSIS'/
     416X,'EQ.2; FIRST STEP-STATIONARY AND OTHER TRANSIENT'///
     411X,'INDICATOR FOR FREE NUMERATION ................ ISNUMER =',I5/
     416X,'EQ.0; NO'/
     416X,'EQ.1; YES')

C==========================================================================
      CALL VISPITA(ACOZ)
      READ(ACOZ,1003) INTEB,INDSC,IFORM,MAXIT,EPSTA,EPSTR,NJRAP
 1003 FORMAT(4I5,2F10.2,I5)
      IF(INTEB.EQ.0) INTEB=1
      IF(NJRAP.EQ.0) NJRAP=1
      IF(MAXIT.EQ.0) MAXIT=15
      IF(DABS(EPSTA).LT.1.D-10.AND.DABS(EPSTR).LT.1.D-10) EPSTR=.001
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2002) INTEB,INDSC,IFORM,MAXIT,EPSTA,EPSTR,NJRAP
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6002) INTEB,INDSC,IFORM,MAXIT,EPSTA,EPSTR,NJRAP
 2002 FORMAT(///
     111X,'PRIMENJEN METOD VREMENSKE INTEGRACIJE .......... INTEB =',I5/
     116X,'EQ.0; POSTAJE "1"'/
     116X,'EQ.1; EULER BACWARD INTEGRACIJA (ALFA = 1)'/
     116X,'EQ.2; NE KORISTI SE'///
     211X,'STAMPANJE REZULTATA U ZELJENIM CVOROVIMA ....... INDSC =',I5/
     216X,'EQ.0; U SVIM CVOROVIMA'/
     216X,'EQ.1; U ZELJENIM CVOROVIMA'///
     311X,'STAMPANJE REZULTATA U ZELJENOM FORMATU ......... IFORM =',I5/
     316X,'EQ.0; U FORMATU  D13.5'/
     316X,'EQ.1; U FORMATU  F10.3'///
     411X,'MAXIMALNI BROJ RAVNOTEZNIH ITERACIJA ........... MAXIT =',I5/
     416X,'EQ.0; POSTAJE "15"'///11X,
     5'APSOLUTNA TACNOST PRI ITERACIJAMA ......... EPSTA =',1PD10.3///
     611X,'RELATIVNA TACNOST PRI ITERACIJAMA ......... EPSTR =',1PD10.3/
     616X,'EQ.0; POSTAJE "1.E-3"'///
     111X,'PRIMENJEN ITERATIVNI METOD ..................... NJRAP =',I5/
     116X,'EQ.0; POSTAJE "1"'/
     116X,'EQ.1; MODIFIKOVAN NJUTN-RAPSONOV METOD'/
     116X,'EQ.2; POTPUN NJUTNOV')
 6002 FORMAT(///
     111X,'TIME INTEGRATION METHOD USED ................... INTEB =',I5/
     116X,'EQ.0; DEFAULT SET "1"'/
     116X,'EQ.1; EULER BACKWARD INTEGRATION (ALPFA = 1)'///
     211X,'PRINT OF RESULTS IN PRESCRIBED NODES ........... INDSC =',I5/
     216X,'EQ.0; IN ALL NODES'/
     216X,'EQ.1; IN PRESCRIBED NODES'///
     311X,'PRINT OF RESULTS IN PRESCRIBED FORMAT .......... IFORM =',I5/
     316X,'EQ.0; IN FORMAT  D13.5'/
     316X,'EQ.1; IN FORMAT  F10.3'///
     411X,'MAXIMUM NUMBER OF ITERATIONS ................... MAXIT =',I5/
     416X,'EQ.0; DEFAULT SET "15"'///11X,
     5'ABSOLUTE ACCURACY AT ITERATIONS ........... EPSTA =',1PD10.3///
     611X,'RELATIVE ACCURACY AT ITERATIONS ........... EPSTR =',1PD10.3/
     616X,'EQ.0; DEFAULT SET "1.E-3"'///
     111X,'ITERATION METHODS EMPLOYED ..................... NJRAP =',I5/
     116X,'EQ.0; DEFAULT SET "1"'/
     116X,'EQ.1; MODIFIED NEWTON'/
     116X,'EQ.2; FULL NEWTON')
C==========================================================================
      CALL VISPITA(ACOZ)
      READ(ACOZ,1004) IREST,ICCGG
      isolver=iccgg
      IF(isolver.eq.-11) NBLOCK=1
 1004 FORMAT(2I5)
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2004)IREST,ICCGG
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6004)IREST,ICCGG
 2004 FORMAT(///
     111X,'INDIKATOR IZVRSENJA PROGRAMA ................... IREST =',I5/
     116X,'EQ.0; IZVRSENJE PROGRAMA'/
     116X,'EQ.1; KONTROLA ULAZNIH PODATAKA'/
     111X,'VRSTA SOLVERA ................................ ICCGG =',I5//)
 6004 FORMAT(///
     111X,'INDICATOR FOR JOB EXECUTION .................... IREST =',I5/
     116X,'EQ.0; JOB EXECUTION'/
     116X,'EQ.1; CHECK OF INPUT DATA'/
     111X,'VRSTA SOLVERA ................................ ICCGG =',I5//)
C==========================================================================
      LVREME=1
      CALL VPERIOD(A(LVREME))
      LMAX=LVREME+NPER*950*IDVA
C==========================================================================
!       LID=LMAX
!       LMAX=LID+1*NPT
!       CALL VDELJIV(LMAX,2,INDL)
!       IF(INDL.EQ.0) LMAX=LMAX+1
!       LCORD=LMAX
!       LMAX=LCORD+3*NPT*IDVA
!       CALL VPROMEM(LMAX)

      if (.not.allocated(ID)) allocate(ID(1,NPT),STAT=istat)
      if (.not.allocated(CORD)) allocate(CORD(3,NPT),STAT=istat) 
      if(ISNUMER.EQ.1) then
        if (.not.allocated(NCVEL)) allocate(NCVEL(NPT),STAT=istat)
      endif
C UCITAVANJE ID(6,*) MATRICE I KOORDINATE CVOROVA-MATRICA CORD(3,*)
Cxxx
C      CALL VCZINIT
C      CALL VCZPROV(NPT,1,34)
Cxxx
!       CALL VULAZT1(ID,CORD,NPT,JEDN)
      CALL VULAZT1(ISNUMER)
      if(ISNUMER.EQ.1) then
        NPM=NPA-NPI+1
        if (.not.allocated(NELCV)) allocate(NELCV(NPM),STAT=istat)
!         write(*,*) 'NPA,NPI', NPA,NPI
        CALL VVEZACV(NPT)
      endif
Cxxx
C      CALL VZASTIT
Cxxx
C==========================================================================
      CALL VISPITA(ACOZ)
      IF(INPT.EQ.1) THEN
       READ(ACOZ,1018) NETIP,NET,INDAX,HH,RIZ,INDTLO
      ELSE
       READ(ACOZ,1008) NETIP,NET,INDAX,HH,RIZ,INDTLO
      ENDIF
 1008 FORMAT(3I5,2F10.3,I5) 
 1018 FORMAT(I5,I10,I5,2F10.3,I5) 
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2011) NGET,NETIP,NGET,NET,NGET,INDAX,HH,RIZ
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6011) NGET,NETIP,NGET,NET,NGET,INDAX,HH,RIZ
 2011 FORMAT(////
     111X,'TIP KONACNOG ELEMENTA GRUPE ELEMENATA ',I6,' ... NETIP =',I5/
     116X,'EQ.2; IZOPARAM. 2D ELEMENT U RAVNI X-Y'/
     116X,'EQ.3; PROSTORNI IZOPARAMETARSKI 3D ELEMENT'///
     211X,'BROJ KONACNIH ELEMENATA GRUPE ELEMENATA ',I6,' .. NET =',I10/
     216X,'GT.0; '///
     311X,'INDIKATOR ZA OSNOSIMETRICNE ELEMENTE',I6,' ..... INDAX =',I5/
     316X,'(SAMO ZA 2D ELEMENTE)'/
     316X,'EQ.0; 2D U RAVNI X-Y'/
     316X,'EQ.1; 2D OSNOSIMETRICNI, Y JE OSA SIMETRIJE'///11X,
     4'VISINA POVRSINE SA ZADATIM FLUKSOM ....... HEIGHT =',1PD12.5///11
     4X,'RASTOJANJE OD OSE SIMETRIJE .............. RADIUS =',1PD12.5//)
 6011 FORMAT(////
     111X,'ELEMENT TYPE FOR GROUP ',I6,' .................. NETIP =',I5/
     116X,'EQ.2; ISOPARAMETRIC 2D ELEMENT IN PLANE X-Y'/
     116X,'EQ.3; ISOPARAMETRIC 3D ELEMENT'///
     211X,'NUMBER OF ELEMENTS IN THE GROUP ',I6,' .......... NET =',I10/
     216X,'GT.0; '///
     311X,'INDICATOR FOR AXISYMMETRIC ELEMENTS',I6,' ...... INDAX =',I5/
     316X,'(ONLY FOR 2D ELEMENTS)'/
     316X,'EQ.0; 2D IN PLANE X-Y'/
     316X,'EQ.1; 2D AXISYMMETRIC, Y IS AXIS OF SYMMETRY '///11X,
     4'HEIGH OF SURFACE WITH PRESCRIBED FLUX .... HEIGHT =',1PD12.5///11
     4X,'DISTANCE FROM AXIS OF SYMMETRY ........... RADIUS =',1PD12.5//)
C
C==========================================================================
      PENALT=0.
      INDJOT=0
      CALL VISPITA(ACOZ)
C      READ(ACOZ,1006) NMAT2D,MAT2D,NP2DMX,PENALT
      READ(ACOZ,1006) NMAT2D,MAT2D,NP2DMX,INDJOT,INDJOP,IMASS
      IMASS=1
!       write(*,*)'INDJOT',INDJOT
C 1006 FORMAT(3I5,D10.3)
 1006 FORMAT(7I5)
C      IF(NMAT2D.GT.1) MAT2D=0
C      IF(NMATT.EQ.1) MAT2D=1
CC      IF(NMATT.LT.NMAT2D) STOP '(NMATT.LT.NMAT2D)'
      IF(NMATT.EQ.1) NMAT2D=1
      IF(NP2DMX.EQ.0) NP2DMX=20
      IF(NETIP.EQ.1) THEN
      IF(ISRPS.EQ.0)
C     *WRITE(IIZLAZ,3022) NMAT2D,NP2DMX,PENALT
     *WRITE(IIZLAZ,3021) NP2DMX,INDJOT,INDJOP,INDJOS,IMASS
      IF(ISRPS.EQ.1)
C     *WRITE(IIZLAZ,6022) NMAT2D,NP2DMX,PENALT
     *WRITE(IIZLAZ,6021) NP2DMX,INDJOT,INDJOP,INDJOS,IMASS
      ENDIF
      IF(NETIP.EQ.2) THEN
      IF(ISRPS.EQ.0)
C     *WRITE(IIZLAZ,3022) NMAT2D,NP2DMX,PENALT
     *WRITE(IIZLAZ,3022) NP2DMX,INDJOT,INDJOP,INDJOS,IMASS
      IF(ISRPS.EQ.1)
C     *WRITE(IIZLAZ,6022) NMAT2D,NP2DMX,PENALT
     *WRITE(IIZLAZ,6022) NP2DMX,INDJOT,INDJOP,INDJOS,IMASS
      ENDIF
      IF(NETIP.EQ.3) THEN
      IF(ISRPS.EQ.0)
C     *WRITE(IIZLAZ,3032) NMAT2D,NP2DMX,PENALT
     *WRITE(IIZLAZ,3032) NP2DMX,INDJOT,INDJOP,INDJOS,IMASS
      IF(ISRPS.EQ.1)
C     *WRITE(IIZLAZ,6032) NMAT2D,NP2DMX,PENALT
     *WRITE(IIZLAZ,6032) NP2DMX,INDJOT,INDJOP,INDJOS,IMASS
      ENDIF
 3021 FORMAT(//
     111X,'ULAZNI PODACI ZA 1D IZOPARAMETARSKI ELEMENT '//
C     111X,'UKUPAN BROJ RAZLICITIH MATERIJALA ............. NMAT2D =',I5/
C     116X,'(ZA NMAT2D.EQ.1 IGNORISE SE UNETI PODATAK KOD ELEMENTA)'/
C     116X,'EQ.1;(ZA NMATT.EQ.1)'/
C     116X,'GT.1; ZA NMATT.GT.1'///
     311X,'MAKSIMALAN BROJ CVOROVA PO ELEMENTU ........... NP2DMX =',I5/
     316X,'EQ.0; POSTAJE "2" (ZA 2 CVORA)'/
     316X,'EQ.3; ZA ELEMENTE OD 3 CVORA'//
     311X,'Broj gausovih tacaka po zapremini ........... INDJOT =',I5/
     311X,'Broj gausovih tacaka po povrsini ............ INDJOP =',I5/
     311X,'Osrednjavanje ............................... INDJOS =',I5/
     311X,'Indikaktor matrice C ........................ IMASS =',I5/
C     311X,'PENALTI FAKTOR ................................ PENALT =',
C     31PD10.3//
     +)
 3022 FORMAT(//
     111X,'ULAZNI PODACI ZA 2D IZOPARAMETARSKI ELEMENT '//
C     111X,'UKUPAN BROJ RAZLICITIH MATERIJALA ............. NMAT2D =',I5/
C     116X,'(ZA NMAT2D.EQ.1 IGNORISE SE UNETI PODATAK KOD ELEMENTA)'/
C     116X,'EQ.1;(ZA NMATT.EQ.1)'/
C     116X,'GT.1; ZA NMATT.GT.1'///
     311X,'MAKSIMALAN BROJ CVOROVA PO ELEMENTU ........... NP2DMX =',I5/
     316X,'EQ.0; POSTAJE "4" (ZA 4 CVORA)'/
     316X,'GT.4; ZA ELEMENTE OD 5 DO 9 CVOROVA'//
     311X,'Broj gausovih tacaka po zapremini ........... INDJOT =',I5/
     311X,'Broj gausovih tacaka po povrsini ............ INDJOP =',I5/
     311X,'Osrednjavanje ............................... INDJOS =',I5/
     311X,'Indikaktor matrice C ........................ IMASS =',I5/
C     311X,'PENALTI FAKTOR ................................ PENALT =',
C     31PD10.3//
     +)
 3032 FORMAT(//
     111X,'ULAZNI PODACI ZA 3D IZOPARAMETARSKI ELEMENT '//
C     111X,'UKUPAN BROJ RAZLICITIH MATERIJALA ............. NMAT3D =',I5/
C     116X,'(ZA NMAT3D.EQ.1 IGNORISE SE UNETI PODATAK KOD ELEMENTA)'/
C     116X,'EQ.1;(ZA NMATT.EQ.1)'/
C     116X,'GT.1; ZA NMATT.GT.1'///
     311X,'MAKSIMALAN BROJ CVOROVA PO ELEMENTU ........... NP3DMX =',I5/
     316X,'EQ.0; POSTAJE "8" (ZA 8 CVORA)'/
     316X,'GT.8; ZA ELEMENTE OD 9 DO 20 CVOROVA'//
     311X,'Broj gausovih tacaka po zapremini ........... INDJOT =',I5/
     311X,'Broj gausovih tacaka po povrsini ............ INDJOP =',I5/
     311X,'Osrednjavanje ............................... INDJOS =',I5/
     311X,'Indikaktor matrice C ........................ IMASS =',I5/
C     311X,'PENALTI FAKTOR ................................ PENALT =',
C     31PD10.3//
     +)
 6021 FORMAT(//
     111X,'INPUT DATA FOR 1D ISOPARAMETRIC ELEMENT'//
C     111X,'TOTAL NUMBER DIFFERENT MATERIALS  ............. NMAT2D =',I5/
C     116X,'EQ.1;(ZA NMATT.EQ.1)'/
C     116X,'GT.1; ZA NMATT.GT.1'///
     311X,'MAKSIMAL NUMBER OF NODES PER ELEMENT .......... NP2DMX =',I5/
     316X,'EQ.0; BECOME "2" (FOR 2 NODES)'/
     316X,'EQ.3; FOR ELEMENT OF 3 NODES'//
     311X,'Broj gausovih tacaka po zapremini ........... INDJOT =',I5/
     311X,'Broj gausovih tacaka po povrsini ............ INDJOP =',I5/
     311X,'Osrednjavanje ............................... INDJOS =',I5/
     311X,'Indikaktor matrice C ........................ IMASS =',I5/
C     311X,'PENALTY FACTOR ................................ PENALT =',
C     31PD10.3//
     +)
 6022 FORMAT(//
     111X,'INPUT DATA FOR 2D ISOPARAMETRIC ELEMENT'//
C     111X,'TOTAL NUMBER DIFFERENT MATERIALS  ............. NMAT2D =',I5/
C     116X,'EQ.1;(ZA NMATT.EQ.1)'/
C     116X,'GT.1; ZA NMATT.GT.1'///
     311X,'MAKSIMAL NUMBER OF NODES PER ELEMENT .......... NP2DMX =',I5/
     316X,'EQ.0; BECOME "4" (FOR 4 NODES)'/
     316X,'GT.4; FOR ELEMENT OF 5 TO 9 NODES'//
     311X,'Broj gausovih tacaka po zapremini ........... INDJOT =',I5/
     311X,'Broj gausovih tacaka po povrsini ............ INDJOP =',I5/
     311X,'Osrednjavanje ............................... INDJOS =',I5/
     311X,'Indikaktor matrice C ........................ IMASS =',I5/
C     311X,'PENALTY FACTOR ................................ PENALT =',
C     31PD10.3//
     +)
 6032 FORMAT(//
     111X,'INPUT DATA FOR 3D ISOPARAMETRIC ELEMENT'//
C     111X,'TOTAL NUMBER DIFFERENT MATERIALS  ............. NMAT3D =',I5/
C     116X,'EQ.1;(ZA NMATT.EQ.1)'/
C     116X,'GT.1; ZA NMATT.GT.1'///
     311X,'MAKSIMAL NUMBER OF NODES PER ELEMENT .......... NP3DMX =',I5/
     316X,'EQ.0; BECOME "8" (FOR 8 NODES)'/
     316X,'GT.8; FOR ELEMENT OF 9 TO 20 NODES'//
     311X,'Broj gausovih tacaka po zapremini ........... INDJOT =',I5/
     311X,'Broj gausovih tacaka po povrsini ............ INDJOP =',I5/
     311X,'Osrednjavanje ............................... INDJOS =',I5/
     311X,'Indikaktor matrice C ........................ IMASS =',I5/
C     311X,'PENALTY FACTOR ................................ PENALT =',
C     31PD10.3//
     +)
C=========================================================================
C DODELJIVANJE ZA NDIM
!       NDIM=NP2DMX
      NDIM=20
      NDIMM=NDIM+2
      NDES=NDIM
C
C==========================================================================
!       LNEL=LMAX
!       LMAX=LNEL+NDIMM*NET
!       CALL VPROMEM(LMAX)
C UCITAVANJE ELEMENATA I SMESTANJE U MATRICU NEL(9,*)
!       CALL VICLEAR(NEL,NDIMM*NET)
      if (.not.allocated(NEL)) allocate(NEL(NDIMM,NET),STAT=istat) 
      if(ISNUMER.EQ.1) then
        if (.not.allocated(MCVEL)) allocate(MCVEL(NET),STAT=istat)
      endif
    
      CALL VULAZR3(ISNUMER)
      IF(ISNUMER.EQ.1) THEN
!       write(*,*) "NMA,NMI", NMA,NMI
!        write(*,*) " EL(1,2)",NEL(1,2)
CE    CONNECT NODE NUMBERS FOR FREE NODE NUMERATION
        CALL VVEZACE(NET,NDIM)
CE      NMA IS MAXIMUM AND NMI IS MINIMUM ELEMENT NUMBERS USED IN FREE
CE      ELEMENT NUMERATION
        LLMEL=NMA-NMI+1
        if (.not.allocated(MELCV)) allocate(MELCV(LLMEL),STAT=istat)
CE      CONNECT ELEMENT NUMBERS FOR FREE ELEMENT NUMERATION
        CALL VVEZAEL(NET)
      ENDIF
Cxxx
C      CALL VCZPROV(NET,1,34)
Cxxx      
      CALL VISPITA(ACOZ)
      IF(INPT.EQ.1)THEN
        READ(ACOZ,1111) NUMZAD,NUMAXISPTSX,NUMAXISPTSY,NUMAXISPTSZ
      ELSE
        READ(ACOZ,1001) NUMZAD,NUMAXISPTSX,NUMAXISPTSY,NUMAXISPTSZ
      ENDIF
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,3007) NUMZAD,NUMAXISPTSX,NUMAXISPTSY,NUMAXISPTSZ
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6007) NUMZAD,NUMAXISPTSX,NUMAXISPTSY,NUMAXISPTSZ
 1111 FORMAT(4I10)
 3007 FORMAT(//
     111X,'UKUPAN BROJ ZADATIH VREDNOSTI ZA POTENCIJAL NUMZAD=  ',I10//
     211X,'BROJ TACAKA NA OSI X ZA INTERPOLACIJU...NUMAXISPTSX=',I10//
     311X,'BROJ TACAKA NA OSI Y ZA INTERPOLACIJU...NUMAXISPTSY=',I10//
     411X,'BROJ TACAKA NA OSI Z ZA INTERPOLACIJU...NUMAXISPTSZ=',I10//)
 6007 FORMAT(//
     111X,'NUMBER OF PRESCRIBED VALUES FOR POTENTIAL ...NUMZAD= ',I10//
     211X,'NUMBER OF POINTS ON X INTERPOLATION AXIS NUMAXISPTSX=',I10//
     311X,'NUMBER OF POINTS ON Y INTERPOLATION AXIS NUMAXISPTSY=',I10//
     411X,'NUMBER OF POINTS ON Z INTERPOLATION AXIS NUMAXISPTSZ=',I10//)
C==========================================================================
      IF(NUMAXISPTSX.GT.0) THEN
        if (.not.allocated(INTAXISPOINTX)) allocate(INTAXISPOINTX
     1                              (2,NUMAXISPTSX),STAT=istat)
        if (.not.allocated(XAXISPTCORD)) allocate(XAXISPTCORD
     1                                  (NUMAXISPTSX),STAT=istat)
      ENDIF
      IF(NUMAXISPTSY.GT.0) THEN
        if (.not.allocated(INTAXISPOINTY)) allocate(INTAXISPOINTY
     1                              (2,NUMAXISPTSY),STAT=istat)
        if (.not.allocated(YAXISPTCORD)) allocate(YAXISPTCORD
     1                                  (NUMAXISPTSY),STAT=istat)
      ENDIF
      IF(NUMAXISPTSZ.GT.0) THEN
        if (.not.allocated(INTAXISPOINTZ)) allocate(INTAXISPOINTZ
     1                              (2,NUMAXISPTSZ),STAT=istat)
        if (.not.allocated(ZAXISPTCORD)) allocate(ZAXISPTCORD
     1                                  (NUMAXISPTSZ),STAT=istat)
      ENDIF
      IF(NUMZAD.GT.0) THEN
       if (.not.allocated(POTZADC)) allocate(POTZADC(NUMZAD),
     1       STAT=istat) 
      ENDIF
      CALL VULAZAXIS()
      LNZAD=LMAX
      LNZADJ=LNZAD+3*NUMZAD
      LMAX=LNZADJ+NUMZAD
      CALL VDELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      LZADVR=LMAX
      LMAX=LZADVR+NUMZAD*IDVA
      CALL VPROMEM(LMAX)
      CALL VULAZT3(A(LNZAD),A(LZADVR),A(LNZADJ),ISNUMER)
      LICUR=LMAX
      LMAX=LICUR+NUMZAD
Cxxx
C      CALL VCZPROV(NUMZAD,1,34)
Cxxx      
C==========================================================================
C UCITAVANJE POCETNIH VREDNOSTI I SPECIFICNE GUSTINE GAMA
      CALL VISPITA(ACOZ)
      IF(IPAKT.EQ.0) THEN
         READ(ACOZ,1014) POCU,POCV,POCP,GAMA,IOSA,ISIL,IFIL,ZASIC,IPOR
         IF(ZASIC.LT.1.D-10) ZASIC=0.999
!       WRITE(*,*)'filtracione sile se racun. u zavisnosti od poroznosti'
!          IPOR=1
         IF(ISRPS.EQ.0)
     *      WRITE(IIZLAZ,3008)GAMA,IOSA,ISIL,IFIL,ZASIC,IPOR
         IF(ISRPS.EQ.1)
     *     WRITE(IIZLAZ,6008)GAMA,IOSA,ISIL,IFIL,ZASIC,IPOR
      ELSEIF(IPAKT.EQ.1) THEN
C ucitavanje pocetnih temperatura
        READ(ACOZ,2233) POCT,INDPT
       ENDIF
 2233 FORMAT(F10.3,I5)
 1014 FORMAT(4D10.3,3I5,F10.3,I5)
 3008 FORMAT(////
C     111X,'POCETNA BRZINA U PRAVCU OSE X ................U= ',F10.3//
C     111X,'POCETNA BRZINA U PRAVCU OSE Y ................V= ',F10.3//
C     111X,'POCETNI PRITISAK .............................P= ',F10.3//
     111X,'SPECIFICNA GUSTINA ........................GAMA= ',F10.3//
     111X,'PRAVAC VERTIKALNE OSE .....................IOSA= ',I5//
     111X,'INDIKATOR ZA SILE IZNAD FS (0-DA, 1-NE) ...ISIL= ',I5//
     111X,'IND. ZA FILTR. SILE(0-UKUPNE,1-FIL,2-POT) .IFIL= ',I5//  
     111X,'KOEFICIJENT ZASICENOSTI ................. ZASIC= ',1PE13.5//
     111X,'KOEFICIJENT POROZNOSTI ................... IPOR= ',I5)
 6008 FORMAT(////
C     111X,'INITIAL VELOCITY IN X DIRECTION ..............U= ',F10.3//
C     111X,'INITIAL VELOCITY IN Y DIRECTION ..............V= ',F10.3//
C     111X,'INITIAL PRESSURE .............................P= ',F10.3//
     111X,'SPECIFIC WEIGHT ...........................GAMA= ',F10.3//
     111X,'DIRECTION OF VERTICAL AXIS ................IOSA= ',I5//
     111X,'IND. FOR FORCES ABOVE FS (0-YES, 1-NO) ....ISIL= ',I5//
     111X,'IND.FOR FILTR.FORCES(0-TOTAL,1-FIL,2-POT) .IFIL= ',I5//  
     111X,'KOEFFICIENT OF SATURATION ............... ZASIC= ',1PE13.5//
     111X,'KOEFFICIENT OF POROSITY................... IPOR= ',I5)
C        
C==========================================================================
C
  130 CALL VISPITA(ACOZ)
      IF(IPAKT.EQ.0)THEN
        IF(INPT.EQ.0) READ(ACOZ,1223) MAXSIL,INDFS,KISA
        IF(INPT.EQ.1) READ(ACOZ,1222) MAXSIL,INDFS,KISA
       IF(ISRPS.EQ.0)
     * WRITE(IIZLAZ,2009) MAXSIL,INDFS,KISA
       IF(ISRPS.EQ.1)
     * WRITE(IIZLAZ,6009) MAXSIL,INDFS,KISA
      ELSEIF(IPAKT.EQ.1) THEN
C ucitavanje broja elemenata na kojima je zadat fluks, prelazanost i zracenje      
      IF(INPT.eq.0) READ(ACOZ,1015) MAXSIL,MAXTQE,MAXER,INDWATER,IFZRAC
      IF(INPT.eq.1) READ(ACOZ,1016) MAXSIL,MAXTQE,MAXER,INDWATER,IFZRAC
        IF(ISRPS.EQ.0)
     * WRITE(IIZLAZ,3088) POCT,INDPT,MAXSIL,MAXTQE,MAXER,INDWATER,IFZRAC
        IF(ISRPS.EQ.1)
     * WRITE(IIZLAZ,6088) POCT,INDPT,MAXSIL,MAXTQE,MAXER,INDWATER,IFZRAC
       ENDIF
 1222 FORMAT(I10,2I5)
 1223 FORMAT(3I5)
 1015 FORMAT(5I5)
 1016 FORMAT(5I10)
 2009 FORMAT(//,
     111X,'PODACI O GRANICNIM USLOVIMA '//
     111X,'UKUPAN BROJ POVRSINA SA ZADATIM FLUKSOM ...... MAXSIL =',I10/
     1//
     111X,'INDIKATOR SLOBODNE POVRSINE ...................  INDFS =',I5/
     116X,'EQ.0; BEZ SLOBODNE POVRSINE'/
     116X,'EQ.1; RACUNANJE SLOBODNE POVRSINE'///
     111X,'INDIKATOR ZADATE INFILTRACIJE (PADANJE KISE) .. INFILT =',I5)
 6009 FORMAT(//,
     111X,'DATA ABOUT BOUNDARY CONDITIONS '//
     111X,'TOTAL NUMBER OF SURFACES WITH PRESCRIBED FLUX .MAXSIL =',I5/
     1//
     111X,'INDICATOR FOR FREE SURFACE ....................  INDFS =',I5/
     116X,'EQ.0; WITHOUT FREE SURFACE'/
     116X,'EQ.1; FREE SURFACE CALCULATION'///
     111X,'INDICATOR FOR PRESCRIBED INFILTRATION (RAIN) .. INFILT =',I5)
 3088 FORMAT(////
     111X,'POCETNA TEMPERATURA.......................... POCT= ',F10.3/
     111X,'INDIKATOR POCETNE TEMPERATURE ................  INDPT =',I5/
     116X,'EQ.-1; POCETNA TEMPERATURA KONSTANTNA I FORMIRA ZIPOCT.IZ'/
     116X,'EQ.0; POCETNA TEMPERATURA KONSTANTNA'/
     116X,'EQ.1; POCETNA TEMPERATURA SE CITA IZ FAJLA ZIPOCT.UL'///
     411X,'UKUPAN BROJ POVRSINA SA ZADATIM FLUKSOM ...... MAXSIL =',I5/
     511X,'UKUPAN BROJ POVRSINA SA ZADATIM PRELAZOM ..... MAXTQE =',I5/
     111X,'UKUPAN BROJ POVRSINA SA ZADATIM ZRACENJA ...... MAXER= ',I5/
     111X,'INDIKATOR ZA VODEE STRUKTURE(BRANE,...)..... INDWATER= ',I5/
     111X,'INDIKATOR ZA ZRACENJE/FLUKS  ..................IFZRAC= ',I5)
 6088 FORMAT(////
     111X,'INITIAL TEMPERATRUE ......................... POCT= ',F10.3/
     111X,'INDIKATOR POCETNE TEMPERATURE ................  INDPT =',I5/
     116X,'EQ.-1; POCETNA TEMPERATURA KONSTANTNA I FORMIRA ZIPOCT.IZ'/
     116X,'EQ.0; POCETNA TEMPERATURA KONSTANTNA'/
     116X,'EQ.1; POCETNA TEMPERATURA SE CITA IZ FAJLA ZIPOCT.UL'///
     111X,'TOTAL NUMBER OF SURFACES WITH PRESCRIBED FLUX .MAXTHP= ',I5/
     111X,'TOTAL NUMBER OF SURFACES WITH PRESCRIBED CONV .MAXTQE= ',I5/
     111X,'TOTAL NUMBER OF SURFACES WITH PRESCRIBED RAD . MAXER= ',I5/
     111X,'INDICATOR FOR WATER STRUCTURES (DAMS,  ...) .INDWATER= ',I5/
     111X,'INDICATOR FOR RADIATION/FLUX  .................IFZRAC= ',I5)
C
      LNGPSI=LMAX
!       IF (NETIP.EQ.2) LMAX=LNGPSI+4*MAXSIL
      LMAX=LNGPSI+12*MAXSIL*2
      CALL VDELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      LPOVSI=LMAX
      LMAX=LPOVSI+MAXSIL*4*IDVA
      CALL VPROMEM(LMAX)
!       IF (NETIP.EQ.2)
!      1CALL VULAZT4(A(LPOVSI),NEL,A(LNGPSI),MAXSIL,NDIM)
      IF(MAXSIL.GT.0) THEN
       if (.not.allocated(IFLUXR)) allocate(IFLUXR(MAXSIL),STAT=istat) 
      ENDIF
      CALL VULA3D4(A(LPOVSI),A(LNGPSI),MAXSIL,ISNUMER,NET)
      IF(IPAKT.EQ.1) THEN
C ucitavanje elemenata na kojima je zadata prelaznost
        IF(MAXTQE.GT.0) THEN
          if (.not.allocated(HFACE)) allocate(HFACE(MAXTQE),STAT=istat) 
      if (.not.allocated(NELTOK)) allocate(NELTOK(12,MAXTQE),STAT=istat) 
          CALL VULA3D4TOK(ISNUMER,NET)
        ENDIF
C ucitavanje elemenata na kojima je zadata radijacija
        IF(MAXER.GT.0) THEN
          if (.not.allocated(NELR)) allocate(NELR(12,MAXER),STAT=istat) 
          CALL VULA3D4R(ISNUMER,NET)
        ENDIF
        IF(INDWATER.EQ.1) THEN
          CALL VISPITA(ACOZ)
          READ(ACOZ,1987) NWATERS,NSENSORS,IBOFANG
          IF(NWATERS.GT.0) THEN
            if(.not.allocated(WATER)) allocate(WATER(3,NWATERS),
     1        STAT=istat)
            CALL VULA3D4VODE() 
          ENDIF
          IF(NSENSORS.GT.0) THEN
            if(.not.allocated(SENSOR)) allocate(SENSOR(3,NSENSORS),
     1        STAT=istat)
            if(.not.allocated(HSENSOR)) allocate(HSENSOR(NSENSORS),
     1        STAT=istat) 
            CALL VULA3D4SENZORI()
          ENDIF 
          IF(IBOFANG.GE.1) THEN
            CALL VISPITA(ACOZ)
            READ(ACOZ,1988) Rd_BOFANG,D_BOFANG,Alpha
          ENDIF
       WRITE(IIZLAZ,6089) Rd_BOFANG,D_BOFANG,Alpha
 6089 FORMAT(////
     111X,'Bofang parameter..................... Rd_BOFANG= ',F10.3/
     111X,'Bofang parameter...................... D_BOFANG= ',F10.3/
     116X,'Bofang parameter......................... Alpha= ',F10.3/)
        ENDIF
C ucitavanje funkcija kojima je definisan vektor normale izvora zracenja
        IF(IFZRAC.EQ.1) THEN
          CALL VISPITA(ACOZ)
          READ(ACOZ,1123) INFX,INFY,INFZ,RKOREKCIJA,PREKIDNAFR
       WRITE(IIZLAZ,6090) INFX,INFY,INFZ,RKOREKCIJA,PREKIDNAFR
 6090 FORMAT(////
     111X,'Function for normal Vn ....................... INFX= ',I10/
     111X,'Function for normal Ve ....................... INFY= ',I10/
     111X,'Function for normal Vz ....................... INFZ= ',I10/
     111X,'Koef. korekcije merengo zracenja ..............rc= ',F10.3/
     111X,'Prekidna funkciaj vode  ..................... Zg= ',F10.3/)
        ENDIF
      ENDIF
C
C KONTURE DUZ KOJIH SE RACUNA PROTOK - cita se samo za pakp
C
      IF(IPAKT.EQ.0)THEN
       CALL VISPITA(ACOZ)
       IF(INPT.EQ.1)THEN
        READ(ACOZ,1122) NKONT,MAXLIN
       ELSE
        READ(ACOZ,1002) NKONT,MAXLIN
       ENDIF
       IF(NKONT.GT.0) THEN
        OPEN(IPROTK,FILE='PROTOK')
CN       IF (INDSC.EQ.0)THEN
          WRITE(IPROTK,3003) NKONT
!         WRITE(IPROTK,3004) (PRTOK(JJ),JJ=1,NKONT)
CN       ENDIF
       ENDIF
       IF(ISRPS.EQ.0)
     * WRITE(IIZLAZ,3217) NKONT
       IF(ISRPS.EQ.1)
     * WRITE(IIZLAZ,6217) NKONT
C
 1122  FORMAT(I5,I10)
 1123  FORMAT(3I10,2F10.4)
 1987  FORMAT(3I10)
 1988  FORMAT(3F10.4)
 3003  FORMAT(I5,'   LIN   ')
 3004  FORMAT('   TIME   ',10A13)
 3217  FORMAT(//
     111X,'UKUPAN BROJ KONTURA PROTOKA.................. NKONT= ',I10//)
 6217 FORMAT(//
     111X,'TOTAL NUMBER OF CONTOUR...................... NKONT= ',I10//)
       LKONT=LMAX
       LMAX=LKONT+9*MAXLIN*NKONT
       LKOJK=LMAX
       LMAX=LKOJK+NPT
       CALL VPROMEM(LMAX)
       CALL VULKONT(A(LKONT),ISNUMER,NET)
       CALL VICLEAR(A(LKOJK),NPT)
      ENDIF
C
      CALL VISPITA(ACOZ)
      READ(ACOZ,1001) NUMMAT
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,3218) NUMMAT
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6218) NUMMAT
 3218 FORMAT(//
     111X,'UKUPAN BROJ MATERIJALA ...................... NUNMAT= ',I5)
 6218 FORMAT(//
     111X,'TOTAL NUMBER OF MATERIALS ................... NUMMAT= ',I5)
      NJUTN=0
      INDOPT=1
      LINDPK=LMAX
      LMAX=LINDPK+NUMMAT
      CALL VDELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      LKONST=LMAX
      LMAX=LKONST+3*5*NUMMAT*IDVA
      CALL VPROMEM(LMAX)
!       IF(NETIP.EQ.1) CALL VULAZKO(A(LKONST),NUMMAT,A(LINDPK))
!       IF(NETIP.EQ.2) CALL VULAZKO(A(LKONST),NUMMAT,A(LINDPK))
      CALL VULKO3D(A(LKONST),NUMMAT,A(LINDPK))
C
C========================================================================
C UCITAVANJE VREMENSKIH FUNKCIJA:
C========================================================================
      CALL VISPITA(ACOZ)
      READ(ACOZ,1002) NTABFT,MAXTFT
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2501) NTABFT,MAXTFT
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6501) NTABFT,MAXTFT
 2501 FORMAT(//
     111X,'UKUPAN BROJ VREMENSKIH FUNKCIJA ............... NTABFT =',I5/
     111X,'MAKSIMALAN BROJ TACAKA ZA VREMENSKE FUNKCIJE .. MAXTFT =',I5
     1//)
 6501 FORMAT(//
     111X,'TOTAL NUMBER OF TIME FUNCTIONS ................ NTABFT =',I5/
     111X,'MAKSIMAL NUMBER OF POINTS FOR FUNCTIONS ....... MAXTFT =',I5
     1//)
      LVRFUN=LMAX
      LMAX=LVRFUN+NTABFT*MAXTFT*2*IDVA
      LITFMA=LMAX
      LMAX=LITFMA+NTABFT
      CALL VPROMEM(LMAX)
      CALL VULTABF(A(LVRFUN),A(LITFMA))

Cxxx
C      CALL VCZPROV(MAXSIL,1,34)
Cxxx      

C==========================================================================
      CALL VISPITA(ACOZ)
      READ(ACOZ,1000) KRAJ
      IF(KRAJ(1:4).NE.'STOP') WRITE(*,*) 'NEMA KRAJA'
      CALL VZATVOR
      CALL VSTAMPA()
      if (.not.allocated(MAXA)) allocate(MAXA(JEDN+1),STAT=istat)
      if (.not.allocated(MHT)) allocate(MHT(JEDN+1),STAT=istat)
!       LMAXA=LMAX
!       LMHT=LMAXA+JEDN+1
!       LMAX=LMHT+JEDN+1
!       CALL VPROMEM(LMAX)
!        write(*,*) "pre MAXATE"
      CALL VMAXATE(MAXA,MHT,NET,NDIM,JEDN,LM2,NWK)
!        write(*,*) "posle MAXATE"
      WRITE(IIZLAZ,*)'NWK= ',NWK
C
!       CALL VDELJIV(LMAX,2,INDL)
!       IF(INDL.EQ.0) LMAX=LMAX+1
      if (.not.allocated(TT1)) allocate(TT1(JEDN),STAT=istat)
      if (.not.allocated(TT10)) allocate(TT10(JEDN),STAT=istat)
      if (.not.allocated(PRIV)) allocate(PRIV(JEDN),STAT=istat)
      if (.not.allocated(PRIV1)) allocate(PRIV1(JEDN),STAT=istat)
      if (.not.allocated(UBRZ0)) allocate(UBRZ0(JEDN),STAT=istat)
      if (.not.allocated(UBRZ)) allocate(UBRZ(JEDN),STAT=istat)
      if (.not.allocated(BRZ0)) allocate(BRZ0(JEDN),STAT=istat)
      if (.not.allocated(BRZ)) allocate(BRZ(JEDN),STAT=istat)
      if (.not.allocated(SKEF)) allocate(SKEF(NDES,NDES),STAT=istat)
      if (.not.allocated(SKEFN)) allocate(SKEFN(NDES*(NDIM+1)/2),
     1STAT=istat)
      if (.not.allocated(AK)) allocate(AK(NDES,NDES),STAT=istat)

!       LTT1=LMAX
!       LTT10=LTT1+JEDN*IDVA
!       LPRIV=LTT10+JEDN*IDVA
!       LPRIV1=LPRIV+JEDN*IDVA
!       LUBRZ0=LPRIV1+JEDN*IDVA
!       LUBRZ=LUBRZ0+JEDN*IDVA
!       LBRZ0=LUBRZ+JEDN*IDVA
!       LBRZ=LBRZ0+JEDN*IDVA
!       LSKEF=LBRZ+JEDN*IDVA
!       LSKEFN=LSKEF+NDES*NDES*IDVA 
!       LAK=LSKEFN+NDIM*(NDIM+1)/2*IDVA
!       LMAX=LAK+NDES*NDES*IDVA
!       CALL VPROMEM(LMAX)


      IDPRIT=NDIM

      PENALT=0.D0
      IF (PENALT.GT.0.D0)
     1LMAX=LMAX+IDPRIT*NET
      CALL VDELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
C
!       LDEFOR=LMAX
!       LMAX=LDEFOR+3*NPT*IDVA
      IBRGT=2
      IBRKT=3
!       LVG=LMAX
!       LMAX=LVG+3*IBRKT*IBRKT*IBRKT*NET*IDVA
!       LGG=LMAX
!       LMAX=LGG+3*IBRKT*IBRKT*IBRKT*NET*IDVA
!       LVECTJ=LMAX
!       LMAX=LVECTJ+NPT*3*IDVA
!       LPOMER=LMAX
!       LMAX=LPOMER+NPT*3*IDVA
!       LIVECT=LMAX
!       LMAX=LIVECT+NPT*IDVA
!       LSILE=LMAX
!       LMAX=LSILE+JEDN*IDVA
      
      if (.not.allocated(DEFOR)) allocate(DEFOR(3,NPT),STAT=istat)
      if (.not.allocated(VG)) allocate(VG(3,NET,IBRKT*IBRKT*IBRKT),
     1STAT=istat)
      if (.not.allocated(GG)) allocate(GG(3,NET,IBRKT*IBRKT*IBRKT),
     1STAT=istat)
      if (.not.allocated(VECTJ)) allocate(VECTJ(3,NPT),STAT=istat)
      if (.not.allocated(POMER)) allocate(POMER(3,NPT),STAT=istat)
      if (.not.allocated(IVECT)) allocate(IVECT(NPT),STAT=istat)
      if (.not.allocated(SILE)) allocate(SILE(JEDN),STAT=istat)

C=======================================================================
C OVO JE ZA DOPUNU ELEMENATA:
      IF (IDODAV.EQ.1) THEN
       LNELG=LMAX
       NCVOR=21000
       NELEM=20000
       LMAX=LNELG+(NDIMM)*NELEM
       LCORDG=LMAX
       LMAX=LCORDG+3*NCVOR*IDVA
       LNNOVI=LMAX
       LMAX=LNNOVI+NCVOR
       CALL VPROMEM(LMAX)
       CALL VDODEL1(NEL,ID,CORD,A(LNELG),A(LCORDG),A(LNNOVI))
      ELSEIF(IREST.EQ.0.AND.NXNASTRAN.EQ.0) THEN
       if(ISOLVER.NE.-11) CALL VBLOCK1()
      ENDIF
C=======================================================================
      CALL VPROMEM(LMAX)
      WRITE(IIZLAZ,3009) MAXVEC,LMAX
      WRITE(*,*) MAXVEC,LMAX
 3009 FORMAT(//
     111X,'MAKSIMALNA DUZINA RADNOG VEKTORA IZNOSI   :',I10/
     111X,'ZA RESAVANJE PROBLEMA POTREBNO JE MEMORIJE:',I10//)

      IF (INDSC.EQ.0) THEN
        CALL VTGRAFC(CORD,NPT,18)
!         IF (NETIP.EQ.2) CALL VTGRAFE(NEL,NDIM,NET,NGE,18)
        IF (NETIP.EQ.3) CALL VTGRAF3(NEL,NDIM,NET,NGE,18)
      ENDIF
      IF (IREST.EQ.1) THEN
        CALL VSTAZP(A(LNZAD),NASLOV,KKORAK)
      ENDIF
!        write(*,*) "pre PIJEZ"
C
      NEPIJ=1
      NEPOR=1
      IPIJEZ=46
!       IF(INDSC.EQ.1.OR.INDSC.EQ.3) THEN
!         IF(IPAKT.EQ.0)THEN
!         OPEN (IPIJEZ,FILE='PIJEZ.DAT',STATUS='UNKNOWN',FORM='FORMATTED',
!      1      ACCESS='SEQUENTIAL')
!         ELSEIF(IPAKT.EQ.1)THEN
!        OPEN (IPIJEZ,FILE='PIJEZT.DAT',STATUS='UNKNOWN',FORM='FORMATTED',
!      1      ACCESS='SEQUENTIAL')
!         ENDIF
! 
!           READ(IPIJEZ,1002,err=9999) NPIJ,MAXNP
!           if (npij.gt.0) then
!              IF(MAXNP.GT.20) STOP 'MAXNP.GT.20 - PAKV1.FOR'
!              DO I=1,NPIJ
!                 READ(IPIJEZ,1002) NPIJI,MAXNPI
!                 IF(NPIJI.GT.100) STOP 'NPIJI.GT.100 - PAKV1.FOR'
!                  NODP(NPIJI)=MAXNPI
!                 DO J=1,MAXNPI
!                    IF(INPT.EQ.1)THEN
!                      READ(IPIJEZ,12220) NPIJEZ(J,I)
!                    ELSE
!                      READ(IPIJEZ,1002) NPIJEZ(J,I)
!                    ENDIF
!                 ENDDO
!              ENDDO
!           endif
!           GO TO 9994
!           
!  9999     write(*,*) ' NE POSTOJI PIJEZOMETRI '         
!           write(3,*) ' NE POSTOJE PIJEZOMETRI '
!           NEPIJ=0
! 12220     FORMAT (I10)     
! C 9997  READ(IPIJEZ,1002,ERR=9996) NPOR
! C      if (npor.gt.0) then
! C         IF(npor.GT.1000) STOP 'npor.GT.1000 - PAKV1.FOR'
! C         DO I=1,npor
! C            READ(IPIJEZ,1002) nporcv(i)
! C         ENDDO
! C      endif
! C      go to 9994
! C 9996 write(*,*) ' NE POSTOJE PORNE CELIJE '         
! C      write(3,*) ' NE POSTOJE PORNE CELIJE '
! C
!  9994     READ(IPIJEZ,1002,ERR=9995) NPORE
!           if (npore.gt.0) then
!              IF(npore.GT.1000) STOP 'npore.GT.1000 - PAKV1.FOR'
!              DO I=1,npore
!                    IF(INPT.EQ.1)THEN
!                     READ(IPIJEZ,13330) nporel(i),(cpor(jj,i),jj=1,3)
!                    ELSE
!                     READ(IPIJEZ,1333) nporel(i),(cpor(jj,i),jj=1,3)
!                    ENDIF
!              ENDDO
!           endif
!           go to 9998
!  9995     write(*,*) ' NE POSTOJE ELEMENTI U KOJIMA SU PORNE CELIJE '  
!           write(3,*) ' NE POSTOJE ELEMENTI U KOJIMA SU PORNE CELIJE '
!           NEPOR=0
!         ENDIF
        IF(INDSC.EQ.2.OR.INDSC.EQ.4) THEN
        IF(IPAKT.EQ.0)THEN
        OPEN (IPIJEZ,FILE='PIJEZ.DAT',STATUS='UNKNOWN',FORM='FORMATTED',
     1      ACCESS='SEQUENTIAL')
          IZPOD1=51
          IMECSV='pijeznivoPAK.CSV'
          OPEN (IZPOD1,FILE=IMECSV,STATUS='UNKNOWN',
     1      FORM='FORMATTED',ACCESS='SEQUENTIAL')
     
        ELSEIF(IPAKT.EQ.1)THEN
       OPEN (IPIJEZ,FILE='PIJEZT.DAT',STATUS='UNKNOWN',FORM='FORMATTED',
     1      ACCESS='SEQUENTIAL')
        ENDIF
     
        READ(IPIJEZ,356) MAX_MPOINTS,MAX_DPOINTS
          if (MAX_MPOINTS.gt.0) then
             if (.not.allocated(MPOINT_ID)) allocate(MPOINT_ID
     1                                  (MAX_MPOINTS),STAT=istat)
             if (.not.allocated(MP_ELEMENT)) allocate(MP_ELEMENT
     1                                  (MAX_MPOINTS),STAT=istat)
             if (.not.allocated(MP_COORDS)) allocate(MP_COORDS
     1                               (6,MAX_MPOINTS),STAT=istat)
             if (.not.allocated(MP_VREME)) allocate(MP_VREME
     1                               (BRKORAKA),STAT=istat)
             if (.not.allocated(MP_RESULTS)) allocate(MP_RESULTS
     1                             (BRKORAKA,MAX_MPOINTS+1),STAT=istat)
             DO I=1,MAX_MPOINTS
                 READ(IPIJEZ,358) MPOINT_ID(I),MP_ELEMENT(I),
!   promenjen format za Grancarevo za nazive termometra na A15
!                 READ(IPIJEZ,358) MPOINT_ID(I),MP_ELEMENT(I),
     1                              (MP_COORDS(jj,I),jj=1,6)
               IF(ISNUMER.EQ.1) THEN
                 NNMP=MP_ELEMENT(I)
                 JJ=NNMP-NMI+1
             IF(NNMP.LT.NMI.OR.NNMP.GT.NMA.OR.MELCV(JJ).EQ.0) THEN
             write(*,*) 'pijez', MPOINT_ID(I),I,NNMP,'van modela'
             STOP 
             ENDIF
                 MP_ELEMENT(I)=MELCV(JJ)
               ENDIF
             ENDDO
          endif
C dodato citanje tacaka u kojima se stampaju vrednosti za prikaz rezultata
          if (MAX_DPOINTS.gt.0) then
             if (.not.allocated(DPOINT_ID)) allocate(DPOINT_ID
     1                                  (MAX_DPOINTS),STAT=istat)
             if (.not.allocated(DP_ELEMENT)) allocate(DP_ELEMENT
     1                                  (MAX_DPOINTS),STAT=istat)
             if (.not.allocated(DP_COORDS)) allocate(DP_COORDS
     1                               (6,MAX_DPOINTS),STAT=istat)
!              if (.not.allocated(DP_VREME)) allocate(MP_VREME
!      1                               (MAX_DPOINTS),STAT=istat)
             if (.not.allocated(DP_RESULTS)) allocate(DP_RESULTS
     1                         (BRKORAKA,2,MAX_DPOINTS+1),STAT=istat)
             DO I=1,MAX_DPOINTS
                READ(IPIJEZ,358) DPOINT_ID(I),DP_ELEMENT(I),
     1                             (DP_COORDS(jj,I),jj=1,6)
!             write(3,*)"dpoint",DPOINT_ID(I),DP_ELEMENT(I)
               IF(ISNUMER.EQ.1) THEN
                 NNMP=DP_ELEMENT(I)
                 JJ=NNMP-NMI+1
             IF(NNMP.LT.NMI.OR.NNMP.GT.NMA.OR.MELCV(JJ).EQ.0) THEN
             write(*,*) 'tacka', DPOINT_ID(I),I,NNMP,'van modela'
             STOP 
             ENDIF
                 DP_ELEMENT(I)=MELCV(JJ)
               ENDIF
             ENDDO
          IZPOD=50
          IB=INDEX(IME,'.')
          IF (IB.EQ.0) THEN
           IMECSV=trim(IME) // '.CSV'
          ELSE
           IMECSV=IME(1:IB-1)//'.CSV'
          ENDIF
          OPEN (IZPOD,FILE=IMECSV,STATUS='UNKNOWN',
     1      FORM='FORMATTED',ACCESS='SEQUENTIAL')
          endif
        ENDIF
        
C KRAJ DELA ZA UCITAVANJE PODATAKA IZ PIJEZO.DAT
C
 9998 IZNEU=49
       IF(IPAKT.EQ.0) THEN
       IZPRO=48
      OPEN (IZPRO,FILE='PROT.NEU',STATUS='UNKNOWN',FORM='FORMATTED',
     1      ACCESS='SEQUENTIAL')
       CALL VTGRMAT(A(LKONST),NUMMAT,NETIP,IZPRO,NDIM,NASLOV,GAMA)
        CALL VTGRAUK(NPT,IZPRO,ISNUMER)      
        IF(NETIP.EQ.3.AND.NKONT.GT.0) 
     *  CALL VTGRAU3K(A(LKONT),NDIM,IZPRO)
       ENDIF
       IZPOT=47
      OPEN (IZPOT,FILE='POT.TXT',STATUS='UNKNOWN',FORM='FORMATTED',
     1      ACCESS='SEQUENTIAL')
       IF (INDSC.EQ.0.OR.INDSC.EQ.3.OR.INDSC.EQ.4) THEN
         CALL VTGRMAT(A(LKONST),NUMMAT,NETIP,IZNEU,NDIM,NASLOV,GAMA)
         CALL VTGRAUK(NPT,IZNEU,ISNUMER)
!          IF (NETIP.EQ.2) CALL VTGRAU2(NEL,NDIM,NET,IZNEU)
         CALL VTGRAU3(NDIM,NET,IZNEU,ISNUMER)
      ENDIF
C      IF (NETIP.EQ.2) CALL VTGRAU2(NEL,NDIM,NET,IZNEU)
C      IF (NETIP.EQ.3) CALL VTGRAU3(NEL,NDIM,NET,IZNEU)
C
CN      IF (INDSC.EQ.0)THEN
!       IF(IPAKT.EQ.0) THEN
       IF(IPAKT.EQ.0)  CALL VNEUTRA(A(LVREME),NPER,NASLOV,1,48)
!       ENDIF
        CALL VNEUTRA(A(LVREME),NPER,NASLOV,1,49)
!  stampanje koeficijenta filtracije u neu fajl        
        IF((NKONT.GT.0).AND.(INDSC.EQ.0.OR.INDSC.EQ.3.OR.INDSC.EQ.4))
     1        CALL VTGRMATK(A(LKONT),A(LKONST),NDIM,ISNUMER)
        IF (INDSC.EQ.0.OR.INDSC.EQ.3.OR.INDSC.EQ.4) 
     1       CALL VTGRMATS(A(LKONST),NDIM,NET,ISNUMER)
CN      ENDIF
!  kada ima preseka 
      IF(ISPRESEK.GT.0) THEN
! ucitavanje preseka (svi preseci za lamelu su u jednom fajlu)      
        CALL VREADPRESEK (ISNUMER)
! otvaranje neu fajlova za stampanje rezultata za preseke
        DO IPR=1,IPRES
          IFILE=300
          IFILE=IFILE+IPR
       SELECT CASE (IPR)
       CASE (1)
       IF ((IDJERDAP.GE.10).AND.(IDJERDAP.LE.14)) THEN
          write(IMECSV2,21) IDJERDAP
          IMEPR='L' // IMECSV2 // '-O-F.NEU'    
       if(IPAKT.EQ.1)IMEPR='L' // IMECSV2 // '-O-T.NEU'
        ELSEIF ((IDJERDAP.LT.10))THEN
          write(IMECSV1,20) IDJERDAP
          IMEPR='L' // IMECSV1 // '-O-F.NEU'    
       if(IPAKT.EQ.1)IMEPR='L' // IMECSV1 // '-O-T.NEU'
        ELSEIF ((IDJERDAP.EQ.15))THEN
          IMEPR='SIII-A6-F.NEU'    
          if(IPAKT.EQ.1)IMEPR='SIII-A6-T.NEU'
        ELSEIF ((IDJERDAP.EQ.16))THEN
          IMEPR='SII-A4-F.NEU'    
        if(IPAKT.EQ.1)IMEPR='SII-A4-T.NEU'
        ELSEIF ((IDJERDAP.EQ.17))THEN
          IMEPR='SI-A2-F.NEU'    
        if(IPAKT.EQ.1)IMEPR='SI-A2-T.NEU'
        ELSEIF ((IDJERDAP.EQ.18))THEN
          IMEPR='MB-O-F.NEU'    
        if(IPAKT.EQ.1)IMEPR='MB-O-T.NEU'
        ENDIF
       CASE (2)
       IF ((IDJERDAP.GE.10).AND.(IDJERDAP.LE.14)) THEN
          write(IMECSV2,21) IDJERDAP
          IMEPR='L' // IMECSV2 // '-OB-F.NEU'    
       if(IPAKT.EQ.1)IMEPR='L' // IMECSV2 // '-OB-T.NEU'
        ELSEIF ((IDJERDAP.LT.10))THEN
          write(IMECSV1,20) IDJERDAP
          IMEPR='L' // IMECSV1 // '-OB-F.NEU'    
       if(IPAKT.EQ.1)IMEPR='L' // IMECSV1 // '-OB-T.NEU'
        ELSEIF ((IDJERDAP.EQ.15))THEN
          IMEPR='SIII-A5-F.NEU'    
          if(IPAKT.EQ.1)IMEPR='SIII-A5-T.NEU'
        ELSEIF ((IDJERDAP.EQ.16))THEN
          IMEPR='SII-A3-F.NEU'    
        if(IPAKT.EQ.1)IMEPR='SII-A3-T.NEU'
        ELSEIF ((IDJERDAP.EQ.17))THEN
          IMEPR='SI-A1-F.NEU'    
        if(IPAKT.EQ.1)IMEPR='SI-A1-T.NEU'
        ELSEIF ((IDJERDAP.EQ.18))THEN
          IMEPR='MB-OA-F.NEU'    
         if(IPAKT.EQ.1)IMEPR='MB-OA-T.NEU'
        ENDIF
       CASE (3)
        IF ((IDJERDAP.EQ.15))THEN
          IMEPR='SIII-OA-F.NEU'    
          if(IPAKT.EQ.1)IMEPR='SIII-OA-T.NEU'
        ELSEIF ((IDJERDAP.EQ.16))THEN
          IMEPR='SII-OA-F.NEU'    
        if(IPAKT.EQ.1)IMEPR='SII-OA-T.NEU'
        ELSEIF ((IDJERDAP.EQ.17))THEN
          IMEPR='SI-OA-F.NEU'    
        if(IPAKT.EQ.1)IMEPR='SI-OA-T.NEU'
         ELSEIF ((IDJERDAP.EQ.18))THEN
           IMEPR='MB-CS-F.NEU'    
         if(IPAKT.EQ.1)IMEPR='MB-CS-T.NEU'
       ENDIF
       CASE (4)
         IF ((IDJERDAP.EQ.15))THEN
           IMEPR='SIII-O-F.NEU'    
           if(IPAKT.EQ.1)IMEPR='SIII-O-T.NEU'
         ELSEIF ((IDJERDAP.EQ.16))THEN
           IMEPR='SII-O-F.NEU'    
         if(IPAKT.EQ.1)IMEPR='SII-O-T.NEU'
         ELSEIF ((IDJERDAP.EQ.17))THEN
           IMEPR='SI-O-F.NEU'    
         if(IPAKT.EQ.1)IMEPR='SI-O-T.NEU'
         ELSEIF ((IDJERDAP.EQ.18))THEN
           IMEPR='MB-O-F.NEU'    
         if(IPAKT.EQ.1)IMEPR='MB-O-T.NEU'
        ENDIF
      CASE DEFAULT
      END SELECT
!          write(*,*),"IFILE,IME ",IFILE,IMEPR
        OPEN (IFILE,FILE=IMEPR,STATUS='UNKNOWN',
     1      FORM='FORMATTED',ACCESS='SEQUENTIAL')
!          write(*,*) "stampanje materijala preseka"
        CALL VTGRMATP(A(LKONST),NUMMAT,IFILE,NASLOV,
     1        GAMA,10)
!  za presek       
!          write(*,*) "stampanje cvorova preseka"
        CALL VTGRAUKP(IFILE,IPR)
!          write(*,*) "stampanje elemenata preseka"
        CALL VTGRAU3P(IFILE,IPR)
!          write(*,*) "stampanje perioda"
        CALL VNEUTRA(A(LVREME),NPER,NASLOV,1,IFILE)
        ENDDO
      ELSEIF(ISPRESEK.LE.-1)THEN
        CALL VREADPRESEKALL (ISNUMER)
! otvaranje neu fajlova za stampanje rezultata za preseke
        DO IPR=1,IPRES
         IFILE=300
         IFILE=IFILE+IPR
       SELECT CASE (IPR)
       CASE (1)
          IF(ISPRESEK.EQ.-1) IMEPR='L1-O-F.NEU'    
          IF(ISPRESEK.EQ.-2) IMEPR='V-OB-F.NEU'    
       CASE (2)
          IF(ISPRESEK.EQ.-1) IMEPR='L1-OB-F.NEU'    
          IF(ISPRESEK.EQ.-2) IMEPR='V-X1-F.NEU'    
       CASE (3)
          IF(ISPRESEK.EQ.-1) IMEPR='L2-O-F.NEU'    
          IF(ISPRESEK.EQ.-2) IMEPR='V-X2-F.NEU'    
       CASE (4)
          IF(ISPRESEK.EQ.-1) IMEPR='L2-OB-F.NEU'    
          IF(ISPRESEK.EQ.-2) IMEPR='V-X3-F.NEU'    
       CASE (5)
          IF(ISPRESEK.EQ.-1) IMEPR='L3-O-F.NEU'    
          IF(ISPRESEK.EQ.-2) IMEPR='V-X4-F.NEU'    
       CASE (6)
          IF(ISPRESEK.EQ.-1) IMEPR='L3-OB-F.NEU'    
          IF(ISPRESEK.EQ.-2) IMEPR='V-X5-F.NEU'    
       CASE (7)
          IF(ISPRESEK.EQ.-1) IMEPR='L4-O-F.NEU'    
          IF(ISPRESEK.EQ.-2) IMEPR='V-X6-F.NEU'    
       CASE (8)
          IF(ISPRESEK.EQ.-1) IMEPR='L4-OB-F.NEU'    
          IF(ISPRESEK.EQ.-2) IMEPR='V-X7-F.NEU'    
       CASE (9)
          IF(ISPRESEK.EQ.-1) IMEPR='L5-O-F.NEU'    
       CASE (10)
          IF(ISPRESEK.EQ.-1) IMEPR='L5-OB-F.NEU'    
       CASE (11)
          IF(ISPRESEK.EQ.-1) IMEPR='L6-O-F.NEU'    
       CASE (12)
          IF(ISPRESEK.EQ.-1) IMEPR='L6-OB-F.NEU'    
       CASE (13)
          IF(ISPRESEK.EQ.-1) IMEPR='L7-O-F.NEU'    
       CASE (14)
          IF(ISPRESEK.EQ.-1) IMEPR='L7-OB-F.NEU'    
       CASE (15)
          IF(ISPRESEK.EQ.-1) IMEPR='L8-O-F.NEU'    
       CASE (16)
          IF(ISPRESEK.EQ.-1) IMEPR='L8-OB-F.NEU'    
       CASE (17)
          IF(ISPRESEK.EQ.-1) IMEPR='L9-O-F.NEU'    
       CASE (18)
          IF(ISPRESEK.EQ.-1) IMEPR='L9-OB-F.NEU'    
       CASE (19)
          IF(ISPRESEK.EQ.-1) IMEPR='L10-O-F.NEU'    
       CASE (20)
          IF(ISPRESEK.EQ.-1) IMEPR='L10-OB-F.NEU'    
       CASE (21)
          IF(ISPRESEK.EQ.-1) IMEPR='L11-O-F.NEU'    
       CASE (22)
          IF(ISPRESEK.EQ.-1) IMEPR='L11-OB-F.NEU'    
       CASE (23)
          IF(ISPRESEK.EQ.-1) IMEPR='L12-O-F.NEU'    
       CASE (24)
          IF(ISPRESEK.EQ.-1) IMEPR='L12-OB-F.NEU'    
       CASE (25)
          IF(ISPRESEK.EQ.-1) IMEPR='L13-O-F.NEU'    
       CASE (26)
          IF(ISPRESEK.EQ.-1) IMEPR='L13-OB-F.NEU'    
       CASE (27)
          IF(ISPRESEK.EQ.-1) IMEPR='L14-O-F.NEU'    
       CASE (28)
          IF(ISPRESEK.EQ.-1) IMEPR='L14-OB-F.NEU'    
       CASE (29)
          IF(ISPRESEK.EQ.-1) IMEPR='PB-OB-F.NEU'    
       CASE (30)
          IF(ISPRESEK.EQ.-1) IMEPR='EL-OA-F.NEU'    
       CASE (31)
          IF(ISPRESEK.EQ.-1) IMEPR='SIII-A6-F.NEU'    
       CASE (32)
          IF(ISPRESEK.EQ.-1) IMEPR='SIII-A5-F.NEU'    
       CASE (33)
          IF(ISPRESEK.EQ.-1) IMEPR='SIII-OA-F.NEU'    
       CASE (34)
          IF(ISPRESEK.EQ.-1) IMEPR='SII-A4-F.NEU'    
       CASE (35)
          IF(ISPRESEK.EQ.-1) IMEPR='SII-A3-F.NEU'    
       CASE (36)
          IF(ISPRESEK.EQ.-1) IMEPR='SII-OA-F.NEU'    
       CASE (37)
          IF(ISPRESEK.EQ.-1) IMEPR='SII-O-F.NEU'    
       CASE (38)
          IF(ISPRESEK.EQ.-1) IMEPR='SI-A2-F.NEU'    
       CASE (39)
          IF(ISPRESEK.EQ.-1) IMEPR='SI-A1-F.NEU'    
       CASE (40)
          IF(ISPRESEK.EQ.-1) IMEPR='SI-OA-F.NEU'    
       CASE (41)
          IF(ISPRESEK.EQ.-1) IMEPR='MB-O-F.NEU'    
       CASE (42)
          IF(ISPRESEK.EQ.-1) IMEPR='MB-OA-F.NEU'    
       CASE (43)
          IF(ISPRESEK.EQ.-1) IMEPR='MB-CS-F.NEU'    
       CASE (44)
          IF(ISPRESEK.EQ.-1) IMEPR='SI-O-F.NEU'    
       CASE (45)
          IF(ISPRESEK.EQ.-1) IMEPR='SIII-O-F.NEU'    
      CASE DEFAULT
      END SELECT
!          write(*,*),"IFILE,IME ",IFILE,IMEPR
        OPEN (IFILE,FILE=IMEPR,STATUS='UNKNOWN',
     1      FORM='FORMATTED',ACCESS='SEQUENTIAL')
!          write(*,*) "stampanje materijala preseka"
        CALL VTGRMATP(A(LKONST),NUMMAT,IFILE,NASLOV,
     1        GAMA,10)
!  za presek       
!          write(*,*) "stampanje cvorova preseka"
        CALL VTGRAUKP(IFILE,IPR)
!          write(*,*) "stampanje elemenata preseka"
        CALL VTGRAU3P(IFILE,IPR)
!          write(*,*) "stampanje perioda"
        CALL VNEUTRA(A(LVREME),NPER,NASLOV,1,IFILE)
        ENDDO
      ENDIF
   20  format (I1)
   21  format (I2)
C
      IF (IREST.EQ.1) THEN 
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2020)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6020)
 1333 FORMAT(I5,3F10.5)
13330 FORMAT(I10,3F10.5)
  356 FORMAT(2I10)
  357 FORMAT(A10,I10,6F10.2)
  358 FORMAT(A15,I10,6F10.2)
 2020 FORMAT(///' ZAVRSENO JE UCITAVANJE I GENERISANJE ULAZNIH PODATAKA'
     1/' ULAZNI PODACI SU BEZ FORMALNIH GRESAKA'/
     2' ZA IZVRSENJE PROGRAMA UNETI:   IREST = 0'//)
 6020 FORMAT(///' READING AND GENERATING OF INPUT DATA IS OVER'/
     1' INPUT DATA ARE WITHOUT FORMAL ERROR'/
     2' FOR EXECUTION OF PROGRAM ENTER:   IREST = 0'//)
      STOP      
      ENDIF

          CALL DATE_AND_TIME(VALUES=Dtime)
          WRITE(*,*) 'vreme pre poziva racun3d', (Dtime(i),i=5,7)
          WRITE(3,*) 'vreme pre poziva racun3d', (Dtime(i),i=5,7)

      NJUTN=0
      IF (NJUTN.NE.0) THEN
 
      CALL VOPTIM1(TT1,SILE,NEL,
     1ID,A(LNZAD),A(LZADVR),A(LNGPSI),MAXA,CORD,SKEF,
     1SKEFN,A(LKONT),DEFOR,A(LVREME),A(LVRFUN),TT10,
     1UBRZ,UBRZ0,BRZ,BRZ0,AK,VECTJ,IVECT,
     1A(LPOVSI),POMER,A(LITFMA),A(LKONST),NASLOV,
     1A(LICUR),PRIV,PRIV1,VG,GG,A(LKOJK))
      ENDIF
C
C      WRITE(IIZLAZ,*) 'PRE USLOVA NETIP'
C      WRITE(IIZLAZ,*) 'NETIP =',NETIP
C      WRITE(IIZLAZ,*) LTT1,LSILE,LNEL,LID,LNZAD,LZADVR,LNGPSI,LMAXA,
C     1LCORD,LSKEF,LSKEFN,LKONT,LDEFOR,LVREME,LVRFUN,LTT10,LUBRZ,LUBRZ0,
C     1LAK,LVECTJ,LIVECT,LPOVSI,LPOMER,LITFMA,LKONST,NASLOV,LICUR
C     1,LVG,LGG,LKOJK
 1234 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(IPAKT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!       IF (NETIP.EQ.2) THEN
!       CALL VRACUN(TT1,SILE,NEL,
!      1ID,A(LNZAD),A(LZADVR),A(LNGPSI),MAXA,CORD,SKEF,
!      1SKEFN,A(LKONT),DEFOR,A(LVREME),A(LVRFUN),TT10,
!      1UBRZ,UBRZ0,AK,VECTJ,IVECT,
!      1A(LPOVSI),POMER,A(LITFMA),A(LKONST),NASLOV,
!      1A(LICUR),VG,GG,A(LKOJK))
!       ELSE
        write(*,*) "pre RACN3D"
      IF(IPAKT.EQ.0) THEN
      CALL VRACN3D(TT1,SILE,
     1A(LNZAD),A(LZADVR),A(LNGPSI),MAXA,SKEF,
     1SKEFN,A(LKONT),DEFOR,A(LVREME),A(LVRFUN),TT10,
     1UBRZ,UBRZ0,AK,VECTJ,IVECT,
     1A(LPOVSI),POMER,A(LITFMA),A(LKONST),NASLOV,
     1A(LICUR),VG,GG,A(LKOJK),ISNUMER)
      ELSEIF(IPAKT.EQ.1) THEN
      CALL VRACN3DT(TT1,SILE,
     1A(LNZAD),A(LZADVR),A(LNGPSI),MAXA,SKEF,
     1SKEFN,DEFOR,A(LVREME),A(LVRFUN),TT10,
     1UBRZ,UBRZ0,AK,VECTJ,IVECT,
     1A(LPOVSI),POMER,A(LITFMA),A(LKONST),NASLOV,
     1A(LICUR),VG,GG,INDPT,ISNUMER)
      ENDIF
      IF ((myid.ne.0).or.(mcm_kojpak.eq.7)) goto 1345
!          WRITE(*,*) "isave",isave
        write(*,*) "posle RACN3D"
!  stampa za pakt temperature u tab fajlove
      IF(IPAKT.EQ.1)THEN
      DO IMP=1,MAX_MPOINTS
        INDMPI=0
        DO JMP=1,15
            IF(lge(MPOINT_ID(IMP)(JMP:JMP),'-').and.lle(MPOINT_ID(IMP)
     1                                          (JMP:JMP),'z')) then
                INDMPI=INDMPI+1
            ENDIF
        ENDDO
!         WRITE(*,*) "INDMPI",INDMPI
!         WRITE(*,*) "MPOINT_ID(",IMP,")",MPOINT_ID(IMP)
!top WINDOWS START      
!          write(*,*) "IMP", IMP
          if (.not.allocated(MP_RESULTS_NIZ)) allocate(MP_RESULTS_NIZ
     1                               (BRKORAKA),STAT=istat)
         do i=1,BRKORAKA
!          write(*,*) "pre save koorak", I
             MP_RESULTS_NIZ(i) = MP_RESULTS(i,IMP+1)
!        write(*,*) "MP_RESULTS(i,IMP+1)",MP_RESULTS(i,IMP+1),MP_VREME(I)
        enddo    
                 CALL save_series(MPOINT_ID(IMP),MP_VREME,
     1            MP_RESULTS_NIZ,BRKORAKA,MAX_MPOINTS,IMP,INDMPI,isave)
!           write(*,*) "posle save koorak"
        if (allocated(MP_RESULTS_NIZ)) deallocate(MP_RESULTS_NIZ)
C      zile asus linux
!top WINDOWS END       
!busarac LINUX START      
c          write(*,*) "pakv1-posle raun"
c             write(*,*) BRKORAKA
c          do ii=1,BRKORAKA
c             write(*,*) "MP_RESULTS",ii,MP_RESULTS(II,IMP+1)
c          enddo
!        CALL save_series(MPOINT_ID(IMP),MP_VREME,
!     1             MP_RESULTS,BRKORAKA,MAX_MPOINTS,IMP,INDMPI,isave)
!busarac LINUX END     
      ENDDO
      ENDIF
!         write(*,*) "posle MPOINT"
!  kada ima preseka yatvaranje neu fajla
      IF((ISPRESEK.GT.0).OR.(ISPRESEK.EQ.-1)) THEN
        DO IPR=1,IPRES
          IFILE=300
          IFILE=IFILE+IPR
          CALL VNEUTRA(VREME,NPER,NASLOV,2,IFILE)
        ENDDO
      ENDIF
!      
      IF (INDSC.NE.2) THEN
      IF(IPAKT.EQ.0) CALL VNEUTRA(A(LVREME),NPER,NASLOV,2,48)
         CALL VNEUTRA(A(LVREME),NPER,NASLOV,2,49)
      ENDIF
C==========================================================================
C     ZAPISIVANJE TEMPERATURA KOJE CE BITI POCETNE
      IF(INDPT.EQ.-1) THEN
         IB=INDEX(IME,'.')
          IF (IB.EQ.0) THEN
           IZIPOCT=trim(IME) // '.UL'
          ELSE
           IZIPOCT=IME(1:IB-1)//'.UL'
          ENDIF
      
         OPEN (99,FILE=IZIPOCT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     1         ACCESS='SEQUENTIAL')
         CALL VPOCETNET(TT1,NPT,99,0)
            CLOSE (99)
      ENDIF
      CALL DATE_AND_TIME(VALUES=Dtime)
      WRITE(*,*) 'vreme na kraju', (Dtime(i),i=5,7)
      WRITE(3,*) 'vreme na kraju', (Dtime(i),i=5,7)
C==========================================================================
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2500)
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,6500)
 2500 FORMAT(//78('*')/,
     1'.............................KRAJ PRORACUNA......................
     1...............'/,78('*')/)
 6500 FORMAT(//78('*')/,
     1'...........................END OF CALCULATION....................
     1...............'/,78('*')/)
     
!          CALL DATE_AND_TIME(VALUES=Dtime)
!          WRITE(*,*) 'kraj proracuna', (Dtime(i),i=5,7)

      WRITE(*,*)'END OF CALCULATION'
 1345 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_FINALIZE(IERR)
      ! AKO JE PAKV EXE ONDA JE STOP, AKO JE LIB ONDA JE RETURN
      !STOP 
      RETURN
      END
C=======================================================================
      SUBROUTINE VSTAMPA()
      USE NODES
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
      COMMON /ULAZNI/ IULAZ,IIZLAZ,IPAKT
      COMMON /PRIKAZ/ INDSC,IZPOT
      COMMON /DODAT/ NDIMM
!       DIMENSION ID(1,*),CORD(3,*)
      
      IF (INDSC.EQ.2) RETURN
      IF(ISRPS.EQ.0)
     *WRITE(IIZLAZ,2000)
      IF(ISRPS.EQ.1)
     *WRITE(IIZLAZ,6000)
 2000 FORMAT(//
     111X,'GENERISANI CVOROVI SA BROJEVIMA JEDNACINA I KOORDINATAMA'/)
 6000 FORMAT(//
     111X,'GENERATED DATA ABOUT NODES, EQUATIONS AND COORDINATES'/)
      IF (INDSC.EQ.0)THEN
       DO 12 I=1,NPT
!         WRITE(IIZLAZ,1005) I,(ID(II,I),II=1,1),(CORD(J,I),J=1,3)
        WRITE(IIZLAZ,1005) I,(ID(1,I)),(CORD(J,I),J=1,3)
   12  CONTINUE
      ENDIF
 1005 FORMAT(I20,2X,I10,3X,3(F10.4,2X))
      RETURN
      END
C==========================================================================
C zapisivanje temperature u fajl *.ZT (ZITEMP)
C=======================================================================
!       SUBROUTINE VRESTEL(VREME,TEMP,NP,IDJSTAMP1,NP3D1,ISNUMER,KKORAK)
      SUBROUTINE VRESTEL(VREME,TEMP,NP,ISNUMER,KKORAK)
      USE NODES
      USE TEMPCVOROVI
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /DJERDAP/ IDJERDAP,ISPRESEK
      COMMON /ICITANJE/ INPT
      DIMENSION TEMP(*)
!       ,IDJSTAMP1(1,*),NP3D1(5)
      IF (ALLOCIRANEtemperature.EQ.(.FALSE.)) THEN
      ALLOCATE (TEMPuCVORU(NP))
      ALLOCIRANEtemperature = .TRUE.
      ENDIF
C
!            IF(IDJERDAP.EQ.1) THEN
!  kad je IDJERDAP > 0 zapisuje temperature za cvorove iz fajla DJ1 ...
          IFILE=105
!           IFILE=IFILE+IDJERDAP+KKORAK
!            IF(IDJERDAP.GE.1) THEN
!              DO IK=1,1
!                IFILE=IFILE+IDJERDAP
!                WRITE(IFILE,5017) NP3D1(IK)
               WRITE(IFILE,5017) NP
               DO NN=1,NP
!  provera da li se za cvor NN pisu temperature               
!                IF(IDJSTAMP1(IK,NN).EQ.IK) THEN
!                   IF(INPT.EQ.1) THEN
                IF(ISNUMER.EQ.0)THEN
                    ! TOPALOVIC OVDE ZAPISUJE TEMPERATURE!!!!!
                  TEMPuCVORU(NN)= TEMP(NN) 
                 WRITE(IFILE,5017) NN,TEMP(NN)
                ELSE
                 WRITE(IFILE,5017) NCVEL(NN),TEMP(NN)
                ENDIF
!                   ELSE
!                 WRITE(IFILE,5007) NN,TEMP(NN)
!                   ENDIF
!                ENDIF 
              ENDDO
!              ENDDO
!            ELSE
!              WRITE(19,5017) NP
!              DO NN=1,NP
! !                IF(INPT.EQ.1) THEN
!                 IF(ISNUMER.EQ.0)THEN
!                  WRITE(19,5017) NN,TEMP(NN)
!                 ELSE
!                  WRITE(19,5017) NCVEL(NN),TEMP(NN)
!                 ENDIF
! !                ELSE
! !                 WRITE(19,5007) NN,TEMP(NN)
! !                ENDIF
!              ENDDO
!            ENDIF
!       IF(IDJERDAP.EQ.1) WRITE(II,5017) NP3D1
!          DO NN=1,NP
!           IF(IDJERDAP.EQ.1) THEN
!              IF(ISTAMP1(NN).EQ.1) THEN
!                 IF(INPT.EQ.1) THEN
!                       WRITE(II,5017) NN, TEMP(NN)
!                 ELSE
!                       WRITE(II,5007) NN,TEMP(NN)
!                 ENDIF
!              ENDIF 
!           ELSE
!              IF(INPT.EQ.1) THEN
!                    WRITE(II,5017) NN,TEMP(NN)
!              ELSE
!                    WRITE(II,5007) NN,TEMP(NN)
!              ENDIF
!           ENDIF
!          ENDDO
C      WRITE (II) VREME,(TEMP(I),I=1,NP)
      close (IFILE)
 5007 FORMAT(I5,1PE13.5)
 5017 FORMAT(I10,1PE13.5)
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE VPOCETNET(TEMP,NP,II,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TEMP(*)
         IF(IND.EQ.0) THEN
         WRITE (II) (TEMP(I),I=1,NP)
      ELSE
         READ (II) (TEMP(I),I=1,NP)
      ENDIF
      RETURN
      END
