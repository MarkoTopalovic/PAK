C=======================================================================
C
C=======================================================================
      program main
      use mcm_database
      USE CVOROVI
      USE WSTAZK
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical newproblem
      INTEGER IVTKCOUNTER
C		  
C ......................................................................
C .
C .                    - P A K N E L -
C .
CE.        FINITE ELEMENT PROGRAM FOR NONLINEAR ANALYSIS
CE.        PROGRAM MAIN
CS.        PROGRAM ZA NELINEARNU ANALIZU KONSTRUKCIJA
CS.        GLAVNI PROGRAM
C .
C ......................................................................
C
CE    NTOTAL - MAXIMUM TOTAL STORAGE AVAILABLE IN BLANK COMMON
CS    NTOTAL - MAKSIMALAN PROSTOR U VEKTORU A
C
      include 'paka.inc'
      INCLUDE 'mpif.h'
      COMMON/VERSION/ IVER
      COMMON /VTKVALUES/ VTKIME,IVTKCOUNTER
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,IOPGL(6),KOSI,NDIN,ITEST
      CHARACTER*6    FIPAKS,FIPAKF
      DIMENSION IA(1)
      EQUIVALENCE(A(1),IA(1))
      DIMENSION NKDT(100),DTDT(100)
      integer ierr, myid
      ALLOCIRANICVOROVI = .FALSE.
      ALLOCIRANIELEMENTI = .FALSE.
      ALLOCIRANAPOMERANJA = .FALSE.
      ALLOCIRANECVORNESILE = .FALSE.
      ALLOCIRANAMATRICA = .FALSE.
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      IF (myid.ne.0) goto 10
        
        IVER=0
        IF(IVER.EQ.1) THEN
          PRINT *, '     - STUDENTSKA VERZIJA PROGRAMA -'
          PRINT *, 'NE SME SE KORISTITI U KOMERCIJALNE SVRHE'
          PRINT *, ' '
        ENDIF
C
CE      ANALYSIS TYPE (=1-PAKS, =2-PAKF, =3-PAKS+PAKF, =4-MCM, =5-PAKS+MCM)
CS      KOJI SE PROGRAM KORISTI (=1-PAKS, =2-PAKF, =3-ZAJEDNO, 4-MCM )
C       WRITE(*,*)'Finite element program for fluid-structure interaction'
C       WRITE(*,*) 'Options: '
C       WRITE(*,*) '1 - Structure analysis only'
C       WRITE(*,*) '2 - Fluid analysis only'
C       WRITE(*,*) '3 - Interaction fluid-structure analysis'
C       WRITE(*,*) '4 - Meshless Continuum Mechanics MCM sph'
C       WRITE(*,*) '5 - PAKS + MCM'
C       READ(*,*) KOJPAK
        KOJPAK=1
        mcm_kojpak = KOJPAK
        IF(KOJPAK.EQ.0) KOJPAK=3
CE      MEMORY INDICATOR (=0-ENOUGH, =1-NOT ENOUGH)
CS      DA LI IMA PROSTORA ZA SVE U MEMORIJI (=0-IMA, =1-NEMA)
        NEMADE=0
CS      PAZNJA - OVAJ USLOV NE MORA DA BUDE UKLJUCEN AKO IMA MEMORIJE
        IF(KOJPAK.EQ.3) NEMADE=1
CE      INDICATOR FOR MORE EXAMPLES IN ONE INPUT FILE (=0-NO, =1-YES)
CS      VISE ULAZNIH DATOTEKA U JEDNOJ (=0-NE, =1-DA)
        IMAJOS=0 
CE      IDVA - IS PRECISION PARAMETER SINGLE(=1)/DOUBLE(=2)
CS      INDIKATOR ZA DVOSTRUKU PRECIZNOST
        IDVA = 2
CE      LMAX IS START POINTER IN ARRAY A(LMAX)
CS      LMAX JE POCETNI REPER U VEKTORU A(LMAX)
        LMAX = 1
        IF(KOJPAK.EQ.3) LMAX = 201
C
CE      PAKF - INPUT DATA 
CS      PAKF - ULAZNI PODACI
C
CE      START POINTER FOR PAKF
CS      POCETNI REPER ZA PAKF
        LPAKF=LMAX
        IF(KOJPAK.EQ.2.OR.KOJPAK.EQ.3) THEN
C
CE         INPUT DATA AND MEMORY ALLOCATION
CS         ULAZNI PODACI FOR PAKS I REZERVISANJE MEMORIJE
C
           CALL INPAKF(A,NTOTAL,LMAX,NKDT,DTDT,NPER)
C
CE         FROM LPAKF TO LSKF ARE IRREMOVABLE DATA FOR PAKF IN ARRAY A(*)   
CS         OD LPAKF DO LSKF SU STALNI PODACI ZA PAKF   
C
           CALL DELJIV(LMAX,2,INDL)
           IF(INDL.EQ.0) LMAX=LMAX+1
           LSKF=LMAX
           KOLKF=0
           IF(KOJPAK.EQ.3.AND.NEMADE.EQ.1) THEN
              FIPAKF='ZIPAKF'
              IPAKF=51
              OPEN (IPAKF,FILE=FIPAKF,STATUS='UNKNOWN',
     1              FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
              KOLKF=(LSKF-LPAKF+1)/IDVA
              IF(KOLKF.GT.0) CALL WRITED(A(LPAKF),KOLKF,IPAKF)
           ENDIF
        ENDIF
CS    LMAXF IS FIRST POINTER AFTER IRREMOVABLE DATA FOR PAKF
CS    LMAXF JE REPER DO KOGA JE UZEO MEMORIJU PAKF
        LMAXF=LMAX
        IF(NEMADE.EQ.1) LMAXF=LPAKF
        KOLKO=0
        IVTKCOUNTER = 0
C
CE    PAKS - INPUT DATA 
CS    PAKS - ULAZNI PODACI
C
        LMAX=LMAXF
CE    START POINTER FOR PAKS
CS    POCETNI REPER ZA PAKS
        CALL DELJIV(LMAX,2,INDL)
        IF(INDL.EQ.0) LMAX=LMAX+1
        LPAKS=LMAX
        LSK=LMAX
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SAMO MCM BEZ PAKA!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        
        IF(KOJPAK.EQ.4) THEN 
! get input file names
            call mcm_startup(newproblem)
            if(newproblem) then
! read input file
                call mcm_getinput
! problem initialisation
                call mcm_initial
                mcm_init_ts = mcm_dt
            else
! read restart file
!call mcm_restart
            endif         
        ENDIF 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SAMO MCM BEZ PAKA!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        
        
 9999   IF(KOJPAK.EQ.1.OR.KOJPAK.EQ.3.OR.KOJPAK.EQ.5) THEN
           IF (myid.ne.0) goto 10
           LMAX=LMAXF
           LPAKS=LMAX
C
CE       INPUT DATA AND MEMORY ALLOCATION
CS       ULAZNI PODACI FOR PAKS I REZERVISANJE MEMORIJE
C
           CALL UPAKS(IMAJOS,NTOTAL,NKDT,DTDT,NPER,LIPODS,LMAX,LSK,LMM,
     +              NGENL,LKAKO6)
           IF(KOJPAK.EQ.3) CALL IJEDN1(A(1),A(LIPODS),200)
           
      IF(KOJPAK.EQ.5) THEN 
! ako je proracuna PAK+MCM ovde se ucitavaju SPH vrednosti
! get input file names
            call mcm_startup(newproblem)
            if(newproblem) then
! read input file
                call mcm_getinput
! problem initialisation
                call mcm_initial
                mcm_init_ts = mcm_dt
            else
! read restart file
!call mcm_restart
            endif     
        ENDIF
           
C
CE       FROM LPAKS TO LSK ARE IRREMOVABLE DATA FOR PAKS IN ARRAY A(*)   
CS       OD LPAKS DO LSK SU STALNI PODACI ZA PAKS   
C
           KOLKO=0
           IF(KOJPAK.EQ.3.AND.NEMADE.EQ.1) THEN
              FIPAKS='ZIPAKS'
              IPAKS=48
              OPEN (IPAKS,FILE=FIPAKS,STATUS='UNKNOWN',
     1              FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
              KOLKO=(LMM-LPAKS+1)/IDVA
              IF(KOLKO.GT.0) CALL WRITED(A(LPAKS),KOLKO,IPAKS)
           ENDIF
        ENDIF
CE    LMAXS IS FIRST POINTER AFTER IRREMOVABLE DATA FOR PAKS
CS    LMAXS JE REPER DO KOGA JE UZEO MEMORIJU PAKS
        LMAXS=LSK
        IF(NEMADE.EQ.1) LMAXS=LPAKS
C
CE    TIME PERIOD LOOP
CS    OSNOVNA PETLJA PO VREMENSKIM PERIODIMA
C     -----------
        IF(KOJPAK.NE.2.AND.KOJPAK.NE.4) CALL VREM(2)
        IF(KOJPAK.EQ.3) CALL IDENTI(A(1),LPAKF,IPAKF,KOLKF)
      IF(KOJPAK.NE.2.AND.KOJPAK.NE.4)
     1 CALL NEUTRA(NKDT,DTDT,1)
      
      IF(KOJPAK.EQ.2.OR.KOJPAK.EQ.3)
     1 CALL NEUTRAF(NKDT,DTDT,NPER,1,69)
10    CALL PERIOD(NKDT,DTDT,NPER,LIPODS,KOJPAK,NEMADE,
     +            LPAKS,KOLKO,IPAKS,NGENL,
     +            LPAKF,KOLKF,IPAKF,LSKF)
      IF (myid.ne.0) goto 20
      IF(KOJPAK.NE.2.AND.KOJPAK.NE.4) CALL VREM(3)
C     -----------
C
CE    PAKS - EIGEN VALUES CALCULATION
CS    PAKS - RACUNANJE SOPSTVENIH VREDNOSTI
C
20    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KOJPAK,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      IF(KOJPAK.EQ.1.OR.KOJPAK.EQ.3.OR.KOJPAK.EQ.5) THEN
         IF (myid.ne.0) goto 30
            IF((KOJPAK.EQ.3.OR.KOJPAK.EQ.5)
     1        .AND.NEMADE.EQ.1.AND.KOLKO.GT.0) THEN
            REWIND IPAKS
            CALL READD(A(LPAKS),KOLKO,IPAKS)
         ENDIF
30       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL SPAKS(IMAJOS)
      ENDIF
40    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(IMAJOS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      IF(IMAJOS.GT.0) GO TO 9999
      IF (myid.ne.0) goto 50
      IF(KOJPAK.NE.2.AND.KOJPAK.NE.4)
     1 CALL NEUTRA(NKDT,DTDT,2)
      IF(KOJPAK.EQ.2.OR.KOJPAK.EQ.3)
     1 CALL NEUTRAF(NKDT,DTDT,NPER,2,69)

C
50    CALL MPI_FINALIZE(IERR)

      
      IF (ALLOCIRANAMATRICA.EQ.(.TRUE.)) THEN
      DEALLOCATE (RTWRITE)
      ENDIF

      STOP
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE PERIOD(NKDT,DTDT,NPER,LIPODS,KOJPAK,NEMADE,
     +                  LPAKS,KOLKO,IPAKS,NGENL,
     +                  LPAKF,KOLKF,IPAKF,LSKF)
      use mcm_database
      USE ELEMENTI
      USE CVOROVI
      USE WSTAZK
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        WITH LOOP OVER TIME PERIODS AND STEPS
CS.    P R O G R A M
CS.        SA PETLJOM PO VREMENSKIM PERIODIMA I KORACIMA
C .
CE.    V A R I A B L E S
CE.        NPER  - TOTAL NUMBER OF PERIODS WITH CONSTANT TIME STEP,
CE.                SEE CARD /3/
C .
C ......................................................................
C
      CHARACTER*20    VTKIME
      INTEGER IVTKCOUNTER
      include 'paka.inc'
       
      COMMON /KOJISR/ KOJSET
      COMMON /AUTSTK/ DELUM,DELFM,KRAJP
      common /crklie/ icrkli(100000)
      INCLUDE 'mpif.h'
      COMMON /EXPLICITNA/ INDEXPL
      COMMON /KRITKOR/ DTCR
      COMMON /EXPLALOK/ IUBRZ,IBRZINA,IBRZINAIPO,IPOMAK,IMASA,IPRIGUSEN
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPERMT,IOPGL(6),KOSI,NDIN,ITEST
      COMMON /VTKVALUES/ VTKIME,IVTKCOUNTER
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA
     1 ,LMXAU,LAPRS
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      DIMENSION NKDT(*),DTDT(*)
      DIMENSION IS(1)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      IF (myid.eq.0) THEN
c
      call iclear(icrkli,100000)
C
      KRAJP=0
      KOJSET=0
      IS(1)=0
      KOCID=0
      VREM0=0.D0
      INDT=0
!!!!!!!!!!!!!!!!!!!!!!!!!! MCM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(KOJPAK.EQ.4) THEN
      call mcm_solution
! Write dynamic relaxation data
      call mcm_write_drelax
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!! MCM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      IF(INDEXPL.EQ.1) THEN
!C
!CE      BASIC LOOP OVER TIME PERIODS OF EXPLICIT INTEGRATION
!CS      OSNOVNA PETLJA PO VREMENSKIM PERIODIMA ZA EKSPLICITNU INTEGRACIJU
!C
!        LPUU=IPOMAK
!        LPUV=IBRZINA
!        LPUA=IUBRZ
!        ISUU=0
!        ISUV=0
!        ISUA=0
!        KORBR = 1
!        BROJAC = 99
!        DT=DTDT(KORBR)
!        VREME=NKDT(KORBR)*DT
!        DTCR=DT*10
!        DTP=DT
!        CALL INTKMM
!C        
!        DO
!            CALL EXPLINTGR(LIPODS,A(LRTDT),A(IUBRZ),A(IBRZINA),
!     1           A(IBRZINAIPO),A(IPOMAK),A(IMASA),A(IPRIGUSEN),
!     2           DT,DTP,DTCR,VREM0,KORBR,KOJPAK)
!            VREM0 = VREM0 + DT
!            DTP=DT
!            IF (DT>DTCR) DT=0.9*DTCR
!                IF(BROJAC.EQ.0) THEN 
!                    CALL STAMP
!                    CALL STAGP
!                    BROJAC = 99
!                ENDIF
!            IF(VREM0.GT.VREME) exit
!            BROJAC = BROJAC - 1
!            KORBR = KORBR + 1
!        ENDDO
!        WRITE(*,*) 'posle petlje po vremenu, trenutno vreme =',VREM0
!        WRITE(*,*) 'posle petlje po vremenu, ukupno vreme =',VREME
!        pause
!          ENDIF ! IF EXPLICIT
      ENDIF ! IF (myid.eq.0)
C
CE    BASIC LOOP OVER TIME PERIODS

CS    OSNOVNA PETLJA PO VREMENSKIM PERIODIMA
C
10    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NPER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      DO 100 IPER=1,NPER
         IF (myid.eq.0) then
            IINDT=NKDT(IPER)
            DT=DTDT(IPER)
            IPDT=INDT+1
            INDT=INDT+IINDT
         endif
C
CE    BASIC LOOP OVER TIME STEPS
CS    OSNOVNA PETLJA PO VREMENSKIM KORACIMA
C
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(IPDT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(INDT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(KRAJP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         DO 500 KORBR=IPDT,INDT
C          END OF ANALYSIS (KRAJP.EQ.1)
	      IF(KRAJP.EQ.1) RETURN
            IF (myid.ne.0) goto 20
            VREM0 = VREM0 + DT
            
C
C          PAKF - PROGRAM
C
             IF(KOJPAK.EQ.2.OR.KOJPAK.EQ.3) THEN
                IF(KOJPAK.EQ.3.AND.NEMADE.EQ.1.AND.KOLKF.GT.0) THEN
                   REWIND IPAKF
                   CALL READD(A(LPAKF),KOLKF,IPAKF)
                ENDIF
                CALL PPAKF(A,LPAKF,IPER,LSKF,IS,VREM0,KORBR)
                IF(KOJPAK.EQ.3.AND.NEMADE.EQ.1.AND.KOLKO.GT.0) THEN
                   REWIND IPAKF
                   CALL WRITED(A(LPAKF),KOLKF,IPAKF)
                ENDIF
             ENDIF
C
C            PAKS - PROGRAM
C
 20          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(KOJPAK,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             IF(KOJPAK.EQ.1.OR.KOJPAK.EQ.3.OR.KOJPAK.EQ.5) THEN
                IF (myid.eq.0) THEN
                    IF((KOJPAK.EQ.3.OR.KOJPAK.EQ.5)
     1            .AND.NEMADE.EQ.1.AND.KOLKO.GT.0) THEN
                        REWIND IPAKS
                        CALL READD(A(LPAKS),KOLKO,IPAKS)
                    ENDIF
                ENDIF
  30            CALL PPAKS(A(LIPODS),KOCID,IPDT,DT,VREM0,KORBR,KOJPAK)
                IF (myid.eq.0) THEN
                    IF((KOJPAK.EQ.3.OR.KOJPAK.EQ.5)
     1               .AND.NEMADE.EQ.1.AND.KOLKO.GT.0) THEN
                        REWIND IPAKS
                        CALL WRITED(A(LPAKS),KOLKO,IPAKS)
                    ENDIF
                ENDIF

                IF(KOJPAK.EQ.1.OR.KOJPAK.EQ.5) THEN     
                      IHELP = IVTKCOUNTER
                      IVTKCOUNTER = IHELP+1
                      CALL VTKPRINT(A(LPAKS))
                      CALL VTKPRINTMODULEELEMENT
                END IF
  40         ENDIF
C
            
!
!!!!!!!!!!!!!!!!!!!!!!!!!! MCM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             IF(KOJPAK.EQ.5) THEN            
             if (pak_initialized.eq.(.false.))then
             pak_initialized = .true.
             mcm_oldtime = 0.0_d
             else
             mcm_oldtime = mcm_endtime
             end if
             mcm_endtime = VREM0
                call mcm_solution
! Write dynamic relaxation data
                call mcm_write_drelax
             ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!! MCM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  500   CONTINUE 
C
  100 CONTINUE
C  
!      ENDIF
C

      
      IF(KOJPAK.EQ.4.OR.KOJPAK.EQ.5) THEN
      DEALLOCATE (VTKELEMENTI)
!      DEALLOCATE (VTKECVOROVI)
!      DEALLOCATE (NODEBROJ)
        CALL mcm_shutdown(1)
      ENDIF
      RETURN
      
      END