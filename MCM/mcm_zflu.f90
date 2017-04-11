!C==========================================================================
!C==========================================================================
      SUBROUTINE mcm_zflu
      use mcm_database
      USE ZFLUID
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IFILE=61
      NSTART = -10
      NPTI = 0

        DO i = 1,mcm_np
        ! second material is the contact material hardoced solution
        if (par(i).mat.eq.2) then
            if (NSTART.LT.0) then
                NSTART = i
            end if
            NPTI =NPTI +1
        end if
        ENDDO
              IF (ALLOCIRANECVORNESILE.EQ.(.FALSE.)) THEN
      ALLOCATE ( NODEBROJ(NPTI))
      ALLOCATE (CVORNESILE(NPTI,3))
      
      
      INODEBROJAC =1
      DO I=NSTART,NPTI+NSTART-1
            NODEBROJ(INODEBROJAC) = I
         CVORNESILE(INODEBROJAC,1) = 0
         CVORNESILE(INODEBROJAC,2) = 0
         CVORNESILE(INODEBROJAC,3) = 0
         INODEBROJAC = INODEBROJAC+1
      ENDDO
      
      
      ALLOCIRANECVORNESILE = .TRUE.
      END IF
     
      !  
      !OPEN(IFILE,FILE='ZFLUID')
      !WRITE(IFILE,200) NPTI
!        INODEBROJAC =1
!      DO I=NSTART,NPTI+NSTART-1
!            NODEBROJ(INODEBROJAC) = I
!        WRITE(IFILE,300) I,((par(I)%h*par(I)%h)*par(I)%sigma(J,J),J=1,3)
!         !WRITE(IFILE,300) I,((par(I)%mass)*par(I)%a(J),J=1,3)
!!         CVORNESILE(INODEBROJAC,1) = (par(I)%mass)*par(I)%a(1)
!!         CVORNESILE(INODEBROJAC,2) = (par(I)%mass)*par(I)%a(2)
!!         CVORNESILE(INODEBROJAC,3) = (par(I)%mass)*par(I)%a(3)
!         CVORNESILE(INODEBROJAC,1) = par(I)%h*par(I)%h*par(I)%sigma(1,1)
!         CVORNESILE(INODEBROJAC,2) = par(I)%h*par(I)%h*par(I)%sigma(2,2)
!         CVORNESILE(INODEBROJAC,3) = par(I)%h*par(I)%h*par(I)%sigma(3,3)
!         INODEBROJAC = INODEBROJAC+1
!      ENDDO

      !CLOSE(IFILE)      

 200  FORMAT (I5)
 300  FORMAT (I5,3D13.5)

      END
!C==========================================================================
!C==========================================================================