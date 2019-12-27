PROGRAM WSMTD_rw_2D
!---------------------------------------------------------------------------------------------------------------------!
!FORTRAN PROGRAM WRITTEN TO COMPUTE UNBIASED DISTRIBUTION (1D AND 2D)FROM TASS SIMULATION ALONG USER DEFINED          !
!COLLECTIVE COORDINATES..                                                                                             !
!ORIGINAL CODE WRITTEN BY Shalini Awasthi (ashalini@iitk.ac.in)                                                       !
!MODIFIED BY Rahul Verma (vrahul@iitk.ac.in)                                                                          !
!                                                                                                                     !
!kt = System Temeprature ; kt0 = Extended CV Temperature ; bias_fact = Biased Factor for MTD ; ct = ct Factor         !
!t_min = Minimum MD steps ; t_max = Maximum MD steps ; narg = Argumets ; ncv = Numer of CV ; v_baias = Total MTD bias !
!cv_mtd = m = MTD CV index ; u_cv = u = Umbrella CV index ; t_cv = t = Temeprature CV index                           !
!UCV = U cv argument ; MTD = MTD cv argument ; Prob_nD = Dimension of Unbiased Probability                            !
!CV_num = Probability is the Dimension of ; pfrqMD = Print Frequency argument ; w_cv Print Frequency in cvmdck File   !
!dtMTD  = Print Frequency argument for MTD bias added ; w_hill = Print Frequency in colvar File                       !
!gridmin = Minimum Grid Size = gridmax = Maximum Grid Size ; griddif = Grid Difference                                !
!width = Hill Width of Gaussian Bias in MTD ; ht = Hill Height of Gaussian Bias in MTD ; ht = MTD CV Displacement     !
!kb = Boltzman Constant in A.U. ; prob = 1D Probability ; prob_2D = 2D Probability                                    !
!---------------------------------------------------------------------------------------------------------------------!

IMPLICIT NONE
REAL*8              :: den,dum,alpha,Ro,kt0,kt,ktb,s1,s2,bias_fact
REAL*8, ALLOCATABLE :: cv1(:),prob(:),fes(:),fes1(:),grid(:)
REAL*8, ALLOCATABLE :: dummy(:,:),cv(:,:),prob_2D(:,:),prob_3D(:,:,:)
REAL*8, ALLOCATABLE :: gridmin(:),gridmax(:),griddif(:),vbias(:),ct(:)
INTEGER,ALLOCATABLE :: nbin(:),indx(:),t(:),t_cv(:)
INTEGER :: md_steps,mtd_steps,dummy1,i,j,index1,k,t_min,t_max,i_s1,i_s2,narg
INTEGER :: i_md,i_mtd,ncv,w_hill,w_cv,n1,n2,n3,n4,prob_nD,cv_mtd,cv_us,cv_num(3)
INTEGER :: ii,jj,kk,l,u,m
LOGICAL :: pmf,inpgrid
CHARACTER*5 :: mtd
CHARACTER*120 :: arg 
REAL*8, PARAMETER :: kb = 1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: au_to_kcal = 627.51 
REAL*8, PARAMETER :: kj_to_kcal = 0.239006 

OPEN(11,FILE='cvmdck_mtd',STATUS='unknown')
OPEN(12,FILE='parvar_mtd',STATUS='unknown')
OPEN(13,FILE='colvar_mtd',STATUS='unknown')
OPEN(14,FILE='cv.dat',STATUS='unknown')
!
CALL get_steps(11,md_steps)
CALL get_steps(13,mtd_steps)
!
kt0       = 300.D0   ; kt      = 300.D0
t_min     = 1        ; t_max   = md_steps
pmf       = .FALSE.  ; inpgrid = .false.
t_max     = md_steps ; narg    = IARGC()
bias_fact  = 1500.D0

DO i=1,narg
  CALL GETARG(i,arg)
  IF(INDEX(arg,'-T0').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt0
  ELSEIF(INDEX(arg,'-T').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt
  ELSEIF(INDEX(arg,'-bias_fact').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)bias_fact
  ELSE IF(INDEX(arg,'-tmin').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_min
  ELSE IF(INDEX(arg,'-tmax').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_max
     IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'
   ELSEIF(INDEX(arg, '-ncv') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)ncv
   ELSEIF(INDEX(arg, '-UCV') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)cv_us
   ELSEIF(INDEX(arg, '-MTD') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)mtd
   ELSEIF(INDEX(arg, '-MCV') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)cv_mtd
   ELSEIF(INDEX(arg, '-Prob_nD') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)prob_nD
       IF (prob_nD .gt. 3) STOP  "MORE THAN 3 DIMENSION IS NOT IMPLIMENTED"
   ELSEIF(INDEX(arg, '-CV_num') .NE. 0)THEN
         DO ii = 1,prob_nD
           CALL GETARG(i+ii,arg)
           READ(arg,*)cv_num(ii)
         ENDDO
  ELSE IF(INDEX(arg,'-grid').NE.0)THEN

ALLOCATE(gridmin(ncv),gridmax(ncv),griddif(ncv)) ; ALLOCATE(indx(ncv))

j = 0 
DO ii = 1,ncv
j = j + 1
      CALL GETARG(i+j,arg)
      READ(arg,*)gridmin(ii)
      CALL GETARG(i+j+1,arg)
      READ(arg,*)gridmax(ii)
      CALL GETARG(i+j+2,arg)
      READ(arg,*)griddif(ii)
j = j + 2 
100 FORMAT (I4,2X,3F16.4)
!WRITE(*,100)ii,gridmin(ii),gridmax(ii),griddif(ii)
ENDDO
      inpgrid=.true.
  ELSE IF(INDEX(arg,'-pfrqMD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_cv
  ELSE IF(INDEX(arg,'-dtMTD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_hill
  END IF
END DO
IF (mtd .eq. 'y' .and. cv_mtd .eq. 0) STOP "***ERROR !! PLEASE SPECIFY METADYNAMICS CV COORDINATE INDEX"
WRITE(*,'(A)')'!------------------------------------------------------------------------------------'
   WRITE(*,'(A,I10)')'! No: of MD  steps                            =',md_steps
   WRITE(*,'(A,I10)')'! No: of max MD  steps                        =',t_max
   WRITE(*,'(A,I9)')'! No: of min MD  steps                        =',t_min
   WRITE(*,'(A,I7)')'! No: of CV                                   =',ncv
   WRITE(*,'(A,I7)')'! Umbrella CV Column                          =',cv_us
IF (mtd .eq. 'n') THEN
   WRITE(*,'(A,6X,A)')'! Metadynamics is Enabled                     =','NO'
ELSEIF(mtd .eq. 'y') THEN
   WRITE(*,'(A,5X,A)')'! Metadynamics is Enabled                     =','YES'
   WRITE(*,'(A,2X,I5)')'! Metadynamics CV Column                      =',cv_mtd
ENDIF
   WRITE(*,'(A,F11.2)')'! Physical System Temperature T0 (K)          =',kt0
   WRITE(*,'(A,F11.2)')'! CV Temperature T (K)                         =',kt
   WRITE(*,'(A,F12.2)')'! Bias Factor (K)                             =',bias_fact
   WRITE(*,'(A,I7)')'! Print Freq. in cvmdck_mtd file              =',w_cv
   WRITE(*,'(A,I8)')'! Freq. of Hill Update                        =',w_hill
   WRITE(*,'(A,I7)')'! Dimension of Probability                    =',prob_nD
   WRITE(*,'(A,3I7)')'! Probability Index                           =',cv_num(1:prob_nD)
WRITE(*,'(A)')'!------------------------------------------------------------------------------------'
ii = cv_num(1); jj = cv_num(2) ; kk = cv_num(3)
u  = cv_us    ; m  = cv_mtd  !! #### m = METADYNAMICS CV INDEX !! #### u = UMBRELLA CV INDEX
!-----------------------------------------------------------------------------------------------------!
WRITE(*,'(A85)')'=========================================================================================='! 
IF (ii .ne. u .and. jj .ne. u .and. kk .ne. u) THEN
WRITE(*,'(10X,A)') "!!SORRY!! PLEASE SPECIFY UMBRELLA CV INDEX PROPERLY" 
WRITE(*,'(A85)')'=========================================================================================='!
STOP ; ENDIF
IF (mtd .eq. 'y' .and. cv_mtd .eq. cv_us) THEN
!IF (ii .ne. m .and. jj .ne. m .and. kk .ne. m) THEN 
WRITE(*,'(10X,A)')"METADYNAMICS ENABLED, CAN'T UNBIAS PROBABILITY WITHOUT UNBIASING 'MTD'"
WRITE(*,'(A85)')'=========================================================================================='!
STOP
ENDIF !; ENDIF 
!-----------------------------------------------------------------------------------------------------!
ALLOCATE(cv(ncv,md_steps)) ; ALLOCATE(dummy(ncv,md_steps))
ALLOCATE(nbin(ncv))        ; ALLOCATE(vbias(md_steps))
ALLOCATE(ct(mtd_steps))    ; ALLOCATE(t(ncv))
ALLOCATE(t_cv(ncv))
!-----------------------------------------------------------------------------------------------------!

101 FORMAT (I5,4F16.6)
DO i_md=1,md_steps
   READ(11,*)dummy1,dummy1,(dummy(j,i_md), j=1,ncv),(cv(j,i_md) ,j=1,ncv)
!   WRITE(*,*)dummy1,(dummy(j,i_md), j=1,ncv),(cv(j,i_md) ,j=1,ncv)
   WRITE(14,101)dummy1,(cv(j,i_md) ,j=1,ncv)
END DO
  DO i = 1,ncv 
   nbin(i) = NINT((gridmax(i)-gridmin(i))/griddif(i)) + 1
  ENDDO

j = 0
WRITE(*,'(9X,4A9)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
WRITE(*,'(A10,3F8.4,I10)')'US   COORD:', gridmin(u),gridmax(u),griddif(u),nbin(u)
IF (mtd .eq. 'y') THEN
WRITE(*,'(A10,3F8.4,I10)')'MTD  COORD:', gridmin(m),gridmax(m),griddif(m),nbin(m)
ENDIF
DO i = 1,ncv ; IF (i .ne. u .and. i .ne. m) THEN
WRITE(*,'(A10,3F8.4,2I10)')'TASS COORD:', gridmin(i),gridmax(i),griddif(i),nbin(i)
j = j + 1 ; t_cv(j) = i
ENDIF ; ENDDO
WRITE(*,'(A85)')'=========================================================================================='!
t(1:j) = t_cv(1:j) !; PRINT*,t(1:j)    !# t_cv TASS CV INDEX
n1 = nbin(1) ;n2 = nbin(2) ;n3 = nbin(3) ;n4 = nbin(4) 

IF(mtd .eq. 'y') THEN
CALL mtd_unbiased(au_to_kcal,bias_fact,kt0,kt,kb,md_steps,mtd_steps,w_cv,w_hill &
                & ,n1,n2,ncv,t_min,t_max,cv,gridmin,gridmax,griddif,vbias,ct,nbin,m,u)
ENDIF
!---------------------------------------------------------------------------------------------------------------------------!
ALLOCATE(prob((nbin(u))))
IF (mtd .ne. 'y') m = t(1)
ALLOCATE(prob_2D(nbin(u),nbin(m)))
prob = 0.d0 ; prob_2D = 0.d0
IF (jj .eq. 0 .and. kk .eq. 0) jj = 1 ; kk = 1
ALLOCATE(prob_3D(nbin(ii),nbin(jj),nbin(kk)))
!---------------------------------------------------------------------------------------------------------------------------!
!Calling subroutine to unbaias the MTD potential
den = 0.d0  ; dum = 0.d0
IF(mtd .eq. 'y') THEN
CALL mtd_pot(md_steps,mtd_steps,w_cv,w_hill,t_min,t_max,gridmin,gridmax,griddif &
            & ,vbias,ct,m,u,ii,jj,kk,ncv,kt,nbin,cv,den,prob_2D)
END IF
!---------------------------------------------------------------------------------------------------------------------------!
!Computing 1-Dimensional Probability along Umbrella Coordinate
IF (Prob_nD .eq. 1) THEN
CALL oneD_prob(ii,u,m,mtd,t_min,t_max,md_steps,den,prob,prob_2D,ncv,cv,nbin,gridmin,gridmax,griddif)
DEALLOCATE(prob)
!---------------------------------------------------------------------------------------------------------------------------!
!Computing 2-Dimensional Probability along Umbrella/MTD (if MTD is enabled) or US/TEMP
ELSEIF (Prob_nd .eq. 2) THEN
IF (mtd .eq. 'y') THEN
CALL twoD_prob(ii,jj,u,m,kt,w_cv,w_hill,ct,vbias,t_min,t_max,md_steps,den,prob_2D,ncv,cv,nbin,gridmin,gridmax,griddif)
ELSEIF (mtd .eq. 'n') THEN
CALL twoD_temp_prob(ii,jj,u,kt,t_min,t_max,md_steps,prob_2D,ncv,cv,nbin,gridmin,gridmax,griddif)
DEALLOCATE(prob_2D)
ENDIF
ENDIF
!---------------------------------------------------------------------------------------------------------------------------!
DEALLOCATE(ct,cv,dummy,vbias,nbin,indx)
END PROGRAM WSMTD_rw_2D
!================================================== US Probability =========================================================!
!Compute 1 Dimensional Probability along Umbrella Coordinate
!
SUBROUTINE oneD_prob(ii,u,m,mtd,t_min,t_max,md_steps,den,prob,prob_2D,ncv,cv,nbin,gridmin,gridmax,griddif)
IMPLICIT NONE
INTEGER :: ii,u,m,i_md,t_min,t_max,index1,md_steps,i_s1,i_s2,ncv,nbin(*)
REAL*8  :: den,dum,s1,prob(*),gridmin(*),gridmax(*),griddif(*),cv(ncv,*),prob_2D(nbin(u),nbin(m))
CHARACTER*5 :: mtd

IF (mtd .eq. 'y') THEN
DO i_s1 = 1,nbin(u)
  dum = 0.0d0
   DO i_s2 = 1,nbin(m)
     dum = dum + prob_2D(i_s1,i_s2)
   ENDDO
prob(i_s1) = dum*griddif(m)
den = den + dum*griddif(m)
ENDDO 
ENDIF

IF (ii .eq. u .and. mtd .eq. 'n') THEN
 DO i_md=1,md_steps
   IF ((i_md.ge.t_min).and.(i_md.le.t_max)) THEN
          index1 = nint((cv(u,i_md) - gridmin(u))/griddif(u)) +1
!----------------------------------------------------------------------------
IF (index1 .gt. nbin(u)) THEN
PRINT*, "            ***************************************************"
PRINT*, "                ERROR !! UMBRELLA CV RANGE IS NOT CORRECT"
PRINT*, "            ***************************************************"
STOP ; ENDIF
!----------------------------------------------------------------------------
    prob(index1) = prob(index1) + 1.d0
   END IF
 END DO
ENDIF 

DO index1 = 1,nbin(u)
    den = den + prob(index1)
END DO   

OPEN(2,FILE='Pu.dat',STATUS='unknown')
DO i_s1 = 1,nbin(u)
     s1 = DFLOAT(i_s1-1)*griddif(u) + gridmin(u)
     IF (mtd .eq. 'y') THEN
     WRITE(2,*)s1,prob(i_s1)/(den*griddif(u)*griddif(m))
     ELSE
     WRITE(2,*)s1,prob(i_s1)/(den*griddif(u))
     ENDIF 
END DO
WRITE(*,'(A)')"Unbiased 1D distribution along US  written in 'Pu.dat'"
END SUBROUTINE 
!================================================== US vs MTD Probability ================================================!
!Compute Unbiased Distribution Along Umbrella Coordinate vs MTD coordinate
!
SUBROUTINE twoD_prob(ii,jj,u,m,kt,w_cv,w_hill,ct,vbias,t_min,t_max,md_steps,den,prob_2D,ncv,cv,nbin,gridmin,gridmax,griddif)

IMPLICIT NONE
INTEGER :: ii,jj,u,m,i_md,i_mtd,t_min,t_max,index1,md_steps,i_s1,i_s2,ncv,nbin(*),w_cv,w_hill,indx(2)
REAL*8  :: dum,den,s1,s2,kt
REAL*8  :: prob_2D(nbin(u),nbin(m)),gridmin(*),gridmax(*),griddif(*),cv(ncv,*),vbias(*),ct(*)

IF (ii .eq. u .and. jj .eq. m ) THEN
DO i_md=1,md_steps
    IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
     indx(u) = nint((cv(u,i_md)-gridmin(u))/griddif(u)) + 1
     indx(m) = nint((cv(m,i_md)-gridmin(m))/griddif(m)) + 1
!----------------------------------------------------------------------------
IF (indx(u) .gt. nbin(u) .or. indx(m) .gt. nbin(m)) THEN
PRINT*, "            ***************************************************"
PRINT*, "                ERROR !! ONE OF THE CV RANGE IS NOT CORRECT"
PRINT*, "            ***************************************************"
STOP ; ENDIF
!----------------------------------------------------------------------------
       i_mtd = i_md*w_cv/w_hill
         dum = vbias(i_md) - ct(i_mtd)
         dum = dexp(dum/kt)
     prob_2D(indx(u),indx(m)) = prob_2D(indx(u),indx(m)) + dum
         den = den + dum
    END IF
END DO
PRINT*,i_mtd
dum = den*griddif(u)*griddif(m) ;  den = 1.D0/dum
OPEN(2,FILE='Pu_2D.dat',STATUS='unknown')
DO i_s1 = 1,nbin(u)
          s1 = DFLOAT(i_s1-1)*griddif(u)+gridmin(u)
  DO i_s2 = 1,nbin(m)
          s2 = DFLOAT(i_s2-1)*griddif(m)+gridmin(m)
          prob_2D(i_s1,i_s2) = prob_2D(i_s1,i_s2)*den
    WRITE(2,'(3E16.8)')s1,s2,prob_2D(i_s1,i_s2)
  END DO
    WRITE(2,*)
END DO
WRITE(*,'(A)')'Unbiased 2D distribution along US vs MTD  written in Pu_2D.dat'
ENDIF 
CLOSE(2) 
END SUBROUTINE
!=================================================== US vs TEMP Probability ================================================!
!Compute Unbiased Distribution along Umbrella Coordinate vs One of the Temperature CV
!
SUBROUTINE twoD_temp_prob(ii,jj,u,kt,t_min,t_max,md_steps,prob_2D,ncv,cv,nbin,gridmin,gridmax,griddif)
IMPLICIT NONE
INTEGER :: ii,jj,u,t,i_md,t_min,t_max,md_steps,i_s1,i_s2,ncv,nbin(*),indx(2),index1,index2
REAL*8  :: dum,den,s1,s2,kt
REAL*8  :: prob_2D(nbin(u),nbin(jj)),gridmin(*),gridmax(*),griddif(*),cv(ncv,*)
den = 0.0d0 ; t = jj 
DO i_md=1,md_steps
    IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
     index1 = nint((cv(u,i_md)-gridmin(u))/griddif(u)) + 1
     index2 = nint((cv(t,i_md)-gridmin(t))/griddif(t)) + 1
!----------------------------------------------------------------------------
IF (index1 .gt. nbin(u) .or. index2 .gt. nbin(t)) THEN
PRINT*, "            ***************************************************" 
PRINT*, "                ERROR !! ONE OF THE CV RANGE IS NOT CORRECT" 
PRINT*, "            ***************************************************" 
STOP ; ENDIF
!----------------------------------------------------------------------------
     prob_2D(index1,index2) = prob_2D(index1,index2) + 1.0
    END IF
END DO

      DO index1=1,nbin(u)
        DO index2=1,nbin(t)
          den=den+prob_2D(index1,index2)
        END DO
      END DO

dum = den*griddif(u)*griddif(t) ;  den = 1.D0/dum
OPEN(2,FILE='Pu_2D.dat',STATUS='unknown')
DO i_s1 = 1,nbin(u)
        s1 = DFLOAT(i_s1-1)*griddif(u)+gridmin(u)
  DO i_s2 = 1,nbin(t)
        s2 = DFLOAT(i_s2-1)*griddif(t)+gridmin(t)
        prob_2D(i_s1,i_s2) = prob_2D(i_s1,i_s2)*den
    WRITE(2,'(3E16.8)')s1,s2,prob_2D(i_s1,i_s2)
  END DO
    WRITE(2,*)
END DO
WRITE(*,'(A)')'Unbiased 2D distribution along US vs TEMP  written in Pu_2D.dat'
CLOSE(2) 

ENDSUBROUTINE 
!============================================= MTD Unbias Potential ========================================================!
!Computing Unbiased Probability from WT-MTD Biased Potential
! 
SUBROUTINE mtd_pot(md_steps,mtd_steps,w_cv,w_hill,t_min,t_max,gridmin,gridmax,griddif &
                 &,vbias,ct,m,u,ii,jj,kk,ncv,kt,nbin,cv,den,prob_2D)
IMPLICIT NONE
REAL*8  :: kt,den,dum,ct(*),gridmin(*),gridmax(*),griddif(*),vbias(*)
INTEGER :: i_md,md_steps,i_mtd,mtd_steps,t_min,t_max,m,u,w_cv,w_hill,ii,jj,kk,ncv,indx(ncv),nbin(*)
REAL*8  :: cv(ncv,*),prob_2D(nbin(u),nbin(m))

den=0.d0 ; prob_2D=0.d0
IF (ii .eq. u .and. jj .eq. m ) THEN
!calculate prob (unbiased from MTD potential)      
DO i_md=1,md_steps
    IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
      indx(u) = nint((cv(u,i_md)-gridmin(u))/griddif(u)) + 1
      indx(m) = nint((cv(m,i_md)-gridmin(m))/griddif(m)) + 1
       i_mtd = i_md*w_cv/w_hill
         dum = vbias(i_md) - ct(i_mtd)
         dum = dexp(dum/kt)
     prob_2D(indx(u),indx(m)) = prob_2D(indx(u),indx(m)) + dum
         den = den + dum
    END IF
END DO
ENDIF
END SUBROUTINE
!---------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE get_steps(iunit,nsteps)
IMPLICIT NONE
INTEGER iunit, nsteps
INTEGER ios
nsteps=0
REWIND(iunit)
Read_Loop: DO
   READ(iunit,*,IOSTAT=ios)
   IF(ios.ne.0)EXIT Read_Loop
   nsteps=nsteps+1
END DO Read_Loop 
REWIND(iunit)
END 
!---------------------------------------------------------------------------------------------------------------------------!
!Compute WT-MTD unbiased potential 
!
SUBROUTINE mtd_unbiased(au_to_kcal,bias_fact,kt0,kt,kb,md_steps,mtd_steps,w_cv,w_hill &
                & ,n1,n2,ncv,t_min,t_max,cv,gridmin,gridmax,griddif,vbias,ct,nbin,m,u)
IMPLICIT NONE
INTEGER :: i,j,i_mtd,mtd_steps,i_md,md_steps,mtd_max,n1,n2,m,u
INTEGER :: i_s1,i_s2,w_cv,w_hill,ncv,index1,index2,t_max,t_min,nbin(ncv)
REAL*8  :: bias_fact,kt0,kt,kb,ktb,au_to_kcal,s1,s2
REAL*8  :: dummy1,diff_s2,ds2,ss,hh,num,den,alpha,dum
REAL*8,ALLOCATABLE             :: width(:),ht(:),hill(:)
REAL*8,DIMENSION(ncv,md_steps) :: cv
REAL*8,DIMENSION(md_steps)     :: vbias
REAL*8,DIMENSION(mtd_steps)    :: ct
REAL*8,DIMENSION(ncv)          :: gridmin,gridmax,griddif
REAL*8,DIMENSION(nbin(m))      :: grid,fes_2D
!REAL*8,DIMENSION(500,500)      :: prob_2D
                                                            !! #### m = METADYNAMICS CV INDEX
!OPEN(12,FILE='parvar_mtd',STATUS='unknown')
!OPEN(13,FILE='colvar_mtd',STATUS='unknown')
!CALL get_steps(13,mtd_steps)
ALLOCATE(width(mtd_steps),ht(mtd_steps),hill(mtd_steps))

IF (MOD(md_steps,mtd_steps) .ne. 0.0) THEN
PRINT*, "*****************************************************************" 
PRINT*, "         ERROR !! MISSING DATA IN colvar OR parvar FILES.        " 
PRINT*, "*****************************************************************" 
STOP ; ENDIF

bias_fact = (kt0 + bias_fact)/kt0
       kt = kb*kt
      ktb = kt*bias_fact
    alpha = bias_fact/(bias_fact-1.D0)

DO i_mtd=1,mtd_steps
  READ(12,*) dummy1,dummy1,width(i_mtd),ht(i_mtd)
  READ(13,*) dummy1,hill(i_mtd)
  ht(i_mtd) = ht(i_mtd)*au_to_kcal
END DO
!
DO i_s1 = 1,nbin(m)                                            ! Metadynamics CV bins
    grid(i_s1) = gridmin(m) + DFLOAT(i_s1 - 1)*griddif(m)
ENDDO
!
OPEN(21,FILE='ct.dat',STATUS='unknown')
DO i_mtd=1,mtd_steps
      ds2 = width(i_mtd)*width(i_mtd)
       ss = hill(i_mtd)
       hh = ht(i_mtd)
num = 0.D0; den = 0.D0 ; fes_2D = 0.0D0
   DO i_s1=1,nbin(m)  ! Metadynamics CV bins
      diff_s2 = grid(i_s1) - ss
!
      diff_s2 = diff_s2*diff_s2*0.5D0
      fes_2D(i_s1) = fes_2D(i_s1) - hh*DEXP(-diff_s2/ds2)
            num = num + DEXP(-fes_2D(i_s1)/kt)
            den = den + DEXP(-fes_2D(i_s1)/ktb)
  END DO
       ct(i_mtd) = kt*DLOG(num/den)
       WRITE(21,'(I10,F16.8)')i_mtd,ct(i_mtd)
END DO
CLOSE(21)
!!@calculate v(s,t)
 DO i_md = 1,md_steps
    mtd_max = (i_md*w_cv/w_hill)
         ss = cv(m,i_md)
dum = 0.d0
    DO i_mtd = 1,mtd_max
         ds2 = width(i_mtd)*width(i_mtd)
          hh = ht(i_mtd)/alpha
     diff_s2 = ss - hill(i_mtd)
!
     diff_s2 = diff_s2*diff_s2*0.5D0
         dum = dum + hh*DEXP(-diff_s2/ds2)
    END DO
     vbias(i_md) = dum
!PRINT*,i_md,vbias(i_md)
 END DO
!
!!!calculate prob (unbiased from MTD potential)       
!!den=0.d0 ; prob_2D=0.d0
!!
!!DO i_md=1,md_steps
!!    IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
!!      index1 = nint((cv(1,i_md)-gridmin(1))/griddif(1)) + 1
!!!      index2 = nint((cv(2,i_md)-gridmin(2))/griddif(2)) + 1
!!       i_mtd = i_md*w_cv/w_hill
!!         dum = vbias(i_md) - ct(i_mtd)
!!         dum = dexp(dum/kt)
!!     prob_2D(index1,index2) = prob_2D(index1,index2) + dum
!!         den = den + dum
!!PRINT*,i_mtd
!!    END IF
!!END DO
!!dum = den*griddif(1)*griddif(2) ;  den = 1.D0/dum
!!OPEN(2,FILE='Pu.dat',STATUS='unknown')
!!DO i_s1 = 1,n1
!!          s1 = DFLOAT(i_s1-1)*griddif(1)+gridmin(1)
!!  DO i_s2 = 1,n2
!!          s2 = DFLOAT(i_s2-1)*griddif(2)+gridmin(2)
!!          prob_2D(i_s1,i_s2) = prob_2D(i_s1,i_s2)*den 
!!    WRITE(2,'(3E16.8)')s1,s2,prob_2D(i_s1,i_s2)
!!  END DO
!!    WRITE(2,*)
!!END DO
!!WRITE(*,'(A)')'Unbiased distribution written in Pu.dat'
!!CLOSE(2)

END SUBROUTINE mtd_unbiased
!---------------------------------------------------------------------------------------------------------------------------!
