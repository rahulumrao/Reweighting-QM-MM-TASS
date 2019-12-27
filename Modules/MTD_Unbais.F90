MODULE MTD_Unbais
CONTAINS
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
!! END DO
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
END
END MODULE MTD_Unbais
!---------------------------------------------------------------------------------------------------------------------------!
