MODULE US_Prob
CONTAINS
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
END MODULE US_Prob
