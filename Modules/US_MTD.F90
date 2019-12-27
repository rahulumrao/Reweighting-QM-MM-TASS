MODULE US_MTD
CONTAINS
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
!PRINT*,i_mtd
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
END MODULE US_MTD
