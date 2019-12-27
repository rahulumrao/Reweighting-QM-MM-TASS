MODULE MTD_Potential
CONTAINS
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
END SUBROUTINE mtd_pot
END MODULE MTD_Potential
