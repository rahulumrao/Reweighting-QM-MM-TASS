MODULE US_TEMP
CONTAINS
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

END SUBROUTINE
END MODULE US_TEMP
