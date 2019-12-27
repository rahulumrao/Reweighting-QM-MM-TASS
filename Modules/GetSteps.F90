MODULE GetSteps
CONTAINS
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
END MODULE GetSteps
