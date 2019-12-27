#Makefile

Probability_analysis.x  : GetSteps.o MTD_Unbais.o MTD_Potential.o US_Prob.o US_MTD.o US_TEMP.o
	gfortran -o Probability_analysis.x Main.F90 GetSteps.o MTD_Unbais.o MTD_Potential.o US_Prob.o US_MTD.o US_TEMP.o

GetSteps.o  :   GetSteps.F90
	gfortran -c GetSteps.F90

MTD_Unbais.o    :   MTD_Unbais.F90
	gfortran -c MTD_Unbais.F90

MTD_Potential.o :   MTD_Potential.F90
	gfortran -c MTD_Potential.F90

US_Prob.o   :   US_Prob.F90
	gfortran -c US_Prob.F90

US_MTD.o    :   US_MTD.F90
	gfortran -c US_MTD.F90

US_TEMP.o   :   US_TEMP.F90
	gfortran -c US_TEMP.F90

clean   :
	rm *.o Probability_analysis.x
