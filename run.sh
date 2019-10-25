#---------------------------------------------------------------------------------------------------------------#
#!/bin/bash
#Run file written by Rahul Verma to Unbias Probability Distribution along 1D and 2D 
#Date : October 2019
#---------------------------------------------------------------------------------------------------------------#
gfortran  -fcheck=all Probability_analysis.F90 -o Probability_analysis.x

#$#-------------------------------------------------------------------------------------------------------------#
./Probability_analysis.x  ,\
 -T0 300                  ,\
 -T 1000                  ,\
 -bias_fact 1500          ,\
 -tmin 1000               ,\
 -ncv 4                   ,\
 -UCV 1                   ,\
 -MTD n                   ,\
 -MCV 0                   ,\
 -Prob_nD 2               ,\
 -CV_num 1 3 3            ,\
 -grid -1.5 3.0 0.02 1.0 10.0 0.02 1.0 9.0 0.02 3.0 5.0 0.02 ,\
 -pfrqMD 1 ,\
 -dtMTD 200

#$#-------------------------------------------------------------------------------------------------------------#
#$./Probability_analysis.x ,\                                      # Executable                                 #                        
#$ -T0 300                  ,\                                     # Physical System Temeprature                #
#$ -T 1000                  ,\                                     # Extended CV Temperature                    #
#$ -bias_fact 1500          ,\                                     # Biased Factor WT-MTD                       #
#$ -tmin 1000               ,\                                     # Minimum MD Steps                           #
#$ -ncv 4                   ,\                                     # Total CV's in Simulation                   #
#$ -UCV 1                   ,\                                     # Umbrella CV                                #
#$ -MTD n                   ,\                                     # Metadynamics Enabling (y/n){Case-sensitive}#
#$ -MCV 0                   ,\                                     # MTD CV Coordinate Index in cvmdck_mtd File #
#$ -Prob_nD 2               ,\                                     # Probability Dimension (Output)             #
#$ -CV_num 1 3 3            ,\                                     # Index of Output Probability                #
#$ -grid -1.5 3.0 0.02 1.0 10.0 0.02 1.0 9.0 0.02 3.0 5.0 0.02 ,\  # Grid Size of All CV's From Simulation      #
#$ -pfrqMD 1 ,\                                                    # Print Frequency in 'cvmdck_mtd' File       #
#$ -dtMTD 200                                                      # Print Frequency in 'colvar_mtd' File       #
#$#-------------------------------------------------------------------------------------------------------------#
