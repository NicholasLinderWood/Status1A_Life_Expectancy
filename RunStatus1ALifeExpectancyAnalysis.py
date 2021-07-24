'''This script runs all the analysis for the status 1A life expectancy paper. The A and B files must have
    been generated before hand.

    Name: Nicholas Wood, PhD
    Institution: USNA
    Email: nwood@usna.edu
'''
from Status1ALifeExpectancyAnalysisFunctions import Add1AInfo, GetStatus1ACohortTable, Bootstrap
from Status1ALifeExpectancyAnalysisFunctions import CalculateLifeExpectancy, Calculate1ASurvival, Generate1ASurvivalFigures


##########################################################################################################
#USER INPUTS
##########################################################################################################
#These can be changed as desired.
##########################################################################################################

#Define the A and B file paths
AFile_path = r'C:/Users/Nicholas/Documents/Survival/Without Transplant Survival Data Generator/Survival Data 01-01-2010 to 01-01-2020/MELD-Na A File.csv'
BFile_path = r'C:/Users/Nicholas/Documents/Survival/Without Transplant Survival Data Generator/Survival Data 01-01-2010 to 01-01-2020/MELD-Na B File.csv'

#Define the paths of the requisite SRTR SAFs
stathist_liin_path = r"C:\Users\Nicholas\Documents\SRTR Data\Liver\stathist_liin.dta"
statjust_li1_path = r"C:\Users\Nicholas\Documents\SRTR Data\Liver\statjust_li1.dta"
cand_liin_path = r'C:/Users/Nicholas/Documents/SRTR Data/Liver/cand_liin.dta'

#Define the start and end dates
start_date = '01-01-2010'
end_date = '01-01-2020'

#Define the number of bootstrap iterations
n = 200

#Define the amount of time to do the survival calculation
T = 90

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################




#Add 1A info to the A and B files
dfA, dfB = Add1AInfo(AFile_path, BFile_path, start_date, stathist_liin_path, statjust_li1_path)

#Create the status 1A cohort table
GetStatus1ACohortTable(dfA, end_date, cand_liin_path)

#Calculate survival for 1A overall and categories
Calculate1ASurvival(dfA, dfB, T = T)

#Calculate Life Expectancy
LE = CalculateLifeExpectancy()

#Generate the 1A survival figures
Generate1ASurvivalFigures()

#Do the bootstrapping!
dfB = Bootstrap(dfA, dfB, LE, n = n, T = T)














