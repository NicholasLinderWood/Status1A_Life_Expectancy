'''This script combines functions from Generate_Master_File and Conver_Master_File_to_A_and_B_Files
    to create the data sets needed for the IPCW survival calculations.

    Name: Nicholas Wood, PhD
    Institution: USNA
    Email: nwood@usna.edu
'''


from Generate_Master_File import GenerateMasterSurvivalDataSet
from Convert_Master_File_to_A_and_B_Files import GenerateFiles
import pandas as pd
import os

##########################################################################################################
#USER INPUTS
##########################################################################################################
#These can be changed as desired.
##########################################################################################################

#Define the start and end dates of the study. The format MUST be %m-%d-%Y. The longer the time frame
#the longer it will take to generate these files. For testing purposes I suggest only using 1 year or
#shorter.
start_date = '01-01-2010'
end_date = '01-01-2020'

#Provide the name of the folder in which you want the resultant data sets to be saved.
#If this folder does not already exist it will be created.
folder = f'Survival Data {start_date} to {end_date}'

#Define which meld columns you want. The more MELD definitions you choose the longer
#it will take to generate the files. For testing purposes I suggest only using 1 MELD definition.
meld_cols = ['Match MELD']

#Define which candidates you want to include/exlcude based on age, exception, and status
age_min = 18
age_max = 120
exc_cands = 'Exclude'
status1A_cands = 'Include'
status1B_cands = 'Exclude'

#Define any extra columns you want included in the A and B files
extra_cols = []

#Define the paths for the cand_liin.dta and stathist_liin.dta SRTR SAFs. The below are where
#my SAFs are, and for you this will be different.
cand_liin_path = r'C:\Users\Nicholas\Documents\SRTR Data\Liver\cand_liin.dta'
stathist_liin_path = r'C:\Users\Nicholas\Documents\SRTR Data\Liver\stathist_liin.dta'

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################







#If the specified folder does not already exist, create it.
try:
    os.mkdir(folder)
    
except FileExistsError:
    pass


#Generate the Master File
df = GenerateMasterSurvivalDataSet(start_date, end_date, folder, meld_cols, cand_liin_path, stathist_liin_path)

#Generate the A and B files
GenerateFiles(df, folder, end_date, meld_cols, extra_cols, age_min = age_min, age_max = age_max, exc_cands = exc_cands, status1A_cands = status1A_cands, status1B_cands = status1B_cands)













