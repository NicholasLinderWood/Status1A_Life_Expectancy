# Status1A_Life_Expectancy
Within this repository will be contained all of the code (described below) I used for the analysis in my "Life expectancy without a transplant of status 1A liver transplant candidates" (Currently undergoing revisions for publication in the American Journal of Transplantation). To do the analysis, I assume the user of this code will have access to the following standard analysis files from the Scientific Registry of Transplant Recipients:  

cand_liin.dta, stathist_liin.dta, statjust_li1.dta


DESCRIPTION OF FILES:

Generate_A_and_B_Files.py - Script that creates two data sets. An "A" file and a "B" file. These files are formatted in such a way to be used for Kaplan-Meier/IPCW survival analysis. See the file for additional details.

Generate_Master_File.py - Used by Generate_A_and_B_Files.py

Convert_Master_File_to_A_and_B_Files.py - Used by Generate_A_and_B_Files.py

RunStatus1ALifeExpectancyAnalysis.py - Script that runs all of the analysis, and creates all tables and figures, for this paper. See the file for additional details.

Status1ALifeExpectancyAnalysisFunctions.py - Script that contains all the functions for the aforementioned analysis. Used by RunStatus1ALifeExpectancyAnalysis.py

IPCW.py - Contains the class definition for the IPCS survival estimator. Used by RunStatus1ALifeExpectancyAnalysis.py

KaplanMeier.py - Contains the class definition for KaplanMeier survival estimator. Used by RunStatus1ALifeExpectancyAnalysis.py and IPCW.py


