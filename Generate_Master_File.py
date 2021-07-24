'''This file contains functions to generate a "MasterFile" that is subsequently used to create the A and B files required
    to do the IPCW survival analysis.

    Name: Nicholas Wood, PhD
    Institution: USNA
    Email: nwood@usna.edu
'''

import pandas as pd
import numpy as np


def CalcMELD(Status):
    '''calculates MELD and return the column containing that score.'''

    #Copy the dataframe
    df = Status.copy()

    #Take the natural log of the relevant laboratory values
    df['ln_creat'] = np.log(df['canhx_creat_bound'].values)
    df['ln_inr'] = np.log(df['canhx_inr_bound'].values)
    df['ln_bili'] = np.log(df['canhx_bili_bound'].values)

    #Calculate untruncated MELD
    df['MELD'] = 10*(0.957*df['ln_creat'].values + 0.378*df['ln_bili'].values + 1.120*df['ln_inr'].values + 0.643).round(decimals = 1)

    #Constrain MELD to be between 6 and 40
    df['MELD'] = np.select([df['MELD'] < 6, df['MELD'] > 40],
                           [6, 40],
                           default = df['MELD'])

    #Return the column of MELD values
    return df['MELD'].values


def CalcMELDNa(Status):
    '''calculates MELD-Na and return the column containing that score'''

    #Copy the dataframe
    df = Status.copy()

    #First calculate basic MELD, which is needed to determine if the MELD-Na
    #definition comes into play
    df['MELD'] = CalcMELD(df)

    #Constrain the serum sodium values
    df['canhx_serum_sodium'] = np.select([df['canhx_serum_sodium'] < 125, df['canhx_serum_sodium'] > 137],
                                         [125, 137],
                                         default = df['canhx_serum_sodium'])

    #Calculate MELD-Na - by construction it is already between 6 and 40
    df['MELD-Na'] = np.where(df['MELD'] > 11, (df['MELD'] + 1.32*(137 - df['canhx_serum_sodium']) - 0.033*df['MELD']*(137 - df['canhx_serum_sodium'])).round(decimals = 0), df['MELD'])

    #Return the column of MELD-Na values
    return df['MELD-Na'].values



def CalcMELD3_With_Albumin(Status):
    '''Calculate MELD3 without Albumin'''

    #Copy the dataframe
    df = Status.copy()

    #Bound creatinine
    df['canhx_creat_bound'] = np.select([df['canhx_creat_bound'] < 1, df['canhx_creat_bound'] > 3],
                                        [1, 3],
                                        default = df['canhx_creat_bound'])

    #Bound serum sodium
    df['canhx_serum_sodium'] = np.select([df['canhx_serum_sodium'] < 125, df['canhx_serum_sodium'] > 137],
                                         [125, 137],
                                         default = df['canhx_serum_sodium'])

    #Bound albumin
    df['canhx_albumin'] = np.select([df['canhx_albumin'] < 1.5, df['canhx_albumin'] > 3.5],
                                    [1.5, 3.5],
                                    default = df['canhx_albumin'])

    #Grab relevant variables/transformation of variables for the calculation of
    #MELD3 with albumin
    df['female'] = (df['can_gender'] == 'F').astype(int)
    df['ln_creat'] = np.log(df['canhx_creat_bound'].values)
    df['ln_inr'] = np.log(df['canhx_inr_bound'].values)
    df['ln_bili'] = np.log(df['canhx_bili_bound'].values)

    #calculate MELD3 with albumin
    df['MELD 3 W Alb'] = (1.33*df['female'] + 4.56*df['ln_bili'] + 0.82*(137 - df['canhx_serum_sodium']) - 0.24*(137 - df['canhx_serum_sodium'])*df['ln_bili']\
                           + 9.09*df['ln_inr'] + 11.14*df['ln_creat'] + 1.85*(3.5 - df['canhx_albumin']) - 1.83*(3.5 - df['canhx_albumin'])*df['ln_creat'] + 6).round(decimals = 0)

    #Bound it between 6 and 40
    df['MELD 3 W Alb'] = np.select([df['MELD 3 W Alb'] < 6, df['MELD 3 W Alb'] > 40],
                                   [6, 40],
                                   default = df['MELD 3 W Alb'])

    #Return the column of MELD3 with albumin values
    return df['MELD 3 W Alb'].values


def CalcMELD3_Without_Albumin(Status):
    '''Calculate MELD3.0 without Albumin'''

    #Copy the dataframe
    df = Status.copy()

    #Bound creatinine
    df['canhx_creat_bound'] = np.select([df['canhx_creat_bound'] < 1, df['canhx_creat_bound'] > 3],
                                        [1, 3],
                                        default = df['canhx_creat_bound'])

    #Bound serum sodium
    df['canhx_serum_sodium'] = np.select([df['canhx_serum_sodium'] < 125, df['canhx_serum_sodium'] > 137],
                                         [125, 137],
                                         default = df['canhx_serum_sodium'])

    #Grab relevant variables/transformation of variables for the calculation of
    #MELD3 with albumin
    df['female'] = (df['can_gender'] == 'F').astype(int)
    df['ln_creat'] = np.log(df['canhx_creat_bound'].values)
    df['ln_inr'] = np.log(df['canhx_inr_bound'].values)
    df['ln_bili'] = np.log(df['canhx_bili_bound'].values)


    #calculate MELD3 without albumin
    df['MELD 3 WO Alb'] = (1.40*df['female'] + 4.85*df['ln_bili'] + 0.88*(137 - df['canhx_serum_sodium']) - 0.25*(137 - df['canhx_serum_sodium'])*df['ln_bili']\
                              + 9.66*df['ln_inr'] + 10.47*df['ln_creat'] + 6).round(decimals = 0)

    #Bound it between 6 and 40
    df['MELD 3 WO Alb'] = np.select([df['MELD 3 WO Alb'] < 6, df['MELD 3 WO Alb'] > 40],
                                    [6, 40],
                                    default = df['MELD 3 WO Alb'])
    
    #Return the column of MELD3 without albumin values
    return df['MELD 3 WO Alb'].values




def CalcBun(Status):
    '''Estimate blood urea nitrogen (BUN) according to the formula in the Asrani (2019 or 2020) supporting information.
        This is needed to calculate GRAIL and MELD-GRAIL-Na'''

    #Copy the dataframe
    df = Status.copy()

    #bound albumin
    df['canhx_albumin'] = np.where(df['canhx_albumin'] < 1, 1, df['canhx_albumin'])

    #Grab relevant variables/transformation of variables
    df['female'] = df['can_gender'] == 'F'
    df['black'] = df['can_race'] == '          16: Black or African American'
    df['ln_age'] = np.log(df['can_age'].values)
    df['ln_creat'] = np.log(df['canhx_creat_bound'].values)
    df['ln_alb'] = np.log(df['canhx_albumin'].values)

    #Calculate BUN
    df['Bun'] = np.exp(2.1329 + 0.0749*df['female'] - 0.1370*df['black'] + 0.1259*df['ln_age'] + 0.9774*df['ln_creat'] + 0.1794*df['ln_alb'])

    #Return the column of BUN
    return df['Bun'].values



def CalcTrigger(Status):
    '''Trigger value used in the calculation of GRAIL and MELD-GRAIL-Na'''

    #Copy the dataframe
    df = Status.copy()

    #Calculate BUN
    df['Bun'] = CalcBun(df)

    #Bound albumin
    df['canhx_albumin'] = np.where(df['canhx_albumin'] < 1, 1, df['canhx_albumin'])

    #Grab relevant variables/transformation of variables
    df['ln_bun'] = np.log(df['Bun'].values)
    df['female'] = df['can_gender'] == 'F'
    df['black'] = df['can_race'] == '          16: Black or African American'
    df['ln_age'] = np.log(df['can_age'].values)
    df['ln_creat'] = np.log(df['canhx_creat_bound'].values)
    df['ln_alb'] = np.log(df['canhx_albumin'].values)

    #Calculate an intermediate value
    df['logit'] = -7.281 + 0.793*df['female'] - 0.730*df['black'] + 0.526*df['ln_age'] + 3.706*df['ln_creat'] + 1.263*df['ln_bun'] - 2.622*df['ln_alb']

    #Get from that intermediate value a probability
    df['P'] = 1/(1+np.exp(-1*df['logit']))

    #Trigger = 0 if P < 0.05 and 1 otherwise
    df['Trigger'] = np.where(df['P'] <= 0.05, 0, 1)

    #Return a column of the trigger values
    return df['Trigger'].values


def CalcGRAIL(Status):
    '''Calculate the GFR asessment in liver disease (GRAIL)'''

    #Copy the dataframe
    df = Status.copy()

    #Calculate BUN
    df['Bun'] = CalcBun(df)

    #Calculate Trigger
    df['Trigger'] = CalcTrigger(df)

    #Grab candidate sex and race
    df['female'] = df['can_gender'] == 'F'
    df['black'] = df['can_race'] == '          16: Black or African American'

    
    #Calculate GRAIL
    df['GRAIL'] = np.where(df['Trigger'] == 0,
            353.435*(~df['female'] + 0.871*df['female'])*(~df['black'] + 1.069*df['black'])*df['can_age']**(-0.318)*df['canhx_creat_bound']**(-0.402)*df['Bun']**(-0.068)*df['canhx_albumin']**(0.123),
            5.253*(~df['female'] + 0.917*df['female'])*(~df['black'] + 0.807*df['black'])*df['can_age']**(0.447)*df['canhx_creat_bound']**(-0.466)*df['Bun']**(-0.042)*df['canhx_albumin']**(0.224))

    #Return a column of GRAIL values
    return df['GRAIL'].values




def CalcMELDGRAILNa(Status):
    '''Calculates MELD-GRAIL-Na'''

    #Copy the dataframe
    df = Status.copy()

    #Calculate GRAIL
    df['GRAIL'] = CalcGRAIL(df)

    #Set eGFR to 15 if on dialysis. Otherwise set it to GRAIL
    df['eGFR'] = np.where(df['canhx_dial_prior_week'] == 'Y', 15, df['GRAIL'])

    #Bound inr
    df['canhx_inr_bound'] = np.select([df['canhx_inr'] < 1, df['canhx_inr'] > 3],
                                      [1, 3],
                                      default = df['canhx_inr'])

    #Bound bilirubin
    df['canhx_bili_bound'] = np.where(df['canhx_bili'] < 1, 1, df['canhx_bili'])

    #bound eGFR
    df['eGFR'] = np.select([df['eGFR'] < 15, df['eGFR'] > 90],
                           [15, 90],
                           default = df['eGFR'])

    #Bound sodium
    df['canhx_serum_sodium'] = np.select([df['canhx_serum_sodium'] < 125, df['canhx_serum_sodium'] > 140],
                                         [125, 140],
                                         default = df['canhx_serum_sodium'])

    #Transform variables
    df['ln_inr'] = np.log(df['canhx_inr_bound'])
    df['ln_bili'] = np.log(df['canhx_bili_bound'])
    df['ln_eGFR'] = np.log(df['eGFR'])
    df['ln_Na'] = np.log(df['canhx_serum_sodium'])

    #Calculate MELD-GRAIL-Na
    df['MELD-GRAIL-Na'] = 29.751 + 10.836*df['ln_inr'] + 3.039*df['ln_bili'] - 5.054*df['ln_eGFR'] - 0.372*df['ln_Na']


    #This step is optional and could be removed or altered. To get MELD-GRAIL-Na to have
    #a similar distribution to MELD-Na, it needs to be scaled. Below is what I did.
    slope = (50 - 6)/(df['MELD-GRAIL-Na'].max() - df['MELD-GRAIL-Na'].min())
    df['MELD-GRAIL-Na'] = (slope*(df['MELD-GRAIL-Na'] - df['MELD-GRAIL-Na'].min())) + 6
                           
    
    #Bound MELD-GRAIL-Na to be between 6 and 40
    df['MELD-GRAIL-Na'] = np.select([df['MELD-GRAIL-Na'] < 6, df['MELD-GRAIL-Na'] > 40],
                                    [6, 40],
                                    default = df['MELD-GRAIL-Na'])



    #Return the column of MELD-GRAIL-Na values
    return df['MELD-GRAIL-Na'].values



def GetMatchMELD(Status):
    '''Gets the match MELD score. I THINK canhx_optn_lab_meld is correct for this
        purpose, but I could be wrong. This is redunant, but for status 1A and 1B
        candidates I set match MELD to 41 (I do this later on in the logic already).'''
    
    #Copy the dataframe
    df = Status.copy()

    #Restructure the canhx_optn_lab_meld
    df['canhx_optn_lab_meld'] = df['canhx_optn_lab_meld'].str.lstrip().str[0:4]

    #Determine which status updates are MELD
    df['meld_status'] = df['canhx_optn_lab_meld'].str[0:2].isin(['61', '62'])

    #Determine which status updates are status 1
    df['status1'] = df['canhx_optn_lab_meld'].isin(['6011', '6012'])

    #only look at meld/status1 updates
    df = df.loc[(df['meld_status']) | (df['status1'])]

    #Now determinine match meld
    df['untrunc_match_meld'] = np.select([(df['status1']), (df['canhx_optn_lab_meld'].str[1] == '1'), (df['canhx_optn_lab_meld'].str[1] == '2')],
                                         [41, df['canhx_optn_lab_meld'].str[2:].astype(int) - 100, df['canhx_optn_lab_meld'].str[2:].astype(int)])

    #Bound match MELD to be between 6 and 40
    df['trunc_match_meld'] = np.select([(df['untrunc_match_meld'] < 6), (df['untrunc_match_meld'] > 40)],
                                       [6, 40],
                                       default = df['untrunc_match_meld'])

    #Incorporate status 1A/1B
    df['match_meld'] = np.where(df['status1'], 41, df['trunc_match_meld'])


    #Return the column of match MELD values
    return df['match_meld']

    

def GenerateMasterSurvivalDataSet(start_date, end_date, folder, meld_cols,
                                  cand_liin_path = r'C:\Users\Nicholas\Documents\SRTR Data\Liver\cand_liin.dta',
                                  stathist_liin_path = r'C:\Users\Nicholas\Documents\SRTR Data\Liver\stathist_liin.dta'):
    
    '''Generates a master survival data set to be used for survival analysis.

        start_date: String, format %m-%d-%Y. Beginning of the study.
        
        end_date: String, format %m-%d-%Y. End of the study.

        folder: String. Folder wherein the resulting master file is saved.

        meld_cols: List. Contains a list of names of MELD definitions you want calculated.

        cand_liin_path: String. Path to the cand_liin.dta SRTR SAF.

        stathist_liin_path: String. Path to the stathist_liin.dta SRTR SAF.


        Returns the resulting dataframe to be used by functions in the Convert_Master_File_to_A_and_B_Files.
    '''


    #Read in the candidate file
    Candidates = pd.read_stata(cand_liin_path, index_col = 'px_id')

    #Read in the status file
    Status = pd.read_stata(stathist_liin_path)

    #Find all status updates that begin between the start and end dates
    Status = Status.loc[(Status['canhx_begin_dt'] >= start_date) & (Status['canhx_begin_dt'] < end_date)]

    #Add information to the Status dataframe that is relevant - i.e. age at listing, sex, etc.
    for col in ['can_age_in_months_at_listing', 'can_gender', 'can_race', 'can_dgn', 'can_rem_cd', 'can_rem_cod', 'can_rem_dt', 'can_listing_dt', 'rec_tx_dt']:
        Status[col] = Status['px_id'].map(Candidates[col])


    #Add Status Information
    Status['stat'] = Status['canhx_stat_cd'].astype(str)
    Status['stat'] = np.select([(Status['stat'].str.contains('Status 1B') | Status['stat'].str.contains('Status 2B')),
                                (Status['stat'].str.contains('Status 1A') | Status['stat'].str.contains('Status 2A'))],
                               ['1B', '1A'],
                               default = '0')


    #Calculate candidate age in years - Specifically needed for GRAIL
    Status['can_age'] = Status['can_age_in_months_at_listing']/12 + (Status['canhx_begin_dt'] - Status['can_listing_dt']).dt.days/365/25


    #Store MELD name to the func that calculates that MELD. Update this as needed. In theory I could
    #Change this to be an input to this function.
    meld2meldfunc = {'MELD': CalcMELD,
                     'MELD-Na': CalcMELDNa,
                     'MELD 3 W Alb': CalcMELD3_With_Albumin,
                     'MELD 3 WO Alb': CalcMELD3_Without_Albumin,
                     'MELD-GRAIL-Na': CalcMELDGRAILNa,
                     'Match MELD': GetMatchMELD}

    #Calculate each definition of MELD
    for meld_col in meld_cols:
        Status[meld_col] = meld2meldfunc[meld_col](Status)

        #For each 1A or 1B candidate, set the meld_col to 41
        Status[meld_col] = np.where(Status['stat'].isin(['1A', '1B']), 41, Status[meld_col])

        


    #Select which columns we want to keep.
    columns = ['px_id', 'can_age_in_months_at_listing', 'can_gender', 'can_race', 'can_dgn', 'can_rem_cd', 'can_rem_cod', 'can_rem_dt', 'rec_tx_dt', 'can_listing_dt',
               'canhx_exc_flg', 'stat', 'canhx_begin_dt', 'canhx_end_dt'] + meld_cols
    


    #Save the resulting Master Data Set
    Status[columns].to_csv(f'{folder}\\Master Data Set {start_date} to {end_date}.csv', index = False)

    #Return the dataframe.
    return Status
    




    
