'''This file contains functions that take the "MasterFile", reduces it according to some exclusion criteria,
    and then uses it to create the A and B files needed to do the IPCW survival analysis.

    Name: Nicholas Wood, PhD
    Institution: USNA
    Email: nwood@usna.edu
'''

import pandas as pd
import numpy as np
from datetime import datetime


def ReduceMasterFile(df, meld_cols, end_date, age_min = 18, age_max = 120, exc_cands = 'Exclude', status1A_cands = 'Exclude', status1B_cands = 'Exclude'):
    '''Take the given master file, and reduce it by only considering candidates within certain criteria.

        df: pandas DataFrame of the Master File created with GenerateMasterSurvivalDataSet from Generate_Master_File.

        meld_cols: List of strings. Each string indicating a column in the Master File that is some definition of MELD.

        end_date: String of the format %m-%d-%Y. Should be identical to the end_date used to generate the master file.

        age_min: Float. The minimum age of candidates (at listing) that are included.

        age_max: Float. The maximum age of candidates (at listing) that are included.

        exc_cands: String. if exc_cands == 'Exclude', any candidate who ever had an exception is excluded.

        status1A_cands: String. if status1A_cands == 'Exclude', any candidate who ever achieved status 1A is excluded.

        status1B_cands: String. if status1B_cands == 'Exclude', any candidate who ever achieved status 1B is excluded.


        Returns the master file after reducing according to these exclusion criteria.
    '''

    #Keep only those candidates within the age range
    df = df.loc[(df['can_age_in_months_at_listing'] >= 12*age_min) & (df['can_age_in_months_at_listing'] < 12*age_max)]

    #Now check regarding exception candidates, status 1A cands, and status 1B cands. Based
    #On the input to the function we either include them or exclude them in our data set.
    if exc_cands == 'Exclude':
        px = df.loc[df['canhx_exc_flg'] == 1]['px_id'].unique()

        df = df.loc[~df['px_id'].isin(px)]

    if status1A_cands == 'Exclude':
        px = df.loc[df['stat'] == '1A']['px_id'].unique()

        df = df.loc[~df['px_id'].isin(px)]


    if status1B_cands == 'Exclude':
        px = df.loc[df['stat'] == '1B']['px_id'].unique()

        df = df.loc[~df['px_id'].isin(px)]


    #Finally, remove any rows for which th canhx_begin_dt == canhx_end_dt == rec_tx_dt. These are status updates
    #that were created when the candidate came in for transplant AFTER allocation occured. For that reason they are not
    #informative
    df = df.loc[~((df['canhx_begin_dt'] == df['canhx_end_dt']) & (df['canhx_begin_dt'] == df['rec_tx_dt']))]


    #for any removal date which is after the end date (or which is null) change it to the given end date
    df['can_rem_dt'] = np.select([pd.isnull(df['can_rem_dt']), df['can_rem_dt'] >= end_date],
                                 [pd.to_datetime(end_date, format ='%m-%d-%Y'), pd.to_datetime(end_date, format ='%m-%d-%Y')],
                                 default = df['can_rem_dt'].dt.date)

    #Similarly for the recipient transplant date - Make it so that those transplanted after the end date
    #(Or not at all) have a tx_dt that is NaT
    df['rec_tx_dt'] = np.select([pd.isnull(df['rec_tx_dt']), df['rec_tx_dt'] >= end_date],
                                 [pd.to_datetime(end_date, format ='%m-%d-%Y'), pd.to_datetime(end_date, format ='%m-%d-%Y')],
                                 default = df['rec_tx_dt'].dt.date)

    df['rec_tx_dt'] = pd.to_datetime(df['rec_tx_dt'], format ='%Y-%m-%d')


    #Likewise for the end date
    df['canhx_end_dt'] = np.select([pd.isnull(df['canhx_end_dt']), df['canhx_end_dt'] >= end_date],
                                 [pd.to_datetime(end_date, format ='%m-%d-%Y'), pd.to_datetime(end_date, format ='%m-%d-%Y')],
                                 default = df['canhx_end_dt'].dt.date)

    df['canhx_end_dt'] = pd.to_datetime(df['canhx_end_dt'], format ='%Y-%m-%d')



    #Let's also create a column for each kind of meld, indicating how those meld scores change for each candidate over time
    #for the purpose of grouping. In the creation of the A and B files, we merge status updates together if the MELD score
    #(for whichever definition) did not change from the first to the second. Furthermore, missing values of MELD are assumed
    #to be the MELD score of the previous status update. If no such MELD score exists, that candidates rows will ultimately be
    #dropped.

    #Group by px_id
    for px_id, tdf in df.groupby('px_id'):

        #Create a dictionary to store MELD groupings
        meldcol2group = {meld_col:0 for meld_col in meld_cols}

        #Sort by begin date of the status
        tdf = tdf.sort_values(by = 'canhx_begin_dt')

        #Store the "last meld". To start that will be the first.
        last_meld = {meld_col: tdf.iloc[0][meld_col] for meld_col in meld_cols}

        #Iterate over each row for this candidate
        for ind, row in tdf.iterrows():

            for meld_col in meld_cols:
                current_meld = row[meld_col]

                #If the current MELD is missing, fill it in with
                #the last known MELD score (if available)
                if current_meld == np.nan:
                    current_meld = last_meld[meld_col]
                    df.at[ind, meld_col] = last_meld[meld_col]
                    
                if current_meld == last_meld[meld_col]:
                    df.at[ind, f'{meld_col} Group'] = meldcol2group[meld_col]

                else:
                    last_meld[meld_col] = current_meld
                    meldcol2group[meld_col] += 1
                    df.at[ind, f'{meld_col} Group'] = meldcol2group[meld_col]
        

    #Return the reduced master file
    return df




def GenerateAFile(Status, meld_col, end_date, folder, extra_cols = []):
    '''Generate the A file for a given MELD score. For each candidate, track through his history.
        For each status update at a new MELD, record the days until this candidates MELD changed. At the
        end of that snippet record if the candidate was transplanted or not. The resulting A file is saved
        as a csv.

        Status: pandas DataFrame. This is the reduced master file returned from ReduceMasterFile

        meld_col: String. The MELD column used to generate the A file (e.g. MELD-Na, MELD-GRAIL-Na, etc.)
                    Must be a MELD column that was included in the creation of the master file.

        end_date: String of the format %m-%d-%Y. Should be identical to the end_date used to generate the master file
                    and the end_date used to reduce the master file.

        folder: String. The folder in which you want the created A file to be saved.

        extra_cols: List of strings. Each element should be a column in the master file. These are added to
                    the resulting A file so stratified survival analysis can be done if desired (e.g. can_gender).



        The resulting A file will be used in the IPCW survival estimator in order to calculate the weights. These
        weights correspond to the probability of remaining untransplanted.
        '''

    #Create a dataframe to store the A file
    df = pd.DataFrame(columns = ['px_id', 'Days', meld_col, 'Transplant'] + extra_cols)

    #Create an index counter variable for df
    ind = 0
    
    #Iterate over each candidate
    for px_id, tdf in Status.groupby('px_id'):

        #Copy the data frame because we will be adding a column, and we don't want
        #that copy warning
        tdf = tdf.copy()

        #Grab the transplant date
        rec_tx_dt = tdf.iloc[0]['rec_tx_dt']

        #Grab the removal code
        can_rem_cd = tdf.iloc[0]['can_rem_cd']

        #Get the candidate remove date
        can_rem_dt = tdf.iloc[0]['can_rem_dt']

        #Grab the extra cols
        extras = [tdf.iloc[0][col] for col in extra_cols]

        #Sort the status updates by the begin date
        tdf = tdf.sort_values(by = 'canhx_begin_dt')
                
        #Grab the last group. Only the last group could possibly be a transplant
        last_group = sorted(list(tdf[f'{meld_col} Group'].unique()))[-1]


        #Group the candidates status updates according to the groups defined when reducing the master file.
        for group, temp_df in tdf.groupby(f'{meld_col} Group'):

            #Ensure the status updates in this group are sorted by begin date.
            temp_df = temp_df.sort_values(by = 'canhx_begin_dt')

            #The number of days spent at this MELD score is the end of the last update minus the
            #beginning of the first, plus 1
            Days = (temp_df.iloc[-1]['canhx_end_dt'] - temp_df.iloc[0]['canhx_begin_dt']).days + 1

            #Determine the meld score for these status updates. It should be the same for all
            #status updates in this group, so simply grab the first.
            meld = temp_df.iloc[0][meld_col]

            #If this is not the last group, then it is not a transplant
            if group != last_group:
                Transplant = 0

            #Otherwise, it might be a transplant
            else:

                #If the recipient was transplanted before the end date, this counts as a transplant
                if rec_tx_dt < pd.to_datetime(end_date, format ='%m-%d-%Y'):
                    Transplant = 1

                #Otherwise that transplant does not count for this study (administrative censoring)
                else:
                    Transplant = 0


                #If this is the last group, then check to see if the next day was the transplant day.
                #If so, we need to add back in one day.
                if (temp_df.iloc[-1]['canhx_end_dt'] + pd.to_timedelta(1, unit = 'days') == rec_tx_dt):
                    Days += 1

                    #However, if the end date of the last status update is the same as the removal update,
                    #we do not need to add that day back in. This is convuluted, and there is probably a
                    #more elegant way to handle this issue than what I have done here. But this works.
                    if temp_df.iloc[-1]['canhx_end_dt'] == can_rem_dt:
                        Days -= 1
            

            #Add the information to our File A dataframe
            df.loc[ind] = [px_id, Days, meld, Transplant] + extras

            #Increment the index
            ind += 1
        
 
    #Finally, drop any rows with missing values.
    df = df.dropna()


    ##Save the A file as a csv to the given folder
    df.to_csv(f'{folder}\\{meld_col} A File.csv', index = False) 




def GenerateBFile(Status, meld_col, end_date, folder, extra_cols = []):
    '''Generate the B file for a given MELD score. For each candidate, track through his history.
        For each status update at a new MELD, record the days until this candidates MELD changed. At the
        end of that snippet record if the candidate died or not. The resulting B file is saved
        as a csv.

        Status: pandas DataFrame. This is the reduced master file returned from ReduceMasterFile

        meld_col: String. The MELD column used to generate the B file (e.g. MELD-Na, MELD-GRAIL-Na, etc.)
                    Must be a MELD column that was included in the creation of the master file.

        end_date: String of the format %m-%d-%Y. Should be identical to the end_date used to generate the master file
                    and the end_date used to reduce the master file.

        folder: String. The folder in which you want the created B file to be saved.

        extra_cols: List of strings. Each element should be a column in the master file. These are added to
                    the resulting B file so stratified survival analysis can be done if desired (e.g. can_gender).



        The resulting B file will be used in the IPCW survival estimator to calculate survival. Without the addition
        of the results from the A file, the B file contains data that could be used with a simple Kaplan-Meier survival
        estimator.
        '''

    #Create dataframe for B file
    df = pd.DataFrame(columns = ['px_id', 'Days', meld_col, 'Death'] + extra_cols)

    #Create an index counter variable for df
    ind = 0

    #Iterate over each candidate
    for px_id, tdf in Status.groupby('px_id'):

        #Grab removal date
        can_rem_dt = tdf.iloc[0]['can_rem_dt']

        #Grab removal code for determining whether this candidate died
        can_rem_cd = tdf.iloc[0]['can_rem_cd']

        #Grab cause of death - we will not consider suicide a death related to liver disease
        can_rem_cod = tdf.iloc[0]['can_rem_cod']

        #Grab the extra cols
        extras = [tdf.iloc[0][col] for col in extra_cols]

        #If the candidates removal date is on or after the end date, then we have administrative censoring
        #and by default the death column is zero
        if can_rem_dt >= pd.to_datetime(end_date, format ='%m-%d-%Y'):
            Death = 0

        #Otherwise, check the removal code. If it is 8, 13, or 5, we count this as a death
        #unless the cod was suicide
        else:
            if can_rem_cd in ['           8: Died',
                                '          13: Candidate condition deteriorated , too sick for tx',
                                '           5: Medically Unsuitable']:

                if can_rem_cod != '        4920: LI:SUICIDE:ATTEMPTED SUICIDE - DIED LATER':
                    Death = 1
                else:
                    Death = 0

            else:
                Death = 0

        #Group the candidates status updates according to the groups defined when reducing the master file.
        for group, temp_df in tdf.groupby(f'{meld_col} Group'):

            #Ensure the status updates in this group are sorted by begin date.
            temp_df = temp_df.sort_values(by = 'canhx_begin_dt')

            #Days are the number of days from the first of these status updates to the removal date
            #plus 1
            Days = (can_rem_dt - temp_df.iloc[0]['canhx_begin_dt']).days + 1

            #Determine the meld score for these status updates. It should be the same for all
            #status updates in this group, so simply grab the first.
            meld = temp_df.iloc[0][meld_col]
            
            #Add the information to our File B dataframe
            df.loc[ind] = [px_id, Days, meld, Death] + extras

            #Increment the index
            ind += 1

    #Drop any missing values
    df = df.dropna()
        
    #Save the B file as a csv to the given folder
    df.to_csv(f'{folder}\\{meld_col} B File.csv', index = False)

    




def GenerateFiles(df, folder, end_date, meld_cols, extra_cols = [], age_min = 18, age_max = 120, exc_cands = 'Exclude', status1A_cands = 'Exclude', status1B_cands = 'Exclude'):
    '''Simply runs all three of the above functions together for each given MELD column.

        df: pandas DataFrame of the Master File created with GenerateMasterSurvivalDataSet from Generate_Master_File.

        folder: String. The folder in which you want the created A and B files to be saved.

        end_date: String of the format %m-%d-%Y. Should be identical to the end_date used to generate the master file.

        meld_cols: List of strings. Each string indicating a column in the Master File that is some definition of MELD.

        extra_cols: List of strings. Each element should be a column in the master file. These are added to
                    the resulting A and B files so stratified survival analysis can be done if desired (e.g. can_gender).

        age_min: Float. The minimum age of candidates (at listing) that are included.

        age_max: Float. The maximum age of candidates (at listing) that are included.

        exc_cands: String. if exc_cands == 'Exclude', any candidate who ever had an exception is excluded.

        status1A_cands: String. if status1A_cands == 'Exclude', any candidate who ever achieved status 1A is excluded.

        status1B_cands: String. if status1B_cands == 'Exclude', any candidate who ever achieved status 1B is excluded.
    '''

    #Reduce the master file
    df = ReduceMasterFile(df, meld_cols, end_date, age_min = age_min, age_max = age_max, exc_cands = exc_cands, status1A_cands = status1A_cands, status1B_cands = status1B_cands)

    #For each MELD definition, generate the A and B files and save them
    for meld_col in meld_cols:
        GenerateAFile(df, meld_col, end_date, folder, extra_cols)
        GenerateBFile(df, meld_col, end_date, folder, extra_cols)
                

        
    

    

                    

