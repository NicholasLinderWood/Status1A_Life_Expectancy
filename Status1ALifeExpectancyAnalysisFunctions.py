'''This script contains all the functions I used to generate the results for the
    "Life expectancy without a transplant for status 1A liver transplant candidates"
    paper".
    
    Name: Nicholas Wood, PhD
    Institution: USNA
    Email: nwood@usna.edu
'''

import pandas as pd
from IPCW import IPCW
from matplotlib import pyplot as plt
from matplotlib import gridspec, rcParams
import io
from PIL import Image




def Add1AInfo(AFile_path, BFile_path, start_date,
              stathist_liin_path = r"C:\Users\Nicholas\Documents\SRTR Data\Liver\stathist_liin.dta",
              statjust_li1_path = r"C:\Users\Nicholas\Documents\SRTR Data\Liver\statjust_li1.dta"):
    '''Add 1A info to the given A and B files. Then save the results in the same
        path but with the added 1A categories.

        AFile_Path: String. Path for the A file

        BFile_Path: String. Path for the B file

        start_date: String. Start date of the study

        stathist_liin_path: String. Path for stathist_liin.dta SRTR SAF

        statjust_li1_path: String. Path for statjust_li1.dta SRTR SAF


    Resulting A and B files with added 1A information are saved in the same location as the original A
    and B files, and they are additionally returned.
    '''

    #Read in all relevant files
    dfA = pd.read_csv(AFile_path)
    dfB = pd.read_csv(BFile_path)
    Status = pd.read_stata(stathist_liin_path)
    StatJust = pd.read_stata(statjust_li1_path)


    #First, find all candidates who ever achieved status 1A
    px_id = dfA.loc[dfA['MELD-Na'] == 41]['px_id'].unique()


    #For each of these candidates, find the first instance they had status 1A in the given time period
    Status = Status.loc[Status['px_id'].isin(px_id)]
    StatJust = StatJust.loc[StatJust['px_id'].isin(px_id)]
    px2date = {}

    for px_id, df in Status.groupby('px_id'):

        df = df.loc[(df['canhx_begin_dt'] >= start_date) & (df['canhx_stat_cd'].str.contains('Status 1A', na = False))]

        df = df.sort_values(by = 'canhx_begin_dt')

        px2date[px_id] = df.iloc[0]['canhx_begin_dt']



    #Next, for each candidate find the row in the StatJust file which corresponds to most recent form
    #submitted and accepted for status 1A
    inds = []
    for px_id, df in StatJust.groupby('px_id'):

        df = df.loc[df['canhx_chg_dt'] <= px2date[px_id]]

        df = df.sort_values(by = 'canhx_chg_dt', ascending = False)

        for ind, row in df.iterrows():
            if row['canhx_crit_not_met'] == False:
                inds.append(ind)
                break


    #Reduce the status justification dataframe to only those found status justification forms
    StatJust = StatJust.loc[StatJust.index.isin(inds)]


    #Create a new dataframe to store results
    Stat1AJust = pd.DataFrame(columns = ['Fulminant Liver Failure', 'Unknown', 'Primary Non Function or HAT', 'Wilsons'])

    #Loop over each row. Find out why the given candidate was given status 1A
    for ind, row in StatJust.iterrows():

        #Acute decompensted wilson's disease flag
        wilsons = row['canhx_wilsons']

        #primary non-function of whole or segmented liver transplant within seven days
        #of transplant, or HAT within seven days of transplant flag
        primary_non_function = row['canhx_prime_nonfunctn']

        #fulminant hepatic failure flag
        fulminant_hepatic_failure = row['canhx_fulminant_fail']

        #Anhepatic could not be determined, so if none of the others apply
        #I choose to not assume it is anhepatic. Instead I call this unknown.
        if sum((wilsons, primary_non_function, fulminant_hepatic_failure)) == 0:
            unknown = 1
        else:
            unknown = 0

        Stat1AJust.loc[row['px_id']] = [fulminant_hepatic_failure, unknown, primary_non_function, wilsons]




    #Create a mapping of each px_id to the status 1A dgn
    px2dgn = {}

    for px_id, row in Stat1AJust.iterrows():

        for col in Stat1AJust.columns:
            if row[col] == 1:
                px2dgn[px_id] = col


    for px_id in dfA.loc[dfA['MELD-Na'] == 41]['px_id'].unique():

        if px_id not in px2dgn:
            px2dgn[px_id] = 'Unknown'

    S = pd.Series(px2dgn)
    

    #Now add this info to the A and B dataframes
    dfA['1A Category'] = dfA['px_id'].map(S)
    dfB['1A Category'] = dfB['px_id'].map(S)


    #Save
    dfA.to_csv(f'{AFile_path[0:-4]} With 1A Categories.csv', index = False)
    dfB.to_csv(f'{BFile_path[0:-4]} With 1A Categories.csv', index = False)

    return dfA, dfB


def GetStatus1ACohortTable(dfA, end_date, cand_liin_path = r'C:/Users/Nicholas/Documents/SRTR Data/Liver/cand_liin.dta'):
    '''Creates a table of information regarding the status 1A candidates. This is
        essentially table 1 in the paper.

        dfA: pandas DataFrame - This is dfA after status 1A info is added. dfB could also be used.

        end_date: String. End date of the study.

        cand_liin_path: String. Indicates the location of the cand_liin.dta SRTR SAF.

        
    The resulting table is saved in the local directory.
    '''

    #Read in the Candidates file
    Candidates = pd.read_stata(cand_liin_path, index_col = 'px_id')

    #Define the removal codes that correspond to transplant
    tx_cds = ['           4: Deceased Donor tx, removed by transplanting center',
              '          18: Deceased Donor Emergency Tx',
              '          21: Patient died during TX procedure',
              '          19: Deceased Donor Multi-Organ Tx']

    #Define the removal codes that correspond to death
    dth_cds = ['           8: Died',
               '          13: Candidate condition deteriorated , too sick for tx',
               '           5: Medically Unsuitable']


    #Create a results dataframe to store all the info
    Results = pd.DataFrame(index = ['Overall', 'Fulminant Liver Failure', 'Primary Non Function or HAT', 'Wilsons', 'Unknown'],
                           columns = ['N', 'Transplanted (%)', 'Died (%)'])



    #Reduce the dataframe to only those candidates who ever achieved status 1A
    dfA = dfA.loc[dfA['MELD-Na'] == 41]

    #Loop over each status 1A category and get the statistics
    for cat, df in dfA.groupby('1A Category'):

        #Get the unique px_ids of all these candidates
        px_id = df['px_id'].unique()

        #Find N
        n = len(px_id)

        #Grab those N candidates
        df = Candidates.loc[px_id]

        #Determine how many were transplanted
        n_Tx = ((df['can_rem_cd'].isin(tx_cds)) & (df['can_rem_dt'] < end_date)).sum()

        #Add percent information
        p = round(100*(n_Tx/n), 1)
        n_Tx = f'{n_Tx} ({p})'

        #Determine how many died
        n_Dth = ((df['can_rem_cd'].isin(dth_cds)) & (df['can_rem_dt'] < end_date)).sum()

        #Add percent information
        p = round(100*(n_Dth/n), 1)
        n_Dth = f'{n_Dth} ({p})'

        #Add this to the results data frame
        Results.loc[cat] = [n, n_Tx, n_Dth]


    #Repeat for overall
    #Get the unique px_ids of all these candidates
    px_id = dfA['px_id'].unique()

    #Find N
    n = len(px_id)

    #Grab those N candidates
    df = Candidates.loc[px_id]

    #Determine how many were transplanted
    n_Tx = ((df['can_rem_cd'].isin(tx_cds)) & (df['can_rem_dt'] < end_date)).sum()

    #Add percent information
    p = round(100*(n_Tx/n), 1)
    n_Tx = f'{n_Tx} ({p})'

    #Determine how many died
    n_Dth = ((df['can_rem_cd'].isin(dth_cds)) & (df['can_rem_dt'] < end_date)).sum()

    #Add percent information
    p = round(100*(n_Dth/n), 1)
    n_Dth = f'{n_Dth} ({p})'

    #Add this to the results data frame
    Results.loc['Overall'] = [n, n_Tx, n_Dth]

    

    #Save the results
    Results.to_csv('Status 1A Cohort Table.csv')


def Bootstrap(dfA, dfB, LE, n = 200, T = 90):
    '''Use bootstrapping to get confidence intervals on the life expectancy estimations.

        dfA: pandas DataFrame - this is dfA after 1A info has been added

        dfB: pandas DataFrame - this is dfB after 1A info has been added

        LE: pandas DataFrame - this is a dataframe created by CalculateLifeExpectancy
                to this dataframe we add confidence intervals

        n: Positive Integer - this is the number of bootstrap iterations. For testing I
                recommend a low number for time's sake.

        T: Positive Integer - The time up to which you want to cap survival estimation.
                In some cases this will set an upper bound on what the 97.5 percentile
                of life expectancies will be.
        
    The resulting LE table is saved. This is table 2 in the paper.
    '''

    #Store all px_id values
    AllPX = dfA['px_id'].unique()

    #Now create the ipcw object
    ipcw = IPCW(dfA, dfB, t = T)

    #Calculate the probability of remaining untransplanted. Assume this remains constant
    #across bootstrap iterations
    ipcw.GetProbabilityOfRemainingUntransplantedByMELD()
    

    #Now, reduce the dfA and dfB files to only those candidates who were ever 1A
    S = pd.Series(dfA.loc[dfA['MELD-Na'] == 41]['px_id'].unique())
    dfA_1A = dfA.loc[dfA['px_id'].isin(S)].copy()
    dfB_1A = dfB.loc[dfB['px_id'].isin(S)].copy()
    dfA_not1A = dfA.loc[~dfA['px_id'].isin(S)].copy()
    dfB_not1A = dfB.loc[~dfB['px_id'].isin(S)].copy()

    #Create a mapping of px_id to rows in the A and B file
    px2dfa = {px:df for px, df in dfA_1A.groupby('px_id')}
    px2dfb = {px:df for px, df in dfB_1A.groupby('px_id')}

    #Create a dataframe to store the bootstrapping results
    BS = pd.DataFrame(index = LE.index, columns = range(n))

    #Go through the bootstrap iterations
    for i in range(n):

        #This process is somewhat slow. So I print to know how quickly it is moving along.
        print(f'Bootstrap Iteration: {i+1}')

        #Create lists that will store the rows to be concatenated for the A and B files
        dfas = []
        dfbs = []

        #Sample the px_ids with replacement
        px_ids = list(S.sample(len(S), replace = True))

        temp_px_id = 1

        for px_id in px_ids:

            #Grab the rows for this candidate
            dfa = px2dfa[px_id].copy()
            dfb = px2dfb[px_id].copy()

            #Need to make sure this new px_id is not already a px_id
            while True:
                if temp_px_id in AllPX:
                    temp_px_id += 1

                else:
                    break

            #Change the px_id
            dfa['px_id'] = temp_px_id
            dfb['px_id'] = temp_px_id

            #Increment temp_px_id
            temp_px_id += 1

            #Add those rows to the list of things being concatenated
            dfas.append(dfa)
            dfbs.append(dfb)


        #Now create the new dfa and dfb
        dfA = pd.concat(dfas + [dfA_not1A], ignore_index = True)
        dfB = pd.concat(dfbs + [dfB_not1A], ignore_index = True)


        #Update A and B files
        ipcw.dfA = dfA
        ipcw.dfB = dfB


        #Create data sets by MELD
        ipcw.CreateDataSetsByMELD()


        #Grab the B file data set for Status 1A candidates
        dfB = ipcw.meld2dfB[41]      

        #Calculate Survival for each 1A category
        for cat, df in dfB.groupby('1A Category'):

        
            ipcw.CalculateSurvival(41, T, df)

            LEFlag = False

            #Determine the life expectancy
            for t in range(1, T+1):
                
                if ipcw.SurvivalAtTime(t, 41) < 0.5:
                    BS.loc[cat, i] = t-1

                    LEFlag = True
                    break


            if not LEFlag:
                BS.loc[cat, i] = T


        ipcw.CalculateSurvival(41, T, dfB)

        LEFlag = False

        #Determine the life expectancy
        for t in range(1, T+1):
            if ipcw.SurvivalAtTime(t, 41) < 0.5:
                BS.loc['Overall', i] = t-1

                LEFlag = True
                break


        if not LEFlag:
            BS.loc['Overall', i] = T


    for cat in LE.index:
        LE.loc[cat, '95% Confidence Interval'] = f'({BS.loc[cat].quantile(0.025)}, {BS.loc[cat].quantile(0.975)})'

    #Save the life expectancy table.
    LE.to_csv('Status 1A Life Expectancy.csv')
     
       

def CalculateLifeExpectancy():
    '''Read in the survival results, calculate the life expectancy, and save that information to a dataframe.
        I assume that Calculate1ASurvival has already been run, which creates and saves the needed file to
        then calculate life expectancy.'''

    Results = pd.read_csv('1A Survival by Category.csv')

    #Create a life expectancy dataframe
    LE = pd.DataFrame(index = Results['1A Category'].unique(), columns = ['Life Expectancy (Days)', '95% Confidence Interval'])


    #Calculate life expectancy for each category
    for cat, df in Results.groupby('1A Category'):

        #Set the index to time
        df.set_index('t', inplace = True)

        #Find all survivals that are above 50%
        df = df.loc[df['Survival'] >= 0.50]

        #Take the last time (index)
        t = df.index[-1]

        LE.loc[cat, 'Life Expectancy (Days)'] = t


    return LE

        

def Calculate1ASurvival(dfA, dfB, T = 90):
    '''Calculates status 1A survival overall and for each 1A category. Saves the results as a csv.

        dfA: pandas DataFrame - this is dfA after 1A info has been added

        dfB: pandas DataFrame - this is dfB after 1A info has been added

        T: Positive Integer - The time up to which you want to cap survival estimation.
                In some cases this will set an upper bound on what the 97.5 percentile
                of life expectancies will be.

    '''

    #Create the IPCW object
    ipcw = IPCW(dfA, dfB, t = T, fast = False)


    #Calculate the probability of remaining untransplanted at each MELD/1A.
    ipcw.GetProbabilityOfRemainingUntransplantedByMELD()


    #Grab the B file data set for Status 1A candidates
    dfB = ipcw.meld2dfB[41]


    #Create a dataframe to store all results
    Results = pd.DataFrame(columns = ['1A Category', 't', 'Survival', 'Standard Error'])

    ind = 0


    #Calculate Survival for each 1A category
    for cat, df in dfB.groupby('1A Category'):
        
        ipcw.CalculateSurvival(41, T, df)

        for t in range(1, T+1):
            Results.loc[ind] = [cat, t, ipcw.SurvivalAtTime(t, 41), ipcw.StandardErrorAtTime(t, 41)]
            ind += 1


    #Calculate Survival for 1A overall
    ipcw.CalculateSurvival(41, T, dfB)

    for t in range(1, T+1):
        Results.loc[ind] = ['Overall', t, ipcw.SurvivalAtTime(t, 41), ipcw.StandardErrorAtTime(t, 41)]
        ind += 1


    #Save the results
    Results.to_csv('1A Survival by Category.csv', index = False)




def Generate1ASurvivalFigures():
    '''Generates the 1A survival over time for 1A overall (Figure 1) and the categories of
        fulminant liver failure, primary non-function or HAT, wilson's and unknown (Figure 2).
        I assume that Calculate1ASurvival has already been run, which creates and saves the needed file to
        then plot survival.'''

    rcParams.update({'figure.autolayout': True})
    
    #Read in the results
    Results = pd.read_csv('1A Survival by Category.csv')



    #First generate a figure for the four subcategories
    fig = plt.figure(figsize = (16, 9), dpi = 150)
    gs = gridspec.GridSpec(2, 2, width_ratios = [0.5, 0.5])

    w = 4

    for (i, cat), letter in zip(enumerate(['Fulminant Liver Failure', 'Primary Non Function or HAT', 'Wilsons', 'Unknown'], start = 1), 'abcd'):

        df = Results.loc[Results['1A Category'] == cat]

        ax = fig.add_subplot(gs[i-1])
        

        if i in [1, 2]:
            ax.axes.xaxis.set_ticks([])

        if i in [2, 4]:
            ax.axes.yaxis.set_ticks([])

        if i in [3, 4]:
            plt.xticks(fontsize = 16)
            ax.set_xlabel('Days', fontsize = 18)

        if i in [1, 3]:
            plt.yticks(fontsize = 16)
            ax.set_ylabel('Without-Transplant Survival', fontsize = 18)


        df = df.set_index('t')

        df.loc[0] = [cat, 1, 0]

        df = df.sort_index()

        x = df.index
        y = df['Survival']

        y_upper = df['Survival'] + 1.96*df['Standard Error']
        y_lower = df['Survival'] - 1.96*df['Standard Error']

        ax.plot(x, y, color = 'k', linewidth = w, zorder = 1)

        ax.plot(x, y_upper, color = 'k', linestyle = '--', zorder = 1)
        ax.plot(x, y_lower, color = 'k', linestyle = '--', zorder = 1)

        xmin, xmax = ax.get_xlim()

        ymin, ymax = -0.05, 1.05

        ax.set_xlim(xmin, xmax)

        ax.set_ylim(ymin, ymax)

        ax.set_aspect((xmax - xmin)/(ymax - ymin))

        ax.set_title(f'({letter})', fontsize = 18)


    plt.tight_layout()
            
    png1 = io.BytesIO()

    fig.savefig(png1, format='png')
            
    png2 = Image.open(png1)

    png2.save(f'Status 1A Survival By Category.tiff')
    png1.close()



    
    #Next generate a figure for status 1A overall
    fig = plt.figure(figsize = (16, 9), dpi = 150)
    gs = gridspec.GridSpec(1, 1, width_ratios = [0.5])

    w = 4

    df = Results.loc[Results['1A Category'] == 'Overall']

    ax = fig.add_subplot(gs[0])
        
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    
    ax.set_xlabel('Days', fontsize = 18)
    ax.set_ylabel('Without-Transplant Survival', fontsize = 18)


    df = df.set_index('t')

    df.loc[0] = [cat, 1, 0]

    df = df.sort_index()

    x = df.index
    y = df['Survival']

    y_upper = df['Survival'] + 1.96*df['Standard Error']
    y_lower = df['Survival'] - 1.96*df['Standard Error']

    ax.plot(x, y, color = 'k', linewidth = w, zorder = 1)

    ax.plot(x, y_upper, color = 'k', linestyle = '--', zorder = 1)
    ax.plot(x, y_lower, color = 'k', linestyle = '--', zorder = 1)

    xmin, xmax = ax.get_xlim()

    ymin, ymax = -0.05, 1.05

    ax.set_xlim(xmin, xmax)

    ax.set_ylim(ymin, ymax)

    ax.set_aspect((xmax - xmin)/(ymax - ymin))


    plt.tight_layout()
            
    png1 = io.BytesIO()

    fig.savefig(png1, format='png')
            
    png2 = Image.open(png1)

    png2.save(f'Status 1A Survival Overall.tiff')
    png1.close()





    
