'''This script defines the KaplanMeier class. This class is specifically designed to be used by the
    IPCW class for the purpose of calculating bias-corrected survival. See the class definition for details.

    Name: Nicholas Wood, PhD
    Institution: USNA
    Email: nwood@usna.edu
'''


import pandas as pd
from matplotlib import pyplot as plt
from bisect import bisect
import numpy as np


class KaplanMeier:

    def __init__(self, df, t = 90, ident_col = 'px_id', time_col = 'Days', meld_col = 'MELD-Na', failure_col = 'Death'):
        '''The KaplanMeier Survival estimator.

            df: pandas DataFrame. Contains survival data.

            t: Non-Negative int. The time up to which survival is calculated.

            ident_col: String. Name of column in the df indicating the identifier for candidates.

            time_col: String. Name of the column in the df indicating the time to event.

            meld_col: String. Name of the column in the df indicating the definition of MELD.

            failure_col: String. Name of the column in the df indicating failure (1 = failure, 0 = not failure).

            See the method definitions for details of their purpose and use. In general, if want wanted to calculate
            survival by MELD (without bias correction) the df would be the generated B file.
        '''

        
        #store all relevant attributes
        self.df = df
        self._t = t
        self._ident_col = ident_col
        self._meld_col = meld_col
        self._time_col = time_col
        self._failure_col = failure_col



        #Create two attibutes that will store, for each meld score, the survival and variance
        self.meld2S = {}
        self.meld2varS = {}


    def CalculateSurvivalByMELD(self):
        '''Estimates survival for each MELD score.'''

        #For each meld score calculate survival
        for meld in list(sorted(self.df[self._meld_col].unique())):
            self.CalculateSurvival(meld)


    def CreateDataSetForMELD(self, m):
        '''Creates a survival data set for the given MELD score. If m == 'All' we ignore MELD.'''

        #Create an empty list to store the indices of the data frane we will use.
        inds = []

        #If m == 'All' we use the entire df. 
        if m == 'All':
            tdf = self.df

        #Otherwise we reduce the data set to only those candidates who ever had the given MELD.
        else:
            tdf = self.df.loc[self.df[self._meld_col] == m]

        #Group by candidates in tdf and grab the first row for each candidate
        for px_id, pdf in tdf.groupby(self._ident_col):

            inds.append(pdf.index[0])

        #Return the resultant data set
        return tdf.loc[tdf.index.isin(inds)]


    def CalculateSurvival(self, meld = 'All'):
        '''Calculates the KaplanMeier based survival and stores the result as a dictionary attribute.'''

        #Grab the time up to which survival is calculated
        t = self._t

        #Grab the dataframe
        df = self.df

        #Reduce the overarching dataset to determine survival at the given meld
        df = self.CreateDataSetForMELD(meld)

        #Determine the unique times at which failure occurs. These are the only times at which we
        #need to calculate survival. Convert to a sorted list
        FailureTimes = sorted(list(df.loc[df[self._failure_col] == True][self._time_col].unique()))

        #only include failure times up to and including t
        if not t is None:
            FailureTimes = [ft for ft in FailureTimes if ft <= t]

        else:
            try:
                FailureTimes = [ft for ft in FailureTimes if ft <= self._t]
            except TypeError:
                pass


        #Now for each unique failure time, we need to calculate the number of individuals who die
        #at that time and the number of individuals at risk (just before) that time. This will be stored
        #as a list of tuples
        DN = [(len(df.loc[(df[self._failure_col] == True) & (df[self._time_col] == t)]), len(df.loc[df[self._time_col] >= t])) for t in FailureTimes]

        #Calculate the survival and standard error
        S = {0:1}
        varS = {0:0}
        s_hat = 1

        #Calculate survival and variance of survival at each failure time
        for (i, (d, n)), t in zip(enumerate(DN, start = 1), FailureTimes):
            S[t] = s_hat*(1 - d/n)

            try:
                varS[t] = s_hat**2*sum(d/(n*(n-d)) for d, n in DN[0:i])
            except ZeroDivisionError:
                varS[t] = np.nan

            #update s_hat
            s_hat = S[t]


        #Store the results as attributes in the class
        self.meld2S[meld] = S
        self.meld2varS[meld] = varS


    
    def SurvivalAtTime(self, t, m):
        '''Determine the estimated survival at the given time for the given meld. If the time
            is not in the dictionary for S, the best estimate is the first time in S
            that is prior to t.'''

        times = sorted(list(self.meld2S[m].keys()))

        #If this is a failure time, return the survival estimate at that survival time
        if t in times:

            return self.meld2S[m][t]

        #If it is not a faiilure time, determine the first failure time before this time
        #and return survival at that time.
        else:

            ind = bisect(times, t) - 1
            t = times[ind]

            return self.meld2S[m][t]



    def StandardErrorAtTime(self, t, m):
        '''Determine the estimated standard error of survival at the given time for the given meld. If the time
            is not in the dictionary for varS, the best estimate is the first time in varS
            that is prior to t.'''

        times = sorted(list(self.meld2S[m].keys()))

        #If this is a failure time, return the standard error of survival estimate at that survival time
        if t in times:

            return self.meld2varS[m][t]**(1/2)

        #If it is not a faiilure time, determine the first failure time before this time
        #and return standard error of survival at that time.
        else:

            ind = bisect(times, t) - 1
            t = times[ind]

            return self.meld2varS[m][t]**(1/2)



    def Plot(self, m, show = True):
        '''Generates a plot for the survival at the given meld. returns the fig and ax objects.
            If show, then display the plot as well.'''

        fig = plt.figure(figsize = (16, 9), dpi = 72)
        ax = fig.add_subplot(111)

        ax.set_ylim([-0.05, 1.05])
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)

        times = list(self.meld2S[m].keys())
        survivals = list(self.meld2S[m].values())
        varsurvivals = list(self.meld2varS[m].values())

        T = [times[0]]
        S = [survivals[0]]

        for t, s in zip(times, survivals):

            T += [t, t]
            S.append(S[-1])
            S.append(s)


        #Plot the survival curve
        w = 4
        ax.plot(T, S, 'k-', linewidth = 1.5*w, zorder = 1)
        ax.plot(T, S, 'r-', linewidth = w, zorder = 2)

        #Plot the 95% confidence interval
        S_Upper = [survivals[0] + 1.96*varsurvivals[0]]
        S_Lower = [survivals[0] - 1.96*varsurvivals[0]]

        for s, var_s in zip(survivals, varsurvivals):
            S_Upper.append(S_Upper[-1])
            S_Lower.append(S_Lower[-1])

            S_Upper.append(s + 1.96*var_s**(1/2))
            S_Lower.append(s - 1.96*var_s**(1/2))

        ax.plot(T, S_Upper, 'r--', zorder = 3)
        ax.plot(T, S_Lower, 'r--', zorder = 3)


        if show:
            plt.show()


        return fig, ax


    def PlotSurvivalByMELD(self, t):
        '''Plots survival at each meld score at the given time.'''

        fig = plt.figure(figsize = (16, 9), dpi = 72)
        ax = fig.add_subplot(111)

        ax.set_ylim([-0.05, 1.05])
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)

        M = [m for m in self.meld2S.keys() if type(m) != str]
        S = [self.SurvivalAtTime(t, m) for m in M]
        SE = [self.StandardErrorAtTime(t, m) for m in M]

        S_Upper = [s + 1.96*se for s, se in zip(S, SE)]
        S_Lower = [s - 1.96*se for s, se in zip(S, SE)]

        #Plot the survival curve
        w = 4
        ax.plot(M, S, 'k-', linewidth = 1.5*w, zorder = 1)
        ax.plot(M, S, 'r-', linewidth = w, zorder = 2)
    

        #Add 95% CI
        ax.plot(M, S_Upper, 'r--', zorder = 3)
        ax.plot(M, S_Lower, 'r--', zorder = 3)

        plt.show()

        return fig, ax


    def SaveSurvivalData(self, filename):
        '''Save all the survival data!'''

        MELDs = [m for m in self.meld2S.keys() if type(m) != str]


        #Determine the max time
        maxT = self._t


        df = pd.DataFrame(index = range(1, maxT + 1))

        for m in MELDs:
            S = []
            seS = []

            for t in range(1, maxT + 1):

                S.append(self.SurvivalAtTime(t, m))
                seS.append(self.StandardErrorAtTime(t, m))

            df[f'{self._meld_col} {int(m)} Survival'] = S
            df[f'{self._meld_col} {int(m)} Standard Error'] = seS



        df.to_csv(filename)


    
    
    
