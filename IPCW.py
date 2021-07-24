'''This script defines the IPCW class. This class is specifically designed to use the A and B
    files to calculate a bias-corrected survival. See the class definition for details.

    Name: Nicholas Wood, PhD
    Institution: USNA
    Email: nwood@usna.edu
'''
import pandas as pd
from KaplanMeier import KaplanMeier as KM
from bisect import bisect
from matplotlib import pyplot as plt
import numpy as np


class IPCW:

    def __init__(self, dfA, dfB, t = 90, ident_col = 'px_id', time_col = 'Days', meld_col = 'MELD-Na', death_col = 'Death', tx_col = 'Transplant', fast = True):
        '''The IPCW Survival estimator.

            dfA: pandas DataFrame. The A file dataframe. Used to calculate the probability of remaining
                    untransplanted.

            dfB: pandas DataFrame. The B file dataframe. Used to calculate bias-corrected survival (bias-
                    correction is done using info from the A file).

            t: Non-Negative int. The time up to which survival is calculated.

            ident_col: String. Name of column in the df indicating the identifier for candidates.

            time_col: String. Name of the column in the df indicating the time to event.

            meld_col: String. Name of the column in the df indicating the definition of MELD.

            death_col: String. Name of the column in the df indicating death (1 = death, 0 = not death).

            tx_col: String. Name of the column in the df indicating transplant (1 = transplant, 0 = not transplant).

            fast: Bool. If true, only uses a maximum of 500 candidates at each MELD score. Useful for quick testing.

            See the method definitions for details of their purpose and use..
        '''


        #Store all the relevant info as attributes of this class
        self.dfA = dfA
        self.dfB = dfB
        self._t = t
        self._ident_col = ident_col
        self._time_col = time_col
        self._meld_col = meld_col
        self._death_col = death_col
        self._tx_col = tx_col

        #Create two attibutes that will store, for each meld score, the survival and variance
        self.meld2S = {}
        self.meld2varS = {}


        #For each MELD score create the data sets from the B files for later use.
        self.CreateDataSetsByMELD()

        #If fast is true, then for each MELD score we will randomly sample a number of
        #candidates equal to the number of candidates there are for the largest MELD
        if fast:
            #Only do a relatively small number of candidates to speed things up
            n = 500

            for meld in self.meld2dfB:
                print(meld, len(self.meld2dfB[meld]))
                if n < len(self.meld2dfB[meld]):
                    self.meld2dfB[meld] = self.meld2dfB[meld].sample(n)


    def CalculateSurvivalByMELD(self):

        #Use the A file to get the probability of being transplanted at each MELD
        self.GetProbabilityOfRemainingUntransplantedByMELD()


        #For each meld score calucalte survival
        for m in self.meld2dfB:
            print(f'{self._meld_col} {m} Survival Calculation in Progress!')
            self.CalculateSurvival(m, self._t, self.meld2dfB[m])


    def CreateDataSetsByMELD(self):
        '''For each MELD score we need to create a data set of all candidates
            who ever had that MELD score. In the even a candidate had a certain MELD score more than once,
            we choose the first.'''

        #Grab the B file
        df = self.dfB

        #Create a dictionary to store the results
        meld2df = {}

        #Group the dataframe by MELD
        for meld, mdf in df.groupby(self._meld_col):

            inds = []

            #Group the dataframe by px_id
            for px_id, pdf in mdf.groupby(self._ident_col):

                #Only use the first
                inds.append(pdf.index[0])

            
            meld2df[meld] = df.loc[df.index.isin(inds)]

        #Store the results
        self.meld2dfB = meld2df



    def GetProbabilityOfRemainingUntransplantedByMELD(self):
        '''For each MELD score, use the appropriate dfA file to get the probability of remaining untransplanted.'''

        #Create the kaplan meier object. Not that "failure" in this case is transplant.
        km = KM(self.dfA, t = self._t, ident_col = self._ident_col, time_col = self._time_col, meld_col = self._meld_col, failure_col = self._tx_col)

        #Within that km object calculate survival at each meld. 
        for meld in sorted(list(self.dfA[self._meld_col].unique())):
            km.CalculateSurvival(meld)

        #Store the km object as an attribute for future use.
        self.TxKM = km


    def CalculateSurvival(self, meld, t, dfB_meld):
        '''Calculate survival for all candidates at the given MELd score up to time t.'''

        #Read in the A File
        dfA = self.dfA

        #Determine all unique failure times. These are the times at which we need to calculate the weights
        FailureTimes = sorted(list(dfB_meld.loc[dfB_meld[self._death_col] == True][self._time_col].unique()))
        
        #only include failure times up to and including t
        if not t is None:
            FailureTimes = [ft for ft in FailureTimes if ft <= t]

        else:
            try:
                FailureTimes = [ft for ft in FailureTimes if ft <= self._t]
            except TypeError:
                pass

        vals = []

        #Get a dictionary of identifiers to A file dfs
        px2dfA = {px: df for px, df in dfA.groupby(self._ident_col)}

        def GetWeights(row, px2dfA, FailureTimes, TxKM):

            #Grab the index
            ind = row.name

            #Determine whether or not this person died
            death = row[self._death_col]

            #Find this candidates identifier
            px_id = row[self._ident_col]

            #find the rows in the dfA file corresponding to this candidate
            pdfA = px2dfA[px_id]

            #Only keep the rows with an index value which is greater than or equal to ind. This only works
            #because of how these files are indexed
            pdfA = pdfA.loc[pdfA.index >= ind]


            #Determine the max time we need
            maxT = pdfA[self._time_col].sum()
            
            
            for ft in FailureTimes:

                if ft > maxT:
                    break

                else:

                    #Create a variable to store the k_hat value
                    k_hat = 1

                    #We need to construct k_hat
                    time = 0
                    
                    for ind1, row1 in pdfA.iterrows():

                        #Determine the MELD score for this row
                        m = row1[self._meld_col]

                        #If the total time plus the time at this row is less than to the failure
                        #time, then then entire time this person was at this MELD should contribute to k_hat
                        if time + row1[self._time_col] < ft:

                            #Update the total time
                            time += row1[self._time_col]

                            #Update k_hat
                            k_hat *= self.TxKM.SurvivalAtTime(row1[self._time_col], m)


                        #If the total time plus the time at this row is greater than or equal to the failure time, then use
                        #the time at this MELD that gets us to the faiure time contributes to k_hat. Then we stop.
                        else:
                            k_hat *= self.TxKM.SurvivalAtTime(ft - time, m)

                            break

                    
                    if ft == row['Days'] and death == True:
                        dth_flag = 1


                    else:
                        dth_flag = 0



                    vals.append([px_id, ft, dth_flag, 1/k_hat])


                
            

        #Iterate over each candidate in the dfB_meld dataframe
        dfB_meld.apply(lambda row: GetWeights(row, px2dfA, FailureTimes, self.TxKM), axis = 1)

        W = pd.DataFrame(data = vals, columns = [self._ident_col, 'Days', self._death_col, 'W_hat'])

        #And now we can calculate survival! (I think)
        #For each unique failure time, calculate the weighted number of individuals who
        #died and are at risk
        DN = {ft: (tdf.loc[tdf[self._death_col] == True]['W_hat'].sum(), tdf['W_hat'].sum()) for ft, tdf in W.groupby(self._time_col)}

        #Calculate the survival and standard error
        S = {0:1}
        varS = {0:0}
        s_hat = 1

        lambd = []
        M = []

        D = []
        N = []

        #Calculate survival and variance of survival at each failure time
        for ft in FailureTimes:

            d, n = DN[ft]
            
            S[ft] = s_hat*(1 - d/n)

            D.append(d)
            N.append(n)

            lambd.append(sum(D)/sum(N))
            M.append(n**2/((W.loc[W[self._time_col] == ft]['W_hat']**2).sum()))
            

            try:
                varS[ft] = s_hat**2*sum(L/(m*(1 - L)) for m, L in zip(M, lambd))
            except ZeroDivisionError:
                varS[ft] = np.nan
                

            #update s_hat
            s_hat = S[ft]

        #Store the results as attributes of the class
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

        M = self.meld2S.keys()
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



    def SaveSurvivalData(self, filename, maxT = 90):
        '''Save all the survival data!'''

        MELDs = self.meld2S.keys()

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
            
                      



        
