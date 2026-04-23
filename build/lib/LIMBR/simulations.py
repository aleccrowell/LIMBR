import numpy as np
import pandas as pd
import pickle
import string
from sklearn.preprocessing import scale
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import random

class simulate:
    """
    Generates a properly formulated simulated circadian dataset for testing.


    This class takes user specifications on dataset size along with the number and frequency of batch effects and levels of background noise and outputs a simulated dataset along with a key showing which rows represent truly circadian data.


    Parameters
    ----------
    tpoints : int
        
    nrows : int

    nreps : int

    tpoint_space : int
        
    pcirc : float

    phase_prop : float

    phase_noise : float

    amp_noise : float

    n_batch_effects : int

    pbatch : float

    Probability of each batch effects appearance in a given peptide.

    effect_size : float

    Average size of batch effect

    p_miss : float

    probability of missing data

    lam_miss : float

    poisson lambda for # of missing columns in row with missing data



        


    Attributes
    ----------

    simdf : dataframe

    Simulated data without noise.

    simndf : dataframe

    Simulated data with noise.

    """

    def __init__(self,tpoints=24,nrows=1000, nreps=3, tpoint_space=2, pcirc=.5, phase_prop=.5, phase_noise=.25, amp_noise=.75, n_batch_effects=3, pbatch=.5, effect_size=2, p_miss=.2,lam_miss=5):
        """
        Simulates circadian data and saves as a properly formatted example .csv file.

        Takes a file from one of two data types protein ('p') which has two index columns or rna ('r') which has only one.  Opens a pickled file matching pooled controls to corresponding samples if data_type = 'p' and opens a picked file matching samples to blocks if designtype = 'b'.

        """

        self.tpoints = int(tpoints)
        self.nreps = int(nreps)
        self.nrows = int(nrows)
        self.tpoint_space = int(tpoint_space)
        self.pcirc = float(pcirc)
        self.phase_prop = float(phase_prop)
        self.phase_noise = float(phase_noise)
        self.amp_noise = float(amp_noise)
        self.n_batch_effects = int(n_batch_effects)
        self.pbatch = float(pbatch)
        self.effect_size = float(effect_size)
        self.p_miss = float(p_miss)

        #procedurally generate column names
        self.cols = []
        for i in range(self.tpoints):
            for j in range(self.nreps):
                self.cols.append("ZT"+("{0:0=2d}".format(self.tpoint_space*i+self.tpoint_space))+'_'+str(j+1))

        #randomly determine which rows are circadian
        self.circ = np.random.binomial(1, self.pcirc, self.nrows)
        #generate a base waveform
        base = np.arange(0,(4*np.pi),(4*np.pi/self.tpoints))
        #simulate data
        self.sim = []
        phases = []
        for i in self.circ:
            if i == 1:
                temp=[]
                p = np.random.binomial(1, self.phase_prop)
                phases.append(p)
                for j in range(self.nreps):
                    temp.append(np.sin(base+np.random.normal(0,self.phase_noise,1)+np.pi*p)+np.random.normal(0,self.amp_noise,self.tpoints))
                self.sim.append([item for sublist in zip(*temp) for item in sublist])
            else:
                phases.append('nan')
                self.sim.append(np.random.normal(0,self.amp_noise,(self.tpoints*self.nreps)))
        #add in batch effects
        batch_effects = []
        for i in range(self.n_batch_effects):
            batch_effects.append(np.random.normal(0.,self.effect_size,(self.tpoints*self.nreps)))
        self.simnoise = []
        for i in self.sim:
            temp = i
            bts = np.random.binomial(1, self.pbatch, self.n_batch_effects)
            for j in range(self.n_batch_effects):
                temp += bts[j]*batch_effects[j]
            self.simnoise.append(temp)
        #p_miss is prob a row is missing any data, generate list of length rows using binomial, then for entries in this list
        #generate lists of length columns with number of missing values drawn from poisson and append these lists to the mask
        #sample from poisson to get num_miss make list of num_miss 1s plus num_cols - num_miss 0s then shuffle with np.shuffle
        miss_rows = np.random.binomial(1, self.p_miss, self.nrows)
        m = []
        for i in miss_rows:
            if i == 1:
                num_miss = np.random.poisson(lam_miss) + 1
                mini_mask = np.asarray(num_miss*[1] + (self.tpoints*self.nreps-num_miss)*[0])
                np.random.shuffle(mini_mask)
            else:
                mini_mask = self.tpoints*self.nreps*[0]
            m.append(mini_mask)

        self.sim_miss = np.ma.masked_array(np.asarray(self.simnoise, dtype=object), mask=m).filled(np.nan)


    def generate_pool_map(self, out_name='pool_map'):
        """

        out_name : str

        output file stem

        """
        self.out_name = str(out_name)
        pool_map = {}
        for i in self.cols:
            pool_map[i] = 1
        pickle.dump(pool_map, open(out_name+'.p', "wb") )


    def write_output(self, out_name='simulated_data'):
        """

        out_name : str

        output file stem

        """

        self.out_name = str(out_name)
        
        self.simndf = pd.DataFrame(self.sim_miss,columns=self.cols).fillna('NULL')
        peps = [''.join([random.choice(string.ascii_uppercase) for j in range(12)]) for i in range(len(self.simndf))]
        prots = [''.join([random.choice(string.ascii_uppercase) for j in range(12)]) for i in range(len(self.simndf))]
        self.simndf.insert(0, 'Protein', prots)
        self.simndf.insert(0, 'Peptide', peps)
        self.simndf.set_index('Peptide',inplace=True)
        self.simndf.insert((len(self.cols)+1), 'pool_01', ['1']*len(self.simndf))
        self.simndf.to_csv(out_name+'_with_noise.txt',sep='\t')

        self.simdf = pd.DataFrame(np.asarray(self.sim),columns=self.cols)
        self.simdf.insert(0, 'Protein', prots)
        self.simdf.insert(0, 'Peptide', peps)
        self.simdf.set_index('Protein',inplace=True)
        self.simdf.groupby(level='Protein').mean()
        self.simdf.index.names = ['#']
        self.simdf.to_csv(out_name+'_baseline.txt',sep='\t')

        pd.DataFrame(self.circ,columns=['Circadian'],index=self.simndf['Protein']).to_csv(out_name+'_true_classes.txt',sep='\t')



class analyze:
    def __init__(self,filename_classes):
        self.true_classes = pd.read_csv(filename_classes,sep='\t')
        self.tags = {}
        self.merged = []
        self.i = 0

    def add_data(self,filename_ejtk,tag,include_missing=True):
        self.tags[tag] = self.i
        ejtk = pd.read_csv(filename_ejtk,sep='\t')
        if include_missing == True:
            self.merged.append(pd.merge(self.true_classes[['Protein','Circadian']], ejtk[['ID','GammaBH']], left_on='Protein', right_on='ID',how='left'))
        else:
            self.merged.append(pd.merge(self.true_classes[['Protein','Circadian']], ejtk[['ID','GammaBH']], left_on='Protein', right_on='ID',how='inner'))
        self.merged[self.i].set_index('Protein',inplace=True)
        self.merged[self.i].drop('ID',axis=1,inplace=True)
        self.merged[self.i].fillna(1.0,inplace=True)
        self.i += 1

    def generate_roc_curve(self):
        for j in self.tags.keys():
            fpr, tpr, thresholds = roc_curve(self.merged[self.tags[j]]['Circadian'].values, (1-self.merged[self.tags[j]]['GammaBH'].values), pos_label=1)
            plt.plot(fpr,tpr,label=j)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic Comparison')
        plt.legend(loc="lower right")
        plt.savefig('ROC.pdf')

    def calculate_auc(self):
        out = {}
        for j in self.tags.keys():
            fpr, tpr, thresholds = roc_curve(self.merged[self.tags[j]]['Circadian'].values, (1-self.merged[self.tags[j]]['GammaBH'].values), pos_label=1)
            out[j] = auc(fpr, tpr)
        self.roc_auc = out 

