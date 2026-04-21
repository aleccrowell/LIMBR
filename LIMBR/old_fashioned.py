import numpy as np
import pandas as pd
from sklearn import preprocessing

class old_fashioned:
    """
    Performs a standard normalization procedure without SVD as a baseline.


    This class performs simple quantile normalization and row scaling along with pool normalization for proteomics experiments using the same methods and interface employed in the sva class.  This provides a baseline comparison point for data processed with LIMBR.


    Parameters
    ----------
    filename : str
        Path to the input dataset.
    data_type : str
        Type of dataset, one of 'p' or 'r'.  'p' indicates proteomic with two index columns specifying peptide and protein.  'r' indicates RNAseq with one index column indicating gene.
    pool : str
        Path to parquet file containing pooled control design for experiment in the case of data_type = 'p'.  Must have a 'pool_number' column indexed by sample column headers.


    Attributes
    ----------
    raw_data : dataframe
        This is where the input data is stored.
    data_type : str
        This is where the data type ('p' or 'r') is stored.
    norm_map : dict
        This is where the assignment of pooled controls to samples are stored if data_type = 'p'.

    """

    def __init__(self, filename,data_type,pool=None):
        """
        Imports data and initializes an old_fashioned object.


        Takes a file from one of two data types protein ('p') which has two index columns or rna ('r') which has only one.  Reads a parquet file matching pooled controls to corresponding samples if data_type = 'p'.

        """

        np.random.seed(4574)
        self.data_type = str(data_type)
        if self.data_type == 'p':
            self.raw_data = pd.read_csv(filename,sep='\t').set_index(['Peptide','Protein'])
        if self.data_type == 'r':
            self.raw_data = pd.read_csv(filename,sep='\t').set_index('#')
        if pool != None:
            self.norm_map = pd.read_parquet(pool)['pool_number'].to_dict()
        self.notdone = True

    def pool_normalize(self):
        """
        Preprocessing normalization.


        Performs pool normalization on an sva object using the raw_data and norm_map if pooled controls were used. Quantile normalization of each column and scaling of each row are then performed.


        Attributes
        ----------
        scaler : sklearn.preprocessing.StandardScaler()
            A fitted scaler from the sklearn preprocessing module.
        data_pnorm : dataframe
            Pool normalized data.

        """

        def pool_norm(df,dmap):
            """
            Pool normalizes samples in a proteomics experiment.


            Peptide abundances of each sample are divided by corresponding pooled control abundances.


            Parameters
            ----------
            df : dataframe
                The dataframe to be pool normalized.
            dmap : dict
                The dictionary connecting each sample to its corresponding pooled control.


            Returns
            -------
            newdf : dataframe
                Dataframe with samples pool normalized and pooled control columns dropped.

            """

            newdf = pd.DataFrame(index=df.index)
            for column in df.columns.values:
                if 'pool' not in column:
                    newdf[column] = df[column].div(df['pool_'+'%02d' % dmap[column]],axis='index')
            nonpool = [i for i in newdf.columns if 'pool' not in i]
            newdf = newdf[nonpool]
            return newdf

        def qnorm(df):
            """
            Quantile normalizes data by columns.


            A reference distribution is generated as the mean across rows of the dataset with all columns sorted by abundance.  Each column is then quantile normalized to this target distribution.


            Parameters
            ----------
            df : dataframe
                The dataframe to be quantile normalized


            Returns
            -------
            newdf : dataframe
                The quantile normalized dataframe.

            """

            ref = pd.concat([df[col].sort_values().reset_index(drop=True) for col in df], axis=1, ignore_index=True).mean(axis=1).values
            for i in range(0,len(df.columns)):
                df = df.sort_values(df.columns[i])
                df[df.columns[i]] = ref
            return df.sort_index()
        if self.data_type == 'r':
            self.data = qnorm(self.raw_data)
            self.scaler = preprocessing.StandardScaler().fit(self.data.values.T)
            self.data = pd.DataFrame(self.scaler.transform(self.data.values.T).T,columns=self.data.columns,index=self.data.index)
        else:
            self.data_pnorm = pool_norm(self.raw_data,self.norm_map)
            self.data_pnorm = self.data_pnorm.replace([np.inf, -np.inf], np.nan)
            self.data_pnorm = self.data_pnorm.dropna()
            self.data_pnorm = self.data_pnorm.sort_index(axis=1)
            self.data_pnorm = qnorm(self.data_pnorm)
            self.scaler = preprocessing.StandardScaler().fit(self.data_pnorm.values.T)
            self.data = pd.DataFrame(self.scaler.transform(self.data_pnorm.values.T).T,columns=self.data_pnorm.columns,index=self.data_pnorm.index)

    def normalize(self,outname):
        """
        Groups peptides by protein and outputs final processed dataset.


        These final results are then written to an output file.


        Parameters
        ----------
        outname : str
            Path to desired output file.

        """

        #self.old_norm = self.scaler.inverse_transform(self.data.values.T).T
        #self.old_norm = pd.DataFrame(self.old_norm,index=self.data.index,columns=self.data.columns)
        if self.data_type == 'p':
            self.data = self.data.groupby(level='Protein').mean()
        self.data.index.names = ['#']
        self.data.to_csv(outname,sep='\t')
