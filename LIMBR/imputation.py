import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm

class imputable:
    """Imputes missing data with K Nearest Neighbors based on user specified parameters.


    This class generates an object from raw data and allows for preprocessing, imputation and output of that data.


    Parameters
    ----------
    filename : str
        This is the path to the file containing the raw data.
    missingness : float
        This is the maximum allowable percentage of missingness expressed as a decimal.  For example a value of 0.25 would mean that all rows for which more than one in four values are missing will be rejected.


    Attributes
    ----------
    data : dataframe
        This is where the raw data read from filename is stored.
    miss : float
        This is where the maximum allowable percentage of missingness expressed as a decimal is stored.
    pats : dict
        This is where the set of missingness patterns found in the dataset will be stored later, initialized here.

    """

    def __init__(self, filename, missingness, neighbors=10):
        """
        Constructor, takes input data and missingness threshold and initializes imputable object.

        This is the initialization function for imputation.  It reads the input file of raw data and sets the user specified value for the missingness threshold and number of nearest neighbors.

        """
        self.data = pd.read_csv(filename,sep='\t')
        self.miss = float(missingness)
        self.pats = {}
        self.notdone = True
        self.NN = neighbors

    def deduplicate(self):
        """
        Removes duplicate peptides.


        Groups rows by peptide, if a peptide appears in more than one row it is removed.

        """
        if (self.data[self.data.columns.values[1]][0][-2] == "T") & (self.data[self.data.columns.values[1]][0][-1].isdigit()):
            self.data[self.data.columns.values[1]] = self.data[self.data.columns.values[1]].apply(lambda x: x.split('T')[0])
        
        self.data = self.data.groupby(['Peptide','Protein']).mean()
        todrop = []
        for name, group in tqdm(self.data.groupby(level='Peptide')):
            if len(group) > 1:
                todrop.append(name)
        self.data = self.data.drop(todrop)

    def drop_missing(self):
        """Removes rows which are missing more data than the user specified missingness threshold."""
        self.miss = np.rint(len(self.data.columns)*self.miss)
        self.data = self.data[self.data.isnull().sum(axis=1)<=self.miss]

    def impute(self,outname):
        """
        Imputes missing data with KNN and outputs the results to the specified file.


        First all of the missingness patterns present in the dataset are identified.  Then those patterns are iterated over and for each pattern, missing values are imputed.  Finally the dataset is reformed with the imputed values and output.


        Parameters
        ----------
        outname : str
            Path to output file.

        """

        def match_pat(l,i):
            """
            finds all missingness patterns present in the dataset


            For each row, if that row has a new missingness pattern, that pattern is added to the list, then whether the missingness pattern is new or not, the index of that row is assigned to the appropriate missingness pattern.


            Parameters
            ----------
            l : list
                A row of data
            i : int
                the index of that row in the original dataset

            """

            l = "".join(np.isnan(l).astype(int).astype(str))
            if l not in self.pats.keys():
                self.pats[l] = [i]
            else:
                self.pats[l].append(i)

        def get_patterns(arr):
            """Calls match_pat on all rows of data"""
            for ind, val in enumerate(arr):
                match_pat(val,ind)

        def sub_imputer(inds,pattern,origarr,comparr):
            """
            single imputation process for a missingness pattern.


            Drops columns missing in a given missingness pattern. Then finds nearest neighbors.  Iterates over rows matching missingness pattern, getting indexes of nearest neighbors, averaging nearest neighbrs and replacing 
missing values with corresponding averages.


            Parameters
            ----------
            inds : list
                indexes of rows sharing the missingness pattern.
            pattern : str
                Binary representation of missingness pattern.
            origarr : arr
                original array of data with missing values
            comparr : arr
                Complete array of only rows with no missing values (complete cases).


            Returns
            -------
            outa : arr
                Imputed array for missingness pattern.

            """

            #drop missing columns given missingness pattern
            newarr = comparr[:,~np.array(list(pattern)).astype(int).astype(bool)]
            #fit nearest neighbors
            nbrs = NearestNeighbors(n_neighbors=self.NN).fit(newarr)
            outa = []
            #iterate over rows matching missingness pattern
            for rowind, row in enumerate(origarr[inds]):
                outl = []
                #get indexes of given rows nearest neighbors
                indexes = nbrs.kneighbors([origarr[inds[rowind],~np.array(list(pattern)).astype(int).astype(bool)]],return_distance=False)
                #get array of nearest neighbors
                means = np.mean(comparr[indexes[0][1:]], axis=0)
                #iterate over entries in each row
                for ind, v in enumerate(row):
                    if not np.isnan(v):
                        outl.append(v)
                    else:
                        outl.append(means[ind])
                outa.append(outl)
            return outa

        def imputer(origarr, comparr):
            """
            Calls sub_imputer on each missingness pattern and outputs the results to a dict.


            Parameters
            ----------
            origarr : arr
                Original array with missing values.
            comparr : arr
                Complete array of only rows with no missing values (complete cases).


            Returns
            -------
            outdict : dict
                Dict of imputed data with index in the original dataset as the key and imputed data as the value.

            """

            outdict = {}
            for k in tqdm(self.pats.keys()):
                temparr = sub_imputer(self.pats[k],k, origarr,comparr)
                for ind, v in enumerate(temparr):
                    outdict[self.pats[k][ind]] = v
            return outdict

        datavals = self.data.values
        #generate array of complete cases
        comparr = datavals[~np.isnan(datavals).any(axis=1)]

        #find missingness patterns
        get_patterns(datavals)

        #impute
        out = imputer(datavals, comparr)

        #reform dataframe with imputed values from outdict
        meld = pd.DataFrame.from_dict(out,orient='index')
        meld.index = meld.index.astype(float)
        meld.sort_index(inplace=True)
        meld.set_index([self.data.index.get_level_values(0),self.data.index.get_level_values(1)], inplace=True)
        meld.columns = self.data.columns
        meld.to_csv(outname,sep='\t')

    def impute_data(self,out_file):
        self.deduplicate()
        self.drop_missing()
        self.impute(out_file)

