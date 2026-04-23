import itertools

import numpy as np
import pandas as pd
from multiprocess import Pool
from sklearn import preprocessing
from sklearn.decomposition import PCA
from scipy.stats import linregress, f_oneway
from statsmodels.nonparametric.smoothers_lowess import lowess
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Module-level helpers for parallel permutation testing
# ---------------------------------------------------------------------------

# Per-worker globals populated once by _init_perm_worker (Pool initializer).
# This avoids re-serialising the large residual matrix for every task.
_perm_res = _perm_tks = _perm_designtype = _perm_tpoints = _perm_block_design = None


def _init_perm_worker(res, tks, designtype, tpoints, block_design):
    """Pool initializer: stores shared read-only data in each worker once."""
    global _perm_res, _perm_tks, _perm_designtype, _perm_tpoints, _perm_block_design
    _perm_res = res
    _perm_tks = tks
    _perm_designtype = designtype
    _perm_tpoints = tpoints
    _perm_block_design = block_design


def _perm_worker(rseed):
    """Single permutation iteration executed inside a pool worker."""
    rng = np.random.default_rng(rseed * 100)
    rstar = rng.permuted(_perm_res.copy(), axis=1)

    if _perm_designtype in ('c', 't'):
        resstar = np.array([row - lowess(row, _perm_tpoints, it=1)[:, 1] for row in rstar])
    else:
        blocks = np.asarray(_perm_block_design)
        ma = np.empty_like(rstar)
        for v in np.unique(blocks):
            idx = blocks == v
            ma[:, idx] = rstar[:, idx].mean(axis=1, keepdims=True)
        resstar = rstar - ma

    pca = PCA(svd_solver='randomized', random_state=rseed * 100)
    pca.fit(resstar)
    tkstar = pca.explained_variance_ratio_

    n = min(len(tkstar), len(_perm_tks))
    out = np.zeros(len(_perm_tks))
    out[:n] = tkstar[:n] > _perm_tks[:n]
    return out

class sva:
    """
    Performs sva based identification and removal of batch effects.


    This class takes a dataset without missing values (raw RNAseq or Proteomics from class imputable).  It performs pool normalization for proteomics datasets and quantile normalization and scaling for all datasets.  The data is then subsetted based on correlation to a primary variable of interest determined by the experimental design.  This subsetted data is used to calculate a residual matrix from which initial estimates of batch effects are produced by SVD.  Significance of these batch effects are estimated by permutation of the residual matrix.  Batch effects deemed significant are regressed against the original dataset and final batch effects are calculated from those rows most correlated with the effect.


    Parameters
    ----------
    filename : str
        Path to the input dataset.
    design : str
        Experimental design, one of: 'c', 't', 'b'.  'c' indicates a circadian timecourse. 't' indicates a general timecourse (lowess fit used).  'b' indicates a block design.
    data_type : str
        Type of dataset, one of 'p' or 'r'.  'p' indicates proteomic with two index columns specifying peptide and protein.  'r' indicates RNAseq with one index column indicating gene.
    blocks : str
        Path to parquet file containing block design in the case of design = 'b'.  Must have a 'block' column whose values (integers) indicate the block assignment for each sample column, in order.
    pool : str
        Path to parquet file containing pooled control design for experiment in the case of data_type = 'p'.  Must have a 'pool_number' column indexed by sample column headers.


    Attributes
    ----------
    raw_data : dataframe
        This is where the input data is stored.
    data_type : str
        This is where the data type ('p' or 'r') is stored.
    designtype : str
        This is where the design type of the experiment ('c','l' or 'b') is stored.
    block_design : list
        This is where the block assignments for each sample are stored if designtype = 'b'.
    norm_map : dict
        This is where the assignment of pooled controls to samples are stored if data_type = 'p'.

    """

    def __init__(self, filename,design,data_type,blocks=None,pool=None):
        """
        Imports data and initializes an sva object.


        Takes a file from one of two data types protein ('p') which has two index columns or rna ('r') which has only one.  Reads a parquet file matching pooled controls to corresponding samples if data_type = 'p' and a parquet file matching samples to blocks if designtype = 'b'.

        """

        np.random.seed(4574)
        self.data_type = str(data_type)
        if self.data_type == 'p':
            self.raw_data = pd.read_csv(filename,sep='\t').set_index(['Peptide','Protein'])
        if self.data_type == 'r':
            self.raw_data = pd.read_csv(filename,sep='\t').set_index('#')
        self.designtype = str(design)
        if self.designtype == 'b':
            self.block_design = pd.read_parquet(blocks)['block'].tolist()
        if pool != None:
            self.norm_map = pd.read_parquet(pool)['pool_number'].to_dict()
        elif pool == None:
            self.norm_map = None
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
            newdf = df.sort_index()
            return newdf

        if (self.data_type == 'r') or (self.norm_map == None):
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


    def get_tpoints(self):
        """
        Extracts timepoints from header of data.


        Splits strings in header based on required syntax and generates timepoints.


        Attributes
        ----------
        tpoints : array
            array of timepoints at which samples were collected.

        """

        tpoints = [i.replace('ZT', '').replace('CT','') for i in self.data.columns.values]
        tpoints = [int(i.split('_')[0]) for i in tpoints]
        #deprecated splitting for alternative header syntax
        #tpoints = [int(i.split('.')[0]) for i in tpoints]
        self.tpoints = np.asarray(tpoints)

    def prim_cor(self):
        """
        calculates correlation for each row against the primary variable of interest.


        The primary variable of interest is defined based on the designtype.  For circadian ('c') this is the difference between the autocorrelation at the expected period and one half that period.  For general timecourse ('l') this is the goodness of fit of a lowess model.  For a block design ('b') this is the ANOVA f statistic.


        Attributes
        ----------
        cors : array
            Array of correlations with primary variable of interest.  The calculation of this variable is determined by the designtype as specified above and the procedure described in the corresponding function.

        """
        def circ_cor():
            """
            Calculates rough estimate of circadianness based on autocorrelation differences.


            Given a period (fixed here at 12 samples) the autocorrelation is calculated at 12 and 6 samples shift and the difference indicates the degree of circadian signal.

            """
            def autocorr_vec(m, shift):
                """Row-wise autocorrelation of a 2-D array at a given shift."""
                shifted = np.roll(m, shift, axis=1)
                return np.einsum('ij,ij->i', m, shifted) / np.einsum('ij,ij->i', m, m)

            per = 12
            tps = np.asarray(self.tpoints)
            data = self.data.values
            # Build (nrows, n_unique_tpoints) matrix of per-timepoint means
            ave = np.column_stack([
                data[:, tps == tp].mean(axis=1) * 1e6
                for tp in np.unique(tps)
            ])
            self.cors = autocorr_vec(ave, per) - autocorr_vec(ave, per // 2)

        def l_cor():
            """
            Calculates how well a row of data fits expectations for a timecourse.


            In a timecourse, rows least affected by batch effects should be those with the best fit from a lowess model.  Goodness of fit in this case is defined as the sum of squared errors for the lowess fit.

            """
            cors = [-np.sum((row - lowess(row, self.tpoints, it=1)[:, 1]) ** 2)
                    for row in tqdm(self.data.values)]
            self.cors = np.asarray(cors)

        def block_cor():
            """
            Calculates how well a row of data fits expectations for a block study design.


            In a blocked study, rows least affected by batch effects should be those for which within group variation is much smaller than between group variation.  In this case the ANOVA f statistic is used as a relative measure of within and between group variability.

            """
            blocks = np.asarray(self.block_design)
            # Pre-compute block indices once rather than rebuilding per row
            block_indices = [np.where(blocks == b)[0] for b in np.unique(blocks)]
            cors = [f_oneway(*[row[idx] for idx in block_indices])[0]
                    for row in tqdm(self.data.values)]
            self.cors = np.asarray(cors)

        if self.designtype == 'c':
            circ_cor()
        elif self.designtype == 'b':
            block_cor()
        elif self.designtype == 't':
            l_cor()


    def reduce(self,perc_red=25):
        """
        Reduces the data based on the correlation to primary variable of interest calculated in prim_cor.


        Rows most correlated with the variable of interest up to the specified percentage are dropped from the dataset.


        Parameters
        ----------
        perc_red : float
            Percentage of data to remove during reduction.


        Attributes
        ----------
        data_reduced : dataframe
            Reduced dataset.

        """
        self.data_reduced = self.data[self.cors < np.percentile(self.cors, float(perc_red))]


    def get_res(self,in_arr):
        """
        Calculates model residuals from which to learn batch effects.


        For time course based designs, this is done with a lowess model.  For block designs, this is done with the block means.


        Parameters
        ----------
        in_arr : arr
            Input array of data from which to calculate residuals.


        Returns
        -------
        res_mat : arr
            Residual matrix.

        """
        def get_l_res(arr):
            """calculates residuals from a lowess model for timecourse designs"""
            return np.array([row - lowess(row, self.tpoints, it=1)[:, 1] for row in arr])

        def get_b_res(arr):
            """calculates residuals from the block mean for block based designs"""
            blocks = np.asarray(self.block_design)
            ma = np.empty_like(arr)
            for v in np.unique(blocks):
                idx = blocks == v
                ma[:, idx] = arr[:, idx].mean(axis=1, keepdims=True)
            return arr - ma

        if self.designtype == 'c':
            res_mat = get_l_res(in_arr)
        elif self.designtype == 'b':
            res_mat = get_b_res(in_arr)
        elif self.designtype == 't':
            res_mat = get_l_res(in_arr)
        return res_mat

    #defined seperately for reuse in perm_test
    def set_res(self):
        """
        Calls get_res()


        Attributes
        ----------
        res : arr
            This is where the residual matrix is stored.

        """

        self.res = self.get_res(self.data_reduced.values)

    def get_tks(self,arr):
        """calculates the fraction of variance for each row of the residual matrix explained by each principle component with PCA"""
        pca = PCA(svd_solver='randomized',random_state=4574)
        pca.fit(arr)
        return pca.explained_variance_ratio_

    #defined seperately for reuse in perm_test
    def set_tks(self):
        """
        Calls get_tks().


        Attributes
        ----------
        tks : list
            This is where the tks (explained variance ratios) are stored.

        """

        self.tks = self.get_tks(self.res)

    def perm_test(self,nperm,npr=1):
        """
        Performs permutation testing on residual matrix SVD.

        The rows of the residual matrix are first permuted.  Then get_tks is called to calculate explained variance ratios and these tks are compared to the values from the actual residual matrix.  A running total is kept for the number of times the explained variance from the permuted matrix exceeds that from the original matrix. Significance is estimated by dividing these totals by the number of permutations.  When npr > 1 a multiprocess Pool is used to parallelise the permutations.
        
        Parameters
        ----------
        nperm : int
            Number of permutations to be tested.
        npr : int
            Number of processors to be used.

        Attributes
        ----------
        sigs : array
            Estimated significances for each batch effect.

        """
        def single_it(rseed):
            rng = np.random.default_rng(rseed * 100)
            rstar = rng.permuted(self.res.copy(), axis=1)
            resstar = self.get_res(rstar)
            tkstar = self.get_tks(resstar)
            n = min(len(tkstar), len(self.tks))
            out = np.zeros(len(self.tks))
            out[:n] = tkstar[:n] > self.tks[:n]
            return out

        seeds = range(int(nperm))
        if int(npr) > 1:
            tpoints = getattr(self, 'tpoints', None)
            block_design = getattr(self, 'block_design', None)
            chunksize = max(1, int(nperm) // (int(npr) * 4))
            init_args = (self.res, self.tks, self.designtype, tpoints, block_design)
            with Pool(int(npr), initializer=_init_perm_worker, initargs=init_args) as pool:
                output = list(tqdm(
                    pool.imap_unordered(_perm_worker, seeds, chunksize=chunksize),
                    total=int(nperm), desc='permuting', smoothing=0,
                ))
        else:
            output = [single_it(s) for s in tqdm(seeds, desc='permuting', smoothing=0)]
        self.sigs = np.sum(np.asarray(output), axis=0) / float(nperm)

    def eig_reg(self,alpha=0.05):
        """
        Regresses eigentrends (batch effects) against the reduced dataset calculating p values for each row being associated with that trend.


        Batch effects are estimated with SVD and only those passing the significance threshold for permutation testing are retained.  These trends are then regressed against each row of the reduced data to estimate a p value for that row's association with that trend.


        Parameters
        ----------
        alpha : float
            Significance threshold for permutation testing, expressed as a decimal.


        Attributes
        ----------
        ps : arr
            P values for association between rows and estimated batch effects.

        """

        alpha = float(alpha)
        U, s, V = np.linalg.svd(self.res)
        sig = V.T[:,:len([i for i in itertools.takewhile(lambda x: x < alpha, self.sigs.copy())])]
        pvals = []
        if len(sig)>0:
            for trend in tqdm(sig.T.copy()):
                temp = []
                for row in self.data_reduced.values.copy():
                    slope, intercept, r_value, p_value, std_err = linregress(row,trend)
                    temp.append(p_value)
                pvals.append(temp)
            self.ps =  pvals
        else:
            print('No Significant Trends')

    def subset_svd(self,lam=0.5):
        """
        Performs SVD on the subset of rows associated with each batch effect.


        The reduced data is first subsetted based on the p values form the regression for each estimated batch effect.  After subsetting the reduced data, SVD is performed on this matrix and the true batch effect is taken as the right singular vector most correlated with the initially estimated batch effect.


        Parameters
        ----------
        lam : float
            P value cutoff above which distribution is assumed to be uniform.  Used to calculate significance cutoff.


        Attributes
        ----------
        ts : arr
            Estimated batch effects (right singular vectors).
        pepts : arr
            Estimates of effect size for each batch effect on each peptide (left singular vectors).

        """

        def est_pi_naught(probs_naught,lam):
            """
            Estimates background distribution of p values.


            Given cutoff lam, calculates ratio of background p values by comparing actual number of p values above cutoff to expected number above cutoff.


            Parameters
            ----------
            probs_naught : list
                List of p values.
            lam : float
                Cutoff above which p values assumed to be drawn from null distribution.


            Returns
            -------
            pi_naught : float
                Estimated ratio of background p values.

            """

            pi_naught = len([i for i in probs_naught if i > lam])/(len(probs_naught)*(1-lam))
            return pi_naught

        def est_pi_sig(probs_sig,l):
            """
            Finds p value pi_sig as cutoff for rows associated with batch effect.


            First estimates ration of background p values.  Then based on this ratio calculates the p value cutoff for rows associated with the batch effect.


            Parameters
            ----------
            probs_sig : list
                List of p values.
            l : float
                Cutoff to be passed to est_pi_naught() as lam.


            Returns
            -------
            pi_sig : float
                Estimated p value cutoff.

            """

            pi_0 = est_pi_naught(probs_sig,l)
            if pi_0 > 1:
                return 'nan'
            sp = np.sort(probs_sig)
            pi_sig = sp[int(np.floor((1-pi_0)*len(probs_sig))-1)]
            return pi_sig

        pt, _, bt = np.linalg.svd(self.res)
        trends = []
        pep_trends = []
        for j, entry in enumerate(tqdm(self.ps)):
            sub = []
            thresh = est_pi_sig(entry,lam)
            if thresh == 'nan':
                self.ts = trends
                self.pepts = pep_trends
                return
            for i in range(len(entry)):
                if entry[i] < thresh:
                    sub.append(self.data_reduced.values[i])
            if len(sub) > 0:
                U, s, V = np.linalg.svd(sub)
                temp = []
                for trend in V:
                    _, _, _, p_value, _ = linregress(bt[j],trend)
                    temp.append(p_value)
                trends.append(V.T[:,np.argmin(temp)])
                pep_trends.append(pt[:,j])
        self.pepts = pep_trends
        self.ts = trends


    def normalize(self,outname):
        """
        Creates diagnostic files, normalizes data based on calculated batch effects, groups peptides by protein and outputs final processed dataset.


        Diagnostic files containing estimated batch effects, explained variance ratios, results of permutation testing and estimated magnitudes of each batch effect on each peptide are generated.  The signal produced by significant batch effects is then estimated and removed from the dataset.  These final results are then written to an output file.


        Parameters
        ----------
        outname : str
            Path to desired output file.


        Attributes
        ----------
        svd_norm : dataframe
            Normalized dataframe with significant batch effects removed.

        """

        pd.DataFrame(self.ts,columns=self.data.columns).to_csv(outname.split('.txt')[0]+'_trends.txt',sep='\t')
        pd.DataFrame(self.sigs).to_csv(outname.split('.txt')[0]+'_perms.txt',sep='\t')
        pd.DataFrame(self.tks).to_csv(outname.split('.txt')[0]+'_tks.txt',sep='\t')
        if len(self.pepts) > 0:
            pd.DataFrame(np.asarray(self.pepts).T,index=self.data_reduced.index).to_csv(outname.split('.txt')[0]+'_pep_bias.txt',sep='\t')
        if len(self.ts) > 0:
            fin_res = np.dot(np.dot(self.data.values,np.linalg.lstsq(self.ts,np.identity(np.shape(self.ts)[0]),rcond=None)[0]),self.ts)
        else:
            fin_res = 0
        self.svd_norm = self.scaler.inverse_transform((self.data.values - fin_res).T).T
        self.svd_norm = pd.DataFrame(self.svd_norm,index=self.data.index,columns=self.data.columns)
        if self.data_type == 'p':
            self.svd_norm = self.svd_norm.groupby(level='Protein').mean()
        self.svd_norm.index.names = ['#']
        self.svd_norm.to_csv(outname,sep='\t')

    def preprocess_default(self):
        self.pool_normalize()
        self.get_tpoints()
        self.prim_cor()
        self.reduce()
        self.set_res()
        self.set_tks()

    def output_default(self,out_file):
        self.eig_reg()
        self.subset_svd()
        self.normalize(out_file)
