LIMBR: Learning and Imputation for Mass-spec Bias Reduction
===========================================================

LIMBR provides a streamlined tool set for imputation of missing data followed by modelling and removal of batch effects.  The software was designed for proteomics datasets, with an emphasis on circadian 
proteomics data, but can be applied to any time course or blocked experiments which produce large amounts of data, such as RNAseq. The two main classes are imputable, which performs missing data imputation, and sva, which performs 
modelling and removal of batch effects.

----------
Motivation
----------

Decreasing costs and increasing ambition are resulting in larger Mass-spec (MS) experiments.  MS experiments have a few limitations which are exacerbated by this increasing scale, namely batch effects and missing data.  Many 
downstream statistical analyses require complete cases for analysis, however, MS produces some missing data at random meaning that as the number of experiments increase the number of peptides rejected due to missing data actually 
*increases*.  This is obviously not good, but fortunately there is a solution!  If the missing data for observations missing only a small number of data points are imputed this issue can be overcome and that's the first thing that 
LIMBR does.  The second issue with large scaler MS experiments is batch effects.  As the number of samples increases, the number of batches necessary for sample processing also increases.  Batch effects from sample processing are 
known to have a large effect on MS data and increasing the number of batches means more batch effects and a higher proportion of observations affected by at least one batch effect.  Here LIMBR capitolizes on the larger amount of 
data and the known correlation structure of the data set to model these batch effects so that they can be removed.

--------
Features
--------

*KNN based imputation of missing data.

*SVD based modelling and removal of batch effects.

*Built for circadian and non-circadian time series as well as block designs

------------
Code Example
------------

------------
Installation
------------

-----------
How to Use?
-----------

-------
Credits
-------

The sva based methods build on work for micro-array datasets by Jeffrey Leek, with particular reliance on his PhD Thesis from the University of Washington.  Other contributions in the field include EigenMS.
