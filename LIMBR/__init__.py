"""
The package LIMBR consists of three modules.  imputation performs imputation of missing data using K nearest neighbors.  batch_fx performs modelling and removal of batch effects. old_fashioned provides normalization methods without
batch effect removal as a baseline for the results of batch_fx.
"""

from LIMBR import imputation
from LIMBR import batch_fx
from LIMBR import old_fashioned
from LIMBR import simulations
