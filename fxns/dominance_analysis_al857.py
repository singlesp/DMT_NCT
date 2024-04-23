#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 12:12:44 2022

@author: al857
"""
# dominance_classification=Dominance(data=df_data,target='Target',objective=0,pseudo_r2="mcfadden",data_format=0)
# cfr https://dominance-analysis.readthedocs.io

# data - Complete Dataset, should be a Pandas DataFrame

# target - Name of the target variable, it should be present in passed dataset.

# top_k - No. of features to choose from all available features. By default, the package will run for top 15 features.

# objective - It can take the values 0 or 1.
# 0 for Classification
# 1 for Regression
# By default, the package will run for Regression.

# pseudo_r2 - It can take one of the Pseudo R-Squared measures - mcfadden, nagelkerke , cox_and_snell or estrella, where default = mcfadden. It is not needed in the case of regression models (i.e. objective = 1 ).

# data_format - It can take the values 0, 1 or 2.
# 0 when raw data is being passed,
# 1 when correlation matrix (correlation of predictors with target variable) is being passed,
# 2 when covariance matrix (covariance of predictors with target variable) is being passed.
# By default, the package will run for raw data (i.e. data_format = 0). This parameter is not needed in case of Classification models.

# Note: While passing a Covariance / Correlation matrix to dominance-analysis , it is advisable to pass the matrix of a dataset having 15 or lesser predictor variables

# Outputs
# Individual Dominance - The individual dominance of a predictor is the R2 of the model
#  between the dependent variable and the predictor. Hence, individual dominance can be
#  interpreted as the variability explained by the predictor alone 

# Interactional Dominance - This is the incremental R2 contribution of the predictor 
# to the complete model. Hence, the Interactional Dominance of a particular predictor 'X'
#  will be the diffrence between the R2 of the complete model and the R2 of the model 
#  with all other predictors except the particular predictor 'X'

# Average Partial Dominance - This is average of average incremental R2 contributions of
#  the predictor to all subset models

# Total Dominance - The last measure of dominance summarizes the additional contributions of each predictor to all subset models



###############################
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
from dominance_analysis import Dominance_Datasets
from dominance_analysis import Dominance
import sys
import os



def dominance_analysis_al857(input_filename, target_var, num_predictors, input_format, output_filename):
    #Example:
    #boston_dataset=Dominance_Datasets.get_boston()
    #corr_data = boston_dataset.corr()
    #corr_data
    #dominance_regression=Dominance(data=corr_data,target='House_Price',data_format=1)
    
    
    myData = pd.read_csv(input_filename)
    
    dominance_regression=Dominance(data=myData,target=target_var, top_k=num_predictors, data_format=input_format)
    
    R2=dominance_regression.incremental_rsquare()
    #dominance_regression.plot_incremental_rsquare()
    DS=dominance_regression.dominance_stats()
    DL=dominance_regression.dominance_level()
    
    #save
    #R2.to_csv(output_filename + "_incrementalR2.csv", index=False)
    DS.to_csv(output_filename + "_DominanceStats.csv")#, index=False)
    DL.to_csv(output_filename + "_DominanceLevels.csv")#, index=False)
    return(dominance_regression)

if __name__ == "__main__":    
    input_filename = str(sys.argv[1])    
    target_var = str(sys.argv[2])    
    num_predictors = int(str(sys.argv[3]))    
    input_format = int(str(sys.argv[4]))    
    output_filename = str(sys.argv[5])
    dominance_analysis_al857(input_filename, target_var, num_predictors, input_format, output_filename)
    
