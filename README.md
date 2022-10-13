# Hope_Graphing_LIET
This repo is dependent on Jacob's LIET repo and is used to graph and analyze LIET output.
*Note on Vocabulary:*
- Runs/Samples: e.g. SRR86293022
- Papers/Experiments: e.g. Core2014Analysis
## Summary Notes
1. I changed the Fitparse class to return a dataframe **additional** to the current dictionaries with the following columns (in ''):
        
        
            Columns: 
        Annotations: 'chrom', 'start', 'stop', 'strand' 
        Sense Main Priors: 'mL_mean', 'mL_stdev', 'sL_mean','sL_stdev', 'tI_mean', 
                           'tI_stdev', 'mT_mean', 'mT_stdev', 'sT_mean','sT_stdev' 
        Sense Weight Priors:'w_LI_mean', 'w_LI_stdev', 'w_E_mean', 'w_E_stdev',
                            'w_T_mean', 'w_T_stdev', 'w_B_mean', 'w_B_stdev' 
        Antisense Main Priors:'mL_a_mean','mL_a_stdev', 'sL_a_mean', 'sL_a_stdev' 
        Antisense Weight Priors:'w_aLI_mean', 'w_aLI_stdev','w_aB_mean', 'w_aB_stdev'
        Log Info: 'pos_cov', 'neg_cov', 'total_cov', 'elbo_lrange', 'elbo_urange', 'fit_time_min' 
        
        
2. All of my graphing uses dataframes since it allows easy merging of different datasets for comparison. This is mediated by the functions in analysis.py. The two of greatest note would be the df_merger which would return the dataframes to ensure the new indexes only include genes shared between the two runs. df_combiner takes dataframes of different samples and gets the columns of interest and adds them to a single dataframe labeled appropriately for comparsion. It uses the merger function to ensure only genes shared by all samples are included.
3. I tried to make the functions modular but if I want to change a tiny nuance in a function, I'll usually just copy and paste a function into the jupyter notebook and use that function rather than the imported one.


## Shared Functions
- analysis_funcs.py
  - Working with dataframes:
    - merging, combining, and filtering dataframes
  - Calculating:
    - calculate strand bias, averages, weighted averages, correlation matrices, coverage

- plotting_funcs.py
  - FitParse class with dataframe added
  - Individual gene plotting among multiple samples
  - Bar & boxplot graphs (use df merger before)
  - Scatter plots with the correlation coefficients (single dataframe or multiple dataframes)
  
## Analysis Jupyter Notebooks
