# Hope_Graphing_LIET
This repo is dependent on Jacob's LIET repo and is used to graph and analyze LIET output.

**Note on Vocabulary:**
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
  
- trashfunctions.py
  - functions or code that are no longer used but in case you wanted it
  - plotting multiple bar plots in one plot
  - calculating weighted average based on stdev AND probability of parameter based on prior distribution
  
## Analysis Jupyter Notebooks (I chose the actually useful ones lol)

### Important Baseline Analysis

**SharedvsSepPriors10.8**
- comparing mT & sigmaT
- calculating strand bias
- exploratory daata analysis of individual runs
- actual analysis comparing posteriors in both runs

**TPMAnalysis2.0**
- seeing coverage, mT_stdev, and mT_mean comparison between subsamples of a bam file
- **LIET_model_plotting-TPM.ipynb** plots genes of interest from this analysis

**Ru's Genes Analysis**
- seeing how Ru's genes compare to our hermit list and my manually curated list

**Iter_Plotting**
- comparing with 10, 30, and 50 thousand iterations
- I also compare w/ and w/o tI_a and sL_a (NOT without mL_a so don't use as actual comparison for shared vs separate priors)

### Comparing Termination

**GRO&PROComparison**
- Compares mT values between GRO & PRO mcf7 samples

**Term_Gene_comparison_10_9 & LIET_model_plotting-TPM**
- Comparing termination across celltypes (within & between experiments)
- This actually uses the padded data lol

**HeatShock1.0**
- comparing heatshock results for Vih 2017 
- This is the one where the heatshock bedgraphs didn't actually look like heatshock graphs.
- Therefore, there is a lot of frantic exploratory data analysis to assess what is happening before I looked at the bedgraphs.

**HeatShock2.0**
- comparing heatshock results for Joe's data for Eric

### Random Other

**RNA_recycling**
- EARLY assessment so excuse the older, less efficient functions
- calculating background coverage for antisense and sense
- seeing if termination peaks correlate with antisense initiation peaks
- Summary: no correlation observed, but seeing if the same RNA pol might be used again to retranscribe the same gene and/or might be transcribing the antisense portion
