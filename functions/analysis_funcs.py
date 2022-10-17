import pandas as pd
import plotly.express as px
from statistics import stdev as statstdev
#import pymc3 as pm
# ============ WORKING WITH DATAFRAMES
def df_merger(df1, df2, name):
    df1 = df1.filter(items = df2.index, axis=0)
    df2 = df2.filter(items = df1.index, axis=0)
    print("After merging " + name + ", # genes= " + str(len(df1)))
    return df1, df2

def df_combiner(col_list, df_list, name_list):
    comb_df = df_list[0]['gene']
    comb_df = pd.DataFrame(comb_df)
    comb_df['strand'] = df_list[0]['strand']
    for y_value in col_list:
        count = 0
        for df, name in zip(df_list, name_list):
            if count != 0:
                # use only same genes
                comb_df, df = df_merger(comb_df, df, name)
            comb_df[str(name)+'_'+str(y_value)] = df[y_value]
            # increase count
            count += 1
    return comb_df

def filter_dfs(df_list, col, value, comparison):
    """
    This function takes a list of dataframes and filters them based 
    on a comparsion. It can also name them.
    Parameters:
    - df_list: list of dataframes
    - col: string of column wanting to filter by
    - value: value by which to filter
    - comparison: string of >, <, ==, !=
    Returns:
    - filtered dataframes in order of df_list
    """
    new_df_list = []
    for df in df_list:
        if comparison == '>':
            df = df.loc[(df[col] > value)]
        elif comparison == '<':
            df = df.loc[(df[col] < value)]
        elif comparison == '==':
            df = df.loc[(df[col] == value)]
        elif comparison == '!=':
            df = df.loc[(df[col] != value)]
        new_df_list.append(df)
    return new_df_list

# ============ CALCULATING. CORRELATION.

def getcov(df,name='', log=False):
    '''This function allows coverage to be calculated for 
    each gene based on if it is on the positive or negative strand. 
        Eg. If it is a positive gene the cov is the pos_cov.
    Parameters:
       - df: FitParse dataframe
    Optional:
       - name: if the columns with the pos_cov & neg_cov have a prefix put as string
           ie: if Stein 22_pos_cov then name='Stein 22'
       - log: If you want the coverage to be log2, 
           True will make a new column called log2_cov
           default: False
     Returns: new dataframe
        '''
    cov_list = []
    for gene in df.index:
        if df.loc[gene]['strand'] == 1:
            cov_list.append(df.loc[gene][name+'pos_cov'])
        elif df.loc[gene]['strand'] == -1:
            cov_list.append(abs(df.loc[gene][name+'neg_cov']))
    df['cov'] = cov_list
    if log == True:
        df['log2_cov'] = np.log2(df['cov'])
    return df

def calc_strand_bias(df):
    '''This function allows calculation of strand bias (favors antisense 0<x<1 favors sense)."    
    Parameters: FitParse dataframe
    Returns: dataframe with the following columns:
       Copied from original: w_LI_mean, w_aLI_mean, pos_cov, neg_cov
       New:
        - sense init: w_LI_mean * pos_cov
        - antisense init: w_aLI_mean * neg_cov
        - total init: sense_init + antisense_init
        - strand bias: sense_init/(total_init)'''
    df = df[['gene','w_LI_mean', 'w_aLI_mean', 'pos_cov', 'neg_cov']]
    df['sense init'] = df['w_LI_mean']*df['pos_cov']
    df['antisense init'] = abs(df['w_aLI_mean']*df['neg_cov'])
    df['total init'] = df['sense init'] + df['antisense init']
    df['strand bias'] = (df['sense init'])/(df['total init'])
    return df

def calc_weighted_average(df_list, factor, average=[]):
    '''This function calculates the weighted average of the parameter of interest
    via its standard deviation.
    Parameters:
    - df_list: list of dataframes
    - factor: "prior" of interest (must have _mean & _stdev in df)
    - average: (OPTIONAL) list of columns you waant to get average of
    Returns:
    - df with the following columns:
        - mean/weight sum (sum(mean/stdev))
        - weight sum (sum(stdev))
        - weighted average ((mean/weight sum)/(weight sum))
        - averaages of values in averagee column
    '''
    # create a new df for calcs
    df2= df_list[0]['gene']
    df2 = pd.DataFrame(df2)
    # iterate through all dataframes in list
    count = 0
    for df in df_list:
        # if starting out
        if count == 0:
            df2['mean/weight sum'] = (df[str(factor + '_mean')]/df[str(factor + '_stdev')])
            df2['weight sum'] = df[str(factor + '_stdev')]    
            for value in average:
                df2[value] = df[value]
        else:
            # need to make sure new one only contains genes that df2 has
            df = df[df.index.isin(df2.index)]
            df2 = df2[df2.index.isin(df.index)]
            # get the proper sums
            df2['mean/weight sum'] += (df[str(factor + '_mean')]/df[str(factor + '_stdev')])
            df2['weight sum'] += df[str(factor + '_stdev')]
            if average != []:
                for value in average:
                    # sum up values
                    df2[value] = df2[value]+df[value]
        count += 1
    # normalize mean/weight sum (* inverse weight sum)
    df2['weighted average'] = df2['mean/weight sum']/((df2['weight sum']))
    # at end get averages
    if average != []:
        for value in average:
            df2[value] = df2[value]/(count)
    return df2



def calc_average(df_list, average=[], stdev=''):
    '''This function calculates the weighted average of the parameter of interest
    via its standard deviation.
    Parameters:
    - df_list: list of dataframes
    - factor: "prior" of interest (must have _mean & _stdev in df)
    - average: (OPTIONAL) list of columns you waant to get average of
    Returns:
    - df with the following columns:
        - mean/weight sum (sum(mean/stdev))
        - weight sum (sum(stdev))
        - weighted average ((mean/weight sum)/(weight sum))
        - averaages of values in averagee column
    '''
    # create a new df for calcs
    df2= df_list[0]['gene']
    df2 = pd.DataFrame(df2)
    # create a new df for stdev calcs (recognize not most efficient)
    df3 = df_list[0]['gene']
    df3 = pd.DataFrame(df3)
    
    # iterate through all dataframes in list
    count = 0
    stdev_list = []
    stdev_name = stdev+'_stdev'
    for df in df_list:
        # if starting out
        if count == 0:   
            for value in average:
                df2[value] = df[value]
        else:
            # need to make sure new one only contains genes that df2 has
            df = df[df.index.isin(df2.index)]
            df2 = df2[df2.index.isin(df.index)]
            # get the proper sums
            if average != []:
                for value in average:
                    # sum up values
                    df2[value] = df2[value]+df[value]
        if stdev != '':
            # make new df with values of each to get stdev of
            df3 = df3[df3.index.isin(df2.index)]
            stdev_name_iter = stdev + '_' + str(count)
            df3[stdev_name_iter] = df[stdev]  
         
        count += 1
    # at end get averages
    if average != []:
        for value in average:
            df2[value] = df2[value]/(count)
    # at end get stdev with stdev list
    if stdev != '':
        df3.pop("gene")
        stdev_calc = df3.std(axis=1)
        stdev_name = stdev+'_stdev'
        df2[stdev_name] = stdev_calc
    return df2

def corr_matrix_provider2(dfs, y_value="", corrtype = 'pearson'):
    '''
    This function takes a list of dfs and performs correlation tests among all combos
    Parameters:
    - dfs: list of dfs
    - y_value: column you want to compare (must be shared among all dfs)
    Optional parameters:
    - corrtype: pearson or spearman
    Returns:
    - correlation df 
    '''
    for x in range(0, len(dfs)):
        name = dfs[x].name
        if x == 0:
            # make initial df
            df = pd.DataFrame(dfs[x][y_value])
            df.columns = [name]
        else:
            # fix df to have only indexes in next df
            df, dfs1 = df_merger(df, dfs[x], name)
            assert df.index.tolist()==dfs1.index.tolist()
            df = pd.DataFrame(df)
            df[name] = dfs1[y_value]
    if corrtype == 'pearson':
        return df.corr(method='pearson')
    elif corrtype == 'spearman':
        return df.corr(method='spearman')
    else:
        print("Improper corrtype provided: please choose pearson or spearman")
        
# Get plotly correlation graph
def corr_matrix_plotter(dfs, y_value=""):
    '''
    This function takes a list of dfs and produces a correlation scatter matrix via plotly and colors by gene.
    Parameters:
    - dfs: list of dfs
    - y_value: column you want to compare (must be shared among all dfs)
    Returns:
    - plots scatterplotmatrix colored by gene
    - returns dataframe that has gene as one column & y_value for other

    '''
    dimensions = []
    for x in range(0, len(dfs)):
        name = dfs[x].name
        dimensions.append(name)
        if x == 0:
            # make initial df
            df = pd.DataFrame(dfs[x][["gene", y_value]])
            df.columns = ["gene", name]
        else:
            # fix df to have only indexes in next df
            df, dfs1 = df_merger(df, dfs[x], name)
            assert df.index.tolist()==dfs1.index.tolist()
            df = pd.DataFrame(df)
            df[name] = dfs1[y_value]
    fig = px.scatter_matrix(df, dimensions = dimensions, color="gene")
    fig.show()
    return df