import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import re
import sys
from collections import OrderedDict

#======================== Ind Gene Plotting (ensure comb,shared) =======
def meanerror_plotter(full_df, gene, par, cats, group_name):
    '''
    This function plots the mean and stdev of parameter par for ONE gene 
    from a df with all applicable info and color codes/labels based on cats.
    In this case, cats should be experiments within a single cell type.
    Pars:
    - full_df: df of all experiments of choice (within single cell type)
    - gene of choice: string
    - par: parameter of interest (must be in cols as par_stdev & par_mean)
        - ex. 'mT'
    - cats: list of categories by which legend & coloring is based (experiments)
        - ex. ['And', 'Li21', 'Li13']
    - group_name: string celltype
    '''
    # make the lists to plot
    error = [[],[],[],[],[],[]]
    name_list = [[],[],[],[],[],[]]
    par_list = [[],[],[],[],[],[]]

    for col in full_df.columns:
        if str(par+'_mean') in col:
            name = col.split('_')[0]
            # name = col.split('_')[0].split(' ')[1]
            for i in range(0,len(cats)):
                if name.startswith(cats[i]):
                    par_list[i].append(full_df.at[gene,col])
                    name_list[i].append(col.split('_')[0].split(' ')[1])
        elif str(par+'_stdev') in col:
            name = col.split('_')[0]
            # name = col.split('_')[0].split(' ')[1]
            for i in range(0,len(cats)):
                if name.startswith(cats[i]):
                    error[i].append(full_df.at[gene,col])
    fig, ax = plt.subplots()
    
    for i in range(0,len(cats)):
        ax.errorbar(x=name_list[i], y=par_list[i], yerr=error[i], fmt='o', linewidth=2, capsize=10, label=cats[i])
    #plt.ylim([-100, 800])
    ax.set_ylabel(par)
    ax.set_title(par+' means of '+gene+' for '+group_name +' experiments with stdev')
    ax.legend()
    
def meanerror_exp_plotter(full_df, gene, par, cats, group_name,
                          expline=True, fontsize=12, coverage=True):
    '''
    This function plots the mean and stdev of parameter par for ONE gene 
    from a df with all applicable info and color codes/labels based on cats.
    In this case, cats should be celltypes. It also adds black lines between 
    experiments if line=true. Default is true.
    Pars:
    - full_df: df of all experiments of choice
        - col names should be "celltype papername #'
    - gene of choice: string
    - par: parameter of interest (must be in cols as par_stdev & par_mean)
        - ex. 'mT'
    - cats: list of categories by which legend & coloring is based (celltypes)
        - ex. ['hct116', 'g401', 'lymph']
    - group_name: string celltypes
    OPTIONAL PARS:
    - expline: DEFAULT is true, whether or not want black lines distinguishing experimnets
    - fontsize: DEFAULT is 12, size of x labels
    - coverage: DEFAULT IS TRUE, whether or not want coverage also plotted on graph
    '''
    # make the lists to plot
    error_list = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    name_list = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    par_list = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    cov_list = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

    for col in full_df.columns:
        # get mean
        if str(par+'_mean') in col:
            # get the number (exp & run #)
            name = col.split('_')[0].split(' ')[1:]
            # for each category (cell type)
            for i in range(0,len(cats)):
                # if the column is part of the cell type
                if col.startswith(cats[i]):
                    # add that mean of the par to the par_list
                    par_list[i].append(full_df.at[gene,col])
                    # add that name (exp & run #) to the name_list
                    name_list[i].append(''.join(name))
        # get standard deviation
        elif str(par+'_stdev') in col:
            name = col.split('_')[0].split(' ')[1]
            for i in range(0,len(cats)):
                if col.startswith(cats[i]):
                    # add the stdev of the par to the error list
                    error_list[i].append(full_df.at[gene,col])
        if 'cov' in col:
            for i in range(0,len(cats)):
                if col.startswith(cats[i]):
                    # add the cov to the error list
                    cov_list[i].append(full_df.at[gene,col])
    fig, ax = plt.subplots()
    if coverage == True:
        ax2 = ax.twinx()
    # make a variable to hold expeeriment
    exp = ''
    # make a variable to hold the # for the black line placement
    lineplace = 0
    colordict = {}
    color_ops = [  'dodgerblue','orange', 'darkseagreen', 'red', 'magenta', 'gold',
                 'turquoise', 'purple', 'brown', 'pink','blue',]
    for i in range(0,len(cats)):
        colordict[cats[i]] = color_ops[i]
    for i in range(0,len(cats)):
        for ername, mean, error, cov in zip(name_list[i], par_list[i], error_list[i], cov_list[i]):
            if ername != []:
                # add 1 to lineplacement
                lineplace += 1
                # graph the error bars for each cell type (color changes based on it)
                ax.errorbar(x=ername, y=mean, yerr=error, fmt='o', linewidth=1.5, capsize=5, 
                            color = colordict[cats[i]], label = cats[i])
                if coverage == True:
                    #ax2.plot(x=ername, y=cov, color = 'gray')
                    ax2.errorbar(x=ername, y=cov, yerr=error*0, fmt='.', capsize=0, 
                            color = 'black')
                if expline == True:
                    # if the experiment changes
                    actualexp = re.sub(r'[^a-zA-Z]','',ername)
                    if actualexp != exp:
                        if exp!='':
                            ax.axvline(lineplace-1.5,color='black',linewidth=0.2)
                        # change exp value
                        exp = actualexp

    ax.set_ylabel(par, color='dodgerblue')
    if coverage == True:
        ax2.set_ylabel('coverage', color='black')
    ax.set_title(par+' means of '+gene +' in '+group_name +' experiments with stdev')
    markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in colordict.values()]
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True)
    ax.legend(markers, colordict.keys(), numpoints=1, 
              bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=fontsize)
    

# def meanerror_exp_plotter(full_df, gene, par, cats, group_name):
#     '''
#     This function plots the mean and stdev of parameter par for ONE gene 
#     from a df with all applicable info and color codes/labels based on cats.
#     In this case, cats should be celltypes.
#     Pars:
#     - full_df: df of all experiments of choice
#     - gene of choice: string
#     - par: parameter of interest (must be in cols as par_stdev & par_mean)
#         - ex. 'mT'
#     - cats: list of categories by which legend & coloring is based (celltypes)
#         - ex. ['hct116', 'g401', 'lymph']
#     - group_name: string celltypes
#     '''
#     # make the lists to plot
#     error = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
#     name_list = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
#     par_list = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

#     for col in full_df.columns:
#         if str(par+'_mean') in col:
#             # get the S81 (exp & run #)
#             name = col.split('_')[0].split(' ')[1]
#             # for each category (cell type)
#             for i in range(0,len(cats)):
#                 # if the column is part of the cell type
#                 if col.startswith(cats[i][0]):
#                     # add that mean of the par to the par_list
#                     par_list[i].append(full_df.at[gene,col])
#                     # add that name (S81) to the name_list
#                     name_list[i].append(col.split('_')[0].split(' ')[1])
#         elif str(par+'_stdev') in col:
#             name = col.split('_')[0].split(' ')[1]
#             for i in range(0,len(cats)):
#                 if col.startswith(cats[i][0]):
#                     # add the stdev of the par to the error list
#                     error[i].append(full_df.at[gene,col])
#     fig, ax = plt.subplots()
#     # make a variable to hold expeeriment
#     exp = ''
#     for i in range(0,len(cats)):
#         # graph the error bars for each cell type (color changes based on it)
#         ax.errorbar(x=name_list[i], y=par_list[i], yerr=error[i], fmt='o', linewidth=2, capsize=10, label=cats[i])
#     ax.set_ylabel(par)
#     ax.set_title(par+' means of '+gene +' in '+group_name +' experiments with stdev')
#     ax.legend()
    
#======================== Plotting (ensure shared before)  ======================
def one_df_bar_plotter(x_df, y_df, x_label, y_label, title=None, sort=True, custom_segment=None):
    '''This function plots a horizontal bargraph of the x_df column against
    the y_df column with options to sort & skip rows.
    Required Parameters:
        - x_df & y_df: dataframe columns that must have the same index size
            that will be plotted (ex. df.index & df['y'])
        - x_label & y_label: String of the desired labels on the graph.
            These do NOT need to be the same labels used in the dfs.
    Optional Parameters:
        - sort: Boolean of whether or not the graph will sort the y_values
           from lowest to highest when graphing.
           - DEFAULT: True
        - custom_segment: string & be in format ::,: or :,:'
        - title: string title of plot
        - y_axis chosen:
    '''
    '''**IMPORTANT: MUST ENSURE only shared genes B4 USE'''
    df = pd.DataFrame({str(x_label): x_df, str(y_label): y_df})
    x_row, y_row, x_col, y_col, skip_rows = None, None, None, None, None
    if sort == True:
        df = df.sort_values(by=y_label)
    if custom_segment is not None:
        seg = custom_segment.split(',')
        seg_row = seg[0].split(':')
        seg_col = seg[1].split(':')
        if seg_row[0]: 
            x_row=int(seg_row[0])
        if seg_row[1]:
            y_row = int(seg_row[1])
        if seg_row[2]:
            skip_rows = int(seg_row[2])
        if seg_col[0]:
            x_col = int(seg_col[0])
        if seg_col[1]:
            y_col = int(seg_col[1])
        df = df.iloc[x_row:y_row:skip_rows, x_col:y_col]
    ax = df.plot.bar(x=x_label, y=y_label, rot=90, figsize=(20,6))
    if title:
        ax.set(title=title)

def box_plotter_all(datadict, title=None, y_label=None, x_label=None, showmeans=True):
    '''This function plots multiple boxplots on a single figure based on the datadict with ALL datapoints.
        Required Parameters:
          - datadict: {"x_label": data} data can be array, or dataframe['column']
        Optional Parameters
          - showmeans: boolean: default=True
          - title: string
          - y_label: string
          - x_label: string
    '''
    fig = plt.figure(figsize=(7.50, 3.50))
    #fig.suptitle('bold figure suptitle', fontsize=14, fontweight='bold')
    ax = fig.add_subplot()
    # get data & plot
    data = pd.DataFrame(datadict)
    ax = data.boxplot(showfliers=False, showmeans=showmeans)
    # get the scatterplots for each dataset
    for i, d in enumerate(data):
       y = data[d]
       x = np.random.normal(i + 1, 0.04, len(y))
       plt.scatter(x, y, alpha=0.4)
    # other parameters
    ax.set_title(title)
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    plt.show()
    
def box_plotter(data, labels, title=None, y_label=None, x_label=None, showmeans=True, meanpointprops=None):
    '''This function will take data in the form of an array and plot boxplots without all datapoints.
    Required parameters:
      - data: form of an array or array of arrays (including dataframe arrays)
      - labels: array of strings for labels (must correspond in index with array of arrays)
    Optional parameters
      - showmeans: boolean: default=True
      - title: string
      - y_label: string
      - x_label: string
      - meanpointprops: dict of properties ex: dict(marker='D', markeredgecolor='black', markerfacecolor='firebrick')
    '''
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot()
    medianprops = dict(linestyle='-', color='black')
    meanpointprops = meanpointprops
    ax.boxplot(data, labels=labels, showmeans=showmeans, patch_artist=True, 
           medianprops= medianprops, meanprops=meanpointprops)
    ax.set_title(title)
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    plt.show()
    
    
#======================== Plotting (checks shared in func)  ================================
def dataframe_plotter(dfx, dfy, y_value, xlabel, ylabel, title, corr=True):
    '''This function plots the y_values (list) between each of the runs. It also shows the NONLINEAR correlation (Spearman's R^2) between the 2 if corr=True'''
    # only compare where have same genes
    dfx = dfx[dfx.index.isin(dfy.index)]
    dfy = dfy[dfy.index.isin(dfx.index)]
    # plot y_value of both runs against each other
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.scatter(dfx[y_value], dfy[y_value], edgecolors='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if corr==True:
        # calculate Spearman's correlation
        corr_value, _ = spearmanr(dfx[y_value], dfy[y_value])
        plt.title(title + '(R^2: ' + str(round(corr_value,4)) + ')'+'p='+str(_))

        
def corr_grapher(x,y, y_value=""):
    '''Same as dataframe_plotter but if more than 2.'''
    compname_list = []
    corr_list = []
    for dfx in x:
        for dfy in y:
            key = str(dfx.name) + ' vs ' + str(dfy.name)
            # only compare where have same genes
            dfx1 = dfx[dfx.index.isin(dfy.index)]
            dfy1 = dfy[dfy.index.isin(dfx.index)]
            # calculate Pearson's correlation
            corr_value, _ = pearsonr(dfx1[y_value], dfy1[y_value])
            compname_list.append(key)
            corr_list.append(corr_value)
    plt.plot(compname_list, corr_list)
    plt.xticks(rotation = 90)

    
    

class FitParse:
    '''
    FitParse loads the results of a LIET fitting run from the standard output 
    file. Results are organzed into several dictionaries and lists for easier 
    parsing. Additionally, the class contains methods for parsing the main 
    results dictionary <fits>.
    '''

    def __init__(self, res_file, log_file=None, antisense=True, ET=True, colon_format=False):

        self.definitions = OrderedDict({
            "mL": "Sense strand loading position (mu)",
            "sL": "Sense strand loading stdev (sigma)",
            "tI": "Sense strand initiation length (tau)",
            "mT": "Sense strand termination position (mu) relative to TCS",
            "sT": "Sense strand termination stdev (sigma)",
            "w": "Sense strand weights [load, elong, terminate, background]",
            "mL_a": "Antisense strand loading position (mu)",
            "sL_a": "Antisense strand loading stdev (sigma)",
            "tI_a": "Antisense strand initiation length (tau)",
            "w_a": "Antisense strand weights [load, background]",
        })

        self.antisense = antisense
        self.ET = ET
        self.genes = []
        self.annotations = OrderedDict()
        self.fits = OrderedDict()
        self.percentiles = None
        self.priors_list = ['mL','sL', 'tI', 'mT','sT','w_LI','w_E', 'w_T', 'w_B', 'mL_a','sL_a', 'w_aLI','w_aB']

        with open(res_file, 'r') as rf:

            for line in rf:
                # Iterate through header
                if line[0] == '#':
                    line_list = line[1:].strip().split()
                    if line_list[0] == "CONFIG":
                        self.config = line_list[1]
                    continue
                else:
                    pass
                
                # Check line (gene) has a fit result
                line_list = line.strip().split('\t')
                if len(line_list) != 6:
                    if "error" in line_list:
                        print(line_list, file=sys.stderr)
                        continue
                    else:
                        print(f"CHECK LINE: {line_list}", file=sys.stderr)
                        continue
                
                # Parse and cast line
                chrom, start, stop, strand, gid, fit = line_list
                start = int(start)
                stop = int(stop)
                strand = int(strand)

                self.genes.append(gid)
                self.annotations[gid] = (chrom, start, stop, strand)

                # Fit line format: param1=val1:err1,param2=val2:err2,...
                # Parse all the best fit parameter values (CSV string)
                temp = OrderedDict()
                for param_val in fit.strip().split(','):
                    # Parameter name and its value
                    p, v = param_val.strip().split('=')
                    # param is Percentiles, need to change data stripping
                    if p == "Percentiles":
                        self.percentiles = True
                        temp["Percentiles"]= v
                    else:
                        # Mean and standard error of value (from posterior dist)
                        v_m, v_s = v.split(':')
                        if p in ['w', 'w_a']:
                            v_m = [float(i) for i in v_m.strip('[]').split()]
                            v_s = [float(i) for i in v_s.strip('[]').split()]
                        else:
                            v_m = float(v_m)
                            v_s = float(v_s)
                        temp[p] = (v_m, v_s)
                
                self.fits[gid] = temp

        # Extract and assign all the variable arrays
        self.mL, self.mL_std = self.param_extract('mL', stdev=True)
        self.sL, self.sL_std = self.param_extract('sL', stdev=True)
        self.tI, self.tI_std = self.param_extract('tI', stdev=True)
        
        # Recalculate mT values so they are relative to end of annotation
        if ET is True:
            absolute_mT, self.mT_std = self.param_extract('mT', stdev=True)
            relative_mT = []
            for i, gene in enumerate(self.genes):
                tss = self.annotations[gene][1]
                tcs = self.annotations[gene][2]
                diff = abs(tcs - tss)
                relative_mT.append(absolute_mT[i] - diff)
            self.mT = relative_mT

            self.sT, self.sT_std = self.param_extract('sT', stdev=True)
        self.w, self.w_std = self.param_extract('w', stdev=True)
        
        if antisense is True:
            self.mL_a, self.mL_a_std = self.param_extract('mL_a', stdev=True)
            self.sL_a, self.sL_a_std = self.param_extract('sL_a', stdev=True)
            self.tI_a, self.tI_a_std = self.param_extract('tI_a', stdev=True)
            self.w_a, self.w_a_std = self.param_extract('w_a', stdev=True)
        
        ## make the unique weighted lists
#         self.priors_list = ['mL','sL', 'tI', 'mT','sT','w_LI','w_E', 'w_T', 'w_B', 'mL_a','sL_a', 'w_aLI','w_aB']
        if ET is True:
            self.w_LI, self.w_E, self.w_T, self.w_B  = self.weight_par_extract(self.w)
            self.w_LI_std, self.w_E_std, self.w_T_std, self.w_B_std = self.weight_par_extract(self.w_std)
        else:
            self.w_LI, self.w_B  = self.weight_par_extract(self.w)
            self.w_LI_std, self.w_B_std = self.weight_par_extract(self.w_std)
        if antisense is True:
            self.w_aLI, self.w_aB = self.weight_par_extract(self.w_a)
            self.w_aLI_std, self.w_aB_std = self.weight_par_extract(self.w_a_std)
#         else:
#             self.w_aLI, self.w_aE, self.w_aT, self.w_aB  = self.weight_par_extract(self.w_a)
#             self.w_aLI_std, self.w_aE_std, self.w_aT_std, self.w_aB_std = self.weight_par_extract(self.w_a_std)
        
        
        
        # Parse strand coverage and min/max elbo values from log file
        if log_file:
            self.log = OrderedDict()
            with open(log_file, 'r') as lf:
                for line in lf:
                    if line[0] == '#':
                        continue
                    elif line[0] == '>':
                        line = line.strip().split(':')
                        if colon_format:
                            gene_id = ":".join([line[0][1:], line[1]])
                        else:
                            gene_id = line[0][1:]
                        self.log[gene_id] = dict()
                    else:
                        field, value = line.strip().split(':')
                        if field == 'strand_cov':
                            value = value.strip('()').split(',')
                            value = tuple(map(int, value))
                        elif field == 'elbo_range':
                            value = value.strip('()').split(',')
                            value = tuple(map(float, value))
                        elif field == 'fit_time_min':
                            value = float(value)
                        else:
                            continue
                        self.log[gene_id].update({field: value})
            
            self.cov_pos = []
            self.cov_neg = []
            print(self.genes)
            for g in self.genes:
                #if (self.log[g]['strand_cov']):
                    self.cov_pos.append(self.log[g]['strand_cov'][0])
                    self.cov_neg.append(abs(self.log[g]['strand_cov'][1]))

        # Get all info in a dataframe
        self.df = self.dataframe_creator()

    
    def param_extract(self, p, stdev=False):
        '''
        Extract param (p) values (ordered based on genes list) from fits 
        dictionary and output them to a list.
        '''
        param_vals = []
        if stdev:
            param_stdev = []

        for g in self.genes:
            param_vals.append(self.fits[g][p][0])
            if stdev:
                param_stdev.append(self.fits[g][p][1])

        if stdev:
            return param_vals, param_stdev
        else:
            return param_vals
    
    def percentile_extract(self):
        """
        Extract the percentiles (ordered baased on genes list) from fits
        dictionary and output as a dictionary with the different percentiles as keys
        and the lists as the lists in gene order.
        """
        percentile_dict = dict()
        # save each percentile as per_pos and per_neg
        first = True
        for g in self.genes:
            if first:
                for perc in self.fits[g]["Percentiles"].split(";"):
                    p_new, pos_res, neg_res = perc.split("-")
                    percentile_dict["_".join([p_new, "perc_pos"])] = [pos_res]
                    percentile_dict["_".join([p_new, "perc_neg"])] = [neg_res]
                    first = False
            else:
                for perc in self.fits[g]["Percentiles"].split(";"):
                    p_new, pos_res, neg_res = perc.split("-")
                    percentile_dict["_".join([p_new, "perc_pos"])] = percentile_dict["_".join([p_new, "perc_pos"])] + [pos_res]
                    percentile_dict["_".join([p_new, "perc_neg"])] = percentile_dict["_".join([p_new, "perc_neg"])] + [neg_res]
        return percentile_dict
            
    def weight_par_extract(self, par_list):
        '''This function will take the par list and separate it, returning new lists of each of the indexes'''
        list1 = []
        list2 = []
        list3 = []
        list4 = []
        if len(par_list[0]) == 4:
            for line in par_list:
                list1 = list1 + [line[0]]
                list2 = list2 + [line[1]]
                list3 = list3 + [line[2]]
                list4 = list4 + [line[3]]
            return list1, list2, list3, list4
        elif len(par_list[0]) == 2:
            for line in par_list:
                list1 = list1 + [line[0]]
                list2 = list2 + [line[1]]
            return list1, list2
        else: print("There are not 2/4 indexes")      
        
    def dataframe_creator(self):
        "This function creates a dataframe where each gene is a column with labeled rows of all fit information"
        '''
            Columns: 
        Annotations: 'chrom', 'start', 'stop', 'strand' 
        Sense Main Priors: 'mL_mean', 'mL_stdev', 'sL_mean','sL_stdev', 'tI_mean', 
                           'tI_stdev', 'mT_mean', 'mT_stdev', 'sT_mean','sT_stdev' 
        Sense Weight Priors:'w_LI_mean', 'w_LI_stdev', 'w_E_mean', 'w_E_stdev',
                            'w_T_mean', 'w_T_stdev', 'w_B_mean', 'w_B_stdev' 
        Antisense Main Priors:'mL_a_mean','mL_a_stdev', 'sL_a_mean', 'sL_a_stdev', 'tI_a_mean', 'tI_a_stdev'
        Antisense Weight Priors:'w_aLI_mean', 'w_aLI_stdev','w_aB_mean', 'w_aB_stdev'
        Log Info: 'pos_cov', 'neg_cov', 'total_cov', 'elbo_lrange', 'elbo_urange', 'fit_time_min' 
        '''
        # Create a dataframe based on annotations
        self.df = pd.DataFrame(data=self.annotations, index=["chrom", "start", "stop", "strand"]).transpose()
        self.df['gene'] = self.df.index
        # Add prior values to the dataframe
        self.df['mL_mean'] = self.mL
        self.df['mL_stdev'] = self.mL_std
        self.df['sL_mean'] = self.sL
        self.df['sL_stdev'] = self.sL_std
        self.df['tI_mean'] = self.tI
        self.df['tI_stdev'] = self.tI_std
        if self.ET is True:
            self.df['mT_mean'] = self.mT
            self.df['mT_stdev'] = self.mT_std
            self.df['sT_mean'] = self.sT
            self.df['sT_stdev'] = self.sT_std
            self.df['w_E_mean'] = self.w_E
            self.df['w_E_stdev'] = self.w_E_std
            self.df['w_T_mean'] = self.w_T
            self.df['w_T_stdev'] = self.w_T_std
        self.df['w_LI_mean'] = self.w_LI
        self.df['w_LI_stdev'] = self.w_LI_std
        self.df['w_B_mean'] = self.w_B
        self.df['w_B_stdev'] = self.w_B_std
        if self.antisense is True:
            self.df['mL_a_mean'] = self.mL_a
            self.df['mL_a_stdev'] = self.mL_a_std
            self.df['sL_a_mean'] = self.sL_a
            self.df['sL_a_stdev'] = self.sL_a_std
            self.df['tI_a_mean'] = self.tI_a
            self.df['tI_a_stdev'] = self.tI_a_std
            self.df['w_aLI_mean'] = self.w_aLI
            self.df['w_aLI_stdev'] = self.w_aLI_std
            self.df['w_aB_mean'] = self.w_aB
            self.df['w_aB_stdev'] = self.w_aB_std
        # get the percentiles if they exist
        if self.percentiles:
            percentile_dict = self.percentile_extract()
            # save each percentile as per_pos and per_neg
            for perc, perc_list in percentile_dict.items():
                self.df[perc] = perc_list
        # Add the + & - strand coverage, and elbow range
        # initiate lists
        pos_cov_list = []
        neg_cov_list = []
        elbo_lrange_list = []
        elbo_urange_list = []
        fit_time_list = []
        # iterate through log items to get values
        for gene, coverage in self.log.items():
            if len(coverage) > 1:
                pos_cov_list = pos_cov_list + [coverage['strand_cov'][0]]
                neg_cov_list = neg_cov_list + [coverage['strand_cov'][1]]
                elbo_lrange_list = elbo_lrange_list + [coverage['elbo_range'][0]]
                elbo_urange_list = elbo_urange_list + [coverage['elbo_range'][1]]
                fit_time_list = fit_time_list + [coverage['fit_time_min']]
        # actually add to dataframe
        self.df['pos_cov'] = pos_cov_list
        self.df['neg_cov'] = neg_cov_list
        self.df['elbo_lrange'] = elbo_lrange_list
        self.df['elbo_urange'] = elbo_urange_list
        self.df['fit_time_min'] = fit_time_list                       
        # Sort the dataframe based on geneid
        self.df.sort_index(inplace=True)
        return self.df
 






 # def dataframe_diff_creator(stringfp1, stringfp2, y_values):
# #     '''This function will take two dataframes, find the difference between desired 
# #     variables, & returns a new dataframe with columns of the differences'''
# #     # Create a new dataframe
# #     test_df = pd.DataFrame((globals()[stringfp1]).genes, index=(globals()[stringfp1]).genes )
# #     # Get the differences & add to data frame
# #     for y_value in y_values: 
# # #         if 'w' in y_value:
            
# #         test_df[y_value] = (globals()[stringfp1]).df[y_value]-(globals()[stringfp2]).df[y_value]
# #     # Get the averages of the fit_coverage
# #     return test_df
    



# #======================== Plotting (1 run) ================================

                
# def cov_std_plotter(stringfp):
#     '''This function plots the standard deviation of the posterior distribution of each prior 
#        against the positive & negative strand coverages.
#         - Input: string of FitParse instance (stringfp)
#         - Output: labeled graphs as described above'''
#     prior_list = class_namer(stringfp, prior_average2)
#     xlim = max((globals()[stringfp]).pos_cov_list)/2.5
#     # Graph the pos strand coverage vs std
#     for i in prior_list:
#         fig = plt.figure(figsize=(20, 6))
#         ax = fig.add_subplot()
#         ax.scatter((globals()[stringfp]).pos_cov_list, (globals()[i]).smean_list, color='red')
#         plt.ylabel(str("standard deviation of "+ str(i)))
#         plt.xlabel("positive strand coverage (rds/kb)")
#         xlabel = str(i).replace(str(stringfp), '')
#         plt.xlim([0, xlim]) 
#         plt.title(str("Standard deviation of posterior distribution of "+str(xlabel)+
#                       " vs positive strand coverage in "+str(stringfp)))
#     # Graph the neg strand coverage vs std
#     for i in prior_list:
#         fig = plt.figure(figsize=(20, 6))
#         ax = fig.add_subplot()
#         ax.scatter((globals()[stringfp]).neg_cov_list, (globals()[i]).smean_list, color='blue')
#         plt.ylabel(str("standard deviation of "+ str(i)))
#         plt.xlabel("negative strand coverage (rds/kb)")
#         xlabel = str(i).replace(str(stringfp), '')
#         plt.title(str("Standard deviation of posterior distribution of "+str(xlabel)+
#                       " vs negative strand coverage in "+str(stringfp)))

      
    
# def len_std_plotter(stringfp):
#     prior_list = class_namer(stringfp, prior_average2)
#     # Graph the length vs std
#     for i in prior_list:
#         fig = plt.figure(figsize=(20, 6))
#         ax = fig.add_subplot()
#         ax.scatter((globals()[stringfp]).len_list, (globals()[i]).smean_list, color='green')
#         plt.ylabel(str("standard deviation of "+ str(i)))
#         plt.xlabel("lengths of genes (bp)")
#         xlabel = str(i).replace(str(stringfp), '')
#         plt.title(str("Standard deviation of posterior distribution of "+str(xlabel)+
#                       " vs gene lengths in "+str(stringfp)))

        
# #======================== Plotting (more than 1 run) ================================
# def dataframe_plotter(stringfp1, stringfp2, y_values, corr=True):
#     '''This function plots the y_values (list) between each of the runs (stringfp1 & stringfp2). It also shows the 
#     LINEAR correlation (Pearson's R^2) between the 2 if corr=True'''
#     # only compare where have same genes
#     # ensure the same length as proxy for same genes
#     assert (len((globals()[stringfp1]).df)==len((globals()[stringfp2]).df))
#     # plot y_value of both runs against each other
#     for y_value in y_values: 
#         fig = plt.figure()
#         ax = fig.add_subplot()
#         ax.scatter((globals()[stringfp1]).df[y_value], (globals()[stringfp2]).df[y_value])
#         plt.xlabel(str(y_value + ' of ' + stringfp1))
#         plt.ylabel(str(y_value + ' of ' + stringfp2))
#         plt.title(str(y_value + ' for each gene of ' + stringfp1 + ' and ' + stringfp2))
#         print("hello")
#         if corr==True:
#             print("hello2")
#             # calculate Pearson's correlation
#             corr_value, _ = pearsonr((globals()[stringfp1]).df[y_value], (globals()[stringfp2]).df[y_value])
#             print('Pearsons correlation: %.3f' % corr)
#             plt.title(str(y_value + ' for each gene of ' + stringfp1 + ' and ' + stringfp2 + '(R^2: ' + str(corr_value) + ')'))
#         else:
#             print('false')

        