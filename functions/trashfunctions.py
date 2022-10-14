##======================================== plotting
## Trying to plot multiple bar graph info onto one plot
#Actually plot the histograms of mT
## individual plots
x_df = hct116_full_df.index
for name in name_list:
    y_df_name = name+'_'+'mT_mean'
    title = 'mT_mean of genes for '+name
    ax = hct116_full_df.plot.bar(x='gene', y=y_df_name, rot=90, figsize=(20,6))
    ax.set(title=title)
## same plot but overlapping and no label
fig = plt.figure(figsize=(20,9))
ax = fig.add_subplot(111)
ax.bar(x=hct116_full_df['gene'], height=hct116_full_df['hct116 S81_mT_mean'], width=1)
ax.bar(x=hct116_full_df['gene'], height=hct116_full_df['hct116 S82_mT_mean'], width=0.8)
ax.bar(x=hct116_full_df['gene'], height=hct116_full_df['hct116 S83_mT_mean'], width=0.6)
ax.bar(x=hct116_full_df['gene'], height=hct116_full_df['hct116 S84_mT_mean'], width=0.4)
## sample plot but overlapping BUT label
# ax.bar(x=ind, height=b,  align='center')
# colors = ['steelblue', 'firebrick', 'darksage', 'goldenrod', 'gray']
# for name, color in zip(name_list, colors):
#     y_df_name = name+'_'+'mT_mean'
#     kwargs = plot_kwargs
#     kwargs.update({'color': color, 'label': label})
#     plt.bar(df['gene'], df[y_df_name], alpha=alpha if i else 1, **kwargs)

plt.xticks(rotation=90)

#plt.xticks(rotation=90)

plt.show()


##======================================== weighted average with a diff method
# def calc_weighted_average2(df_list, factor, prior_dist='exp'):
#     '''This function calculates the weighted average of the parameter of interest
#     via its standard deviation AND probability that parameter is that value based on     prior distirbution.
#     Parameters:
#     - df_list: list of dataframes
#     - factor: "prior" of interest (must have _mean & _stdev in df)
#     - average: (OPTIONAL) list of columns you waant to get average of
#     Returns:
#     - df with the following columns:
#         - mean/weight sum (sum(mean/stdev))
#         - weight sum (sum(stdev))
#         - weighted average ((mean/weight sum)/(weight sum))
#         - averaages of values in averagee column
#     '''
#     # create a new df for calcs
#     df2= df_list[0]['gene']
#     df2 = pd.DataFrame(df2)
#     # iterate through all dataframes in list
#     count = 0
#     for df in df_list:
#         # get the prob of mean existing based on prior dist
#         # save values in list
#         exp_list = []
#         for par in df[str(factor + '_mean')]:
#             # find the chance that in exponential distribution
#             prob = 0.0
#             # use tau from config file
#             tau=10000
#             exp_lam = 1/tau
#             # 0 offset for mT & use cdf --> wouldn't that be getting =<par?
#             prob_log = pm.Exponential(var_name='mT', lam = exp_lam).logcdf(par)
#             # to get actual prob, unlog?
#             prob = exp(prob_log)
#             # add chance to list
#             exp_list.append(float(prob))
#         # make the list a dataframe column
#         df2['exp_prob'] = exp_list
#         # if starting out
        
#         if count == 0:
#             df2['mean/weight sum'] = (df[str(factor + '_mean')]/df[str(factor[0] + '_stdev')]*df['exp_prob'])
#             df2['weight sum'] = df[str(factor + '_stdev')] + df['exp_prob'] 
#         else:
#             # need to make sure new one only contains genes that df2 has
#             df = df[df.index.isin(df2.index)]
#             df2 = df2[df2.index.isin(df.index)]
#             # get the proper sums
#             df2['mean/weight sum'] += (df[str(factor + '_mean')]/df[str(factor[0] + '_stdev')]*df['exp_prob'])
#             df2['weight sum'] += df[str(factor + '_stdev')] + df['exp_prob'] 
#         count += 1
#     # normalize mean/weight sum (* inverse weight sum)
#     df2['weighted average'] = df2['mean/weight sum']/((df2['weight sum']))
#     return df2
