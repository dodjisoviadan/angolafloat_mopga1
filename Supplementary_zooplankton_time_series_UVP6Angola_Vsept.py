
#Script for plot time series based on K-means clustering of detritus (Emilia's method)

        # Setup
# import needed libraries
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings # these two lines to remove the annoying warnings from pandas
from matplotlib.backends.backend_pdf import PdfPages
warnings.simplefilter(action='ignore', category=FutureWarning)
import itertools
import matplotlib.cm as cm
from matplotlib.patches import Patch
from matplotlib import patches
tp_lista = [[20210720, 20210815],[20210819, 20210920], [20211015, 20211110], [20211120, 20211215], [20220130, 20220225],[20220315, 20220410]]  # 0

tp_listA = ['2021-08-03','2021-09-04','2021-10-24','2021-12-02','2022-02-11','2022-03-29'] # 0
#tp_listA_char = pd.to_datetime(tp_listA, format="%Y%m%d", errors='coerce')  #

# plot time series
cluster_conc =pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/angola_indiv_angola_rough.csv')
cluster_conc = cluster_conc.rename(columns={'angola_rough': 'category'})
# and the sampled volume
volumes = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/volumes_angola.csv')
indiv_binned =cluster_conc
# Compute abundance and biovolume per taxa and bin for either the rough or medium taxonomic definiton
indiv_binned = indiv_binned.groupby(['category', 'depth_bin', 'profile'], group_keys=False).agg(n = pd.NamedAgg(column = 'category', aggfunc = 'count'),
                                                                                         vol_sph = pd.NamedAgg(column = "vol_sph", aggfunc = 'sum'),
                                                                                         perim = pd.NamedAgg(column = "perim.", aggfunc = 'mean'),
                                                                                         circ = pd.NamedAgg(column = "circ.", aggfunc = 'mean'),
                                                                                         mean = pd.NamedAgg(column = "mean", aggfunc = 'mean'),
                                                                                         kurt = pd.NamedAgg(column = "kurt", aggfunc = 'mean'),
                                                                                         esd = pd.NamedAgg(column = "esd", aggfunc = 'mean'),
                                                                                         fractal = pd.NamedAgg(column = "fractal", aggfunc = 'mean'))

indiv_binned.reset_index(inplace = True) # to keep a column with exports_groups and depth_bin values
indiv_binned.columns = indiv_binned.columns.str.replace('_bin', '')

# add the volume and compute the concentrations
obs = pd.merge(indiv_binned, volumes, how="left", on=['profile', 'depth'])
obs["watervolume"] = obs["volume_L"] / 1000          # volume in m3
obs["conc"] = obs["n"] / obs["watervolume"]          # concentration in n/m3
#obs["vol_sph"] = obs["vol_sph"] / obs["watervolume"] # biovolume concentration in mm3/m3

# keep only the 5m depth bins which have a watervolume <= 1000L
obs = obs[obs["watervolume"] <= 1]
obs = obs.dropna(subset=['watervolume'])

volumes.head()

# compute all taxon x bins combinations
depth_bins = np.unique(volumes.depth)

# and all the profiles
profile_list = np.unique(volumes.profile)

category_list = indiv_binned.category.unique()

# regroup them with the list of taxa
list_ = [profile_list, depth_bins, category_list]
# compute the unique combinations of the 3 lists
combination  = [p for p in itertools.product(*list_)]

# Convert it into a dataframe
column_names = ["profile", "depth", "category"]
all = pd.DataFrame(combination, columns =["profile", "depth", "category"])

# Add the data from obs
full = pd.merge(all, obs, how="left", on=['profile', 'depth', "category"])
full.head()

full = pd.DataFrame(full)

# remove the volume_L column. We'll keep only watervolume (in m3)
full.drop('volume_L', inplace = True, axis = 1)

# consider the observations we don't have as zeroes
cols_to_fill = ['n', 'conc', 'vol_sph', 'watervolume']
full[cols_to_fill]=full[cols_to_fill].fillna(0)

# save it
full['category'].replace('other<living', 'otherliving', inplace=True)
full.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/category_concentrations_Mai2023.csv', index=False)

category_newconc = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/category_concentrations_Mai2023.csv')
category_newconc['category'].replace('t001', 'Acantharia', inplace=True)
category_newconc['category'].replace('t004', 'Acantharia', inplace=True)
category_newconc['category'].replace('t019', 'Collodaria', inplace=True)
category_newconc['category'].replace('t018', 'Other_Rhizaria', inplace=True)


category_newconc['category'].replace('Acantharia', 'Other_Rhizaria', inplace=True)
category_newconc['category'].replace('Foraminifera ', 'Other_Rhizaria', inplace=True)
category_newconc['category'].replace('Phaeodaria', 'Other_Rhizaria', inplace=True)
category_newconc['category'].replace('Copepoda', 'Crustacea', inplace=True)
category_newconc['category'].replace('Other_Crustacea', 'Other_Rhizaria', inplace=True)

category_newconc['category'].replace('Actinopterygii', 'Other_living', inplace=True)
category_newconc['category'].replace('Annelida', 'Other_living', inplace=True)
category_newconc['category'].replace('Appendicularia', 'Other_living', inplace=True)
category_newconc['category'].replace('Chaetognatha', 'Other_living', inplace=True)
category_newconc['category'].replace('Cnidaria', 'Other_living', inplace=True)
category_newconc['category'].replace('Ctenophora', 'Other_living', inplace=True)
category_newconc['category'].replace('Echinodermata', 'Other_living', inplace=True)
category_newconc['category'].replace('Mollusca', 'Other_living', inplace=True)
category_newconc['category'].replace('Thaliacea', 'Other_living', inplace=True)
category_newconc['category'].replace('otherliving', 'Other_living', inplace=True)

category_newconc.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/categoryshort_concentrations_Sept2023.csv', index=False)

# load for time series
          #cluster_conc = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/category_concentrations_Mai2023.csv')
cluster_conc = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/categoryshort_concentrations_Sept2023.csv')
taxa_list = category_newconc.category.unique()
taxashort_list=taxa_list[[3,5,0,4,1]]
# add date
profiles_date = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/list_of_profiles_angola.csv', sep=';', encoding='unicode_escape')
profiles_date.rename(columns={'date': 'date_ancien'}, inplace=True)
profiles_date.rename(columns={'dateok': 'date'}, inplace=True)
profiles_date.rename(columns = {'Profile':'profile'}, inplace = True)
# select only ascent profiles :
profiles_date = profiles_date[profiles_date['profile'].str.contains('a')]
profiles_date['CYCLE_NUMBER'] = (range(1,119)) # profiles_date['CYCLE_NUMBER'] = list(range(1,119))
cluster_conc = pd.merge(cluster_conc, profiles_date, on = 'profile')
list
# convert it into date format
# convert 'date' column to datetime format
AAA=round(cluster_conc['date'],0)
cluster_conc['date'] = pd.to_datetime(AAA, format="%Y%m%d", errors='coerce')#pd.to_datetime(cluster_conc['date'])

# extract date in 'YYYY-MM-DD' format from datetime column
cluster_conc['date'] = cluster_conc['date'].dt.strftime('%Y-%m-%d')

# add a column with depth levels
cluster_conc['layer'] = np.where(cluster_conc.depth < 100, "0-100m",
                      np.where(np.logical_and(cluster_conc.depth >= 100, cluster_conc.depth < 300), "100-300m",
                      np.where(np.logical_and(cluster_conc.depth >= 300, cluster_conc.depth < 600), "300-600m",
                      np.where(np.logical_and(cluster_conc.depth >= 600, cluster_conc.depth < 1000), "600-1000m",
                      np.where(cluster_conc.depth > 1000, "1000-2000m",
                      "else")))))

cluster_conc = cluster_conc[cluster_conc['depth'] <= 1000]
# add physical cluster to it

cluster_conc = cluster_conc.drop("CYCLE_NUMBER", axis = 1)
cluster_conc['date'] = pd.to_datetime(cluster_conc['date'])
# group by cluster, date and layer -> sum, mean and std
cluster_conc_depth = cluster_conc
cluster_conc_depth = cluster_conc_depth.groupby(['category', 'date', 'layer']).agg({'conc': ['mean', 'std', 'sum'], 'vol_sph': ['mean', 'std', 'sum']}).reset_index()
# flatten the MultiIndex column names to a single level
cluster_conc_depth.columns = ['{}_{}'.format(col[0], col[1]) for col in cluster_conc_depth.columns]
cluster_conc_depth.columns = cluster_conc_depth.columns.str.replace('/', '')
cluster_conc_depth = cluster_conc_depth.rename(columns={'category_': 'category', 'date_': 'date', 'layer_' : 'layer'})

cluster_conc_depth['date'] = pd.to_datetime(cluster_conc_depth['date'])

# define cluster of interest

t = "Copepoda"
#t = "t004"
#t = "t001"
t = "Acantharia"
#t = "t019"
t = "Collodaria"
#t = "t018"
t = "Other_Rhizaria"

# Plot the sum concentration of cluster of interest in each depth layer
# Define the plot
#fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3)
'''fig = plt.figure()
gs = fig.add_gridspec(2, 3, hspace=0, wspace=0)
(ax1,ax2,ax3),(ax4,ax5,ax6) = gs.subplots(sharex='col', sharey='row')
fig = plt.figure(figsize=(5, 5))
gs = fig.add_gridspec(5, hspace=0)
(ax1, ax2, ax3, ax4, ax5) = gs.subplots(sharex=True, sharey=True)
#fig.suptitle('Sharing both axes')'''

path_to_figures = Path(
    'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/category_or_taxa time_series_sumVsept').expanduser()

count_t=1

pp = PdfPages(str(path_to_figures) +'/Time series Zooplankton' + '.pdf')
for t in taxashort_list:
    #fig, ax = plt.subplots()
    fig = plt.figure(figsize=(5, 3.5))  # (19.20, 9.83)
    plt.rcParams.update({'font.size': 8})
    ax = plt.subplot(1, 1, 1)
    #ax= eval('ax'+ str(count_t))

    # Loop over the depth layers
    plt.yscale("log")  # np.log10+1
    for layer in ['0-100m']: #, '100-300m', '300-600m', '600-1000m'
        layer_data = cluster_conc_depth.loc[cluster_conc_depth['layer'] == layer]
        if len(layer_data) > 0:
            # Filter the data to only include the current depth layer and taxon
            taxon_data = layer_data.loc[layer_data['category'] == t]
            # Plot the concentration vs. date for the current taxon and depth layer

            ax.plot(taxon_data['date'], (taxon_data['conc_sum']), marker='o', linestyle='--', label=layer, markersize=3,
                    markeredgecolor='white', markeredgewidth=0.25)
            #0/0
            high = np.max(np.log10(taxon_data['conc_sum']+1))
            tt=taxon_data['date']
            #tp_listA   pour le premier groupe (31),
            # ax.plot(taxon_data['date'][taxon_data.index[[30,40,55,67,87,102]]], [high,high,high,high,high,high] , marker='o',         
            #linestyle='None',label='bloom', markersize=8, markeredgecolor='green', markeredgewidth=0.5)
            # add shadded area arround the bloom
            high=10**(5.5)
            high1=10**(5.5-0.3) #high-0.3
            '''ax.fill_between(taxon_data['date'][taxon_data.index[[27,34]]], [1,1],[high,high], alpha=0.2, color='gray')
            ax.fill_between(taxon_data['date'][taxon_data.index[[36, 46]]], [1, 1], [high, high], alpha=0.2, color='gray')
            ax.fill_between(taxon_data['date'][taxon_data.index[[53, 61]]], [1, 1], [high, high], alpha=0.2, color='gray')
            ax.fill_between(taxon_data['date'][taxon_data.index[[64, 72]]], [1, 1], [high, high], alpha=0.2, color='gray')
            ax.fill_between(taxon_data['date'][taxon_data.index[[84, 92]]], [1, 1], [high, high], alpha=0.2, color='gray')
            ax.fill_between(taxon_data['date'][taxon_data.index[[98, 108]]], [1, 1], [high, high], alpha=0.2, color='gray')
            '''
            ax.fill_between(taxon_data['date'][taxon_data.index[[25, 35]]], [1, 1], [high, high], alpha=0.2,color='gray')
            ax.fill_between(taxon_data['date'][taxon_data.index[[36, 46]]], [1, 1], [high, high], alpha=0.2,color='gray')
            ax.fill_between(taxon_data['date'][taxon_data.index[[50, 61]]], [1, 1], [high, high], alpha=0.2,color='gray')
            ax.fill_between(taxon_data['date'][taxon_data.index[[64, 75]]], [1, 1], [high, high], alpha=0.2,color='gray')
            ax.fill_between(taxon_data['date'][taxon_data.index[[83, 95]]], [1, 1], [high, high], alpha=0.2,color='gray')
            ax.fill_between(taxon_data['date'][taxon_data.index[[96, 108]]], [1, 1], [high, high], alpha=0.2,color='gray')

            ax.plot([], marker="s", markersize=11, linestyle="", color='gray', label="Bloom Column")

            plt.text(taxon_data['date'][taxon_data.index[30]], high1, '1', fontsize=10, color='gray')
            ax.text(taxon_data['date'][taxon_data.index[41]], high1, '2', fontsize=10, color='gray')
            ax.text(taxon_data['date'][taxon_data.index[55]], high1, '3', fontsize=10, color='gray')
            ax.text(taxon_data['date'][taxon_data.index[68]], high1, '4', fontsize=10, color='gray')
            ax.text(taxon_data['date'][taxon_data.index[89]], high1, '5', fontsize=10, color='gray')
            ax.text(taxon_data['date'][taxon_data.index[102]], high1, '6', fontsize=10, color='gray')

            # Add the layer label to the legend
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='lower right')
            plt.ylim([1, 10**5.5])
            #plt.tight_layout()
    count_t +=1
        # Set the y axis label for the first depth layer  'Log10 (+ 1)
    ax.set_ylabel('Integrated concentration [#.$m^{-2}$]')
    ax.set_xlabel('Date')
    #ax.set_yticks()

    plt.title(t + '')
    # Show the plot
    plt.show()
    
    '''for ax in fig.get_axes():
        ax.label_outer()'''

    path_to_figures = Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/category_or_taxa time_series_sumVsept').expanduser()
    plt.savefig(str(path_to_figures) + '/FEV2024category_sum_concentration' + t + 'okVFone.png', dpi=1200)
    pp.savefig(fig, dpi=1200)  # Save each figure in the pdf #plt.gcf()
    # close the pdf document
    pp.close()
    0/0

    # Plot the mean concentration of taxon of interest in each depth layer with two different axis
    # (because too much differences between layers) and confidence interval (std)
    # Define the plot
    fig, ax1 = plt.subplots()

    # Loop over the depth layers

    for layer in ['0-100m']: #, '100-300m', '300-600m', '600-1000m'
        layer_data = cluster_conc_depth.loc[cluster_conc_depth['layer'] == layer]
        if len(layer_data) > 0:
            # Filter the data to only include the current depth layer and taxon
            taxon_data = layer_data.loc[layer_data['category'] == t]
            # Plot the concentration vs. date for the current taxon and depth layer
            plt.yscale("log")
            ax1.plot(taxon_data['date'], taxon_data['conc_mean'], marker='o', linestyle='--', label=layer,
                     markersize=3.5,
                     markeredgecolor='white', markeredgewidth=0.5)
            # Add the layer label to the legend
            handles, labels = ax1.get_legend_handles_labels()
            ax1.legend(handles, labels, loc='upper right')
            # Clip the negative values of the error bars to zero
            error_bar_lower = np.clip(taxon_data['conc_mean'] - taxon_data['conc_std'], 0, None)
            # Add confidence intervals to the plot with clipped error bars
            ax1.fill_between(taxon_data['date'], error_bar_lower,
                             taxon_data['conc_mean'] + taxon_data['conc_std'], alpha=0.2)

            # Set the y axis label for the first depth layer
    ax1.set_ylabel('Mean concentration (ind/m3)')

    # Set the y axis label for the other depth layers
    plt.title(t + ' mean concentration (ind/m3)')
    # Show the plot
    plt.show()

    path_to_figures = Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/category_or_taxa time_series_meanVsept').expanduser()
    plt.savefig(str(path_to_figures) + '/category_mean_concentration' + t + '.png', dpi=1200)

plt.close('all')


