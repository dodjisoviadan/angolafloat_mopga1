
# import needed libraries
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings # these two lines to remove the annoying warnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)
import itertools
# load for time series
          #cluster_conc = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/category_concentrations_Mai2023.csv')
cluster_conc = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/categoryshort_concentrations_Sept2023.csv')
taxa_list = cluster_conc.category.unique()
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
cluster_conc_depth.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/zoo_heatmap_Nov2023.csv', index=False)
