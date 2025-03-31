
#Script for plot time series based on K-means clustering of detritus (Emilia's method)

        # Setup
# import needed libraries
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings # these two lines to remove the annoying warnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)
import itertools
import matplotlib.cm as cm
from matplotlib.patches import Patch
from matplotlib import patches

# import csv file
# import TEMP SAL POC
ecopartf = pd.read_table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/Ecopart_diagnostics_data_647.tsv", sep='\t', encoding='unicode_escape')
ecopartf = ecopartf[ecopartf['Profile'].str.contains('a')]
mask=ecopartf["Chlorophyll-a [mg/m3]"]<0
colonne="Chlorophyll-a [mg/m3]"
ecopartf.loc[mask,colonne]=0
ecopartf['colonne']=ecopartf["bbp POC [mgC/m3]"]/ecopartf["Chlorophyll-a [mg/m3]"]
ecopartf.colonne[ecopartf.colonne.isin([-np.inf,np.inf])]=np.nan
ecopartf = ecopartf.iloc[:,[1, 4, 5, 7, 11, 12, 17, 18, 37,38,39,40,42,41,44,36]] # data # date,PRESSION,flux,fluxmip,fluxmap,mip,map,density,tempera,salinity,Doxy,bbppoc,Chla,bbp/CHLA,depth
ecopartf.rename(columns = {'colonne':'bbpPOC/Chla'}, inplace = True)
ecopartf['DATE'] = ecopartf['Date_Time']
ecopartf['DEPTH'] = ecopartf['Depth [m]']
cluster=ecopartf# = ecopartf.iloc[:,[7, 9, 10]] # data
cluster['depth_min'] = cluster['DEPTH']
# create a column layer with depth levels
cluster['layer'] = np.where(cluster.depth_min < 100, "0-100m",
                      np.where(np.logical_and(cluster.depth_min >= 100, cluster.depth_min < 200), "100-200m",
                      np.where(np.logical_and(cluster.depth_min >= 200, cluster.depth_min < 300), "200-300m",
                      np.where(np.logical_and(cluster.depth_min >= 300, cluster.depth_min < 400), "300-400m",
                      np.where(np.logical_and(cluster.depth_min >= 400, cluster.depth_min < 500), "400-500m",
                      np.where(np.logical_and(cluster.depth_min >= 500, cluster.depth_min < 600), "500-600m",
                      np.where(np.logical_and(cluster.depth_min >= 600, cluster.depth_min < 800), "600-800m",
                      np.where(np.logical_and(cluster.depth_min >= 800, cluster.depth_min < 1000), "800-1000m",
                      np.where(cluster.depth_min > 1000, "1000-2000m","else")))))))))

cluster = cluster.drop('depth_min', axis = 1)
cluster = cluster.drop('Date_Time', axis = 1)
cluster = cluster.drop('Pressure [dbar]', axis = 1)
cluster = cluster.drop('Depth [m]', axis = 1)
cluster_depth = cluster.groupby(['Profile','DATE', 'DEPTH','layer']).mean().reset_index()

cluster_depth = cluster_depth[cluster_depth['layer'] != '1000-2000m']

# group the dataframe by 'layer'
grouped = cluster_depth.groupby('layer')
# extract date in 'YYYY-MM-DD' format from datetime column
cluster_conc=cluster_depth
cluster_conc = cluster_conc[cluster_conc['DEPTH'] <= 500]
#cluster_conc = cluster_conc.drop('DEPTH', axis = 1)
# add physical cluster to it
cluster_conc.rename(columns = {'Profile':'profile'}, inplace = True)
# convert 'date' column to datetime format
cluster_conc_depth = cluster_conc
cluster_conc_depth = cluster_conc_depth.groupby(['profile','DATE', 'layer']).agg('mean').reset_index()
#cluster_conc_depth.columns = ['{}_{}'.format(col[0], col[1]) for col in cluster_conc_depth.columns]
#cluster_conc_depth.columns = cluster_conc_depth.columns.str.replace('/', '')
#cluster_conc_depth = cluster_conc_depth.rename(columns={'DATE_': 'DATE', 'layer_' : 'layer'})
0/0
cluster_conc_depth.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/MDS_Environment' + '.csv',index=True)
