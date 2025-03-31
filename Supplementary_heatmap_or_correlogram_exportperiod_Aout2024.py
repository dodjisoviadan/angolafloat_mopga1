from gettext import install

import matplotlib
#################################################################################
# heat map with all environmental variables and lagrandian diagnostics with
# particles concentrations #
#################################################################################

# libraries

import pandas as pd
from pathlib import Path
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'
pd.options.display.max_columns = 100
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import warnings # these two lines to remove the annoying warnings from pandas
from matplotlib.backends.backend_pdf import PdfPages
warnings.simplefilter(action='ignore', category=FutureWarning)

from datetime import date

# Import lagrangian and eulerian data
lagrangian_eulerian =pd.read_table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/final_data_an35an.csv", sep=',')
# import environmental data
bgc_angola_clean0 = pd.read_csv("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_angola_clean.csv")
# import particles cluster concentration (in this case only exclusif members)
#cluster_conc = pd.read_csv(r"C:\Users\lemoi\Desktop\GIT\k-means\clusters_concentrations_q25.csv")
cluster_conc0 = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/clusters_concentrations.csv')
zoo_heatmap = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/zoo_heatmap_NOV2023.csv')
profiles_date = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/list_of_profiles_angola.csv', sep=';', encoding='unicode_escape')
profiles_date['CYCLE_NUMBER'] = (range(1,119))
prof=profiles_date[['profile','dateok','CYCLE_NUMBER']]
prof.rename(columns={'profile': 'Profile'}, inplace=True)
prof['date'] = [str(int2date(d)) for d in prof['dateok']]
zoo_conc = pd.merge(zoo_heatmap, prof, on = 'date')
# clycle num and date
cycleIDate = pd.read_csv('C:\MOPGA 2022-2023 Dodji\ARGO float Angola\To be send to Alberto/albertoeddyokTable.csv')
cycleIDate['CYCLE_NUMBER'] = list(range(1, 119))
# import MiP, MaP abundance and flux data
ecopart = pd.read_table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/Ecopart_diagnostics_data_647.tsv", sep='\t', encoding='unicode_escape')
ecopart = ecopart[ecopart['Profile'].str.contains('a')]

ecopart_profiles = ecopart["Profile"].unique()
ecopart_profiles = pd.DataFrame(ecopart_profiles, columns=['Profile'])
ecopart_profiles = ecopart_profiles.sort_values('Profile')
ecopart_profiles['CYCLE_NUMBER'] = list(range(1, 119))
bbpdecomp0 = pd.read_table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bbp_decomposed.csv", sep=',', encoding='unicode_escape')
layers=[[0,100],[100,300],[300,600],[600,1000],[100,1000]]
layerstr=['0-100m','100-300m','300-600m','600-1000m','100-1000m']
tp_lista = [[20210505, 20220426],[20210819, 20210925],[20211001, 20211113],[20220120, 20220306]] # 0
count_tp_list=1
tp_lista = [[20210505, 20220426]]
for tp_list in tp_lista:
    tp_list_char = pd.to_datetime(tp_list, format="%Y%m%d", errors='coerce')  # pd.to_datetime(cluster_conc['date'])
    # select only the EXPORT period
    ecopart['Date_Time'] = pd.to_datetime(ecopart['Date_Time'])
    date_mask1 = (cycleIDate['date'] >= tp_list[0]) & (cycleIDate['date'] <= tp_list[1])
    cyclo = cycleIDate[date_mask1]
    period = cyclo['CYCLE_NUMBER']
    date_mask = (ecopart['Date_Time'] >= tp_list_char[0]) & (ecopart['Date_Time'] <= tp_list_char[1])
    ecopartmk = ecopart[date_mask]

    path_to = Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/Heat_Map/Period' + str(
        count_tp_list)).expanduser()
    pp = PdfPages(str(path_to) + '/Heatmaps' + '.pdf')

    count_layer=0
    for layer in layers[0:-4]:
        # _100 is not really surface but the layer
        # try other layer depth
        ecopart_100 = ecopartmk
        depth_mask = (ecopart_100['Depth [m]'] >= layer[0]) & (ecopart_100['Depth [m]'] <= layer[1])
        ecopart_100 = ecopart_100[depth_mask]
        ecopart_100 = ecopart_100[['Profile', 'Flux_mgC_m2', 'MiP_abun', 'MaP_abun']]
        ecopart_100 = ecopart_100.groupby(['Profile']).agg({'MaP_abun': 'mean',
                                                            'MiP_abun': 'mean',
                                                            'Flux_mgC_m2': 'mean', }).reset_index()

        # select divergence data
        divergence = lagrangian_eulerian[['Profile',
                                          'Eul_div_GlobEkmanDt_000daysBackward_mean',
                                          'Eul_div_GlobEkmanDt_005daysBackward_mean',
                                          'Eul_div_GlobEkmanDt_010daysBackward_mean',
                                          'Eul_div_GlobEkmanDt_015daysBackward_mean',
                                          'Eul_div_GlobEkmanDt_020daysBackward_mean',
                                          'Eul_div_GlobEkmanDt_030daysBackward_mean',
                                          'Eul_div_GlobEkmanDt_045daysBackward_mean']]

        # select vorticity data
        vorticity = lagrangian_eulerian[['Profile',
                                         'Eul_vort_GlobEkmanDt_000daysBackward_mean',
                                         'Eul_vort_GlobEkmanDt_005daysBackward_mean',
                                         'Eul_vort_GlobEkmanDt_010daysBackward_mean',
                                         'Eul_vort_GlobEkmanDt_015daysBackward_mean',
                                         'Eul_vort_GlobEkmanDt_020daysBackward_mean',
                                         'Eul_vort_GlobEkmanDt_030daysBackward_mean',
                                         'Eul_vort_GlobEkmanDt_045daysBackward_mean']]

        # select front data
        front = lagrangian_eulerian[['Profile',
                                     'Ftle_GlobEkmanDt_005daysBackward_mean_delta0ftle010',
                                     'Ftle_GlobEkmanDt_010daysBackward_mean_delta0ftle010',
                                     'Ftle_GlobEkmanDt_015daysBackward_mean_delta0ftle010',
                                     'Ftle_GlobEkmanDt_020daysBackward_mean_delta0ftle010',
                                     'Ftle_GlobEkmanDt_030daysBackward_mean_delta0ftle010',
                                     'Ftle_GlobEkmanDt_045daysBackward_mean_delta0ftle010']]

        # select only lagrangian chlorophyll data
        lagrangian_chla = lagrangian_eulerian[['Profile',
                                               'Lagr_Chl_GlobEkmanDt_Oceancolour_000daysBackward_mean',
                                               'Lagr_Chl_GlobEkmanDt_Oceancolour_005daysBackward_mean',
                                               'Lagr_Chl_GlobEkmanDt_Oceancolour_010daysBackward_mean',
                                               'Lagr_Chl_GlobEkmanDt_Oceancolour_015daysBackward_mean',
                                               'Lagr_Chl_GlobEkmanDt_Oceancolour_020daysBackward_mean',
                                               'Lagr_Chl_GlobEkmanDt_Oceancolour_025daysBackward_mean',
                                               'Lagr_Chl_GlobEkmanDt_Oceancolour_030daysBackward_mean',
                                               'Lagr_Chl_GlobEkmanDt_Oceancolour_035daysBackward_mean',
                                               'Lagr_Chl_GlobEkmanDt_Oceancolour_040daysBackward_mean',
                                               'Lagr_Chl_GlobEkmanDt_Oceancolour_045daysBackward_mean']]

        # select only the EXPORT period
        bgc_angola_clean = pd.merge(bgc_angola_clean0, ecopart_profiles, on='CYCLE_NUMBER')

        date_mask = (bgc_angola_clean['CYCLE_NUMBER'] >= min(period)) & (
                    bgc_angola_clean['CYCLE_NUMBER'] <= max(period))
        bgc_angola_cleanmk = bgc_angola_clean[date_mask]

        bbpdecomp = pd.merge(bbpdecomp0, ecopart_profiles, on='CYCLE_NUMBER')
        date_mask2 = (bbpdecomp['CYCLE_NUMBER'] >= min(period)) & (
                bbpdecomp['CYCLE_NUMBER'] <= max(period))
        bbpdecompmk = bbpdecomp[date_mask2]

        # first try surface layer
        # bgc_100 = bgc_angola_clean[bgc_angola_clean['binned_depth'] < 100]

        # try other layer depth
        bgc_100 = bgc_angola_cleanmk
        depth_mask = (bgc_100['binned_depth'] >= layer[0]) & (bgc_100['binned_depth'] <= layer[1])
        bgc_100 = bgc_100[depth_mask]
        bgc_100 = bgc_100.groupby(['CYCLE_NUMBER']).agg({'TEMP': 'mean',
                                                         'PSAL': 'mean',
                                                         'CHLA_ADJUSTED': 'mean',
                                                         'DOXY_ADJUSTED': 'mean',
                                                         'BBP700': 'mean'}).reset_index()

        bgc_100 = pd.merge(bgc_100, ecopart_profiles, on='CYCLE_NUMBER')
        bgc_100 = bgc_100.drop('CYCLE_NUMBER', axis=1)

        bbp_100 = bbpdecompmk
        depth_maski = (bbp_100['PRES'] >= layer[0]) & (bbp_100['PRES'] <= layer[1])
        bbp_100 = bbp_100[depth_maski]
        bbp_100 = bbp_100.groupby(['CYCLE_NUMBER']).agg({'bbsr': 'mean',
                                                         'bbl': 'mean'}).reset_index()
        bbp_100 = pd.merge(bbp_100, ecopart_profiles, on='CYCLE_NUMBER')
        bbp_100 = bbp_100.drop('CYCLE_NUMBER', axis=1)

        cluster_conc0.rename(columns={'profile': 'Profile'}, inplace=True)
        cluster_conc = pd.merge(cluster_conc0, ecopart_profiles, on="Profile")
        # select only the EXPORT period
        date_mask = (cluster_conc['CYCLE_NUMBER'] >= min(period)) & (cluster_conc['CYCLE_NUMBER'] <= max(period))
        cluster_concmk = cluster_conc[date_mask]

        # try sufrace layer first
        # cluster_conc = cluster_conc[cluster_conc['depth'] < 100]

        # try other layer depth
        depth_mask = (cluster_concmk['depth'] >= layer[0]) & (cluster_concmk['depth'] <= layer[1])
        cluster_concd = cluster_concmk[depth_mask]
        cluster_concd = cluster_concd[['Profile', 'Cluster', 'conc']]
        cluster_concd = cluster_concd.groupby(['Profile', 'Cluster']).agg({'conc': 'mean'}).reset_index()
        cluster_conc_reshape = pd.pivot_table(cluster_concd, values='conc', index='Profile',
                                              columns='Cluster').reset_index()

        date_mask = (zoo_conc['CYCLE_NUMBER'] >= min(period)) & (zoo_conc['CYCLE_NUMBER'] <= max(period))
        zoo_concmk = zoo_conc[date_mask]
        depth_mask = (zoo_concmk['layer'] == layerstr[count_layer])
        zoo_concd = zoo_concmk[depth_mask]
        zoo_concd = zoo_concd[['Profile', 'category', 'conc_mean']]
        zoo_conc_reshapes = pd.pivot_table(zoo_concd, values='conc_mean', index='Profile',
                                              columns='category').reset_index()
        zoo_conc_reshape=zoo_conc_reshapes[['Profile','Collodaria','Crustacea','Cyanobacteria','Other_Rhizaria','Other_living']]

        cluster_conc_reshape.rename(columns={'cluster 1':'Flakes'}, inplace=True)
        cluster_conc_reshape.rename(columns={'cluster 2': 'Agglomerates'}, inplace=True)
        cluster_conc_reshape.rename(columns={'cluster 3': 'Strings'}, inplace=True)
        cluster_conc_reshape.rename(columns={'cluster 4': 'Spheres'}, inplace=True)

        # merge all dataset in the same one
        heatmap = pd.merge(ecopart_100, lagrangian_chla, on='Profile')
        heatmap = pd.merge(heatmap, bgc_100, on='Profile')
        heatmap = pd.merge(heatmap, bbp_100, on='Profile')
        heatmap = pd.merge(heatmap, front, on='Profile')
        heatmap = pd.merge(heatmap, divergence, on='Profile')
        heatmap = pd.merge(heatmap, vorticity, on='Profile')
        heatmap = pd.merge(heatmap, cluster_conc_reshape, on='Profile')
        heatmap = pd.merge(heatmap, zoo_conc_reshape, on='Profile')
        heatmap.rename(columns={'Chlorophyll-a [mg/m3]': 'chla_100_mean'}, inplace=True)
        heatmap = heatmap.drop('Profile', axis=1)
        aaa=list(heatmap.columns)
        aaa1=heatmap[aaa[26:40]]
        heatmap[aaa[26:40]]=aaa1.abs()



        # Create mask for significant correlations
        # Calculate Spearman correlation and p-values for all columns
        corr, pval = spearmanr(heatmap)

        # Filter correlation matrix for significant correlations only
        corr = pd.DataFrame(corr, columns=heatmap.columns, index=heatmap.columns)
        corr = corr.mask(pval >= 0.05)

        # Filter correlation matrix for specific columns only
        corr = corr[['MaP_abun', 'MiP_abun', 'Flux_mgC_m2', 'Flakes', 'Agglomerates', 'Strings', 'Spheres','Collodaria','Crustacea','Cyanobacteria','Other_Rhizaria','Other_living']]
        corr['index'] = range(1, 50)
        #0/0
        corr = corr[corr['index'] != 1]
        corr = corr[corr['index'] != 2]
        corr = corr[corr['index'] != 3]
        corr = corr[corr['index'] != 4]
        corr = corr[corr['index'] != 5]
        corr = corr[corr['index'] != 7]
        corr = corr[corr['index'] != 9]
        corr = corr[corr['index'] != 10]
        corr = corr[corr['index'] != 11]
        corr = corr[corr['index'] != 12]
        corr = corr[corr['index'] != 13]

        corr = corr[corr['index'] != 19]
        corr = corr[corr['index'] != 20]

        corr = corr[corr['index'] != 21]
        corr = corr[corr['index'] != 23]

        corr = corr[corr['index'] != 25]
        corr = corr[corr['index'] != 26]
        corr = corr[corr['index'] != 27]
        corr = corr[corr['index'] != 28]

        corr = corr[corr['index'] != 30]
        corr = corr[corr['index'] != 32]
        corr = corr[corr['index'] != 33]
        corr = corr[corr['index'] != 34]
        corr = corr[corr['index'] != 35]

        corr = corr[corr['index'] != 37]
        corr = corr[corr['index'] != 39]
        corr = corr[corr['index'] != 40]

        '''corr = corr[corr['index'] != 41]
        corr = corr[corr['index'] != 42]
        corr = corr[corr['index'] != 43]
        corr = corr[corr['index'] != 44]
        corr = corr[corr['index'] != 45]
        corr = corr[corr['index'] != 46]
        corr = corr[corr['index'] != 47]
        corr = corr[corr['index'] != 48]
        corr = corr[corr['index'] != 49]'''

        corr = corr.drop('index', axis=1)

        '''ax1 = plt.subplot(1, 1, 1)'''

        fig, ax = plt.subplots(figsize=(10, 10))
        fig.tight_layout()
        plt.subplot(1, 1, 1)
        plt.clf()
        plt.subplots_adjust(hspace=0.2, wspace=0.2, top=0.88, bottom=0.117, left=0.42, right=1)
        # fig, ax = plt.subplots()
        s1=sns.heatmap(corr,vmin=-1, vmax=1, center=0,cmap='bwr',square=False,annot=True,fmt='.2f') #
        fig.set_tight_layout(True)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45,horizontalalignment='right', size=8)

        '''ax.set_yticklabels(ax.get_yticklabels(),rotation=0,horizontalalignment='right',size = 4)'''

        # ax.xaxis(tickangle = 45,automargin = True)???
        plt.title("Spearman correlation (" + str(layer[0]) + '-' + str(layer[1]) + 'm )' + ';' + "mean:" + str(
            tp_list_char.date[0]) + ' to ' + str(
            tp_list_char.date[1]))
        #plt.get_current_fig_manager().window.state('zoomed')
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.show()
        plt.pause(3)

        0/0

        path_to = Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/Heat_Map/Period'+str(count_tp_list)).expanduser()
        plt.savefig(str(path_to) + '/NEW2024-HEATMAP_ANGOLA_DEC2023 (' + str(layer[0]) + '-' + str(layer[1]) + 'm)' + '-' + str(
            tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.png',dpi=1200)
        pp.savefig(fig, dpi=1200)  # Save each figure in the pdf # plt.gcf()

        plt.close('all')
        count_layer += 1
    pp.close()
    #count_tp_list=+1