
# # # 1. Setup
# import needed libraries
import os
from pathlib import Path

# for data management
import pandas as pd
import warnings # these two lines to remove the annoying warnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.set_option('display.max_rows', 500, 'display.max_columns', 500) # to be able to see a whole table
import numpy as np
import datetime
from datetime import date
from datetime import datetime
from datetime import timedelta
# Settling speed computations
from scipy.optimize import curve_fit
import scipy.stats

# for figures
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors, cm

# for interpolation
from plotting_funcs import contour_levels_func, gridding_func
from matplotlib import ticker

# set up working space
path_to_data = Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/QuartileMiriamfromSpeedpython').expanduser()
os.chdir(str(path_to_data))

tp_list = [20211001, 20211113]  # 2
tp_list = [20210713, 20210819] # 0
tp_list = [20210819, 20210925]  # 1
tp_list = [20211113, 20211215]  # 3
tp_list = [20211113, 20220119]  # 3b
tp_list = [20220120, 20220306]  # 4
tp_list = [20220306, 20220426]  # 5

path_to_SPEED= Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/ClusterMiriamdoc/SPEED').expanduser()

tp_lista = [[20210713, 20210819], [20210819, 20210925],[20211001, 20211113],[20211113, 20211215],[20211113, 20220119],[20220120, 20220306],[20220306, 20220426]] # 0
#tp_lista = [[20211001, 20211113]] # 0
tp_lista = [[20210819, 20210925],[20211001, 20211113],[20220120, 20220306]] # 0
tp_lista = [[20210819, 20210925],[20211001, 20211113],[20211113, 20211215],[20220120, 20220306],[20220310, 20220416]] # 0
ecopart_bins = [0.064, 0.128, 0.256, 0.512, 1.02, 2.05, 4.1, 8.19, 16.4]  # bins
# cluster from kmeans or umap ( check the script : cluster_time_series......)
# the differente method of Miriam save using the below names ( save files) and read from the csv files (csvname)
savefile = ['Miriam_wosizeKmeans4', 'Miriam_sizeUMAPhclust4', 'Miriam_wosizeUMAPKmeans4',
            'Miriam_wosizeUMAPhclust5']
csvname = ['Angola_detritus_forMiriam_mod', 'umapclusters', 'umapclusters_sizerm', 'umapclusters_sizerm']
colomnname_incsv = ['cluster_reduced_wosize_kmeans4', 'clust_h4', 'kmeans_4', 'clust_h5']
Q = range(0, 4)
for q in Q:

    for tp_list in tp_lista:

        tp_list_char = pd.to_datetime(tp_list, format="%Y%m%d", errors='coerce')  # pd.to_datetime(cluster_conc['date'])

        for cl in range(1, 7):
            # cl=2

            path_to_figures = Path(
                'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/ClusterMiriamdoc/SPEED/SPEED' + str(
                    cl)).expanduser()

            arg1 = 0  # For the defitition choice of taxo definition : either 0 for EXPORTS_rough or 1 for EXPORTS_medium
            argb = 0

            # # # 2. Compute settling speed for Ecotaxa images
            if q < 4: # Just to have 3 choice and not modified to much my script
                # first I need to import the csv indiv_angola_rough that contains object id in order to compute the concentration of each cluster considered as species
                indiv0 = pd.read_csv(
                    'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/angola_indiv_angola_rough.csv',
                    low_memory=False, index_col=False)
                indiv0.rename(columns={'angola_rough': 'EXPORTS_rough'}, inplace=True)
                indiv0.rename(columns={'angola_medium': 'EXPORTS_medium'}, inplace=True)
                # cluster kmeans and object id without size
                cluster00 = pd.read_csv(
                    'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/aMiriamAout2023/' + csvname[q] + '.csv', sep=';',
                    encoding='unicode_escape')
                cluster0 = cluster00[['id', colomnname_incsv[q]]]
                cluster0.rename(columns={'id': 'object_id'}, inplace=True)
                cluster0.rename(columns={colomnname_incsv[q]: 'Cluster'}, inplace=True)
                # define a lambda function to add a hyphen to the start of each value
                if q > 0:
                    add_hyphen = lambda x: 'cluster ' + str('%.0f' % (x + 1))
                else:
                    add_hyphen = lambda x: 'cluster ' + str('%.0f' % (x))
                # apply the lambda function to the 'my_column' column using the .apply() function
                cluster0['Cluster'] = cluster0['Cluster'].apply(add_hyphen)
                # CHECK if luster 5 exist
                CL0 = np.unique(cluster0['Cluster'].str.find('cluster 5'))

                # add the column cluster
                indiv0 = pd.merge(indiv0, cluster0, on='object_id')

                if cl == 6:  # all detritus without cluster
                    indiv = indiv0
                else:
                    # for cl in range(1,5): the four OR FIVE clusters
                    indiv = indiv0[indiv0['Cluster'] == 'cluster ' + str(cl)]
                if np.size(indiv) != 0:

                    # Choose the definition either as rough or medium
                    definition_needed = ['EXPORTS_rough', 'EXPORTS_medium']
                    number = arg1
                    def_chosen = definition_needed[number]
                    def_chosen
                    # For each size class, compute the abundance per layer
                    # volumes = pd.read_csv(str(path_to_data) + '/Data_created/2.volumes_EXPORTS_NA_0_500m.csv')
                    volumes = pd.read_csv(
                        'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/volumes_angola.csv')
                    detritus0 = indiv[indiv['EXPORTS_rough'] == 'detritus']

                    count_tp = 0
                    count = 1
                    fig = plt.figure(figsize=(19.20, 9.83))  # (5, 3.5)
                    plt.rcParams.update({'font.size': 12})
                    ax1 = fig.add_axes([0, 0, 0, 0])
                    for tp in tp_list[0:-1]:
                        detritus = detritus0[detritus0['date'] > tp_list[count_tp]].reset_index(drop=True)
                        detritus = detritus[detritus['date'] < tp_list[count_tp + 1]].reset_index(drop=True)
                        detritus = detritus[detritus['depth'] < 500].reset_index(drop=True)

                        detritus.shape


                        # bin the data per 100m
                        def rounditup(x, precision, method):
                            import math
                            if method == "floor":
                                return math.floor(x / precision) * precision
                            elif method == "ceiling":
                                return math.ceil(x / precision) * precision
                            else:
                                return "give the parameter floor or ceiling"


                        volumes.loc[:, 'depth_part'] = 99999
                        volumes['depth_part'] = volumes["depth"].apply(
                            lambda x: rounditup(x, 100, "floor") + 50 if x < 600 else rounditup(x, 200, "floor") + 100)
                        detritus.loc[:, 'depth_part'] = 99999
                        detritus['depth_part'] = detritus["depth"].apply(
                            lambda x: rounditup(x, 100, "floor") + 50 if x < 600 else rounditup(x, 200, "floor") + 100)
                        # Create a function to determine the bin_min and bin_max of each particle based on the ecopart bins
                        # ecopart_bins = [0.645, 1.020,2.050,4.100,6.500,8.190,10.300,16.400,26.000, 10e7]
                        # ecopart_bins = [0.645, 0.813, 1.020,1.290,1.630,2.050,2.580,3.250,4.100,5.160,6.500,8.190
                        #        ,10.300,13.000,16.400,20.600,26.000, 10e7]
                        # Not working for 0.645-0.813µm

                        ecopart_bins = [0.645, 0.813, 1.020, 1.290, 1.630, 2.050, 2.580, 3.250, 4.100, 5.160, 6.500,
                                        8.190
                            , 10.300, 13.000, 16.400, 20.600, 26.000, 10e7]
                        ecopart_binsb = [0.645, 1.020, 2.050, 4.100, 6.500, 8.190, 10.300, 16.400, 26.000, 10e7]


                        def b_min_max(x, lst):
                            bin_max = next(b_max for b_max in lst if x < b_max)
                            bin_min = lst[lst.index(bin_max) - 1]
                            return [bin_min, bin_max]


                        # Add two columns to the table det_1 with bin_max and bin_min
                        if argb == 0:
                            detritus['bin_min'] = detritus['esd'].apply(lambda x: b_min_max(x, ecopart_bins)[0])
                            detritus['bin_max'] = detritus['esd'].apply(lambda x: b_min_max(x, ecopart_bins)[1])
                        else:
                            detritus['bin_min'] = detritus['esd'].apply(lambda x: b_min_max(x, ecopart_binsb)[0])
                            detritus['bin_max'] = detritus['esd'].apply(lambda x: b_min_max(x, ecopart_binsb)[1])

                        # do this only for aggregate for the moment
                        if arg1 == 0:
                            detritus.rename(columns={'EXPORTS_rough': 'taxa'}, inplace=True)
                        else:
                            detritus.rename(columns={'EXPORTS_medium': 'taxa'}, inplace=True)
                        taxa_list = np.unique(detritus.taxa)
                        taxa_list

                        # Now compute number of particles per profile, size class, taxo group and depth bin
                        df = detritus.groupby(['profile', 'depth_part', 'taxa', 'bin_min', 'bin_max']) \
                            .agg({'esd': 'mean', 'vol_sph': 'sum', 'object_id': 'size'}) \
                            .rename(columns={'esd': 'esd_mean', 'vol_sph': 'vol_sph_tot', 'object_id': 'n'}) \
                            .reset_index()

                        # add the volume and compute the concentrations
                        df = pd.merge(df, volumes, how="left", on=['profile', 'depth_part'])
                        df["watervolume"] = df["volume_L"] / 1000  # volume in m3
                        df["conc"] = df["n"] / df["watervolume"]  # concentration in n/m3
                        df["vol_sph"] = df["vol_sph_tot"] / df["watervolume"]  # biovolume concentration in mm3/m3

                        # Add the date
                        list_prof_500 = pd.read_csv(
                            'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/list_of_profiles_angola.csv',
                            sep=';',
                            encoding='unicode_escape')
                        list_prof_500.rename(columns={'date': 'date_ancien'}, inplace=True)
                        list_prof_500.rename(columns={'dateok': 'date'}, inplace=True)
                        df = df.merge(list_prof_500, on=['profile'])

                        ## change the date format from string 20210504.0 to date format 20210504
                        # df['date'] = df['date'].apply(lambda x: str(x))
                        ## get the year, month, and day
                        # df['date'] = df['date'].apply(lambda x: datetime.datetime(int(x[0:4]), int(x[4:6]), int(x[6:8])))
                        ## change it to the 20210504 date format
                        # df['date'] = df['date'].apply(lambda x: x.strftime('%Y%m%d'))
                        ref_date = datetime(2021, 5, 5)

                        df.head()


                        # Plot the distribution of abundance (nb m-3) according to the depth layer, group and time
                        # to compute the max abundance

                        # Gaussian distribution
                        def gaus(x, amp, mu, s):
                            return amp * np.exp(-(x - mu) ** 2 / (2 * s ** 2))


                        # Do it in a loop for all sizes and depth
                        arr = np.empty((0, 7), int)

                        for t in taxa_list:
                            df_taxa = df[df['taxa'] == t]
                            bin_min_list = ecopart_bins[: -1]
                            depth_part_list = np.unique(df_taxa.depth_part)

                            # Add a counter for the b_max selection
                            count_b = 0

                            # Open a pdf to save the figures rough definition
                            if cl == 6:
                                pp = PdfPages(
                                    str(path_to_figures) + '\Ecotaxa_maximum_concentrations' + str(def_chosen) + '_'
                                                                                                                 ' All detritus ' + str(
                                        tp_list_char.date[count_tp]) + ' to ' + str(
                                        tp_list_char.date[count_tp + 1]) + '.pdf')
                            else:
                                pp = PdfPages(
                                    str(path_to_figures) + '\Ecotaxa_maximum_concentrations' + str(
                                        def_chosen) + '_' + str(
                                        t) + ' Cluster' + str(cl) + '; ' + str(
                                        tp_list_char.date[count_tp]) + ' to ' + str(
                                        tp_list_char.date[count_tp + 1]) + '.pdf')

                            for b in bin_min_list[0:-1]:
                                day0 = -1
                                # b = bin_min_list[1]
                                df_temp = df_taxa[df_taxa["bin_min"] == b].reset_index(drop=True)
                                bmax = bin_min_list[count_b + 1]
                                ndepth = len(np.unique(df_temp.depth_part))
                                # For each depth and size cobination, compute the maximum of abundance
                                count_d = 0
                                for d in depth_part_list:
                                    plt.clf()
                                    # d = depth_part_list[0]
                                    df_temp_d = df_temp[df_temp["depth_part"] == d].reset_index(drop=True)
                                    df_temp_d = df_temp_d.groupby(['date', 'depth_part', 'bin_min', 'bin_max']) \
                                        .agg({'conc': 'median'}) \
                                        .reset_index()

                                    if df_temp_d.shape[
                                        0] > 4:  # if we have a median value for at least 9 out of the 20 days
                                        xo = pd.to_datetime(df_temp_d.date, format="%Y%m%d", errors='coerce')
                                        # x = df_temp_d.date
                                        x = (xo - ref_date).dt.days
                                        y = df_temp_d.conc
                                        try:
                                            popt, pcov = curve_fit(gaus, x, y, p0=[max(y), np.mean(x), np.std(x)])
                                            # print(popt)

                                            # better resolution for fit curve representation
                                            curvex = np.linspace(np.min(x), np.max(x), 1000)
                                            curvey = gaus(curvex, *popt)
                                            # get the maximum of abundance
                                            max_curvey = np.max(curvey)
                                            if max_curvey == 0:
                                                max_curvey = np.nan
                                            # the day at which we have this maximum
                                            curvex_df = pd.DataFrame(curvex, columns=['curvex_date']).reset_index()
                                            curvey_df = pd.DataFrame(curvey, columns=['curvey_abund']).reset_index()
                                            curve = pd.merge(curvex_df, curvey_df, on=['index'])
                                            # plt.plot(curve.curvex_date, curve.curvey_abund, 'r', label='gaussian fit')
                                            # plt.show()
                                            max_val = curve[
                                                curve['curvey_abund'] == np.max(curve.curvey_abund)].reset_index(
                                                drop=True)
                                            day = int(max_val.curvex_date[0])
                                            day0 = day

                                            # the std to have the error bar
                                            std_curvex = np.std(x)

                                            # Plot and save it
                                            plt.plot(x, y, 'o', label='Daily median concentration')
                                            plt.plot(curvex, curvey, 'r', label='Gaussian fit' + ', day:' + str(day))
                                            if d < 600:
                                                plt.title(
                                                    str(t) + ' Cluster' + str(cl) + ":" + str(b) + "-" + str(
                                                        bmax) + " mm: " + str(
                                                        d - 50) + "-" + str(d + 50) + " m")
                                            else:
                                                plt.title(
                                                    str(t) + ' Cluster' + str(cl) + ":" + str(b) + "-" + str(
                                                        bmax) + " mm: " + str(
                                                        d - 100) + "-" + str(d + 100) + " m")

                                            plt.legend()
                                            pp.savefig(plt.gcf())  # Save each figure in the pdf

                                        except RuntimeError:
                                            max_val = np.percentile(y, 75)  # Q3
                                            # dayc=x[count_d]-1
                                            # dayarray = np.array([x[count_d],np.min(x[x>dayc]),np.min(x[x>day0])])
                                            dayarray = np.min(x[x > day0])
                                            if np.isnan(dayarray) == False:
                                                day = dayarray.max()
                                                day0 = day
                                                plt.plot(x, y, 'o', label='Daily median concentration')
                                                plt.plot(day, max_val, 'sr',
                                                         label='maximum point' + ', day:' + str(day))
                                                if d < 600:
                                                    plt.title(
                                                        "ERROR ON THE METHOD-" + str(t) + ' Cluster' + str(
                                                            cl) + ":" + str(
                                                            b) + "-" + str(bmax) + " mm: " + str(d - 50) + "-" + str(
                                                            d + 50) + " m")
                                                else:
                                                    plt.title(
                                                        "ERROR ON THE METHOD-" + str(t) + ' Cluster' + str(
                                                            cl) + ":" + str(
                                                            b) + "-" + str(bmax) + " mm: " + str(d - 100) + "-" + str(
                                                            d + 100) + " m")

                                                plt.legend()
                                                pp.savefig(plt.gcf())  # Save each figure in the pdf
                                                max_curvey = np.nan
                                                std_curvex = np.nan
                                                # max_curvey = max_val
                                                # std_curvex = np.nan
                                            else:
                                                day = day0
                                                plt.plot(x, y, 'o', label='Daily median concentration')
                                                plt.plot(day, max_val, 'sr',
                                                         label='maximum point' + ', day:' + str(day))
                                                if d < 600:
                                                    plt.title(
                                                        "ERROR ON THE METHOD-" + str(t) + ' Cluster' + str(
                                                            cl) + ":" + str(
                                                            b) + "-" + str(bmax) + " mm: " + str(d - 50) + "-" + str(
                                                            d + 50) + " m")
                                                else:
                                                    plt.title(
                                                        "ERROR ON THE METHOD-" + str(t) + ' Cluster' + str(
                                                            cl) + ":" + str(
                                                            b) + "-" + str(bmax) + " mm: " + str(d - 100) + "-" + str(
                                                            d + 100) + " m")

                                                plt.legend()
                                                pp.savefig(plt.gcf())  # Save each figure in the pdf
                                                max_curvey = np.nan
                                                std_curvex = np.nan
                                                # max_curvey = max_val
                                                # std_curvex = np.nan

                                        arr = np.append(arr, np.array([[t, b, bmax, d, max_curvey, std_curvex, day]]),
                                                        axis=0)
                                        day0
                                        # print(b,d)
                                    count_d += 1

                                # Add one to the bmax counter
                                count_b += 1

                            # close the pdf document
                            pp.close()

                        arr = pd.DataFrame(arr,
                                           columns=['taxa', 'bin_min', 'bin_max', 'depth', 'max_abun', 'std_curvex',
                                                    'day'])

                        # keep only the depth values below 100m for the linear regression
                        # arr = arr[arr['depth'].astype(int) > 100]
                        np.unique(arr.day)
                        # keep only the data where the maximum is found before the 29th to avoid having detected maximum
                        # which are not true maximums
                        ########################################################################################################################arr = arr[arr['day'].astype(int) < 20210529]
                        # Compute the linear regression on the previously computed points
                        sinking_speed = np.empty((0, 14), int)
                        # Convert the bin_min and bin_max columns in

                        for t in taxa_list:
                            df_taxa = df[df['taxa'] == t]
                            bin_min_list = ecopart_bins[: -1]
                            depth_part_list = np.unique(df_taxa.depth_part)

                            # Add a counter for the b_max selection
                            count_b = 0

                            # Open a pdf to save the figures
                            '''pp = PdfPages(
                                str(path_to_figures) + '\Ecotaxa_sinking_speeds' + str(def_chosen) + '_' + str(
                                    t) + ' Cluster' + str(
                                    cl) +
                                '; ' + str(tp_list_char.date[count_tp]) + ' to ' + str(tp_list_char.date[count_tp + 1]) + '.pdf') '''

                            for b in bin_min_list[0:-1]:
                                # Get the data for the size b and taxa_t in concentration and from the previous computations
                                df_temp = df_taxa[df_taxa["bin_min"] == b].reset_index(drop=True)
                                AAA = round(df_temp['date'], 0)
                                # full['date'] = pd.to_datetime(full['date'], format="%Y%m%d")
                                df_temp['date'] = pd.to_datetime(AAA, format="%Y%m%d", errors='coerce')
                                df_temp['days'] = (df_temp['date'] - ref_date).dt.days
                                arr_taxa = arr[arr["taxa"] == t].reset_index(drop=True)
                                arr_temp = arr_taxa[arr_taxa["bin_min"] == str(b)].reset_index(drop=True)
                                # change format for linear fit
                                arr_temp = arr_temp[arr_temp['day'] != "nan"]  # Remove the nan rows
                                arr_temp['day'] = arr_temp['day'].astype(int)
                                arr_temp['depth'] = arr_temp['depth'].astype(int)
                                arr_temp = arr_temp[arr_temp['max_abun'] != "nan"]  # Remove the nan rows
                                bmax = bin_min_list[count_b + 1]

                                if arr_temp.shape[0] > 3:  # if we have more than 3 points to do the linear regression
                                    if len(np.unique(
                                            arr_temp.day)) > 1:  # if the values give different days for the maximums accross the depth layers
                                        linear_fit = scipy.stats.linregress(arr_temp.day, -arr_temp.depth)
                                        slope = linear_fit.slope
                                        tinv = lambda p, df: abs(scipy.stats.t.ppf(p / 2, df))
                                        ts = tinv(0.05, len(x) - 2)
                                        slope_ci_inf = linear_fit.slope - linear_fit.stderr * ts  # 95% condidence intervall inferior value
                                        slope_ci_sup = linear_fit.slope + linear_fit.stderr * ts  # 95% condidence intervall superior value
                                        intercept = linear_fit.intercept
                                        R2 = linear_fit.rvalue ** 2

                                        linear_fit = scipy.stats.linregress(np.log(arr_temp.depth / 50),
                                                                            np.log(arr_temp.max_abun.astype(float)))
                                        slopeb = linear_fit.slope
                                        tinv = lambda p, df: abs(scipy.stats.t.ppf(p / 2, df))
                                        ts = tinv(0.05, len(x) - 2)
                                        slope_ci_infb = linear_fit.slope - linear_fit.stderr * ts  # 95% condidence intervall inferior value
                                        slope_ci_supb = linear_fit.slope + linear_fit.stderr * ts  # 95% condidence intervall superior value
                                        interceptb = linear_fit.intercept
                                        R2b = linear_fit.rvalue ** 2
                                        F50 = arr_temp.max_abun[arr_temp.depth == 50]

                                        # save this data in a table
                                        sinking_speed = np.append(sinking_speed, np.array([[def_chosen, t, b, bmax,
                                                                                            slope, slope_ci_inf,
                                                                                            slope_ci_sup,
                                                                                            intercept, R2, slopeb,
                                                                                            slope_ci_infb,
                                                                                            slope_ci_supb, F50[0],
                                                                                            R2b]]),
                                                                  axis=0)

                                        # select the values of max of abundance computed above
                                        iday = arr_temp.day
                                        idays = iday.reset_index(drop=True)
                                        date_min = idays[0]
                                        date_max = arr_temp.day.iloc[-1]

                                        # 3. Interpolate the concentration field
                                        # dat = np.array(df_temp['date'])
                                        dat = np.array(df_temp['days'])
                                        unique_pos_dat = np.unique(dat)
                                        depth = np.array(df_temp['depth'] * (-1))
                                        param = np.array(df_temp['conc'])
                                        # choose the limits in depth and time
                                        depth_min_max = [0, -1000]
                                        # pos_min_max = [min(df_temp.date), max(df_temp.date)]  # Date transect
                                        pos_min_max = [min(df_temp.days), max(df_temp.days)]  # Date transect
                                        # interpolate
                                        xi, yi, zi = gridding_func(pos_min_max, depth_min_max, dat, depth, param)
                                        levels = 15
                                        min_contour_level = min(df_temp.conc)
                                        max_contour_level = max(df_temp.conc) / 10
                                        contour_levels = contour_levels_func(min_contour_level, max_contour_level,
                                                                             levels)
                                        date_lab = [datetime(2021, 5, 5) + timedelta(days=i) for i in
                                                    xi]  # changer ma date dedans
                                        # Format each date in the list as a string in the format "YYYY-MM-DD"
                                        formatted_dates = [d.strftime("%Y-%m-%d") for d in date_lab]
                                        a_lab = [datetime(2021, 5, 5) + timedelta(days=i) for i in
                                                 arr_temp.day]  # changer ma date dedans
                                        # Format each date in the list as a string in the format "YYYY-MM-DD"
                                        formatted_dates_arr = [d.strftime("%Y-%m-%d") for d in a_lab]

                                        ''''# 4. Plot
                                        plt.clf()
                                        ax1 = plt.subplot(1, 1, count)
                                        # Concentration field
                                        # p1 = plt.contourf(xi, yi, zi, contour_levels, cmap=cm.viridis, alpha=1, extend="both",norm=colors.LogNorm(vmin=np.nanmin(zi), vmax=np.nanmax(zi) / 5))
                                        p1 = plt.contourf(formatted_dates, yi, zi, contour_levels, cmap='viridis', alpha=1,
                                                          extend="both",
                                                          norm=colors.LogNorm(vmin=np.nanmin(zi), vmax=np.nanmax(zi) / 5))
                                        # Maximum of concentration
                                        p2 = plt.scatter(x=formatted_dates_arr, y=- arr_temp.depth, s=15, c="black")
                                        # Linear regression
                                        p3 = plt.plot(np.linspace(date_min, date_max, 20),
                                                      np.linspace(date_min * slope + intercept, date_max * slope + intercept, 20),
                                                      c='black')
                                        # ax1.hlines(y=300, xmin=min(formatted_dates), xmax=max(formatted_dates), colors='r')
                                        cb = plt.colorbar(p1, orientation='vertical')
                                        tick_locator = ticker.MaxNLocator(nbins=5)
                                        cb.locator = tick_locator
                                        cb.ax.set_ylabel("concentration [nb m-3]", fontsize=10)
                                        plt.ylabel("Depth [m]")
                                        # plt.xlim(min(df_temp.date), max(df_temp.date))
                                        plt.xlim(min(formatted_dates), max(formatted_dates))
                                        ax1.set_xticks(formatted_dates[::60])
                                        plt.ylim(-1000, 0)
                                        plt.xlabel("Time [days]")
                                        plt.xticks(rotation=70)
                                        plt.title(str(t) + ' Cluster' + str(cl) + ":" + str(b) + "-" + str(
                                            bmax) + " mm ; sinking speed=" + str(
                                            - round(slope, 1)) + " m d-1 (R2=" + str(round(R2, 2)) + ")")
                                        pp.savefig(plt.gcf())  # Save each figure in the pdf '''

                                # Add one to the bmax counter
                                count_b += 1

                            # close the pdf document
                            '''pp.close() '''

                        sinking_speed = pd.DataFrame(sinking_speed, columns=['definition', 'taxa', 'bin_min', 'bin_max',
                                                                             'slope', 'slope_ci_inf', 'slope_ci_sup',
                                                                             'intercept', 'R2',
                                                                             'attenuation', 'attenuation_ci_inf',
                                                                             'attenuation_ci_sup',
                                                                             'Abundance50m', 'attenuationR2'])
                        # save it
                        if cl == 6:
                            arr.to_csv(
                                str(path_to_SPEED) + '\Info_sinking_' + savefile[q] + def_chosen + ' ' + ' All detritus' + str(
                                    tp_list_char.date[count_tp]) + ' to ' + str(
                                    tp_list_char.date[count_tp + 1]) + '.csv',
                                index=False)
                            sinking_speed.to_csv(
                                str(path_to_figures) + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' All detritus' + str(
                                    tp_list_char.date[count_tp]) + ' to ' + str(
                                    tp_list_char.date[count_tp + 1]) + '.csv',
                                index=False)
                            sinking_speed.to_csv(
                                str(path_to_SPEED) + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' All detritus' + str(
                                    tp_list_char.date[count_tp]) + ' to ' + str(
                                    tp_list_char.date[count_tp + 1]) + '.csv',
                                index=False)
                        else:
                            arr.to_csv(
                                str(path_to_SPEED) + '\Info_sinking_' + savefile[q] + def_chosen + ' ' + ' Cluster' + str(cl) + str(
                                    tp_list_char.date[count_tp]) + ' to ' + str(
                                    tp_list_char.date[count_tp + 1]) + '.csv',
                                index=False)
                            sinking_speed.to_csv(
                                str(path_to_figures) + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' Cluster' + str(
                                    cl) + str(
                                    tp_list_char.date[count_tp]) + ' to ' + str(
                                    tp_list_char.date[count_tp + 1]) + '.csv',
                                index=False)
                            sinking_speed.to_csv(
                                str(path_to_SPEED) + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' Cluster' + str(
                                    cl) + str(
                                    tp_list_char.date[count_tp]) + ' to ' + str(
                                    tp_list_char.date[count_tp + 1]) + '.csv',
                                index=False)

                        # Dodji added speed plot

                        if cl == 6:
                            vitesses = pd.read_csv(
                                str(path_to_figures) + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' All detritus' +
                                str(tp_list_char.date[count_tp]) + ' to ' + str(
                                    tp_list_char.date[count_tp + 1]) + '.csv')
                        else:
                            vitesses = pd.read_csv(
                                str(path_to_figures) + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' Cluster' + str(
                                    cl) +
                                str(tp_list_char.date[count_tp]) + ' to ' + str(
                                    tp_list_char.date[count_tp + 1]) + '.csv')
                        # x_vit = (vitesses.bin_min + vitesses.bin_max) / 2

                        st1 = vitesses.bin_min
                        st2 = vitesses.bin_max
                        st1s = st1.astype("string")
                        st2s = st2.astype("string")
                        bin_label = st1s + '-' + st2s
                        x_vit = bin_label
                        y_vit = -vitesses.slope
                        y_vit[y_vit < 0] = np.NaN

                        if cl == 6:
                            pp = PdfPages(str(path_to_figures) + '\ vitesse_graphique_allsize' + savefile[q] + str(
                                def_chosen) + '_' + ' All detritus' + '; ' + str(
                                tp_list_char.date[count_tp]) + ' to ' + str(tp_list_char.date[count_tp + 1]) + '.pdf')
                        else:
                            pp = PdfPages(
                                str(path_to_figures) + '\ vitesse_graphique_allsize' + savefile[q] + str(def_chosen) + '_' + str(
                                    t) + ' Cluster' + str(cl) + '; ' + str(tp_list_char.date[count_tp]) + ' to ' + str(
                                    tp_list_char.date[count_tp + 1]) + '.pdf')

                        plt.clf()
                        ax1 = plt.subplot(1, 1, count)
                        plt.plot(x_vit, y_vit, 'b')
                        plt.plot(x_vit, y_vit, 'ob', label='Speeds(m/d)', linestyle='--', markersize=3.5)
                        plt.ylabel("Speed [m/d]")
                        plt.xlabel("Size classes (mm)")
                        plt.title("Speeds variation with size classes " + ' Cluster' + str(cl) + str(
                            tp_list_char.date[count_tp]) + ' to ' + str(tp_list_char.date[count_tp + 1]))
                        plt.legend()
                        # ax1.set_xticklabels(bin_label[0:4])
                        plt.xticks(rotation=20)
                        pp.savefig(plt.gcf(), bbox_inches='tight')  # Save each figure in the pdf
                        pp.close()
                        plt.close()

                        count_tp += 1

        path_to_fig = Path(
            'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/ClusterMiriamdoc/SPEED/SPEED').expanduser()

        vitesses = pd.read_csv(
            str(path_to_fig) + '1' + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' Cluster1' +
            str(tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.csv')

        st1 = vitesses.bin_min
        st2 = vitesses.bin_max
        st1s = st1.astype("string")
        st2s = st2.astype("string")
        bin_label = st1s + '-' + st2s
        x_vit1 = bin_label
        y_vit1 = -vitesses.slope
        y_vit1[y_vit1 < 0] = np.NaN

        vitesses = pd.read_csv(
            str(path_to_fig) + '2' + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' Cluster2' +
            str(tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.csv')

        st1 = vitesses.bin_min
        st2 = vitesses.bin_max
        st1s = st1.astype("string")
        st2s = st2.astype("string")
        bin_label = st1s + '-' + st2s
        x_vit2 = bin_label
        y_vit2 = -vitesses.slope
        y_vit2[y_vit2 < 0] = np.NaN

        vitesses = pd.read_csv(
            str(path_to_fig) + '3' + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' Cluster3' +
            str(tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.csv')

        st1 = vitesses.bin_min
        st2 = vitesses.bin_max
        st1s = st1.astype("string")
        st2s = st2.astype("string")
        bin_label = st1s + '-' + st2s
        x_vit3 = bin_label
        y_vit3 = -vitesses.slope
        y_vit3[y_vit3 < 0] = np.NaN

        vitesses = pd.read_csv(
            str(path_to_fig) + '4' + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' Cluster4' +
            str(tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.csv')

        st1 = vitesses.bin_min
        st2 = vitesses.bin_max
        st1s = st1.astype("string")
        st2s = st2.astype("string")
        bin_label = st1s + '-' + st2s
        x_vit4 = bin_label
        y_vit4 = -vitesses.slope
        y_vit4[y_vit4 < 0] = np.NaN

        if np.size(CL0) > 1:
            vitesses = pd.read_csv(
                str(path_to_fig) + '5' + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' Cluster5' +
                str(tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.csv')

            st1 = vitesses.bin_min
            st2 = vitesses.bin_max
            st1s = st1.astype("string")
            st2s = st2.astype("string")
            bin_label = st1s + '-' + st2s
            x_vit5 = bin_label
            y_vit5 = -vitesses.slope
            y_vit5[y_vit5 < 0] = np.NaN

        vitesses = pd.read_csv(
            str(path_to_fig) + '6' + '\Ecotaxa_sinking_speed_' + savefile[q] + def_chosen + ' ' + ' All detritus' +
            str(tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.csv')
        st1 = vitesses.bin_min
        st2 = vitesses.bin_max
        st1s = st1.astype("string")
        st2s = st2.astype("string")
        bin_label = st1s + '-' + st2s
        x_vit6 = bin_label
        y_vit6 = -vitesses.slope
        y_vit6[y_vit6 < 0] = np.NaN

        pp = PdfPages(str(path_to_fig) + '\ vitesse_graphique_allsize' + savefile[q] + '_' + 'Detritus ' + ' Clusters ' + ' ' + str(
            tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.pdf')

        plt.clf()
        ax1 = plt.subplot(1, 1, count)
        p1 = plt.plot(x_vit1, y_vit1, 'ob', label='Cluster1 Speeds(m/d)', linestyle='--', markersize=3.5)
        p2 = plt.plot(x_vit2, y_vit2, 'or', label='Cluster2 Speeds(m/d)', linestyle='--', markersize=3.5)
        p3 = plt.plot(x_vit3, y_vit3, 'og', label='Cluster3 Speeds(m/d)', linestyle='--', markersize=3.5)
        p4 = plt.plot(x_vit4, y_vit4, 'ok', label='Cluster4 Speeds(m/d)', linestyle='--', markersize=3.5)
        if np.size(CL0) > 1:
            p5 = plt.plot(x_vit5, y_vit5, 'oc', label='All detritus Speeds(m/d)', linestyle='--', markersize=3.5)

        p6 = plt.plot(x_vit6, y_vit6, 'om', label='All detritus Speeds(m/d)', linestyle='--', markersize=3.5)

        plt.ylabel("Speed [m/d]")
        plt.xlabel("Size classes (mm)")
        plt.title("Speeds variation with size classes " + ' Clusters ' + str(tp_list_char.date[0]) + ' to ' + str(
            tp_list_char.date[1]))
        plt.legend()
        # ax1.set_xticklabels(bin_label[0:4])
        plt.xticks(rotation=20)
        pp.savefig(plt.gcf(), bbox_inches='tight')  # Save each figure in the pdf
        plt.rcParams.update({'font.size': 20})
        pp.close()
    plt.close('all')

''' detritus = detritus0[detritus0['date']> 20210715].reset_index(drop = True)
    detritus = detritus0[detritus0['date']< 20210819].reset_index(drop = True)

    detritus = detritus0[detritus0['date']> 20210819].reset_index(drop = True)
    detritus = detritus0[detritus0['date']< 20211001].reset_index(drop = True)

    detritus = detritus0[detritus0['date']> 20211001].reset_index(drop = True)
    detritus = detritus0[detritus0['date']< 20211030].reset_index(drop = True)

    detritus = detritus0[detritus0['date']> 20211113].reset_index(drop = True)
    detritus = detritus0[detritus0['date']< 20211215].reset_index(drop = True)

    detritus = detritus0[detritus0['date']> 20211215].reset_index(drop = True)
    detritus = detritus0[detritus0['date']< 20220110].reset_index(drop = True)

    detritus = detritus0[detritus0['date']> 20220116].reset_index(drop = True)
    detritus = detritus0[detritus0['date']< 20220225].reset_index(drop = True)

    detritus = detritus0[detritus0['date']> 20220321].reset_index(drop = True)
    detritus = detritus0[detritus0['date']< 20220424].reset_index(drop = True) '''

    # Script Alberto
    # import scipy.stats
    # import numpy as np
    #
    # def lin_fit(x,y):
    #     result = scipy.stats.linregress(x, y)
    #     tinv = lambda p, df: abs(scipy.stats.t.ppf(p/2, df))
    #     ts = tinv(0.05, len(x)-2)
    #     slpe_ci=np.array([0.0,0.0]);intrcpt_ci=np.array([0.0,0.0])
    #     slpe_ci[0] = result.slope - result.stderr * ts
    #     slpe_ci[1] = result.slope + result.stderr * ts
    #     intrcpt_ci[0] = result.intercept - result.intercept_stderr * ts
    #     intrcpt_ci[1] = result.intercept + result.intercept_stderr * ts
    #     R2=result.rvalue**2
    #     t = np.sqrt(R2 * (len(x) - 2)) / np.sqrt(1 - R2)
    #     t1=scipy.stats.t.ppf(0.975, len(x) - 2)
    #     t2 = scipy.stats.t.ppf(0.995, len(x) - 2)
    #     t3 = scipy.stats.t.ppf(0.9995, len(x) - 2)
    #
    #     if t < t1:
    #         signif = 0;signif_label=''
    #     elif t < t2:
    #         signif = 1;signif_label='*'
    #     elif t < t3:
    #         signif = 2;signif_label='**'
    #     elif t >= t3:
    #         signif = 3;signif_label='***'
    #     return result,slpe_ci,intrcpt_ci,signif,signif_label
    #
    # #if __name__ == "__main__":
    # #    x = np.array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    # #    y = np.array([2, 1, 4, 5, 8, 12, 18, 25, 96, 48])
    # #    (result,slpe_ci,intrcpt_ci,signif, signif_label) = lin_fit(x, y)
    # #    print(result)
    # #    print('significativity of fit is %s (%s)' % signif,signif_label)











