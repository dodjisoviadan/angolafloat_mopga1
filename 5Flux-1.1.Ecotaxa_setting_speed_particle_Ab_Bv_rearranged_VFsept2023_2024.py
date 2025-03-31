
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
path_to_data = Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/DEPTH_rearranged/ECOPART').expanduser()
os.chdir(str(path_to_data))

# import abundance biovolume of ecopart
ecopart = pd.read_table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/ECOPART_2023-05-16-15h_03min_PAR_Aggregated.tsv", sep='\t', encoding='unicode_escape')
ecopart = ecopart[ecopart['Profile'].str.contains('a')]
listsize0 = list(ecopart.columns)
ecopartAB = ecopart.iloc[:,[0, 2, 4, 12,13,14,15,16,17,18,19]] # abundance /l
ecopartBV = ecopart.iloc[:,[0, 2, 4, 27,28,29,30,31,32,33,34]] # biovolume /l
taxa_list=['OK']

ecopart_bins = [0.064, 0.128, 0.256, 0.512, 1.02, 2.05, 4.1, 8.19, 16.4]  # bins

# bin the data per 100m
def rounditup(x, precision, method):
    import math
    if method == "floor":
        return math.floor(x / precision) * precision
    elif method == "ceiling":
        return math.ceil(x / precision) * precision
    else:
        return "give the parameter floor or ceiling"


ecopart_profiles = ecopart["Profile"].unique()
ecopart_profiles = pd.DataFrame(ecopart_profiles, columns=['Profile'])
ecopart_profiles = ecopart_profiles.sort_values('Profile')
ecopart_profiles['CYCLE_NUMBER'] = list(range(1, 119))

tp_list = [20211001, 20211113]  # 2
tp_list = [20210713, 20210819] # 0
tp_list = [20210819, 20210925]  # 1
tp_list = [20211113, 20211215]  # 3
tp_list = [20211113, 20220119]  # 3b
tp_list = [20220120, 20220306]  # 4
tp_list = [20220306, 20220426]  # 5

path_to_SPEED= Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/DEPTH_rearranged/ECOPART').expanduser()
path_to_SPEED= Path('C:\MOPGA 2022-2023 Dodji\ARGO float Angola\EcotaxaPart_BCG_UVP6_Plot\ClusterCenterQselection\SPEED\AB_or_BV\A2024-data').expanduser()
tp_lista = [[20210713, 20210819],[20210819, 20210925],[20211001, 20211113],[20211113, 20211215],[20211113, 20220119],[20220120, 20220306],[20220306, 20220426]] # 0
#tp_lista = [[20211001, 20211113]] # 0
tp_lista = [[20210819, 20210925],[20211001, 20211113],[20220120, 20220306]] # 0
tp_lista = [[20210819, 20210925],[20211001, 20211113],[20211113, 20211215],[20220120, 20220306],[20220310, 20220416]] # 0

tp_lista = [[20210720, 20210815],[20210819, 20210920], [20211015, 20211110], [20211120, 20211215], [20220130, 20220225],[20220315, 20220410]]  # 0
tp_lista = [[20210717, 20210819],[20210822, 20210924], [20211008, 20211113], [20211123, 20211231], [20220128, 20220309],[20220313, 20220412]]  # 0

PartABV=['AB','BV']#

for tp_list in tp_lista:

    tp_list_char = pd.to_datetime(tp_list, format="%Y%m%d", errors='coerce')  #

    for cl in PartABV:

        ecopartchoose=eval('ecopart' + cl)

        ecopartchoose.loc[:, 'depth_part'] = 99999
        ecopartchoose['depth_part'] = ecopartchoose["Depth [m]"].apply(
            lambda x: rounditup(x, 100, "floor") + 50 if x < 600 else rounditup(x, 200, "floor") + 100)

        '''path_to_figures = Path(
            'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/DEPTH_rearranged/ECOPART ' + str(
                cl)).expanduser()'''

        if cl=='AB':
            path_to_figures = Path(
                'C:\MOPGA 2022-2023 Dodji\ARGO float Angola\EcotaxaPart_BCG_UVP6_Plot\ClusterCenterQselection\SPEED\AB_or_BV\AB\A2024').expanduser()
        else:
            path_to_figures = Path(
                'C:\MOPGA 2022-2023 Dodji\ARGO float Angola\EcotaxaPart_BCG_UVP6_Plot\ClusterCenterQselection\SPEED\AB_or_BV\BV\A2024').expanduser()



        arg1 = 0  # For the defitition choice of taxo definition : either 0 for EXPORTS_rough or 1 for EXPORTS_medium
        argb = 0

        count_tp = 0
        count = 1
        fig = plt.figure(figsize=(19.20, 9.83))  # (5, 3.5)
        plt.rcParams.update({'font.size': 12})
        ax1 = fig.add_axes([0, 0, 0, 0])

        for tp in tp_list[0:-1]:
            ecopartchoose['yyyy-mm-dd hh:mm'] = pd.to_datetime(ecopartchoose['yyyy-mm-dd hh:mm']) #.loc[row_indexer,col_indexer] = value instead
            #ecopartchoose.drop(['yyyy-mm-dd hh:mm'], inplace=True, axis=1)
            #ecopartchoose.rename(columns={'yyyy-mm-dd hh:mmx': 'yyyy-mm-dd hh:mm'}, inplace=True)
            detritus = ecopartchoose[ecopartchoose['yyyy-mm-dd hh:mm'] > tp_list_char[count_tp]].reset_index(drop=True)
            detritus = detritus[detritus['yyyy-mm-dd hh:mm']< tp_list_char[count_tp + 1]].reset_index(drop=True)
            detritus = detritus[detritus['Depth [m]'] < 500].reset_index(drop=True)

            detritus.shape


            # Now compute number and volumes of particles per profile, size class, taxo group and depth bin
            df = detritus

            ref_date = datetime(2021, 5, 5)

            df.head()


            # Plot the distribution of abundance (nb m-3) according to the depth layer, group and time
            # to compute the max abundance

            # Gaussian distribution
            def gaus(x, amp, mu, s):
                return amp * np.exp(-(x - mu) ** 2 / (2 * s ** 2))


            # Do it in a loop for all sizes and depth
            arr = np.empty((0, 8), int)

            for t in taxa_list:
                df_taxa = df
                bin_min_list = ecopart_bins
                depth_part_list = np.unique(df_taxa.depth_part)

                # Open a pdf to save the figures rough definition
                pp = PdfPages(str(path_to_figures) + '\Ecotaxa_maximum_concentrations' +  '_' + str(
                    t) + ' Particle' + str(cl) +
                              '; ' + str(tp_list_char.date[count_tp]) + ' to ' + str(
                    tp_list_char.date[count_tp + 1]) + '.pdf')
                # Add a counter for the b_max selection
                count_b = 0

                for b in bin_min_list[0:-1]:
                    countb=count_b+3
                    day0 = -1
                    # b = bin_min_list[1]
                    listsize = list(df_taxa.columns)
                    df_temp = df_taxa # .iloc[:,countb]*0.001
                    df_temp[listsize[countb]] = df_taxa[listsize[countb]]*1000  # abundance or biovolume / m3
                    bmax = bin_min_list[count_b + 1]
                    # For each depth and size cobination, compute the maximum of abundance
                    count_d = 0
                    for d in depth_part_list:
                        plt.clf()
                        # d = depth_part_list[0]
                        df_temp_d = df_temp[df_temp["depth_part"] == d].reset_index(drop=True)
                        df_temp_d = df_temp_d.groupby(['yyyy-mm-dd hh:mm', 'depth_part']) \
                            .agg({listsize[countb]: 'median'}) \
                            .reset_index()
                        df_temp_d['date'] = df_temp_d['yyyy-mm-dd hh:mm']
                        df_temp_d['varr'] = df_temp_d[listsize[countb]]
                        if df_temp_d.shape[0] > 4:  # about time series point or median values but can vary: if we have a median value for at least 9 out of the 20 days
                            xo = pd.to_datetime(df_temp_d.date, format="%Y%m%d", errors='coerce')
                            # x = df_temp_d.date
                            x = (xo - ref_date).dt.days
                            y = df_temp_d.varr
                            try:
                                popt, pcov = curve_fit(gaus, x, y, p0=[max(y), np.mean(x), np.std(x)])
                                # print(popt)
                                # now compute de R2 of the gaussian fit
                                modelPredictions = gaus(x, *popt)
                                absError = modelPredictions - y
                                SE = np.square(absError)  # squared errors
                                MSE = np.mean(SE)  # mean squared errors
                                RMSE = np.sqrt(MSE)  # Root Mean Squared Error, RMSE
                                if np.var(y) == 0:
                                    Rsquared = np.nan
                                else:
                                    Rsquared = 1.0 - (np.var(absError) / np.var(y))

                                # better resolution for fit curve representation
                                curvex = np.linspace(np.min(x), np.max(x), 1000)
                                curvey = gaus(curvex, *popt)
                                # get the maximum of abundance
                                max_curvey = np.max(curvey)
                                if max_curvey==0:
                                    max_curvey=np.nan
                                    Rsquared = np.nan

                                # the day at which we have this maximum
                                curvex_df = pd.DataFrame(curvex, columns=['curvex_date']).reset_index()
                                curvey_df = pd.DataFrame(curvey, columns=['curvey_abund']).reset_index()
                                curve = pd.merge(curvex_df, curvey_df, on=['index'])
                                # plt.plot(curve.curvex_date, curve.curvey_abund, 'r', label='gaussian fit')
                                # plt.show()
                                max_val = curve[curve['curvey_abund'] == np.max(curve.curvey_abund)].reset_index(
                                    drop=True)
                                day = int(max_val.curvex_date[0])
                                day0 = day

                                # the std to have the error bar
                                std_curvex = np.std(x)

                                # Plot and save it
                                plt.plot(x, y, 'o', label='Daily median concentration')
                                plt.plot(curvex, curvey, 'r', label='Gaussian fit' + ', day:' + str(day))
                                if d<600:
                                    plt.title(
                                        str(t) + ' Particle' + str(cl) + ":" + str(b) + "-" + str(bmax) + " mm: " + str(
                                            d - 50) + "-" + str(d + 50) + " m")
                                else:
                                    plt.title(
                                        str(t) + ' Particle' + str(cl) + ":" + str(b) + "-" + str(bmax) + " mm: " + str(
                                            d - 100) + "-" + str(d + 100) + " m")

                                plt.legend()
                                pp.savefig(plt.gcf())  # Save each figure in the pdf

                            except RuntimeError:
                                max_val = np.percentile(y, 75)  # Q3
                                # dayc=x[count_d]-1
                                # dayarray = np.array([x[count_d],np.min(x[x>dayc]),np.min(x[x>day0])])
                                dayarray = np.min(x[x > day0])
                                if np.isnan(dayarray) == False:
                                    day = np.max(dayarray) #dayarray.max()
                                    day0 = day
                                    plt.plot(x, y, 'o', label='Daily median concentration')
                                    plt.plot(day, max_val, 'sr', label='maximum point' + ', day:' + str(day))
                                    if d<600:
                                        plt.title(
                                            "ERROR ON THE METHOD-" + str(t) + ' Particle' + str(cl) + ":" + str(
                                                b) + "-" + str(bmax) + " mm: " + str(d - 50) + "-" + str(d + 50) + " m")
                                    else:
                                        plt.title(
                                            "ERROR ON THE METHOD-" + str(t) + ' Particle' + str(cl) + ":" + str(
                                                b) + "-" + str(bmax) + " mm: " + str(d - 100) + "-" + str(d + 100) + " m")

                                    plt.legend()
                                    pp.savefig(plt.gcf())  # Save each figure in the pdf
                                    max_curvey = np.nan
                                    std_curvex = np.nan
                                    Rsquared = np.nan
                                    # max_curvey = max_val
                                    # std_curvex = np.nan
                                else:
                                    day = day0
                                    plt.plot(x, y, 'o', label='Daily median concentration')
                                    plt.plot(day, max_val, 'sr', label='maximum point' + ', day:' + str(day))
                                    if d<600:
                                        plt.title(
                                            "ERROR ON THE METHOD-" + str(t) + ' Particle' + str(cl) + ":" + str(
                                                b) + "-" + str(bmax) + " mm: " + str(d - 50) + "-" + str(d + 50) + " m")
                                    else:
                                        plt.title(
                                            "ERROR ON THE METHOD-" + str(t) + ' Particle' + str(cl) + ":" + str(
                                                b) + "-" + str(bmax) + " mm: " + str(d - 100) + "-" + str(d + 100) + " m")

                                    plt.legend()
                                    pp.savefig(plt.gcf())  # Save each figure in the pdf
                                    max_curvey = np.nan
                                    std_curvex = np.nan
                                    Rsquared = np.nan
                                    # max_curvey = max_val
                                    # std_curvex = np.nan

                            arr = np.append(arr, np.array([[t, b, bmax, d, max_curvey, std_curvex, Rsquared, day]]), axis=0)
                            day0
                            # print(b,d)
                        count_d += 1

                    # Add one to the bmax counter
                    count_b += 1

                # close the pdf document
                pp.close()

            arr = pd.DataFrame(arr, columns=['taxa', 'bin_min', 'bin_max', 'depth', 'max_abun', 'std_curvex', 'R2gauss', 'day'])

            # respect here the period (6 days) of the float and the order in time space (t1<t2)
            '''arr.insert(1, 'daycorrect', 'ans')
            BN = np.unique(arr.bin_min)
            NPa = []
            count_bn = 0
            for bn in BN:
                arr_bn = arr[arr["bin_min"] == bn].reset_index(drop=True)
                NP = np.size(arr_bn, 0)
                NPa.append(NP)
                count_da = 0
                for da in arr_bn['day']:
                    if count_da == 0:
                        dx = da
                    else:
                        nexto = int(arr_bn['day'][0]) + count_da * 100
                        if int(da) <= nexto:
                            if int(da) > int(arr_bn['day'][count_da - 1]):
                                dx = da
                            else:
                                dx = np.nan
                        else:
                            dx = np.nan
                    arr['daycorrect'][count_da + np.sum(NPa[0:-1])] = str(dx)
                    count_da += 1
                count_bn += 1

            arr.drop(['day'], inplace=True, axis=1)
            arr.rename(columns={'daycorrect': 'day'}, inplace=True) '''

            # keep only the depth values below 100m for the linear regression
            # arr = arr[arr['depth'].astype(int) > 100]
            np.unique(arr.day)
            # keep only the data where the maximum is found before the 29th to avoid having detected maximum
            # which are not true maximums
            ########################################################################################################################arr = arr[arr['day'].astype(int) < 20210529]
            # Compute the linear regression on the previously computed points
            sinking_speed = np.empty((0, 13), int)
            # Convert the bin_min and bin_max columns in

            for t in taxa_list:
                df_taxa = df
                bin_min_list = ecopart_bins[: -1]
                depth_part_list = np.unique(df_taxa.depth_part)

                # Add a counter for the b_max selection
                count_b = 0

                # Open a pdf to save the figures
                '''pp = PdfPages(
                    str(path_to_figures) + '\Ecotaxa_sinking_speeds'  + '_' + str(
                        t) + ' Particle' + str(
                        cl) +
                    '; ' + str(tp_list_char.date[count_tp]) + ' to ' + str(tp_list_char.date[count_tp + 1]) + '.pdf')'''

                for b in bin_min_list[0:-1]:
                    # Get the data for the size b and taxa_t in concentration and from the previous computations
                    df_temp = df_taxa # .reset_index(drop=True)  # abundance or biovolume / m3
                    AAA = df_temp['yyyy-mm-dd hh:mm']
                    df_temp['date'] = pd.to_datetime(AAA, format="%Y%m%d", errors='coerce')
                    df_temp['days'] = (df_temp['date'] - ref_date).dt.days
                    arr_taxa = arr
                    arr_temp = arr_taxa[arr_taxa.bin_min == str(b)].reset_index(drop=True)

                    # change format for linear fit
                    arr_temp = arr_temp[arr_temp['day'] != "nan"]  # Remove the nan rows
                    arr_temp['day'] = arr_temp['day'].astype(int)
                    arr_temp['depth'] = arr_temp['depth'].astype(int)
                    arr_temp = arr_temp[arr_temp['max_abun'] != "nan"]  # Remove the nan rows
                    bmax = bin_min_list[count_b + 1]

                    if arr_temp.shape[0] >= 2:  # if we have more than 3 points to do the linear regression
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

                            linear_fit = scipy.stats.linregress(np.log(arr_temp.depth/50), np.log(arr_temp.max_abun.astype(float)))
                            slopeb = linear_fit.slope
                            tinv = lambda p, df: abs(scipy.stats.t.ppf(p / 2, df))
                            ts = tinv(0.05, len(x) - 2)
                            slope_ci_infb = linear_fit.slope - linear_fit.stderr * ts  # 95% condidence intervall inferior value
                            slope_ci_supb = linear_fit.slope + linear_fit.stderr * ts  # 95% condidence intervall superior value
                            interceptb = linear_fit.intercept
                            R2b = linear_fit.rvalue ** 2
                            #F50=arr_temp.max_abun[arr_temp.depth==50]
                            F50 = arr_temp.max_abun[arr_temp.depth == np.min(arr_temp['depth'])]

                            # save this data in a table
                            sinking_speed = np.append(sinking_speed, np.array([[t, b, bmax,
                                                                                slope, slope_ci_inf, slope_ci_sup,
                                                                                intercept, R2, slopeb, slope_ci_infb, slope_ci_supb, np.max(F50), R2b]]), axis=0)

                            ''''# select the values of max of abundance computed above
                            iday = arr_temp.day
                            idays = iday.reset_index(drop=True)
                            date_min = np.min(idays)
                            date_max = np.max(idays)

                            # 3. Interpolate the concentration field
                            # dat = np.array(df_temp['date'])
                            dat = np.array(df_temp['days'])
                            unique_pos_dat = np.unique(dat)
                            depth = np.array(df_temp['Depth [m]'] * (-1))
                            param = np.array(df_temp[listsize[countb]])
                            # choose the limits in depth and time
                            depth_min_max = [0, -1000]
                            # pos_min_max = [min(df_temp.date), max(df_temp.date)]  # Date transect
                            pos_min_max = [min(df_temp.days), max(df_temp.days)]  # Date transect
                            # interpolate
                            xi, yi, zi = gridding_func(pos_min_max, depth_min_max, dat, depth, param)
                            levels = 15
                            min_contour_level = min(df_temp[listsize[countb]])
                            max_contour_level = max(df_temp[listsize[countb]]) / 10
                            contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)
                            date_lab = [datetime(2021, 5, 5) + timedelta(days=i) for i in xi]  # changer ma date dedans
                            # Format each date in the list as a string in the format "YYYY-MM-DD"
                            formatted_dates = [d.strftime("%Y-%m-%d") for d in date_lab]
                            a_lab = [datetime(2021, 5, 5) + timedelta(days=i) for i in
                                     arr_temp.day]  # changer ma date dedans
                            # Format each date in the list as a string in the format "YYYY-MM-DD"
                            formatted_dates_arr = [d.strftime("%Y-%m-%d") for d in a_lab]

                            # 4. Plot
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
                            plt.title(str(t) + ' Particle' + str(cl) + ":" + str(b) + "-" + str(
                                bmax) + " mm ; sinking speed=" + str(
                                - round(slope, 1)) + " m d-1 (R2=" + str(round(R2, 2)) + ")")
                            pp.savefig(plt.gcf())  # Save each figure in the pdf '''

                    # Add one to the bmax counter
                    count_b += 1

                # close the pdf document
                #pp.close()

            sinking_speed = pd.DataFrame(sinking_speed, columns=['taxa', 'bin_min', 'bin_max',
                                                                 'slope', 'slope_ci_inf', 'slope_ci_sup',
                                                                 'intercept', 'R2',
                                                                 'attenuation', 'attenuation_ci_inf', 'attenuation_ci_sup',
                                                                 'Abundance50m', 'attenuationR2'])
            # save it
            arr.to_csv(
                str(path_to_SPEED) + '\Info_sinking_' + 'Particle' + str(cl) + str(
                    tp_list_char.date[count_tp]) + ' to ' + str(tp_list_char.date[count_tp + 1]) + '.csv', index=False)
            sinking_speed.to_csv(
                str(path_to_figures) + '\Ecotaxa_sinking_speed_'  + 'Particle' + str(cl) + str(
                    tp_list_char.date[count_tp]) + ' to ' + str(tp_list_char.date[count_tp + 1]) + '.csv', index=False)
            sinking_speed.to_csv(
                str(path_to_SPEED) + '\Ecotaxa_sinking_speed_'  + 'Particle' + str(cl) + str(
                    tp_list_char.date[count_tp]) + ' to ' + str(tp_list_char.date[count_tp + 1]) + '.csv', index=False)
            sinking_speed;

            # Dodji added speed plot

            vitesses = pd.read_csv(
                str(path_to_figures) + '\Ecotaxa_sinking_speed_' + 'Particle' + str(cl) +
                str(tp_list_char.date[count_tp]) + ' to ' + str(tp_list_char.date[count_tp + 1]) + '.csv')
            # x_vit = (vitesses.bin_min + vitesses.bin_max) / 2

            st1 = vitesses.bin_min
            st2 = vitesses.bin_max
            st1s = st1.astype("string")
            st2s = st2.astype("string")
            bin_label = st1s + '-' + st2s
            x_vit = bin_label
            y_vit = -vitesses.slope
            y_vit[y_vit <0]=np.NaN


            pp = PdfPages(
                str(path_to_figures) + '\ vitesse_graphique_allsize'  + '_' + str(
                    t) + ' Particle' + str(
                    cl) + '; ' + str(
                    tp_list_char.date[count_tp]) + ' to ' + str(tp_list_char.date[count_tp + 1]) + '.pdf')

            plt.clf()
            ax1 = plt.subplot(1, 1, count)
            plt.plot(x_vit, y_vit, 'b')
            plt.plot(x_vit, y_vit, 'ob', label='Speeds(m/d)', linestyle='--', markersize=3.5)
            plt.ylabel("Speed [m/d]")
            plt.xlabel("Size classes (mm)")
            plt.title("Speeds variation with size classes " + ' Particle' + str(cl) + str(
                tp_list_char.date[count_tp]) + ' to ' + str(tp_list_char.date[count_tp + 1]))
            plt.legend()
            # ax1.set_xticklabels(bin_label[0:4])
            plt.xticks(rotation=20)
            pp.savefig(plt.gcf(), bbox_inches='tight')  # Save each figure in the pdf
            pp.close()

            count_tp += 1

    path_to_fig=path_to_SPEED

    vitesses = pd.read_csv(
        str(path_to_fig) + '\Ecotaxa_sinking_speed_'  +  'ParticleAB' +
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
        str(path_to_fig) + '\Ecotaxa_sinking_speed_' + 'ParticleBV' +
        str(tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.csv')

    st1 = vitesses.bin_min
    st2 = vitesses.bin_max
    st1s = st1.astype("string")
    st2s = st2.astype("string")
    bin_label = st1s + '-' + st2s
    x_vit2 = bin_label
    y_vit2 = -vitesses.slope
    y_vit2[y_vit2 < 0] = np.NaN

    '''vitesses = pd.read_csv(
        str(path_to_fig) + ' 3' + '\Ecotaxa_sinking_speed_' + ' ' + ' Particle3' +
        str(tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.csv')

    st1 = vitesses.bin_min
    st2 = vitesses.bin_max
    st1s = st1.astype("string")
    st2s = st2.astype("string")
    bin_label = st1s + '-' + st2s
    x_vit3 = bin_label
    y_vit3 = -vitesses.slope

    vitesses = pd.read_csv(
        str(path_to_fig) + ' 4' + '\Ecotaxa_sinking_speed_' + ' ' + ' Particle4' +
        str(tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.csv')

    st1 = vitesses.bin_min
    st2 = vitesses.bin_max
    st1s = st1.astype("string")
    st2s = st2.astype("string")
    bin_label = st1s + '-' + st2s
    x_vit4 = bin_label
    y_vit4 = -vitesses.slope '''

    pp = PdfPages(str(path_to_fig) + '\ vitesse_graphique_allsize' + '_' + ' Particles ' + ' ' + str(
        tp_list_char.date[0]) + ' to ' + str(tp_list_char.date[1]) + '.pdf')

    plt.clf()
    ax1 = plt.subplot(1, 1, count)
    p1 = plt.plot(x_vit1, y_vit1, 'ob', label='ParticleAB Speeds(m/d)', linestyle='--', markersize=3.5)
    p2 = plt.plot(x_vit2, y_vit2, 'or', label='ParticleBV Speeds(m/d)', linestyle='--', markersize=3.5)
    '''p3 = plt.plot(x_vit3, y_vit3, 'og', label='Particle3 Speeds(m/d)', linestyle='--', markersize=3.5)
    p4 = plt.plot(x_vit4, y_vit4, 'ok', label='Particle4 Speeds(m/d)', linestyle='--', markersize=3.5) '''

    plt.ylabel("Speed [m/d]")
    plt.xlabel("Size classes (mm)")
    plt.title("Speeds variation with size classes " + ' Particles ' + str(tp_list_char.date[0]) + ' to ' + str(
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







