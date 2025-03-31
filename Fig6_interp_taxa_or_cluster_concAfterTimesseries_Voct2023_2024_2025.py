
#************************* READ FIRST interpolation....._V_Janv2024.py****#
#***********************************************************************#

# import needed libraries
import os
from pathlib import Path
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
# for data management
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
pd.options.display.max_columns = 100
import numpy as np
from matplotlib.legend import Legend

# for plots

import matplotlib.pyplot as plt
import datetime
from datetime import datetime
from datetime import timedelta
import warnings # these two lines to remove the annoying warnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)

# set up working space
path_to_data = Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6').expanduser()
os.chdir(str(path_to_data))

arg2 = 1 # log(n+1) figure for the interpolation or not
# Choose the definition either as rough or medium
definition_needed = ['angola_rough', 'angola_medium']
number = 0
def_chosen = definition_needed[number]
def_chosen

# import physical cluster
station_clusterphysic = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/station_cluster_kmeans_600.csv')
xiZ = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/xiZ.csv')
xii = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/xi.csv')

# import date, lat, lon of each profile
#profile = pd.read_csv("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/list_of_profiles_angola.csv")
profile = pd.read_csv(str(path_to_data) + '/list_of_profiles_angola.csv', sep=';', encoding='unicode_escape')
profile.rename(columns={'date':'date_ancien'}, inplace=True)
profile.rename(columns={'dateok':'date'}, inplace=True)
profile['CYCLE_NUMBER'] = (range(1, 119)) # create a list with station number (in my case I have 118 stations)


# There are two options here: you can import taxa concentrations or detritus cluster concentrations
# choose one option to create the full data frame
Q=[5,10,20,25,30,40,50,60,70,75, 80,90,95]
Q=[100]
for q in Q:
    # Open the corresponding table containing taxa concentrations
    #full = pd.read_csv(str(path_to_data) + '/angola_indiv_' + def_chosen + '.csv')
    # or open a table containing detritus cluster concentrations
    full = pd.read_csv(
        'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/clusters_concentrationsQq' + str(q) + '.csv')
    if q < 100:
        full = full.replace('cluster 1', 'cluster 11')
        full = full.replace('cluster 2', 'cluster 22')
        full = full.replace('cluster 3', 'cluster 33')
        full = full.replace('cluster 4', 'cluster 44')

        full = full.replace('cluster 11', 'cluster 4')
        full = full.replace('cluster 22', 'cluster 1')
        full = full.replace('cluster 33', 'cluster 3')
        full = full.replace('cluster 44', 'cluster 2')

    # merge
    full = pd.merge(full, profile, on="profile")
    # convert the date column to datetime objects
    AAA = round(full['date'], 0)
    # full['date'] = pd.to_datetime(full['date'], format="%Y%m%d")
    full['date'] = pd.to_datetime(AAA, format="%Y%m%d", errors='coerce')

    # add cluster physic
    full = pd.merge(full, station_clusterphysic, on="CYCLE_NUMBER")

    # define the base date
    ref_date = datetime(2021, 5, 5)  # METTRE LE DEBUT DE LA DATE DE MON FLOTTEUR

    # calculate the number of days between each date and the base date
    full['days'] = (full['date'] - ref_date).dt.days

    # Keep only 0-1000m
    full = full[full['depth'] < 1000]
    # choose the taxa that you want to interpolate
    # full = full[full['taxon'] == 'Copepoda']
    fulll = full
    fulll.rename(columns={'Cluster': 'taxon'}, inplace=True)

    path_to_figures = Path(
        'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/ClusterCenterQselection').expanduser()
    pp = PdfPages(str(path_to_figures) + '/interp clusters fev2024/' + str(t) + ' ARGO-ECO bbpPOCChla Ratio' +
                  str(q) + '.pdf')

    for cl in range(1, 5):
        # cl=2
        full = fulll[fulll['taxon'] == 'cluster ' + str(cl)]
        # full = full[full['taxon'] == 'conc']

        taxa_list = np.unique(full.taxon)
        taxa_list


        # plotting and interpolate functions
        def contour_levels_func(min_contour_level, max_contour_level, levels):
            """Function to define contour levels for contourf"""
            import numpy as np
            distance_levels = max_contour_level / levels
            contour_levels = np.arange(min_contour_level, max_contour_level, distance_levels)
            return contour_levels


        def gridding_func(pos_min_max, depth_min_max, pos_array, depth, param):
            import numpy as np
            from scipy.interpolate import griddata
            grid_method = "linear"  # choose the gridding method here
            # method to do the regridding, can be nearest (for nearest neighbour) or linear

            xi = np.linspace(min(pos_min_max), max(pos_min_max), 1000)
            yi = np.linspace(min(depth_min_max), max(depth_min_max), 200)
            zi = griddata((pos_array, depth), param, (xi[None, :], yi[:, None]), method=grid_method)
            return xi, yi, zi


        def nonlinear_colormap():
            import pylab as pyl
            # import numpy as np
            levels1 = [0, 1, 2]


        from matplotlib import ticker

        cm = nonlinear_colormap()

        log10_or_not = arg2  # if 0 then not log10 ; if 1 then log10
        if arg2 == 1:
            correction = "log"
        else:
            correction = "normal"

        # sort the groups
        taxa_list = sorted(taxa_list)

        # Run the entire for loop to make the cotour plot
        for t in taxa_list:
            df = full[full.taxon == t].reset_index(drop=True)

            # put the concentration in log10 if arg2 is 1
            if arg2 == 1:
                df['conc_log'] = df["conc"].apply(lambda x: np.log1p(x))
            else:
                df = df

            # put the date in float. It doesn't work with dates in int
            df.days = df.days.apply(lambda x: float(x))

            dat = np.array(df['days'])
            '''date_labbaty = [datetime(2021, 5, 5) + timedelta(days=int(idt)) for idt in dat]  # changer ma date dedans
            # Format each date in the list as a string in the format "YYYY-MM-DD"
            formatted_datesbaty = [d.strftime("%Y-%m-%d") for d in date_labbaty]
            df_temp_baty = df
            df_temp_baty['datbaty'] = date_labbaty
            df_temp_batyshort = df_temp_baty.groupby(['profile', 'datbaty']).agg({'depth': 'max'}).reset_index()
            # df_temp_baty['Maxbaty']=df_temp_baty['DEPTH']
            # df_temp_baty.drop('DEPTH', inplace=True, axis=1)
            df_temp_batyh = pd.merge(df_temp_baty, df_temp_batyshort, how="left", on=['profile'])
            AA, AAindices = np.unique(formatted_datesbaty, return_index=True)
            maxz = df_temp_batyh['depth_y'][AAindices]
            xiZ = pd.DataFrame(xi)
            count_ia = 0
            for ia in AA:
                xiZ[pd.DataFrame(formatted_dates) == ia] = maxz[maxz.index[count_ia]]
                count_ia += 1
            xiZ[xiZ <= 355] = np.nan
            '''
            pos_array = np.array(df['lon'])
            unique_pos = np.unique(pos_array)

            unique_pos_dat = np.unique(dat)
            depth = np.array(df['depth'] * (-1))

            fig = plt.figure(figsize=(19.20, 9.83))  # (5, 3.5)
            plt.rcParams.update({'font.size': 13.5})
            ax1 = fig.add_axes([0, 0, 0, 0])
            count = 1

            # limits of concentration for the figure
            if arg2 == 1:
                parameter_dict = {"conc": [min(df.conc_log), max(df.conc_log)]}
                legend_dict = {"conc": "log Concentration [#.$m^{-3}$]"}
                legend_dict = {"conc": "Concentration [#.$m^{-3}$]"} #, log$_{e}$ scale
            else:
                parameter_dict = {"conc": [0, (max(df.conc) / 6)]}
                legend_dict = {"conc": "Concentration [#.$m^{-3}$]"}

            for p in parameter_dict.keys():
                #fig.tight_layout()
                ax1 = plt.subplot(1, 1, count)
                #plt.subplots_adjust(hspace=0.2, wspace=0.2, top=0.92, bottom=0.15, left=0.125, right=0.9)
                plt.subplots_adjust(hspace=0.2, wspace=0.2, top=0.88, bottom=0.19, left=0.125, right=0.9)
                param = np.array(df[p])  # select the concentration

                # choose the limits in depth and time
                depth_min_max = [0, -1000]
                pos_min_max = [min(df.days), max(df.days)]  # Date transect

                # interpolate
                xi, yi, zi = gridding_func(pos_min_max, depth_min_max, dat[~np.isnan(dat)], depth[~np.isnan(dat)],
                                           param[~np.isnan(dat)])
                # xi = [datetime.datetime(2021, 5, 5) + datetime.timedelta(days=i) for i in xi]
                parameter = p

                # set the parameters of contour levels
                levels = 15
                min_contour_level = parameter_dict[parameter][0]
                max_contour_level = parameter_dict[parameter][1]
                contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)

                # Convert each element of xi to a datetime object by adding the number of days to the start date
                # date_list = [ref_date + timedelta(days=x) for x in xi]
                date_list = [datetime(2021, 5, 5) + timedelta(days=i) for i in xi]  # changer ma date dedans
                # Format each date in the list as a string in the format "YYYY-MM-DD"
                formatted_dates = [d.strftime("%Y-%m-%d") for d in date_list]
                formatted_datesym = [d.strftime("%Y-%m") for d in date_list]

                # plot it
                if arg2 == 1:
                    p1 = plt.contour(formatted_dates, yid, zid, levels=[0, 1025.7], colors='black',linewidths=2.5)
                    h1,_=p1.legend_elements()
                    #plt.legend(h1[1],['Isopycnal 1025.7 Kg/$m^3$'], loc="lower left")
                    p2 = plt.plot(MLD, color='hotpink',  linewidth=2.5,)#label='Mixed Depth Layer', linestyle='dashed',
                    h2=p2[:2]
                    p3 = plt.contourf(formatted_dates, yi, np.log1p(zi), contour_levels,
                                      cmap='viridis', alpha=1, extend="both")


                else:
                    p1 = plt.contour(formatted_dates, yid, zid, levels=[0, 1025.7], colors='black',linewidths=2.5)
                    p2 = plt.plot(MLD, color='hotpink', linewidth=2.5,
                                  label='Mixed Layer Depth') # , linestyle='dashed'

                    p3 = plt.contourf(formatted_dates, yi, zi, contour_levels,
                                      cmap='viridis', alpha=1, extend="both")
                    plt.legend(loc="lower left",fontsize=15)


                ax1.fill_between(np.array(xii[np.isfinite(xiZ.iloc[:,1])].iloc[:,1]), np.array(-xiZ[np.isfinite(xiZ.iloc[:,1])].iloc[:,1]),
                                 -1000, alpha=1, color='white')
                high = -1000
                ha = 15
                ha1 = 2.5 #3.5
                '''ax1.fill_between([formatted_dates[214], formatted_dates[289]], [0, 0], [high, high], alpha=0.2,color='cyan')
                ax1.fill_between([formatted_dates[299], formatted_dates[389]], [0, 0], [high, high], alpha=0.2,color='cyan')
                ax1.fill_between([formatted_dates[459], formatted_dates[534]], [0, 0], [high, high], alpha=0.2,color='cyan')
                ax1.fill_between([formatted_dates[561], formatted_dates[633]], [0, 0], [high, high], alpha=0.2,color='cyan')
                ax1.fill_between([formatted_dates[760], formatted_dates[835]], [0, 0], [high, high], alpha=0.2,color='cyan')
                ax1.fill_between([formatted_dates[884], formatted_dates[959]], [0, 0], [high, high], alpha=0.2,color='cyan')
                '''

                r = patches.Rectangle((xi[206], 0 + ha1), xi[301] - xi[206], high - ha, lw=2, ls='--',edgecolor='#ED7D31', facecolor='none', clip_on=False)
                ax1.add_patch(r)
                r = patches.Rectangle((xi[307], 0 + ha1), xi[402] - xi[307], high - ha, lw=2, ls='--',edgecolor='#ED7D31', facecolor='none', clip_on=False)
                ax1.add_patch(r)
                r = patches.Rectangle((xi[439], 0 + ha1), xi[543] - xi[439], high - ha, lw=2, ls='--',edgecolor='#ED7D31', facecolor='none', clip_on=False)
                ax1.add_patch(r)
                r = patches.Rectangle((xi[569], 0 + ha1), xi[678] - xi[569], high - ha, lw=2, ls='--',edgecolor='#ED7D31', facecolor='none', clip_on=False)
                ax1.add_patch(r)
                r = patches.Rectangle((xi[755], 0 + ha1), xi[869] - xi[755], high - ha, lw=2, ls='--',edgecolor='#ED7D31', facecolor='none', clip_on=False)
                ax1.add_patch(r)
                r = patches.Rectangle((xi[884], 0 + ha1), xi[970] - xi[884], high - ha, lw=2, ls='--',edgecolor='#ED7D31', facecolor='none', clip_on=False)
                ax1.add_patch(r)

                '''ax1.legend([h1[1]], ['iso'], loc="lower left")
                leg2 = Legend(ax1, h2, ['mld'], loc='lower right')
                ax1.add_artist(leg2)'''

                b1, = ax1.plot([], marker="", markersize=15, linestyle="-", color='hotpink', label="Mixed layer Depth")
                b2, = ax1.plot([], marker="", markersize=15, linestyle="-", color="black", label="Isopycnal 1025.7 Kg/$m^3$")
                b3, = ax1.plot([], marker="", markersize=15, linestyle="--", color='#ED7D31', label="Export Event")
                ax1.legend(handles=[b1, b2,b3])
                handles, labels = ax1.get_legend_handles_labels()
                ax1.legend(handles, labels, loc='lower left',fontsize=20)
                HIG=-2.7#HIG=500
                ssi=25
                plt.text(xi[250], -HIG, '1', fontsize = ssi, color='#ED7D31')
                ax1.text(xi[350], -HIG, '2', fontsize = ssi, color='#ED7D31')
                ax1.text(xi[490], -HIG, '3', fontsize = ssi, color='#ED7D31')
                ax1.text(xi[595], -HIG, '4', fontsize = ssi, color='#ED7D31')
                ax1.text(xi[800], -HIG, '5', fontsize = ssi, color='#ED7D31')
                ax1.text(xi[922], -HIG, '6', fontsize = ssi, color='#ED7D31')

                alphabet = ['','A', 'B', 'C', 'D']
                ax1.text(-55, 4, alphabet[cl], fontsize=25, color='black', weight='bold')

                # plotting parameters
                ax1.hlines(y=300, xmin=min(formatted_dates), xmax=max(formatted_dates), colors='r')
                plt.ylabel("Depth [m]", fontsize=22) # 21
                plt.xlim(min(formatted_dates), max(formatted_dates))
                plt.ylim(-1000, 0)
                plt.xlabel("Date", fontsize=22) #21
                plt.xticks(rotation=25)
                #ax1.tick_params(axis='y', labelsize=13.5)
                #ax1.tick_params(axis='x', labelsize=13.5)
                plt.title(t, size=12)

                # Set the x-axis tick frequency to every one month
                #ax1.set_xticks(formatted_dates[::60])
                #ax1.set_xticklabels(formatted_datesym[::60])
                formatted_datesq = []
                for ii in [0, 88, 172, 259, 347, 431, 518, 603, 690, 777, 856, 943]:
                    formatted_datesq.append(formatted_dates[ii])

                ax1.set_xticks(formatted_datesq)

                ax1.tick_params(axis='x', labelsize=22) # 18
                ax1.tick_params(axis='y', labelsize=22)
                cb = plt.colorbar(p3, orientation='vertical')
                cb.ax.set_ylabel(legend_dict[parameter], fontsize=22)
                cb.ax.tick_params(labelsize=18)

                tick_locator = ticker.MaxNLocator(nbins=5)
                cb.locator = tick_locator

                ticks = [float(t.get_text().replace('−', '-')) for t in cb.ax.get_yticklabels()]
                expticks = np.exp(ticks) - 1
                cb.ax.set_yticklabels(np.round(expticks, 1))  # vertically oriented colorbar

                count += 1
                #plt.title('cluster' + str(cl))
                strclust=['','Flakes','Agglomerates','Strings','Spheres']
                plt.title(strclust[cl],fontsize=25, weight='bold', x=0.5,y=1.06)

                # Show the plot
                plt.show()

                '''path_to_figures = Path(
                    'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/ClusterCenterQselection').expanduser()
                plt.savefig(str(path_to_figures) + '/interp clusters fev2024/cluster' + str(cl) + 'Q' + str(q) + 'BOXFEV2024.png', dpi=1200)
                pp.savefig(fig,dpi=1200)  # Save each figure in the pdf # plt.gcf()
                '''

                path_to_figures2024 = Path(
                    'C:/Users/syawo/Downloads/Figures_pdf_Angola_Biogeosciences_submission2025/Interp_Cluster').expanduser()
                plt.savefig(str(path_to_figures2024) + '/' + str(cl) + 'Q' + str(q) + 'April2025.png', dpi=300)
                pp.savefig(fig, dpi=1200)  # Save each figure in the pdf #plt.gcf()

                plt.close('all')
    pp.close()

