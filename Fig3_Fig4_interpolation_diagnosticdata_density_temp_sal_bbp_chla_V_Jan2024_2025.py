#VJan2024
# # # 1. Setup
# import needed libraries
import os
from pathlib import Path
import matplotlib.patches as patches


# for data management
import pandas as pd
import warnings  # these two lines to remove the annoying warnings from pandas

warnings.simplefilter(action='ignore', category=FutureWarning)
pd.set_option('display.max_rows', 500, 'display.max_columns', 500)  # to be able to see a whole table
import numpy as np
import datetime
from datetime import date
from datetime import datetime
from datetime import timedelta
# Settling speed computations
from scipy.optimize import curve_fit
import scipy.stats
from scipy import interpolate

# for figures
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors, cm

# for interpolation
from plotting_funcs import contour_levels_func, gridding_func
from matplotlib import ticker

'''import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(3, 2)'''



# set up working space
path_to_data = Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/INTERP_ECOPART/').expanduser()
os.chdir(str(path_to_data))

# import TEMP SAL POC
ecopartf = pd.read_table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/Ecopart_diagnostics_data_647.tsv", sep='\t', encoding='unicode_escape')

ecopartf = ecopartf[ecopartf['Profile'].str.contains('a')]
mask = ecopartf["Chlorophyll-a [mg/m3]"] < 0
colonne = "Chlorophyll-a [mg/m3]"
ecopartf.loc[mask, colonne] = 0
maskb = ecopartf["bbp POC [mgC/m3]"] < 0
colonneb = "bbp POC [mgC/m3]"
ecopartf.loc[maskb, colonneb] = 0
ecopartf['colonne'] = ecopartf["bbp POC [mgC/m3]"] / ecopartf["Chlorophyll-a [mg/m3]"]
#ecopartf.colonne[ecopartf.colonne.isin([-np.inf, np.inf])] = np.nan
ecopartff = ecopartf
'''ecopartf = ecopartff.iloc[:,[1, 4, 5,38,39,42,41,44,36]] # data
ecopartf = ecopartf.iloc[:,[1, 4, 5, 7, 11, 12, 17, 18, 37,38,39,40,42,41,44,36]] # data # date,PRESSION,flux,fluxmip,fluxmap,mip,map,density,tempera,salinity,Doxy,bbppoc,Chla,bbp/CHLA,depth
ecopartf = ecopartf.iloc[:,[1, 4, 5, 7, 11, 12, 17, 18, 37,38,39,40,42,41,44,36]] # data # date,PRESSION,flux,fluxmip,fluxmap,mip,map,density,tempera,salinity,Doxy,bbppoc,Chla,bbp/CHLA,depth
'''


#0/0

ecopartf = ecopartff.iloc[:, [1, 4, 5, 37, 38, 39, 40, 42, 41, 11,12, 44,36]]  # data # date,PRESSION,flux,fluxmip,fluxmap,mip,map,density,tempera,salinity,Doxy,bbppoc,Chla,bbp/CHLA,depth
ecopartf.rename(columns={'colonne': 'bbpPOC/Chla'}, inplace=True)
ecopartf.rename(columns={'TEMPERATURE': 'Temperature'}, inplace=True)
ecopartf.rename(columns={'SALINITY': 'Salinity'}, inplace=True)
ecopartf.rename(columns={'Chla': 'Chlorophyll a'}, inplace=True)
ecopartf.rename(columns={'bbp POC': 'Particulate Organic Carbon'}, inplace=True)

ecopartf['DATE'] = ecopartf['Date_Time']
ecopartf['DEPTH'] = ecopartf['Depth [m]']
taxa_list = ['']
ecopart_binsf = ['Density', 'TEMPERATURE', 'SALINITY', 'Dissolved Oxygen', 'bbp POC', 'Chla','Flux Mip','Flux Map','']  # ,'bbpPOC/Chla'  bins
ecopart_binsf = ['Density', 'Temperature', 'Salinity', 'Dissolved Oxygen', 'Particulate Organic Carbon', 'Chlorophyll a','Flux Mip','Flux Map','']  # ,'bbpPOC/Chla'  bins

# I used already saved vector (array) from interopalation see below
xii = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/xi.csv')
xii=xii['0']
# bin the data per 100m
def rounditup(x, precision, method):
    import math
    if method == "floor":
        return math.floor(x / precision) * precision
    elif method == "ceiling":
        return math.ceil(x / precision) * precision
    else:
        return "give the parameter floor or ceiling"


tp_list = [20211001, 20211113]  # 2
tp_list = [20210713, 20210819]  # 0
tp_list = [20210819, 20210925]  # 1
tp_list = [20211113, 20211215]  # 3
tp_list = [20211113, 20220119]  # 3b
tp_list = [20220120, 20220306]  # 4
tp_list = [20220306, 20220426]  # 5

path_to_SPEED = Path(
    'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/INTERP_ECOPART/').expanduser()

tp_listavf = [[20210720, 20210815], [20210819, 20210920], [20211015, 20211110], [20211120, 20211215],
              [20220130, 20220225], [20220315, 20220410]]  # 0
tp_lista = [[20210713, 20210819], [20210819, 20210925], [20211001, 20211113], [20211113, 20211215],
            [20211113, 20220119], [20220120, 20220306], [20220306, 20220426]]  # 0
tp_lista = [[20210505, 20220426]]  # 0
PartABV = ['f']
clist = ["Potential density [kg/$m^3$]", "Temperature [°C]", "Salinity", "Dissolved Oxygen [micromol/kg] ",
         "Bbp POC[mgC/$m^3$]", "Chla[mg/$m^3$]", "Flux Mip [mgC/$m^2$/d]", "Flux Map [mgC/$m^2$/d]", "Bbp POC[mgC/$m^3$]", "bbp POC/Chla [mgC/mg]"]

clist = ["Potential density [kg/$m^3$]", "Temperature [°C]", "Salinity", "Dissolved Oxygen [micromol/kg] ",
         "Bbp POC[mgC/$m^3$]", "Chla [mg/$m^3$]", "Flux Mip [mgC/$m^2$/d]", "Flux Map [mgC/$m^2$/d]", "POC[mgC/$m^3$]", "bbp POC/Chla [mgC/mg]"]

for tp_list in tp_lista:

    tp_list_char = pd.to_datetime(tp_list, format="%Y%m%d", errors='coerce')  #

    for cl in PartABV:

        ecopartchoose = eval('ecopart' + cl)

        ecopartchoose.loc[:, 'depth_part'] = 99999
        ecopartchoose['depth_part'] = ecopartchoose["DEPTH"].apply(
            lambda x: rounditup(x, 100, "floor") + 50 if x < 600 else rounditup(x, 200, "floor") + 100)


        path_to_figures = Path(
            'C:/Users/syawo/Downloads/Figures_pdf_Angola_Biogeosciences_submission2025/Interp_PhysqiueEnv').expanduser()

        arg1 = 0  # For the defitition choice of taxo definition : either 0 for EXPORTS_rough or 1 for EXPORTS_medium
        argb = 0

        count_tp = 0
        count = 1
        fig = plt.figure(figsize=(19.20,9.83))  # (19.20,9.83) # (5, 3.5)
        plt.rcParams.update({'font.size': 13.5})
        plt.rcParams['contour.linewidth'] = 1.8
        ax1 = fig.add_axes([0, 0, 0, 0])

        for tp in tp_list[0:-1]:
            ecopartchoose['DATE'] = pd.to_datetime(ecopartchoose['DATE'])  # .loc[row_indexer,col_indexer] = value instead
            detritus = ecopartchoose[ecopartchoose['DATE'] > tp_list_char[count_tp]].reset_index(drop=True)
            detritus = detritus[detritus['DATE'] < tp_list_char[count_tp + 1]].reset_index(drop=True)
            detritus = detritus[detritus['DEPTH'] < 1000].reset_index(drop=True)

            detritus.shape

            # Now compute number and volumes of particles per profile, size class, taxo group and depth bin
            df = detritus

            ref_date = datetime(2021, 5, 5)

            df.head()

            # Plot the distribution according to the depth and time

            for t in taxa_list:
                df_taxa = df
                if cl == 'f':
                    bin_min_list = ecopart_binsf
                else:
                    bin_min_list = ecopart_binsf

                depth_part_list = np.unique(df_taxa.depth_part)

                # Open a pdf to save the figures rough definition
                pp = PdfPages(str(path_to_figures) + '/' + str(t) + ' ARGO-ECO ENV Aout2024' + str(tp_list_char.date[count_tp]) + ' to ' + str(
                    tp_list_char.date[count_tp + 1]) + '.pdf')
                # Add a counter for the b_max selection
                count_b = 0
                for b in bin_min_list[0:-1]:

                    if cl == 'f':
                        bmax = ''
                    else:
                        bmax = bin_min_list[count_b + 1]

                    countb = count_b + 3

                    day0 = -1
                    # b = bin_min_list[1]
                    listsize = list(df_taxa.columns)
                    df_temp = df_taxa  # .iloc[:,countb]*0.001
                    if cl == 'f':
                        if listsize[countb] == listsize[3]:
                            df_temp[listsize[countb]] = df_taxa[listsize[countb]]  #
                        else:
                            df_temp[listsize[countb]] = df_taxa[listsize[countb]]  #

                    else:
                        df_temp[listsize[countb]] = df_taxa[listsize[countb]]  #

                    depth_part_list = np.unique(df_taxa.depth_part)

                    # Get the data for the size b and taxa_t in concentration and from the previous computations
                    AAA = df_temp['DATE']
                    df_temp['date'] = pd.to_datetime(AAA, format="%Y%m%d", errors='coerce')
                    df_temp['days'] = (df_temp['date'] - ref_date).dt.days
                    # select the values for interpolation
                    iday = df_temp.days
                    idays = iday.reset_index(drop=True)
                    date_min = np.min(idays)
                    date_max = np.max(idays)

                    # 3. Interpolate the concentration field
                    # dat = np.array(df_temp['date'])
                    dat = np.array(df_temp['days'])
                    date_labbaty = [datetime(2021, 5, 5) + timedelta(days=int(idt)) for idt in
                                    dat]  # changer ma date dedans
                    # Format each date in the list as a string in the format "YYYY-MM-DD"
                    formatted_datesbaty = [d.strftime("%Y-%m-%d") for d in date_labbaty]
                    df_temp_baty = df_temp
                    df_temp_baty['datbaty'] = date_labbaty
                    df_temp_batyshort = df_temp_baty.groupby(['Profile', 'datbaty']).agg({'DEPTH': 'max'}).reset_index()
                    # df_temp_baty['Maxbaty']=df_temp_baty['DEPTH']
                    # df_temp_baty.drop('DEPTH', inplace=True, axis=1)
                    df_temp_batyh = pd.merge(df_temp_baty, df_temp_batyshort, how="left", on=['Profile'])
                    AA, AAindices = np.unique(formatted_datesbaty, return_index=True)
                    maxz = df_temp_batyh['DEPTH_y'][AAindices]
                    '''xiZ=pd.DataFrame(xi)
                    count_ia= 0
                    for ia in AA:
                        xiZ[pd.DataFrame(formatted_dates) == ia] = maxz[maxz.index[count_ia]]
                        count_ia += 1 
                    xiZ[xiZ<=355]=np.nan

                    '''
                    '''count_di=0
                    for di in df_temp_batyshort['datbaty']:
                        df_temp_baty['Maxbaty'][df_temp_baty['datbaty']== di]= df_temp_batyshort['DEPTH'][count_di]
                        count_di += 1'''

                    unique_pos_dat = np.unique(dat)
                    depth = np.array(df_temp['DEPTH'] * (-1))
                    param = np.array(df_temp[listsize[countb]])

                    # choose the limits in depth and time
                    depth_min_max = [0, -1000]
                    if b == 'bbp POC [0-1000m]':
                        depth_min_max = [0, -1000]
                    if b == 'Chlorophyll a':
                        depth_min_max = [0, -1000]

                    # pos_min_max = [min(df_temp.date), max(df_temp.date)]  # Date transect
                    pos_min_max = [min(df_temp.days), max(df_temp.days)]  # Date transect
                    # interpolate *********************************************************************
                          #****************************************************************************
                    xi, yi, zi = gridding_func(pos_min_max, depth_min_max, dat, depth, param)


                    '''# check anf fix nan for the batymetry
                    Zibaty=np.empty((200,1000))
                    YIbaty=np.empty((1000,1))
                    count_dii=0
                    for dii in formatted_dates:
                        Maxbaty=np.max(df_temp_batyshort['DEPTH'][df_temp_batyshort['datbaty']== dii])
                        #zi[yi>Maxbaty,count_dii]=np.nan
                        Zibaty[yi>= -Maxbaty,count_dii]=np.nan
                        Zibaty[yi < -Maxbaty, count_dii] = -Maxbaty
                        YIbaty[count_dii]=-Maxbaty
                        count_dii += 1'''

                    levels = 15

                    if b == 'Density':
                        levels = 1017
                        zi[161:170, 402:439] = np.mean(zi[161:170, 430:439])
                        # zi[164, 402:439] = np.mean(zi[145:159, 409:439])

                        # MLD COMPUTATION
                        z197=zi[197,] # THE SURFACE DENSITY AT 10m
                        zidiff = np.abs(zi - z197) # THE DIFFERENCE
                        mldKEEP=[]
                        for ie in range (0,1000): # LOOP FOR EACH PROFIL
                            zimaskmld = np.where(zidiff[:, ie] > 0.03)  # THE CUTOFF
                            zimaskmld10=zimaskmld[0][zimaskmld[0]<197]
                            #ziy = zimaskmld[0][len(zimaskmld[0]) - 1]  # THE CLOSEST POSITION TO THE SURFACE
                            ziy = zimaskmld10[len(zimaskmld10) - 1]  # THE CLOSEST POSITION TO THE SURFACE
                            fx=interpolate.interp1d(zidiff[[ziy, ziy+1],ie],yi[[ziy,ziy+1]])
                            mldKEEP.append(fx(0.03))
                            #mldKEEP.append(yi[ziy])  # DEPTH MLD
                            MLD=pd.DataFrame(mldKEEP)
                        yi[zimaskmld]

                    if b == 'Chlorophyll a':
                        levels = 10
                        # 0/0
                        dcm = np.nanmax(zi, axis=0)
                        dcmlg = zi == dcm
                        dcmdepth = [np.mean(yi[dcmlg[:, ii]]) for ii in range(0, 1000)]
                        dcmok = -np.asarray(dcmdepth)
                        '''dcm=zi==np.nanmax(zi,axis=0)
                        x = [1, 2, 3]
                        a = ones(3, 1) * x
                        YYI=yi*np.full((200, 1000), 1, dtype=int)
                        DCM=np.mean(ziZ(dcm))'''

                    if b == 'Salinity':
                        # 0/0
                        zi[161:170, 402:439] = np.mean(zi[161:170, 430:439])

                        # B=zi[153:176,350:450]
                        # B[B<35]=np.nan
                        # zi[153:176,350:450]=B
                        min_contour_level = 32
                        max_contour_level = max(df_temp[listsize[countb]])
                        max_contour_level = 37.5
                        levels = 33


                    elif b == 'Temperature':
                        min_contour_level = 10
                        max_contour_level = max(df_temp[listsize[countb]])
                        levels = 10
                        z198 = zi[198,]  # THE SURFACE TEMPERATURE AT 5m
                        zidiffT = np.abs(zi - z198)  # THE DIFFERENCE
                        mldKEEPT = []
                        for ie in range(0, 1000):  # LOOP FOR EACH PROFIL
                            zimaskmldT = np.where(zidiffT[:, ie] > 0.2)  # THE CUTOFF
                            zimaskmld5T = zimaskmldT[0][zimaskmldT[0] < 198]
                            # ziy = zimaskmld[0][len(zimaskmld[0]) - 1]  # THE CLOSEST POSITION TO THE SURFACE
                            ziyT = zimaskmld5T[len(zimaskmld5T) - 1]  # THE CLOSEST POSITION TO THE SURFACE
                            fx = interpolate.interp1d(zidiffT[[ziyT, ziyT + 1], ie], yi[[ziyT, ziyT + 1]])
                            mldKEEPT.append(fx(0.2))
                            #mldKEEPT.append(yi[ziyT])  # DEPTH MLD
                            MLDT = pd.DataFrame(mldKEEPT)
                    elif b == 'bbpPOC/Chla':
                        df_temp['conc_log'] = df_temp[listsize[countb]].apply(lambda x: np.log1p(x))
                        df_temp['conc_log'] = df_temp[listsize[countb]].apply(lambda x: np.log10(x))
                        tour = [min(df_temp.conc_log[~df_temp.conc_log.isin([-np.inf])]), max(df_temp.conc_log)]
                        min_contour_level = tour[0]
                        max_contour_level = tour[1] / 3 + tour[0] / 3
                        # max_contour_level = tour[1]
                        levels = 15
                    elif b == 'Particulate Organic Carbon' : # 'bbp POC': #
                        df_temp['conc_log'] = df_temp[listsize[countb]].apply(lambda x: np.log1p(x))
                        # df_temp['conc_log'] = df_temp[listsize[countb]].apply(lambda x: np.log10(x))
                        tour = [min(df_temp.conc_log[~df_temp.conc_log.isin([-np.inf])]), max(df_temp.conc_log)]
                        min_contour_level = tour[0]
                        max_contour_level = tour[1] / 3 + tour[0] / 3
                        # max_contour_level = tour[1]
                        levels = 10

                    elif b == 'Chlorophyll a':
                        df_temp['conc_log'] = df_temp[listsize[countb]].apply(lambda x: np.log1p(x))
                        # df_temp['conc_log'] = df_temp[listsize[countb]].apply(lambda x: np.log10(x))
                        tour = [min(df_temp.conc_log[~df_temp.conc_log.isin([-np.inf])]), max(df_temp.conc_log)]
                        min_contour_level = tour[0]
                        max_contour_level = tour[1] / 3 + tour[0] / 3
                        # max_contour_level = tour[1]
                        levels = 100
                    elif b == 'Flux Mip':
                        levels = 100
                    elif b == 'Flux Map':
                        levels = 100
                    else:
                        min_contour_level = min(df_temp[listsize[countb]])
                        max_contour_level = max(df_temp[listsize[countb]])

                    contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)

                    date_lab = [datetime(2021, 5, 5) + timedelta(days=i) for i in xi]  # changer ma date dedans
                    # Format each date in the list as a string in the format "YYYY-MM-DD"
                    formatted_dates = [d.strftime("%Y-%m-%d") for d in date_lab]
                    formatted_datesym = [d.strftime("%Y-%m") for d in date_lab]
                    a_lab = [datetime(2021, 5, 5) + timedelta(days=i) for i in
                             df_temp.days]  # changer ma date dedans
                    # Format each date in the list as a string in the format "YYYY-MM-DD"
                    formatted_dates_arr = [d.strftime("%Y-%m-%d") for d in a_lab]

                    # for batymetry plot
                    xiZ = pd.DataFrame(xii)
                    count_ia = 0
                    for ia in AA:
                        loo = np.where((pd.DataFrame(formatted_dates) == ia) == True)
                        xiZ['0'][loo[0]] = maxz[maxz.index[count_ia]]
                        count_ia += 1
                    xiZ[xiZ <= 355] = np.nan

                    # 4. Plot
                    #
                    #
                    #

                    plt.clf()
                    ax1 = plt.subplot(1,1,count)
                    '''if count_b == 0:
                        ax1 = plt.subplot(gs[0, 0])
                    if count_b == 1:
                        ax2 = plt.subplot(gs[0, 1])
                    if count_b == 2:
                        ax3 = plt.subplot(gs[1, 0])
                    if count_b == 3:
                        ax4 = plt.subplot(gs[1, 1])
                    if count_b == 4:
                        ax5 = plt.subplot(gs[2, 0])
                    if count_b == 5:
                        ax6 = plt.subplot(gs[2, 1])'''


                    #plt.subplots_adjust(hspace=0.2, wspace=0.2, top=0.92, bottom=0.15, left=0.125, right=0.9)
                    plt.subplots_adjust(hspace=0.2, wspace=0.2, top=0.885   , bottom=0.19, left=0.125, right=0.9)

                    # Concentration field
                    # p1 = plt.contourf(xi, yi, zi, contour_levels, cmap=cm.viridis, alpha=1, extend="both",norm=colors.LogNorm(vmin=np.nanmin(zi), vmax=np.nanmax(zi) / 5))


                    if b == 'Density':#contour_levels cmap="RdBu_r"
                        zid=zi
                        yid=yi
                        p1 = plt.contour(formatted_dates, yi, zi, levels=[0, 1025.7],colors='black',linewidths=2.5)
                        h1, _ = p1.legend_elements()
                        #plt.legend([h1[1]], ['isopycnal 1025.7 Kg/$m^3$'])
                        #plt.legend(p1, ["","isopycnal 1025.7 Kg/$m^3$"])
                        #plt.clabel(p1, inline=1, fontsize=10)
                        #labels = ['', 'Isopycnal 1025.7 Kg/m3']
                        p1 = plt.plot(MLD,color='hotpink',linewidth=2.5)#,linestyle='dashed'
                        #p1 = plt.plot(MLDT, color='orange', linestyle='dashed', linewidth=1.5,label='Mixed Layer Depth from temperature at 5m')
                        p1 = plt.contourf(formatted_dates, yi, zi,levels=[1021,1021.5,1022,1022.5,1023,1023.5,1024,1024.5,1025,1025.5,1026,1026.5,1027,1027.5,1028,1028.5] ,cmap='viridis' , alpha=1,
                                          extend="both", vmin=1022, vmax=1027)

                        #plt.legend(loc="lower left")

                    elif b == 'Salinity':#
                        p1 = plt.contour(formatted_dates, yid, zid, levels=[0, 1025.7], colors='black',linewidths=2.5)
                        p1 = plt.plot(MLD, color='hotpink', linewidth=2.5)#, linestyle='dashed'

                        p1 = plt.contourf(formatted_dates, yi, zi,levels=[34,34.1,34.25,34.4,34.5,34.6,34.75,34.9,35,35.10,35.25,35.4, 35.5,35.6, 35.75,35.9, 36,36.1, 36.25,36.4, 36.5,36.6, 36.75,36.9, 37] , cmap='viridis' , alpha=1,
                                          extend="both", vmin=34, vmax=36.3)

                        #plt.legend(loc="lower left")

                    elif b == 'Temperature':

                        #
                        p1 = plt.contour(formatted_dates, yid, zid, levels=[0, 1025.7], colors='black',linewidths=2.5)
                        #p1 = plt.plot(MLDT, color='green', linestyle='dashed', linewidth=1.5, label='Mixed Layer Depth from temperature at 5m')
                        p1 = plt.plot(MLD, color='hotpink', linewidth=2.5)#, linestyle='dashed'

                        #plt.legend(loc="lower left")

                        p1 = plt.contourf(formatted_dates, yi, zi,
                                          levels=[8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28],
                                          cmap='viridis', alpha=1,extend="both", vmin=8, vmax=28)

                        #plt.legend(loc="lower left")
                        # 0/0
                    elif b == 'Particulate Organic Carbon':
                        p1 = plt.contour(formatted_dates, yid, zid, levels=[0, 1025.7], colors='black',linewidths=2.5)
                        p1 = plt.plot(MLD, color='hotpink', linewidth=2.5)#, linestyle='dashed'
                        p1 = plt.contourf(formatted_dates, yi, np.log1p(zi), levels=[0, 0.20, 0.40, 0.60, 0.80, 1.0, 1.10, 1.20, 1.40, 1.60, 1.80, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4,3.6], cmap='viridis', alpha=1,
                                          extend="both")
                        #plt.legend(loc="lower left")

                    elif b == 'Chlorophyll a':
                        p1 = plt.contour(formatted_dates, yid, zid, levels=[0, 1025.7], colors='black',linewidths=2.5)
                        p1 = plt.plot(MLD, color='hotpink', linewidth=2.5)#, linestyle='dashed'

                        p1 = plt.contourf(formatted_dates, yi, np.log1p(zi), levels=[0,0.05, 0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1], cmap='viridis', alpha=1, extend="both", vmin=0, vmax=1)
                        #plt.legend(loc="lower left")
                    elif b == 'Dissolved Oxygen':  # contour_levels
                        p1 = plt.contour(formatted_dates, yid, zid, levels=[0, 1025.7], colors='black',linewidths=2.5)
                        p1 = plt.plot(MLD, color='hotpink', linewidth=2.5) #, linestyle='dashed'
                        #cmap = 'RdBu_r'
                        p1 = plt.contourf(formatted_dates, yi, zi,
                                          levels=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,110, 120, 130, 140, 150, 160, 170, 180, 190, 200,210, 220, 230, 240, 250, 260],
                                          cmap='viridis', alpha=1, extend="both", vmin=0, vmax=260)
                        '''[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105,
                                                  110, 115, 120, 125, 130,135, 140, 145,150, 155,160,165, 170, 175,180, 185,190,195, 200,
                                                  205,210, 215,220, 225,230, 235,240, 245,250, 255,260]'''
                        #plt.legend(loc="lower left")
                    elif b == 'Flux Mip':  # contour_levels
                        p1 = plt.contour(formatted_dates, yid, zid, levels=[0, 1025.7], colors='black',linewidths=2.5)
                        p1 = plt.plot(MLD, color='hotpink', linewidth=2.5)#, linestyle='dashed'
                        #cmap = 'RdBu_r'
                        p1 = plt.contourf(formatted_dates, yi, np.log1p(zi),
                                          levels=[0, 0.40, 0.80, 1.20, 1.60, 2.0, 2.40, 2.80, 3.20, 3.60, 4, 4.40, 4.80, 5.20],
                                          cmap='viridis', alpha=1, extend="both", vmin=0, vmax=4.5)

                        #plt.legend(loc="lower left")
                        #0/0
                    elif b == 'Flux Map':  # contour_levels
                        p1 = plt.contour(formatted_dates, yid, zid, levels=[0, 1025.7], colors='black',linewidths=2.5)
                        p1 = plt.plot(MLD, color='hotpink', linewidth=2.5)#, linestyle='dashed'
                        #cmap = 'RdBu_r'
                        p1 = plt.contourf(formatted_dates, yi, np.log1p(zi),
                                          levels=[0, 0.40, 0.80, 1.20, 1.60, 2.0, 2.40, 2.80, 3.20, 3.60, 4, 4.40, 4.80, 5.20, 5.60, 6.0],
                                          cmap='viridis', alpha=1, extend="both", vmin=0, vmax=5)

                        #plt.legend(loc="lower left")
                    else:
                        p1 = plt.contourf(formatted_dates, yi, zi, contour_levels, cmap='viridis', alpha=1,
                                          extend="both")

                    # pgon=plt.Polygon(pd.DataFrame(xi)[~np.nan(xiZ[0])], -xiZ[~np.nan(xiZ)])
                    # ax1.add_patch(pgon)

                    #
                    ax1.fill_between(xii[np.isfinite(xiZ)['0']], np.array(np.transpose(-xiZ[np.isfinite(xiZ)['0']]))[0], -1000, alpha=1, color='white')

                    # plt.plot(pd.DataFrame(xi)[np.isfinite(xiZ)[0]],np.array(np.transpose(-xiZ[np.isfinite(xiZ)[0]]))[0], marker='o', linestyle='-',markersize=3.5,markeredgecolor='white', markeredgewidth=0.5)

                    # r = patches.Polygon(xi,np.array(np.transpose(-xiZ))[0], lw=2, ls='--',edgecolor='#ED7D31', facecolor='none', clip_on=False)
                    #
                    HIG = 500
                    if b == 'Particulate Organic Carbon': #'bbp POC':
                        '''high = -100
                        HIG = -0.65#5# 70
                        ha = 1.25
                        ha1 = 0.5'''

                        
                        high = -1000
                        HIG = -2.65#-1
                        ha = 4.25
                        ha1 = 2.5

                    elif b == 'Flux Mip':
                        high = -1000
                        ha = 4.25
                        ha1 = 2.5
                        HIG = -2.65#-1
                    elif b == 'Flux Map':
                        high = -1000
                        ha = 4.25
                        ha1 = 2.5
                        HIG = -2.65#-1
                    else:
                        # ax1.add_patch(r)#
                        high = -100
                        ha = 1.25  # 4.25
                        ha1 = 0.5  # 2.5
                        HIG = -0.65#5# 70



                    # plt.plot([formatted_dates[214],formatted_dates[214], formatted_dates[289],formatted_dates[289]], [0+3,high-10,high-10,0+3], linestyle='dashed', alpha=1,color='red')

                    '''ax1.fill_between([formatted_dates[214], formatted_dates[289]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    ax1.fill_between([formatted_dates[299], formatted_dates[389]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    ax1.fill_between([formatted_dates[459], formatted_dates[534]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    ax1.fill_between([formatted_dates[561], formatted_dates[633]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    ax1.fill_between([formatted_dates[760], formatted_dates[835]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    ax1.fill_between([formatted_dates[884], formatted_dates[959]], [0, 0], [high, high], alpha=0.2,color='cyan')'''

                    '''ax1.fill_between([formatted_dates[206], formatted_dates[301]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    ax1.fill_between([formatted_dates[307], formatted_dates[402]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    ax1.fill_between([formatted_dates[439], formatted_dates[543]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    ax1.fill_between([formatted_dates[569], formatted_dates[678]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    ax1.fill_between([formatted_dates[755], formatted_dates[869]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    ax1.fill_between([formatted_dates[884], formatted_dates[970]], [0, 0], [high, high], alpha=0.2,color='cyan')
                    '''
                    ###############
                    r = patches.Rectangle((xii[206], 0 + ha1), xii[301] - xii[206], high - ha, lw=2, ls='--', edgecolor='#ED7D31', facecolor='none', clip_on=False)
                    ax1.add_patch(r)
                    r = patches.Rectangle((xii[307], 0 + ha1), xii[402] - xii[307], high - ha, lw=2, ls='--', edgecolor='#ED7D31', facecolor='none', clip_on=False)
                    ax1.add_patch(r)
                    r = patches.Rectangle((xii[439], 0 + ha1), xii[543] - xii[439], high - ha, lw=2, ls='--', edgecolor='#ED7D31', facecolor='none', clip_on=False)
                    ax1.add_patch(r)
                    r = patches.Rectangle((xii[569], 0 + ha1), xii[678] - xii[569], high - ha, lw=2, ls='--', edgecolor='#ED7D31', facecolor='none', clip_on=False)
                    ax1.add_patch(r)
                    r = patches.Rectangle((xii[755], 0 + ha1), xii[869] - xii[755], high - ha, lw=2, ls='--', edgecolor='#ED7D31', facecolor='none', clip_on=False)
                    ax1.add_patch(r)
                    r = patches.Rectangle((xii[884], 0 + ha1), xii[970] - xii[884], high - ha, lw=2, ls='--', edgecolor='#ED7D31', facecolor='none', clip_on=False)
                    ax1.add_patch(r)

                    b1, = ax1.plot([], marker="", markersize=20, linestyle="-", color='hotpink',label="Mixed layer Depth")
                    b2, = ax1.plot([], marker="", markersize=20, linestyle="-", color="black",label="Isopycnal 1025.7 Kg/$m^3$")
                    b3, = ax1.plot([], marker="", markersize=20, linestyle="--", color='#ED7D31', label="Export Event")
                    ax1.legend(handles=[b1, b2, b3])#,loc="lower left"
                    #ax1.legend(loc="lower left")
                    handles, labels = ax1.get_legend_handles_labels()
                    ax1.legend(handles, labels, loc='lower left',fontsize=25)
                    ssi=25
                    plt.text(xi[250], -HIG, '1', fontsize=ssi, color='#ED7D31')
                    ax1.text(xi[350], -HIG, '2', fontsize=ssi, color='#ED7D31')
                    ax1.text(xi[490], -HIG, '3', fontsize=ssi, color='#ED7D31')
                    ax1.text(xi[595], -HIG, '4', fontsize=ssi, color='#ED7D31')
                    ax1.text(xi[800], -HIG, '5', fontsize=ssi, color='#ED7D31')
                    ax1.text(xi[922], -HIG, '6', fontsize=ssi, color='#ED7D31')
                    alphabet = ['C', 'A', 'B', 'D', 'F', 'E', 'B', 'C', '']
                    #alphabet = ['C', 'A', 'B', 'D', 'A', 'E', 'B', 'C', '']
                    ax1.text(-55, 4, alphabet[count_b], fontsize=25, color='black', weight='bold')


                    # ax1.hlines(y=300, xmin=min(formatted_dates), xmax=max(formatted_dates), colors='r')
                    cb = plt.colorbar(p1, orientation='vertical')
                    tick_locator = ticker.MaxNLocator(nbins=5)
                    cb.locator = tick_locator
                    #cb.set_label (fontsize=20)



                    if cl == 'AB':
                        cb.ax.set_ylabel(" []", fontsize=21)
                    elif cl == 'BV':
                        cb.ax.set_ylabel(" []", fontsize=21)
                    else:
                        cb.ax.set_ylabel(clist[count_b], fontsize=20)
                        if b == 'bbpPOC/Chla':
                            cb.ax.set_ylabel('log10 ' + clist[count_b], fontsize=21)
                        if b == 'Particulate Organic Carbon': # 'bbp POC':
                            # cb.ax.set_ylabel('log ' + clist[count_b], fontsize=15)
                            cb.ax.set_ylabel(clist[count_b], fontsize=21)
                            ticks = [float(t.get_text().replace('−', '-')) for t in cb.ax.get_yticklabels()]
                            expticks = np.exp(ticks) - 1
                            #cb.ax.set_yticklabels(np.round(expticks, 1))  # vertically oriented colorbar
                            cb.ax.set_yticklabels(np.round(expticks, 0))

                        if b == 'Flux Mip':
                            #cb.ax.set_ylabel('log ' + clist[count_b], fontsize=15)
                            cb.ax.set_ylabel(clist[count_b], fontsize=21)
                            ticks = [float(t.get_text().replace('−', '-')) for t in cb.ax.get_yticklabels()]
                            expticks = np.exp(ticks) - 1
                            #cb.ax.set_yticklabels(np.round(expticks, 1))  # vertically oriented colorbar
                            cb.ax.set_yticklabels(np.round(expticks, 0))
                        if b == 'Flux Map':
                            cb.ax.set_ylabel(clist[count_b], fontsize=21)
                            ticks = [float(t.get_text().replace('−', '-')) for t in cb.ax.get_yticklabels()]
                            expticks = np.exp(ticks) - 1
                            #cb.ax.set_yticklabels(np.round(expticks, 1))  # vertically oriented colorbar
                            cb.ax.set_yticklabels(np.round(expticks, 0))
                        if b == 'Chlorophyll a': #
                            cb.ax.set_ylabel(clist[count_b], fontsize=21)
                            ticks = [float(t.get_text().replace('−', '-')) for t in cb.ax.get_yticklabels()]
                            expticks = np.exp(ticks) - 1
                            cb.ax.set_yticklabels(np.round(expticks, 1))  # vertically oriented colorbar
                            #cb.ax.set_yticklabels(np.round(expticks, 0))

                        cb.ax.tick_params(labelsize=18)


                    plt.ylabel("Depth [m]", fontsize=22)
                    # plt.xlim(min(df_temp.date), max(df_temp.date))
                    plt.xlim(min(formatted_dates), max(formatted_dates))
                    #ax1.set_xticks(formatted_dates[::60])

                    formatted_datesq=[]
                    if b == 'Chlorophyll aA': #
                        for ii in [0,88 ,259, 431, 603 ,777, 943]:
                            formatted_datesq.append(formatted_dates[ii])
                    elif b == 'Particulate Organic CarbonA': # 'bbp POC':
                        for ii in [0,88 ,259, 431, 603 ,777, 943]:
                            formatted_datesq.append(formatted_dates[ii])
                    else:
                        for ii in [0,88 ,172 ,259 ,347, 431, 518, 603 ,690 ,777, 856 ,943]:
                            formatted_datesq.append(formatted_dates[ii])

                    ax1.set_xticks(formatted_datesq)

                    if high==100:
                        ax1.set_yticks(np.arange(high, 10, 10))
                    elif high==1000:
                        ax1.set_yticks(np.arange(high, 10, 100))
                        #ax1.set_yticks(np.arange(high, 0, 10))


                    #ax1.set_xticklabels(formatted_datesym[::60])
                    ax1.set_xticks(formatted_datesq)
                    plt.ylim(high, 0)
                    ax1.tick_params(axis='x', labelsize=22) #18
                    ax1.tick_params(axis='y', labelsize=22)



                    #plt.xlabel("Time [days]")
                    plt.xlabel("Date", fontsize=22) #17
                    plt.xticks(rotation=25)
                    if cl == 'f':
                        plt.title(str(t) + str(b), fontsize=25, x=0.5,y=1.06)
                    else:
                        plt.title(str(t) + ":" + str(b), fontsize=25, x=0.5,y=1.06)
                    bb=1





                    if bb == 1:
                        plt.savefig(str(path_to_figures) + '/' + 'Aout2025lARS' + str(
                            -high) + 'mResoluEEE' + str(b) + str(tp_list_char.date[count_tp]) +
                                    ' to ' + str(tp_list_char.date[count_tp + 1]) + 'Year2024.png', dpi=300)

                    pp.savefig(fig, dpi=1200)  # Save each figure in the pdf #plt.gcf()


                    # Add one to the bmax counter
                    count_b += 1

                # close the pdf document
                pp.close()
                plt.close('all')
# plt.plot(pd.DataFrame(xi),-xiZ,marker='o', linestyle='--',markersize=3.5,markeredgecolor='white', markeredgewidth=0.5)

1 + 1



''' for j, lab in enumerate(['$0$','$1$','$2$','$>3$']):
    cbar.ax.text(.5, (2 * j + 1) / 8.0, lab, ha='center', va='center')
    tick_locator '''