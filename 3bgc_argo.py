
#Script for treat bgc argo data

        # Setup

import netCDF4 as nc
import numpy as np
import pandas as pd
#import seawater as sw
#import gsw
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import xarray as xr
#import mplcursors
import seaborn as sns
import math
import matplotlib.cm as cm
from matplotlib.patches import Patch
from matplotlib import patches
import matplotlib.dates as mdates
import warnings # these two lines to remove the annoying warnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)
########################################################################################################################
                                        # Import the netcdf file #
########################################################################################################################

# Open the netCDF file as an xarray dataset
data = xr.open_dataset('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/6903096_Sprof.nc')
dac = xr.open_dataset('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/6903096_Sprof.nc')
data = data.rename({'CYCLE_NUMBER':'PROF_NUM'}).swap_dims({'N_PROF':'PROF_NUM'})

data2025_check= xr.open_dataset('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/nc files dowload '
                                'from the fleet/GL_PR_PF_6903096.nc')
#data2025_check = data2025_check.rename({'POSITION':'PROF_NUM'}).swap_dims({'POSITION':'PROF_NUM'})
# Create empty DataFrame
dfcheck = pd.DataFrame()
var_listcheck = ['PRES', 'TEMP', 'PSAL', 'BBP700']  # Add any other variables you want to extract

# Loop over each profile and extract variables
for i in range(len(data2025_check['POSITION'])):
    row_dict = {'POSITION': data2025_check['POSITION'][i].values}
    for var_name in var_listcheck:
        row_dict[var_name] = data2025_check[var_name].isel(TIME=i).values
    row = pd.DataFrame(row_dict)
    dfcheck = pd.concat([dfcheck, row], ignore_index=True)
# Get list of variable names to extract
var_list = ['PRES', 'JULD', 'TEMP', 'PSAL', 'CHLA_ADJUSTED', 'DOXY_ADJUSTED', 'BBP700', 'LONGITUDE', 'LATITUDE']  # Add any other variables you want to extract

# Create empty DataFrame
df = pd.DataFrame()

# Loop over each profile and extract variables
for i in range(len(data['PROF_NUM'])):
    row_dict = {'PROF_NUM': data['PROF_NUM'][i].values}
    for var_name in var_list:
        row_dict[var_name] = data[var_name].isel(PROF_NUM=i).values
    row = pd.DataFrame(row_dict)
    df = pd.concat([df, row], ignore_index=True)

bgc_angola = df

# rename PROF_NUMBER column into CYCLE_NUMBER
bgc_angola = bgc_angola.rename(columns={'PROF_NUM': 'CYCLE_NUMBER'})

# rename JULD column into Date_Time
bgc_angola = bgc_angola.rename(columns={'JULD': 'Date_Time'})

# Delete the station number 157 because the CTD was not working during this profile (profile name : 0156a_WMO6903095)
bgc_angola = bgc_angola.drop(bgc_angola[bgc_angola['CYCLE_NUMBER'] == 157].index)

# Delete the station number 185 because the UVP was not working during this profile (profile name : 0184a_WMO6903095)
bgc_angola = bgc_angola.drop(bgc_angola[bgc_angola['CYCLE_NUMBER'] == 185].index)

# drop rows where 'PRES' column has NaN values
bgc_angola.dropna(subset=['PRES'], inplace=True)

# delete rows with negative values in the 'PRES' column
bgc_angola = bgc_angola[bgc_angola['PRES'] >= 0]

# delete negative values in chlorophylle a
bgc_angola['CHLA_ADJUSTED'] = bgc_angola['CHLA_ADJUSTED'].apply(lambda x: 0 if x < 0 else x)

# convert pressure into depth
depth = sw.eos80.dpth(bgc_angola['PRES'], bgc_angola['LATITUDE'])
bgc_angola['depth'] = depth

# I compute the potential density: for that, I need absolute salinity and conservative temperature, so I transform
# salinity and temperature first
abs_sal = gsw.SA_from_SP(bgc_angola['PSAL'], bgc_angola['PRES'], bgc_angola['LONGITUDE'], bgc_angola['LATITUDE'])
bgc_angola['ABS_SAL'] = abs_sal

cons_temp = gsw.CT_from_t(bgc_angola['PSAL'], bgc_angola['TEMP'], bgc_angola['PRES'])
bgc_angola['CONS_TEMP'] = cons_temp

dens = gsw.density.sigma0(abs_sal, cons_temp)
bgc_angola['dens'] = dens + 1000

# Compute potential temperature based on absolute salinity, in situ temperature and pressure
POT_TEMP = gsw.conversions.pt0_from_t(bgc_angola['PSAL'], bgc_angola['TEMP'], bgc_angola['PRES'])
bgc_angola['POT_TEMP'] = POT_TEMP
#save it in order to clean ctd cast and compute clines with the R script : /Users/Alexandre/GIT/ecotaxa/Scripts/clines.R
bgc_angola.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_angola.csv', index = False)

# Create a df with station and date
Profiles_dates = bgc_angola.loc[:, ['CYCLE_NUMBER', 'Date_Time']]
Profiles_dates = Profiles_dates.drop_duplicates()

# Clean ctd casts and compute clines and MLD with JO package castr

station_cline_depths = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/station_clines.csv') # data set compute
# with r script /Users/Alexandre/GIT/ecotaxa/Scripts/clines.R
station_cline_depths = pd.merge(station_cline_depths, Profiles_dates, on='CYCLE_NUMBER')

########################################################################################################################
                                                  # Plot #
########################################################################################################################


# Interactive plot of all temperature profiles
# Create a figure and axes
fig, ax = plt.subplots()

# Set up the color gradient
cmap = cm.get_cmap('rainbow')
norm = plt.Normalize(vmin=bgc_angola['CYCLE_NUMBER'].min(), vmax=bgc_angola['CYCLE_NUMBER'].max())

# Loop over each profile
for station in bgc_angola['CYCLE_NUMBER'].unique():
    # Select data for this profile
    profile_data = bgc_angola[bgc_angola['CYCLE_NUMBER'] == station]

    # Plot a line for this profile with a color gradient
    ax.plot(profile_data['TEMP'], profile_data['depth'],
            color=cmap(norm(station)), label=station)

# Add axis labels and a legend
ax.set_xlabel('Temperature [degrees Celsius]')
ax.xaxis.set_ticks_position('top')  # put x-axis on top
ax.xaxis.set_label_position('top')
ax.set_ylabel('Depth [m]')
#ax.legend()

# Reverse the y-axis
ax.invert_yaxis()

# Add annotations
#annotations = ax.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
 #                         bbox=dict(boxstyle="round", fc="w"),
   #                       arrowprops=dict(arrowstyle="->"))
#annotations.set_visible(False)

#def update_annot(ind):
    #pos = ax.get_yticklabels()[ind['ind'][0]].get_position()
    #annotations.xy = pos
    #text = bgc_angola['Profile'].unique()[ind['ind'][0]]
    #annotations.set_text(text)

#def hover(event):
    #vis = annotations.get_visible()
    #if event.inaxes == ax:
     #   cont, ind = ax.contains(event)
      #  if cont:
       #     update_annot(ind)
        #    annotations.set_visible(True)
         #   fig.canvas.draw_idle()
        #else:
         #   if vis:
          #      annotations.set_visible(False)
           #     fig.canvas.draw_idle()

#mplcursors.cursor(ax, hover=True)

# Show the plot
plt.show()


# Interactive plot of all density profiles
# Create a figure and axes
fig, ax = plt.subplots()

# Set up the color gradient
cmap = cm.get_cmap('rainbow')
norm = plt.Normalize(vmin=bgc_angola['CYCLE_NUMBER'].min(), vmax=bgc_angola['CYCLE_NUMBER'].max())

# Loop over each profile
for station in bgc_angola['CYCLE_NUMBER'].unique():
    # Select data for this profile
    profile_data = bgc_angola[bgc_angola['CYCLE_NUMBER'] == station]

    # Plot a line for this profile with a color gradient
    ax.plot(profile_data['dens'], profile_data['depth'],
            color=cmap(norm(station)), label=station)

# Add axis labels and a legend
ax.set_xlabel('Potential density [kg/m3]')
ax.xaxis.set_ticks_position('top')  # put x-axis on top
ax.xaxis.set_label_position('top')
ax.set_ylabel('Depth [m]')

# Reverse the y-axis
ax.invert_yaxis()

# Show the plot
plt.show()

# Interactive plot of all salinity profiles
# Create a figure and axes
fig, ax = plt.subplots()
# Set up the color gradient
cmap = cm.get_cmap('rainbow')
norm = plt.Normalize(vmin=bgc_angola['CYCLE_NUMBER'].min(), vmax=bgc_angola['CYCLE_NUMBER'].max())
# Loop over each profile
for station in bgc_angola['CYCLE_NUMBER'].unique():
    # Select data for this profile
    profile_data = bgc_angola[bgc_angola['CYCLE_NUMBER'] == station]

    # Plot a line for this profile
    ax.plot(profile_data['PSAL'], profile_data['depth'],
            color=cmap(norm(station)), label=station)

# Add axis labels and a legend
ax.set_xlabel('Practical salinity [psu]')
ax.xaxis.set_ticks_position('top')  # put x-axis on top
ax.xaxis.set_label_position('top')
ax.set_ylabel('Depth [m]')

# Reverse the y-axis
ax.invert_yaxis()

# Show the plot
plt.show()

# Interactive plot of all chla profiles
# Create a figure and axes
fig, ax = plt.subplots()
# Set up the color gradient
cmap = cm.get_cmap('rainbow')
norm = plt.Normalize(vmin=bgc_angola['CYCLE_NUMBER'].min(), vmax=bgc_angola['CYCLE_NUMBER'].max())
# Loop over each profile
for station in bgc_angola['CYCLE_NUMBER'].unique():
    # Select data for this profile and drop NaN values in CHLA_ADJUSTED column
    profile_data = bgc_angola[bgc_angola['CYCLE_NUMBER'] == station].dropna(subset=['CHLA_ADJUSTED'])

    # Plot a line for this profile
    ax.plot(profile_data['CHLA_ADJUSTED'], profile_data['depth'],
            color=cmap(norm(station)), label=station, linestyle='-', marker=None)

# Add axis labels and a legend
ax.set_xlabel('Chlorophyll-a [mg/m3]')
ax.xaxis.set_ticks_position('top')  # put x-axis on top
ax.xaxis.set_label_position('top')
ax.set_ylabel('Depth [m]')

# Reverse the y-axis
ax.invert_yaxis()

# Show the plot
plt.show()

# Interactive plot of all bbp700 profiles
# Create a figure and axes
fig, ax = plt.subplots()
# Set up the color gradient
cmap = cm.get_cmap('rainbow')
norm = plt.Normalize(vmin=bgc_angola['CYCLE_NUMBER'].min(), vmax=bgc_angola['CYCLE_NUMBER'].max())
# Loop over each profile
for station in bgc_angola['CYCLE_NUMBER'].unique():
    # Select data for this profile
    profile_data = bgc_angola[bgc_angola['CYCLE_NUMBER'] == station].dropna(subset=['BBP700'])

    # Plot a line for this profile
    ax.plot(profile_data['BBP700'], profile_data['depth'],
            color=cmap(norm(station)), label=station, linestyle='-', marker=None)

# Add axis labels and a legend
ax.set_xlabel('BBP700') # add units of bbp
ax.xaxis.set_ticks_position('top')  # put x-axis on top
ax.xaxis.set_label_position('top')
ax.set_ylabel('Depth [m]')

# Reverse the y-axis
ax.invert_yaxis()

# Show the plot
plt.show()

# Interactive plot of all oxygen profiles
# Create a figure and axes
fig, ax = plt.subplots()
# Set up the color gradient
cmap = cm.get_cmap('rainbow')
norm = plt.Normalize(vmin=bgc_angola['CYCLE_NUMBER'].min(), vmax=bgc_angola['CYCLE_NUMBER'].max())
# Loop over each profile
for station in bgc_angola['CYCLE_NUMBER'].unique():
    # Select data for this profile
    profile_data = bgc_angola[bgc_angola['CYCLE_NUMBER'] == station].dropna(subset=['DOXY_ADJUSTED'])

    # Plot a line for this profile
    ax.plot(profile_data['DOXY_ADJUSTED'], profile_data['depth'],
            color=cmap(norm(station)), label=station, linestyle='-', marker=None)

# Add axis labels and a legend
ax.set_xlabel('Doxy [micromol/kg]') # add units of bbp
ax.xaxis.set_ticks_position('top')  # put x-axis on top
ax.xaxis.set_label_position('top')
ax.set_ylabel('Depth [m]')

# Reverse the y-axis
ax.invert_yaxis()

# Show the plot
plt.show()

########################################################################################################################
                                        # Time series #
########################################################################################################################


# Define depth intervals
depth_intervals = [(0, 100), (100, 300), (300, 600), (600, 1000), (1000, 2000)]
depth_interval_names = ['0-100m', '100-300m', '300-600m', '600-2000m', '1000-2000m']

# Group data by depth interval and time, and calculate mean temperature
mean_temp = []
std_temp = []
for interval, name in zip(depth_intervals, depth_interval_names):
    depth_mask = (bgc_angola['depth'] >= interval[0]) & (bgc_angola['depth'] < interval[1])
    depth_data = bgc_angola[depth_mask].groupby('Date_Time')['TEMP'].agg(['mean', 'std'])
    depth_data.columns = [name + '_mean', name + '_std']
    mean_temp.append(depth_data[name + '_mean'])
    std_temp.append(depth_data[name + '_std'])
mean_temp = pd.concat(mean_temp, axis=1)
std_temp = pd.concat(std_temp, axis=1)

# Create a figure and axes
fig, ax = plt.subplots()

# Plot the data
sns.lineplot(data=mean_temp, ax=ax, marker = 'o', markersize=3.5,
                 markeredgecolor='white', markeredgewidth = 0.5)
for line in ax.lines:
    line.set_linestyle('--')  # add dashed lines

# Add error bars representing the standard deviation
for i, col in enumerate(mean_temp.columns):
    ax.fill_between(mean_temp.index, mean_temp[col] - std_temp.iloc[:, i], mean_temp[col] + std_temp.iloc[:, i],
                    alpha=.1)

# Add axis labels and a legend
ax.set_xlabel('Date')
ax.set_ylabel('Temperature [degrees Celsius]')
ax.legend()

# Show the plot
plt.show()


# Group data by depth interval and time, and calculate mean potential density
mean_dens = []
std_dens = []
for interval, name in zip(depth_intervals, depth_interval_names):
    depth_mask = (bgc_angola['depth'] >= interval[0]) & (bgc_angola['depth'] < interval[1])
    depth_data = bgc_angola[depth_mask].groupby('Date_Time')['dens'].agg(['mean', 'std'])
    depth_data.columns = [name + '_mean', name + '_std']
    mean_dens.append(depth_data[name + '_mean'])
    std_dens.append(depth_data[name + '_std'])
mean_dens = pd.concat(mean_dens, axis=1)
std_dens = pd.concat(std_dens, axis=1)

# Create a figure and axes
fig, ax = plt.subplots()

# Plot the data
sns.lineplot(data=mean_dens, ax=ax, marker = 'o', markersize=3.5,
                 markeredgecolor='white', markeredgewidth = 0.5)
for line in ax.lines:
    line.set_linestyle('--')  # add dashed lines

# Add error bars representing the standard deviation
for i, col in enumerate(mean_dens.columns):
    ax.fill_between(mean_dens.index, mean_dens[col] - std_dens.iloc[:, i], mean_dens[col] + std_dens.iloc[:, i],
                    alpha=.1)

# Add axis labels and a legend
ax.set_xlabel('Date')
ax.set_ylabel('Potential density [kg/m3]')
ax.legend()

# Show the plot
plt.show()

# Group data by depth interval and time, and calculate mean practical salinity
mean_sal = []
std_sal = []
for interval, name in zip(depth_intervals, depth_interval_names):
    depth_mask = (bgc_angola['depth'] >= interval[0]) & (bgc_angola['depth'] < interval[1])
    depth_data = bgc_angola[depth_mask].groupby('Date_Time')['PSAL'].agg(['mean', 'std'])
    depth_data.columns = [name + '_mean', name + '_std']
    mean_sal.append(depth_data[name + '_mean'])
    std_sal.append(depth_data[name + '_std'])
mean_sal = pd.concat(mean_sal, axis=1)
std_sal = pd.concat(std_sal, axis=1)

# Create a figure and axes
fig, ax = plt.subplots()

# Plot the data
sns.lineplot(data=mean_sal, ax=ax, marker = 'o', markersize=3.5,
                 markeredgecolor='white', markeredgewidth = 0.5)
for line in ax.lines:
    line.set_linestyle('--')  # add dashed lines

# Add error bars representing the standard deviation
for i, col in enumerate(mean_sal.columns):
    ax.fill_between(mean_sal.index, mean_sal[col] - std_sal.iloc[:, i], mean_sal[col] + std_sal.iloc[:, i],
                    alpha=.1)

# Add axis labels and a legend
ax.set_xlabel('Date')
ax.set_ylabel('Practical salinity [psu]')
ax.legend()

# Show the plot
plt.show()


# Group data by depth interval and time, and calculate mean oxygen
mean_oxy = []
std_oxy = []
for interval, name in zip(depth_intervals, depth_interval_names):
    depth_mask = (bgc_angola['depth'] >= interval[0]) & (bgc_angola['depth'] < interval[1])
    depth_data = bgc_angola[depth_mask].groupby('Date_Time')['DOXY_ADJUSTED'].agg(['mean', 'std'])
    depth_data.columns = [name + '_mean', name + '_std']
    mean_oxy.append(depth_data[name + '_mean'])
    std_oxy.append(depth_data[name + '_std'])
mean_oxy = pd.concat(mean_oxy, axis=1)
std_oxy = pd.concat(std_oxy, axis=1)

# Create a figure and axes
fig, ax = plt.subplots()

# Plot the data
sns.lineplot(data=mean_oxy, ax=ax, marker = 'o', markersize=3.5,
                 markeredgecolor='white', markeredgewidth = 0.5)
for line in ax.lines:
    line.set_linestyle('--')  # add dashed lines

# Add error bars representing the standard deviation
for i, col in enumerate(mean_oxy.columns):
    ax.fill_between(mean_oxy.index, mean_oxy[col] - std_oxy.iloc[:, i], mean_oxy[col] + std_oxy.iloc[:, i],
                    alpha=.1)

# Add axis labels and a legend
ax.set_xlabel('Date')
ax.set_ylabel('Doxy [micromol/kg]')
ax.legend()

# Show the plot
plt.show()


# Group data by date and calculate mean chla

# select only 300 first meters (because no more chla deeper)
bgc_angola_chla = bgc_angola[bgc_angola['depth'] <= 300]

# Group data by date and calculate mean and std
mean_chla = bgc_angola_chla.groupby('Date_Time')['CHLA_ADJUSTED'].agg(['mean', 'std'])

# Create a figure and axes
fig, ax = plt.subplots()

# Plot the data
sns.lineplot(data=mean_chla['mean'], ax=ax, marker = 'o', markersize=3.5,
                 markeredgecolor='white', markeredgewidth = 0.5, linestyle = '--')

lower_err = np.clip(mean_chla['mean'] - mean_chla['std'], 0, np.inf)
ax.fill_between(mean_chla.index, lower_err, mean_chla['mean'] + mean_chla['std'], alpha=.1)

# Add axis labels and a legend
ax.set_xlabel('Date')
ax.set_ylabel('Chlorophylle-a [mg/m3]')
plt.title('Mean Chlorophyll-a (0-300m)')

# Show the plot
plt.show()

# Group data by depth interval and time, and calculate mean BBP700
mean_bbp = []
std_bbp = []
for interval, name in zip(depth_intervals, depth_interval_names):
    depth_mask = (bgc_angola['depth'] >= interval[0]) & (bgc_angola['depth'] < interval[1])
    depth_data = bgc_angola[depth_mask].groupby('Date_Time')['BBP700'].agg(['mean', 'std'])
    depth_data.columns = [name + '_mean', name + '_std']
    mean_bbp.append(depth_data[name + '_mean'])
    std_bbp.append(depth_data[name + '_std'])
mean_bbp = pd.concat(mean_bbp, axis=1)
std_bbp = pd.concat(std_bbp, axis=1)

# Create a figure and axes
fig, ax = plt.subplots()

# Plot the data
sns.lineplot(data=mean_bbp, ax=ax, marker = 'o', markersize=3.5,
                 markeredgecolor='white', markeredgewidth = 0.5)
for line in ax.lines:
    line.set_linestyle('--')  # add dashed lines

# Add error bars representing the standard deviation
for i, col in enumerate(mean_bbp.columns):
    lower_error = np.maximum(mean_bbp[col] - std_bbp.iloc[:, i], 0)  # clip negative values to zero
    ax.fill_between(mean_bbp.index, lower_error, mean_bbp[col] + std_bbp.iloc[:, i], alpha=.1)

# Add axis labels and a legend
ax.set_xlabel('Date')
ax.set_ylabel('BBP700') #add units of BBP700
ax.legend()

# Show the plot
plt.show()


########################################################################################################################
                                        # Compute 5m bins #
########################################################################################################################

# create a function to round to the closest 5 value below the number (2 options)
def round_to_multiple(number, multiple):
    from math import floor
    return multiple * floor(number / multiple)
def rounditup(x, precision, method):
    if method == "floor":
        return math.floor(x / precision) * precision
    elif method == "ceiling":
        return math.ceil(x / precision) * precision
    else:
        return "give the parameter floor or ceiling"

# use list comprehension to create the column depth_bin using the round_to_multiple function
# m93["depth_bin"] = [round_to_multiple(m93["depth"][r], 5) + 2.5 for r in range(len(m93))]

bgc_binned = bgc_angola
bgc_binned.loc[:,'depth_bin'] = 99999
bgc_binned['depth_bin'] = bgc_binned["depth"].apply(lambda x: rounditup(x, 5, "floor")+ 2.5)
bgc_binned = bgc_binned.groupby(['CYCLE_NUMBER', 'Date_Time', 'depth_bin']).mean().reset_index()
bgc_binned = bgc_binned.drop('depth', axis = 1)
bgc_binned.columns = bgc_binned.columns.str.replace('_bin', '')

# add a column with depth levels
bgc_binned['layer'] = np.where(bgc_binned.depth < 100, "0-100m",
                      np.where(np.logical_and(bgc_binned.depth >= 100, bgc_binned.depth < 300), "100-300m",
                      np.where(np.logical_and(bgc_binned.depth >= 300, bgc_binned.depth < 600), "300-600m",
                      np.where(np.logical_and(bgc_binned.depth >= 600, bgc_binned.depth < 1000), "600-1000m",
                      np.where(bgc_binned.depth > 1000, "1000-2000m","else")))))

bgc_binned = bgc_binned.dropna()
bgc_binned_1000m = bgc_binned[bgc_binned['layer'] != '1000-2000m']
bgc_binned_600m = bgc_binned[bgc_binned['layer'] != '600-1000m']
bgc_binned_600m = bgc_binned_600m[bgc_binned_600m['layer'] != '1000-2000m']
bgc_binned_100m = bgc_binned[bgc_binned['layer'] == '0-100m']
bgc_binned.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_ben_binned.csv', index=False)
bgc_binned_1000m.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_ben_binned_1000m.csv', index=False)
bgc_binned_600m.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_ben_binned_600m.csv', index=False)
bgc_binned_100m.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_ben_binned_100m.csv', index=False)
# compute the mean of the 100 first meters to use it in the PCA script
bgc_binned_100 = bgc_binned.groupby(['CYCLE_NUMBER', 'layer']).mean().reset_index()
# select only 100 first meters
bgc_binned_100 = bgc_binned_100[bgc_binned_100['layer'] == '0-100m']
bgc_binned_100 = bgc_binned_100.drop(columns = ['depth', 'LATITUDE', 'LONGITUDE', 'CYCLE_NUMBER', 'layer', 'PRES'])

# save it
bgc_binned_100.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_ben_binned_100m.csv', index = False)



########################################################################################################################
                          # plot time series according to clustering results #
########################################################################################################################

# import bgc data clean, despike and binned (5m)

bgc_angola_clean = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_angola_clean.csv')

# select only 600 first meters
bgc_angola_clean = bgc_angola_clean[bgc_angola_clean['binned_depth'] <= 600]


# import station with cluster

station_cluster = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/station_cluster_600.csv')
station_cluster_kmeans = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/station_cluster_kmeans_600.csv')
bgc_binned_600m = pd.merge(bgc_angola_clean, station_cluster_kmeans, on='CYCLE_NUMBER')


# represent profile group by cluster

# Create a figure and axes
fig, ax = plt.subplots()
# Set up the color gradient
cmap = cm.get_cmap('rainbow', len(bgc_binned_600m['cluster'].unique()))
norm = plt.Normalize(vmin=bgc_binned_600m['cluster'].min(), vmax=bgc_binned_600m['cluster'].max())

# Create a dictionary to store the colors for each cluster
colors = {}
# Loop over each profile
for station in bgc_binned_600m['CYCLE_NUMBER'].unique():
    # Select data for this profile and drop NaN values in CHLA_ADJUSTED column
    profile_data = bgc_binned_600m[bgc_binned_600m['CYCLE_NUMBER'] == station].dropna(subset=['TEMP'])

    # Get the cluster for this profile
    cluster = profile_data['cluster'].iloc[0]

    # Plot a line for this profile with the color associated with its cluster
    ax.plot(profile_data['TEMP'], profile_data['binned_depth'],
            color=cmap(norm(cluster)), label=f'Cluster {cluster}', linestyle='-', marker=None)

    # Store the color for this cluster in the colors dictionary
    if cluster not in colors:
        colors[cluster] = cmap(norm(cluster))
# Add axis labels and a legend
#ax.set_xlabel('Potential density [kg/m3]')
ax.set_xlabel('Temperature [degrees Celsius]')
ax.xaxis.set_ticks_position('top')  # put x-axis on top
ax.xaxis.set_label_position('top')
ax.set_ylabel('Depth [m]')

# Reverse the y-axis
ax.invert_yaxis()
# Create the legend based on the colors dictionary
legend_elements = [Patch(facecolor=colors[cluster], label=f'Cluster {cluster}') for cluster in colors]
ax.legend(handles=legend_elements, loc='lower right')
# Show the plot
plt.show()


bgc_binned_600m = pd.merge(bgc_binned_600m, Profiles_dates, on='CYCLE_NUMBER')

# Define depth intervals
depth_intervals = [(0, 100), (100, 300), (300, 600)]
depth_interval_names = ['0-100m', '100-300m', '300-600m']

# Group data by depth interval and time, and calculate mean temperature
mean_temp = []
std_temp = []
for interval, name in zip(depth_intervals, depth_interval_names):
    depth_mask = (bgc_binned_600m['binned_depth'] >= interval[0]) & (bgc_binned_600m['binned_depth'] < interval[1])
    depth_data = bgc_binned_600m[depth_mask].groupby('Date_Time')['TEMP'].agg(['mean', 'std'])
    depth_data.columns = [name + '_mean', name + '_std']
    mean_temp.append(depth_data[name + '_mean'])
    std_temp.append(depth_data[name + '_std'])
mean_temp = pd.concat(mean_temp, axis=1)
std_temp = pd.concat(std_temp, axis=1)


# Create a figure and axes
fig, ax = plt.subplots()

# Plot the data
sns.lineplot(data=mean_temp, ax=ax, marker = 'o', markersize=3.5,
                 markeredgecolor='white', markeredgewidth = 0.5)
for line in ax.lines:
    line.set_linestyle('--')  # add dashed lines

# Add error bars representing the standard deviation
for i, col in enumerate(mean_temp.columns):
    ax.fill_between(mean_temp.index, mean_temp[col] - std_temp.iloc[:, i], mean_temp[col] + std_temp.iloc[:, i],
                    alpha=.1)

# Add background colors based on the 'cluster' column
cmap = cm.get_cmap('rainbow', len(bgc_binned_600m['cluster'].unique()))
norm = plt.Normalize(vmin=bgc_binned_600m['cluster'].min(), vmax=bgc_binned_600m['cluster'].max())
for i, (start, end, cluster) in enumerate(zip(bgc_binned_600m['Date_Time'].iloc[:-1],
                                              bgc_binned_600m['Date_Time'].iloc[1:],
                                              bgc_binned_600m['cluster'].iloc[:-1])):
    rect = patches.Rectangle((start, -100), end-start, 2000, linewidth=0, alpha=0.2, facecolor=cmap(norm(cluster)))
    ax.add_patch(rect)

# Convert dates to numerical values
date1 = mdates.datestr2num('2021-08-01')
date2 = mdates.datestr2num('2021-08-11')

# Add vertical dashed lines
''' ax.axvline(date1, linestyle='--', color='white')
ax.axvline(date2, linestyle='--', color='white')
'''

# Add axis labels and a legend
ax.set_xlabel('Date')
ax.set_ylabel('Temperature [°C]')
ax.set_ylim(0, 25)  # Set y-axis limits to 0-30 degrees Celsius
ax.legend(loc='upper right')

# Show the plot
plt.show()


# Plot the MLD time series

station_cline_depths = pd.merge(station_cline_depths, station_cluster, on='CYCLE_NUMBER')


# Create a figure and axes
fig, ax = plt.subplots()
# Plot the data
sns.lineplot(data=station_cline_depths, x = station_cline_depths['Date_Time'], y = station_cline_depths['MLD'], ax=ax, marker = 'o', markersize=3.5,
                 markeredgecolor='white', markeredgewidth = 0.5, linestyle = '--')
# Add background colors based on the 'cluster' column
cmap = cm.get_cmap('rainbow', len(station_cline_depths['cluster'].unique()))
norm = plt.Normalize(vmin=station_cline_depths['cluster'].min(), vmax=station_cline_depths['cluster'].max())
for i, (start, end, cluster) in enumerate(zip(station_cline_depths['Date_Time'].iloc[:-1],
                                              station_cline_depths['Date_Time'].iloc[1:],
                                              station_cline_depths['cluster'].iloc[:-1])):
    rect = patches.Rectangle((start, -100), end-start, 2000, linewidth=0, alpha=0.2, facecolor=cmap(norm(cluster)))
    ax.add_patch(rect)
# Convert dates to numerical values
date1 = mdates.datestr2num('2021-08-01')
date2 = mdates.datestr2num('2021-08-11')
# Add vertical dashed lines
''' ax.axvline(date1, linestyle='--', color='white')
ax.axvline(date2, linestyle='--', color='white')
'''
# Add annotations
annotations = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                          bbox=dict(boxstyle="round", fc="w"),
                          arrowprops=dict(arrowstyle="->"))
annotations.set_visible(False)
def update_annot(ind):
    pos = ax.get_yticklabels()[ind['ind'][0]].get_position()
    text = '\n'.join(station_cline_depths.iloc[ind['ind'][0]][['CYCLE_NUMBER']])
    annotations.xy = pos
    annotations.set_text(text)
def hover(event):
    vis = annotations.get_visible()
    if event.inaxes == ax:
        cont, ind = ax.contains(event)
        if cont:
            update_annot(ind)
            annotations.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annotations.set_visible(False)
                fig.canvas.draw_idle()
mplcursors.cursor(ax, hover=True)
# Add axis labels and a legend
ax.set_xlabel('Date')
ax.set_ylabel('MLD [m]')
plt.title('Mixed Layer Depth [m]')
ax.set_ylim(-5, 360)
# Show the plot
plt.show()
# Reverse the y-axis
ax.invert_yaxis()


# plot clines time series

# Create a figure and axes
fig, ax = plt.subplots()

# Plot the data for pycno
#sns.lineplot(data=station_cline_depths, x='Date_Time', y='pycnocline', ax=ax, marker='o', markersize=3.5,
 #            markeredgecolor='white', markeredgewidth=0.5, linestyle='--', label='pycnocline', color='red')

# Plot the data for halo
#sns.lineplot(data=station_cline_depths, x='Date_Time', y='halocline', ax=ax, marker='o', markersize=3.5,
  #           markeredgecolor='white', markeredgewidth=0.5, linestyle='--', label='halo', color='green')

# Plot the data for thermo
sns.lineplot(data=station_cline_depths, x='Date_Time', y='thermocline', ax=ax, marker='o', markersize=3.5,
             markeredgecolor='white', markeredgewidth=0.5, linestyle='--', label='thermocline', color = 'orange')

# Add background colors based on the 'cluster' column
cmap = cm.get_cmap('rainbow', len(station_cline_depths['cluster'].unique()))
norm = plt.Normalize(vmin=station_cline_depths['cluster'].min(), vmax=station_cline_depths['cluster'].max())
for i, (start, end, cluster) in enumerate(zip(station_cline_depths['Date_Time'].iloc[:-1],
                                              station_cline_depths['Date_Time'].iloc[1:],
                                              station_cline_depths['cluster'].iloc[:-1])):
    rect = patches.Rectangle((start, -100), end-start, 2000, linewidth=0, alpha=0.2, facecolor=cmap(norm(cluster)))
    ax.add_patch(rect)

# Convert dates to numerical values
date1 = mdates.datestr2num('2021-08-01')
date2 = mdates.datestr2num('2021-08-11')
# Add vertical dashed lines
''' ax.axvline(date1, linestyle='--', color='white')
ax.axvline(date2, linestyle='--', color='white')
'''
# Add annotations
#annotations = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
 #                         bbox=dict(boxstyle="round", fc="w"),
  #                        arrowprops=dict(arrowstyle="->"))
#annotations.set_visible(False)
#def update_annot(ind):
 #   pos = ax.get_yticklabels()[ind['ind'][0]].get_position()
  #  text = '\n'.join(station_cline_depths.iloc[ind['ind'][0]][['CYCLE_NUMBER']])
   # annotations.xy = pos
    #annotations.set_text(text)
#def hover(event):
 #   vis = annotations.get_visible()
  #  if event.inaxes == ax:
   #     cont, ind = ax.contains(event)
    #    if cont:
     #       update_annot(ind)
      #      annotations.set_visible(True)
       #     fig.canvas.draw_idle()
        #else:
         #   if vis:
          #      annotations.set_visible(False)
           #     fig.canvas.draw_idle()
# mplcursors.cursor(ax, hover=True)

# Add axis labels and a legend
ax.set_xlabel('Date')
ax.set_ylabel('Depth [m]')
plt.title(' Thermocline depth [m]')
#ax.legend()
ax.set_ylim(-100, 1400)
# Show the plot
plt.show()
# Reverse the y-axis
ax.invert_yaxis()

# compute a new data set in the same format that the one produce during the PCA on biology concentration in order to perform a mantel test on it


bgc_angola_mantel = bgc_angola

bgc_angola_mantel['layer'] = np.where(bgc_angola_mantel.depth < 100, "0-100m",
                      np.where(np.logical_and(bgc_angola_mantel.depth >= 100, bgc_angola_mantel.depth < 300), "100-300m",
                      np.where(np.logical_and(bgc_angola_mantel.depth >= 300, bgc_angola_mantel.depth < 600), "300-600m",
                      np.where(np.logical_and(bgc_angola_mantel.depth >= 600, bgc_angola_mantel.depth < 1000), "600-1000m",
                      np.where(bgc_angola_mantel.depth > 1000, "1000-2000m","else")))))

bgc_angola_mantel = bgc_angola_mantel[['CYCLE_NUMBER', 'layer', 'TEMP', 'PSAL', 'CHLA_ADJUSTED', 'DOXY_ADJUSTED', 'BBP700', 'dens']]

bgc_angola_mantel = bgc_angola_mantel.groupby(['CYCLE_NUMBER', 'layer']).agg({'TEMP': 'mean', 'PSAL': 'mean', 'CHLA_ADJUSTED': 'mean',
                                                                                  'DOXY_ADJUSTED' : 'mean', 'BBP700':'mean', 'dens' : 'mean'}).reset_index()

bgc_angola_mantel = bgc_angola_mantel[bgc_angola_mantel['layer'] == '0-100m']
bgc_angola_mantel = bgc_angola_mantel[['CYCLE_NUMBER', 'TEMP', 'PSAL', 'CHLA_ADJUSTED', 'DOXY_ADJUSTED', 'BBP700', 'dens']]
bgc_angola_mantel.rename(columns={'TEMP':'TEMP_100m',
                                     'PSAL':'PSAL_100m',
                                     'CHLA_ADJUSTED':'CHLA_ADJUSTED_100m',
                                     'DOXY_ADJUSTED':'DOXY_ADJUSTED_100m',
                                     'BBP700':'BBP700_100m',
                                     'dens':'dens_100m'}, inplace=True)

# save it
bgc_angola_mantel.to_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/physical_data_100m.csv', index=False)