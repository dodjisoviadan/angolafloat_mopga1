# Script to create a table from the Ecotaxa project uvp6_WMO6904139_WMO6903095_merged

# 1. Setup

# import needed libraries
import os
from pathlib import Path
# for data management
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
pd.options.display.max_columns = 100
import numpy as np
import math
from numpy import unique
import itertools
from math import floor
from tabulate import tabulate
import timer
import time
import pyarrow
# for plots
import seaborn as sns
import plotly.express as px
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# set up working space
path_to_data = Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6').expanduser()
os.chdir(str(path_to_data))
path_to_figures = Path('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot').expanduser()

# 2. Data preparation

# 2.1 Ecopart

# Angola :
ecopart = pd.read_table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/Ecopart_diagnostics_data_647.tsv", sep='\t', encoding='unicode_escape')
# print(ecopart)
print(ecopart.columns)

# uvp was not working during the last profile (profile name : 0184a_WMO6903095) and
# ctd was not working during the 157th profile (profile name : 0156a_WMO6903095)

# rename some of the columns
ecopart.rename(columns={'Depth [m]': 'depth', 'Vol [L] (sampled for this depth bin)': 'volume_L', 'Profile': 'profile'}, inplace=True)

# select only ascent profiles :
ecopart = ecopart[ecopart['profile'].str.contains('a')]

# compute max and min depth for the profiles
profiles_selected = ecopart.groupby("profile").agg({"depth": [np.max, np.min]})
profiles_selected.columns = ['max_depth', "min_depth"]
profiles_selected = profiles_selected.reset_index()
profiles_selected["range"] = profiles_selected["max_depth"] - profiles_selected["min_depth"]  # depth range of the layer

# look at their distribution
# max depth
plt.figure(figsize=(10, 10), dpi=600)
fig_depth_max = sns.displot(profiles_selected, x="max_depth", kde=True, binwidth=50)
plt.show()
# fig_depth_max.savefig(str(path_to_figures) + '/Distribution_max_depth.pdf')

# min depth
plt.figure(figsize=(10, 10), dpi=600)
fig_depth_min = sns.displot(profiles_selected, x="min_depth", kde=True, binwidth=50)
plt.show()
# fig_depth_min.savefig(str(path_to_figures) + '/Distribution_min_depth.pdf')

# add the lon lat and date to it
profiles_selected_list = list(profiles_selected["profile"].unique())
len(profiles_selected_list)

# Keep in the ecopart table only the selected profiles
ecopart_selected= ecopart.merge(profiles_selected, on=['profile'])
ecopart_selected1 = ecopart_selected.filter(items=['profile', 'Pressure [dbar]', 'volume_L'])
ecopart_selected1.rename(columns={'Pressure [dbar]': 'depth'}, inplace=True)
print(ecopart_selected1)
print(len(np.unique(ecopart_selected1['profile'])))

# save the table
volumes = ecopart_selected1
volumes.to_csv(str(path_to_data) + '/volumes_angola.csv', index=False)
volumes = pd.read_csv(str(path_to_data) + '/volumes_angola.csv')

# 2.2 Ecotaxa

# Import ecotaxa table
# Open and prepare the data
start = time.time()
ecotaxa_angola = pd.read_csv('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/ecotaxa_export_Angola_20230323.tsv',sep='\t')
end = time.time()
print(end - start)

ecotaxa_angola.head()
ecotaxa_angola.shape
ecotaxa_angola.object_depth_min.describe()

print("Max lat:", np.max(ecotaxa_angola['object_lat']))
print("min lat:", np.min(ecotaxa_angola['object_lat']))

print("Max lon:", np.max(ecotaxa_angola['object_lon']))
print("min lon:", np.min(ecotaxa_angola['object_lon']))

# remove "object_" from the columns names
ecotaxa_angola.columns = ecotaxa_angola.columns.str.replace('object_', '')

# keep only relevant information :
ecotaxa_angola_df = ecotaxa_angola[['id', 'lat', 'lon', 'date', 'time', 'depth_min', 'depth_max', 'sample_id', 'annotation_category',
     'annotation_hierarchy', 'area', 'minor', 'major', 'acq_pixel', 'perim.', 'circ.', 'mean', 'kurt', 'fractal']]
print(ecotaxa_angola_df)

# Add 2.5m to the depth columns to be able to compare the ecotaxa data with the ecopart and keep data above 500m
ecotaxa_angola_df['depth_min'] = ecotaxa_angola_df['depth_min'] + 2.5
ecotaxa_angola_df['depth_max'] = ecotaxa_angola_df['depth_max'] + 2.5

# rename columns
ecotaxa_angola_df = ecotaxa_angola_df.rename(columns={"sample_id": "profile", "id": "object_id", "annotation_category": "taxon", "acq_pixel": "px_size",
             "area": "area_px"})

# select only ascent profiles :
ecotaxa_angola_df = ecotaxa_angola_df[ecotaxa_angola_df['profile'].str.contains('a')]

# get the metadata of the profiles
coordinates_profiles = ecotaxa_angola_df.filter(items=['profile', 'lat', 'lon', 'date'])

# 2.2.2 Add the ecopart data and keep only the selected profiles

# add the ecotaxa lon and lat to the profiles we selected earlier
profiles = pd.DataFrame(np.unique(ecopart_selected1.profile), columns=['profile'])
# add the coordinates and date of each profile
profiles = profiles.merge(coordinates_profiles, on=['profile'])

profiles.loc[profiles['date'] == 20210731, 'date'] = 20210801 # to figured out the problem in the next agg function because of different date for the same profil
# which create a gap in the time serie
# for angola data the problem

# keep only one line per profile and per date
list_prof = profiles.groupby(['profile']).agg({'lat': 'mean', 'lon': 'mean', 'date': 'mean'}).reset_index() #j'ai enlever le groupby date

print(len(np.unique(list_prof['profile'])))

# check that it was done properly
print(len(np.unique(ecopart_selected1.profile)))
print(len(np.unique(list_prof.profile)))
ecopart_selected_unique = ecopart_selected1.groupby(['profile']).size().reset_index(name='counts')
list_prof_unique = list_prof.groupby(['profile']).size().reset_index(name='counts')

# check this selection to see if the same profiles are found here
df_1notin2 = ecopart_selected_unique[~(ecopart_selected_unique['profile'].isin(list_prof_unique['profile']))].reset_index(drop=True)
# this shows profiles which are in ecopart_selected_unique (filtered ecopart profiles)
# but not in list_prof_unique (profiles from ecotaxa which are among the filtered ecopart profiles)
len(df_1notin2) == len(ecopart_selected_unique) - len(list_prof_unique)  # True
print(df_1notin2)  # it's empty so we have the same profiles in both tables. Great !

# save the table for later
list_prof.to_csv(str(path_to_data) + '/list_of_profiles_angola.csv', index=False)

# select in ecotaxa only the profiles which are present in list_prof
list_prof_reduced = list_prof.filter(items=['profile'])  # keep only the profile column

# this will allow us to not have two date, lon and lat columns as they are not exactly the same
# between CTD and UVP data
ecotaxa_angola_before_taxo = ecotaxa_angola_df.merge(list_prof_reduced, on=['profile'])


# save
ecotaxa_angola_df.reset_index().to_csv(str(path_to_data) + '/ecotaxa_angola_df.csv')


# next step do the UVP regrouped script (Rstudio)

UVP5_taxo_regrouped = pd.read_csv(str(path_to_data) + "/Angola_level_taxo.csv", sep = ';')
#print(UVP5_taxo_regrouped)
UVP5_taxo_regrouped = UVP5_taxo_regrouped.filter(items=['angola_rough', 'angola_medium', 'taxon'])

# remove the lines of elements which are not regrouped by angola
UVP5_taxo_regrouped = UVP5_taxo_regrouped[UVP5_taxo_regrouped['taxon'].notna()]
UVP5_taxo_regrouped = UVP5_taxo_regrouped[UVP5_taxo_regrouped['angola_rough'].notna()]
UVP5_taxo_regrouped.head()

# merge the ecotaxa dataframe with the regrouping taxa table
indiv = pd.merge(ecotaxa_angola_before_taxo, UVP5_taxo_regrouped, how="left", on="taxon")

# remove the elements which were not regrouped
indiv = indiv[indiv['angola_rough'].notna()]

print(indiv.shape)

# count the images in each category
count = indiv.groupby(['angola_rough']).agg(n = pd.NamedAgg(column = "angola_rough", aggfunc = 'count')).reset_index()

# keep only the profiles from "profiles_selected"
indiv = indiv.merge(profiles_selected, on=['profile'])

# remove columns that we don't need anymore
indiv.drop(['max_depth', 'min_depth', 'range'], inplace=True, axis=1)

indiv.head(5)

# Plot the distribution of taxonomic groups

    # With barplots

## Plot
# the counts of each group
plt.clf()
plt.figure(figsize=(25, 15), dpi = 600)
sns.countplot(y="taxon", data=indiv, order = indiv['taxon'].value_counts().index)
plt.tick_params(axis='both', which='major', labelsize=6)
# save it
#plt.savefig(str(path_to_figures) + '/taxon.pdf', dpi= 600)
# let's see it
plt.show()

# angola rough regrouping
plt.figure(figsize=(4, 2), dpi = 600)
sns.countplot(y="angola_rough", data=indiv, order = indiv['angola_rough'].value_counts().index)
plt.tick_params(axis='both', which='major', labelsize=3)
# save it
#plt.savefig(str(path_to_figures) + '/exports_groups_rough_regrouping.pdf', dpi= 600)
# let's see it
plt.show()

# angola medium regrouping
plt.figure(figsize=(4, 2), dpi = 600)
plt.clf()
sns.countplot(y="angola_medium", data=indiv, order = indiv['angola_medium'].value_counts().index)
plt.tick_params(axis='both', which='major', labelsize=2.5)
# save it
#plt.savefig(str(path_to_figures) + '/exports_groups_medium_regrouping.pdf', dpi= 600)
# let's see it
plt.show()

    # With lon/lat distribution scatter plots
indiv.head()

# The distribution of the counts according to longitude and latitude
cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
indiv = indiv.dropna()
taxa_list_rough = indiv.angola_rough.unique()
taxa_list_medium = indiv.angola_medium.unique()

# remove the artefacts
taxa_list_rough = taxa_list_rough[taxa_list_rough != 'Artefact']
taxa_list_medium = taxa_list_medium[taxa_list_medium != 'Artefact']

taxa_list_medium
taxa_list_rough

# Choose the definition either as rough or medium
definition_needed = ['angola_rough', 'angola_medium']
number = 0
def_chosen = definition_needed[number]
def_chosen

if def_chosen == 'angola_rough':
    taxa_list = taxa_list_rough
else:
    taxa_list = taxa_list_medium
taxa_list

# Open a pdf to save the figures rough definition
pp = PdfPages(str(path_to_figures) + '/distribution_number_taxa_' + def_chosen + '.pdf')

# sort the groups
taxa_list = sorted(taxa_list)

#Run dans la console puis run les lignes "df_temp = indiv[indiv[def_chosen] == t]" à "plt.show()"

 ''' for t in taxa_list:
    # select the taxa of interest
    df_temp = indiv[indiv[def_chosen] == t]
    # for every location, count the number of organisms from this taxa
    df_temp = df_temp.groupby(['lon', 'lat', 'date']).agg(n=pd.NamedAgg(column=def_chosen, aggfunc='count'))
    df_temp.columns = ['n']
    # reset the index so that it won't be an issue later
    df_temp = df_temp.reset_index()
    # put the date in the day of May
    #df_temp.date = df_temp.date - 20210500
    # set the style
    sns.set_style("white")

# plot
    plt.figure(figsize=(10, 12), dpi = 600)
    g = sns.relplot(
        data=df_temp,
        x="lon", y="lat", size="n", hue="date",
        palette=cmap, sizes=(10, 200), alpha = 0.7
    ).set(title=t)
    # make the title visible
    g.fig.subplots_adjust(top=.90)
    # limits of the axes
    plt.xlim(min(indiv.lon)-0.5, max(indiv.lon)+0.5)
    plt.ylim(min(indiv.lat)-0.5, max(indiv.lat)+0.5)
    plt.show()
    pp.savefig(plt.gcf())            #Save each figure in the pdf

# close the pdf document
pp.close() '''

# make some conversions to have everything in metrics (mm instead of pixels)
indiv["area"] = indiv["area_px"]*indiv["px_size"]**2
indiv["esd"] = 2*np.sqrt(indiv["area"]/np.pi)
indiv["minor"] = indiv["minor"]*indiv["px_size"]
indiv["major"] = indiv["major"]*indiv["px_size"]
indiv["vol_sph"] = 4/3*np.pi*(indiv["esd"]/2)**3

indiv.head()

# compute bins of 5m depth
indiv["depth"] = indiv.loc[:, ["depth_min","depth_max"]].mean(axis= 1)
# remove depth_min and depth_max
indiv.drop(['depth_min', 'depth_max'], inplace=True, axis=1)

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

indiv.head()

# use list comprehension to create the column depth_bin using the round_to_multiple function
# m93["depth_bin"] = [round_to_multiple(m93["depth"][r], 5) + 2.5 for r in range(len(m93))]
indiv.loc[:,'depth_bin'] = 99999
indiv['depth_bin'] = indiv["depth"].apply(lambda x: rounditup(x, 5, "floor")+ 2.5)

indiv.head()

# Compute abundance and biovolume per taxa and bin for either the rough or medium taxonomic definiton
indiv_binned = indiv.groupby([def_chosen, 'depth_bin', 'profile'], group_keys=False).agg(n = pd.NamedAgg(column = def_chosen, aggfunc = 'count'),
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
obs = pd.merge(indiv_binned, ecopart_selected1, how="left", on=['profile', 'depth'])
obs["watervolume"] = obs["volume_L"] / 1000          # volume in m3
obs["conc"] = obs["n"] / obs["watervolume"]          # concentration in n/m3
obs["vol_sph"] = obs["vol_sph"] / obs["watervolume"] # biovolume concentration in mm3/m3

# keep only the 5m depth bins which have a watervolume <= 1000L
obs = obs[obs["watervolume"] <= 1]
obs = obs.dropna(subset=['watervolume'])

volumes.head()

# compute all taxon x bins combinations
#taxon_list = np.unique(obs['EXPORTS_rough'])
#volume = ecopart.filter(items = ["profile", "Project", "depth", "volume_L"])
#volume.columns = volume.columns.str.replace('depth', 'mid_depth_bin')
# create lists containing all the depth and profiles possibilities
depth_bins = np.unique(volumes.depth)

# and all the profiles
profile_list = np.unique(volumes.profile)
# regroup them with the list of taxa
list = [profile_list, depth_bins, taxa_list]
# compute the unique combinations of the 3 lists
combination  = [p for p in itertools.product(*list)]

# Convert it into a dataframe
column_names = ["profile", "depth", def_chosen]
all = pd.DataFrame(combination, columns =["profile", "depth", def_chosen])

# Add the data from obs
full = pd.merge(all, obs, how="left", on=['profile', 'depth', def_chosen])
full.head()

full = pd.DataFrame(full)

# remove the volume_L column. We'll keep only watervolume (in m3)
full.drop('volume_L', inplace = True, axis = 1)

# rename the taxon column
full.rename(columns = {def_chosen:'taxon'}, inplace = True)

# consider the observations we don't have as zeroes
cols_to_fill = ['n', 'conc', 'vol_sph', 'watervolume']
full[cols_to_fill]=full[cols_to_fill].fillna(0)

# Check the histograms for anomalies
plt.hist(full['vol_sph'], bins=100, edgecolor='black')
plt.show()

# do some stats
stats = [['0 in %', sum(i < pow(10,-3) for i in full.vol_sph)/ len(full) * 100],
         ['q25',np.quantile(full.vol_sph, 0.25)],
         ['med',np.quantile(full.vol_sph, 0.5)],
         ['q75',np.quantile(full.vol_sph, 0.75)],
         ['q95',np.quantile(full.vol_sph, 0.95)],
         ['mean', np.mean(full.vol_sph)]]
col_names = ["stat", "value"]
print(tabulate(stats, headers=col_names))

# and look at the number of vignettes per category
toto = full[full['n'] > 0]
toto = toto.groupby('taxon').agg({'n': 'sum'}) #, 'vol_sph': 'sum'})
toto['percent'] = toto['n'] / sum(toto['n'] )* 100
toto = toto.sort_values('percent', ascending=False).reset_index()

# plot count
plt.clf()
plt.figure(figsize=(2, 1), dpi = 600)
sns.barplot(data = toto, y = 'taxon', x = 'n')
plt.tick_params(axis='both', which='major', labelsize=2)
#plt.savefig(str(path_to_figures) + '/barplots_taxa_numbers_without_aggregates_' + def_chosen + '.pdf')
plt.show()

# plot percentage
sns.barplot(data = toto, y = 'taxon', x = 'percent')
plt.tick_params(axis='both', which='major', labelsize=8)
plt.show()

# Reduce to volume and number per sample, taxon and layer
reduced = full.merge(list_prof, on=['profile'])

# Plot for each taxa a map with size of point = size of volume or number of images
# create a pdf
pp = PdfPages(str(path_to_figures) + '/distribution_number_taxa_' + def_chosen + '.pdf')

taxa_list = sorted(taxa_list)

for t in taxa_list:
    df_temp = reduced[reduced['taxon'] == t]
    df_temp = df_temp[df_temp['n'] > 0]
    # reset the index so that it won't be an issue later
    df_temp = df_temp.reset_index()
    # for every location, count the number of organisms from this taxa
    df_temp = df_temp.groupby(['lon', 'lat']).agg(n=pd.NamedAgg(column="taxon", aggfunc='count'))
    df_temp.columns = ['n']
    df_temp = df_temp.reset_index()

    # plot
    plt.figure(figsize=(10, 10), dpi=600)
    g = sns.relplot(
        data=df_temp,
        x="lon", y="lat", size="n", sizes=(10, 200),
        palette=cmap).set(title=t)
    # make the title visible
    g.fig.subplots_adjust(top=.90)
    # limits of the axes
    plt.xlim(min(indiv.lon) - 0.5, max(indiv.lon) + 0.5)
    plt.ylim(min(indiv.lat) - 0.5, max(indiv.lat) + 0.5)
    plt.show()
    pp.savefig(plt.gcf())  # Save each figure in the pdf

pp.close()

# do the same for the concentration
pp = PdfPages(str(path_to_figures) + '/distribution_concentration_n_per_m3_taxa_' + def_chosen + '.pdf')

taxa_list = sorted(taxa_list)

for t in taxa_list:
    # select the taxa of interest
    df_temp = reduced[reduced['taxon'] == t]
    df_temp = df_temp[df_temp['n'] > 0]
    # reset the index so that it won't be an issue later
    df_temp = df_temp.reset_index()
    # for every location, count the number of organisms from this taxa
    df_temp = df_temp.groupby(['lon', 'lat']).agg(n=pd.NamedAgg(column="taxon", aggfunc='count'))
    df_temp.columns = ['conc']
    df_temp = df_temp.reset_index()

    # plot
    plt.figure(figsize=(10, 10), dpi=600)
    g = sns.relplot(
        data=df_temp,
        x="lon", y="lat", size="conc", sizes=(10, 200),
        palette=cmap).set(title=t)  # , sizes=(10, 200),
    # make the title visible
    g.fig.subplots_adjust(top=.90)
    # limits of the axes
    plt.xlim(min(reduced.lon) - 0.5, max(reduced.lon) + 0.5)
    plt.ylim(min(reduced.lat) - 0.5, max(reduced.lat) + 0.5)
    plt.show()
    pp.savefig(plt.gcf())  # Save each figure in pdf

pp.close()

# save table

indiv.to_csv(str(path_to_data) + '/angola_indiv_' + def_chosen + '.csv', index=False)

full.to_csv(str(path_to_data) + '/angola_full_' + def_chosen + '.csv', index=False)

reduced.to_csv(str(path_to_data) + '/angola_treated_' + def_chosen +'.csv', index=False)

plt.close('all')