# Script to make an animated map with satellite data and float trajectory


import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime

from matplotlib import pyplot as plt, animation



# define the output directory
output_dir = 'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/CHLA_satellite_data'

# import profiles location
profile = pd.read_csv("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/list_of_profiles_angola.csv", sep=';', encoding='unicode_escape')
profile.rename(columns={'date':'date_ancien'}, inplace=True)
profile.rename(columns={'dateok':'date'}, inplace=True)
datetime64_list = []
for float_date in profile['date']:
    # Convert to datetime object
    date_str = str(int(float_date))
    date_obj = datetime.strptime(date_str, '%Y%m%d')

    # Convert to datetime64[s] format

    datetime64 = np.datetime64(date_obj)

    datetime64_list.append(datetime64)

# Add the 'datetime64' values to the DataFrame
profile['date'] = datetime64_list
profile['date'] = profile['date'].dt.date

date_list = list(profile['date'].unique())

# open the netCDF dataset with chlorophyll values
data = xr.open_dataset('https://my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D',
                       engine='pydap').sel(
    lat=slice(-7, -17),  # first the highest latitude and then the lowest
    lon=slice(7, 17),  # opposite for longitude
    time=date_list  # period of interest
)

# iterate over the time dimension and save each time step to a separate file




dirname ="C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6"


# extract time variable as np array
time_var = np.array(data.variables['time'])

# convert the time values to a list of formatted date strings
time_list = [dt.astype('datetime64[s]').astype('O').strftime('%Y-%m-%d') for dt in time_var]


def create_movie(data, profile, dirname):
    # Get a handle on the figure and the axes
    fig, ax = plt.subplots(figsize=(9,12))

    # Plot the initial frame.
    cax = data.CHL.isel(time=0).plot(
        add_colorbar=True,
        vmin=0, vmax=data.CHL.quantile(0.98),
        cbar_kwargs={
            'extend': 'max'
        }
    )

    # Plot the initial position of the boat.
    boat_pos = profile.iloc[0]
    boat_marker, = ax.plot(
        boat_pos['lon'], boat_pos['lat'], 'kd', markersize=10
    )

    # Create a line for the boat trajectory.
    boat_traj, = ax.plot([], [], 'k--')

    # Next we need to create a function that updates the values for the colormesh, as well as the title.
    def animate(frame):
        # Update the chlorophyll concentration.
        cax.set_array(data.CHL.sel(time=data.time.values[frame]).values.flatten())
        ax.set_title("Time = " + str(data.time.values[frame])[:10])

        # Update the boat position and trajectory.
        boat_pos = profile.iloc[frame]
        boat_marker.set_data(boat_pos['lon'], boat_pos['lat'])
        boat_traj.set_data(profile['lon'][:frame+1], profile['lat'][:frame+1])

    # Finally, we use the animation module to create the animation.
    ani = animation.FuncAnimation(
        fig,             # figure
        animate,         # name of the function above
        frames=range(1, 117),  # Number of frames       # NOMBRE DE PAGES D'IMAGES
        interval=200     # ms between frames
    )

    # Save the animation as mp4 video.
    writer = animation.FFMpegWriter(fps=1)
    ani.save("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/animation_chla_anim.mp4", writer=writer)
create_movie(data, profile, dirname)


'''
# A partir d'ici ça ne marche pas pour l'instant

# open the netCDF dataset with SST values
data = xr.open_dataset('https://nrt.cmems-du.eu/thredds/dodsC/METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2',
                       engine='pydap').sel(
    lat=slice(-30, -45),  # first the highest latitude and then the lowest
    lon=slice(5, 20)).sel( # opposite for longitude
    time=date_list,
    method = 'nearest'# period of interest
)

# iterate over the time dimension and save each time step to a separate file





dirname ="C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6"


# extract time variable as np array
time_var = np.array(data.variables['time'])

# convert the time values to a list of formatted date strings
time_list = [dt.astype('datetime64[s]').astype('O').strftime('%Y-%m-%d') for dt in time_var]


def create_movie(data, profile, dirname):
    # Get a handle on the figure and the axes
    fig, ax = plt.subplots(figsize=(9,12))

    # Plot the initial frame.
    cax = data.analysed_sst.isel(time=0).plot(
        add_colorbar=True,
        vmin=0, vmax=data.analysed_sst.max(),
        cbar_kwargs={
            'extend': 'max'
        }
    )

    # Plot the initial position of the boat.
    boat_pos = profile.iloc[0]
    boat_marker, = ax.plot(
        boat_pos['lon'], boat_pos['lat'], 'kd', markersize=10
    )

    # Create a line for the boat trajectory.
    boat_traj, = ax.plot([], [], 'k--')

    # Next we need to create a function that updates the values for the colormesh, as well as the title.
    def animate(frame):
        # Update the chlorophyll concentration.
        cax.set_array(data.analysed_sst.sel(time=data.time.values[frame]).values.flatten())
        ax.set_title("Time = " + str(data.time.values[frame])[:10])

        # Update the boat position and trajectory.
        boat_pos = profile.iloc[frame]
        boat_marker.set_data(boat_pos['lon'], boat_pos['lat'])
        boat_traj.set_data(profile['lon'][:frame+1], profile['lat'][:frame+1])

    # Finally, we use the animation module to create the animation.
    ani = animation.FuncAnimation(
        fig,             # figure
        animate,         # name of the function above
        frames=range(1, 180),  # Number of frames
        interval=200     # ms between frames
    )

    # Save the animation as mp4 video.
    writer = animation.FFMpegWriter(fps=1)
    ani.save("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/SST_anim.mp4", writer=writer)
create_movie(data, profile, dirname)                  '''
