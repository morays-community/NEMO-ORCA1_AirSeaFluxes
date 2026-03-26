import os
import argparse
import numpy as np
import xarray as xr
import cmocean

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import matplotlib
matplotlib.use('Agg')

def make_plot(data,lon,lat,infos,output):
    # args
    title, cmap, norm, tfs = infos
    data = tfs(data)
    # figure
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.EqualEarth())
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
    # color map
    pcm = ax.pcolormesh(lon, lat, data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    cbar = plt.colorbar(pcm, ax=ax, orientation='vertical', pad=0.05, shrink=0.5)
    plt.title(title)
    # write fig
    plt.savefig(output, bbox_inches='tight')
    plt.close()


def main(filepath_ref, filepath, var_name, fig_name, infos, freq, time=-1):

    # read files
    try:
        ds = xr.open_dataset(filepath)
        ds_ref = xr.open_dataset(filepath_ref)
    except:
        return

    # coordinates
    lon = ds.nav_lon.values
    lat = ds.nav_lat.values

    # get fields
    if time == -1:
        var_ref = getattr(ds_ref,var_name).values
        var = getattr(ds,var_name).values
    else:
        var_ref = getattr(ds_ref,var_name).values[time,:,:]
        var = getattr(ds,var_name).values[time,:,:]
    diff_var = var - var_ref

    # plot
    plotpath = fig_name + '_diff_' + config +'_' + freq + '.png'
    make_plot(diff_var,lon,lat,infos,plotpath)



if __name__=="__main__":

    # Config name
    # -----------
    try:
        namelist = nml.read('namelist_cfg')
        config = namelist['namrun']['cn_exp']
    except:
        config = 'ORCA1_W25_DET'

    # SST difference
    infos = [ 'SST (ºC): ORCA1.W25_DET - ORCA1' , cmocean.cm.balance , colors.Normalize(vmin=-2.0, vmax=2.0), lambda x: x ]
    main( filepath_ref='../../../../ORCA1_REF/EXP00/OUT/RUN1/ORCA1_REF_40y_ocebudget.nc' , filepath=config+'_40y_ocebudget.nc' , var_name='tos', fig_name='SST', infos=infos , freq='40y' )

    infos = [ 'SST (ºC): ORCA1.W25_DET - ORCA1' , cmocean.cm.balance , colors.Normalize(vmin=-2.0, vmax=2.0), lambda x: x ]
    main( filepath_ref='../../../../ORCA1_REF/EXP00/OUT/RUN1/ORCA1_REF_10y_ocebudget.nc' , filepath=config+'_10y_ocebudget.nc' , var_name='tos', fig_name='SST', infos=infos , freq='10y' )

    # QNS difference
    infos = [ 'Non-solar heat flux (W/m2): ORCA1.W25_DET - ORCA1' , cmocean.cm.balance , colors.Normalize(vmin=-40.0, vmax=40.0), lambda x: x ]
    main( filepath_ref='../../../../ORCA1_REF/EXP00/OUT/RUN1/ORCA1_REF_40y_ocebudget.nc' , filepath=config+'_40y_ocebudget.nc' , var_name='qns', fig_name='QNS', infos=infos , freq='40y' )

    infos = [ 'Non-solar heat flux (W/m2): ORCA1.W25_DET - ORCA1' , cmocean.cm.balance , colors.Normalize(vmin=-40.0, vmax=40.0), lambda x: x ]
    main( filepath_ref='../../../../ORCA1_REF/EXP00/OUT/RUN1/ORCA1_REF_10y_ocebudget.nc' , filepath=config+'_10y_ocebudget.nc' , var_name='qns', fig_name='QNS', infos=infos , freq='10y' )

    # MLD difference - March
    infos = [ 'MLD_03 (m): ORCA1.W25_DET - ORCA1' , cmocean.cm.balance , colors.Normalize(vmin=-50.0, vmax=50.0), lambda x: x ]
    main( filepath_ref='../../../../ORCA1_REF/EXP00/OUT/RUN1/ORCA1_REF_1m40y_averaged_grid_T.nc' , filepath=config+'_1m40y_averaged_grid_T.nc' , var_name='mldr10_1', fig_name='MLD03', infos=infos , freq='40y' , time=2 )

    infos = [ 'MLD_03 (m): ORCA1.W25_DET - ORCA1' , cmocean.cm.balance , colors.Normalize(vmin=-50.0, vmax=50.0), lambda x: x ]
    main( filepath_ref='../../../../ORCA1_REF/EXP00/OUT/RUN1/ORCA1_REF_1m10y_averaged_grid_T.nc' , filepath=config+'_1m10y_averaged_grid_T.nc' , var_name='mldr10_1', fig_name='MLD03', infos=infos , freq='10y' , time=2 )

    # MLD difference - september
    infos = [ 'MLD_09 (m): ORCA1.W25_DET - ORCA1' , cmocean.cm.balance , colors.Normalize(vmin=-50.0, vmax=50.0), lambda x: x ]
    main( filepath_ref='../../../../ORCA1_REF/EXP00/OUT/RUN1/ORCA1_REF_1m40y_averaged_grid_T.nc' , filepath=config+'_1m40y_averaged_grid_T.nc' , var_name='mldr10_1', fig_name='MLD09', infos=infos , freq='40y' , time=8 )

    infos = [ 'MLD_09 (m): ORCA1.W25_DET - ORCA1' , cmocean.cm.balance , colors.Normalize(vmin=-50.0, vmax=50.0), lambda x: x ]
    main( filepath_ref='../../../../ORCA1_REF/EXP00/OUT/RUN1/ORCA1_REF_1m10y_averaged_grid_T.nc' , filepath=config+'_1m10y_averaged_grid_T.nc' , var_name='mldr10_1', fig_name='MLD09', infos=infos , freq='10y' , time=8 )
