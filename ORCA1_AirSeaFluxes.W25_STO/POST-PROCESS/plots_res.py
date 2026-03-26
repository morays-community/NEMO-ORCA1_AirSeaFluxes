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


def main(filepath, var_name, fig_name, infos, freq):

    # read files
    try:
        ds = xr.open_dataset(filepath)
    except:
        return

    print(f'Plotting {var_name}')

    # get fields
    lon = ds.nav_lon.values
    lat = ds.nav_lat.values
    var_val = getattr(ds,var_name).values

    # plot
    plotpath = fig_name + '_' + config +'_' + freq + '.png'
    make_plot(var_val,lon,lat,infos,plotpath)



if __name__=="__main__":

    # Config name
    # -----------
    try:
        namelist = nml.read('namelist_cfg')
        config = namelist['namrun']['cn_exp']
    except:
        config = 'ORCA1_W25_STO'

    # snapshots
    # ---------
    # SST
    infos = [ 'SST (ºC)' , cmocean.cm.thermal , colors.Normalize(vmin=-5.0, vmax=35.0), lambda x: x ]
    main( filepath=config+'_40y_ocebudget.nc' , var_name='tos' , fig_name='SST', infos=infos , freq='40y' )

    # Heat Fluxes
    infos = [ 'Non-solar Heat flux (W/m2)' , cmocean.cm.balance , colors.Normalize(vmin=-300.0, vmax=50.0), lambda x: x ]
    main( filepath=config+'_40y_ocebudget.nc' , var_name='qns' , fig_name='QNS', infos=infos , freq='40y' )
