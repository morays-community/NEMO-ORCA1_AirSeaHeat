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


def main(filepath_ref, filepath, var_name, fig_name, infos, freq):

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
    var_ref = getattr(ds_ref,var_name).values[-1]
    var = getattr(ds,var_name).values[-1]
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
        config = 'ORCA1-AirSeaHeat.S24'


    # QNS difference
    infos = [ 'Non-solar heat flux : ORCA1.S24 - ORCA1 (W/m2)' , cmocean.cm.balance , colors.Normalize(vmin=-100.0, vmax=100.0), lambda x: x ]
    main( filepath_ref='ORCA1-AirSeaHeat_5d_20100101_20100130_grid_T.nc' , filepath=config+'_5d_20100101_20100130_grid_T.nc' , var_name='qns_oce', fig_name='qns_20100126', infos=infos , freq='5d' )

