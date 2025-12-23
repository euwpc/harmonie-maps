import requests
import xml.etree.ElementTree as ET
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap, Normalize, LogNorm
import matplotlib
import datetime
import os
import glob
from PIL import Image
import pandas as pd

matplotlib.use('Agg')

# --- Step 1: Latest model run ---
wfs_url = "https://opendata.fmi.fi/wfs?service=WFS&version=2.0.0&request=getFeature&storedquery_id=fmi::forecast::harmonie::surface::grid"
response = requests.get(wfs_url, timeout=60)
response.raise_for_status()
tree = ET.fromstring(response.content)

ns = {'gml': 'http://www.opengis.net/gml/3.2', 'omso': 'http://inspire.ec.europa.eu/schemas/omso/3.0'}
origintimes = [elem.text for elem in tree.findall('.//omso:phenomenonTime//gml:beginPosition', ns)] or \
              [elem.text for elem in tree.findall('.//gml:beginPosition', ns)]
latest_origintime = max(origintimes)
run_time_str = datetime.datetime.strptime(latest_origintime, "%Y-%m-%dT%H:%M:%SZ").strftime("%Y-%m-%d %H:%M UTC")

# --- Step 2: Download with wide bbox ---
download_url = (
    "https://opendata.fmi.fi/download?"
    "producer=harmonie_scandinavia_surface&"
    "param=temperature,Dewpoint,Pressure,CAPE,WindGust,Precipitation1h&"
    "format=netcdf&"
    "bbox=10,53,35,71&"
    "projection=EPSG:4326"
)
response = requests.get(download_url, timeout=300)
response.raise_for_status()
nc_path = "harmonie.nc"
with open(nc_path, "wb") as f:
    f.write(response.content)

# --- Step 3: Load data ---
ds = xr.open_dataset(nc_path)

temp_c = ds['air_temperature_4'] - 273.15
dewpoint_c = ds['dew_point_temperature_10'] - 273.15
pressure_hpa = ds['air_pressure_at_sea_level_1'] / 100
cape = ds['atmosphere_specific_convective_available_potential_energy_59']
windgust_ms = ds['wind_speed_of_gust_417']
precip_mm = ds['precipitation_amount_353']

# --- Step 4: High-res temperature colormap ---
tree = ET.parse("temperature_color_table_high.qml")
root = tree.getroot()
items = []
for item in root.findall(".//item"):
    value = float(item.get('value'))
    color_hex = item.get('color').lstrip('#')
    r = int(color_hex[0:2], 16) / 255.0
    g = int(color_hex[2:4], 16) / 255.0
    b = int(color_hex[4:6], 16) / 255.0
    items.append((value, (r, g, b, 1.0)))
items.sort(key=lambda x: x[0])
temp_colors = [i[1] for i in items]
temp_cmap = ListedColormap(temp_colors)
temp_norm = Normalize(vmin=-40, vmax=50)

# Colormaps for other variables
dewpoint_cmap = temp_cmap
dewpoint_norm = Normalize(vmin=-40, vmax=30)

pressure_cmap = plt.cm.viridis_r
pressure_norm = Normalize(vmin=950, vmax=1050)

cape_cmap = plt.cm.YlOrRd
cape_norm = Normalize(vmin=0, vmax=2000)

windgust_cmap = plt.cm.plasma
windgust_norm = Normalize(vmin=0, vmax=25)

precip_cmap = plt.cm.Blues
precip_norm = Normalize(vmin=0, vmax=10)

# --- Step 5: Helper ---
def get_analysis(var):
    if 'time' in var.dims:
        return var.isel(time=0)
    elif 'time_h' in var.dims:
        return var.isel(time_h=0)
    return var

# --- Step 6: Views ---
views = {
    'focused': {'extent': [19, 30, 56, 61], 'suffix': ''},
    'wide':    {'extent': [10, 35, 53, 71], 'suffix': '_wide'}
}

variables = {
    'temperature': {'var': temp_c, 'cmap': temp_cmap, 'norm': temp_norm, 'unit': '°C', 'title': '2m Temperature (°C)', 'levels': range(-40, 51, 2)},
    'dewpoint':    {'var': dewpoint_c, 'cmap': dewpoint_cmap, 'norm': dewpoint_norm, 'unit': '°C', 'title': '2m Dew Point (°C)', 'levels': range(-40, 31, 2)},
    'pressure':    {'var': pressure_hpa, 'cmap': pressure_cmap, 'norm': pressure_norm, 'unit': 'hPa', 'title': 'MSLP (hPa)', 'levels': range(950, 1051, 4)},
    'cape':        {'var': cape, 'cmap': cape_cmap, 'norm': cape_norm, 'unit': 'J/kg', 'title': 'CAPE (J/kg)', 'levels': range(0, 2001, 200)},
    'windgust':    {'var': windgust_ms, 'cmap': windgust_cmap, 'norm': windgust_norm, 'unit': 'm/s', 'title': 'Wind Gust (m/s)', 'levels': range(0, 26, 2)},
    'precip':      {'var': precip_mm, 'cmap': precip_cmap, 'norm': precip_norm, 'unit': 'mm/h', 'title': 'Precipitation (1h)', 'levels': [0, 0.1, 0.5, 1, 2, 5, 10]}
}

# --- Generate for each view ---
for view_key, view_conf in views.items():
    extent = view_conf['extent']
    suffix = view_conf['suffix']

    for var_key, conf in variables.items():
        # Analysis map
        data = get_analysis(conf['var'])
        min_val = float(data.min())
        max_val = float(data.max())
        
        fig = plt.figure(figsize=(14 if view_key == 'wide' else 12, 10))
        ax = plt.axes(projection=ccrs.PlateCarree())
        data.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap=conf['cmap'], norm=conf['norm'], levels=100,
                           cbar_kwargs={'label': conf['unit'], 'shrink': 0.8})
        cl = data.plot.contour(ax=ax, transform=ccrs.PlateCarree(), colors='black', linewidths=0.5, levels=conf['levels'])
        ax.clabel(cl, inline=True, fontsize=8, fmt="%.1f" if var_key == 'precip' else "%d")
        
        ax.coastlines(resolution='10m')
        ax.gridlines(draw_labels=True)
        ax.set_extent(extent)
        plt.title(f"HARMONIE {conf['title']}\nModel run: {run_time_str} | Analysis\nMin: {min_val:.1f} {conf['unit']} | Max: {max_val:.1f} {conf['unit']}")
        plt.savefig(f"{var_key}{suffix}.png", dpi=200, bbox_inches='tight')
        plt.close()

        # Animation with 1h steps up to +48h, then 3h steps after
        frames = []
        time_dim = 'time' if 'time' in conf['var'].dims else 'time_h'
        time_values = ds[time_dim].values
        
        for i in range(len(time_values)):
            # Skip frames after +48h that are not on 3-hour steps
            if i >= 48 and (i - 48) % 3 != 0:
                continue

            fig = plt.figure(figsize=(12 if view_key == 'wide' else 10, 8))
            ax = plt.axes(projection=ccrs.PlateCarree())
            slice_data = conf['var'].isel(**{time_dim: i})
            hour_offset = i

            slice_data.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap=conf['cmap'], norm=conf['norm'], levels=100)
            cl = slice_data.plot.contour(ax=ax, transform=ccrs.PlateCarree(), colors='black', linewidths=0.5, levels=conf['levels'])
            ax.clabel(cl, inline=True, fontsize=8, fmt="%.1f" if var_key == 'precip' else "%d")

            ax.coastlines(resolution='10m')
            ax.gridlines(draw_labels=True)
            ax.set_extent(extent)
            
            valid_dt = pd.to_datetime(time_values[i])
            valid_dt_eet = valid_dt + pd.Timedelta(hours=2)
            valid_str = valid_dt_eet.strftime("%a %d %b %H:%M EET")
            
            plt.title(f"HARMONIE {conf['title']}\nValid: {valid_str} | +{hour_offset}h from run {run_time_str}")

            frame_path = f"frame_{var_key}{suffix}_{i:03d}.png"
            plt.savefig(frame_path, dpi=110, bbox_inches='tight')
            plt.close()
            frames.append(Image.open(frame_path))

        frames[0].save(f"{var_key}{suffix}_animation.gif", save_all=True, append_images=frames[1:], duration=500, loop=0)

        # Clean frames
        for f in glob.glob(f"frame_{var_key}{suffix}_*.png"):
            os.remove(f)

# --- Cleanup large file ---
if os.path.exists("harmonie.nc"):
    os.remove("harmonie.nc")
    print("Removed large harmonie.nc file")

print("All analysis maps + animations for focused and wide views generated and cleanup complete")
