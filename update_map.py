import requests
import xml.etree.ElementTree as ET
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap, Normalize
import matplotlib
import datetime
import os
from PIL import Image

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

# --- Step 2: Download all variables (fixed param names) ---
download_url = (
    "https://opendata.fmi.fi/download?"
    "producer=harmonie_scandinavia_surface&"
    "param=temperature,Dewpoint,Pressure,CAPE&"  # Fixed: Dewpoint, Pressure, CAPE
    "format=netcdf&"
    "bbox=19,56,30,61&"
    "projection=EPSG:4326"
)
response = requests.get(download_url, timeout=300)
response.raise_for_status()
nc_path = "harmonie.nc"
with open(nc_path, "wb") as f:
    f.write(response.content)

# --- Step 3: Load data (fixed variable names) ---
ds = xr.open_dataset(nc_path)
temp_c = ds['air_temperature_4'] - 273.15
dewpoint_c = ds['dew_point_temperature_10'] - 273.15  # Fixed name
pressure_hpa = ds['air_pressure_at_sea_level_1'] / 100  # Pa → hPa
cape = ds['atmosphere_specific_convective_available_potential_energy_59']  # J/kg

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

# Colormaps for others (can customize later)
dewpoint_cmap = temp_cmap
dewpoint_norm = Normalize(vmin=-40, vmax=30)

pressure_cmap = plt.cm.viridis_r
pressure_norm = Normalize(vmin=950, vmax=1050)

cape_cmap = plt.cm.YlOrRd
cape_norm = Normalize(vmin=0, vmax=2000)

# --- Step 5: Generate analysis maps ---
variables = {
    'temperature': {'data': temp_c.isel(time=0), 'cmap': temp_cmap, 'norm': temp_norm, 'unit': '°C', 'title': '2m Temperature (°C)', 'levels': range(-40, 51, 2), 'file': 'temperature.png'},
    'dewpoint':    {'data': dewpoint_c.isel(time=0), 'cmap': dewpoint_cmap, 'norm': dewpoint_norm, 'unit': '°C', 'title': '2m Dew Point (°C)', 'levels': range(-40, 31, 2), 'file': 'dewpoint.png'},
    'pressure':    {'data': pressure_hpa.isel(time=0), 'cmap': pressure_cmap, 'norm': pressure_norm, 'unit': 'hPa', 'title': 'MSLP (hPa)', 'levels': range(950, 1051, 4), 'file': 'pressure.png'},
    'cape':        {'data': cape.isel(time=0), 'cmap': cape_cmap, 'norm': cape_norm, 'unit': 'J/kg', 'title': 'CAPE (J/kg)', 'levels': range(0, 2001, 200), 'file': 'cape.png'}
}

for key, conf in variables.items():
    fig = plt.figure(figsize=(12, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    data = conf['data']
    min_val = float(data.min())
    max_val = float(data.max())

    data.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap=conf['cmap'], norm=conf['norm'], levels=100,
                       cbar_kwargs={'label': conf['unit'], 'shrink': 0.8})
    cl = data.plot.contour(ax=ax, transform=ccrs.PlateCarree(), colors='black', linewidths=0.5, levels=conf['levels'])
    ax.clabel(cl, inline=True, fontsize=8, fmt="%d")

    ax.coastlines(resolution='10m')
    ax.gridlines(draw_labels=True)
    ax.set_extent([19, 30, 56, 61])
    plt.title(f"HARMONIE {conf['title']}\nModel run: {run_time_str} | Analysis\nMin: {min_val:.1f} {conf['unit']} | Max: {max_val:.1f} {conf['unit']}")
    plt.savefig(conf['file'], dpi=200, bbox_inches='tight')
    plt.close()

# --- Step 6: Temperature animation (unchanged) ---
frames = []
for i in range(len(temp_c.time)):
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    temp_slice = temp_c.isel(time=i)
    hour_offset = i

    temp_slice.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap=temp_cmap, norm=temp_norm, levels=100)
    cl = temp_slice.plot.contour(ax=ax, transform=ccrs.PlateCarree(), colors='black', linewidths=0.5, levels=range(-40, 51, 2))
    ax.clabel(cl, inline=True, fontsize=8, fmt="%d")

    ax.coastlines(resolution='10m')
    ax.gridlines(draw_labels=True)
    ax.set_extent([19, 30, 56, 61])
    plt.title(f"HARMONIE 2m Temperature (°C)\n+{hour_offset}h | Run: {run_time_str}")

    frame_path = f"frame_{i:03d}.png"
    plt.savefig(frame_path, dpi=150, bbox_inches='tight')
    plt.close()
    frames.append(Image.open(frame_path))

frames[0].save("animation.gif", save_all=True, append_images=frames[1:], duration=500, loop=0)

print("All maps + temperature animation generated")
