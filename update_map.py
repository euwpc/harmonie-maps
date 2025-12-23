import requests
import xml.etree.ElementTree as ET
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap, Normalize
import matplotlib
import datetime
import os
from PIL import Image  # For creating animated GIF

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

# --- Step 2: Download data ---
download_url = "https://opendata.fmi.fi/download?producer=harmonie_scandinavia_surface&param=temperature&format=netcdf&bbox=19,56,30,61&projection=EPSG:4326"
response = requests.get(download_url, timeout=300)
response.raise_for_status()
nc_path = "harmonie.nc"
with open(nc_path, "wb") as f:
    f.write(response.content)

# --- Step 3: Load data ---
ds = xr.open_dataset(nc_path)
temp_k = ds['air_temperature_4']
temp_c = temp_k - 273.15

# --- Step 4: Parse high-res color ramp ---
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
colors = [i[1] for i in items]
cmap = ListedColormap(colors)
norm = Normalize(vmin=-40, vmax=50)

# --- Step 5: Generate main map (analysis) + contour labels ---
fig = plt.figure(figsize=(12, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
analysis = temp_c.isel(time=0)
min_temp = float(analysis.min())
max_temp = float(analysis.max())

# Filled contour plot
cf = analysis.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm, levels=100, add_colorbar=True,
                            cbar_kwargs={'label': 'Temperature (°C)', 'shrink': 0.8})

# Contour lines + labels (the numbers you want)
cl = analysis.plot.contour(ax=ax, transform=ccrs.PlateCarree(), colors='black', linewidths=0.5, levels=range(-40, 51, 2))
ax.clabel(cl, inline=True, fontsize=8, fmt="%d")  # Labels every 2°C

ax.coastlines(resolution='10m')
ax.gridlines(draw_labels=True)
ax.set_extent([19, 30, 56, 61])
plt.title(f"HARMONIE 2m Temperature (°C)\nModel run: {run_time_str} | Analysis\nMin: {min_temp:.1f}°C | Max: {max_temp:.1f}°C")
plt.savefig("map.png", dpi=200, bbox_inches='tight')
plt.close()

# --- Step 6: Generate animation frames (one per forecast hour) ---
frames = []
forecast_hours = len(temp_c.time)  # Usually ~67
for i in range(forecast_hours):
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    hour_offset = i  # hours from analysis
    temp_slice = temp_c.isel(time=i)
    
    temp_slice.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm, levels=100)
    temp_slice.plot.contour(ax=ax, transform=ccrs.PlateCarree(), colors='black', linewidths=0.5, levels=range(-40, 51, 4))
    
    ax.coastlines(resolution='10m')
    ax.set_extent([19, 30, 56, 61])
    plt.title(f"+{hour_offset}h | Run: {run_time_str}")
    
    frame_path = f"frame_{i:03d}.png"
    plt.savefig(frame_path, dpi=150, bbox_inches='tight')
    plt.close()
    frames.append(Image.open(frame_path))

# Save animated GIF (loop forever, 500ms per frame)
frames[0].save("animation.gif", save_all=True, append_images=frames[1:], duration=500, loop=0)

print("Main map + contour labels + animation generated")
