import requests
import xml.etree.ElementTree as ET
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap, Normalize
import matplotlib
import datetime
import os

matplotlib.use('Agg')  # Headless mode for GitHub Actions

# --- Step 1: Find the latest available model run (origintime) ---
wfs_url = (
    "https://opendata.fmi.fi/wfs?service=WFS&version=2.0.0&request=getFeature&"
    "storedquery_id=fmi::forecast::harmonie::surface::grid"
)

response = requests.get(wfs_url, timeout=60)
response.raise_for_status()
tree = ET.fromstring(response.content)

ns = {
    'gml': 'http://www.opengis.net/gml/3.2',
    'omso': 'http://inspire.ec.europa.eu/schemas/omso/3.0',
}

origintimes = []
for elem in tree.findall('.//omso:phenomenonTime//gml:beginPosition', ns):
    text = elem.text
    if text:
        origintimes.append(text)

if not origintimes:
    for elem in tree.findall('.//gml:beginPosition', ns):
        text = elem.text
        if text:
            origintimes.append(text)

if not origintimes:
    raise Exception("No model runs found – FMI service may have issues")

latest_origintime = max(origintimes)
print(f"Latest model run (origintime): {latest_origintime}")

# --- Step 2: Build the download URL exactly like your original (defaults to latest run) ---
# Your original URL (with lowercase 'temperature' and no extra params)
download_url = (
    "https://opendata.fmi.fi/download?"
    "producer=harmonie_scandinavia_surface&"
    "param=temperature&"  # lowercase as in your original
    "format=netcdf&"
    "bbox=19,56,30,61"
    # No origintime/starttime/timesteps/levels/projection – let it default to latest full run
)

print(f"Downloading from: {download_url}")
response = requests.get(download_url, timeout=300)  # Longer timeout for larger file
response.raise_for_status()

nc_path = "harmonie.nc"
with open(nc_path, "wb") as f:
    f.write(response.content)

file_size_mb = os.path.getsize(nc_path) / 1024 / 1024
print(f"Downloaded NetCDF ({file_size_mb:.1f} MB)")

# --- Step 3: Load and convert temperature to °C ---
ds = xr.open_dataset(nc_path)
print("Available variables:", list(ds.data_vars))

# Temperature is usually 'temperature' (lowercase) in the NetCDF from this producer
if 'temperature' in ds:
    temp_k = ds['temperature']
elif 'Temperature' in ds:
    temp_k = ds['Temperature']
elif 't2m' in ds:
    temp_k = ds['t2m']
else:
    raise Exception("Temperature variable not found")

temp_c = temp_k - 273.15

# --- Step 4: Load your custom color ramp from the .qml file ---
qml_path = "temperature_color_table.qml"
if not os.path.exists(qml_path):
    raise FileNotFoundError("temperature_style.qml not found – add your QGIS style file")

tree = ET.parse(qml_path)
root = tree.getroot()

items = []
for item in root.findall(".//colorrampshader/colorrampitem"):
    value = float(item.get('value'))
    color_str = item.get('color')
    if color_str.startswith('#'):
        color_str = color_str.lstrip('#')
        lv = len(color_str)
        rgb = tuple(int(color_str[i:i + lv // 3], 16) / 255.0 for i in range(0, lv, lv // 3))
        if len(rgb) == 3:
            rgb += (1.0,)
    else:
        rgb = tuple(int(x) / 255.0 for x in color_str.split(','))
    items.append((value, rgb))

if not items:
    raise Exception("Could not parse color ramp from .qml – check format or export as text from QGIS")

items.sort(key=lambda x: x[0])
values = [i[0] for i in items]
colors = [i[1] for i in items]

cmap = ListedColormap(colors)
norm = Normalize(vmin=min(values), vmax=max(values))

# --- Step 5: Plot the analysis timestep (first time step) ---
fig = plt.figure(figsize=(12, 10))
ax = plt.axes(projection=ccrs.PlateCarree())

temp_c.isel(time=0).plot(
    ax=ax,
    transform=ccrs.PlateCarree(),
    cmap=cmap,
    norm=norm,
    cbar_kwargs={'label': 'Temperature (°C)', 'shrink': 0.8}
)

ax.coastlines(resolution='10m')
ax.gridlines(draw_labels=True)
ax.set_extent([19, 30, 56, 61])
plt.title(f"Latest HARMONIE 2m Temperature (°C)\nAnalysis from latest run")

output_path = "map.png"
plt.savefig(output_path, dpi=200, bbox_inches='tight')
plt.close()

print(f"Map successfully generated: {output_path}")
