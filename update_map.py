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

# --- Step 2: Download using your original URL + projection for standard lat/lon ---
download_url = (
    "https://opendata.fmi.fi/download?"
    "producer=harmonie_scandinavia_surface&"
    "param=temperature&"
    "format=netcdf&"
    "bbox=19,56,30,61&"
    "projection=EPSG:4326"  # Add this for WGS84 lat/lon grid
)

print(f"Downloading from: {download_url}")
response = requests.get(download_url, timeout=300)
response.raise_for_status()

nc_path = "harmonie.nc"
with open(nc_path, "wb") as f:
    f.write(response.content)

file_size_mb = os.path.getsize(nc_path) / 1024 / 1024
print(f"Downloaded NetCDF ({file_size_mb:.1f} MB)")

# --- Step 3: Load and convert temperature to °C ---
ds = xr.open_dataset(nc_path)
print("Available variables:", list(ds.data_vars))

temp_k = ds['air_temperature_4']
temp_c = temp_k - 273.15

# --- Step 4: Parse your exact color ramp from temperature_style.qml ---
qml_path = "temperature_color_table.qml"
if not os.path.exists(qml_path):
    raise FileNotFoundError("temperature_style.qml not found – add your QGIS style file to the repo")

tree = ET.parse(qml_path)
root = tree.getroot()

items = []
for item in root.findall(".//item"):
    value = float(item.get('value'))
    color_hex = item.get('color')  # e.g., "#ff6eff"
    alpha = int(item.get('alpha', 255)) / 255.0

    color_hex = color_hex.lstrip('#')
    r = int(color_hex[0:2], 16) / 255.0
    g = int(color_hex[2:4], 16) / 255.0
    b = int(color_hex[4:6], 16) / 255.0

    items.append((value, (r, g, b, alpha)))

if not items:
    raise Exception("No color items found in .qml")

items.sort(key=lambda x: x[0])  # Sort by temperature value
values = [i[0] for i in items]
colors = [i[1] for i in items]

cmap = ListedColormap(colors)
norm = Normalize(vmin=-40, vmax=50)  # Matches your QML classificationMin/Max

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

print(f"Map successfully generated with your custom colors: {output_path}")
