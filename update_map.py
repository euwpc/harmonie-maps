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
now = datetime.datetime.utcnow()
# Look back up to 12 hours to catch the most recent run
start_time = (now - datetime.timedelta(hours=12)).strftime('%Y-%m-%dT%H:00:00Z')
end_time = now.strftime('%Y-%m-%dT%H:00:00Z')

wfs_url = (
    "https://opendata.fmi.fi/wfs?service=WFS&version=2.0.0&request=GetFeature&"
    f"storedquery_id=fmi::forecast::harmonie::surface::grid&"
    f"starttime={start_time}&endtime={end_time}"
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
    raise Exception("No model runs found in the last 12 hours")

latest_origintime = max(origintimes)
print(f"Latest model run (origintime): {latest_origintime}")

# --- Step 2: Download the NetCDF for the latest run ---
download_url = (
    "https://opendata.fmi.fi/download?"
    "producer=harmonie_scandinavia_surface&"
    "param=Temperature&"
    "format=netcdf&"
    "bbox=19,56,30,61&"
    f"origintime={latest_origintime}&"
    f"starttime={latest_origintime}&"
    "timestep=60&"
    "timesteps=67&"  # Full forecast run
    "projection=EPSG:4326&"
    "levels=2m"
)

response = requests.get(download_url, timeout=120)
response.raise_for_status()

nc_path = "harmonie.nc"
with open(nc_path, "wb") as f:
    f.write(response.content)

# --- Step 3: Load and convert temperature to °C ---
ds = xr.open_dataset(nc_path)
# The temperature variable is usually called 'Temperature' – confirm with print(ds) if needed
temp_k = ds['Temperature']  # change to ds['t2m'] if that's the name
temp_c = temp_k - 273.15

# --- Step 4: Load your custom color ramp from the .qml file ---
qml_path = "temperature_color_table.qml"
if not os.path.exists(qml_path):
    raise FileNotFoundError("temperature_style.qml not found – add your QML file to the repo")

tree = ET.parse(qml_path)
root = tree.getroot()

# Find colorramp items (works for most singleband pseudocolor or paletted styles)
items = []
for item in root.findall(".//colorrampshader/colorrampitem"):
    value = float(item.get('value'))
    color_str = item.get('color')  # format "#rrggbb,aa"
    if color_str.startswith('#'):
        # Convert #RRGGBBAA or #RRGGBB to rgba 0-1
        color_str = color_str.lstrip('#')
        lv = len(color_str)
        rgb = tuple(int(color_str[i:i + lv // 3], 16) / 255 for i in range(0, lv, lv // 3))
        if len(rgb) == 3:
            rgb = rgb + (1.0,)  # add full opacity if missing
    else:
        # fallback: sometimes it's "r,g,b,a"
        rgb = tuple(int(x)/255 for x in color_str.split(','))
    items.append((value, rgb))

# Fallback for gradient stops format
if not items:
    stops_elem = root.find(".//prop[@k='stops']")
    if stops_elem is not None:
        stops_text = stops_elem.get('v')
        for part in stops_text.split(':'):
            frac, color_part = part.split(';')
            r,g,b,a = map(int, color_part.split(','))
            # Approximate value from fraction (0-1) – you'll need to set min/max manually
            # This is rough; better to use discrete items above
            pass

if not items:
    raise Exception("Could not extract color ramp from .qml – check file or use exported txt instead")

items.sort(key=lambda x: x[0])
values = [i[0] for i in items]
colors = [i[1] for i in items]

cmap = ListedColormap(colors)
norm = Normalize(vmin=min(values), vmax=max(values))  # Use your QGIS min/max

# --- Step 5: Plot the analysis (first timestep) ---
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
plt.title(f"HARMONIE 2m Temperature (°C)\nRun: {latest_origintime} UTC (analysis)")

output_path = "map.png"
plt.savefig(output_path, dpi=200, bbox_inches='tight')
plt.close()

print(f"Map successfully saved to {output_path}")
