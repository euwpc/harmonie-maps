import requests
import xml.etree.ElementTree as ET
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap, Normalize
import matplotlib
import datetime
import os
import glob
from PIL import Image
import pandas as pd

matplotlib.use('Agg')

# --- Helper to parse QML color ramp ---
def parse_qml_colormap(qml_file, vmin, vmax):
    tree = ET.parse(qml_file)
    root = tree.getroot()
    items = []
    for item in root.findall(".//colorrampshader/item"):
        value = float(item.get('value'))
        color_hex = item.get('color').lstrip('#')
        r = int(color_hex[0:2], 16) / 255.0
        g = int(color_hex[2:4], 16) / 255.0
        b = int(color_hex[4:6], 16) / 255.0
        items.append((value, (r, g, b, 1.0)))
    items.sort(key=lambda x: x[0])
    colors = [i[1] for i in items]
    return ListedColormap(colors), Normalize(vmin=vmin, vmax=vmax)

# --- Step 1-3: Model run, download, load data (unchanged) ---
# ... (your existing code for wfs_url, download_url, ds load, variables)

# --- Step 4: Load custom colormaps ---
temp_cmap, temp_norm = parse_qml_colormap("temperature_color_table_high.qml", vmin=-40, vmax=50)
cape_cmap, cape_norm = parse_qml_colormap("cape_color_table.qml", vmin=0, vmax=5000)
pressure_cmap, pressure_norm = parse_qml_colormap("pressure_color_table.qml", vmin=890, vmax=1064)
windgust_cmap, windgust_norm = parse_qml_colormap("wind_gust_color_table.qml", vmin=0, vmax=50)
dewpoint_cmap = temp_cmap
dewpoint_norm = Normalize(vmin=-40, vmax=30)
precip_cmap = plt.cm.Blues
precip_norm = Normalize(vmin=0, vmax=10)

# --- Step 5-6: Helper and views (unchanged) ---
# ... (get_analysis, views, variables dict)

# --- Generate for each view ---
for view_key, view_conf in views.items():
    extent = view_conf['extent']
    suffix = view_conf['suffix']

    for var_key, conf in variables.items():
        # Analysis map
        data = get_analysis(conf['var'])
        
        # Accurate min/max only for visible area (analysis only)
        lon_min, lon_max, lat_min, lat_max = extent
        try:
            data_cropped = data.sel(
                lon=slice(lon_min, lon_max),
                lat=slice(lat_max, lat_min),
                method='nearest'
            )
            if data_cropped.size == 0:
                raise ValueError("Empty crop")
            min_val = float(data_cropped.min())
            max_val = float(data_cropped.max())
        except:
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
        
        # Watermark: © tormiinfo.ee — no box, clean
        fig.text(0.99, 0.01, '© tormiinfo.ee', fontsize=10, color='white',
                 ha='right', va='bottom', alpha=0.9,
                 path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground='black')])
        
        plt.title(f"HARMONIE {conf['title']}\nModel run: {run_time_str} | Analysis\nMin: {min_val:.1f} {conf['unit']} | Max: {max_val:.1f} {conf['unit']}")
        plt.savefig(f"{var_key}{suffix}.png", dpi=200, bbox_inches='tight')
        plt.close()

        # Animation
        frames = []
        time_dim = 'time' if 'time' in conf['var'].dims else 'time_h'
        time_values = ds[time_dim].values
        
        for i in range(len(time_values)):
            # After +48h, use 3-hour steps
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

            # Watermark: © tormiinfo.ee — no box
            fig.text(0.99, 0.01, '© tormiinfo.ee', fontsize=10, color='white',
                     ha='right', va='bottom', alpha=0.9,
                     path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground='black')])
            
            valid_dt = pd.to_datetime(time_values[i])
            valid_dt_eet = valid_dt + pd.Timedelta(hours=2)
            valid_str = valid_dt_eet.strftime("%a %d %b %H:%M EET")
            
            plt.title(f"HARMONIE {conf['title']}\nValid: {valid_str} | +{hour_offset}h from run {run_time_str}")

            frame_path = f"frame_{var_key}{suffix}_{i:03d}.png"
            plt.savefig(frame_path, dpi=110, bbox_inches='tight')
            plt.close()
            frames.append(Image.open(frame_path))

        frames[0].save(f"{var_key}{suffix}_animation.gif", save_all=True, append_images=frames[1:], duration=500, loop=0)

        for f in glob.glob(f"frame_{var_key}{suffix}_*.png"):
            os.remove(f)

# --- Cleanup ---
if os.path.exists("harmonie.nc"):
    os.remove("harmonie.nc")

print("All maps + animations generated with custom colormaps")
