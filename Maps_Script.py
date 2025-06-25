### Plots on `OpenStreetMap`
## Credits to Myriam Prasow-Émond, PhD candidate in Earth Science & Engineering, Imperial College London (m.prasow-emond22@imperial.ac.uk)

pip install --upgrade osmnx
pip install --upgrade earthengine-api
pip install --upgrade geemap
pip install --upgrade rasterio

from io import StringIO
import ee # this is the 'earthengine-api' package
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import osmnx as ox # OpenStreetMap API package
import shapely # for geometry operations
import geopandas as gpd
import pandas as pd
import geemap
import requests
import os
import zipfile
import xarray as xr
import rasterio
%matplotlib inline 
# or qt for interactive plots

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'



## Define the area of interest -----------------------------------------------------
# There are various methods to query OSM data for a specific location, such as using coordinates (longitude, latitude), a polygon, or a bounding box.
# Can also input the name of a place, such as Pakistan

roi = 'Pakistan'
# Extent of Pakistan
gdf_Pakistan = ox.geocode_to_gdf(roi)

# Define provinces
provinces = ['Khyber Pakhtunkhwa, Pakistan', 'Punjab, Pakistan', 'Sindh, Pakistan', 'Balochistan, Pakistan', 'Gilgit-Baltistan, Pakistan','Kaschmir, Pakistan']

gdfs = []
for province in provinces:
    gdf = ox.geocode_to_gdf(province)
    gdfs.append(gdf)

gdf_Pakistan = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True))



## Values of e.g. MAF in the order of provinces -----------------------------------------------------
participants = gpd.GeoDataFrame(np.array([0.431,0.432,0.0885,0.0411,0.000719,0.00511]))
gdf_Pakistan['Participants'] = participants



## plot the data -----------------------------------------------------
fig, ax = plt.subplots(figsize=(15, 15))

gdf_Pakistan.plot(ax=ax, column='Participants', legend=True, cmap='RdPu', legend_kwds={
        'shrink': 0.4,         # Shrink the size of the colorbar
        'aspect': 20,          # Adjust the aspect ratio
        'orientation': 'vertical',  # Orientation of the colorbar
    }
)

# Add title and province labels
ax.set_title("Participants by provincial origin", fontsize=16, fontweight='bold')

# Aesthetics
ax.axis('off')

# Save the plot to a file
output_path = "participants_by_province.png"  # Specify your desired file path and name
plt.savefig(output_path, dpi=300, bbox_inches='tight')  # Adjust dpi for resolution and bbox for padding

# Show the plot
plt.show()
