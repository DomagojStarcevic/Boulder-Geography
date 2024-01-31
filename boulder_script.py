import rasterio
from rasterstats import zonal_stats
import geopandas as gpd
from shapely.geometry import LineString
import matplotlib.pyplot as plt
from rasterio.plot import show
import openpyxl

#user defined block name
block_attribute="02"
#path to GeoTiff file
file_path_tif='MBES grid/Test_Encoded_Depths_File.tif'
#path to SHP file
file_path_shp="Boulder polygons/Test_Manually_Picked_Boulders.shp"
#path to output file
file_path_output_shp="Boulder polygons/Test_Manually_Picked_Boulders_Centroids.shp"

mbes_raster = rasterio.open(file_path_tif)
mbes_array = mbes_raster.read(1)
mbes_affine = mbes_raster.transform

boulders = gpd.read_file(file_path_shp)

#fig, ax = plt.subplots(figsize=(10,10))
#show(mbes_raster,1,ax=ax)
#boulders.plot(ax=ax,facecolor='None', linewidth=2)
#boulders.centroid.plot(ax=ax,color='black', markersize=10)

boulders['Block']=block_attribute
boulders['Poly_ID']=boulders.index
boulders['Target ID']="MBES_"+block_attribute+"_"+boulders['Poly_ID'].apply(str)
boulders['Easting']=boulders.centroid.x
boulders['Northing']=boulders.centroid.y


depths=[]
for centroid in boulders.centroid:
    depths.append(mbes_array[mbes_raster.index(centroid.x,centroid.y)])

boulders['WaterDepth']=depths

rects=boulders.minimum_rotated_rectangle()
#rects.plot(ax=ax, facecolor='None', linewidth=2)

def get_width_length(rect):
    
    rect_geos = gpd.GeoSeries([rect])
    
    rect_geometry = rect_geos.geometry.iloc[0]
    rect_coords = list(rect_geometry.exterior.coords)

    side_lengths = [
        LineString([rect_coords[i], rect_coords[i + 1]]).length
        for i in range(len(rect_coords) - 1)
    ]
    width = min(side_lengths)
    length = max(side_lengths)

    return width,length

rects_widths, rects_lengths = [],[]
for rect in rects:
    width, length = get_width_length(rect)
    rects_widths.append(width)
    rects_lengths.append(length)
    
boulders['Width']=rects_widths
boulders['Length']=rects_lengths

boulders_zonal_stats = zonal_stats(boulders, mbes_array, affine=mbes_affine, stats=['min', 'max'])
heights=[]
for zone_stats in boulders_zonal_stats:
    heights.append(float(zone_stats['max']-zone_stats['min']))
boulders['Height']=heights

output_gdf = gpd.GeoDataFrame(boulders[['geometry','Poly_ID','Target ID', 'Block','Easting','Northing','WaterDepth','Length','Width','Height']], geometry=boulders.centroid)
output_gdf.to_file(file_path_output_shp)
#output_gdf.to_excel("Boulder polygons/Test_Manually_Picked_Boulders_Centroids.xlsx")