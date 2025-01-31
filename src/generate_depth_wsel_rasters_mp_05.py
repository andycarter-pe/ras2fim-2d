# Script 05 - From the HEC-RAS output, create floodplain rasters of the water
# surface elevation (WSEL) and depth for the simulation.
#
# This requires as HEC-RAS 2D model from which the following inputs are provided:
#  (1) plan file where the run has been executed
#  (2) terrain that will be used to compute the flood inundation raster stack
#  (3) hydraulic results geopackage (from previous steps)
#  (4) directory to write the output rasters
#
# Created by: Andy Carter, PE
# Created - 2024.07.08

# ************************************************************
import h5py
from shapely.geometry import Point, LineString, Polygon, MultiPolygon, box, mapping, shape
import geopandas as gpd
import pandas as pd
import numpy as np
import os

import rasterio
from rasterio.transform import from_origin
from rasterio.transform import Affine
from rasterio.enums import Resampling
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject
from rasterio.features import shapes
import rasterio.mask

from scipy.interpolate import griddata
from scipy.spatial import ConvexHull
import configparser
import multiprocessing as mp

from osgeo import gdal, gdalconst

#import concurrent.futures
import argparse
import time
import datetime
import warnings
# ************************************************************


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist" % arg)
    else:
        # File exists so return the directory
        return arg
        return open(arg, 'r')  # return an open file handle
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# ------------------------
def fn_get_group_names(hdf5_file_path, group_path):
    """
    Retrieve the names of groups within a specified HDF5 file under a given group path.

    Parameters:
    hdf5_file_path (str): The file path to the HDF5 file.
    group_path (str): The path to the group whose subgroups' names are to be retrieved.

    Returns:
    list or None: A list containing the names of groups found under the specified group path. 
                  Returns None if the group path does not exist in the HDF5 file.
    """
    try:
        with h5py.File(hdf5_file_path, 'r') as hdf_file:
            # Check if the specified group path exists
            if group_path in hdf_file:
                group = hdf_file[group_path]

                # Extract names of HDF5 Group objects
                group_names = [name for name in group if isinstance(group[name], h5py.Group)]

                return group_names
            else:
                print(f"Group '{group_path}' not found in the HDF5 file.")
                return None
    except Exception as e:
        print(f"An error occurred: {e}")
# ------------------------


# ~~~~~~~~~~~~~~~~~~~~
def fn_gdf_cell_centerpoints(hdf_file_path,list_unique_indices_sorted, gdf_cells_wsel):

    # Specify the HDF5 file path and group path
    group_path = '/Geometry/2D Flow Areas/'

    # Get names of HDF5 Group objects in the specified group
    list_group_names = fn_get_group_names(hdf_file_path, group_path)

    # for now... lets assume that there is only one 2D area
    str_cell_info_path = '/Geometry/2D Flow Areas/Cell Info'

    str_cell_center_point_path = group_path + list_group_names[0] + '/' + 'Cells Center Coordinate'

    # From the plan HDF, get the cells center coordinates (not the centroid)
    with h5py.File(hdf_file_path, 'r') as hdf_file:
        arr_cell_center_coords = hdf_file[str_cell_center_point_path][:]

    # Filter the array to those cells that are wet
    arr_center_coords_wet_cells = arr_cell_center_coords[list_unique_indices_sorted]

    geometry = [Point(x, y) for x, y in arr_center_coords_wet_cells]

    # Create a GeoDataFrame
    gdf_center_points = gpd.GeoDataFrame(geometry=geometry, columns=['geometry'], crs=gdf_cells_wsel.crs)

    # create a new coloumn that contains the cell index
    gdf_center_points['cell_idx'] = list_unique_indices_sorted
    
    return (gdf_center_points)
# ~~~~~~~~~~~~~~~~~~~~


# -------------------------
def fn_get_raster_resolution(str_dem_path):
    
    # Get the resolution of the ground DEM
    # Open the raster file
    dataset = gdal.Open(str_dem_path)

    # Get raster resolution
    flt_pixel_width = dataset.GetGeoTransform()[1]

    # Close the dataset raster file
    dataset = None  # Release the reference to the dataset
    
    return(flt_pixel_width)
# -------------------------


# .....................
def fn_wet_cells_per_run(str_wsel_col_name, gdf_cells_wsel, gdf_center_points):

    # from geodataframe, create dataframe from only the 'cell_idx' and 'wsel' coloumns
    df_wsel = gdf_cells_wsel[['cell_idx', str_wsel_col_name]]

    # left join the WSEL per cell to populate the geodataframe of cells with WSEL
    gdf_wsel_wet_cells_pnt = pd.merge(gdf_center_points, df_wsel, on='cell_idx', how='left')

    # Drop rows where "wsel_max_1" is NaN
    gdf_wsel_wet_cells_cleaned_pnt = gdf_wsel_wet_cells_pnt.dropna(subset=[str_wsel_col_name]).copy()

    # Drop rows where "wsel_max_1" is NaN
    gdf_wsel_wet_cells_cleaned_ar = gdf_cells_wsel.dropna(subset=[str_wsel_col_name]).copy()

    # set constant arbitrary value used to dissolve geometry
    gdf_wsel_wet_cells_cleaned_ar.loc[:, 'val'] = 1

    # create a geodataframe of the merged cells
    gdf_dissolved = gdf_wsel_wet_cells_cleaned_ar.dissolve(by='val')

    # TODO - 2024.06.12 - Is the first polygon correct?
    # get shapely geometry of the first merged polygon
    shp_wet_cells = gdf_dissolved.geometry.iloc[0]

    # Convert Shapely polygon to GeoJSON-like format
    geom = mapping(shp_wet_cells)
    
    return (geom, gdf_dissolved, gdf_wsel_wet_cells_cleaned_pnt, gdf_wsel_wet_cells_cleaned_ar)
# .....................


# -------------
def fn_compute_convex_hull(points):
    hull = ConvexHull(points)
    return [points[vertex] for vertex in hull.vertices]

def fn_calculate_vertex_points(polygon):
    list_points = []
    coords = polygon.exterior.coords
    num_coords = len(coords)
    for i in range(num_coords):  # Iterate over vertex
        x, y = coords[i]
        list_points.append(Point(x, y))
    return list_points

def fn_is_point_inside_convex_hull(point, hull_points):
    hull_polygon = Polygon(hull_points)
    return hull_polygon.contains(point)

# Define a function to merge points within a certain distance threshold
def fn_merge_points_within_threshold(gdf, threshold, str_wsel_col_name):
    merged_points = []
    grouped = gdf.groupby(gdf.geometry.apply(lambda x: (round(x.x, 2), round(x.y, 2))))
    
    for _, group in grouped:
        if len(group) > 1:
            merged_point = Point(group.geometry.x.mean(), group.geometry.y.mean())
            avg_wsel = group[str_wsel_col_name].mean()  # Calculate the average of 'wsel' column
            merged_points.append({'geometry': merged_point, str_wsel_col_name: avg_wsel})
        else:
            merged_points.append({'geometry': group.iloc[0].geometry, str_wsel_col_name: group.iloc[0][str_wsel_col_name]})
    
    return gpd.GeoDataFrame(merged_points)
# -------------

# --------------
def fn_nudge_raster(str_raster_path, str_output_path, str_crs, flt_res):
    flt_res = float(flt_res)  # Ensure the resolution is a float

    # Open the original raster
    with rasterio.open(str_raster_path) as src:
        # Calculate the transformation for the desired CRS
        
        # TODO -- there is a Linux PROJ issue with "str_crs" -- 2025.01.30
        
        transform, width, height = calculate_default_transform(
            src.crs, str_crs, src.width, src.height, *src.bounds
        )
        
        # Debugging - 2025.01.30
        #print(f"Source CRS: {src.crs}")
        #print(f"Target CRS: {str_crs}")
        
        # Update transform to match the desired resolution
        transform = from_origin(
            round(transform.c / flt_res) * flt_res,  # Nudge X to nearest multiple of flt_res
            round(transform.f / flt_res) * flt_res,  # Nudge Y to nearest multiple of flt_res
            flt_res, flt_res
        )

        # Prepare metadata for the output raster
        profile = src.profile
        profile.update({
            'crs': str_crs,
            'transform': transform,
            'width': width,
            'height': height,
            'dtype': 'float32',  # Ensure output raster is float32 for NaNs
            'compress': 'lzw',    # Compression if needed
            'nodata': np.nan      # Set nodata value to NaN
        })

        # Create the output raster
        with rasterio.open(str_output_path, 'w', **profile) as dst:
            # Loop through each band
            for i in range(1, src.count + 1):
                # Create an empty array for the output band with NaN
                output_band = np.full((height, width), np.nan, dtype='float32')

                # Reproject source data to the output band
                reproject(
                    source=rasterio.band(src, i),
                    destination=output_band,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=str_crs,
                    resampling=Resampling.nearest,
                    src_nodata=src.nodata  # Treat src_nodata as NaN
                )

                # Update the output raster band only with valid data
                dst.write(output_band, i)

# --------------


# ~~~~~~~~~~~~~
def fn_set_nan(str_original_raster_path):

    str_intermediate_path = str_original_raster_path[:-4] + '_temp.tif'

    # Open the original raster file
    with rasterio.open(str_original_raster_path) as src:
        profile = src.profile
        profile['nodata'] = np.nan  # Set NoData value to NaN
        data = src.read(1)  # Read the first band

        # Convert NaN values to the specified NoData value
        data_with_nodata = np.where(np.isnan(data), np.nan, data)

        # Create a new raster file with NoData set to NaN
        with rasterio.open(str_intermediate_path, 'w', **profile) as dst:
            dst.write(data_with_nodata, 1)

        #print("New raster with NoData set to NaN has been created.")

        return(str_intermediate_path)
# ~~~~~~~~~~~~~


# --------------------------------
def fn_create_raster_products(gdf_wsel_wet_cells_cleaned_pnt,
                              gdf_wsel_wet_cells_cleaned_ar,
                              str_wsel_col_name,
                              gdf_dissolved,
                              flt_pixel_width,
                              str_dem_output_folder,
                              geom,
                              str_dem_path):

    # List to store created points
    created_points = []
    list_cell_wsel = []

    # Wrap Convex Hull
    convex_hull_points = fn_compute_convex_hull(gdf_wsel_wet_cells_cleaned_pnt.geometry.apply(lambda point: (point.x, point.y)).tolist())

    # Iterate through cell polygons and perform checks
    for _, cell_polygon in gdf_wsel_wet_cells_cleaned_ar.iterrows():
        vertex_points = fn_calculate_vertex_points(cell_polygon.geometry)
        for v_point in vertex_points:
            if not fn_is_point_inside_convex_hull(v_point, convex_hull_points):
                created_points.append(v_point)
                list_cell_wsel.append(cell_polygon[str_wsel_col_name])

    # Create a new GeoDataFrame from created points
    gdf_additional_points = gpd.GeoDataFrame(geometry=created_points)

    gdf_additional_points[str_wsel_col_name] = list_cell_wsel
    gdf_additional_points.crs = gdf_wsel_wet_cells_cleaned_pnt.crs

    # Remove the duplicate points... that sit on two cells
    # Set the distance threshold
    threshold_distance = 0.01  # in feet

    # Merge points within the threshold distance
    gdf_additional_points_duplicates_removed = fn_merge_points_within_threshold(gdf_additional_points, threshold_distance, str_wsel_col_name)

    gdf_additional_points_duplicates_removed.crs = gdf_wsel_wet_cells_cleaned_pnt.crs

    points1 = gdf_wsel_wet_cells_cleaned_pnt[['geometry', str_wsel_col_name]].copy()

    if len(gdf_additional_points_duplicates_removed) > 0:
        points2 = gdf_additional_points_duplicates_removed[['geometry', str_wsel_col_name]].copy()
        points = pd.concat([points1, points2], axis=0)
    else:
        points = points1

    # create the points of computed wsel
    x = points.geometry.x.values
    y = points.geometry.y.values
    z = points[str_wsel_col_name].values

    # Define the extent of the cell point raster
    xmin, ymin, xmax, ymax = gdf_dissolved.total_bounds
    x_res = flt_pixel_width  # resolution in x direction
    y_res = flt_pixel_width  # resolution in y direction

    # Create the grid
    xi = np.arange(xmin, xmax, x_res)
    yi = np.arange(ymin, ymax, y_res)
    xi, yi = np.meshgrid(xi, yi)

    # Perform bilinear interpolation
    zi = griddata((x, y), z, (xi, yi), method='linear')

    # Define the affine transformation
    transform = Affine.translation(xmin, ymin) * Affine.scale(x_res, y_res)

    # ~~~~~~~~~~~~~~~~~~~~~~~~
    # Begin the creation of the raster products
    list_raster_to_purge = []

    # ~~~~~~~
    # 01 - the wsel raster that is not masked

    # Define output raster file path
    output_raster_path_01 = os.path.join(str_dem_output_folder,'01_' + str_wsel_col_name + '_interp.tif')

    # Write the interpolated raster to a GeoTIFF file
    with rasterio.open(output_raster_path_01, 'w', driver='GTiff', height=zi.shape[0], width=zi.shape[1], count=1,
                       dtype=zi.dtype,
                       crs=gdf_wsel_wet_cells_cleaned_ar.crs,
                       transform=transform) as dst:
        dst.write(zi, 1)

    list_raster_to_purge.append(output_raster_path_01)
    
    # ~~~~~~~
    # 02 - the masked wsel

    # Open the raster dataset in read mode
    with rasterio.open(output_raster_path_01, 'r') as dst:

        # Mask the interpolated raster with the polygon
        zi_masked, transform_masked = rasterio.mask.mask(dst, [geom], crop=True, nodata=np.nan)  # Specify nodata value as np.nan

    # Define output masked raster file path
    output_masked_raster_path_02 = os.path.join(str_dem_output_folder,'02_' + str_wsel_col_name + '_interp_masked.tif')

    # Write the masked interpolated raster to a GeoTIFF file
    with rasterio.open(output_masked_raster_path_02, 'w', driver='GTiff', height=zi_masked.shape[1], width=zi_masked.shape[2], count=1,
                       dtype=zi_masked.dtype,
                       crs=gdf_wsel_wet_cells_cleaned_ar.crs,
                       transform=transform_masked,
                       nodata=np.nan) as dst_masked:  # Specify nodata value as np.nan
        dst_masked.write(zi_masked[0], 1)

    list_raster_to_purge.append(output_masked_raster_path_02)
    
    # ~~~~~~~
    # 03 and 04 - aligning the wsel and ground rasters

    try:
        # Open and read the first raster
        with rasterio.open(output_masked_raster_path_02) as src1:
            raster1 = src1.read(1)
            transform1 = src1.transform
            crs1 = src1.crs

            try:
                # Open and read the second raster (ground dem raster)
                with rasterio.open(str_dem_path) as src2:
                    raster2 = src2.read(1)
                    transform2 = src2.transform
                    crs2 = src2.crs

                    # TODO - I want all rasters to be alligned to the limits of the terrain.  This should
                    # preserve the same pixels on all the 'runs'. - 2024.07.08
                    
                    bbox = gpd.GeoDataFrame(geometry=[box(*src2.bounds), box(*src2.bounds)], crs=crs1)
                    common_extent = bbox.unary_union.bounds

                    # Resample and align rasters
                    dst_crs = crs1
                    transform, width, height = calculate_default_transform(crs2, dst_crs, src2.width, src2.height, *src2.bounds)
                    kwargs = src2.meta.copy()
                    kwargs.update({
                        'crs': dst_crs,
                        'transform': transform,
                        'width': width,
                        'height': height
                    })

                    raster2_aligned = np.zeros_like(raster1)
                    reproject(
                        source=raster2,
                        destination=raster2_aligned,
                        src_transform=src2.transform,
                        src_crs=src2.crs,
                        dst_transform=transform,
                        dst_crs=dst_crs,
                        resampling=Resampling.nearest)

                    # Clip rasters
                    output_raster_path_03 = os.path.join(str_dem_output_folder,'03_' + str_wsel_col_name + '_interp_mask_align.tif')
                    output_raster_path_04 = os.path.join(str_dem_output_folder,'04_' + str_wsel_col_name + '_ground_align.tif')

                    with rasterio.open(output_raster_path_03, 'w', **src1.meta) as dst1:
                        out_img1, out_transform1 = mask(src1, [bbox.iloc[0]['geometry']], crop=True)
                        dst1.write(out_img1)

                    with rasterio.open(output_raster_path_04, 'w', **src2.meta) as dst2:
                        out_img2, out_transform2 = mask(src2, [bbox.iloc[1]['geometry']], crop=True)
                        dst2.write(out_img2)

                    list_raster_to_purge.append(output_raster_path_03)
                    list_raster_to_purge.append(output_raster_path_04)

            except rasterio.errors.RasterioIOError as e:
                print("Error opening or processing raster file:", e)

    except rasterio.errors.RasterioIOError as e:
        print("Error opening raster file:", e)
        
    return list_raster_to_purge
# --------------------------------


# -------------------
def fn_raster_subtraction(raster1_path, raster2_path, output_folder, str_wsel_col_name):
    
    list_created_rasters = []
    
    # Open the first raster
    raster1 = gdal.Open(raster1_path, gdalconst.GA_ReadOnly)
    if raster1 is None:
        raise FileNotFoundError(f"Raster not found at {raster1_path}")

    # Open the second raster
    raster2 = gdal.Open(raster2_path, gdalconst.GA_ReadOnly)
    if raster2 is None:
        raise FileNotFoundError(f"Raster not found at {raster2_path}")

    # Get the extent and resolution of the first raster
    x_min, pixel_width, _, y_max, _, pixel_height = raster1.GetGeoTransform()
    x_max = x_min + raster1.RasterXSize * pixel_width
    y_min = y_max + raster1.RasterYSize * pixel_height

    # Resample the second raster to match the extent and resolution of the first raster
    #resampled_raster2_path = r'E:\working\05_ground_resampled_align.tif'
    resampled_raster2_path = os.path.join(output_folder,'05_' + str_wsel_col_name + '_ground_resampled_align.tif')
    list_created_rasters.append(resampled_raster2_path)
    
    gdal.Warp(resampled_raster2_path, raster2, format='GTiff', outputBounds=[x_min, y_min, x_max, y_max], xRes=pixel_width, yRes=pixel_height)

    # Perform raster subtraction
    raster1_data = raster1.ReadAsArray()
    resampled_raster2 = gdal.Open(resampled_raster2_path, gdalconst.GA_ReadOnly)
    resampled_raster2_data = resampled_raster2.ReadAsArray()

    result_data = raster1_data - resampled_raster2_data

    # Create output raster
    driver = gdal.GetDriverByName('GTiff')
    
    output_raster_path_06 = os.path.join(output_folder,'06_' + str_wsel_col_name + '_depth.tif')
    output_raster = driver.Create(output_raster_path_06, raster1.RasterXSize, raster1.RasterYSize, 1, gdalconst.GDT_Float32)

    # Set projection and transform
    output_raster.SetProjection(raster1.GetProjection())
    output_raster.SetGeoTransform(raster1.GetGeoTransform())

    # TODO - 2024.02.23 -- likely writing two raster bands
    # Write result data to output raster
    output_band = output_raster.GetRasterBand(1)
    output_band.WriteArray(result_data)

    # Close rasters
    output_band = None
    output_raster = None
    raster1 = None
    resampled_raster2 = None
    
    list_created_rasters.append(output_raster_path_06)
    
    return(list_created_rasters)
# -------------------


# --------------------
def fn_fix_depth_raster(output_raster_path_06, str_dem_output_folder, str_wsel_col_name):

    list_depth_raster_path = []
    
    # Build the raster file paths
    output_raster_path_07 = os.path.join(str_dem_output_folder,'07_' + str_wsel_col_name + '_depth_with_nan.tif')

    with rasterio.open(output_raster_path_06, 'r') as src:
        # Read raster data
        data = src.read(1)

        # Set values less than zero to NaN
        data[data < 0] = float('nan')

        # Copy metadata
        kwargs = src.meta.copy()

        # Update metadata for the new file
        kwargs.update(dtype=rasterio.float32, compress='lzw')

        # Write the modified data to a new GeoTIFF file
        with rasterio.open(output_raster_path_07, 'w', **kwargs) as dst:
            dst.write(data, 1)
            
        list_depth_raster_path.append(output_raster_path_07)
        
        return(list_depth_raster_path)
# --------------------


# -------------------
def fn_raster_addition_wsel(ground_raster_path, depth_raster_path, output_folder,str_wsel_col_name):

    list_wsel_raster_path = []
    
    # Open the GeoTIFF rasters
    ground_raster = gdal.Open(ground_raster_path)
    depth_raster = gdal.Open(depth_raster_path)

    # Read raster data as NumPy arrays
    ground_data = ground_raster.ReadAsArray()
    depth_data = depth_raster.ReadAsArray()

    # Convert depth data to float, as NaN is a floating-point concept
    depth_data = depth_data.astype(float)

    # Mask NaN values in depth data
    depth_data[np.isnan(depth_data)] = np.nan

    # Add the two rasters where depth is not NaN
    result = np.where(~np.isnan(depth_data), ground_data + depth_data, np.nan)

    # Save the result to a new GeoTIFF
    
    output_path = os.path.join(output_folder, '08_' + str_wsel_col_name + '_wsel_with_nan.tif')
    
    driver = gdal.GetDriverByName("GTiff")
    #output_raster = driver.Create(output_path, ground_raster.RasterXSize, ground_raster.RasterYSize, 1, gdal.GDT_Float32)
    output_raster = driver.Create(output_path, ground_raster.RasterXSize, ground_raster.RasterYSize, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])
    output_raster.SetProjection(ground_raster.GetProjection())
    output_raster.SetGeoTransform(ground_raster.GetGeoTransform())
    output_band = output_raster.GetRasterBand(1)
    output_band.WriteArray(result)
    output_band.FlushCache()

    # Close the rasters
    ground_raster = None
    depth_raster = None
    output_raster = None
    
    list_wsel_raster_path.append(output_path)
    
    return(list_wsel_raster_path)
# -------------------


# ----------------------
def fn_purge_rasters(list_raster_to_purge):
    for file_path in list_raster_to_purge:
        try:
            os.remove(file_path)
            #print(f"File '{file_path}' deleted successfully.")
        except OSError as e:
            print(f"Error deleting file '{file_path}': {e}")
# ----------------------


# ------------------------------------
def fn_process_raster_stack(str_wsel_col_name,
                            gdf_cells_wsel,
                            gdf_center_points,
                            flt_pixel_width,
                            str_output_dir):
    
    geom, gdf_dissolved, gdf_wsel_wet_cells_cleaned_pnt, gdf_wsel_wet_cells_cleaned_ar = fn_wet_cells_per_run(str_wsel_col_name,
                                                                                                              gdf_cells_wsel,
                                                                                                              gdf_center_points)

    #print(f'Creating Rasters for {str_wsel_col_name}...')
    list_raster_to_purge = fn_create_raster_products(gdf_wsel_wet_cells_cleaned_pnt,
                                                     gdf_wsel_wet_cells_cleaned_ar,
                                                     str_wsel_col_name,
                                                     gdf_dissolved,
                                                     flt_pixel_width,
                                                     str_output_dir,
                                                     geom)

    print(f'Computing Depth and WSEL Raster for {str_wsel_col_name}...')
    list_created_rasters = fn_raster_subtraction(list_raster_to_purge[2],
                                                 list_raster_to_purge[3],
                                                 str_output_dir,
                                                 str_wsel_col_name)

    list_depth_raster_path = fn_fix_depth_raster(list_created_rasters[1],
                                                 str_output_dir,
                                                 str_wsel_col_name)

    list_wsel_raster_path = fn_raster_addition_wsel(list_created_rasters[0],
                                                    list_depth_raster_path[0],
                                                    str_output_dir,
                                                    str_wsel_col_name)

    # Purge the created rasters
    # Keeping #7(depth) & #8(wsel)
    fn_purge_rasters(list_raster_to_purge)
    fn_purge_rasters(list_created_rasters)
    
    return str_wsel_col_name
# ------------------------------------


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
def fn_make_floodplain_gdf(str_input_raster_filepath, flt_minvalue, flt_max_value):

    """
    Converts a raster file representing flooded areas to a GeoDataFrame of a single MultiPolygon representing the floodplain.

    Parameters:
    - str_input_raster_filepath (str): File path of the input raster file.
    - flt_minvalue (float): minimum value of depth allowed
    - flt_max_value (float): maximum value of depth allowed

    Returns:
    - gdf_single_floodplain (geopandas.GeoDataFrame): GeoDataFrame containing a merged MultiPolygon representing the floodplain.
    """
    
    # ----------------------
    # open the raster file and convert to 'binary' of flooded=1 and null = 255
    with rasterio.open(str_input_raster_filepath) as src:
        # Read the raster data as a numpy array
        data = src.read(1)

        # Set values greater than flt_minvalue and <= flt_max_value to 1, others to 255
        data = np.where((data > flt_minvalue) & (data <= flt_max_value), 1, 255).astype('uint8')

        # Create a new raster file with the updated data type and nodata value
        profile = src.profile
        profile.update(dtype=rasterio.uint8, nodata=255)

        # Use in-memory storage
        output_raster_memory = rasterio.MemoryFile()
        with output_raster_memory.open(**profile) as dst:
            dst.write(data, 1)

        #print("Conversion completed. Output raster stored in memory.")

    # ----------------------
    # convert the binary raster to geodataframe of polygons
    # Open the input GeoTIFF file
    with output_raster_memory.open() as src:
        # Read the raster data
        data = src.read(1, masked=True)  # using masked=True to handle NoData as a mask

        # Set all non-null values to 1
        data = np.where(data.mask, 255, 1).astype('uint8')

        # Extract shapes from the raster data
        shapes_gen = shapes(data, mask=data != 255, transform=src.transform)

        # Convert shapes to Shapely geometries and create a GeoDataFrame
        geometries = [shape(geom) for geom, _ in shapes_gen]
        gdf = gpd.GeoDataFrame(geometry=geometries, crs=src.crs)

    # Merge all polygons into a single MultiPolygon
    merged_polygon = gdf['geometry'].unary_union

    # Create a new GeoDataFrame with the merged MultiPolygon
    gdf_single_floodplain = gpd.GeoDataFrame(geometry=[merged_polygon], crs=gdf.crs)

    # Close output_raster_memory to release the memory
    output_raster_memory.close()

    return(gdf_single_floodplain)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~


# -------------------------------
def fn_process_wsel(str_wsel_col_name, gdf_cells_wsel,
                    gdf_center_points, flt_pixel_width,
                    str_output_crs, flt_desired_res,
                    str_output_dir, str_dem_path):
    try:
        # Initialization
        #print('--------')
        #print(f'Processing {str_wsel_col_name}')
        
        # Wet points per run (wsel) and polygon boundary
        geom, gdf_dissolved, gdf_wsel_wet_cells_cleaned_pnt, gdf_wsel_wet_cells_cleaned_ar = fn_wet_cells_per_run(
            str_wsel_col_name, gdf_cells_wsel, gdf_center_points)
        
       # print('Creating Rasters...')
        list_raster_to_purge = fn_create_raster_products(
            gdf_wsel_wet_cells_cleaned_pnt, gdf_wsel_wet_cells_cleaned_ar, 
            str_wsel_col_name, gdf_dissolved, flt_pixel_width, str_output_dir, geom, str_dem_path)
        
        #print('Computing WSEL Raster...')
        list_created_rasters = fn_raster_subtraction(
            list_raster_to_purge[2], list_raster_to_purge[3], str_output_dir, str_wsel_col_name)
        
        list_depth_raster_path = fn_fix_depth_raster(
            list_created_rasters[1], str_output_dir, str_wsel_col_name)
        
        list_wsel_raster_path = fn_raster_addition_wsel(
            list_created_rasters[0], list_depth_raster_path[0], str_output_dir, str_wsel_col_name)
        
        # Nudge depth raster
        str_itermediate_raster_filepath = fn_set_nan(list_depth_raster_path[0])
        str_file_dir, str_file_name = os.path.split(str_itermediate_raster_filepath)
        str_output_filename = "09" + str_file_name[2:-18] + '_nudge.tif'
        str_output_path_depth = os.path.join(str_file_dir, str_output_filename)
        fn_nudge_raster(str_itermediate_raster_filepath, str_output_path_depth, str_output_crs, flt_desired_res)
        
        #print('Begin nudge wsel: 10...')
        # Nudge wsel raster
        str_itermediate_raster_filepath = fn_set_nan(list_wsel_raster_path[0])
        str_file_dir, str_file_name = os.path.split(str_itermediate_raster_filepath)
        str_output_filename = "10" + str_file_name[2:-18] + '_nudge.tif'
        str_output_path = os.path.join(str_file_dir, str_output_filename)
        fn_nudge_raster(str_itermediate_raster_filepath, str_output_path, str_output_crs, flt_desired_res)
        
        # Purge the created rasters
        list_created_rasters.append(list_depth_raster_path[0])
        list_created_rasters.append(list_wsel_raster_path[0])
        fn_purge_rasters(list_raster_to_purge)
        fn_purge_rasters(list_created_rasters)
        
        # Create a gdf of the floodplain from the depth raster
        flt_minvalue = 0
        flt_max_value = 10000
        gdf_floodplain = fn_make_floodplain_gdf(str_output_path_depth, flt_minvalue, flt_max_value)
        
        # Save the gdf_floodplain as a GeoPackage
        str_flood_limits_gpkg = os.path.join(str_output_dir, '11_' + str_wsel_col_name + '_floodplain_ar.gpkg')
        gdf_floodplain.to_file(str_flood_limits_gpkg)
        
        #print(f'Completed processing for {str_wsel_col_name}')
    except Exception as e:
        print(f'Error processing {str_wsel_col_name}: {e}')
# -------------------------------


# ----------------
def fn_str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {'true', 't', '1'}:
        return True
    elif value.lower() in {'false', 'f', '0'}:
        return False
    else:
        raise argparse.ArgumentTypeError(f"Boolean value expected. Got '{value}'.")
# ----------------



# ---------------------------------------------
def fn_generate_depth_wsel_rasters(str_config_file_path,
                                   str_hdf_file_path,
                                   str_dem_path,
                                   str_hydraulic_results_gpkg,
                                   str_output_dir,
                                   b_print_output):
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    if b_print_output:
        print(" ")
        print("+=================================================================+")
        print("|       GENERATE DEPTH AND WATER SURFACE ELEVATION RASTERS        |")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        
        print("  ---(c) INPUT GLOBAL CONFIGUTATION FILE: " + str_config_file_path)
        print("  ---(i) INPUT HEC-RAS PLAN (HDF): " + str_hdf_file_path)
        print("  ---(t) INPUT TERRAIN (TIF): " + str_dem_path)
        print("  ---(g) INPUT HYDRAUNIC RESULTS (GPKG): " + str_hydraulic_results_gpkg)
        print("  ---(o) OUTPUT PATH: " + str_output_dir)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    
    # Read the GeoPackage file
    gdf_cells_wsel = gpd.read_file(str_hydraulic_results_gpkg, layer='02_cells_wsel_ar')
    
    # Get unique values of 'cell_idx' and sort them
    list_unique_indices_sorted = sorted(gdf_cells_wsel['cell_idx'].unique())
    
    list_wsel_col_names = [col for col in gdf_cells_wsel.columns if col.startswith('wsel')]
    
    # Create a gdf of all cells that are "wet"
    gdf_center_points = fn_gdf_cell_centerpoints(str_hdf_file_path,
                                                 list_unique_indices_sorted,
                                                 gdf_cells_wsel)
    
    # get the resolution of the ground DEM
    flt_pixel_width = fn_get_raster_resolution(str_dem_path)
    
    # --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    
    # ----
    # Get the variables in the [04_time_gradient] section
    if '04_time_gradient' in config:
        section = config['06_output_products']
        
        # Read variables in the section
        str_output_crs = section.get('str_output_crs', '')
        flt_desired_res = section.getfloat('flt_desired_res', 0)
    else:
        print("[06_output_products] section not found in the config file.")
    # ----
    
    # **********
    # HARED CODED MAX NUMBER OF PROCESSORS
    int_max_processors = 16
    # **********
    num_processors = (mp.cpu_count() - 1)
    
    if num_processors > int_max_processors:
        num_processors = int_max_processors
    
    # TODO -- 2025.01.09 - add tqdm progress bar to this multiprocessing.
    # Create a pool of workers
    with mp.Pool(processes=num_processors) as pool:
        # Use partial to pass additional arguments to process_wsel
        from functools import partial
        process_func = partial(
            fn_process_wsel, 
            gdf_cells_wsel=gdf_cells_wsel,
            gdf_center_points=gdf_center_points,
            flt_pixel_width=flt_pixel_width,
            str_output_crs=str_output_crs,
            flt_desired_res=flt_desired_res,
            str_output_dir=str_output_dir,
            str_dem_path=str_dem_path
        )
        pool.map(process_func, list_wsel_col_names)
# ---------------------------------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='======= GENERATE DEPTH AND WATER SURFACE ELEVATION RASTERS =======')
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=True,
                        default=r'C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-i',
                        dest = "str_hdf_file_path",
                        help=r'REQUIRED: path to plan HDF file (HDF) Example: E:\ras2fim_test_20250107\03_run_hdf\1884650_wb-2410416_wb-2410422_21-hr_7518-cfs_to_16893-cfs_step_500-cfs.p01.hdf',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-t',
                        dest = "str_dem_path",
                        help=r'REQUIRED: path to terrain dem (TIF) Example: E:\ras2fim_test_20250107\processing\02_model_copies\source_terrain',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-g',
                        dest = "str_hydraulic_results_gpkg",
                        help=r'REQUIRED: hydraulic results geopackage (geopackage) Example: E:\ras2fim_test_20250107\processing\05_large_terrain',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-o',
                        dest = "str_output_dir",
                        help=r'REQUIRED: directory to write model raster output Example: E:\ras2fim_test_20250107\processing\05_large_terrain',
                        required=True,
                        metavar='DIR',
                        type=str)
    
    parser.add_argument('-r',
                        dest = "b_print_output",
                        help=r'OPTIONAL: Print output messages Default: True',
                        required=False,
                        default=True,
                        metavar='T/F',
                        type=fn_str_to_bool)
    

    args = vars(parser.parse_args())
    
    str_config_file_path = args['str_config_file_path']
    str_hdf_file_path = args['str_hdf_file_path']
    str_dem_path = args['str_dem_path']
    str_hydraulic_results_gpkg = args['str_hydraulic_results_gpkg']
    str_output_dir = args['str_output_dir']
    b_print_output = args['b_print_output']

    fn_generate_depth_wsel_rasters(str_config_file_path,
                                   str_hdf_file_path,
                                   str_dem_path,
                                   str_hydraulic_results_gpkg,
                                   str_output_dir,
                                   b_print_output)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~