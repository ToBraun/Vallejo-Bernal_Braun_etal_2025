# Copyright (C) 2025 by
# Sara M. Vallejo-Bernal <vallejo.bernal@pik-potsdam.de>
# Tobias Braun <tobraun@pik-potsdam.de> 
# All rights reserved.
# GNU General Public License v3.0.

'''
Utility functions for the atmospheric river detection tool underlying the 
PIKART version 1 catalog.
'''

# %% IMPORT MODULES

import numpy as np
import math
import pandas as pd
from netCDF4 import Dataset, stringtochar
from geopy.distance import geodesic
from skimage.util import view_as_windows
from skimage import morphology, measure
from scipy import ndimage
import networkx as nx
import matplotlib.pyplot as plt

# %% CREATE AND SAVE NC FILE

def create_netcdf(time, lat, lon, data, file_name, var, lon_var, lat_var, long_name):

    # Create nc file
    ds = Dataset(file_name, 'w', format='NETCDF4')
    ds.conventions = 'CF-1.12'

    # Create dimensions
    ds.createDimension('time', len(time))
    ds.createDimension(lat_var, len(lat))
    ds.createDimension(lon_var, len(lon))

    # Create variables

    # Longitude
    lon_nc = ds.createVariable(lon_var, np.float32, (lon_var,))
    lon_nc.units = 'degrees_east'
    lon_nc.long_name = 'Longitude'

    # Latitude
    lat_nc = ds.createVariable(lat_var, np.float32, (lat_var,))
    lat_nc.units = 'degrees_north'
    lat_nc.long_name = 'Latitude'

    # Time
    time_nc = ds.createVariable('time', np.int32, ('time',))
    time_nc.units = 'hours since 1900-01-01 00:00:00.0'
    time_nc.long_name = 'Time'
    time_nc.calendar = 'proleptic_gregorian'

    # IVT-related variable
    var_nc = ds.createVariable(var, np.float32, ('time', lat_var, lon_var,),
                               fill_value=-32767., compression='zlib',
                               complevel=6, least_significant_digit=4)
    var_nc.units = 'kg m**-1 s**-1'
    var_nc.long_name = long_name

    # Write data
    lon_nc[:] = lon
    lat_nc[:] = lat
    time_nc[:] = ((time - np.datetime64('1900-01-01T00:00:00')) /
                  np.timedelta64(1, 'h')).astype(int)
    var_nc[:,:,:] = data

    # Close nc file
    ds.close()
    
# %% CHECK THE SPATIAL GRID OF A NC FILE

def check_grid(nc_dataset, var, lon_var, lat_var, spatial_res, lon_lims, lat_lims):

    # Initialize variable
    grid_correct = True

    # Load spatial grid
    lon = nc_dataset[lon_var].values
    lat = nc_dataset[lat_var].values

    # Create reference grid
    if lon_lims[0] < lon_lims[-1]:
        lon_ref = np.arange(
            lon_lims[0], lon_lims[-1] + spatial_res[1], spatial_res[1])
    else:
        lon_ref = np.flip(np.arange(
            lon_lims[-1], lon_lims[0] + spatial_res[1], spatial_res[1]))
        
    if lat_lims[0] < lat_lims[-1]:
        lat_ref = np.arange(
            lat_lims[0], lat_lims[-1] + spatial_res[0], spatial_res[0])
    else:
        lat_ref = np.flip(np.arange(
            lat_lims[-1], lat_lims[0] + spatial_res[0], spatial_res[0]))

    # Check spatial grid
    if not np.array_equal(lon.data, lon_ref):
        grid_correct = False

    if not np.array_equal(lat.data, lat_ref):
        grid_correct = False

    # Print report
    if grid_correct:
        print('%s grid is ok \n' %var.upper())
    else:
        print('GridError: %s must be on a regularly spaced grid with longitude \
              from %s° to %s° and latitude from %s° to %s° \n' %(var.upper(),
              lon_lims[0], lon_lims[-1], lat_lims[0], lat_lims[-1]))

    return grid_correct

# %% GEODESIC DISTANCE IN KM

def geo_dist(lons, lats):

    geodist = np.zeros(np.size(lats,0)-1)

    for i in range(np.size(geodist,0)):

        point0 = (lats[i], lons[i])
        point1 = (lats[i+1], lons[i+1])

        geodist[i] = geodesic(point0, point1).km

    return geodist

# %% AREA OF PIXELS

#-------------------------------------------------------------------------------
# Name: area_grid (similar to areaquad from Matlab)
# Purpose: Surface area of latitude-longitude spheroid rectangle
# Produces a 2D array in geographic coordinate system
# with cell size determined by the user. Value stored in the cell
# is the area of the cell. 
#
# Author: Eugenio Arima, UT Austin
#
# For more details, please check reference:
# Snyder, J. P. (1987). Map projections--A working manual (Vol. 1395):
# US Government Printing Office.
#
# Created: 06/12/2016
# Copyright: (c) ea9267 2016
# Modifications: The spatial grid matches the one of the PIKART catalog
#-------------------------------------------------------------------------------

def area_grid(datum, spatial_res, latitudes):
    '''Creates a 2D array with the area of each grid cell in km2
    Usage: >>> test = area_grid('wgs84, 0.5)'''

    # Ellipsoid major and minor axis.
    datum_lib = {'wgs84': [6378137.0, 6356752.3]}
    a,b = datum_lib.get(datum)

    # Eccentricity
    e2 = (a**2 - b**2)/(a**2)
    e = math.sqrt(e2)

    # Authalic radius
    r2 = math.sqrt((a**2/2)+(b**2/2)*(math.atanh(e)/e))  # in m

    # Output array dimension
    ncol = int(360/spatial_res[1])

    # Borders of grid cells
    lats1 = latitudes - spatial_res[0]/2
    lats2 = latitudes + spatial_res[0]/2

    # Convert to radians
    latin1 = np.radians(lats1)
    latin2 = np.radians(lats2)

    # Convert latitudes to authalic latitudes. See Snyder (1987, p.16 eq. 3-18)
    # (want to find beta, not theta, that's why you subtract the series)
    # factor expansion series
    fact1 = e**2 /3 + 31*e**4 /180 + 517*e**6 /5040 # expansion series 1
    fact2 = 23*e**4 /360 + 251*e**6 /3780 # expansion series 2
    fact3 = 761*e**6 /45360 # expansion series 3
    latout1 = latin1 - fact1*np.sin(2*latin1) + fact2*np.sin(4*latin1) + \
        fact3*np.sin(6*latin1)
    latout2 = latin2 - fact1*np.sin(2*latin2) + fact2*np.sin(4*latin2) + \
        fact3*np.sin(6*latin2)
        
    # Convert zonal resolution to radians
    lon_res = (np.pi/180)*spatial_res[1]

    # Calculate area of square on spherical surface
    cst = r2**2 # just a constant; see Synder 1987.
    area = cst*lon_res*(np.absolute(np.sin(latout1)-np.sin(latout2)))

    # Replicate column over Earth's extent
    grid = np.tile(area, (ncol,1)).T  # in m2
    grid = grid/1000000  # in km2

    return grid
    
# %% EROSION

def erosion(data, kernel):  
        
    # Calculate temporal minima
    min_time = np.nanmin(data, axis=0)
    
    # Extract sliding windows
    windows = view_as_windows(min_time, (2*kernel[1]+1,2*kernel[2]+1))

    # Reshape windows into a 2D array where each row is a flattened window
    reshaped_windows = windows.reshape(-1, np.prod((2*kernel[1]+1,2*kernel[2]+1)))

    # Calculate spatial percentile across each flattened window
    min_space = np.nanmin(reshaped_windows, axis=1)
    
    # Calculate shape of the output    
    valid_output_shape = tuple(np.array((np.size(data,1), np.size(data,2))) - 
                               np.array((2*kernel[1]+1,2*kernel[2]+1)) + 1)

    # Reshape back to match the valid output shape
    result = min_space.reshape(valid_output_shape)
        
    return result

# %% MASK BASED ON THE PERCENTILE OF SLIDING WINDOWS

def mask_perc(data, percentile, kernel):  
        
    # Calculate temporal percentiles
    perc_time = np.nanpercentile(data, percentile, axis=0)
    
    # Extract sliding windows
    windows = view_as_windows(perc_time, (2*kernel[1]+1,2*kernel[2]+1))

    # Reshape windows into a 2D array where each row is a flattened window
    reshaped_windows = windows.reshape(-1, np.prod((2*kernel[1]+1,2*kernel[2]+1)))

    # Calculate spatial percentile across each flattened window
    perc_space = np.nanpercentile(reshaped_windows, percentile, axis=1)
    
    # Calculate shape of the output    
    valid_output_shape = tuple(np.array((np.size(data,1), np.size(data,2))) - 
                               np.array((2*kernel[1]+1,2*kernel[2]+1)) + 1)

    # Reshape back to match the valid output shape
    result = perc_space.reshape(valid_output_shape)
        
    return result

# %% CHECK CYCLICITY OF A MASK

def check_cyclicity(mask):

    left=mask[:,0]
    right=mask[:,-1]
    
    if np.sum(left)>0 and np.sum(right)>0:
        return True
    else:
        return False
    
# %% HOLE CLOSSING

def mask_clossing(mask, ele, pad):
    
    # Pad the mask
    pad_mask = np.pad(mask, (pad,pad), mode='constant', constant_values=0)
    
    # Perform morphological closing
    pad_mask = morphology.closing(pad_mask, footprint=ele)
    
    # Crop the mask
    pad_mask = pad_mask[pad:-pad, pad:-pad]
    
    return pad_mask

# %% RELABEL MASK OF AR SHAPES

def relabel_mask(mask):
    
    # re-assign labels so that there is not gap in label numbers
    unique_labels = list(np.unique(mask))
    unique_labels.remove(0)
    new_mask = np.zeros(np.shape(mask))

    for i, li in enumerate(unique_labels):
        new_mask += np.where(mask == li, i+1, 0)

    new_mask = np.array(new_mask, dtype='int')
    
    return new_mask

# %% ZONALLY-CYCLIC LABELS FOR MASK OF AR CHAPES
   
def cyclic_labels(mask):
    
    # Find connected regions and assign a label
    labels = measure.label(mask, connectivity=2)
    #plot_map(time[t], lon, lat, labels, cm.batlowW_r)
    
    # Repeat the labeling process with a rolled mask
    mask_roll = np.roll(mask, np.size(mask,1)//2, axis=1)
    labels_roll = measure.label(mask_roll, connectivity=2)
            
    # Roll back
    labels_roll = np.roll(labels_roll, -np.size(mask,1)//2, axis=1)
    #plot_map(time[t], lon, lat, labels_roll, cm.batlowWS)
    
    #plot_ivt_map(lon, lat, np.ma.masked_less(ivt_anom, 0), cm.vik, 0, 1000, 'IVT')

    # Identify border crossings
    left = labels[:, 0]
    right = labels[:, -1]
    left_roll = labels_roll[:, 0]
    right_roll = labels_roll[:, -1]
    
    # Find divided labels
    labels_split = np.unique([left, right])
    
    # Find and isolate connected labels
    labels_correct = np.unique([left_roll, right_roll])
    labels_correct = labels_correct[labels_correct>0]
    labels_replace = np.where(~np.isin(labels_roll, labels_correct), 0, 
                              labels_roll + 2*np.max(labels))

    # Erase divided labels
    cyclic_mask = np.where(np.isin(labels, labels_split), 0, labels)

    # Replace with connected labels
    cyclic_mask = cyclic_mask + labels_replace    
    
    # Relabel mask to avoig gaps
    cyclic_mask = relabel_mask(cyclic_mask)

    return cyclic_mask

# %% IDENTIFY AR SHAPES FROM IVT ANOMALIES

def identify_ARs(ivt_anom, areas, min_area, max_area_hard, fill_radius, fill_disk):
    
    # Initialize variables
    shapes_list = []
    shapes_mask = np.zeros(np.shape(ivt_anom))

    ##### AREA FILTERING 
    
    # Isolate positive anomalies
    ar_shapes = np.ma.where((ivt_anom > 0), 1, 0)

    # Generate labels for AR shapes
    labels = cyclic_labels(ar_shapes)
        
    # Total number of labeled AR shapes
    n = labels.max() + 1
    
    # Areas of AR shapes
    ar_areas = ndimage.sum(areas, labels, np.arange(n))
    
    # Filter min and max area
    sel = np.ones(n, dtype='bool')
    sel = np.where(ar_areas < min_area, False, sel)
    sel = np.where(ar_areas > max_area_hard, False, sel)
    
    # Remove background area
    sel[0] = False   # Why this? Because labels start in 1. Label=0 is the ocean
    
    # Remove AR shapes that don't fulfill area requirements
    keep_labels = np.arange(np.size(sel))[sel]
    ar_shapes = labels*np.isin(labels, keep_labels)
    
    # If no connected regions were indentified, return empty masks
    if ar_shapes.max()==0:
        return shapes_list, shapes_mask
    # Else, re-assign labels
    else:        
        ar_shapes = relabel_mask(ar_shapes)
        
    ##### FILL HOLES IN AR SHAPE

    # For each identified AR shape
    for i in range(ar_shapes.max()):
        
        # Isolate AR
        mask_i = np.ma.where(ar_shapes == i+1,1,0)
            
        # Check zonal cyclicity
        if check_cyclicity(mask_i):                        
            
            # Roll mask
            mask_i = np.roll(mask_i, np.size(ar_shapes,1)//2, axis=1)
            
            # Fill holes 
            mask_i = mask_clossing(mask_i, fill_disk, fill_radius)
            
            # Roll back
            mask_i = np.roll(mask_i, -np.size(ar_shapes,1)//2, axis=1)            

        else:
            # Fill holes
            mask_i = mask_clossing(mask_i, fill_disk, fill_radius)

        shapes_list.append(mask_i)
        shapes_mask = shapes_mask + mask_i

    return shapes_list, shapes_mask

# %% ADD WEIGHTED EDGES TO A GRAPH 

def weighted_edges(graph, nodes, ivt, edge_flux, edge_eps, d):

    # Identify edges along which IVT is flowing
    ratio = np.where(np.isclose(ivt, 0., atol=1e-5), 0, edge_flux/ivt)
    idx = np.where(ratio >= edge_eps, 1, 0)*nodes
    idx = zip(*np.where(idx>0))

    for i, (yi,xi) in enumerate(idx):
        
        # ni0: start, nif: end
        yi = int(yi)
        xi = int(xi)
        ni0 = (yi,xi)
        
        if d=='r' and xi < np.size(ivt,1)-1:
            nif=(yi,xi+1)        
        elif d=='l' and xi > 0:
            nif=(yi,xi-1)
        elif d=='u' and yi < np.size(ivt,0)-1:
            nif=(yi+1,xi)
        elif d=='d' and yi > 0:
            nif=(yi-1,xi)
        elif d=='tr' and yi < np.size(ivt,0)-1 and xi < np.size(ivt,1)-1:
            nif=(yi+1,xi+1)
        elif d=='br' and yi > 0 and xi < np.size(ivt,1)-1:
            nif=(yi-1,xi+1)
        elif d=='tl' and yi < np.size(ivt,0)-1 and xi > 0:
            nif=(yi+1,xi-1)
        elif d=='bl' and yi > 0 and xi > 0:
            nif=(yi-1,xi-1)
        else:
            continue

        mean_flux = edge_flux[yi,xi]

        graph.add_edge(ni0, nif,
                weight = np.exp(-mean_flux/1e2).astype(np.float32),
                ivt = mean_flux.astype(np.float32))
    
    return graph

# %% CREATE A DIRECTED GRAPH FROM AN AR SHAPE MASK

def mask_graph(mask, ivt, ivtu, ivtv, cos_ang, sin_ang, edge_eps):

    # Initialize graph
    graph = nx.DiGraph()
    
    ##### 1 CONNECTIVITY EDGES

    # Prepare nodes
    right = np.roll(mask, -1, axis=1)*mask
    left = np.roll(mask, 1, axis=1)*mask
    up = np.roll(mask, -1, axis=0)*mask
    down = np.roll(mask, 1, axis=0)*mask

    # Add edges
    graph = weighted_edges(graph, right, ivt, ivtu, edge_eps, 'r')
    graph = weighted_edges(graph, left, ivt, -ivtu, edge_eps, 'l')
    graph = weighted_edges(graph, up, ivt, ivtv, edge_eps, 'u')
    graph = weighted_edges(graph, down, ivt, -ivtv, edge_eps, 'd')

    # 2 CONNECTIVITY EDGES
        
    # Prepare nodes
    tr = np.roll(np.roll(mask, -1, axis=0), -1, axis=1)*mask
    br = np.roll(np.roll(mask, 1, axis=0), -1, axis=1)*mask
    tl = np.roll(np.roll(mask, -1, axis=0), 1, axis=1)*mask
    bl = np.roll(np.roll(mask, 1, axis=0), 1, axis=1)*mask

    # Add edges
    graph = weighted_edges(
        graph, tr, ivt, ivtu*cos_ang + ivtv*sin_ang, edge_eps, 'tr')
    graph = weighted_edges(
        graph, br, ivt, ivtu*cos_ang - ivtv*sin_ang, edge_eps, 'br')
    graph = weighted_edges(
        graph, tl, ivt, -ivtu*cos_ang + ivtv*sin_ang, edge_eps, 'tl')
    graph = weighted_edges(
        graph, bl, ivt, -ivtu*cos_ang - ivtv*sin_ang, edge_eps, 'bl')

    return graph

# %% CALCULATE EDGE DISTANCES ALONG A PATH IN A GRAPH

def path_dists(path, attr, graph):
    
    path_sum = 0
    for i in range(len(path)-1):
        
        si = graph.adj[path[i]][path[i+1]][attr] 

        path_sum += si

    return path_sum

# %% IDENTIFY AR AXIS FROM A DIRECTED GRAPH

def graph_axis(graph, ivt, ivtu, ivtv, mask, contour_mask=None):
    
    # Extract nodes
    nodes = list(graph.nodes())

    # Find boundary nodes
    if contour_mask is None:
        contour_mask = mask - morphology.binary_erosion(mask)

    # Calculate the gradients of the mask
    gy, gx = np.gradient(np.array(mask))
        
    # Calculate the transboundary flux
    trans_flux = (gx*ivtu + gy*ivtv)*contour_mask
    
    # Find boundary pixels with net input moisture fluxes
    in_node_coords = np.where(trans_flux > 0)
    in_node_coords = zip(in_node_coords[0], in_node_coords[1])
    in_node_coords = list(set(in_node_coords).intersection(nodes))
    n_in=len(in_node_coords)
    
    # Find boundary pixels with net output moisture fluxes
    out_node_coords = np.where(trans_flux < 0)
    out_node_coords = zip(out_node_coords[0], out_node_coords[1])
    out_node_coords = list(set(out_node_coords).intersection(nodes))
    n_out=len(out_node_coords)

    # When mask is at edge of the map. Rarely happens.
    if n_in==0:
        in_node_coords = nodes
        n_in = len(in_node_coords)
    if n_out==0:
        out_node_coords = nodes
        n_out = len(out_node_coords)

    # Initialize variable
    dists = np.zeros((n_in,n_out))

    ##### FIND HEAVIEST PATH
    
    for i in range(n_in):
        
        # Starting income node
        ei = in_node_coords[i]
        
        # Find all possible paths
        paths_i = nx.single_source_dijkstra_path(graph, ei, weight='weight')
        
        # Extract the paths ending in an output node
        paths_i = dict([(k,v) for k,v in paths_i.items() if k in out_node_coords])
        
        # Calculate the weight of each path and choose the heaviest one
        if len(paths_i)>0:
            
            dist_dict = dict([(k, path_dists(v,'ivt', graph)) for k,v in paths_i.items()])
            node_i = sorted(dist_dict, key=dist_dict.get)[-1]
            dist_i = dist_dict[node_i]
            dists[i, out_node_coords.index(node_i)] = dist_i

    if np.max(dists)==0:
        # this may happen when a mask is touching the map edges, and in_node_coords
        # out_node_coords can't be linked by a path. Very rarely happen, but damn
        # annoying. A fallback solution is to use an undirected graph linking
        # the most inward and most outward pixels.
        most_in = np.unravel_index(np.argmax(trans_flux), np.shape(mask))
        most_out = np.unravel_index(np.argmin(trans_flux), np.shape(mask))
        graph_und = graph.to_undirected()
        
        try:
            path = nx.dijkstra_path(graph_und, most_in, most_out, weight='weight')
        except:
            # if it still can't find a path, make a full connected network
            graph_full = mask_graph(mask, ivt, ivtu, ivtv, np.ones(np.shape(mask)), 
                                    np.ones(np.shape(mask)), -np.inf)

            path = nx.dijkstra_path(graph_full, most_in, most_out, weight='weight')
    else:
        # Choose the heaviest path of the network
        max_idx = np.argmax(dists)
        yi_dx, xi_dx = np.unravel_index(max_idx,(n_in,n_out))
        path = nx.dijkstra_path(graph, in_node_coords[yi_dx], out_node_coords[xi_dx], weight='weight')

    # Get a mask for the AR axis
    axis_mask = np.zeros(np.shape(mask))
    
    for (y,x) in path:
        axis_mask[y,x] = 1

    path = np.array(path)

    return path, axis_mask

# %% IDENTIFY AXES OF AR SHAPES

def identify_AR_axis(shapes_list, ivt, ivtu, ivtv, 
                     ivt_roll, ivtu_roll, ivtv_roll, 
                     cos_ang, sin_ang, cos_ang_roll, sin_ang_roll, 
                     edge_eps, down_sampling):

    # Initialize variables
    axis_list=[]
    axis_mask = np.zeros(np.shape(ivtu))

    # Find axes
    for i in range(np.size(shapes_list,0)):
        
        # Extract mask with AR shape
        mask_i = shapes_list[i]

        # Check if AR shape is zonally cyclic
        if check_cyclicity(mask_i):
            roll_i = True
            ivt_i = ivt_roll
            ivtu_i = ivtu_roll
            ivtv_i = ivtv_roll
            cos_ang_i = cos_ang_roll
            sin_ang_i = sin_ang_roll
            mask_i = np.roll(mask_i, np.size(mask_i,1)//2, axis=1)
            
        else:
            roll_i = False
            ivt_i = ivt
            ivtu_i = ivtu
            ivtv_i = ivtv
            cos_ang_i = cos_ang
            sin_ang_i = sin_ang

        ds_path_success = None
        
        if down_sampling > 1:
            # Try quicker graph path search first
            
            # Down sample the AR shape
            mask_i_ds = mask_i[::down_sampling, ::down_sampling]
            ivt_i_ds = ivt_i[::down_sampling, ::down_sampling]
            ivtu_i_ds = ivtu_i[::down_sampling, ::down_sampling]
            ivtv_i_ds = ivtv_i[::down_sampling, ::down_sampling]
            cos_ang_i_ds = cos_ang_i[::down_sampling, ::down_sampling]
            sin_ang_i_ds = sin_ang_i[::down_sampling, ::down_sampling]

            # Find axis for down sampled AR shape
            try:
                # Create flux graph
                graph_i_ds = mask_graph(mask_i_ds, ivt_i_ds, ivtu_i_ds, ivtv_i_ds, 
                                        cos_ang_i_ds, sin_ang_i_ds, edge_eps)
                
                # Find AR axes
                axis_arr_i_ds, axis_mask_i_ds = graph_axis(
                    graph_i_ds, ivt_i_ds, ivtu_i_ds, ivtv_i_ds, mask_i_ds, None)
                
            except:
                ds_path_success = False
            else:
                ds_path_success = True

                # Get AR axis end points
                p0 = axis_arr_i_ds[0]*down_sampling
                pf = axis_arr_i_ds[-1]*down_sampling

                # Contour points to search from
                contour_mask = np.zeros_like(mask_i)
                
                # Enlarge the search region a bit: down_sampling*3
                contour_mask[max(0, p0[0]-down_sampling*3):min(np.size(mask_i,0), p0[0]+down_sampling*3),
                          max(0, p0[1]-down_sampling*3):min(np.size(mask_i,1), p0[1]+down_sampling*3)] = 1
                contour_mask[max(0, pf[0]-down_sampling*3):min(np.size(mask_i,0), pf[0]+down_sampling*3),
                          max(0, pf[1]-down_sampling*3):min(np.size(mask_i,1), pf[1]+down_sampling*3)] = 1

                # Isolate the contour
                contour_mask = contour_mask * (mask_i - morphology.binary_erosion(mask_i))

                # Create flux graph
                graph_i = mask_graph(mask_i, ivt_i, ivtu_i, ivtv_i, 
                                     cos_ang_i, sin_ang_i, edge_eps)

                # Get AR axis from graph
                axis_arr_i, axis_mask_i = graph_axis(
                    graph_i, ivt_i, ivtu_i, ivtv_i, mask_i, contour_mask)


        # If quicker graph search failed, or downsampling was not done:
        if down_sampling == 1 or ds_path_success == False:
            
            # Convert mask to directed graph
            graph_i = mask_graph(
                mask_i, ivt_i, ivtu_i, ivtv_i, cos_ang_i, sin_ang_i, edge_eps)
            
            # Get AR axis from graph
            axis_arr_i, axis_mask_i = graph_axis(
                graph_i, ivt_i, ivtu_i, ivtv_i, mask_i, None)

        # Roll back
        if roll_i:
            axis_mask_i = np.roll(axis_mask_i, -np.size(mask_i,1)//2, axis=1)    
            new_lon = (axis_arr_i[:,1] - np.size(mask_i,1)//2)%np.size(axis_mask_i,1)
            axis_arr_i[:,1] = new_lon

        axis_list.append(axis_arr_i)
        axis_mask = axis_mask + axis_mask_i

    return axis_list, axis_mask

# %% UV DECOMPOSITION OF BACKGROUND AND ANOMALOUS IVT

def uv_decomp(ivt, ivtu, ivtv, ivt_rec, ivt_anom):
    
    # U components
    ivtu_rec = ivtu * ivt_rec / ivt
    ivtu_anom = ivtu * ivt_anom / ivt
    
    # V components
    ivtv_rec = ivtv * ivt_rec / ivt
    ivtv_anom = ivtv * ivt_anom / ivt
    
    return ivtu_rec, ivtu_anom, ivtv_rec, ivtv_anom

# %% GET AR CONTOUR FROM AR SHAPE MASK

def get_contour(mask, lon, lat):
        
    # Extract contours
    cs = plt.contourf(lon, lat, mask, [0.9,1.1]).allsegs
    plt.close()
    
    # Choose the longest contour
    ar_cont = max((segment for level in cs for segment in level), key=len, default=None)

    return ar_cont

# %% EXTRACT GEOMETRICAL PROPERTIES OF ARs

def get_AR_data(ivt, ivt_roll, areas, lat, lon, lon_roll, shapes_list, axis_list, 
                time_str, min_area):

    # Initialize variables
    labels = np.zeros(np.shape(ivt))
    ar_prop = {}

    # For each identified AR shape
    for i in range(len(shapes_list)):
        
        # Get AR mask
        mask_i = shapes_list[i]
        
        # Calculate AR area
        area_i = (mask_i*areas).sum()  # km^2, insensitive to zonal cyclicity 

        # Delete too small ARs
        if area_i < min_area:
            continue
        
        # Extract the axis coordinates
        skel_i = axis_list[i]
        ax_lat_i = lat[skel_i[:,0]]
        ax_lon_i = lon[skel_i[:,1]]
        axis_i = np.c_[ax_lat_i, ax_lon_i]

        # Compute AR axis length
        lens = geo_dist(axis_i[:,1], axis_i[:,0])
        len_i = np.sum(lens) #km

        # Delete if too short (this is a first filter to delete noisy identification)
        if len_i < 1200:
            continue

        # Check if AR shape is zonally cyclic
        if check_cyclicity(mask_i):   
            
            # Roll mask
            roll_i = True
            lon_i = lon_roll
            ivt_i = ivt_roll
            mask_i = np.roll(mask_i, np.size(ivt,1)//2, axis=1)
            
        else:
            roll_i = False
            lon_i = lon
            ivt_i = ivt

        # Calculate AR properties, in pixel units
        arp_i = measure.regionprops(mask_i, intensity_image=np.array(ivt_i))[0]
        
        # Extract AR centroid
        centroid_lat, centroid_lon = arp_i.weighted_centroid
        centroid_lat = lat[int(centroid_lat)]
        centroid_lon = lon_i[int(centroid_lon)]        

        # Extract AR contour
        cont_i = get_contour(mask_i, lon_i, lat)
        
        # Roll back to a [-180,180] longitudinal grid
        if roll_i:
            mask_i = np.roll(mask_i, -np.size(ivt,1)//2, axis=1) 
            centroid_lon = ((centroid_lon + 180) % 360) - 180
            cont_i[:,0] = ((cont_i[:,0] + 180) % 360) - 180
        
        # Label AR shape
        labels = labels + mask_i*(i+1)

        # Create pandas dataframe with properties
        ar_prop[i+1]={
            'id': i+1,
            'time': time_str,
            'contour_lon': cont_i[:,0],
            'contour_lat': cont_i[:,1],
            'centroid_lon': centroid_lon,
            'centroid_lat': centroid_lat,
            'axis_lon': axis_i[:,1],
            'axis_lat': axis_i[:,0],
            }

    return labels, ar_prop

# %% READ CSV RECORD

def read_csv_record(file_name):

    def conv_array(text):
        '''Convert array texts to ndarray'''
        text = text.replace('[','').replace(']','')
        array = np.array(text.split()).astype('float')
        return array

    conv_keys = ['contour_lon', 'contour_lat', 'axis_lon', 'axis_lat']

    converters = dict([(key_i, conv_array) for key_i in conv_keys])

    dtypes={'id': 'int', 
            'time': 'str',
            'deleted': 'str',
            'is_relaxed': 'bool',
            'why_relaxed': 'str'}

    ardf = pd.read_csv(file_name, dtype=dtypes, converters=converters)

    return ardf

# %% AREA AVERAGE OF VARIABLE

def area_average(flux_masked, ar_mask, areas):
    '''Compute area-weighted average.'''

    result = np.sum(flux_masked * areas)/np.sum(ar_mask * areas)

    return result

# %% COMPUTE AR RANK BASED ON THE AR SCALE BY RALPH ET AL., 2019

def compute_ar_rank(dur, intens):
    """
    Assigns AR rank based on duration and strength.
    
    Parameters
    ----------
    dur: int
    AR duration in hours
    intens: float
    AR intensity/max. IVT
    Returns
    ------- 
    rank: AR rank
    """
    ## Define masks
    l_mask_cat = []
    # CATEGORY 1
    tmp_cond1 = (((dur > 24) & (dur <= 48)) & ((intens > 250) & (intens <= 500)))
    tmp_cond2 = (((dur > 0) & (dur <= 24)) & ((intens > 500) & (intens <= 750)))
    l_mask_cat.append(tmp_cond1 | tmp_cond2)
    # CATEGORY 2
    tmp_cond1 = (((dur > 48)) & ((intens > 250) & (intens <= 500))) 
    tmp_cond2 = (((dur > 24) & (dur <= 48)) & ((intens > 500) & (intens <= 750))) 
    tmp_cond3 = (((dur > 0) & (dur <= 24)) & ((intens > 750) & (intens <= 1000))) 
    l_mask_cat.append(tmp_cond1 | tmp_cond2 | tmp_cond3)
    # CATEGORY 3
    tmp_cond1 = (((dur > 48)) & ((intens > 500) & (intens <= 750))) 
    tmp_cond2 = (((dur > 24) & (dur <= 48)) & ((intens > 750) & (intens <= 1000))) 
    tmp_cond3 = (((dur > 0) & (dur <= 24)) & ((intens > 1000) & (intens <= 1250))) 
    l_mask_cat.append(tmp_cond1 | tmp_cond2 | tmp_cond3)
    # CATEGORY 4
    tmp_cond1 = (((dur > 48)) & ((intens > 750) & (intens <= 1000))) 
    tmp_cond2 = (((dur > 24) & (dur <= 48)) & ((intens > 1000) & (intens <= 1200))) 
    tmp_cond3 = (((dur > 0) & (dur <= 24)) & (intens > 1250)) 
    l_mask_cat.append(tmp_cond1 | tmp_cond2 | tmp_cond3)
    # CATEGORY 5
    tmp_cond1 = (((dur > 48)) & ((intens > 1000) & (intens <= 1250))) 
    tmp_cond2 = (((dur > 24) & (dur <= 48)) & (intens > 1250)) 
    tmp_cond3 = ((dur > 48) & (intens > 1250)) 
    l_mask_cat.append(tmp_cond1 | tmp_cond2 | tmp_cond3)

    # depending on criterion, dimensions need to be modified
    #if np.ndim(l_mask_cat[0]) == 0:
    #    for i in range(5):
    #        l_mask_cat[i] = np.hstack([l_mask_cat[i]])
    
    # match category
    idx = np.where([l_mask_cat[i] for i in range(5)])[0]
    if len(idx) > 0:
        rank = idx[0] + 1
    else:
        rank=0
    
    return rank

# %% FIND SEQUENCES OF ONES IN A BINARY TIME SERIES

def rle(array):
    """
    Run-length encoding for an array.
    
        Parameters
    ----------
    array : array-like
        Input array to be encoded.
    
    Returns
    -------
    values : ndarray
        Array containing unique values in the input array.
    counts : ndarray
        Array containing run lengths for each unique value.
    """
    array = np.asarray(array)
    changes = np.where(array[:-1] != array[1:])[0] + 1
    
    indices = np.concatenate(([0], changes, [len(array)]))
    counts = np.diff(indices)
    values = array[indices[:-1]]
    
    return values, counts

# %% CLASSIFY AR TRAJECTORY AS INLAND-PENETRATING

a_continents = np.array(['north_america', 'south_america', 'europe', 'africa', 
                         'asia', 'australia', 'oceania', 'antarctica'])
l_contlats = ['conti_na_lat', 'conti_sa_lat', 'conti_eu_lat', 'conti_af_lat', 
              'conti_as_lat', 'conti_au_lat', 'conti_oc_lat', 'conti_an_lat']

def inland_penetration(ARcont):
    """
    Determines if an Atmospheric River (AR) is penetrating inland.
    
    Parameters
    ----------
    ARcont : DataFrame
        DataFrame containing information about the atmospheric river.
    
    Returns
    ------- 
    IP : bool
        True if the AR is penetrating inland, False otherwise.
    """
    l = ARcont.time.size
    # isolate continental from insular component
    tmp_mask = np.zeros(l, dtype=bool)
    for j in range(l):
        # what's the current continent?
        idx = np.where(a_continents == ARcont.continent.values[j])[0][0]
        # only grab non-insular li-latitudes and check if they're all NAN
        lat = ARcont[l_contlats[idx]].values[j]
        # if they're not, AR intersects continent
        tmp_mask[j] = ~np.isnan(lat)
    
    # compute landfalling fraction
    a_frac = ARcont['land'].values[tmp_mask]
    #ARcont[a_continents].max(axis=1).values[tmp_mask]

    # trajectory can only be inland penetrating if the section after landfall is still longer than two instances
    n = a_frac.size
    if  n >= 2:
        # is it increasing?
        slope, _ = np.polyfit(np.arange(n), a_frac, 1)
        if slope > 0:
            IP = True
        else:
            IP = False
    else:
        IP = False
    return IP

# %% POPULATE VARIABLES FOR NETCDF FORMAT

def populate_var_4D(ardf, time, var, var_specif):
    
    # Extract track IDs
    tracks = ardf.drop_duplicates('trackid')
    track_id = np.array(tracks['trackid'])
    
    # Intialize variable
    variable = np.full(var_specif[0], np.nan, dtype=var_specif[1])
    assert np.ndim(variable)==4, '%s has an incorrect shape' %var.upper()
    
    # Populate variable
    for i in range(len(ardf)):
        
        # Time mask
        t = (time == ardf['time'].iloc[i])
        
        # Track ID mask
        j = (track_id == ardf['trackid'].iloc[i])
            
        # Populate variable
        variable[t, j, 1, :ardf['%s_len' %var].iloc[i]] = np.array(
            ardf['%s_lon' %var].iloc[i])
        variable[t, j, 0, :ardf['%s_len' %var].iloc[i]] = np.array(
            ardf['%s_lat' %var].iloc[i])
        
    return variable
        
def populate_var_3D_coord(ardf, time, var, var_specif):
    
    # Extract track IDs
    tracks = ardf.drop_duplicates('trackid')
    track_id = np.array(tracks['trackid'])
    
    # Intialize variable
    variable = np.full(var_specif[0], np.nan, dtype=var_specif[1])
    assert np.ndim(variable)==3, '%s has an incorrect shape' %var.upper()
    
    # Populate variable
    for i in range(len(ardf)):
        
        # Time mask
        t = (time == ardf['time'].iloc[i])
        
        # Track ID mask
        j = (track_id == ardf['trackid'].iloc[i])
            
        # Populate variable
        variable[t, j, 1] = ardf['%s_lon' %var].iloc[i]
        variable[t, j, 0] = ardf['%s_lat' %var].iloc[i]
        
    return variable        

def populate_var_3D_comp(ardf, time, var, var_specif):
    
    # Extract track IDs
    tracks = ardf.drop_duplicates('trackid')
    track_id = np.array(tracks['trackid'])
    
    # Intialize variable
    variable = np.full(var_specif[0], np.nan, dtype=var_specif[1])
    assert np.ndim(variable)==3, '%s has an incorrect shape' %var.upper()
    
    # Populate variable
    for i in range(len(ardf)):
        
        # Time mask
        t = (time == ardf['time'].iloc[i])
        
        # Track ID mask
        j = (track_id == ardf['trackid'].iloc[i])
            
        # Populate variable
        variable[t,j,0] = ardf[var].iloc[i] 
        variable[t,j,1] = ardf['%su' %var].iloc[i] 
        variable[t,j,2] = ardf['%sv' %var].iloc[i] 
        
    return variable 

def populate_var_2D(ardf, time, var, var_specif, fill_value=np.nan):
    
    # Extract track IDs
    tracks = ardf.drop_duplicates('trackid')
    track_id = np.array(tracks['trackid'])
    
    # Intialize variable
    variable = np.full(var_specif[0], fill_value, dtype=var_specif[1])
    assert np.ndim(variable)==2, '%s has an incorrect shape' %var.upper()
    
    # Populate variable
    for i in range(len(ardf)):
        
        # Time mask
        t = (time == ardf['time'].iloc[i])
        
        # Track ID mask
        j = (track_id == ardf['trackid'].iloc[i])
            
        # Populate variable
        variable[t, j] = ardf[var].iloc[i]
        
    return variable

def populate_var_1D(ardf, var, var_specif, fill_value=np.nan):
    
    # Extract track IDs
    tracks = ardf.drop_duplicates('trackid')
    track_id = np.array(tracks['trackid'])
    
    # Intialize variable
    variable = np.full(var_specif[0], fill_value, dtype=var_specif[1])
    assert np.ndim(variable)==1, '%s has an incorrect shape' %var.upper()
    
    # Populate variable
    for j in range(len(track_id)):
        
        # Extract track
        track = ardf.loc[ardf['trackid'] == track_id[j], var]
        
        # Get unique values
        unique_value = np.unique(track)
        
        # Check that all values are the same
        assert np.ndim(unique_value)==1, \
            '%s has different values for track %s' %(var.upper(), track_id)
            
        # Populate variable
        variable[j] = unique_value[0]
        
    return variable 

def populate_var_string(ardf, time, var, var_specif, max_string_size):
    
    # Extract track IDs
    tracks = ardf.drop_duplicates('trackid')
    track_id = np.array(tracks['trackid'])
    
    # Intialize variable
    variable = np.full(var_specif[0], '', dtype=var_specif[1])
    assert np.ndim(variable)==2, '%s has an incorrect shape' %var.upper()
    
    # Populate variable
    for i in range(len(ardf)):
        
        # Time mask
        t = (time == ardf['time'].iloc[i])
        
        # Track ID mask
        j = (track_id == ardf['trackid'].iloc[i])
            
        # Populate variable
        variable[t, j] = ardf[var].iloc[i]
        
    # Convert to numpy fixed-width unicode array
    variable = np.array(variable, dtype='S%s' %max_string_size)
    
    # Convert to a char array (adds a third dimension: character index)    
    variable = stringtochar(variable)
            
    return variable
