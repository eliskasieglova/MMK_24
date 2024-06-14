import arcpy
import numpy as np
from math import *
import geopandas as gpd
import pandas as pd
import os
from pathlib import Path
from map_projections import projs
import math


if not Path('results').exists():
    os.mkdir('results')

if not Path('temp').exists():
    os.mkdir('temp')


# Set workspace - SET YOUR OWN WORKSPACE
arcpy.env.workspace = 'PATH TO YOUR WORKSPACE'
arcpy.env.overwriteOutput = 1

# Create Project - MAKE PROJECT IN ARCGIS PRO CALLED "project.aprx" AND SAVE IN CURRENT DIRECTORY
aprx = arcpy.mp.ArcGISProject("project.aprx")

# Set page width and page height (mm)
pw = 300
ph = 430

pw = 300
ph = 430

# Create Layout
layout = aprx.createLayout(pw, ph, 'MILLIMETER')
layout.name = "Layout"

# antarctica, arctic, africa, australia, australia2, oceania, south america, africa, asia, japan, pacific, america
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
rotace = [-180, -216, 144, -288, 0, -72, -144, -216, -144, -72, 0, 72]
posuny_x = [271, 202, 242.8, 185, 190.1, 251, 283.65, 247.5, 116, 30.4, 108.65, 242.7]
posuny_y = [370, 159, 283.8, 307.9, 370.2, 384.5, 331.3, 229.5, 277, 167, 51.5, 90.2]


def mapframe_vertices(x_center, y_center, edge_length, theta_0):
    """
    Create vertices of map frame.

    :param x_center: x-coord of center point of pentagon
    :param y_center: y-coord of center point of pentagon
    :param edge_length: length of one edge
    :param theta_0: rotation

    :return: list of tuples (x,y)-coordinates of each vertex
    """

    # Calculate the circumradius R
    R = edge_length / (2 * math.sin(math.pi / 5))

    # Calculate the vertices
    vertices = []
    for i in range(5):
        theta = theta_0 + i * 2 * math.pi / 5
        x = x_center + R * math.cos(theta)
        y = y_center + R * math.sin(theta)
        vertices.append((x, y))

    return vertices


def graticule(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, R, uk, vk, u0, proj):
    # Create meridians
    XM = []
    YM = []

    for v in np.arange(v_min, v_max + D_v, D_v):
        # Meridian
        um = np.arange(u_min, u_max + d_u, d_u)
        m = len(um)
        vm = np.repeat(v, m)

        # Convert to oblique aspect
        sm, dm = uv_to_sd(um, vm, uk, vk)

        # Project a meridian
        xm, ym = proj(R, sm, dm)

        # Add meridian
        XM.extend(xm)
        YM.extend(ym)

    # Create parallels
    XP = []
    YP = []
    for u in np.arange(u_min, u_max + D_u, D_u):
        # Longitude of parallel points
        vp = np.arange(v_min, v_max + d_v, d_v)
        n = len(vp)
        up = np.repeat(u, n)

        # Convert to oblique aspect
        sp, dp = uv_to_sd(up, vp, uk, vk)

        # Project a parallel
        xp, yp = proj(R, sp, dp)

        # Add parallel
        XP.extend(xp)
        YP.extend(yp)

    return XM, YM, XP, YP


def uv_to_sd(u, v, uk, vk):

    dv = vk - v
    # Transformed latitude
    s = np.arcsin(np.sin(u) * np.sin(uk) + np.cos(u) * np.cos(uk) * np.cos(dv))

    # Transformed longitude
    d = -1 * np.arctan2(np.cos(u) * np.sin(dv), np.cos(u) * np.sin(uk) * np.cos(dv) - np.sin(u) * np.cos(uk))

    return s, d


def gnom(R, s, d):
    x = R * np.tan(np.pi/2 - s) * np.cos(d)
    y = R * np.tan(np.pi/2 - s) * np.sin(d)
    return x, y


def boundary(u, v, R, uk, vk, u0, proj):
    XB, YB = [], []
    for u, v in zip(u, v):
        s, d = uv_to_sd(u, v, uk, vk)
        dv = vk - v
        XB_, YB_ = proj(R, s, d)

        XB.append(XB_)
        YB.append(YB_)

    return XB, YB

def rotate(X, Y, angle):
    # Convert angle to radians
    angle_rad = np.radians(angle)
    
    # Define the rotation matrix
    rotation_matrix = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],
                                 [np.sin(angle_rad), np.cos(angle_rad)]])
    
    # Convert boundary points to numpy array for easier manipulation
    boundary_points = np.array([X, Y])
    
    # Apply the rotation matrix to the boundary points
    rotated_boundary_points = np.dot(rotation_matrix, boundary_points)
    
    # Extract rotated X and Y coordinates
    rotated_X = rotated_boundary_points[0]
    rotated_Y = rotated_boundary_points[1]
    
    return rotated_X, rotated_Y


def splitMeridians(XM, YM):

    meridians = []

    temp_meridian = []
    xm_current = XM[0]
    for i in range(0, len(XM)):
        if (xm_current != XM[i]) & (i != 0):
            meridians.append(temp_meridian)
            temp_meridian = []
            temp_meridian.append((XM[i], YM[i]))

        else:
            temp_meridian.append((XM[i], YM[i]))
            xm_current = XM[i]

    return meridians

def splitParallels(XP, YP):

    parallels = []

    temp_parallel = []
    yp_current = YP[0]
    for i in range(0, len(YP)):
        if (yp_current != YP[i]) & (i != 0):
            parallels.append(temp_parallel)
            temp_parallel = []
            temp_parallel.append((XP[i], YP[i]))

        else:
            temp_parallel.append((XP[i], YP[i]))
            yp_current = YP[i]

    return parallels

us1 = 52.6226
uj1 = -us1
us2 = 10.8123
uj2 = -us2
ualpha = 26.5651

u_min = [-90, 30, -60, -60, -60, -60, -60, -20, -20, -20, -20, -20]
u_max = [-30, 90, 20,   20, 20,   20, 20,   60, 60, 60, 60, 60]
v_min = [-180, -180, -10, 60, 130, 210, 270,-40, 30, 100, 170, 240]
v_max = [180, 180, 80, 150, 225, 290,   370, 40, 120, 190, 260, 330]
 
u_k = [-90, 90, (uj1+us2)/2, (uj1+us2)/2, (uj1+us2)/2, (uj1+us2)/2, (uj1+us2)/2,
      (us1-uj2)/2, (us1-uj2)/2, (us1-uj2)/2, (us1-uj2)/2, (us1-uj2)/2]

v_k = [0, 0, 36, (144-72)/2+72, 180, 216+36, 288+36, 0, 72, 144, 216, 288]
u_b = [[uj1, uj1, uj1, uj1, uj1, uj1],
      [us1, us1, us1, us1, us1, us1],
      [uj1, uj1, uj2, us2, uj2, uj1],
      [uj1, uj1, uj2, us2, uj2, uj1],
      [uj1, uj1, uj2, us2, uj2, uj1],
      [uj1, uj1, uj2, us2, uj2, uj1],
      [uj1, uj1, uj2, us2, uj2, uj1],
      [uj2, us2, us1, us1, us2, uj2],
      [uj2, us2, us1, us1, us2, uj2],
      [uj2, us2, us1, us1, us2, uj2],
      [uj2, us2, us1, us1, us2, uj2],
      [uj2, us2, us1, us1, us2, uj2]]
v_b = [[0, 72, 144, 216, 288, 0],
      [36, 108, 180, 252, 324, 36],
      [0, 72, 72, 36, 0, 0],
      [72, 144, 144, 108, 72, 72],
      [144, 216, 216, 180, 144, 144],
      [216, 288, 288, 252, 216, 216],
      [288, 360, 360, 324, 288, 288],
      [0, 36, 36, 324, 324, 0],
      [72, 108, 108, 36, 36, 72],
      [144, 180, 180, 108, 108, 144],
      [216, 252, 252, 180, 180, 216],
      [288, 324, 324, 252, 252, 288]]
ylim = [-5000, 5000]

for i in range(12):
    print(i)
    umin = u_min[i] * np.pi/180
    umax = u_max[i] * np.pi/180
    vmin = v_min[i] * np.pi/180
    vmax = v_max[i] * np.pi/180
    Du = 10 * np.pi / 180
    Dv = Du
    du = np.pi / 180
    dv = du
    R = 6378000
    uk = u_k[i] * np.pi/180
    vk = v_k[i] * np.pi/180
    proj = gnom
    u0 = 0
    ub = [num * np.pi/180 for num in u_b[i]]
    vb = [num * np.pi/180 for num in v_b[i]]

    # create boundary face shapefile if it does not exist
    if not Path(f'temp/rotated_boundary_face{i+1}.shp').is_file():
        XB, YB = boundary(ub, vb, R, uk, vk, u0, proj)
        XM, YM, XP, YP = graticule(umin, umax, vmin, vmax, Du, Dv, du, dv, R, uk, vk, u0, proj)

        # Rotate boundary points by -90 degrees
        rotated_XB, rotated_YB = rotate(XB, YB, -90)

        # Create polyline from rotated boundary points
        boundary_pts = arcpy.Array([arcpy.Point(x, y) for x, y in zip(rotated_XB, rotated_YB)])

        # Save the rotated boundary as a shapefile
        output = arcpy.Polyline(boundary_pts)
        arcpy.CopyFeatures_management(output, f'temp/rotated_boundary_face{i+1}.shp')

        rotated_XM, rotated_YM = rotate(XM, YM, -90)
        rotated_XP, rotated_YP = rotate(XP, YP, -90)


    # create shapefile of meridians if it does not exist
    if not Path(f'temp/meridians_face{i + 1}_split.shp').is_file():
        # generate meridians
        meridians = arcpy.Array()
        for x in range(len(rotated_XM)):
            meridians.add(arcpy.Point(rotated_XM[x], rotated_YM[x]))

        output = arcpy.Polyline(meridians)
        arcpy.CopyFeatures_management(output, f'temp/meridians_face{i + 1}_notsplit.shp')

        arcpy.management.SplitLine(f'temp/meridians_face{i+1}_notsplit.shp', f'temp/meridians_face{i+1}_split.shp')

        with arcpy.da.UpdateCursor(f"temp/meridians_face{i+1}_split.shp", ["FID", "SHAPE@LENGTH"]) as cursor:
            for row in cursor:
                if row[1] > 1000000:
                    cursor.deleteRow()

    # create shapefile of parallels if it does not exist
    if not Path(f'temp/parallels_face{i+1}_split.shp').is_file():

        # generate parallels
        parallels = arcpy.Array()
        for y in range(len(rotated_XP)):
            parallels.add(arcpy.Point(rotated_XP[y], rotated_YP[y]))

        output = arcpy.Polyline(parallels)
        arcpy.CopyFeatures_management(output, f'temp/parallels_face{i + 1}_notsplit.shp')

        arcpy.management.SplitLine(f'temp/parallels_face{i+1}_notsplit.shp', f'temp/parallels_face{i+1}_split.shp')

        with arcpy.da.UpdateCursor(f"temp/parallels_face{i+1}_split.shp", ["FID", "SHAPE@LENGTH"]) as cursor:
            for row in cursor:
                if row[1] > 1000000:
                    cursor.deleteRow()

    # Create Map
    new_map = aprx.createMap(f"Face{i+1}")

    # Add basemap territory (CHOOSE ANY BASEMAP FROM https://www.arcgis.com/home/group.html?id=702026e41f6641fb85da88efe79dc166&view=list#content)
    new_map.addBasemap("Charted Territory Map")

    # set target spatial reference
    spatial_reference = projs[i]

    # reproject meridians to given spatial_reference
    if not Path(f'results/meridians_face{i+1}.shp').is_file():
        # open and project shapefile with geopandas - arcpy gave us fatalerror on arcpy.Project_management()
        m_gdf = gpd.read_file(f"temp/meridians_face{i+1}_split.shp").set_crs(spatial_reference)
        m_gdf.to_file(f'results/meridians_face{i + 1}.shp')

    # reproject parallels to given spatial_reference
    if not Path(f'results/parallels_face{i+1}.shp').is_file():
        p_gdf = gpd.read_file(f"temp/parallels_face{i+1}_split.shp").set_crs(spatial_reference)
        p_gdf.to_file(f'results/parallels_face{i+1}.shp')

    # reproject boundary to given spatial_reference
    if not Path(f'results/boundary_face{i + 1}.shp').is_file():
        p_gdf = gpd.read_file(f"temp/rotated_boundary_face{i + 1}.shp").set_crs(spatial_reference)
        p_gdf.to_file(f'results/boundary_face{i + 1}.shp')

    # make feature layers from projected shapefiles
    m_path = f"results/meridians_face{i+1}.shp"
    p_path = f"results/parallels_face{i+1}.shp"
    b_path = f"results/boundary_face{i+1}.shp"

    meridians_shp = arcpy.management.MakeFeatureLayer(m_path, "meridians")[0]
    parallels_shp = arcpy.management.MakeFeatureLayer(p_path, "parallels")[0]
    boundary_shp = arcpy.management.MakeFeatureLayer(b_path, "boundary")[0]

    # add meridians and parallels to map
    new_map.addLayer(meridians_shp)
    new_map.addLayer(parallels_shp)

    # input parameters for counting vertices of map frame
    dx = posuny_x[i]  # x-coordinate of center point (set at the beginning of code)
    dy = posuny_y[i]  # y-coordinate of center point (set at the beginning of code)
    rot = rotace[i]

    # count vertices
    edge_length = 50  # edge length
    center = (dx, dy)
    theta_0 = math.pi/10  # rotation of map frames for northern hemisphere
    if i > 6:
        theta_0 = 19*math.pi/10  # rotation of map frames for southern hemisphere

    # count vertices of map frame
    vertices = mapframe_vertices(center[0], center[1], edge_length, theta_0)

    p1 = (vertices[0])
    p2 = (vertices[1])
    p3 = (vertices[2])
    p4 = (vertices[3])
    p5 = (vertices[4])

    # Set boundary points of map frame
    map_boundary_points = arcpy.Array([
        arcpy.Point(p1[0], p1[1]),
        arcpy.Point(p2[0], p2[1]),
        arcpy.Point(p3[0], p3[1]),
        arcpy.Point(p4[0], p4[1]),
        arcpy.Point(p5[0], p5[1]),
        arcpy.Point(p1[0], p1[1])
    ])

    # Create the polygon of map frame boundary
    map_boundary = arcpy.Polygon(map_boundary_points)

    # Create the map frame
    map_frame = layout.createMapFrame(map_boundary, new_map)
    map_frame.map = new_map
    map_frame.setAnchor = "CENTER_POINT"

    # Edit symbology
    layer1 = new_map.listLayers('meridians')[0]
    layer2 = new_map.listLayers('parallels')[0]

    symbology1 = layer1.symbology
    symbology2 = layer2.symbology

    # Set the symbol color to black
    symbology1.renderer.symbol.color = {'RGB': [0, 0, 0, 20]}
    symbology1.renderer.symbol.width = 1
    symbology2.renderer.symbol.color = {'RGB': [0, 0, 0, 20]}
    symbology2.renderer.symbol.width = 1

    # Update the layer symbology
    layer1.symbology = symbology1
    layer2.symbology = symbology2

    # set projection of background map
    spatial_ref = arcpy.SpatialReference(text=projs[i])
    map_frame.map.spatialReference = spatial_ref

    # set extent to fit the boundary face
    desc = arcpy.Describe(f'results/boundary_face{i+1}.shp')
    extent = desc.extent
    map_frame.camera.setExtent(extent)
    map_frame.camera.heading = rot
    map_frame.elementRotation = rot

# Save the changes to the project
#aprx.save()

# Export the layout to a PDF
print('exporting')
layout.exportToPDF("globe_faces.pdf", resolution=600)

