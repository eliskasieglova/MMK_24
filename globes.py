import arcpy
import numpy as np
from pathlib import Path
import os

# set workspace
arcpy.env.workspace = 'C:/Users/eliss/documents/pgis/MATKARTO'
arcpy.env.overwriteOutput = 1

# create folder with the name 'temp' for storing temporary data and 'results' for storing results
tempfolder = Path('temp/')
if not tempfolder.exists():
    os.mkdir(tempfolder)
resultsfolder = Path('results/')
if not resultsfolder.exists():
    os.mkdir(resultsfolder)

# functions
def graticule(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, R, uk, vk, u0, proj):
    """
    Creates graticule consisting of meridians and parallels based on input parameters.
    Returns lists of X and Y coordinates of meridians and parallels (4 lists).
    """

    # initiate empty lists for meridians
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

    # initiate empty lists for paralells
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
    """
    Transforms (u, v) to (s, d).
    """

    # Transformed latitude
    dv = vk - v
    s = np.arcsin(np.sin(u) * np.sin(uk) + np.cos(u) * np.cos(uk) * np.cos(dv))

    # Transformed longitude
    d = -1 * np.arctan2(np.cos(u) * np.sin(dv), np.cos(u) * np.sin(uk) * np.cos(dv) - np.sin(u) * np.cos(uk))

    return s, d


def gnom(R, s, d):
    """
    Gnomonic projection. Returns x, y coordinates from input s, d.
    """
    x = R * np.tan(np.pi/2 - s) * np.cos(d)
    y = R * np.tan(np.pi/2 - s) * np.sin(d)
    return x, y


def boundary(u, v, R, uk, vk, u0, proj):
    """
    Creates boundary of globe face based on input parameters.
    Returns list of X-coordinates and list of Y-coordinates for the boundary points.
    """
    XB, YB = [], []
    for u, v in zip(u, v):
        s, d = uv_to_sd(u, v, uk, vk)
        dv = vk - v
        XB_, YB_ = proj(R, s, d)

        XB.append(XB_)
        YB.append(YB_)

    return XB, YB


def rotate(X, Y, angle):
    """
    Rotates input points by rotation matrix (90 degrees).
    """
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

# set input parameters
us1 = 52.6226
uj1 = -us1
us2 = 10.8123
uj2 = -us2
ualpha = 26.5651

u_min = [-90, 30, -60, -60, -60, -60, -60, -20, -20, -20, -20, -20]
u_max = [-30, 90, 20,   20, 20,   20, 20,   60, 60, 60, 60, 60]
v_min = [-180, -180, -10, 60, 130, 200, -80, -50, 20, 90, 170, 240]
v_max = [180, 180, 80, 150, 230, 300, 10, 50, 120, 180, 260, 330]
 
u_k = [-90, 90, (uj1+us2)/2, (uj1+us2)/2, (uj1+us2)/2, (uj1+us2)/2, (uj1+us2)/2,
      (us1-uj2)/2, (us1-uj2)/2, (us1-uj2)/2, (us1-uj2)/2, (us1-uj2)/2]
v_k = [0, 0, 36, (144-72)/2+72, 144+36, 216+36, 288+36, 0, 72, 144, 216, 288]

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

num_faces = 12
# loop through faces with given input parameters
for i in range(num_faces):
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

    # create boundary
    XB, YB = boundary(ub, vb, R, uk, vk, u0, proj)

    # create graticule
    XM, YM, XP, YP = graticule(umin, umax, vmin, vmax, Du, Dv, du, dv, R, uk, vk, u0, proj)
    
    # rotate boundary and graticule by -90 degrees (for visualizing in ArcGIS)
    rotated_XB, rotated_YB = rotate(XB, YB, -90)
    rotated_XM, rotated_YM = rotate(XM, YM, -90)
    rotated_XP, rotated_YP = rotate(XP, YP, -90)

    # create polyline from rotated boundary points
    boundary_pts = arcpy.Array([arcpy.Point(x, y) for x, y in zip(rotated_XB, rotated_YB)])
    
    # save the rotated boundary as a shapefile to temporary folder
    output = arcpy.Polyline(boundary_pts)
    arcpy.CopyFeatures_management(output, f'results/boundary_face{i+1}.shp')

    # convert meridians from lists to array
    meridians = arcpy.Array()
    for x in range(len(rotated_XM)):
        meridians.add(arcpy.Point(rotated_XM[x], rotated_YM[x]))

    # save meridians as polyline to temporary folder
    output = arcpy.Polyline(meridians)
    arcpy.CopyFeatures_management(output, f'temp/meridians_face{i+1}_notsplit.shp')

    # delete longest lines from meridians (lines connecting end point of M1 and starting point of M2)
    arcpy.management.SplitLine(f'temp/meridians_face{i+1}_notsplit.shp', f'results/meridians_face{i+1}.shp')
    with arcpy.da.UpdateCursor(f"results/meridians_face{i+1}.shp", ["FID", "SHAPE@LENGTH"]) as cursor:
        for row in cursor:
            if row[1] > 1000000:
                cursor.deleteRow()

    # convert parallels from lists to array
    parallels = arcpy.Array()
    for y in range(len(rotated_XP)):
        parallels.add(arcpy.Point(rotated_XP[y], rotated_YP[y]))

    # save parallels as polyline to temporary folder
    output = arcpy.Polyline(parallels)
    arcpy.CopyFeatures_management(output, f'temp/parallels_face{i+1}_notsplit.shp')

    # delete longest lines from parallels (lines connecting end point of P1 and starting point of P2)
    arcpy.management.SplitLine(f'temp/parallels_face{i+1}_notsplit.shp', f'results/parallels_face{i+1}.shp')
    with arcpy.da.UpdateCursor(f"results/parallels_face{i+1}.shp", ["FID", "SHAPE@LENGTH"]) as cursor:
        for row in cursor:
            if row[1] > 1000000:
                cursor.deleteRow()

    print(f'face{i+1} done!')

print('done')
