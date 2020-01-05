import random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from itertools import product, combinations

import sphericalgeometry

# Create array of random points on the sphere
# See https://mathworld.wolfram.com/SpherePointPicking.html
# Part of the code has been inspired by https://github.com/VictorDavis/GeoConvexHull/blob/master/GeoConvexHull.py by Victor Davis
lons = np.array([])
lats = np.array([])
number_of_points = 20

for jj in range(0, number_of_points):
    u = random.random()
    v = random.random()

    # convert to spherical coordinates
    theta = 2 * np.pi * u # range [0,2pi)
    phi = np.arccos(2 * v - 1) * 0.5 # range (0, pi/2] - points limited to Northern hemisphere

    # convert to lng,lat
    lng = theta / (2 * np.pi) * 360 - 180
    lat =   phi / (2 * np.pi) * 360 

    lons = np.append(lons, lng)
    lats = np.append(lats, lat)

# Compute convex hull
chull_idx = sphericalgeometry.sphere_convhull(lats, lons)

# Display results
fig = plt.figure()
ax = fig.gca(projection='3d')
# ax.set_aspect("equal")

earthRadius =1 # 6378.137;

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = earthRadius * np.cos(u)*np.sin(v)
y = earthRadius * np.sin(u)*np.sin(v)
z = earthRadius * np.cos(v)
ax.plot_wireframe(x, y, z, color="k")
#TODO: use geog2cart function
#The formulaes used do not fit the usual geographic lat/lon convention and instead are for colatitude. Since this is a sphere mesh it does not matter.

pts_x, pts_y, pts_z = sphericalgeometry.geog2cart(lats, lons)

ax.scatter(pts_x, pts_y, pts_z, color="g", s=100)

#close polygon
chull_lons = np.append(lons[chull_idx], lons[chull_idx[0]])
chull_lats = np.append(lats[chull_idx], lats[chull_idx[0]])
# chull_pts_x, chull_pts_y, chull_pts_z = sphericalgeometry.geog2cart(chull_lats, chull_lons)
# ax.plot(chull_pts_x, chull_pts_y, chull_pts_z, color="r")

#interpolate so that spherical polygon is correctly displayed
chull_pts_x = np.array([])
chull_pts_y = np.array([])
chull_pts_z = np.array([])

for ii in range(len(chull_lons)-1):
    print(ii)
    lat0, lon0 = chull_lats[ii], chull_lons[ii]
    lat1, lon1 = chull_lats[ii+1], chull_lons[ii+1]
    x_slerp, y_slerp, z_slerp=sphericalgeometry.geodesic_waypoints(lat0, lon0, lat1, lon1)
    chull_pts_x = np.append(chull_pts_x, x_slerp)
    chull_pts_y = np.append(chull_pts_y, y_slerp)
    chull_pts_z = np.append(chull_pts_z, z_slerp)


ax.plot(chull_pts_x, chull_pts_y, chull_pts_z, color='r')

plt.show()