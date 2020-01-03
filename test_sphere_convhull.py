import random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from itertools import product, combinations

import sphere_convhull

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
chull_idx = sphere_convhull.sphere_convhull(lats, lons)



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

lonr = np.deg2rad(lons)
latr = np.deg2rad(lats)
pts_x = earthRadius * np.cos(latr)*np.cos(lonr)
pts_y = earthRadius * np.cos(latr)*np.sin(lonr)
pts_z = earthRadius * np.sin(latr)

ax.scatter(pts_x, pts_y, pts_z, color="g", s=100)

chull_lonr = np.deg2rad(lons[chull_idx])
chull_latr = np.deg2rad(lats[chull_idx])
#close polygon
chull_lonr = np.append(chull_lonr, chull_lonr[1])
chull_latr = np.append(chull_latr, chull_latr[1])
chull_pts_x = earthRadius * np.cos(chull_latr)*np.cos(chull_lonr)
chull_pts_y = earthRadius * np.cos(chull_latr)*np.sin(chull_lonr)
chull_pts_z = earthRadius * np.sin(chull_latr)

#TODO: interpolate so that a correct polygon on the sphere is displayed

ax.plot(chull_pts_x, chull_pts_y, chull_pts_z, color="r")

plt.show()