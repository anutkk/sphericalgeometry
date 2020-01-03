import random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from itertools import product, combinations

#TODO

lons = np.array([])
lats = np.array([])
number_of_points = 20

for jj in range(0, number_of_points):
    u = random.random()
    v = random.random()

    # convert to spherical coordinates
    theta = 2 * np.pi * u # range [0,2pi)
    phi = np.arccos(2 * v - 1) * 0.5 # range (0, pi/2]

    # convert to lng,lat
    lng = theta / (2 * np.pi) * 360 - 180
    lat =   phi / (2 * np.pi) * 360 

    lons = np.append(lons, lng)
    lats = np.append(lats, lat)


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_wireframe(x, y, z, color="r")

earthRadius =1 # 6378.137;
lonr = np.deg2rad(lons)
latr = np.deg2rad(lats)
pts_x = earthRadius * np.cos(latr)*np.cos(lonr)
pts_y = earthRadius * np.cos(latr)*np.sin(lonr)
pts_z = earthRadius * np.sin(latr)

ax.scatter(pts_x, pts_y, pts_z, color="g", s=100)

plt.show()