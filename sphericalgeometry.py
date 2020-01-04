# -*- coding: utf-8 -*-
""" A module of convex hull on the unit sphere.
    A work in progress.

By Samuel Londner

TODO:
    * For module TODOs
    * Add functions

"""
import numpy as np
from scipy.spatial import ConvexHull

def sphere_convhull(lats, lons):
    """Computes the spherical convex hull of a set of points on the unit sphere.
    For now works only on points on North hemisphere

    Args:
        lats (double): latitudes of the points in degrees
        lons (double): longitudes of the points in degrees

    Returns:
        Indexes of the points forming the convex hull.

    References:
        - Frank Weller, Carsten Kirstein, "Computing the Convex Hull of a Simple Polygon on the Sphere" (1996) [https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.43.987]
        -  Lin-Lin Chen, T. C. Woo, "Computational Geometry on the Sphere With Application to Automated Machining" (1992)
        - C. Grima and A. Marquez, "Computational Geometry on Surfaces: Performing Computational Geometry on the Cylinder, the Sphere, the Torus, and the Cone", Springer, 2002. 

    .. _PEP 484:
        https://www.python.org/dev/peps/pep-0484/

    """
    #TODO: validate input
    #TODO: improve doc, correct doc (from "array" to list/NumPy array), specify lat/lon convention
    #TODO: change output format to be consistent with base module?
    #TODO: check hemispherity and cope with it
    #TODO: add support for south hemisphere
    #TODO: add support for any hemisphere

    phis = np.deg2rad(lats)
    lambdas = np.deg2rad(lons)
    
    #define central point of gnomonic projection
    assert np.all(phis>0) #check all poits on the open north hemisphere
    phi0 = np.deg2rad(90)
    lambda0 = np.deg2rad(0)

    #gnomonic projection - see http://mathworld.wolfram.com/GnomonicProjection.html
    cosc = np.sin(phi0) * np.sin(phis) + np.cos(phi0) * np.cos(phis) * np.cos(lambdas-lambda0)
    x = (np.cos(phis) * np.sin(lambdas-lambda0)) / cosc
    y = (np.cos(phi0)*np.sin(phis) - np.sin(phi0)*np.cos(phis)*np.cos(lambdas-lambda0) ) / cosc
    pts = np.column_stack((x,y))

    #compute planar convex hull
    # planar_set = shapely.geometry.MultiPoint(pts)
    # planar_convhull = planar_set.convex_hull
    # chull_x,chull_y = planar_convhull.exterior.xy
    planar_convhull_idx = ConvexHull(pts).vertices

    #return indexes
    return planar_convhull_idx
        
def geodesic_waypoints(lat0:float, lon0:float, lat1:float, lon1:float, pts_number=20):
    """Computes waypoints of the geodesic between points on the unit sphere.
    Assumes points are not antipodal. Uses SLERP interpolation

    Args:
        Geographic latitudes and longitudes of points in degrees.
        pts_number - the number of waypoints required.

    Returns:
        X, Y, Z of the waypoints

    References:
        - https://en.wikipedia.org/wiki/Slerp

    .. _PEP 484:
        https://www.python.org/dev/peps/pep-0484/
    """
#TODO: return lat/lons instead of XYZ

    #compute cartesian coordinates
    x0, y0, z0 = geog2cart(lat0, lon0)
    x1, y1, z1 = geog2cart(lat1, lon1)

    p0 = [x0, y0, z0]
    p1 = [x1, y1, z1]

    omega = np.arccos(np.dot(p0,p1))

    t = np.linspace(0, 1, pts_number)

    x_slerp = x0 * np.sin((1-t)*omega) / np.sin(omega) + x1 * np.sin(t*omega)/np.sin(omega)
    y_slerp = y0 * np.sin((1-t)*omega) / np.sin(omega) + y1 * np.sin(t*omega)/np.sin(omega)
    z_slerp = z0 * np.sin((1-t)*omega) / np.sin(omega) + z1 * np.sin(t*omega)/np.sin(omega)

    return x_slerp, y_slerp, z_slerp




def geog2cart(lats, lons, R=1):
    """Converts geographic coordinates on sphere to 3D Cartesian points.

    Args:
       Latitude(s) and longitude(s) in degrees

    Returns:
        X,Y,Z



    .. _PEP 484:
        https://www.python.org/dev/peps/pep-0484/
    """
    lonr = np.deg2rad(lons)
    latr = np.deg2rad(lats)
    pts_x = R * np.cos(latr)*np.cos(lonr)
    pts_y = R * np.cos(latr)*np.sin(lonr)
    pts_z = R * np.sin(latr)

    return pts_x, pts_y, pts_z
    