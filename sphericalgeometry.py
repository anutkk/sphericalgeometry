# -*- coding: utf-8 -*-
""" Spherical geometry functions.

By Samuel Londner

TODO:
    * For module TODOs

"""
import numpy as np
from scipy.spatial import ConvexHull


def sphere_convhull(lats, lons):
    """Computes the spherical convex hull of a set of points on the unit sphere.
    For now works only when no point is one the Equator (lat=0).

    Args:
        lats (double): geographic latitudes of the points in degrees
        lons (double): geographic longitudes of the points in degrees

    Returns:
        Indexes of the points forming the convex hull.

    References:
        -  Lin-Lin Chen, T. C. Woo, "Computational Geometry on the Sphere With Application to Automated Machining" (1992)
        - O'Rourke  J., Chien  C., Olson T. and Naddor, D.,"A New Linear Algorithm for Intersecting Convex Polygons" (1982)

    """
    # TODO: improve doc, correct doc (from "array" to list/NumPy array)
    # TODO: change output format to be consistent with base module?
    # TODO: check hemispherity and cope with it
    # TODO: add support for any hemisphere
    # TODO: maybe use Grima-Marquez version of Graham's algorithm?

    # validate input
    assert len(lats)==len(lons), "Latitudes and longitudes must have the same length."
    assert len(lats)>2, "There must be at least 3 points."
    assert np.all(lats >= -90) and np.all(lats <=90) ,  "Latitudes should be in the range [-90, 90]"
    assert np.all(lons >= -180) and np.all(lats <= 180) , "Longitudes should be in the range [-180, 180]"

    # define central point of gnomonic projection
    # assert np.all(lats > 0)  # check all points on the open north hemisphere
    lat0 = 90
    lon0 = 0

    # ensure no points are on the equator
    if np.any(lats == 0):
        raise NotImplementedError('Points on the Equator are currently not supported.')

    # check points are on the same hemisphere by intersecting their convex hulls
    kk_north = lats > 0
    kk_south = lats < 0

    # if all points are in the Northen or Southern hemisphere the solution is simple
    if np.count_nonzero(kk_north) == 0 or np.count_nonzero(kk_south) == 0: 
        # gnomonic projection
        pts = gnomonic_proj(lat0, lon0, lats, lons)
        # compute planar convex hull
        planar_convhull_idx = ConvexHull(pts).vertices
        # return indexes
        return planar_convhull_idx

    # no clear hemisphere - check if there is another hemisphere by 
    # intersecting hulls of the two sets of points
    else:
        raise NotImplementedError('Points accros the Equator are currently not supported.')
        north_lats = lats[lats > 0]
        north_lons = lons[lats > 0]
        south_lats = lats[lats < 0]
        south_lons = lons[lats < 0]
        north_xypts = gnomonic_proj(lat0, lon0, north_lats, north_lons)
        south_xypts = gnomonic_proj(lat0, lon0, south_lats, south_lons)
        # north_cvhull = np.array(ConvexHull(north_xypts).points)
        # south_cvhull = np.array(ConvexHull(south_xypts).points)
        north_cvhull_idx = np.array(ConvexHull(north_xypts).vertices)
        south_cvhull_idx = np.array(ConvexHull(south_xypts).vertices)
        north_cvhull = north_xypts[north_cvhull_idx,:]
        south_cvhull = south_xypts[south_cvhull_idx,:]

        #check if they intersect by passing over every edge and check whether it forms a separating axis
        #currently the method used is O(n^2), based on https://stackoverflow.com/questions/753140/how-do-i-determine-if-two-convex-polygons-intersect
        #TODO: consider logarithmic time algorithm using https://epubs.siam.org/doi/pdf/10.1137/1.9781611973730.109
        #TODO: values are computed multiple times, optimize
        #TODO: code is duplicated, clean up
        jj = 0
        intersecting = True
        while jj<len(north_cvhull):
            jj1 = (jj+1) % len(north_cvhull)
            m = (north_cvhull[jj1,1]-north_cvhull[jj,1])/(north_cvhull[jj1,0]-north_cvhull[jj,0])
            b = north_cvhull[jj,1] - m*north_cvhull[jj,0]
            north_points_signs = np.sign(north_cvhull[:,1] - m*north_cvhull[:,0] - b)
            # northPointsSigns = np.delete(northPointsSigns, northPointsSigns==0)
            south_points_signs = np.sign(south_cvhull[:,1] - m*south_cvhull[:,0] - b)
            if (np.all(north_points_signs>=0) and np.all(south_points_signs<0)) or (np.all(north_points_signs<=0) and np.all(south_points_signs>0)):
                intersecting = False
                distances_from_line = np.abs(south_cvhull[:,1] - m*south_cvhull[:,0] - b) #/(m**2+1)
                min_dist_idx = np.argmin(distances_from_line)
                min_dist = distances_from_line[min_dist_idx]
                b_delta = np.sign(south_cvhull[min_dist_idx,1] - m*south_cvhull[min_dist_idx,0] - b) * (min_dist) * 0.5
                xs = np.array([np.min(south_cvhull[:,0]), np.max(north_cvhull[:,0])])
                separating_cart = np.column_stack((xs, m*xs+b+b_delta,[1,1]))

                # lats_separating = north_lats[north_cvhull_idx[[jj, jj1]]]
                # lons_separating = north_lons[north_cvhull_idx[[jj, jj1]]]
                break
            jj+=1
        if intersecting: #keep checking the edges of the other polygon
            jj=0
            while jj<len(south_cvhull):
                jj1 = (jj+1) % len(south_cvhull)
                m = (south_cvhull[jj1,1]-south_cvhull[jj,1])/(south_cvhull[jj1,0]-south_cvhull[jj,0])
                b = south_cvhull[jj,1] - m*south_cvhull[jj,0]
                north_points_signs = np.sign(north_cvhull[:,1] - m*north_cvhull[:, 0] - b)
                # northPointsSigns = np.delete(northPointsSigns, northPointsSigns==0)
                south_points_signs = np.sign(south_cvhull[:,1] - m*south_cvhull[:,0] - b)
                if (np.all(north_points_signs>0) and np.all(south_points_signs<=0)) or (np.all(north_points_signs<0) and np.all(south_points_signs>=0)):
                    intersecting = False
                    distances_from_line = np.abs(north_cvhull[:,1] - m*north_cvhull[:,0] - b) #/(m**2+1)
                    min_dist_idx = np.argmin(distances_from_line)
                    min_dist = distances_from_line[min_dist_idx]
                    b_delta = np.sign(north_cvhull[min_dist_idx,1] - m*north_cvhull[min_dist_idx,0] - b) * (min_dist) * 0.5
                    separating_cart = np.array([
                        [south_cvhull[jj,0], south_cvhull[jj,1]+b_delta, 1 ],
                        [south_cvhull[jj1,0], south_cvhull[jj1,1]+b_delta, 1 ]
                    ])
                    break
                jj += 1
        
        if intersecting:
            raise ValueError('Points are not hemispheric, the convex hull is the whole sphere.')
        
        #the points are on the same hemisphere, so find the appropriate point for central projection
        else:
            #find normal vector of plane of great circle (sign does not matter, it will work either way)
            # separating_cart = np.column_stack(geog2cart(lats_separating, lons_separating))
            normal_vector = np.cross(separating_cart[1,:], separating_cart[0,:])
            normal_vector = normal_vector / np.linalg.norm(normal_vector)
            lat0, lon0 = cart2geog(normal_vector[0], normal_vector[1], normal_vector[2])
            
            # gnomonic projection
            pts = gnomonic_proj(lat0, lon0, lats, lons)
            # compute planar convex hull
            planar_convhull_idx = ConvexHull(pts).vertices
            # return indexes
            return planar_convhull_idx

            #TODO: use already computed convex hulls!


def gnomonic_proj(lat0: float, lon0: float, lats, lons):
    """Computes gnomonic projection (also called central projections) of points
    on the sphere.

    Args:
        lat0/lon0 - reference/tangent point of projection
        lats/lons - Geographic latitudes and longitudes of points in degrees to be projected.

    Returns:
        2D array of the projected points on the projective plane.

    References:
        - http://mathworld.wolfram.com/GnomonicProjection.html
    """
    phi0 = np.deg2rad(lat0)
    lambda0 = np.deg2rad(lon0)
    phis = np.deg2rad(lats)
    lambdas = np.deg2rad(lons)

    cosc = np.sin(phi0) * np.sin(phis) + np.cos(phi0) * \
        np.cos(phis) * np.cos(lambdas-lambda0)
    x = (np.cos(phis) * np.sin(lambdas-lambda0)) / cosc
    y = (np.cos(phi0)*np.sin(phis) - np.sin(phi0) *
         np.cos(phis)*np.cos(lambdas-lambda0)) / cosc
    pts = np.column_stack((x, y))

    return pts


def geodesic_waypoints(lat0: float, lon0: float, lat1: float, lon1: float, pts_number=20):
    """Computes waypoints of the geodesic between points on the unit sphere.
    Assumes points are not antipodal. Uses SLERP interpolation

    Args:
        Geographic latitudes and longitudes of points in degrees.
        pts_number - the number of waypoints required.

    Returns:
        X, Y, Z of the waypoints.

    References:
        - https://en.wikipedia.org/wiki/Slerp
    """
    # TODO: return lat/lons instead of XYZ ?

    # compute cartesian coordinates
    x0, y0, z0 = geog2cart(lat0, lon0)
    x1, y1, z1 = geog2cart(lat1, lon1)

    p0 = [x0, y0, z0]
    p1 = [x1, y1, z1]

    omega = np.arccos(np.dot(p0, p1))

    t = np.linspace(0, 1, pts_number)

    x_slerp = x0 * np.sin((1-t)*omega) / np.sin(omega) + \
        x1 * np.sin(t*omega)/np.sin(omega)
    y_slerp = y0 * np.sin((1-t)*omega) / np.sin(omega) + \
        y1 * np.sin(t*omega)/np.sin(omega)
    z_slerp = z0 * np.sin((1-t)*omega) / np.sin(omega) + \
        z1 * np.sin(t*omega)/np.sin(omega)

    return x_slerp, y_slerp, z_slerp


def geog2cart(lats, lons, R=1):
    """Converts geographic coordinates on sphere to 3D Cartesian points.

    Args:
       Latitude(s) and longitude(s) in degrees.

    Returns:
        X,Y,Z
    """
    lonr = np.deg2rad(lons)
    latr = np.deg2rad(lats)
    pts_x = R * np.cos(latr)*np.cos(lonr)
    pts_y = R * np.cos(latr)*np.sin(lonr)
    pts_z = R * np.sin(latr)

    return pts_x, pts_y, pts_z

def cart2geog(x,y,z, R=1):
    """Converts 3D Cartesian points to geographic coordinates on sphere.

    Args:
       3D coordinates of points..

    Returns:
        Latitudes and longitudes of points in degrees.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    # assert np.all(np.abs(r-R) < np.finfo(float).eps) , "Points are not on sphere!"
    assert np.all(np.abs(r-R) < R/100) , "Points are not on sphere!"

    lons = np.arctan2(y, x)
    lats = np.arcsin(z/R)

    return np.rad2deg(lats), np.rad2deg(lons)