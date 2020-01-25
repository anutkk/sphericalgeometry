# -*- coding: utf-8 -*-
""" Spherical geometry functions.

By Samuel Londner


"""
import numpy as np
from scipy.spatial import ConvexHull
from .lowlevel import *

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
    # TODO: add support for points on equator
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
    # intersecting hulls of the two sets of points (see in reference)
    else:
        #project points in each hemisphere separately
        north_lats = lats[kk_north]
        north_lons = lons[kk_north]
        south_lats = lats[kk_south]
        south_lons = lons[kk_south]
        north_xypts = gnomonic_proj(lat0, lon0, north_lats, north_lons)
        south_xypts = gnomonic_proj(lat0, lon0, south_lats, south_lons)
        #compute planar convex hulls
        north_cvhull_idx = np.array(ConvexHull(north_xypts).vertices)
        south_cvhull_idx = np.array(ConvexHull(south_xypts).vertices)
        north_cvhull = north_xypts[north_cvhull_idx,:]
        south_cvhull = south_xypts[south_cvhull_idx,:]

        #check if the convex hulls intersect by passing over every edge and check whether it forms a separating axis
        #TODO: currently the method used is O(n^2), based on https://stackoverflow.com/questions/753140/how-do-i-determine-if-two-convex-polygons-intersect
        #consider logarithmic time algorithm using https://epubs.siam.org/doi/pdf/10.1137/1.9781611973730.109
        #TODO: values are computed multiple times, optimize
        #TODO: code is duplicated, clean up
        jj = 0
        intersecting = True

        while jj<len(north_cvhull): #first iterate over edges of "northern" convex hull
            jj1 = (jj+1) % len(north_cvhull)
            # compute linear equation describing the edge
            m = (north_cvhull[jj1,1]-north_cvhull[jj,1])/(north_cvhull[jj1,0]-north_cvhull[jj,0])
            b = north_cvhull[jj,1] - m*north_cvhull[jj,0]
            # check which side is every vertex of the hulls
            north_points_signs = np.sign(north_cvhull[:,1] - m*north_cvhull[:,0] - b)
            south_points_signs = np.sign(south_cvhull[:,1] - m*south_cvhull[:,0] - b)

            if (np.all(north_points_signs>=0) and np.all(south_points_signs<0)) or (np.all(north_points_signs<=0) and np.all(south_points_signs>0)): #if the edge separates the points...
                intersecting = False #mark as non-intersecting
                # move slightly the separating line so that it does not touch any of the hulls
                distances_from_line = np.abs(south_cvhull[:,1] - m*south_cvhull[:,0] - b)
                min_dist_idx = np.argmin(distances_from_line)
                min_dist = distances_from_line[min_dist_idx]
                b_delta = np.sign(south_cvhull[min_dist_idx,1] - m*south_cvhull[min_dist_idx,0] - b) * (min_dist) * 0.5
                xs = np.array([np.min(south_cvhull[:,0]), np.max(north_cvhull[:,0])])
                separating_cart = np.column_stack((xs, m*xs+b+b_delta,[1,1]))
                break
            jj+=1
        if intersecting: #if there was no success looking for in the northern polygon edges, search in the southern polygon edges
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
        
        if intersecting: #no separating edge was found
            raise ValueError('Points are not hemispheric, the convex hull is the whole sphere.')
        
        #else the points are on the same hemisphere, so find the appropriate point for central projection
        else:
            #find normal vector of plane of great circle
            normal_vector = np.cross(separating_cart[1,:], separating_cart[0,:])
            normal_vector = normal_vector / np.linalg.norm(normal_vector)
            lat0, lon0 = cart2geog(normal_vector[0], normal_vector[1], normal_vector[2])
            
            # project
            pts = gnomonic_proj(lat0, lon0, lats, lons)
            # compute planar convex hull
            planar_convhull_idx = ConvexHull(pts).vertices
            # return indexes
            return planar_convhull_idx

            #TODO: use already computed convex hulls!