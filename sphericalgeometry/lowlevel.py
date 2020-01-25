import numpy as np

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
    y = (np.cos(phis) * np.sin(lambdas-lambda0)) / cosc
    x = (np.cos(phi0)*np.sin(phis) - np.sin(phi0) *
         np.cos(phis)*np.cos(lambdas-lambda0)) / cosc
    #x and y are chosen so that they fit the 3D cartesian coordinates convention
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