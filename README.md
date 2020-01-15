# sphericalgeometry
A Python module of various spherical computational geometry related functions.

This module contains various spherical geometry functions that I needed and could not find in well-known libraries.
The functions work on sets of points on the unit sphere.

## High-level functions
Functions having high-level functionality:

* `sphere_convhull` - Computes the spherical convex hull of a set of points on the unit sphere.
<!--
* `pspolydist` - Computes the minimum distance from a point to a polyline on a sphere (unfinished). This is an adaptation to Python of the original MATLAB function [`p_spoly_dist`](https://www.mathworks.com/matlabcentral/fileexchange/52734-p_spoly_dist) by Michael Yoshpe (in-work).

* `insphericalpolygon` - Returns points inside of a spherical polygon (in-work).
-->
## Service functions
Basic non-original functions, which were included mainly to be used by high-level functions or scripts (some do not exist yet):

* `gnomonic_proj` - Computes the gnomonic projection (also called central projection) of points on the sphere.

* `geog2cart` - Converts geographic coordinates (lat/lon) to 3D Cartesian coordinates.

* `geodesic_waypoints` - Computes waypoints along a geodesic ("great circle") between two points.

* `cart2geog` - Convert 3D Cartesian coordinates to geographic coordinates.

## Some basic notions of spherical geometry
### Great circles and geodesics

In "regular" Euclidean geometry, the shortest path between two points on the plane is simply the straight segment between them. But when the points are on a sphere, the shortest path on the sphere is not a straight 3D segment - since this segment crosses the sphere but is not part of it. Instead, this would be a "line on the sphere". The shortest between two points over a surface is sometimes called a geodesic <sup>[1](#myfootnote1)</sup>.

Given two "close" points on the sphere, they define a unique circle on the sphere which passes through them. We call this circle a "great circle". The great circle is separated in two paths, both of them connecting the two points. The shortest path on the sphere between the points is the smallest part of the corresponding great circle.

<!-- TODO: add figure -->

Now consider the two poles of the sphere (the North pole and the South Pole, by analogy to the Earth, which is not a sphere). What is the shortest path between them? In fact, since the sphere is symmetric around its center, there are an infinity of shortest paths.

<!--
### Convexity on the plane and on the sphere

TODO

Convexity of antipodal points

Quasi-convexity

### Is Convex Hull always the solution?

TODO

See [this question](https://mathoverflow.net/questions/76875/convex-hull-on-a-riemannian-manifold) on Math Overflow.

### Minimum Enclosing Polygon

Minimum-perimeter bounding spherical polygon


TODOs
-->

## References:
1. C. Grima and A. Marquez, ["Computational Geometry on Surfaces: Performing Computational Geometry on the Cylinder, the Sphere, the Torus, and the Cone"](https://www.springer.com/gp/book/9781402002021), Springer (2002). 
1. Lin-Lin Chen, T. C. Woo, ["Computational Geometry on the Sphere With Application to Automated Machining"](https://asmedigitalcollection.asme.org/mechanicaldesign/article-abstract/114/2/288/431533/Computational-Geometry-on-the-Sphere-With) (1992).
1. Frank Weller, Carsten Kirstein, ["Computing the Convex Hull of a Simple Polygon on the Sphere"](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.43.987) (1996).
1. Chamberlain, Robert G.; Duquette, William H., ["Some algorithms for polygons on a sphere"](https://trs.jpl.nasa.gov/handle/2014/40409), NASA Jet Propulsion Laboratory (2007).
1. Ning Wang, ["An efficient search algorithm for minimum covering polygons on the sphere"](https://esrl.noaa.gov/gsd/nim/references/efficient_search_algorithm_for_minimum_covering_polygons_on_the_sphere.pdf), SIAM J. SCI. COMPUT (2013).

## Notes

<a name="myfootnote1">1</a>: This is the original meaning of the word "geodesic". In Riemann geometry the geodesic is the whole great circle - see on [Wikipedia](https://en.wikipedia.org/wiki/Geodesic#Metric_geometry). Here I use the word as its meaning in geodesy, in order to differentiate from the whole great circle.