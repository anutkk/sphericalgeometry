# sphericalgeometry
A Python package for spherical computational geometry.

This package contains various spherical geometry functions that I needed and could not find in well-known libraries.
The functions work on sets of points on the sphere.

The package's API is focused on simplicity, casuality and easiness of use. 
Inputs are arrays of geographic coordinates and outputs are scalars or arrays of geographic coordinates.
Specifically, the package is **not** object-oriented.

__This package is under active development. Expect API changes.__

## High-level functions

* `sphere_convhull`: Computes the spherical convex hull of a set of points on the unit sphere.
* `polygon_area`: Computes the area of a polygon on the sphere.

<!--

* `minareapolygon`: Computes the minimum-area bounding polygon of a set of points on the sphere.

* `pspolydist` - Computes the minimum distance from a point to a polyline on a sphere. 
This is an adaptation to Python of the original MATLAB function [`p_spoly_dist`](https://www.mathworks.com/matlabcentral/fileexchange/52734-p_spoly_dist) by Michael Yoshpe (in-work).

* `insphericalpolygon` - Returns points inside of a spherical polygon.
-->

## Low-level functions

Those functions were included mainly to be used by high-level functions or example scripts.
Their syntax may change in the future.

* `gnomonic_proj` - Computes the gnomonic projection (also called central projection) of points on the sphere.

* `geog2cart` - Converts geographic coordinates (lat/lon) to 3D Cartesian coordinates.

* `cart2geog` - Convert 3D Cartesian coordinates to geographic coordinates.

* `geodesic_waypoints` - Computes waypoints along a geodesic ("great circle") between two points.

## Coordinates system

The package was designed so that the user needs to deal exclusively with geographic coordinates (latitude and longitude).
The units are degrees.
However, some of the low-level functions return and use 3D Cartesian coordinates (X, Y, Z).

<!--
## The inside of a spherical polygon


-->

## Requirements

* Python >= 3.6
* Numpy
* Scipy

## Limitations

There is no support for disjoint polygons.

## Alternatives

Check out the [spherical_geometry](https://github.com/spacetelescope/spherical_geometry) package by spacetelescope. This is a stable and strongly object-oriented package, but I am unsure if it is still maintained.

## Some basic notions of spherical geometry
### Great circles and geodesics

In "regular" Euclidean geometry, the shortest path between two points on the plane is simply the straight segment between them. But when the points are on a sphere, the shortest path on the sphere is not a straight 3D segment - since this segment crosses the sphere but is not part of it. Instead, this would be a "line on the sphere". The shortest between two points over a surface is sometimes called a geodesic <sup>[1](#myfootnote1)</sup>.

Given two "close" points on the sphere, they define a unique circle on the sphere which passes through them. We call this circle a "great circle". The great circle is separated in two paths, both of them connecting the two points. The shortest path on the sphere between the points is the smallest part of the corresponding great circle.

<!-- TODO: add figure -->

Now consider the two poles of the sphere (the North pole and the South Pole, by analogy to the Earth, which is not a sphere). What is the shortest path between them? In fact, since the sphere is symmetric around its center, there are an infinity of shortest paths. These points are called _antipodal_. 

### Convexity on the plane and on the sphere

> A planar polygon is convex if it contains all the line segments connecting any pair of its points.
> -- [Wolfram MathWorld](http://mathworld.wolfram.com/ConvexPolygon.html)

As such, the convex hull of a set of points is the smallest convex polygon that contains the set. This definition is simple and intuitive in the (Euclidean) plane, since the "line segment connecting two points" is well-defined and always unique<sup>[2](#myfootnote2)</sup>.

However, if we adapt the definition to a set of points on a sphere, the convex hull problem may be ill-defined if the set contains antipodal points. In this case, the shortests paths connecting the points cover the whole sphere, and the convex hull, which must contain all shortest paths is the whole sphere. More generally, this would happen if the set of points are not contained within an open hemisphere. This makes the definition of "convex hull" not very useful, since for many simple cases it would be the entire surface.

In order to define the spherical convex hull problem without encountering this problem, two approachs exist:

* Spherical convexity is defined for any set of points. A set of points on the sphere is said to be spherically convex, if for any two points, the geodesic joining them lies entirely in the set.  If the points are not contained within an open hemisphere, the convex hull is defined to be the whole sphere [Grima and Marquez, p. 47].
* We define a new concept, _strict_ convexity, which is defined only for points contained within an open hemisphere. If the points are not contained within an open hemisphere, there is no strictly convex hull, and the non-strictly convex hull is the whole sphere [Chen and Woo].

The implementation of the (strictly) convex hull in this library follows the second approach. If the points are not contained within an open hemisphere, the function throws an exception.

### Is Convex Hull always the solution? An example

Consider the following configuration of points (from [Math Overflow](https://mathoverflow.net/questions/76875/convex-hull-on-a-riemannian-manifold)), which do not fit in a hemisphere.

![](images/convex_hull_degenerate_case.jpg?raw=true)

 Suppose the sphere is Earth, and think of the points as locations of hurricanes during last year. We would like to find the area on Earth which is most prone to hurricanes. Intuitively, the points can be enclosed by a small shape - the blue polygon. However, the convex hull is the whole sphere, since the points do not fit in a hemisphere.

**In other words, the spherical convex hull fails to represent the actual geometry of physical problems**. This makes it less useful in real-world applications.

<!--
### Minimum Enclosing Polygon

Minimum-perimeter bounding spherical polygon

-->

## References:
1. C. Grima and A. Marquez, ["Computational Geometry on Surfaces: Performing Computational Geometry on the Cylinder, the Sphere, the Torus, and the Cone"](https://www.springer.com/gp/book/9781402002021), Springer (2002). 
1. Lin-Lin Chen, T. C. Woo, ["Computational Geometry on the Sphere With Application to Automated Machining"](https://asmedigitalcollection.asme.org/mechanicaldesign/article-abstract/114/2/288/431533/Computational-Geometry-on-the-Sphere-With) (1992).
1. Frank Weller, Carsten Kirstein, ["Computing the Convex Hull of a Simple Polygon on the Sphere"](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.43.987) (1996).
1. Chamberlain, Robert G.; Duquette, William H., ["Some algorithms for polygons on a sphere"](https://trs.jpl.nasa.gov/handle/2014/40409), NASA Jet Propulsion Laboratory (2007).
1. Ning Wang, ["An efficient search algorithm for minimum covering polygons on the sphere"](https://esrl.noaa.gov/gsd/nim/references/efficient_search_algorithm_for_minimum_covering_polygons_on_the_sphere.pdf), SIAM J. SCI. COMPUT (2013).

## Notes

<a name="myfootnote1">1</a>: This is the original meaning of the word "geodesic". In Riemann geometry the geodesic is the whole great circle - see on [Wikipedia](https://en.wikipedia.org/wiki/Geodesic#Metric_geometry). Here I use the word as its meaning in geodesy, in order to differentiate from the whole great circle.

<a name="myfootnote2">1</a>: Other definitions of convexity exist, which are equivalent in the planar case but not necessarily in the spherical case [Grima and Marquez, pp 31-36].
