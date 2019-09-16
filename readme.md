# Visibility

## based on a javascript library visibility_polygon.js version 1.9
### forked from [Github](https://github.com/byronknoll/visibility-polygon-js "Byron Knoll's Github repo")

This library can be used to construct a visibility polygon for a set of line segments.  Additional functions allow for intersection and union of polygons (requires Boost::Geometry).  Time permitting, cleanup is ongoing to convert some of the functions to templates to reduce the duplicated code. Welcome to Boost I guess...

The time complexity of the **compute()** implementation is O(N log N) (where N is the total number of line segments). This is the optimal time complexity for this problem.

The following functions should be useful:

1. Visibility::compute(position, segments)
    : Computes a visibility polygon. O(N log N) time complexity (where N is the number of line segments).
    * Arguments:
        :   position - The location of the observer. If the observer is not completely surrounded by line segments, an outer bounding-box will be automatically created (so that the visibility polygon does not extend to infinity).
        :   segments - A list of line segments. Each line segment should be a list of two points. Each point should be a list of two coordinates. Line segments can not intersect each other. Overlapping vertices are OK, but it is not OK if a vertex is touching the middle of a line segment. Use the "breakIntersections" function to fix intersecting line segments.
    * Returns: The visibility polygon (in clockwise vertex order).

2. Visibility::computeViewport(position, segments, viewportMinCorner, viewportMaxCorner)
    : Computes a visibility polygon within the given viewport. This can be faster than the "compute" function if there are many segments outside of the viewport.
    * Arguments:
        :   position - The location of the observer. Must be within the viewport.
        :   segments - A list of line segments. Line segments can not intersect each other. It is OK if line segments intersect the viewport.
        :   viewportMinCorner - The minimum X and Y coordinates of the viewport.
        :   viewportMaxCorner - The maximum X and Y coordinates of the viewport.
    * Returns: The visibility polygon within the viewport (in clockwise vertex order).

3) Visibility::inPolygon(position, polygon)
    * Calculates whether a point is within a polygon. O(N) time complexity (where N is the number of points in the polygon).
    * Arguments:
        :   position - The point to check: a list of two coordinates.
        :   polygon - The polygon to check: a list of points. The polygon can be specified in either clockwise or counterclockwise vertex order.
    * Returns: True if "position" is within the polygon.

4) Visibility::convertToSegments(polygons)
    : Converts the given polygons to list of line segments. O(N) time complexity (where N is the number of polygons).
    * Arguments: a list of polygons (in either clockwise or counterclockwise vertex order). Each polygon should be a list of points. Each point should be a list of two coordinates.
    * Returns: a list of line segments.

5) Visibility::breakIntersections(segments)
    : Breaks apart line segments so that none of them intersect. O(N^2) time complexity (where N is the number of line segments).
    * Arguments: a list of line segments. Each line segment should be a list of two points. Each point should be a list of two coordinates.
    * Returns: a list of line segments.

### Example code:

~~~~
#include <vector>
#include "VisibilityPolygon.h"

namespace V = Visibility;

vector<V::Polygon> polygons {
    {{
        {-1,-1}, 
        {501,-1},
        {501,501},
        {-1,501}
    }},
    {{
        {250,100},
        {260,140},
        {240,140}
    }},
};

auto segments = V::convertToSegments(polygons);
segments = V::breakIntersections(segments);
V::Point position(60, 60);
if( V::inPolygon(position, polygons[0]) ) {
  var visibility = V::compute(position, segments);
}
auto viewportVisibility = V::computeViewport(position, segments, {50, 50}, {450, 450});
~~~~
