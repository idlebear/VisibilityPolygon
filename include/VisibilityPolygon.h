//
// Created by bjgilhul on 19/06/19.
//

#ifndef OB1_VISIBILITYPOLYGON_H
#define OB1_VISIBILITYPOLYGON_H

#include <limits>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <memory>
#include <stdexcept>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/core/coordinate_system.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

#include <boost/geometry/algorithms/intersection.hpp>


using namespace std;
namespace bg = boost::geometry;

namespace Visibility {

    const double VISIBILITY_POLYGON_EPSILON = 0.0000001;

    typedef bg::model::d2::point_xy<double, bg::cs::cartesian> Point;
    typedef bg::model::polygon<Point, false, true> Polygon;
    typedef bg::model::segment<Point> Segment;

    inline bool
    inViewport(const Point &position, const Point &viewportMinCorner, const Point &viewportMaxCorner) {
        if (position.x() < viewportMinCorner.x() - VISIBILITY_POLYGON_EPSILON ||
            position.y() < viewportMinCorner.y() - VISIBILITY_POLYGON_EPSILON ||
            position.x() > viewportMaxCorner.x() + VISIBILITY_POLYGON_EPSILON ||
            position.y() > viewportMaxCorner.y() + VISIBILITY_POLYGON_EPSILON) {
            return false;
        }
        return true;
    }

    inline bool
    operator==( const Point& lhs, const Point& rhs ) {
        return bg::distance( lhs, rhs ) < VISIBILITY_POLYGON_EPSILON;
    }

    inline bool
    operator!=( const Point& lhs, const Point& rhs ) {
        return !( lhs == rhs );
    }

    vector<Point>
    convertToPoints(const Segment &segment);

    inline bool
    doSegmentsIntersect( const Segment& s1, const Segment& s2 ) {
        return bg::intersects( s1, s2 );
    }

    inline bool
    intersectSegments( const Segment& s1, const Segment& s2, vector<Point>& res ) {
        bg::intersection( s1, s2, res );
        return !res.empty();
    }

    inline bool
    intersectLines( const Segment& a, const Segment& b, Point& res ) {
        auto aPts = convertToPoints(a);
        auto bPts = convertToPoints(b);
        auto dax = aPts[1].x() - aPts[0].x();
        auto day = aPts[1].y() - aPts[0].y();
        auto dbx = bPts[1].x() - bPts[0].x();
        auto dby = bPts[1].y() - bPts[0].y();

        auto u_b = dby * dax - dbx * day;
        if( u_b == 0) {
            return false;
        }
        auto dabx = aPts[1].x() - bPts[1].x();
        auto daby = aPts[1].y() - bPts[1].y();
        auto ua = (dby * dabx - dbx * daby)/u_b;
        dax *= ua;
        day *= ua;
        res = Point( aPts[1].x() + dax, aPts[1].y() + day );
        return true;
    };

    struct PointIndex {
        int i;
        int j;
        double a;
    };

    vector <PointIndex> sortPoints(const Point &position, const vector <Segment> &segments);

    inline double
    angle(const Point &a, const Point &b) {
      return atan2(b.y() - a.y(), b.x() - a.x());
    };

    inline double
    angle2(const Point &a, const Point &b, const Point &c) {
      auto res = angle(a, b) - angle(b, c);
      if (res < 0) {
        res += 2 * M_PI;
      } else if (res > 2 * M_PI) {
        res -= 2 * M_PI;
      }
      return res;
    };

    inline bool
    contains( const Polygon& poly, const Point& pt ) {
        return bg::within( pt, poly );
    }

    vector<Segment>
    convertToSegments( const vector<Polygon>& polygons );

    inline void
    addPoint( Polygon& poly, const Point& pt ) {
        bg::append( poly, pt );
    }

    Polygon
    compute(const Point &position, const vector<Segment> &segments);

    Polygon
    computeViewport(const Point &position, const vector<Segment> &segments, const Point &viewportMinCorner,
                    const Point &viewportMaxCorner);

    vector<Segment>
    convertToSegments(const Polygon &poly);

    vector<Segment>
    convertToSegments(const vector<Polygon> &poly);

    vector<Point>
    convertToPoints(const Polygon &poly);

    vector<Segment>
    breakIntersections(const vector<Segment> &segments);

    inline
    int heapParent(int index) {
      return floor((index - 1) / 2);
    };

    inline
    int heapChild(int index) {
      return 2 * index + 1;
    };

    void remove(int index, vector<int>& heap, const Point &position, const vector<Segment> &segments,
                const Point &destination, vector<int>& map);

    void
    insert(int index, vector<int>& heap, const Point &position, const vector<Segment> &segments, const Point &destination,
           vector<int>& map);

    bool lessThan(int i1, int i2, const Point &position, const vector<Segment> &segments, const Point &destination);

};



#endif //OB1_VISIBILITYPOLYGON_H
