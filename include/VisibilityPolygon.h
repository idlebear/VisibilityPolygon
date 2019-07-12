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
#include <boost/geometry/multi/geometries/multi_polygon.hpp>

#include <boost/geometry/algorithms/intersection.hpp>

#include <boost/math/special_functions/sign.hpp>
#include <boost/assign.hpp>

using namespace std;
namespace bg = boost::geometry;
// Bring "+=" for a vector into scope
using namespace boost::assign;

namespace Visibility {

    // define a margin of error (horseshoes and hand grenades...)
    const double VISIBILITY_POLYGON_EPSILON = 0.0000001;

    typedef bg::model::d2::point_xy<double, bg::cs::cartesian> Point;
    typedef bg::model::polygon<Point, false, true> Polygon;
    typedef bg::model::segment<Point> Segment;
    typedef bg::model::linestring<Point> PolyLine;
    typedef bg::model::multi_polygon<Polygon> MultiPolygon;
    typedef bg::model::ring<Point> Ring;
    typedef bg::model::box<Point> Box;

    /////////////////////////////////
    //
    // Point operations
    //
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


    // Ref: https://en.wikipedia.org/wiki/Quickhull
    Polygon
    quickhull( const vector<Point>& pts );


    /////////////////////////////////
    //
    // Segment operations
    //
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

    inline double
    length( const Segment& segment ) {
        return bg::distance( segment.first, segment.second );
    }

    inline bool
    intersectLines( const Segment& a, const Segment& b, Point& res ) {
        auto dax = a.second.x() - a.first.x();
        auto day = a.second.y() - a.first.y();
        auto dbx = b.second.x() - b.first.x();
        auto dby = b.second.y() - b.first.y();

        auto u_b = dby * dax - dbx * day;
        if( u_b == 0) {
            return false;
        }
        auto dabx = a.second.x() - b.second.x();
        auto daby = a.second.y() - b.second.y();
        auto ua = (dbx * daby - dby * dabx)/u_b;
        res = Point( a.second.x() + dax * ua, a.second.y() + day * ua );
        return true;
    };

    vector<Segment>
    breakIntersections(const vector<Segment> &segments);


    inline int
    side( const Segment& s, const Point& p ) {
        auto dsx = s.second.x() - s.first.x();
        auto dsy = s.second.y() - s.first.y();
        auto dpx = p.x() - s.first.x();
        auto dpy = p.y() - s.first.y();

        return boost::math::sign( dpy * dsx - dpx * dsy );
    }

    inline Point
    nearestPoint( const Segment& s, const Point& p ) {
        auto dsx = s.second.x() - s.first.x();
        auto dsy = s.second.y() - s.first.y();
        auto dpx = p.x() - s.first.x();
        auto dpy = p.y() - s.first.y();

        auto ab_bb = (dpx * dsx + dpy * dsy) / (dsx * dsx + dsy * dsy);
        return { s.first.x() + dsx * ab_bb, s.first.y() + dsy * ab_bb };
    }


    /////////////////////////////////
    //
    // Polyline operations
    //
    inline void
    addPoint( PolyLine& line, const Point& point ) {
        line += point;
    }

    vector<Segment>
    convertToSegments( const PolyLine& line );

    MultiPolygon
    expand(const PolyLine& line, double distance, int pointsPerCircle = 36 );

    /////////////////////////////////
    //
    // Box operations
    //

    inline bool
    contains( const Polygon& poly, const Box& box ) {
        auto w = box.max_corner().x() - box.min_corner().x();
        auto h = box.max_corner().x() - box.min_corner().x();

        if( bg::within(box.min_corner(), poly) ||
            bg::within(Point(box.min_corner().x(), box.min_corner().y() + h), poly) ||
            bg::within(Point(box.min_corner().x() + w, box.min_corner().y()), poly) ||
            bg::within(box.max_corner(), poly) ) {
            return true;
        }
        return false;
    }


    /////////////////////////////////
    //
    // Polygon operations
    //
    struct PointIndex {
        int i;
        int j;
        double a;
    };

    vector <PointIndex> sortPoints(const Point &position, const vector <Segment> &segments);


    inline int
    numPoints( const Polygon& poly )  {
        return bg::num_points( poly );
    }

    inline bool
    contains( const Polygon& poly, const Point& pt ) {
        return bg::within( pt, poly );
    }

    inline bool
    contains( const Polygon& container, const Polygon& object ) {
        return bg::within( object, container );
    }

    vector<Segment>
    convertToSegments(const Polygon &poly);

    vector<Point>
    convertToExteriorPoints(const Polygon &poly);

    inline void
    addPoint( Polygon& poly, const Point& pt ) {
        bg::append( poly, pt );
    }

    inline void
    closePolygon( Polygon& poly ) {
        bg::correct(poly);
    }

    Polygon
    compute(const Point &position, const vector<Segment> &segments);

    Polygon
    computeViewport(const Point &position, const vector<Segment> &segments, const Point &viewportMinCorner,
                    const Point &viewportMaxCorner);

    //
    // Intersect to polygons, returning the result (which may be empty...).  We're limiting the Boost interface
    // a bit here by only accepting polygons, but for now, that's all we need...

    inline MultiPolygon
    intersectPolygons( const Polygon& lhs, const Polygon& rhs ) {
        MultiPolygon res;
        bg::intersection( lhs, rhs, res );
        return res;
    }

    inline MultiPolygon
    intersectPolygons( const MultiPolygon& lhs, const Polygon& rhs ) {
        MultiPolygon res;
        bg::intersection( lhs, rhs, res );
        return res;
    }

    inline MultiPolygon
    differencePolygons( const MultiPolygon& lhs, const Polygon& rhs ) {
        MultiPolygon res;
        bg::difference( lhs, rhs, res );
        return res;
    }

    inline MultiPolygon
    unionPolygons( const MultiPolygon& lhs, const Polygon& rhs ) {
        MultiPolygon res;
        bg::union_( lhs, rhs, res );
        return res;
    }


    inline MultiPolygon
    unionPolygons( const MultiPolygon& lhs, const MultiPolygon& rhs ) {
        MultiPolygon res;
        bg::union_( lhs, rhs, res );
        return res;
    }

    // TODO: fix these so we don't have a proliferation for every bloody type -- fine for prototyping, but ugly...
    MultiPolygon
    expand(const Polygon& line, double distance, int pointsPerCircle = 36);

    inline bool
    contains( const MultiPolygon& poly, const Point& pt ) {
        return bg::within( pt, poly );
    }

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
