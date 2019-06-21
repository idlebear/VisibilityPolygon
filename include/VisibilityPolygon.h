//
// Created by bjgilhul on 19/06/19.
//

#ifndef OB1_VISIBILITYPOLYGON_H
#define OB1_VISIBILITYPOLYGON_H

#include <Eigen/Dense>

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <memory>

using namespace std;
using namespace Eigen;

namespace Visibility {

    const double VISIBILITYPOLYGON_EPSILON = 0.0000001;

    class Point {
        friend Point operator+( Point lhs, const Point& rhs );
        friend Point operator-( Point lhs, const Point& rhs );

    public:
        explicit Point( const Vector2d &pt) : pt(pt) {};
        Point( const double x, const double y ) : pt{ x, y } {};
        virtual ~Point() {};

        const double
        operator[](int i) const {
          assert(i == 0 || i == 1);
          return pt[i];
        };

        double &
        operator[](int i) {
          assert(i == 0 || i == 1);
          return pt[i];
        };

        Point &
        operator-=(const Point &rhs) {
          pt -= rhs.pt;
          return *this;
        };

        Point &
        operator+=(const Point &rhs) {
          pt -= rhs.pt;
          return *this;
        };

        double squaredDistanceTo(const Point &other) const {
          double dx = pt[0] - other[0];
          double dy = pt[1] - other[1];
          return dx * dx + dy * dy;
        };

        double distanceTo(const Point &other) const {
          return (pt - other.pt).norm();
        };

    private:
        Vector2d pt;
    };

    inline
    bool operator==(const Point &p1, const Point &p2) {
      return p1.distanceTo(p2) < VISIBILITYPOLYGON_EPSILON;
    };

    inline
    bool operator!=(const Point &p1, const Point &p2) {
      return !(p1 == p2);
    };

    inline
    Point
    operator-(Point lhs, const Point &rhs) {
      lhs -= rhs;
      return lhs;
    };

    inline
    Point
    operator+(Point lhs, const Point &rhs) {
      lhs += rhs;
      return lhs;
    }

    class Line {
    public:
        Line(const Point &p1, const Point &p2) : pts{p1, p2} {};

        Point &
        operator[](std::size_t i) {
          assert(i == 0 || i == 1);
          return pts[i];
        };

        const Point &
        operator[](std::size_t i) const {
          assert(i == 0 || i == 1);
          return pts[i];
        };

        Point
        intersectWith( const Line& other ) const throw( int ) {
            auto db = other[1] - other[0];
            auto da = pts[1] - pts[0];

            auto u_b = db[1] * da[0] - db[0] * da[1];
            if( u_b == 0) {
              throw -1;
            }
            Point dab = pts[1] - other[1];
            auto ua = -(db[0] * dab[1] - db[1] * dab[0]) / u_b;
            da[0] *= ua;
            da[1] *= ua;
            return pts[1] + da;
        }

    private:
        vector <Point> pts;

    };

    inline
    bool operator==(const Line &l1, const Line &l2) {
      return (l1[0] == l2[0]) && (l1[1] == l2[1]);
    };

    inline
    bool operator!=(const Line &l1, const Line &l2) {
      return !(l1 == l2);
    };

    inline bool
    inViewport(const Point &position, const Point &viewportMinCorner, const Point &viewportMaxCorner) {
      if (position[0] < viewportMinCorner[0] - VISIBILITYPOLYGON_EPSILON ||
          position[1] < viewportMinCorner[1] - VISIBILITYPOLYGON_EPSILON ||
          position[0] > viewportMaxCorner[0] + VISIBILITYPOLYGON_EPSILON ||
          position[1] > viewportMaxCorner[1] + VISIBILITYPOLYGON_EPSILON) {
        return false;
      }
      return true;
    }


    typedef struct _PointIndex {
        int i;
        int j;
        double a;
    } PointIndex;

    vector <PointIndex> sortPoints(const Point &position, const vector <Line> &segments);

    Point intersectLines(const Line &a, const Line &b);
    bool doLineSegmentsIntersect(const Line &l1, const Line &l2);

    inline double angle(const Point &a, const Point &b) {
      return atan2(b[1] - a[1], b[0] - a[0]);
    };

    inline
    double angle2(const Point &a, const Point &b, const Point &c) {
      auto res = angle(a, b) - angle(b, c);
      if (res < 0) {
        res += 2 * M_PI;
      } else if (res > 2 * M_PI) {
        res -= 2 * M_PI;
      }
      return res;
    };


    class Polygon {
    public:
        explicit Polygon(vector<Point> &pts) : pts(pts) {};

        Polygon() = default;

        Polygon(Point &p1, Point &p2, Point &p3, Point &p4) : pts{p1, p2, p3, p4} {};

        Point &
        operator[](std::size_t i) {
          assert(i < pts.size());
          return pts[i];
        };

        const Point &
        operator[](std::size_t i) const {
          assert(i < pts.size());
          return pts[i];
        };


        Polygon &
        operator+=(const Point &pt) {
          // TODO: After adding the point, we have to decide whether to sort the list (again)
          pts.emplace_back(pt);
        };

        size_t
        size() const {
          return pts.size();
        }

        bool contains(const Point &position) const;
        vector<Line> toSegments() const;

    private:
        vector<Point> pts;
    };


    Polygon
    compute(Point &position, vector<Line> &segments);

    Polygon
    computeViewport(const Point &position, const vector<Line> &segments, const Point &viewportMinCorner,
                    const Point &viewportMaxCorner);

    vector<Line>
    breakIntersections(const vector<Line> &segments);

    inline
    int heapParent(int index) {
      return floor((index - 1) / 2);
    };

    inline
    int heapChild(int index) {
      return 2 * index + 1;
    };

    void remove(int index, vector<int> heap, const Point &position, const vector<Line> &segments,
                const Point &destination, vector<int>& map);

    void
    insert(int index, vector<int> heap, const Point &position, const vector<Line> &segments, const Point &destination,
           vector<int>& map);

    bool lessThan(int i1, int i2, const Point &position, const vector<Line> &segments, const Point &destination);

};



#endif //OB1_VISIBILITYPOLYGON_H
