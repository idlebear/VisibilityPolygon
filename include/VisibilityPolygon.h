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

    bool operator==(const Point &p1, const Point &p2) {
      return p1.distanceTo(p2) < VISIBILITYPOLYGON_EPSILON;
    };

    bool operator!=(const Point &p1, const Point &p2) {
      return !(p1 == p2);
    };

    Point
    operator-(Point lhs, const Point &rhs) {
      lhs -= rhs;
      return lhs;
    };

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

    private:
        vector <Point> pts;

    };

    bool operator==(const Line &l1, const Line &l2) {
      return (l1[0] == l2[0]) && (l1[1] == l2[1]);
    };

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
          // TODO: Add the point to the polygon, then sort the list (again)
          pts.emplace_back(pt);
        };

        size_t
        size() const {
          return pts.size();
        }

        bool
        contains(const Point &position) {

          double val = pts[0][0];
          for (auto const &pt: pts) {
            val = min(pt[0], val);
            val = min(pt[1], val);
          }
          Point edge(val - 1, val - 1);
          int parity = 0;
          for (int i = 0; i < pts.size(); ++i) {
            int j = (i + 1) % pts.size();

            auto l1 = Line(edge, position);
            auto l2 = Line(pts[i], pts[j]);
            if (doLineSegmentsIntersect(l1, l2)) {
              Point intersect = intersectLines(l1, l2);
              if (position == intersect) {
                return true;
              }
              if ((intersect == pts[i] &&
                   angle2(position, edge, pts[j]) >= M_PI) ||
                  (intersect == pts[j] &&
                   angle2(position, edge, pts[i]) >= M_PI)) {
                continue;
              }
              ++parity;
            }
          }
          return parity % 2 != 0;
        }

    private:
        vector<Point> pts;
    };


    Polygon
    compute(Point &position, vector<Line> &segments);

    Polygon
    computeViewport(const Point &position, const vector<Line> &segments, const Point &viewportMinCorner,
                    const Point &viewportMaxCorner);

    vector<Line>
    convertToSegments(const vector<Polygon> &polygons);

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
                const Point &destination, int map[]);

    void
    insert(int index, vector<int> heap, const Point &position, const vector<Line> &segments, const Point &destination,
           int map[]);

    bool lessThan(int i1, int i2, const Point &position, const vector<Line> &segments, const Point &destination);

};



#endif //OB1_VISIBILITYPOLYGON_H
