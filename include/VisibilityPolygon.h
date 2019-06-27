//
// Created by bjgilhul on 19/06/19.
//

#ifndef OB1_VISIBILITYPOLYGON_H
#define OB1_VISIBILITYPOLYGON_H

#include <limits>
#include <Eigen/Dense>

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <memory>
#include <stdexcept>

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

        Point &
        operator*=( double scale ) {
            pt *= scale;
            return *this;
        };

        Point &
        operator/=( double scale ) {
            pt /= scale ;
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

    inline
    Point
    operator*(Point lhs, double rhs) {
        lhs *= rhs;
        return lhs;
    }

    inline
    Point
    operator/(Point lhs, double rhs) {
        lhs /= rhs;
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
        intersectWith( const Line& other ) const {
            auto db = other[1] - other[0];
            auto da = pts[1] - pts[0];

            auto u_b = db[1] * da[0] - db[0] * da[1];
            if( u_b == 0) {
              throw domain_error( "No intersection" );
            }
            Point dab = pts[1] - other[1];
            auto ua = (db[1] * dab[0] - db[0] * dab[1]) / u_b;
            da *= ua;
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
        explicit Polygon(vector<Point> &pts) : pts(pts), centre(findCenter()) {
        }

        Polygon() : pts(), centre(Point(0, 0)) {
        }

        Polygon(const Point &p1, const Point &p2, const Point &p3, const Point &p4) : pts{p1, p2, p3, p4},
                                                                                      centre(findCenter()) {
        }

        virtual ~Polygon() = default;

        Point &
        operator[](std::size_t i) {
            assert(i < pts.size());
            return pts[i];
        }

        const Point &
        operator[](std::size_t i) const {
            assert(i < pts.size());
            return pts[i];
        }

        Polygon &
        operator+=(const Point &pt) {
            // TODO: After adding the point, we have to decide whether to sort the list (again)
            pts.emplace_back(pt);
            return *this;
        }

        size_t
        size() const {
            return pts.size();
        }

        Point
        findCenter() const {
            Point centre(0, 0);
            for (auto const &pt: pts) {
                centre += pt;
            }
            centre /= pts.size();

            return centre;
        }

        Point
        getCentre() const {
            return centre;
        }

        Point
        intersectWithSegment( const Line& segment ) const {
            auto edges = toSegments();
            for( auto const& edge: edges ) {

                Point s1 = segment[1] - segment[0];
                Point s2 = edge[1] - edge[0];

                double denom = s1[0] * s2[1] - s2[0] * s1[1];
                if (denom == 0.0) {
                    continue;
                }

                Point tmp = segment[0] - edge[0];
                double s = (s1[0] * tmp[1] - s1[1] * tmp[0]) / denom;
                if (s < 0.0 || s > 1.0) {
                    continue;
                }
                double t = (s2[0] * tmp[1] - s2[1] * tmp[0]) / denom;
                if( t < 0.0 || t > 1.0) {
                    continue;
                }
                return edge[1] + s1 * t;
            }
            throw domain_error( "No intersection" );  // no intersection...
        }

        double
        distanceTo(const Point &pt) const {
            auto edges = toSegments();
            double minDist = std::numeric_limits<double>::infinity();

            for (auto const &edge: edges) {
                auto vec = edge[1] - edge[0];

                // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
                auto dist = (abs(vec[1] * pt[0] - vec[0] * pt[1] + edge[1][0] * edge[0][1]
                                 - edge[1][1] * edge[0][0])
                            ) / (edge[1].distanceTo(edge[0]));
                if (dist < minDist) {
                    minDist = dist;
                }
            }
            return minDist;
        }

        bool contains(const Point &position) const;
        vector<Line> toSegments() const;

    private:
        vector<Point> pts;
        Point centre;
    };


    Polygon
    compute(const Point &position, const vector<Line> &segments);

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

    void remove(int index, vector<int>& heap, const Point &position, const vector<Line> &segments,
                const Point &destination, vector<int>& map);

    void
    insert(int index, vector<int>& heap, const Point &position, const vector<Line> &segments, const Point &destination,
           vector<int>& map);

    bool lessThan(int i1, int i2, const Point &position, const vector<Line> &segments, const Point &destination);

    vector<Line>
    convertToSegments( const vector<Polygon>& polygons );
};



#endif //OB1_VISIBILITYPOLYGON_H
