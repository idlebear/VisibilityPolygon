//
// Created by bjgilhul on 19/06/19.
//

#include "../include/VisibilityPolygon.h"


namespace Visibility {


    Polygon
    compute( const Point& position, const vector<Line>& segments) {
      vector <Line> bounded;
      auto minX = position[0];
      auto minY = position[1];
      auto maxX = position[0];
      auto maxY = position[1];
      for (auto const &segment : segments) {
        for ( int i = 0; i < 2; i++ ) {
          minX = min(minX, segment[i][0]);
          minY = min(minY, segment[i][1]);
          maxX = max(maxX, segment[i][0]);
          maxY = max(maxY, segment[i][1]);
        }
        bounded.emplace_back(segment);
      }

      --minX;
      --minY;
      ++maxX;
      ++maxY;

      bounded.emplace_back(Line(Point(minX, minY), Point(maxX, minY)));
      bounded.emplace_back(Line(Point(maxX, minY), Point(maxX, maxY)));
      bounded.emplace_back(Line(Point(maxX, maxY), Point(minX, maxY)));
      bounded.emplace_back(Line(Point(minX, maxY), Point(minX, minY)));

      Polygon polygon;

      auto sorted = sortPoints(position, bounded);

      // TODO: At some point this should be converted to use the STL heap
      // make_heap( heap.begin(), heap.end(), [&position,&vertex]( const Line& lhs, const Line& rhs ){
      //     return LessThan( lhs, rhs, position, start );
      // } );
      vector <int> heap;
      auto vertex = Point(position[0] + 1, position[1]);
      int map[bounded.size()] = {-1};

      int i = 0;
      for (auto const &segment: bounded) {
        auto a1 = angle(segment[0], position);
        auto a2 = angle(segment[1], position);
        bool active = false;
        if (a1 > -M_PI && a1 <= 0 && a2 <= M_PI && a2 >= 0 && a2 - a1 > M_PI) {
          active = true;
        }
        if (a2 > -M_PI && a2 <= 0 && a1 <= M_PI && a1 >= 0 && a1 - a2 > M_PI) {
          active = true;
        }
        if (active) {
          insert(i, heap, position, bounded, vertex, map);
        }
        i++;
      }

      for (int i = 0; i < sorted.size();) {
        bool extend = false;
        bool shorten = false;
        auto orig = i;
        vertex = bounded[sorted[i].i][sorted[i].j];
        int old_segment = heap[0];
        do {
          if (map[sorted[i].i] != -1) {
            if (sorted[i].i == old_segment) {
              extend = true;
              vertex = bounded[sorted[i].i][sorted[i].j];
            }

            remove(map[sorted[i].i], heap, position, bounded, vertex, map);
          } else {
            insert(sorted[i].i, heap, position, bounded, vertex, map);
            if (heap[0] != old_segment) {
              shorten = true;
            }
          }
          ++i;
          if (i == sorted.size()) {
            break;
          }
        } while (sorted[i].a < sorted[orig].a + VISIBILITYPOLYGON_EPSILON);

        if (extend) {
          polygon += vertex;
          Point cur = intersectLines(bounded[heap[0]], Line(position, vertex));
          if (cur != vertex) {
            polygon += cur;
          }
        } else if (shorten) {
          polygon += intersectLines(bounded[old_segment], Line(position, vertex));
          polygon += intersectLines(bounded[heap[0]], Line(position, vertex));
        }
      }
      return polygon;
    }

    Polygon
    computeViewport(const Point &position, const vector <Line> &segments, const Point &viewportMinCorner,
                    const Point &viewportMaxCorner) {
      vector <Line> brokenSegments;
      Point viewport[4] = {Point(viewportMinCorner[0], viewportMinCorner[1]),
                           Point(viewportMaxCorner[0], viewportMinCorner[1]),
                           Point(viewportMaxCorner[0], viewportMaxCorner[1]),
                           Point(viewportMinCorner[0], viewportMaxCorner[1])};

      for (auto const &segment: segments) {
        if (segment[0][0] < viewportMinCorner[0] && segment[1][0] < viewportMinCorner[0]) continue;
        if (segment[0][1] < viewportMinCorner[1] && segment[1][1] < viewportMinCorner[1]) continue;
        if (segment[0][0] > viewportMaxCorner[0] && segment[1][0] > viewportMaxCorner[0]) continue;
        if (segment[0][1] > viewportMaxCorner[1] && segment[1][1] > viewportMaxCorner[1]) continue;

        vector <Point> intersections;
        for (int j = 0; j < 4; ++j) {
          int k = (j + 1) % 4;
          if( doLineSegmentsIntersect( segment, Line(viewport[j], viewport[k] ) ) ) {
            Point intersect = intersectLines(segment, Line(viewport[j], viewport[k]));
            if (intersect[0] == NAN) {
              // TODO: another bogus check for failure -- implement exceptions and be done with it
              continue;
            }
            if (intersect == segment[0] || intersect == segment[1]) {
              continue;
            }
            intersections.push_back(intersect);
          }
        }

        Point start = segment[0];
        while (intersections.size() > 0) {
          int endIndex = 0;
          double endDis = start.distanceTo(intersections[0]);
          for (int j = 1; j < intersections.size(); ++j) {
            double dis = start.distanceTo(intersections[j]);
            if (dis < endDis) {
              endDis = dis;
              endIndex = j;
            }
          }
          brokenSegments.emplace_back(Line(start, intersections[endIndex]));
          start = intersections[endIndex];
          intersections.erase(intersections.begin() + endIndex);
        }
        brokenSegments.emplace_back(Line(start, segment[1]));
      }

      vector <Line> viewportSegments;
      for (auto const &segment : brokenSegments) {
        if (inViewport(segment[0], viewportMinCorner, viewportMaxCorner) &&
            inViewport(segment[1], viewportMinCorner, viewportMaxCorner)) {
          viewportSegments.push_back(segment);
        }
      }
      int eps = VISIBILITYPOLYGON_EPSILON * 10;
      viewportSegments.push_back(Line(Point(viewportMinCorner[0] - eps, viewportMinCorner[1] - eps),
                                      Point(viewportMaxCorner[0] + eps, viewportMinCorner[1] - eps)));
      viewportSegments.push_back(Line(Point(viewportMaxCorner[0] + eps, viewportMinCorner[1] - eps),
                                      Point(viewportMaxCorner[0] + eps, viewportMaxCorner[1] + eps)));
      viewportSegments.push_back(Line(Point(viewportMaxCorner[0] + eps, viewportMaxCorner[1] + eps),
                                      Point(viewportMinCorner[0] - eps, viewportMaxCorner[1] + eps)));
      viewportSegments.push_back(Line(Point(viewportMinCorner[0] - eps, viewportMaxCorner[1] + eps),
                                      Point(viewportMinCorner[0] - eps, viewportMinCorner[1] - eps)));
      return compute(position, viewportSegments);
    }

    vector <Line> convertToSegments(const vector <Polygon> &polygons) {
      vector <Line> segments;
      for (auto const &poly: polygons) {
        for (int j = 0; j < poly.size(); ++j) {
          int k = (j + 1) % poly.size();
          segments.emplace_back(Line(Point(poly[j][0], poly[j][1]), Point(poly[k][0], poly[k][1])));
        }
      }
      return segments;
    }

    vector<Line> breakIntersections( const vector<Line>& segments ) {
      vector<Line> output;
      for (int i = 0; i < segments.size(); ++i) {
        vector<Point> intersections;
        for (int j = 0; j < segments.size(); ++j) {
          if (i == j) {
            continue;
          }
          if ( doLineSegmentsIntersect( segments[i], segments[j] ) ) {
            Point intersect = intersectLines(segments[i], segments[j] );
            if( intersect[0] == NAN ) {
              // TODO: Intersect lines should never fail here -- we just checked for intersection in the
              //  line above...
              continue;
            }
            if( intersect == segments[i][0] || intersect == segments[i][1] ) {
              continue;
            }
            intersections.emplace_back(intersect);
          }
        }
        Point start = segments[i][0];
        while (intersections.size() > 0) {
          int endIndex = 0;
          double endDis = start.distanceTo( intersections[0] );
          for( int  j = 1; j < intersections.size(); ++j ) {
            double dis = start.distanceTo( intersections[j] );
            if( dis < endDis ) {
              endDis = dis;
              endIndex = j;
            }
          }
          output.emplace_back( Line( start, intersections[endIndex] ) );
          start = intersections[endIndex];
          intersections.erase(intersections.begin() + endIndex );
        }
        output.emplace_back( Line( start, segments[i][1] ) );
      }
      return output;
    };

    void remove( int index, vector<int> heap, const Point& position, const vector<Line>& segments, const Point& destination, int map[] ) {
      map[heap[index]] = -1;
      if (index == heap.size() - 1) {
        heap.pop_back();
        return;
      }
      heap[index] = heap.back();
      heap.pop_back();
      map[heap[index]] = index;
      int cur = index;
      int parent = heapParent(cur);
      if (cur != 0 && lessThan(heap[cur], heap[parent], position, segments, destination)) {
        while (cur > 0) {
          int parent = heapParent(cur);
          if (!lessThan(heap[cur], heap[parent], position, segments, destination)) {
            break;
          }
          map[heap[parent]] = cur;
          map[heap[cur]] = parent;
          int temp = heap[cur];
          heap[cur] = heap[parent];
          heap[parent] = temp;
          cur = parent;
        }
      } else {
        while (true) {
          int left = heapChild(cur);
          int right = left + 1;
          if (left < heap.size() && lessThan(heap[left], heap[cur], position, segments, destination) &&
              (right == heap.size() || lessThan(heap[left], heap[right], position, segments, destination))) {
            map[heap[left]] = cur;
            map[heap[cur]] = left;
            int temp = heap[left];
            heap[left] = heap[cur];
            heap[cur] = temp;
            cur = left;
          } else if (right < heap.size() && lessThan(heap[right], heap[cur], position, segments, destination)) {
            map[heap[right]] = cur;
            map[heap[cur]] = right;
            int temp = heap[right];
            heap[right] = heap[cur];
            heap[cur] = temp;
            cur = right;
          } else {
            break;
          }
        }
      }
    }

    void
    insert(int index, vector<int> heap, const Point &position, const vector <Line> &segments, const Point &destination,
           int map[]) {
      Point intersect = intersectLines( segments[index], Line( position, destination ) );
      if (intersect[0] == NAN) {
        return;
      }
      int cur = heap.size();
      heap.push_back(index);
      map[index] = cur;
      while (cur > 0) {
        int parent = heapParent(cur);
        if( !lessThan(heap[cur], heap[parent], position, segments, destination ) ) {
          break;
        }
        map[heap[parent]] = cur;
        map[heap[cur]] = parent;
        int temp = heap[cur];
        heap[cur] = heap[parent];
        heap[parent] = temp;
        cur = parent;
      }
    }

    bool lessThan(int i1, int i2, const Point &position, const vector <Line> &segments, const Point &destination) {

      Line s1 = segments[i1];
      Line s2 = segments[i2];

      Line destLine = Line(position, destination);
      auto inter1 = intersectLines(s1, destLine);
      auto inter2 = intersectLines(s2, destLine);
      if (inter1 != inter2) {
        return inter1.distanceTo(position) < inter2.distanceTo(position);
      }
      int end1 = 0;
      if (inter1 == s1[0]) {
        end1 = 1;
      }
      int end2 = 0;
      if (inter2 == s2[0]) {
        end2 = 1;
      }
      double a1 = angle2(s1[end1], inter1, position);
      double a2 = angle2(s2[end2], inter2, position);
      if (a1 < M_PI) {
        if (a2 > M_PI) {
          return true;
        }
        return a2 < a1;
      }
      return a1 < a2;
    }


    vector <PointIndex> sortPoints(const Point &position, const vector <Line> &segments) {
      vector < PointIndex > points(segments.size() * 2);
      for (int i = 0; i < segments.size(); ++i) {
        for (int j = 0; j < 2; ++j) {
          double a = angle(segments[i][j], position);
          points[2 * i + j] = PointIndex{i, j, a};
        }
      }
      sort(points.begin(), points.end(), [](const PointIndex &a, const PointIndex &b) {
          return a.a - b.a;
      });
      return points;
    };

    Point intersectLines(const Line &a, const Line &b) {
      auto db = b[1] - b[0];
      auto da = a[1] - a[0];

      auto u_b = db[1] * da[0] - db[0] * da[1];
      if (u_b != 0) {
        Point dab = a[1] - b[1];
        auto ua = -(db[0] * dab[1] - db[1] * dab[0]) / u_b;
        da[0] *= ua;
        da[1] *= ua;
        return a[1] + da;
      }
      return Point(NAN, NAN);
    };

    bool doLineSegmentsIntersect(const Line &l1, const Line &l2) {
      Point s1 = l1[1] - l1[2];
      Point s2 = l2[1] - l2[2];

      double denom = (-s2[0] * s1[1] + s1[0] * s2[1]);
      if (denom == 0.0) {
        return false;
      }

      Point tmp = l1[0] - l2[0];
      double s = (-s1[1] * (tmp[0]) + s1[0] * (tmp[1])) / denom;
      if (s < 0.0 || s > 1.0) {
        return false;
      }
      int t = (s2[0] * (tmp[1]) - s2[1] * (tmp[0])) / denom;
      return (t >= 0.0 && t <= 1.0);
    }


}