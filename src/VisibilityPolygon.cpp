//
// Created by bjgilhul on 19/06/19.
//

#include "../include/VisibilityPolygon.h"
#include <memory>

namespace Visibility {
    bool
    Polygon::contains(const Point &position) const {

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
          try {
            Point intersect = l1.intersectWith(l2);
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
          } catch ( ... ) {
            // TODO -- this case should never happen given that we check for an intersection
            // before calling the actual calculation...
            assert( false );
          }

        }
      }
      return parity % 2 != 0;
    }

    vector <Line>
    Polygon::toSegments() const {
      vector <Line> segments;

      auto n = pts.size();
      for (int j = 0; j < n; ++j) {
          int k = (j + 1) % n;
          segments.emplace_back(Line(Point(pts[j][0], pts[j][1]), Point(pts[k][0], pts[k][1])));
      }
      return segments;
    }


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

      // TODO: At some point this should be converted to use the STL heap (?)
      vector <int> heap;
      auto vertex = Point(position[0] + 1, position[1]);
      vector<int> map(bounded.size(), -1);

      int c = 0;
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
          insert(c, heap, position, bounded, vertex, map);
        }
        c++;
      }

      try {
        for (int i = 0; i < sorted.size();) {
          bool extend = false;
          bool shorten = false;
          auto orig = i;
          vertex = bounded[sorted[i].i][sorted[i].j];
          int old_segment = heap[0];
          do {
            if( map[sorted[i].i] != -1) {
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
            Point cur = bounded[heap[0]].intersectWith(Line(position, vertex));
            if (cur != vertex) {
              polygon += cur;
            }
          } else if (shorten) {
            polygon += bounded[old_segment].intersectWith(Line(position, vertex));
            polygon += bounded[heap[0]].intersectWith(Line(position, vertex));
          }
        }
      } catch (...) {
        // TODO: Should only be one throw here -- should have had an intersection that didn't exist
        // for some reason.  Figure out proper recovery...
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
            try {
              Point intersect = segment.intersectWith(Line(viewport[j], viewport[k]));
              if (intersect == segment[0] || intersect == segment[1]) {
                continue;
              }
              intersections.push_back(intersect);
            } catch( ... ) {
              // No intersection found -- should never happen given the previous check -- just
              // ignore it and continue
            }
          }
        }

        Point start = segment[0];
        while (!intersections.empty()) {
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
      auto eps = VISIBILITYPOLYGON_EPSILON * 10.0;
      viewportSegments.emplace_back(Line(Point(viewportMinCorner[0] - eps, viewportMinCorner[1] - eps),
                                      Point(viewportMaxCorner[0] + eps, viewportMinCorner[1] - eps)));
      viewportSegments.emplace_back(Line(Point(viewportMaxCorner[0] + eps, viewportMinCorner[1] - eps),
                                      Point(viewportMaxCorner[0] + eps, viewportMaxCorner[1] + eps)));
      viewportSegments.emplace_back(Line(Point(viewportMaxCorner[0] + eps, viewportMaxCorner[1] + eps),
                                      Point(viewportMinCorner[0] - eps, viewportMaxCorner[1] + eps)));
      viewportSegments.emplace_back(Line(Point(viewportMinCorner[0] - eps, viewportMaxCorner[1] + eps),
                                      Point(viewportMinCorner[0] - eps, viewportMinCorner[1] - eps)));
      return compute(position, viewportSegments);
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
            try {
              Point intersect = segments[i].intersectWith( segments[j] );
              if( intersect[0] == NAN ) {
                // TODO: Intersect lines should never fail here -- we just checked for intersection in the
                //  line above...
                continue;
              }
              if( intersect == segments[i][0] || intersect == segments[i][1] ) {
                continue;
              }
              intersections.emplace_back(intersect);
            } catch ( ... ) {
              // Nothing to clean up -- just ignore the bad recommendation...
            }
          }
        }
        Point start = segments[i][0];
        while (!intersections.empty()) {
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

    void remove( int index, vector<int>& heap, const Point& position, const vector<Line>& segments,
                 const Point& destination, vector<int>& map ) {
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
          parent = heapParent(cur);
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
    insert(int index, vector<int>& heap, const Point &position, const vector <Line> &segments, const Point &destination,
           vector<int>& map) {
      try {
        Point intersect = segments[index].intersectWith( Line( position, destination ) );
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
      } catch (...) {
        // bad insertion request -- ignored
      }
    }

    bool lessThan(int i1, int i2, const Point &position, const vector <Line> &segments, const Point &destination) {

      Line s1 = segments[i1];
      Line s2 = segments[i2];

      Line destLine = Line(position, destination);
      try{
        auto inter1 = s1.intersectWith(destLine);
        auto inter2 = s2.intersectWith(destLine);
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
      } catch ( ... ) {
        // TODO: Bad intersection request/calculation -- unable to compare the two line segments
        return false;
      }
    }


    vector <PointIndex> sortPoints(const Point &position, const vector <Line> &segments) {
      vector < PointIndex > points;
      points.reserve(segments.size() * 2);

      for (int i = 0; i < segments.size(); ++i) {
        for (int j = 0; j < 2; ++j) {
          double a = angle(segments[i][j], position);
          points.emplace_back( PointIndex{i, j, a} );
        }
      }
      sort(points.begin(), points.end(), [](const PointIndex& a, const PointIndex& b) {
          return a.a < b.a;
      });
      return points;
    };

    bool doLineSegmentsIntersect(const Line &l1, const Line &l2) {
      Point s1 = l1[1] - l1[0];
      Point s2 = l2[1] - l2[0];

      double denom = s1[0] * s2[1] - s2[0] * s1[1];
      if (denom == 0.0) {
        return false;
      }

      Point tmp = l1[0] - l2[0];
      double s = (s1[0] * tmp[1] - s1[1] * tmp[0]) / denom;
      if (s < 0.0 || s > 1.0) {
        return false;
      }
      double t = (s2[0] * tmp[1] - s2[1] * tmp[0]) / denom;
      return (t >= 0.0 && t <= 1.0);
    }

    vector<Line>
    convertToSegments( const vector<Polygon>& polygons ) {
        vector<Line> combinedSegments;

        for( auto const& poly: polygons ) {
            auto segments = poly.toSegments();
            combinedSegments.insert(
                combinedSegments.end(),
                std::make_move_iterator(segments.begin()),
                std::make_move_iterator(segments.end())
            );
        }
        return combinedSegments;
    }

}