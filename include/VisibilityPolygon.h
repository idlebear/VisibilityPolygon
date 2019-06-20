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

typedef vector<Vector2d> segment2d;

namespace Visibility {

    const VISIBILITYPOLYGON_EPSILON = 0.0000001;

    class VisibilityPolygon {

        VisibilityPolygon() {

        }

        virtual ~VisibilityPolygon() {

        }


        segment2d compute( Vector2d position, vector<segment2d>segments) {
          vector<segment2d> bounded();
          auto minX = position[0];
          auto minY = position[1];
          auto maxX = position[0];
          auto maxY = position[1];
          for( auto&& segment : segments ) {
            for( auto&& vertex : segment ) {
              minX = min(minX, vertex[0]);
              minY = min(minY, vertex[1]);
              maxX = max(maxX, vertex[0]);
              maxY = max(maxY, vertex[1]);
            }
            bounded.push(segment);
          }

          --minX;
          --minY;
          ++maxX;
          ++maxY;

          segment2d v(2);
          v[0] = Vector2d(minX, minY);
          v[1] = Vector2d(maxX, minY);
          bounded.push( v );
          v[0] = Vector2d(maxX, minY);
          v[1] = Vector2d(maxX, maxY);
          bounded.push( v );
          v[0] = Vector2d(maxX, maxY);
          v[1] = Vector2d(minX, maxY);
          bounded.push( v );
          v[0] = Vector2d(minX, maxY);
          v[1] = Vector2d(minX, minY);
          bounded.push( v );


          segment2d polygon();
          auto sorted = sortPoints( position, bounded );
          auto map = Array( bounded.size() );

          for ( auto it = map.begin(); it != map.end(); ++it ) {
            *it = -1;
          }


          segment2d heap();
          make_heap( heap.begin(), heap.end() );
          auto start = Vector2d( position[0] + 1, position[1] );

          int i = 0;
          for ( auto&& segment: bounded ) {
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
              insert(i, heap, position, bounded, start, map);
            }
            i++;
          }

          for( int i = 0; i < sorted.size(); ) {
            bool extend = false;
            bool shorten = false;
            auto orig = i;
            Vector2d vertex = bounded[sorted[i][0]][sorted[i][1]];
            vector<Vector2d> old_segment = heap[0];
            do {
              if (map[sorted[i][0]] != -1) {
                if (sorted[i][0] == old_segment) {
                  extend = true;
                  vertex = bounded[sorted[i][0]][sorted[i][1]];
                }
                VisibilityPolygon.remove(map[sorted[i][0]], heap, position, bounded, vertex, map);
              } else {
                VisibilityPolygon.insert(sorted[i][0], heap, position, bounded, vertex, map);
                if (heap[0] != old_segment) {
                  shorten = true;
                }
              }
              ++i;
              if (i == sorted.length) break;
            } while (sorted[i][2] < sorted[orig][2] + VISIBILITYPOLYGON_EPSILON);

            if (extend) {
              polygon.push(vertex);
              var cur = VisibilityPolygon.intersectLines(bounded[heap[0]][0], bounded[heap[0]][1], position, vertex);
              if (!VisibilityPolygon.equal(cur, vertex)) polygon.push(cur);
            } else if (shorten) {
              polygon.push(VisibilityPolygon.intersectLines(bounded[old_segment][0], bounded[old_segment][1], position, vertex));
              polygon.push(VisibilityPolygon.intersectLines(bounded[heap[0]][0], bounded[heap[0]][1], position, vertex));
            }
          }
          return polygon;
        };

        VisibilityPolygon.computeViewport = function(position, segments, viewportMinCorner, viewportMaxCorner) {
          var brokenSegments = [];
          var viewport = [[viewportMinCorner[0],viewportMinCorner[1]],[viewportMaxCorner[0],viewportMinCorner[1]],[viewportMaxCorner[0],viewportMaxCorner[1]],[viewportMinCorner[0],viewportMaxCorner[1]]];
          for (var i = 0; i < segments.length; ++i) {
            if (segments[i][0][0] < viewportMinCorner[0] && segments[i][1][0] < viewportMinCorner[0]) continue;
            if (segments[i][0][1] < viewportMinCorner[1] && segments[i][1][1] < viewportMinCorner[1]) continue;
            if (segments[i][0][0] > viewportMaxCorner[0] && segments[i][1][0] > viewportMaxCorner[0]) continue;
            if (segments[i][0][1] > viewportMaxCorner[1] && segments[i][1][1] > viewportMaxCorner[1]) continue;
            var intersections = [];
            for (var j = 0; j < viewport.length; ++j) {
              var k = j + 1;
              if (k == viewport.length) k = 0;
              if (VisibilityPolygon.doLineSegmentsIntersect(segments[i][0][0], segments[i][0][1], segments[i][1][0], segments[i][1][1], viewport[j][0], viewport[j][1], viewport[k][0], viewport[k][1])) {
                var intersect = VisibilityPolygon.intersectLines(segments[i][0], segments[i][1], viewport[j], viewport[k]);
                if (intersect.length != 2) continue;
                if (VisibilityPolygon.equal(intersect, segments[i][0]) || VisibilityPolygon.equal(intersect, segments[i][1])) continue;
                intersections.push(intersect);
              }
            }
            var start = [segments[i][0][0], segments[i][0][1]];
            while (intersections.length > 0) {
              var endIndex = 0;
              var endDis = VisibilityPolygon.distance(start, intersections[0]);
              for (var j = 1; j < intersections.length; ++j) {
                var dis = VisibilityPolygon.distance(start, intersections[j]);
                if (dis < endDis) {
                  endDis = dis;
                  endIndex = j;
                }
              }
              brokenSegments.push([[start[0], start[1]], [intersections[endIndex][0], intersections[endIndex][1]]]);
              start[0] = intersections[endIndex][0];
              start[1] = intersections[endIndex][1];
              intersections.splice(endIndex, 1);
            }
            brokenSegments.push([start, [segments[i][1][0], segments[i][1][1]]]);
          }

          var viewportSegments = [];
          for (var i = 0; i < brokenSegments.length; ++i) {
            if (VisibilityPolygon.inViewport(brokenSegments[i][0], viewportMinCorner, viewportMaxCorner) && VisibilityPolygon.inViewport(brokenSegments[i][1], viewportMinCorner, viewportMaxCorner)) {
              viewportSegments.push([[brokenSegments[i][0][0], brokenSegments[i][0][1]], [brokenSegments[i][1][0], brokenSegments[i][1][1]]]);
            }
          }
          var eps = VISIBILITYPOLYGON_EPSILON * 10;
          viewportSegments.push([[viewportMinCorner[0]-eps,viewportMinCorner[1]-eps],[viewportMaxCorner[0]+eps,viewportMinCorner[1]-eps]]);
          viewportSegments.push([[viewportMaxCorner[0]+eps,viewportMinCorner[1]-eps],[viewportMaxCorner[0]+eps,viewportMaxCorner[1]+eps]]);
          viewportSegments.push([[viewportMaxCorner[0]+eps,viewportMaxCorner[1]+eps],[viewportMinCorner[0]-eps,viewportMaxCorner[1]+eps]]);
          viewportSegments.push([[viewportMinCorner[0]-eps,viewportMaxCorner[1]+eps],[viewportMinCorner[0]-eps,viewportMinCorner[1]-eps]]);
          return VisibilityPolygon.compute(position, viewportSegments);
        }

        VisibilityPolygon.inViewport = function(position, viewportMinCorner, viewportMaxCorner) {
          if (position[0] < viewportMinCorner[0] - VISIBILITYPOLYGON_EPSILON) return false;
          if (position[1] < viewportMinCorner[1] - VISIBILITYPOLYGON_EPSILON) return false;
          if (position[0] > viewportMaxCorner[0] + VISIBILITYPOLYGON_EPSILON) return false;
          if (position[1] > viewportMaxCorner[1] + VISIBILITYPOLYGON_EPSILON) return false;
          return true;
        }

        VisibilityPolygon.inPolygon = function(position, polygon) {
          var val = polygon[0][0];
          for (var i = 0; i < polygon.length; ++i) {
            val = Math.min(polygon[i][0], val);
            val = Math.min(polygon[i][1], val);
          }
          var edge = [val-1, val-1];
          var parity = 0;
          for (var i = 0; i < polygon.length; ++i) {
            var j = i + 1;
            if (j == polygon.length) j = 0;
            if (VisibilityPolygon.doLineSegmentsIntersect(edge[0], edge[1], position[0], position[1], polygon[i][0], polygon[i][1], polygon[j][0], polygon[j][1])) {
              var intersect = VisibilityPolygon.intersectLines(edge, position, polygon[i], polygon[j]);
              if (VisibilityPolygon.equal(position, intersect)) return true;
              if (VisibilityPolygon.equal(intersect, polygon[i])) {
                if (VisibilityPolygon.angle2(position, edge, polygon[j]) < 180) ++parity;
              } else if (VisibilityPolygon.equal(intersect, polygon[j])) {
                if (VisibilityPolygon.angle2(position, edge, polygon[i]) < 180) ++parity;
              } else {
                ++parity;
              }
            }
          }
          return (parity%2)!=0;
        };

        vector<segment2d> convertToSegments( const vector<segment2d>&  polygons ) {
          vector<segment2d> segments();
          for( auto&& poly: polygons ) {
            for( int j = 0; j < poly.size(); ++j ) {
              int k = (j+1)%poly.size();
              segments.push( segment2d( Vector2d( poly[j][0], poly[j][1] ), Vector2d( poly[k][0], poly[k][1] ) ) );
            }
          }
          return segments;
        };

        vector<segment2d> breakIntersections( const vector<segment2d>& segments ) {
          vector<segment2d> output();
          for (int i = 0; i < segments.size(); ++i) {
            segment2d intersections();
            for (int j = 0; j < segments.size(); ++j) {
              if (i == j) {
                continue;
              }
              if ( doLineSegmentsIntersect( segments[i], segments[j] ) {
                Vector2d intersect = intersectLines(segments[i], segments[j] );
                if (intersect.size() != 2) {
                  continue;
                }
                if ( equal(intersect, segments[i][0]) || equal(intersect, segments[i][1])) {
                  continue;
                }
                intersections.push(intersect);
              }
            }
            Vector2d start = segments[i][0];
            while (intersections.size() > 0) {
              var endIndex = 0;
              var endDis = (start - intersections[0]).norm();
              for (var j = 1; j < intersections.length; ++j) {
                var dis = (start - intersections[j]).norm();
                if (dis < endDis) {
                  endDis = dis;
                  endIndex = j;
                }
              }
              output.push( segment2d( start, intersections[endIndex] ) );
              start = intersections[endIndex];
              intersections.splice(endIndex, 1);
            }
            output.push([start, [segments[i][1][0], segments[i][1][1]]]);
          }
          return output;
        };

        bool equal( const Vector2d& a, const Vector2d& b) {
          if( (a - b).norm() < VISIBILITYPOLYGON_EPSILON ) {
            return true;
          }
          return false;
        };

        void
        remove( int index, vector<segment2d>& heap, const Vector2d& position, vector<vector<Vector2d>>& segments, destination, map) {
          map[heap[index]] = -1;
          if (index == heap.length - 1) {
            heap.pop();
            return;
          }
          heap[index] = heap.pop();
          map[heap[index]] = index;
          var cur = index;
          var parent = VisibilityPolygon.parent(cur);
          if (cur != 0 && VisibilityPolygon.lessThan(heap[cur], heap[parent], position, segments, destination)) {
            while (cur > 0) {
              var parent = VisibilityPolygon.parent(cur);
              if (!VisibilityPolygon.lessThan(heap[cur], heap[parent], position, segments, destination)) {
                break;
              }
              map[heap[parent]] = cur;
              map[heap[cur]] = parent;
              var temp = heap[cur];
              heap[cur] = heap[parent];
              heap[parent] = temp;
              cur = parent;
            }
          } else {
            while (true) {
              var left = VisibilityPolygon.child(cur);
              var right = left + 1;
              if (left < heap.length && VisibilityPolygon.lessThan(heap[left], heap[cur], position, segments, destination) &&
                  (right == heap.length || VisibilityPolygon.lessThan(heap[left], heap[right], position, segments, destination))) {
                map[heap[left]] = cur;
                map[heap[cur]] = left;
                var temp = heap[left];
                heap[left] = heap[cur];
                heap[cur] = temp;
                cur = left;
              } else if (right < heap.length && VisibilityPolygon.lessThan(heap[right], heap[cur], position, segments, destination)) {
                map[heap[right]] = cur;
                map[heap[cur]] = right;
                var temp = heap[right];
                heap[right] = heap[cur];
                heap[cur] = temp;
                cur = right;
              } else break;
            }
          }
        };

        VisibilityPolygon.insert = function(index, heap, position, segments, destination, map) {
          var intersect = VisibilityPolygon.intersectLines(segments[index][0], segments[index][1], position, destination);
          if (intersect.length == 0) return;
          var cur = heap.length;
          heap.push(index);
          map[index] = cur;
          while (cur > 0) {
            var parent = VisibilityPolygon.parent(cur);
            if (!VisibilityPolygon.lessThan(heap[cur], heap[parent], position, segments, destination)) {
              break;
            }
            map[heap[parent]] = cur;
            map[heap[cur]] = parent;
            var temp = heap[cur];
            heap[cur] = heap[parent];
            heap[parent] = temp;
            cur = parent;
          }
        };

        bool lessThan( const segment2d& s1, const segment2d& s2, position, destination) {
          auto inter1 = intersectLines( s1, position, destination);
          auto inter2 = intersectLines( s2, position, destination);
          if (!equal(inter1, inter2)) {
            return (inter1 - position).norm() < (inter2 - position).norm();
          }
          int end1 = 0;
          if (equal(inter1, s1[0])) {
            end1 = 1;
          }
          var end2 = 0;
          if (equal(inter2, s2[0])) {
            end2 = 1;
          }
          var a1 = angle2(s1[end1], inter1, position);
          var a2 = angle2(s2[end2], inter2, position);
          if (a1 < M_PI) {
            if (a2 > M_PI) {
              return true;
            }
            return a2 < a1;
          }
          return a1 < a2;
        };

        double angle2( const Vector2d& a, const Vector2d& b, const Vector2d& c) {
          auto res = angle(a,b) - angle(b,c);
          if( res < 0 ) {
            res += 2 * M_PI;
          } else if( res > 2 * M_PI ) {
            res -= 2 * M_PI;
          }
          return a;
        };

        vector<Vector2d> sortPoints( const Vector2d& position, const vector<segment2d> segments) {
          auto points = vector<vector<double>>( segments.size() * 2 );
          for( int i = 0; i < segments.size(); ++i ) {
            for( int j = 0; j < 2; ++j) {
              auto a = angle(segments[i][j], position);
              points[2*i+j] = [i, j, a];
            }
          }
          sort( points.begin(), points.end(), []( const vector<double>& a, const vector<double>&b ) {
            return a[2]-b[2];
          });
          return points;
        };

        double angle( const Vector2d& a, const Vector2d& b) {
          return atan2( b[1]-a[1], b[0]-a[0]);
        };

        Vector2d intersectLines( const segment2d& a, const segment2d& b ) {
          auto db = b[1] - b[0];
          auto da = a[1] - a[0];

          auto u_b  = db[1] * da[0] - db[0] * da[1];
          if (u_b != 0) {
            Vector2d dab = a[1] - b[1];
            auto ua = (db[0] * dab[1] - db[1] * dab[0]) / u_b;
            return a[1] - (ua * -da);
          }
          return Vector2d( NAN, NAN );
        };

        VisibilityPolygon.distance = function(a, b) {
          var dx = a[0]-b[0];
          var dy = a[1]-b[1];
          return dx*dx + dy*dy;
        };

        VisibilityPolygon.doLineSegmentsIntersect = function(x1, y1, x2, y2, x3, y3, x4, y4) {
          var s1_x = x2 - x1;
          var s1_y = y2 - y1;
          var s2_x = x4 - x3;
          var s2_y = y4 - y3;

          var denom = (-s2_x * s1_y + s1_x * s2_y);
          if(denom === 0.0) {
            return false;
          }

          var s = (-s1_y * (x1 - x3) + s1_x * (y1 - y3)) / denom;
          if(s < 0 || s > 1) {
            return false;
          }

          var t = (s2_x * (y1 - y3) - s2_y * (x1 - x3)) / denom;
          return (t >= 0.0 && t <= 1.0);
        }


    };

}



#endif //OB1_VISIBILITYPOLYGON_H
