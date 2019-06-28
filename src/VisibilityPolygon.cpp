//
// Created by bjgilhul on 19/06/19.
//

#include "../include/VisibilityPolygon.h"
#include <memory>

namespace Visibility {


    vector<Segment>
    convertToSegments(const Polygon &poly) {
        vector<Segment> segments;
        vector<Point> points;

        for( auto const& point : bg::exterior_ring(poly) ) {
            points.emplace_back(point);
        }

        auto n = points.size();
        for (int j = 0; j < n; ++j) {
            int k = (j + 1) % n;
            segments.emplace_back(Segment(points[j], points[k]));
        }
        return segments;
    }

    vector<Point>
    convertToPoints( const Segment& segment ) {
        vector<Point> res;
        res.emplace_back( Point( bg::get<0,0>(segment), bg::get<0,1>(segment) ) );
        res.emplace_back( Point( bg::get<1,0>(segment), bg::get<1,1>(segment) ) );
        return res;
    }

    vector<Point>
    convertToPoints( const Polygon& poly ) {
        vector<Point> res;

        for(auto it = boost::begin(bg::exterior_ring(poly)); it != boost::end(bg::exterior_ring(poly)); ++it) {
            double x = bg::get<0>(*it);
            double y = bg::get<1>(*it);
            //use the coordinates...
            res.emplace_back( Point( x, y ) );
        }
        return res;
    }

    vector<Segment>
    convertToSegments(const vector<Polygon> &polygons) {
        vector<Segment> combinedSegments;

        for (auto const &poly: polygons) {
            auto segments = convertToSegments( poly );
            combinedSegments.insert(
                    combinedSegments.end(),
                    std::make_move_iterator(segments.begin()),
                    std::make_move_iterator(segments.end())
            );
        }
        return combinedSegments;
    }


    Polygon
    compute(const Point &position, const vector<Segment> &segments) {
        vector<Segment> bounded;
        auto minX = position.x();
        auto minY = position.y();
        auto maxX = position.x();
        auto maxY = position.y();
        for (auto const &segment : segments) {

            auto pts = convertToPoints( segment );
            for( auto const &pt : pts ) {
                minX = min(minX, pt.x());
                minY = min(minY, pt.y());
                maxX = max(maxX, pt.x());
                maxY = max(maxY, pt.y());
            }
            bounded.emplace_back(segment);
        }

        --minX;
        --minY;
        ++maxX;
        ++maxY;

        bounded.emplace_back(Segment(Point(minX, minY), Point(maxX, minY)));
        bounded.emplace_back(Segment(Point(maxX, minY), Point(maxX, maxY)));
        bounded.emplace_back(Segment(Point(maxX, maxY), Point(minX, maxY)));
        bounded.emplace_back(Segment(Point(minX, maxY), Point(minX, minY)));

        Polygon polygon;

        auto sorted = sortPoints(position, bounded);

        // TODO: At some point this should be converted to use the STL heap (?)
        vector<int> heap;
        auto vertex = Point(position.x() + 1, position.y());
        vector<int> map(bounded.size(), -1);

        int c = 0;
        for (auto const &segment: bounded) {
            auto pts = convertToPoints( segment );
            auto a1 = angle(pts[0], position);
            auto a2 = angle(pts[1], position);
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

        for (int i = 0; i < sorted.size();) {
            bool extend = false;
            bool shorten = false;
            auto orig = i;
            vertex = convertToPoints(bounded[sorted[i].i])[sorted[i].j];
            int old_segment = heap[0];
            do {
                if (map[sorted[i].i] != -1) {
                    if (sorted[i].i == old_segment) {
                        extend = true;
                        vertex = convertToPoints(bounded[sorted[i].i])[sorted[i].j];
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
            } while (sorted[i].a < sorted[orig].a + VISIBILITY_POLYGON_EPSILON);

            if (extend) {
                bg::append(polygon, vertex);
                Point cur;
                if (intersectLines(bounded[heap[0]], Segment(position, vertex), cur)) {
                    if( cur != vertex ) {
                        bg::append(polygon, cur);
                    }
                }
            } else if (shorten) {
                Point pt;
                if (intersectLines(bounded[old_segment], Segment(position, vertex), pt)) {
                    bg::append( polygon, pt );
                }
                if (intersectLines(bounded[heap[0]], Segment(position, vertex), pt)) {
                    bg::append( polygon, pt );
                }
            }
        }
        bg::correct(polygon);
        return polygon;
    }

    Polygon
    computeViewport(const Point &position, const vector <Segment> &segments, const Point &viewportMinCorner,
                    const Point &viewportMaxCorner) {
        vector <Segment> brokenSegments;
        Point viewport[4] = {Point(viewportMinCorner.x(), viewportMinCorner.y()),
                           Point(viewportMaxCorner.x(), viewportMinCorner.y()),
                           Point(viewportMaxCorner.x(), viewportMaxCorner.y()),
                           Point(viewportMinCorner.x(), viewportMaxCorner.y())};

        for (auto const &segment: segments) {
            auto pts = convertToPoints( segment );
            if (pts[0].x() < viewportMinCorner.x() && pts[1].x() < viewportMinCorner.x()) { continue; }
            if (pts[0].y() < viewportMinCorner.y() && pts[1].y() < viewportMinCorner.y()) { continue; }
            if (pts[0].x() > viewportMaxCorner.x() && pts[1].x() > viewportMaxCorner.x()) { continue; }
            if (pts[0].y() > viewportMaxCorner.y() && pts[1].y() > viewportMaxCorner.y()) { continue; }

            vector<Point> intersections;
            for (int j = 0; j < 4; ++j) {
                int k = (j + 1) % 4;
                vector<Point> newInts;
                if( intersectSegments(segment, Segment(viewport[j], viewport[k]), newInts)) {
                    for (auto const &intersect: newInts) {
                        if (intersect != pts[0] && intersect != pts[1]) {
                            intersections.push_back(intersect);
                        }
                    }
                }
            }

            Point start( pts[0] );
            while (!intersections.empty()) {
                int endIndex = 0;
                double endDis = bg::distance(start, intersections[0]);
                for (int j = 1; j < intersections.size(); ++j) {
                    double dis = bg::distance(start, intersections[j]);
                    if (dis < endDis) {
                        endDis = dis;
                        endIndex = j;
                    }
                }
                brokenSegments.emplace_back(Segment(start, intersections[endIndex]));
                start = intersections[endIndex];
                intersections.erase(intersections.begin() + endIndex);
            }
            brokenSegments.emplace_back(Segment(start, pts[1]));
        }

        vector<Segment> viewportSegments;
        for (auto const &segment : brokenSegments) {
            auto pts = convertToPoints( segment );
            if (inViewport(pts[0], viewportMinCorner, viewportMaxCorner) &&
                inViewport(pts[1], viewportMinCorner, viewportMaxCorner)) {
                viewportSegments.push_back(segment);
            }
        }
        auto eps = VISIBILITY_POLYGON_EPSILON * 10.0;
        viewportSegments.emplace_back(Segment(Point(viewportMinCorner.x() - eps, viewportMinCorner.y() - eps),
                                              Point(viewportMaxCorner.x() + eps, viewportMinCorner.y() - eps)));
        viewportSegments.emplace_back(Segment(Point(viewportMaxCorner.x() + eps, viewportMinCorner.y() - eps),
                                              Point(viewportMaxCorner.x() + eps, viewportMaxCorner.y() + eps)));
        viewportSegments.emplace_back(Segment(Point(viewportMaxCorner.x() + eps, viewportMaxCorner.y() + eps),
                                              Point(viewportMinCorner.x() - eps, viewportMaxCorner.y() + eps)));
        viewportSegments.emplace_back(Segment(Point(viewportMinCorner.x() - eps, viewportMaxCorner.y() + eps),
                                              Point(viewportMinCorner.x() - eps, viewportMinCorner.y() - eps)));
        return compute(position, viewportSegments);
    }

    vector<Segment> breakIntersections(const vector<Segment> &segments) {
        vector<Segment> output;
        for (int i = 0; i < segments.size(); ++i) {
            vector<Point> pts = convertToPoints(segments[i]);
            vector<Point> intersections;
            for (int j = 0; j < segments.size(); ++j) {
                if (i == j) {
                    continue;
                }
                vector<Point> newIntersections;
                if (intersectSegments(segments[i], segments[j], newIntersections)) {
                    for (auto const &intersect: newIntersections) {
                        if (intersect != pts[0] && intersect != pts[1]) {
                            intersections.emplace_back(intersect);
                        }
                    }
                }
            }
            Point start = pts[0];
            while (!intersections.empty()) {
                int endIndex = 0;
                double endDis = bg::distance(start, intersections[0]);
                for (int j = 1; j < intersections.size(); ++j) {
                    double dis = bg::distance(start, intersections[j]);
                    if (dis < endDis) {
                        endDis = dis;
                        endIndex = j;
                    }
                }
                output.emplace_back(Segment(start, intersections[endIndex]));
                start = intersections[endIndex];
                intersections.erase(intersections.begin() + endIndex);
            }
            output.emplace_back(Segment(start, pts[1]));
        }
        return output;
    };

    void remove(int index, vector<int> &heap, const Point &position, const vector<Segment> &segments,
                const Point &destination, vector<int> &map) {
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
    insert(int index, vector<int> &heap, const Point &position, const vector<Segment> &segments,
           const Point &destination, vector<int> &map) {
        Point intersection;
        if (intersectLines(segments[index], Segment(position, destination), intersection)) {
            int cur = heap.size();
            heap.push_back(index);
            map[index] = cur;
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
        }
    }

    bool lessThan(int i1, int i2, const Point &position, const vector<Segment> &segments, const Point &destination) {

        Segment s1 = segments[i1];
        Segment s2 = segments[i2];

        Segment destSegment = Segment(position, destination);
        Point inter1;
        if (!intersectLines(s1, destSegment, inter1)) {
            return false;
        }
        Point inter2;
        if (!intersectLines(s2, destSegment, inter2)) {
            return true;
        }
        if (inter1 != inter2) {
            return bg::distance(inter1, position) < bg::distance(inter2, position);
        }
        int end1 = 0;
        auto s1pts = convertToPoints(s1);
        if (inter1 == s1pts[0]) {
            end1 = 1;
        }
        int end2 = 0;
        auto s2pts = convertToPoints(s2);
        if (inter2 == s2pts[0]) {
            end2 = 1;
        }
        double a1 = angle2(s1pts[end1], inter1, position);
        double a2 = angle2(s2pts[end2], inter2, position);
        if (a1 < M_PI) {
            if (a2 > M_PI) {
                return true;
            }
            return a2 < a1;
        }
        return a1 < a2;
    }


    vector<PointIndex> sortPoints(const Point &position, const vector<Segment> &segments) {
        vector<PointIndex> points;
        points.reserve(segments.size() * 2);

        for (int i = 0; i < segments.size(); ++i) {
            auto pts = convertToPoints(segments[i]);
            for (int j = 0; j < 2; ++j) {
                double a = angle(pts[j], position);
                points.emplace_back(PointIndex{i, j, a});
            }
        }
        sort(points.begin(), points.end(), [](const PointIndex &a, const PointIndex &b) {
            return a.a < b.a;
        });
        return points;
    };


}