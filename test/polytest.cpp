#include <iostream>
#include <boost/geometry.hpp>

#include "../include/VisibilityPolygon.h"

using namespace std;

namespace V = Visibility;

int main() {

    vector<V::Polygon> polys{
            {{
                    {893.542, 1338},
                    {893.542, 1295},
                    {1004.54, 1273.5},
                    {1078.54, 1402.5},
                    {1004.54, 1402.5},
                    {967.542, 1381},
                    {893.542, 1338}

            }},
            {{
                     {0, 0},
                     {7, 0},
                     {2, 3},
                     {1, 3},
                     {-1, 1},
                     {0,    0}
             }},

            {{
                     {0, 0},
                     {1, 0},
                     {0, 1},
                     {0,    0}
             }},
            {{
                     {0, 0},
                     {1, 0},
                     {1, 1},
                     {0, 1},
                     {0,    0}
             }},
            {{
                     {0, 0},
                     {1, 0},
                     {1, 1},
                     {0, 1},
                     {-0.5, 0.5},
                     {0,   0}
             }},
            {{
                     {0, 2},
                     {3, 1},
                     {5, 1},
                     {7, 4},
                     {5, 6},
                     {0.5, 9},
                     {-1, 8},
                     {-1, 7},
                     {0                                                                                                              , 2}
             }}
    };


    cout << "--- Find Antipods --------------" << endl;
    auto i = 0;
    for(auto const& p : polys ) {
        auto res = V::findAntipodals( p );
        cout << "** Test " << ++i << endl;
        for( auto r : res ) {
            cout << r.first << ", " << r.second << endl;
        }
    }
    cout << "--- Find Heights --------------" << endl;
    i = 0;
    for(auto const& p : polys ) {
        auto res = V::findHeights( p );
        int a, b, c;
        double height;
        cout << "** Test " << ++i << endl;
        for( auto const& r : res ) {
            tie( a, b, c, height ) = r;
            cout << a << ", "<< b << ", "<< c << " -- " << height << endl;
        }
    }
    cout << "--- Find Coverage --------------" << endl;
    i = 0;
    for(auto const& p : polys ) {
        auto res = V::planMinHeightCoverage( p, 50);
        cout << "** Test " << ++i << endl;
        for( auto const& r : res ) {
            cout << "( " << r.first.x() << ", "<< r.first.y() << ") -- ( "<< r.second.x() << ", " <<  r.second.y() << " )" << endl;
        }
    }

    cout << "Segment tests" << endl;
    {
        V::Segment a( {0,0}, {5,1} );
        V::Segment before( {-3, 5}, {-3, -4} );
        V::Segment thru( {3, 3}, {3,-3} );
        V::Segment after( {7, 2}, {7, -2} );

        V::Point inter;
        double aOff, bOff;
        V::intersectSegments( a, before, inter, aOff, bOff );
        cout << "Before: (" << inter.x() << ", " << inter.y() << ") -- Offset A:" << aOff << " Offset B: " << bOff << endl;
        assert( aOff < 0 );

        V::intersectSegments( a, thru, inter, aOff, bOff );
        cout << "Thru: (" << inter.x() << ", " << inter.y() << ") -- Offset A:" << aOff << " Offset B: " << bOff << endl;
        assert( aOff >= 0 && aOff <= 1 );

        V::intersectSegments( a, after, inter, aOff, bOff );
        cout << "After: (" << inter.x() << ", " << inter.y() << ") -- Offset A:" << aOff << " Offset B: " << bOff << endl;
        assert( aOff > 1 );

    }


}

