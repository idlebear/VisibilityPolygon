#include <iostream>
#include <boost/geometry.hpp>

#include "../include/VisibilityPolygon.h"

using namespace std;

namespace V = Visibility;

int main() {

    vector<V::Polygon> polys{
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
        auto res = V::planMinHeightCoverage( p, 0.25 );
        cout << "** Test " << ++i << endl;
        for( auto const& r : res ) {
            cout << "( " << r.first.x() << ", "<< r.first.y() << ") -- ( "<< r.second.x() << ", " <<  r.second.y() << " )" << endl;
        }
    }
}

