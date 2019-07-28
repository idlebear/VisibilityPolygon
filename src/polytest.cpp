#include <iostream>
#include <boost/geometry.hpp>

#include "../include/VisibilityPolygon.h"

using namespace std;

namespace V = Visibility;

int main() {

    V::Polygon p { { {0,0}, {1,0}, {1,1}, {0,1}}};

    auto res = findAntipodals( p );
    for( auto r : res ) {
        cout << r.first << ", ", << r.second << endl;
    }
}