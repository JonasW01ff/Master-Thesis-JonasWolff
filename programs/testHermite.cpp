#include <iostream>

#include <string>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string_view>
#include <random>
#include "bond.h"
#include "date.h"
#include "rates.h"
#include "interpolation.h"

using namespace std;
using namespace chrono;


// data fetching
vector<bond> getBondData(){

}


// data fetching
vector<double> getBondPriceToday(){

}


int main() {

    vector<double> xVec = {0,1,2,3,4};
    vector<double> yVec = {1,2,1.5,1.7,5};
    vector<double> sVec = {1,1,1,1,1};

    hermiteMCspline<double> myHermite;
    myHermite.calibrate(xVec.cbegin(), xVec.cend(), yVec.cbegin(),
                        sVec.begin());

    for (double x = -3.; x<8.1; x+=0.1)
        cout << to_string(myHermite.getNode(x)) << endl;

    return 0;
}
