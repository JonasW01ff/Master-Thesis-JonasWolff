//
// Created by Jonas Wolff on 31/10/2024.
//
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
#include "numbers.h"
#include "optimizer.h"

int main() {
    vector <dual<double>> xVec = {0, 1, 2, 3, 4};
    vector <dual<double>> yVec = {1, 2, {1.5, 1}, 1.7, 5};
    vector <dual<double>> sVec = {1, 1, 1, 1, 1};
    hermiteMCspline <dual<double>> myHermite;
    myHermite.calibrate(xVec.cbegin(), xVec.cend(), yVec.cbegin(),
                        sVec.begin());

    dual<double> myDual(4., 1.);
    cout << pow(myDual, 2.) * exp(4. * abs(myDual)) << endl;

    for (double x = -3.; x < 8.1; x += 0.1)
        cout << to_string(myHermite.getNode(x)) << endl;

    return 0;
}