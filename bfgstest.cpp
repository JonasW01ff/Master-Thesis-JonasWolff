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

double myfunc(const mat<double> &x){

    mat<double> myX = x;
    return cos(x[0,0]) * sin(x[1,0]*3-cos(x[0,0]*x[1,0]));
}

int main(){
    mat<double> x(2,1);
    mat<double> xsol(2,1);
    x[0][0] = 1.0;
    x[1][0] = 6.33;
    double ysol;
    bfgs<double> mySolver;
    function<double(const mat<double>&)> myWrppedFunc = myfunc;
    mySolver.minimize(myWrppedFunc, x, 1e-6);
    x.print(); //cout << x << endl;
    cout << myfunc(x) << endl;
    return 0;
}