//
// Created by Jonas Wolff on 29/10/2024.
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

int main(){

    mat<double> x(2,1);
    mat<double> xsol(2,1);
    x[0][0] = 2.0;
    x[1][0] = 2.0;
    mat<double> jacobian(2,1);
    mat<double> direction (2,1);
    direction[0,0] = -1.;
    direction[1,0] = -1.;
    double y = myfunc(x);
    mat<double> yMat(1,1);
    double ysol;
    yMat[0,0] = y;
    func_fd<double, double>(
            [](const mat<double> &x, mat<double> &y){y[0,0] = myfunc(x);},
            x,
            yMat,
            jacobian);
    bool check = false;
    lnsrch<double>(x, y, jacobian, direction, xsol, ysol,0.1, check, myfunc);
    jacobian.print();
    cout << jacobian.getIsTransposed() << endl;

    return 0;

}