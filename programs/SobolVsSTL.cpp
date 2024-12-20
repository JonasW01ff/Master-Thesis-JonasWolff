//
// Created by Jonas Wolff on 23/11/2024.
//

#include "../sobol.h"
#include "matplot/matplot.h"
#include <random>

using namespace matplot;

int main(){
    sobol1D rng = sobol1D();
    vector<double> myU1(1000);
    vector<double> myU2(1000);
    for (int i=0; i<myU1.size(); i++) {
        myU1[i] = rng.getU<0>();
        myU2[i] = rng.getU<1>();
    }
    subplot(1,2,1);
    title("Sobol: N=1000");
    xlabel("x");
    auto h1 = scatter(myU1,myU2);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::random_device rd2;
    std::mt19937 gen2(rd2());
    uniform_real_distribution<> dis(0,1);
    for (int i=0; i<myU1.size(); i++) {
        myU1[i] = dis(gen);
        myU2[i] = dis(gen2);
    }
    subplot(1,2,2);
    title("STL Mersenne Twister: N=1000");
    xlabel("x");
    ylabel("y");
    auto h2 = scatter(myU1,myU2);

    auto fig = gcf();
    array<float,3> mycolors = {1,1,1};
    fig->color(mycolors);
    fig->size(800, 400);

    show();
}