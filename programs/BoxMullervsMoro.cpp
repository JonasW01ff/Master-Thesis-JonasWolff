//
// Created by Jonas Wolff on 23/11/2024.
//
#include "../sobol.h"
#include "matplot/matplot.h"
#include <chrono>

using namespace matplot;
using namespace chrono;

int main(){
    sobol1D rng = sobol1D();
    vector<double> myNBM(100000000);
    auto t_start = high_resolution_clock::now();
    for (int i=0; i<myNBM.size(); i++)
        myNBM[i] = rng.getNBoxMuller();
    auto t_end = high_resolution_clock::now();
    double box_muller_time = duration<double>(t_end - t_start).count();
    vector<double> myNBMNOSIMD(100000000);
    t_start = high_resolution_clock::now();
    for (int i=0; i<myNBM.size(); i++)
        myNBMNOSIMD[i] = rng.getNBoxMullerNOSIMD();
    t_end = high_resolution_clock::now();
    double box_muller_timeNOSIMD = duration<double>(t_end - t_start).count();
    vector<double> myNMoro(100000000);
    t_start = high_resolution_clock::now();
    for (int i=0; i<myNMoro.size(); i++)
        myNMoro[i] = rng.getNMoro();
    t_end = high_resolution_clock::now();
    double moro_time = duration<double>(t_end - t_start).count();

    cout <<"Box-Muller SIMD: " << box_muller_time << endl;
    cout <<"Box-Muller NO SIMD: " << box_muller_timeNOSIMD << endl;
    cout <<"Moro: " << moro_time << endl;
    subplot(1,2,1);
    title("Box-Muller");
    xlabel("x");
    auto h1 = matplot::hist(myNBM,100);
    text(0.0, 6e6*0.95, "RNG: Sobol");
    text(0.0, 6e6*0.9, "Simuls:100000000");
    text(0.0, 6e6*0.85, "Time:"+ to_string(box_muller_time));
    subplot(1,2,2);
    title("Moro");
    xlabel("x");
    ylabel("y");
    auto h2 = matplot::hist(myNMoro,100);
    text(0.0, 6e6*0.95, "RNG: Sobol");
    text(0.0, 6e6*0.9, "Simuls:100000000");
    text(0.0, 6e6*0.85, "Time:"+ to_string(moro_time));

    auto fig = gcf();
    array<float,3> mycolors = {1,1,1};
    fig->color(mycolors);
    fig->size(800, 400);
    auto h3 = matplot::hist(myNBMNOSIMD,100);
    show();
    return 1;
}