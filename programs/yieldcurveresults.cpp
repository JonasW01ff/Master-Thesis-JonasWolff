//
// Created by Jonas Wolff on 12/12/2024.
//
//
// Created by Jonas Wolff on 09/12/2024.
//
#include <iostream>
#include "settings.h"
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
#include "termstructureHermite.h"
#include "termstructureForwardHermite.h"
#include "adjointNumeric.h"
#include "quadrature.h"
#include "shortrateHWEV.h"
#include "sobol.h"
#include <matplot/matplot.h>
#include "plot.h"
#include <random>
#include "ThreadPool.h"
#include <string>
#include "ois.h"

#ifndef GNUPLOT_EXECUTABLE
#define GNUPLOT_EXECUTABLE "gnuplot" // Default fallback
#endif

using namespace std;
using namespace chrono;

// data fetching
vector<bond> getBondData(){

}

// data fetching
vector<num> getBondPriceToday(){

}

// data fetching
std::string fetchPrice(const std::string& filePath, const std::string& year, const std::string& month, const std::string& day, const std::string& isin) {

}

// Swaption 1Y10Y, Swaprate today. data fetching
pair<double, double> getSwaptionAtDate(const date &myDate){

}


// data fetching
vector<num> getBondPriceAtDate(const date &myDate){

}

using namespace matplot;

pair<num,num> getHWEVvolAtDate(const date &monitordate){

    vector<bond> myBonds = getBondData();
    vector<num> myPrices = getBondPriceAtDate(monitordate);
    cout << num(monitordate)<< endl;

    // add accrued interest.
    for (int i=0; i< myPrices.size(); i++){
        myPrices[i] += myBonds[i].getAcrruedIntrest(monitordate);
    }

    // Discounted cashflows are dirty prices
    //myRatesLagrange.setLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, priceDate);
    //myRatesNelson.setNSSYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));
    //myRatesSparse.setSparseLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));
    vector<pair<bond, num>> myBondData;
    for (int i=0; i< myBonds.size(); i++) {
        if (myPrices[i] < 10.)
            continue;
        myBondData.push_back(make_pair(myBonds[i], myPrices[i]));
    }
    if (myBondData.size() < 10)
        throw runtime_error("Too few bonds");

    // Start the timer
    auto start_time = std::chrono::high_resolution_clock::now();
    termstructureForwardHermite<num> myHermite("FRG",myBondData, monitordate);
    auto end_hermite_time = std::chrono::high_resolution_clock::now();


    shortrateHWEV<num, termstructureForwardHermite<num>> myHWEV(myHermite);
    cout << "CALLIBRATING" << endl;
    bool useAAD = false;
    num error;
    pair<double, double> swaptiondata = getSwaptionAtDate(monitordate);
    const date swaptionmaturity = monitordate + years(1) +days(2);
    const date swapmaturity = monitordate + years(11) + days(2);
    int floatfreq = 1;
    num swaptionprice = swaptiondata.first;
    num swaprate = swaptiondata.second;
    if (useAAD)
        error = myHWEV.calibrateAAD(myBondData,monitordate);
    else {
        error = myHWEV.calibrate(myBondData, swaptionprice, swaprate, monitordate, swaptionmaturity, swapmaturity,floatfreq);
    }
    num myVol = myHWEV.getVol(0.1);
    num mya = myHWEV.getMeanReversion(0.1);
    cout << "a=" << mya << endl;
    cout << "v=" << myVol << endl;
    return pair<num,num>(mya, myVol);

}

void plotYielCurve (const date& monitorDate) {
    vector<bond> myBonds = getBondData();
    vector<num> myPrices = getBondPriceAtDate(monitorDate);
    pair<double, double> swaptiondata = getSwaptionAtDate(monitorDate);
    // Start the timer
    auto start_time = std::chrono::high_resolution_clock::now();


    cout << num(monitorDate)<< endl;

    // add accrued interest.
    for (int i=0; i< myPrices.size(); i++){
        myPrices[i] += myBonds[i].getAcrruedIntrest(monitorDate);
    }

    // Find delivery date
    date fut_delivery_date(monitorDate.year(), monitorDate.month(), 10d);
    while (fut_delivery_date.month() != March and fut_delivery_date.month() != June and fut_delivery_date.month() != September and fut_delivery_date.month() != December)
        fut_delivery_date += months(1);

    // Discounted cashflows are dirty prices
    vector<pair<bond, num>> myBondData;
    for (int i=0; i< myBonds.size(); i++) {
        if (myPrices[i] < 10.)
            continue;
        // Check if bond is matures at delivery date
        if (myBonds[i].getMaturityDate() < fut_delivery_date +  months(1))
            continue;
        myBondData.push_back(make_pair(myBonds[i], myPrices[i]));
    }
    if (myBondData.size() < 10)
        throw runtime_error("Too few bonds");


    // Calibrate yield curve
    termstructureForwardHermite<num> myHermite("FRG",myBondData, monitorDate);
    auto end_hermite_time = std::chrono::high_resolution_clock::now();

    // Calibrate short rate model
    shortrateHWEV<num, termstructureForwardHermite<num>> myHWEV(myHermite);
    cout << "CALLIBRATING" << endl;
    bool useAAD = false;
    num error;
    const date swaptionmaturity = monitorDate + years(1) +days(2);
    const date swapmaturity = monitorDate + years(11) + days(2);
    int floatfreq = 1;
    num swaptionprice = swaptiondata.first;
    num swaprate = swaptiondata.second;
    error = myHWEV.calibrate(myBondData, swaptionprice, swaprate, monitorDate, swaptionmaturity, swapmaturity,floatfreq);

    // Plot bond price errors
    vector<num> hermitePrice(myPrices.size());
    vector<int> bondIDX(myPrices.size());
    num avgErr = 0;
    for (int i=0; i< myPrices.size(); i++){
        bondIDX[i] = i;
        hermitePrice[i] = myHermite.getBondPrice(myBonds[i], monitorDate);
        avgErr += abs(hermitePrice[i] - myPrices[i]);
    }
    avgErr /= myPrices.size();
    ::matplot::cla();
    plot(bondIDX, hermitePrice, bondIDX, myPrices);
    xlabel("Bond IDX");
    ylabel("price");
    title("Hermite Bond Pricing");
    ::matplot::legend({"Hermite", "Market"});
    save("bondPricingHermite.png");
    ::matplot::cla();

    // Print yield curve
    vector<num> yield(0);
    vector<num> forward(0);
    vector<num> forwardT(0);
    vector<num> tenors;
    tenors.clear();
    for (int i=0; i< myPrices.size(); i++){
        bondIDX[i] = i;
        hermitePrice[i] = myHermite.getBondPrice(myBonds[i], monitorDate);
    }
    cout << "YIELD CURVE" << endl;
    cout << "Date;Yield;Forward;ForwardDerivative" << endl;
    for (num i= 0.01; i< 30; i+=0.01){
        tenors.push_back(i);
        yield.push_back(myHermite.interpolateYield(i));
        forward.push_back(myHermite.interpolateForward(i));
        forwardT.push_back(myHermite.interpolateForwardDerivative(i));
    }
    plot(tenors, yield, tenors, forward, tenors, forwardT);
    xlabel("tenor");
    ylabel("%");
    title("Hermite Yield Curve");
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_hermite_time - start_time).count();
    text(24, -0.015, "Time: " + to_string(elapsed_time) + "ms");
    text(24, -0.02, "Error: " + to_string(avgErr).substr(0,6));
    text(24, -0.025, "Precision: "+to_string(numeric_limits<num>::digits10));
    ylim({-0.03,0.05});
    ::matplot::legend({"Yield", "Forward", "Forward Derivative"});
    save("hermiteYieldCurve"+to_string(numeric_limits<num>::digits10)+".png");
}

int main() {
    ThreadPool* mypool = &ThreadPool::getThreadPool();
    //mypool->addMaxThreads();

    date priceDate(2024y/September/4d);
    date startDate(2024y/September/2d);
    date year0 = date(0,1,1);

    plotYielCurve(startDate);

    auto endtime = std::chrono::high_resolution_clock::now();




}
