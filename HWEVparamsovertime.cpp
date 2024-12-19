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

    // add accrued interest.
    for (int i=0; i< myPrices.size(); i++){
        if (myPrices[i] < 0.)
            continue;
        myPrices[i] += myBonds[i].getAcrruedIntrest(monitordate);
    }

    // Discounted cashflows are dirty prices
    //myRatesLagrange.setLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, priceDate);
    //myRatesNelson.setNSSYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));
    //myRatesSparse.setSparseLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));
    vector<pair<bond, num>> myBondData;
    for (int i=0; i< myBonds.size(); i++) {
        if (myPrices[i] < 0.)
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
    cout << "a=" << myHWEV.getMeanReversion(0.1) << endl;
    cout << "v=" << myHWEV.getVol(0.1) << endl;
    return pair<num,num>(myHWEV.getMeanReversion(0.1), myHWEV.getVol(0.1));

}

int main() {
    ThreadPool* mypool = &ThreadPool::getThreadPool();
    mypool->addMaxThreads();

    date priceDate(2024y/September/2d);
    date startDate(2023y/January/2d);
    date year0 = date(0,1,1);

    vector<num> mydates((priceDate-startDate).count(),0.);
    vector<num> myvols((priceDate-startDate).count(),0.);
    vector<num> myMeanReversion((priceDate-startDate).count(),0.);
    vector<future<void>> myfuts;
    auto starttime = std::chrono::high_resolution_clock::now();
    for (date monitordate = startDate; monitordate<priceDate; monitordate= monitordate + days(1)){
        int idx = (monitordate-startDate).count();
        auto myf = [monitordate,idx, year0, &mydates, &myvols,&myMeanReversion]() {

            try {
                auto wd = chrono::weekday(monitordate);
                if (wd == Saturday or wd == Sunday)
                    return;
                num myYear = year0.yearsuntil(monitordate);
                pair<num, num> myres = getHWEVvolAtDate(monitordate);
                num mya = myres.first;
                num myVol = myres.second;
                mydates[idx] = myYear;
                myvols[idx] = abs(myVol);
                myMeanReversion[idx] = abs(mya);
            } catch (const runtime_error& e) {
                cout << "Callable not run:" << monitordate << e.what() << endl;
            } catch (const invalid_argument& e){
                cout << "Callable not run:" << monitordate << e.what() << endl;
            } catch (const exception &e){
                cout << "Callable not run:" << monitordate << e.what() << endl;
            };

        };

        myfuts.push_back(mypool->spawnTask(myf));

    }

    for (auto& fut : myfuts)
        mypool->helpQueueWhile(fut);

    auto dit = std::remove(mydates.begin(), mydates.end(),0.);
    mydates.resize(distance(mydates.begin(),dit));
    auto vit = std::remove(myvols.begin(), myvols.end(),0.);
    mydates.resize(distance(myvols.begin(),vit));
    auto mit = std::remove(myMeanReversion.begin(), myMeanReversion.end(),0.);
    myMeanReversion.resize(distance(myMeanReversion.begin(),mit));


    auto endtime = std::chrono::high_resolution_clock::now();

    plot(mydates,myvols);
    //ylim({0.0,0.06});
    ylabel("\u03c3");
    xlabel("year");
    title("Hull-White Extended Vasicek Volatility");
    text(2024.3, 0.1, "Time="+ to_string(chrono::duration_cast<chrono::seconds>(endtime-starttime).count())+"s");
    text(2024.3,0.09, "Precision="+to_string(numeric_limits<num>::digits10));
    auto fig = gcf();
    //xticks({2022, 2022.5, 2023, 2023.5, 2024, 2024.5});
    xticks({2023.1,2023.4,2023.7,2024, 2024.3, 2024.6});
    array<float,3> mycolors = {1,1,1};
    fig->color(mycolors);
    ylim({0,0.11});
    if (is_same_v<num,double>)
        save("HWEVvolOverTimeDouble.png");
    else
        save("HWEVvolOverTimeSingle.png");
    ::matplot::cla();

    plot(mydates,myMeanReversion);
    ylabel("a");
    xlabel("year");
    title("Hull-White Extended Vasicek Mean Reversion");
    text(2024.3, 0.55, "Time="+ to_string(chrono::duration_cast<chrono::seconds>(endtime-starttime).count())+"s");
    text(2024.3,0.5, "Precision="+to_string(numeric_limits<num>::digits10));
    //xticks({2022, 2022.5, 2023, 2023.5, 2024, 2024.5});
    xticks({2023.1,2023.4,2023.7,2024.0, 2024.3, 2024.6});
    ylim({0,0.6});
    if (is_same_v<num,double>)
        save("HWEVaOverTimeDouble.png");
    else
        save("HWEVaOverTimeSingle.png");


}
