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
    cout << "a=" << myHWEV.getMeanReversion(0.1) << endl;
    cout << "v=" << myHWEV.getVol(0.1) << endl;
    return pair<num,num>(myHWEV.getMeanReversion(0.1), myHWEV.getVol(0.1));

}

pair<num, num> getfutprice(termstructureForwardHermite<num>& oldTS, bool newTS, const date& monitorDate, const num volbump, const num abump, const num splinebump,
                           const size_t knotidx, const size_t knotsubidx, const size_t bondidx=0, const num bondpricebump=0) {
    vector<bond> myBonds = getBondData();
    vector<num> myPrices = getBondPriceAtDate(monitorDate);
    pair<double, double> swaptiondata = getSwaptionAtDate(monitorDate);
    // Start the timer
    auto start_time = std::chrono::high_resolution_clock::now();

    // add accrued interest.
    for (int i=0; i< myPrices.size(); i++){
        myPrices[i] += myBonds[i].getAcrruedIntrest(monitorDate);
        if (i==bondidx) {
            //cout << "Bond bumped by " << bondpricebump << endl;
            myPrices[i] += bondpricebump;
        }
    }

    // Find delivery date
    date fut_delivery_date(monitorDate.year(), monitorDate.month(), 10d);
    while (fut_delivery_date.month() != March and fut_delivery_date.month() != June and fut_delivery_date.month() != September and fut_delivery_date.month() != December)
        fut_delivery_date += months(1);

    // Discounted cashflows are dirty prices
    //myRatesLagrange.setLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, priceDate);
    //myRatesNelson.setNSSYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));
    //myRatesSparse.setSparseLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));
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

    if(newTS)
        oldTS = termstructureForwardHermite<num>("FRG",myBondData, monitorDate);
    auto end_hermite_time = std::chrono::high_resolution_clock::now();
    /*if (!myHermite.bumpParam(splinebump, knotidx, knotsubidx))
        throw runtime_error("Cant bump spline param");*/
    shortrateHWEV<num, termstructureForwardHermite<num>> myHWEV(oldTS);
    //cout << "CALLIBRATING" << endl;
    bool useAAD = false;
    num error;
    const date swaptionmaturity = monitorDate + years(1) +days(2);
    const date swapmaturity = monitorDate + years(11) + days(2);
    int floatfreq = 1;
    num swaptionprice = swaptiondata.first;
    num swaprate = swaptiondata.second;
    static num sigma(-1);
    static num a;
    if (useAAD)
        error = myHWEV.calibrateAAD(myBondData,monitorDate);
    else {
        //auto temp_time = std::chrono::high_resolution_clock::now();
        if (sigma==-1 or volbump>0){
            error = myHWEV.calibrate(myBondData, swaptionprice, swaprate, monitorDate, swaptionmaturity, swapmaturity,floatfreq, volbump>0);
            sigma = myHWEV.getVol(.1);
            a = myHWEV.getMeanReversion(0.1);
        } else{
            myHWEV.setVol(sigma);
            myHWEV.setMeanReversion(a);
        }

    }

    myHWEV.setVol(myHWEV.getVol(.1)+volbump);
    myHWEV.setMeanReversion(myHWEV.getMeanReversion(.1)+abump);
    termstructureForwardHermite<num>& myBumpHermite = myHWEV.getTermStructureRef();
    if (!myBumpHermite.bumpParam(splinebump, knotidx, knotsubidx))
        throw runtime_error("Cant bump spline param");

    vector<bond> basket;
    for (auto myBond : myBonds){
        bool inc = true;
        num maturity = monitorDate.yearsuntil(myBond.getMaturityDate());
        num startmaturity = myBond.getIssuanceDate().yearsuntil(myBond.getMaturityDate());
        if (startmaturity>11.)
            inc = false;
        if (maturity>10.5 or maturity<8.5)
            inc = false;
        if (inc)
            basket.push_back(myBond);
    }

    num test = 10000;
    mat<num> CF = getCashFlowMatrix(basket, monitorDate);
    vector<num> myTenors(CF.nCols());
    CF.getRow(0, myTenors);
    CF.removeRow(0);
    mat<num> df(CF.nCols(),1);
    for (int i=0;i<myTenors.size();i++)
        df[i,0] = myBumpHermite.interpolateDiscount(myTenors[i]);
    vector<num> conversionfactors(basket.size());
    for (int i=0; i<conversionfactors.size(); i++){
        conversionfactors[i] = basket[i].conversionfactor(fut_delivery_date,6.);
    }
    vector<num> simulBONDS(CF.nRows());
    //cout << "Basket has size=" << simulBONDS.size() << endl;
    num simuldate = monitorDate.yearsuntil(fut_delivery_date);
    num priceSum = 0;
    num maxPrice;
    num tempprice;

    auto calcerr = [&](const num fut){
        int i=0;
        priceSum = 0;
        while (true) {
            i++;
            myHWEV.getBONDsimulSobol(CF, myTenors, simuldate, simulBONDS);
            maxPrice = fut * conversionfactors[0] - simulBONDS[0];
            for (int idx = 1; idx < basket.size(); idx++) {
                tempprice = fut * conversionfactors[idx] - simulBONDS[idx];
                maxPrice = (tempprice > maxPrice) ? tempprice : maxPrice;
            }
            priceSum += maxPrice;
            if (abs(test - priceSum / i) < numSQRTEPS)
                break;
            else
                test = priceSum / i;
        }
        //cout << "IDX=" << i << "; err=" << priceSum / i << "; fut=" << fut << endl;
        return priceSum/i;
    };

    num futlow = 1;
    num errlow = calcerr(futlow);
    num futhigh = 300;
    num errhigh = calcerr(futhigh);
    num err = 1000;
    num futmid = (futlow+futhigh)/2.;
    num errmid;
    while (abs(futhigh-futlow) > numSQRTEPS) {
        futmid = (futlow+futhigh)/2.;
        errmid = calcerr(futmid);
        if (errmid > 0){
            futhigh = futmid;
            errhigh = errmid;
        } else if (errmid<0) {
            futlow = futmid;
            errlow = errmid;
        } else {
            cout << "JACKPOT" << endl;
            break;
        }

        if (abs(errhigh)/(errhigh) == abs(errlow)/(errlow))
            throw runtime_error("???");

    }

    // End the timer
    auto end_time = std::chrono::high_resolution_clock::now();

    // Calculate the elapsed time in milliseconds
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();


    //cout<< "Fut Price " << futmid << " with error " << errmid <<" date="<< monitorDate <<  endl;
    if (volbump)
        return make_pair(futmid, myHWEV.getdVoldBach());
    return make_pair(futmid, elapsed_time);


    std::cout << "The program took " << elapsed_time << " milliseconds to execute." << std::endl;


}

int main() {
    ThreadPool* mypool = &ThreadPool::getThreadPool();
    mypool->addMaxThreads();

    date priceDate(2024y/September/3d);
    date startDate(2024y/September/2d);
    date year0 = date(0,1,1);

    vector<num> mydates((priceDate-startDate).count(),0.);
    vector<num> myFutPrices((priceDate-startDate).count(),0.);
    vector<future<void>> myfuts;
    vector<num> myruntime((priceDate-startDate).count(),0.);
    vector<num> myknots = termstructureForwardHermite<num>("Dummy").getTimePoints();
    size_t nknots = myknots.size();
    auto starttime = std::chrono::high_resolution_clock::now();
    for (date monitordate = startDate; monitordate<priceDate; monitordate= monitordate + days(1)){
        int idx = (monitordate-startDate).count();
        auto myf = [monitordate,idx, year0, nknots, &mydates, &myFutPrices, &myruntime]() {

            try {
                cout << "begin" << endl;
                termstructureForwardHermite<num> myTS("Dummy");
                num myYear = year0.yearsuntil(monitordate);
                auto starttime = std::chrono::high_resolution_clock::now();
                pair<num,num> myres = getfutprice(myTS, true, monitordate,0,0, 0, 0, 0);
                pair<num,num> myresvolbum = getfutprice(myTS, false, monitordate,numSQRTSQRTEPS,0, 0, 0, 0);
                pair<num,num> myresabum = getfutprice(myTS, false, monitordate,0,numSQRTSQRTEPS, 0,0,0);
                cout << "dFut/dVol=" << (myresvolbum.first- myres.first)/(numSQRTSQRTEPS) << endl;
                cout << "dFut/dBach=" <<(myresvolbum.first- myres.first)/(numSQRTSQRTEPS)*myresvolbum.second << endl;
                cout << "dFut/da=" << (myresabum.first- myres.first)/numSQRTSQRTEPS << endl;
                for (int i =0; i< nknots-1;i++){
                    pair<num, num> myressplinebump1 = getfutprice(myTS, false, monitordate,0,0, numSQRTSQRTEPS, i, 1);
                    cout << "dFut/da_"<< i << "=" << (myressplinebump1.first- myres.first)/numSQRTSQRTEPS << endl;
                    pair<num, num> myressplinebump2 = getfutprice(myTS, false,monitordate,0,0, numSQRTSQRTEPS, i, 2);
                    cout << "dFut/db_"<< i << "=" << (myressplinebump2.first- myres.first)/numSQRTSQRTEPS << endl;
                    pair<num, num> myressplinebump3 = getfutprice(myTS, false,monitordate,0,0, numSQRTSQRTEPS, i, 3);
                    cout << "dFut/dc_"<< i << "=" << (myressplinebump3.first- myres.first)/numSQRTSQRTEPS << endl;
                    pair<num, num> myressplinebump4= getfutprice(myTS, false,monitordate,0,0, numSQRTSQRTEPS, i, 4);
                    cout << "dFut/dd_"<< i << "=" << (myressplinebump4.first- myres.first)/numSQRTSQRTEPS << endl;
                }
                /*for (int i=0; i<46; i++){
                    pair<num, num> myresbondbump = getfutprice(monitordate,0,0, 0., 0., 0., i, numSQRTSQRTEPS);
                    cout << "dFut/dbond_"<< i << "=" << (myresbondbump.first- myres.first)/numSQRTSQRTEPS << endl;
                }*/

                auto endtime = std::chrono::high_resolution_clock::now();
                cout << "Risk Report Time=" << chrono::duration_cast<chrono::milliseconds>(endtime-starttime) << endl;
                mydates[idx] = myYear;
                myFutPrices[idx] = myres.first;
                myruntime[idx] = myres.second;
            } catch (exception &e) {
                cout << monitordate << e.what() << endl;
            };

        };

        myfuts.push_back(mypool->spawnTask(myf));

    }

    for (auto& fut : myfuts)
        mypool->helpQueueWhile(fut);



    // total data read time
    num mytotalruntime = reduce(myruntime.begin(), myruntime.end())/myruntime.size();

    auto dit = std::remove(mydates.begin(), mydates.end(),0.);
    mydates.resize(distance(mydates.begin(),dit));
    auto vit = std::remove(myFutPrices.begin(), myFutPrices.end(),0.);
    mydates.resize(distance(myFutPrices.begin(),vit));


    auto endtime = std::chrono::high_resolution_clock::now();
    return 0;
    plot(mydates,myFutPrices);
    ylim({125.,140.});
    ylabel("€");
    xlabel("year");
    title("Hull-White Extended Vasicek Bund-Future Price");
    text(2024.41, 139, "Time="+ to_string(chrono::duration_cast<chrono::seconds>(endtime-starttime).count())+"s");
    //text(2024.45, 138.5, "Read Time="+ to_string(int(datareadtime))+"s");
    text(2024.41, 138, "Avg Calc="+ to_string(int(mytotalruntime))+"ms");
    text(2024.41,137, "Precision="+to_string(numeric_limits<num>::digits10));
    auto fig = gcf();
    //xticks({2022, 2022.5, 2023, 2023.5, 2024, 2024.5});
    xticks({2023.7, 2023.9,2024.1, 2024.3, 2024.5});
    array<float,3> mycolors = {1,1,1};
    fig->color(mycolors);
    //ylim({0,0.11});
    if (is_same_v<num,double>)
        save("HWEVpriceOverTimeDouble.png");
    else
        save("HWEVpriceOverTimeSingle.png");
    ::matplot::cla();


}
