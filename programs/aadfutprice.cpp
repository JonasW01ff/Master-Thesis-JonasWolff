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
#include <ctime>
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

pair<num, num> getfutprice(const date& monitorDate) {
    vector<bond> myBonds = getBondData();
    vector<num> myPrices = getBondPriceAtDate(monitorDate);
    pair<double, double> swaptiondata = getSwaptionAtDate(monitorDate);
    // Start the timer
    auto start_time = std::chrono::high_resolution_clock::now();
    std::clock_t start_clock = std::clock();


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
    error = myHWEV.calibrate(myBondData, swaptionprice, swaprate, monitorDate, swaptionmaturity, swapmaturity,floatfreq, true);

    // Convert to adjoint for backprop. Seems like nodes are not being set up
    shortrateHWEV<adjointIntegral, termstructureForwardHermite<adjointIntegral>> myHWEVAdjoint(myHWEV,myBondData, monitorDate);
    termstructureForwardHermite<adjointIntegral>& myHermiteAdjoint = myHWEVAdjoint.getTermStructureRef();

    // Find the basket
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

    // Get cashflow
    num test = 10000;
    mat<num> CFnum = getCashFlowMatrix(basket, monitorDate);
    mat<adjointIntegral> CF(CFnum);
    //CF[0]->setPause();
    vector<adjointIntegral> myTenorsAdj(CF.nCols());
    vector<num> myTenors(CF.nCols());
    CF.getRow(0, myTenorsAdj);
    transform(myTenorsAdj.cbegin(), myTenorsAdj.cend(), myTenors.begin(),
              [](const adjointIntegral &adj){return adj.getVal();});
    CF.removeRow(0);
    CFnum.removeRow(0);

    // Find discount factors
    mat<adjointIntegral> df(CF.nCols(),1);
    for (int i=0;i<myTenors.size();i++)
        df[i,0] = myHermiteAdjoint.interpolateDiscount(adjointIntegral(myTenors[i]));

    // Find convetsion factors
    vector<num> conversionfactors(basket.size());
    for (int i=0; i<conversionfactors.size(); i++){
        conversionfactors[i] = basket[i].conversionfactor(fut_delivery_date,6.);
    }

    // Allocate memory for simulation
    vector<adjointIntegral> simulBONDS(CF.nRows());
    vector<num> simulBONDSnum(CF.nRows());
    cout << "Basket has size=" << simulBONDS.size() << endl;
    num simuldate = monitorDate.yearsuntil(fut_delivery_date);

    int calcerradjIts = 0;
    auto calcerradj = [&](const adjointIntegral &fut){
        calcerradjIts=0;
        num priceSum(0);
        adjointIntegral maxPrice(0.);
        adjointIntegral tempprice(0.);
        while (true) {
            calcerradjIts++;
            adjointIntegral checkpoint(0.);
            checkpoint.setPause();
            //cout << "simulate start" << endl;
            myHWEVAdjoint.getBONDsimulSobol(CF, myTenors, simuldate, simulBONDS);
            //cout << "simulate end" << endl;
            maxPrice = fut * conversionfactors[0] - simulBONDS[0];
            tempprice = 0.;
            for (int idx = 1; idx < basket.size(); idx++) {
                tempprice = fut * conversionfactors[idx] - simulBONDS[idx];
                maxPrice = (tempprice > maxPrice) ? tempprice : maxPrice;
            }
            priceSum += maxPrice.getVal();
            if (abs(test - priceSum / calcerradjIts) < numSQRTEPS){
                maxPrice.setAdjoint(1.);
                maxPrice.propagateAll();
                break;
            } else {
                maxPrice.setAdjoint(1.);
                maxPrice.propagateToPause();
                test = (priceSum / num(calcerradjIts));
            }
        }

        num err = priceSum / calcerradjIts;

        //err.setAdjoint(1.);
        //err.propagateAll();
        cout << "IDX=" << calcerradjIts << "; err=" << err << "; fut=" << fut << endl;

        return err;
    };

    auto calcerr = [&](const num fut)->num{
        int i=0;
        num priceSum = 0;
        num maxPrice = 0;
        num tempprice;
        while (true) {
            i++;
            myHWEV.getBONDsimulSobol(CFnum, myTenors, simuldate, simulBONDSnum);
            maxPrice = fut * conversionfactors[0] - simulBONDSnum[0];
            for (int idx = 1; idx < basket.size(); idx++) {
                tempprice = fut * conversionfactors[idx] - simulBONDSnum[idx];
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


    adjointIntegral finalfut(futmid);
    num finalerr = calcerradj(finalfut);
    adjointIntegral& myVol =  myHWEVAdjoint.getVolRef(0.1);
    adjointIntegral& mya = myHWEVAdjoint.getMeanReversionRef(0.1);
    num partialToFut = finalfut.getAdjoint();
    num partialToVol = myVol.getAdjoint();
    num partialToa = mya.getAdjoint();
    auto& mySpline = myHermiteAdjoint.getSplineRef();
    auto mySplineParams = mySpline.getParams();
    mat<num> dFutdY(1,3*4,0.);
    cout << "dFut/dVol=" << -partialToVol/partialToFut << endl;
    cout << "dFut/dBach=" << -partialToVol/partialToFut*myHWEV.getdVoldBach() << endl;
    cout << "dFut/da=" << -partialToa/partialToFut << endl;
    for (int i =0; i<3; i++){
        for (int j=1;j<5;j++){
            dFutdY[0,i*4+(j-1)] = -mySplineParams[i][j].getAdjoint()/partialToFut;
        }
    }
    cout << "dFut/da_0=" << -mySplineParams[0][1].getAdjoint()/partialToFut << endl;
    cout << "dFut/db_0=" << -mySplineParams[0][2].getAdjoint()/partialToFut << endl;
    cout << "dFut/dc_0=" << -mySplineParams[0][3].getAdjoint()/partialToFut << endl;
    cout << "dFut/dd_0=" << -mySplineParams[0][4].getAdjoint()/partialToFut << endl;
    cout << "dFut/da_1=" << -mySplineParams[1][1].getAdjoint()/partialToFut << endl;
    cout << "dFut/db_1=" << -mySplineParams[1][2].getAdjoint()/partialToFut << endl;
    cout << "dFut/dc_1=" << -mySplineParams[1][3].getAdjoint()/partialToFut << endl;
    cout << "dFut/dd_1=" << -mySplineParams[1][4].getAdjoint()/partialToFut << endl;
    cout << "dFut/da_2=" << -mySplineParams[2][1].getAdjoint()/partialToFut << endl;
    cout << "dFut/db_2=" << -mySplineParams[2][2].getAdjoint()/partialToFut << endl;
    cout << "dFut/dc_2=" << -mySplineParams[2][3].getAdjoint()/partialToFut << endl;
    cout << "dFut/dd_2=" << -mySplineParams[2][4].getAdjoint()/partialToFut << endl;

    // Market Sensitivity
    const auto& TSparamsAdjErr = myHermiteAdjoint.getParamsAdjErr();
    const auto& TSBondPriceAdjErr = myHermiteAdjoint.getBondAdjErr();
    mat<num> dYdBond(3*4,TSBondPriceAdjErr.size(), 0);
    for (int i=0; i< TSBondPriceAdjErr.size();i++) {
        dYdBond[0,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][0];
        dYdBond[1,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][1];
        dYdBond[2,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][2];
        dYdBond[3,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][3];

        dYdBond[4,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][0];
        dYdBond[5,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][1];
        dYdBond[6,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][2];
        dYdBond[7,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][3];

        dYdBond[8,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][0];
        dYdBond[9,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][1];
        dYdBond[10,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][2];
        dYdBond[11,i] = -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][3];

        cout << "da_0/" << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][0] << endl;
        cout << "db_0/"  << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][1] << endl;
        cout << "dc_0/"  << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][2] << endl;
        cout << "dd_0/"  << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[0][3] << endl;

        cout << "da_1/" << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[1][0] << endl;
        cout << "db_1/"  << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[1][1] << endl;
        cout << "dc_1/"  << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[1][2] << endl;
        cout << "dd_1/"  << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[1][3] << endl;

        cout << "da_2/"  << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[2][0] << endl;
        cout << "db_2/"  << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[2][1] << endl;
        cout << "dc_2/"  << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[2][2] << endl;
        cout << "dd_2/"  << "dBond_" << i << "=" << -TSBondPriceAdjErr[i] / TSparamsAdjErr[2][3] << endl;
    }

    mat<num> dFutDBond = dFutdY.matmul(move(dYdBond));
    for (int i=0; i< dFutDBond.size();i++)
        cout << "dFut/dBond_" << i<<"=" << dFutDBond[0,i] << endl;
    //dFutDBond.print();

    // End the timer
    std::clock_t end_clock = std::clock();
    auto end_time = std::chrono::high_resolution_clock::now();

    // Calculate the elapsed time in milliseconds
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    cout<< "Fut Price " << futmid << " with error " << errmid <<" date="<< monitorDate <<  endl;
    std::cout << "The program took " << elapsed_time << " milliseconds to execute." << std::endl;
    std::cout << "The program took " << int(num(end_clock-start_clock)/CLOCKS_PER_SEC*1000) << " milliseconds CPU time." << std::endl;
    return make_pair(futmid, elapsed_time);




}

int main() {
    ThreadPool* mypool = &ThreadPool::getThreadPool();
    //mypool->addMaxThreads();

    date priceDate(2024y/September/4d);
    date startDate(2024y/September/2d);
    date year0 = date(0,1,1);

    vector<num> mydates((priceDate-startDate).count(),0.);
    vector<num> myFutPrices((priceDate-startDate).count(),0.);
    vector<future<void>> myfuts;
    vector<num> myruntime((priceDate-startDate).count(),0.);
    for (date monitordate = startDate; monitordate<priceDate; monitordate= monitordate + days(1)){
        int idx = (monitordate-startDate).count();
        auto myf = [monitordate,idx, year0, &mydates, &myFutPrices, &myruntime]() {

            try {
                num myYear = year0.yearsuntil(monitordate);

                auto starttime = std::chrono::high_resolution_clock::now();
                pair<num,num> myres = getfutprice(monitordate);
                mydates[idx] = myYear;
                myFutPrices[idx] = myres.first;
                myruntime[idx] = myres.second;
                auto endtime = std::chrono::high_resolution_clock::now();
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




}
