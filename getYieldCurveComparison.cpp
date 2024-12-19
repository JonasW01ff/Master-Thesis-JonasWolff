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

using namespace std;
using namespace chrono;

// data fetching
vector<bond> getBondData(){

}


// data fetching
vector<double> getBondPriceToday(){

}


int main() {

    rates myRatesLagrange, myRatesNelson, myRatesSparse;
    vector<bond> myBonds = getBondData();
    vector<double> myPrices = getBondPriceToday();
    date priceDate(2024y/September/2d);
    cout << double(priceDate)<< endl;

    // add accrued interest.
    for (int i=0; i< myPrices.size(); i++){
        myPrices[i] += myBonds[i].getAcrruedIntrest(priceDate);
    }

    // Discounted cashflows are dirty prices
    myRatesLagrange.setLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, priceDate);
    myRatesNelson.setNSSYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));
    myRatesSparse.setSparseLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));

    // Set intial parameters
    //mat initParams(6,1,0.5);
    /*initParams.setEntry(1,0,3.15);
    initParams.setEntry(2,0,-1.13);
    initParams.setEntry(3,0,6.28);
    initParams.setEntry(4,0,1.2);
    initParams.setEntry(5,0,12.7);
    myRates.setNelsonSiegelParams(initParams);*/

    // Print bond price errors
    cout << "BOND ERRORS" << endl;
    cout << "Date;Price;Lagrange;Nelson;Sparse" << endl;
    for (int i=0; i< myPrices.size(); i++){
        double myPrice = myPrices[i];
        double myEst = myRatesLagrange.getBondPrice(myBonds[i], priceDate);
        double myEstNSS = myRatesNelson.getBondPrice(myBonds[i], priceDate);
        double myEstSparse = myRatesSparse.getBondPrice(myBonds[i], priceDate);
        string sPrice = to_string(myPrice);
        replace(sPrice.begin(), sPrice.end(), '.', ',');
        string sEst = to_string(myEst), sEstNSS = to_string(myEstNSS),
            sEstSparse = to_string(myEstSparse);
        replace(sEst.begin(), sEst.end(), '.', ',');
        replace(sEstNSS.begin(), sEstNSS.end(), '.', ',');
        replace(sEstSparse.begin(), sEstSparse.end(), '.', ',');

        cout << myBonds[i].getMaturityDate() << ";" << sPrice << ";" << sEst << ";" << sEstNSS << ";" << sEstSparse << endl;
    }

    cout << endl;

    // Print yield curve
    cout << "YIELD CURVE" << endl;
    cout << "Date;Lagrange;Nelson;Sparse" << endl;
    for (int i= 1; i< 30*4; i++){
        double lagrY = myRatesLagrange.getYieldCurve(.25*i);
        double nelY = myRatesNelson.getYieldCurve(.25*i);
        double sparY = myRatesSparse.getYieldCurve(.25*i);
        string slagrY = to_string(lagrY), snelY = to_string(nelY),
                ssparY = to_string(sparY);
        replace(slagrY.begin(), slagrY.end(), '.', ',');
        replace(snelY.begin(), snelY.end(), '.', ',');
        replace(ssparY.begin(), ssparY.end(), '.', ',');
        string myMaturity = to_string(.25*i);
        replace(myMaturity.begin(), myMaturity.end(), '.', ',');
        cout  << myMaturity << ";" << slagrY << ";";
        cout << snelY << ";";
        cout <<ssparY << ";";
        cout << endl;

    }

    return 0;
}
