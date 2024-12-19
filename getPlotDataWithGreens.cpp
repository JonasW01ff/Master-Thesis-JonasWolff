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

    rates myRates;
    vector<bond> myBonds = getBondData();
    vector<double> myPrices = getBondPriceToday();
    date priceDate(2024y/September/2d);
    cout << double(priceDate)<< endl;
    // Train on clean prices since it tries to discount future payments
    myRates.setSparseLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, priceDate);
    //myRates.setNSSYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));


    // Set intial parameters
    //mat initParams(6,1,0.5);
    /*initParams.setEntry(1,0,3.15);
    initParams.setEntry(2,0,-1.13);
    initParams.setEntry(3,0,6.28);
    initParams.setEntry(4,0,1.2);
    initParams.setEntry(5,0,12.7);
    myRates.setNelsonSiegelParams(initParams);*/

    // add accrued interest.
    for (int i=0; i< myPrices.size(); i++){
        myPrices[i] += myBonds[i].getAcrruedIntrest(priceDate);

    }

    for (int i=0; i< myPrices.size(); i++){
        double myPrice = myPrices[i];
        double myEst = myRates.getBondPrice(myBonds[i], priceDate);
        string sPrice = to_string(myPrice);
        replace(sPrice.begin(), sPrice.end(), '.', ',');
        string sEst = to_string(myEst);
        replace(sEst.begin(), sEst.end(), '.', ',');
        cout << myBonds[i].getMaturityDate() << ";" << sPrice << ";" << sEst << endl;
    }
    return 0;
}
