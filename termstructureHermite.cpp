//
// Created by Jonas Wolff on 31/10/2024.
//

#include "termstructureHermite.h"

template<class T>
termstructureHermite<T>::termstructureHermite(const std::string &name, vector<pair<bond, T>> &bondData, date priceDate)
            : termstructureStore<T, termstructureHermite<T>>(name) {

    /*for (int i = 1; i < 10; i+=1)
        timePoints.push_back(i);
    for (int i = 15; i<31; i+=1)
        timePoints.push_back(i);*/

    function<T(const mat<T>&)> errFunc = [&bondData, &priceDate, this]
            (const mat<T> &myTermStruct){

        // Check dimensions
        if (timepoints.size() != myTermStruct.size()/3.)
            throw runtime_error("The provided term structure does not follow the provided time line");

        // Clear last parameters
        termstructureStore<T, termstructureHermite<T>>::clearAllParameters();

        // Set yields and forwards
        for (int i = 0; i<timepoints.size() ; i++){

            T myTimePoint = myTermStruct[timepoints.size()*2+i,0];
            //myTimePoint = (myTimePoint < 1.) ?  1. : myTimePoint;
            termstructureStore<T, termstructureHermite<T>>::addYield(myTimePoint, myTermStruct[i,0]);
            termstructureStore<T, termstructureHermite<T>>::addForward(T(0.),myTimePoint, myTermStruct[timepoints.size()+i,0]);
            //termstructureStore<T>::addForward(0.,timePoints[i], 0.);

        }

        // Calibrate hermite interpolation
        calibrate();

        // Calculate error
        T err(0.);
        for (auto& [myBond, myPrice] : bondData)
            err += pow(myPrice - termstructureStore<T, termstructureHermite<T>>::getBondPrice(myBond, priceDate),2.);
        mat<T> addjdiff(myTermStruct.nRows(), 1);
        // Try to punish osilation harder by increasing cost for slope change
        // Use differnce between slop and total change in the inteval
        for (num i = .1; i< 30; i+=.1){
            // Follow linear approximation
            err += pow(interpolateYieldImp(T(i+numSQRTSQRTEPS))-2.* interpolateYieldImp(T(i))+ interpolateYieldImp(T(i-numSQRTSQRTEPS)),2.)/numSQRTEPS;

        }
        return sqrt(err);
    };

    // Get initial guess by fitting Nelson-Siegel-Svenson model
    rates<T> myNSS;
    vector<bond> myBonds;
    vector<T> myPrices;
    for (auto& [myBond, myPrice] : bondData){
        myBonds.push_back(myBond);
        myPrices.push_back(myPrice);
    }
    myNSS.setNSSYieldCurveFrgBonds(myBonds, myPrices, priceDate);


    mat<T> x0(3*timepoints.size(),1);
    for (int i = 0; i< timepoints.size(); i++) {
        x0[i, 0] = myNSS.getYieldCurve(timepoints[i]);
        x0[timepoints.size()+i,0] = -(x0[i,0] -myNSS.getYieldCurve(timepoints[i]+numSQRTEPS) )/numSQRTEPS*timepoints[i] + x0[i,0];
        x0[timepoints.size()*2+i,0] = timepoints[i];
    }


    bfgs<T> myBFGS;
    myBFGS.minimize(move(errFunc), x0, 1e-15);

}

template class termstructureHermite<float>;
template class termstructureHermite<double>;
template class termstructureHermite<adjointIntegral>;