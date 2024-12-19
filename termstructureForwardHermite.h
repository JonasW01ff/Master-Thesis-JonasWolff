//
// Created by Jonas Wolff on 06/11/2024.
//

#ifndef CODE_TERMSTRUCTUREFORWARDHERMITE_H
#define CODE_TERMSTRUCTUREFORWARDHERMITE_H

#include "termstructureStore.h"
#include "interpolation.h"
#include "optimizer.h"
#include "rates.h"
#include "ThreadPool.h"

template <typename T>
class termstructureForwardHermite : public termstructureStore<T, termstructureForwardHermite<T>> {
private:

    hermiteMCspline<T> mySpline;
    vector<T> timepoints = {T(1.),T(5.),T(15.),T(30.)};
    //hermitespline<T> mySpline;

    vector<num> bondAdjErr;
    vector<array<num, 4>> paramAdjErr;

    T parallelshift = T(0);

public:
    termstructureForwardHermite(const string &name) : termstructureStore<T, termstructureForwardHermite<T>>(name) {};

    hermiteMCspline<T> getSpline() const{
        return mySpline;
    }

    const hermiteMCspline<T>& getSplineRef() const {
        return mySpline;
    }

    vector<T> getTimePoints() const{
        return timepoints;
    }

    bool bumpParam(T bump, size_t knotidx, size_t knotsubidx){
        return mySpline.bumpParam(bump, knotidx, knotsubidx);
    }

    template<typename = enable_if<is_same_v<T, adjointIntegral>>>
    explicit termstructureForwardHermite(const termstructureForwardHermite<num> &rhsTS) : mySpline(rhsTS.getSpline()){
        static_assert(is_same_v<T, adjointIntegral>);

        parallelshift = rhsTS.getParallelShift();
    }

    template<typename = enable_if<is_same_v<T, adjointIntegral>>>
    explicit termstructureForwardHermite(const termstructureForwardHermite<num> &rhsTS, vector<pair<bond, num>> bondData, date priceDate) : mySpline(rhsTS.getSpline()){
        static_assert(is_same_v<T, adjointIntegral>);

        parallelshift = rhsTS.getParallelShift();

        adjointIntegral checkpoint(0.);
        checkpoint.setPause();

        // Make bond price adjoint nodes
        vector<pair<bond, T>> bondDataAdj;
        bondDataAdj.reserve(bondData.size());
        for (int i=0; i<bondData.size();i++){
            pair<bond, T> myPair(bondData[i].first,adjointIntegral(bondData[i].second));
            myPair.second.setAdjoint(0.);
            bondDataAdj.push_back(myPair);
        }

        auto & temp = mySpline.getParamsRef();
        for (int i=0; i < temp.size(); i++){
            temp[i][0].setAdjoint(0.);
            temp[i][1].setAdjoint(0.);
            temp[i][2].setAdjoint(0.);
            temp[i][3].setAdjoint(0.);
            temp[i][4].setAdjoint(0.);
        }

        // Calculate error
        T err(0.);
        T myErr;
        bondAdjErr.resize(bondData.size());
        int idx=0;
        for (auto& [myBond, myPrice] : bondDataAdj){
            myErr = myPrice - this->getBondPrice(myBond, priceDate);
            bondAdjErr[idx] = 2.*(to_num(myErr));
            idx++;
            err += myErr*myErr;
        }
        err = err;

        err.setAdjoint(1.);
        err.propagateToPause();

        auto& myParams = mySpline.getParamsRef();
        paramAdjErr.resize(myParams.size());
        for (int i=0; i<paramAdjErr.size(); i++){
            // Get Adjoint
            paramAdjErr[i][0] = get<1>(myParams[i]).getAdjoint();
            paramAdjErr[i][1] = get<2>(myParams[i]).getAdjoint();
            paramAdjErr[i][2] = get<3>(myParams[i]).getAdjoint();
            paramAdjErr[i][3] = get<4>(myParams[i]).getAdjoint();

            // Recreate adjointNode
            myParams[i][0] = myParams[i][0].getVal();
            myParams[i][1] = myParams[i][1].getVal();
            myParams[i][2] = myParams[i][2].getVal();
            myParams[i][3] = myParams[i][3].getVal();
            myParams[i][4] = myParams[i][4].getVal();
        }
    }

    termstructureForwardHermite(const string &name, vector<pair<bond, T>> &bondData, date priceDate): termstructureStore<T, termstructureForwardHermite<T>>(name) {

        /*for (int i = 1; i < 10; i+=1)
            timePoints.push_back(i);
        for (int i = 15; i<31; i+=1)
            timePoints.push_back(i);*/

        termstructureForwardHermite<adjointIntegral> myAdjHermite(*this);

        function<void(mat<num>&,const num&, mat<num>&)> errFuncGrad = [&myAdjHermite,&bondData, &priceDate, this]
                (mat<num> &myTermStruct, const num& val, mat<num>& grad){
            // Check dimensions
            if (timepoints.size() != myTermStruct.size()/3.)
                throw runtime_error("The provided term structure does not follow the provided time line");

            // Clear last parameters
            myAdjHermite.clearAllParameters();
            // Set yields and forwards
            for (int i = 0; i<timepoints.size() ; i++){

                T myTimePoint = myTermStruct[timepoints.size()*2+i,0];
                //myTimePoint = (myTimePoint < 1.) ?  1. : myTimePoint;
                myAdjHermite.addForward(adjointIntegral(0.),adjointIntegral(myTimePoint), adjointIntegral(myTermStruct[i,0]));
                myAdjHermite.addForwardDerivative(adjointIntegral(0.),adjointIntegral(myTimePoint), adjointIntegral(myTermStruct[timepoints.size()+i,0]));
                //termstructureStore<T>::addForward(0.,timePoints[i], 0.);

            }

            // Calibrate hermite interpolation
            myAdjHermite.calibrate();

            // Calculate error
            adjointIntegral err(0.);
            adjointIntegral myErr;
            for (auto [myBond, myPrice] : bondData){
                myErr = myPrice - this->getBondPrice(myBond, priceDate);
                err += myErr*myErr;
            }

            adjointIntegral res = sqrt(err);
            res.setAdjoint(1.);
            res.propagateAll();
            auto& myParams = myAdjHermite.getAllForwardsRef();
            for (int i=0; i<myParams.size();i++){
                grad[i,0] = get<2>(myParams[i]).getAdjoint();
            }
            auto& myParams2 = myAdjHermite.getAllForwardsDerivativeRef();
            for (int i=0; i<myParams.size();i++){
                grad[myParams.size()+i,0] = get<2>(myParams[i]).getAdjoint();
            }
            myParams = myAdjHermite.getAllForwards();
            for (int i=0; i<myParams.size();i++){
                grad[myParams.size()*2+i,0] = get<1>(myParams[i]).getAdjoint()+get<1>(myParams2[i]).getAdjoint();
            }
            grad.print();
            return grad;
        };

        function<T(const mat<T>&)> errFunc = [&bondData, &priceDate, this]
                (const mat<T> &myTermStruct){

            // Check dimensions
            if (timepoints.size() != myTermStruct.size()/3.)
                throw runtime_error("The provided term structure does not follow the provided time line");

            // Clear last parameters
            termstructureStore<T, termstructureForwardHermite<T>>::clearAllParameters();
            // Set yields and forwards
            for (int i = 0; i<timepoints.size() ; i++){

                T myTimePoint = myTermStruct[timepoints.size()*2+i,0];
                //myTimePoint = (myTimePoint < 1.) ?  1. : myTimePoint;
                termstructureStore<T, termstructureForwardHermite<T>>::addForward(0.,myTimePoint, myTermStruct[i,0]);
                termstructureStore<T, termstructureForwardHermite<T>>::addForwardDerivative(0.,myTimePoint, myTermStruct[timepoints.size()+i,0]);
                //termstructureStore<T>::addForward(0.,timePoints[i], 0.);

            }

            // Calibrate hermite interpolation
            calibrate();

            // Calculate error
            T err(0.);
            T myErr;
            for (auto [myBond, myPrice] : bondData){
                myErr = myPrice - this->getBondPrice(myBond, priceDate);
                err += myErr*myErr;
            }

            return sqrt(err);
        };

        // Get initial guess by fitting Nelson-Siegel-Svenson model
        rates<T> myNSS;
        vector<bond> myBonds;
        vector<T> myPrices;
        for (auto [myBond, myPrice] : bondData){
            myBonds.push_back(myBond);
            myPrices.push_back(myPrice);
        }

        myNSS.setNSSYieldCurveFrgBonds(myBonds, myPrices, priceDate);

        mat<T> x0(3*timepoints.size(),1);
        T yield_p8, yield_m8,yield_p4, yield_m4, yield, yield_T, yield_TT;
        dual<num> yield_dual, yield_m8_dual, yield_p8_dual;
        constexpr T TEPS = (is_same_v<T, adjointIntegral>)? numEPS : numeric_limits<T>::epsilon();
        T TSQRTEPS = sqrt(TEPS);
        T TSQRTSQRTEPS = sqrt(sqrt(TEPS));
        for (int i = 0; i< timepoints.size(); i++) {

            yield = myNSS.getYieldCurve(timepoints[i]);
            dual<num> mytimepoint(to_num(timepoints[i]), 1.);
            yield_dual = myNSS.getNSSYieldCurveDual(mytimepoint);
            yield_p8_dual = myNSS.getNSSYieldCurveDual(mytimepoint+num(TSQRTEPS));
            yield_m8_dual = myNSS.getNSSYieldCurveDual(mytimepoint-num(TSQRTEPS));
            yield_p8 = myNSS.getYieldCurve(timepoints[i]+num(TSQRTEPS));
            yield_m8 = myNSS.getYieldCurve(timepoints[i]-num(TSQRTEPS));
            yield_p4 = myNSS.getYieldCurve(timepoints[i]+TSQRTSQRTEPS);
            yield_m4 = myNSS.getYieldCurve(timepoints[i]-TSQRTSQRTEPS);
            yield_T = yield_dual.getEpsilon();
            //yield_T = (yield_p8- yield_m8)/TSQRTEPS;
            yield_TT = (yield_p8_dual.getEpsilon() - yield_m8_dual.getEpsilon())/TSQRTEPS;
            //yield_TT = (yield_p4- 2.*yield + yield_m4)/TSQRTEPS;
            x0[i, 0] = yield_T*timepoints[i] + yield;
            x0[timepoints.size()+i,0] = timepoints[i]*yield_TT + 2.*yield_T;
            x0[timepoints.size()*2+i,0] = timepoints[i];
        }

        bfgs<T> myBFGS;
        //x0.print();
        /*int its;
        T myErr = errFunc(x0);

        dfpmin(x0, TEPS*2*2,its,myErr,move(temperrFunc),move(errFuncGrad));
        myErr = errFunc(x0);*/
        auto temperrFunc = errFunc;
        myBFGS.minimize(move(errFunc), x0, TEPS*2*2);
        cout << "Ending Hermite with ERROR=" << pow(temperrFunc(x0),2.)/bondData.size() << endl;

    }

    const decltype(paramAdjErr)& getParamsAdjErr() const{
        return paramAdjErr;
    }

    const decltype(bondAdjErr)& getBondAdjErr() const{
        return bondAdjErr;
    }

    void calibrate(){

        // container to keep data in
        vector<tuple<T, T, T>> myData;

        // Find yields and forwards where forwards exists with start = 0
        for (auto [myStart, myEnd, myForwardRate] : termstructureStore<T, termstructureForwardHermite<T>>::getAllForwards()){

            // Check if the forward is from today
            if (myStart != 0.)
                continue;

            // Get yield, this will throw error if doesnt exists
            T myForward = termstructureStore<T, termstructureForwardHermite<T>>::getForward(T(0.),myEnd);

            // Calculate slope
            T mySlope = termstructureStore<T, termstructureForwardHermite<T>>::getForwardDerivative(T(0.),myEnd);

            // Make tuple
            tuple<T, T, T> myPoint(myEnd, myForward, mySlope);

            // Add point
            myData.push_back(myPoint);

        };

        // Check if any data points were added
        if (myData.size() <= 0)
            throw runtime_error("No data was provided for hermite interpolation");

        // Sort my data
        sort(myData.begin(), myData.end(), []
                (const tuple<T, T, T> &first, const tuple<T, T, T> &second){
            return get<0>(first) < get<0>(second);
        });

        // Calibrate spline
        mySpline.calibrate(myData);

    }

    void setParallelShift(const T &shiftsize){
        parallelshift = shiftsize;
    }

    T getParallelShift() const{
        return parallelshift;
    }

    // Override interpolateYield with Hermite spline interpolation logic
    T interpolateYieldImp(const T &time) const {

        return mySpline.getIntegral(T(0.),time)/time + parallelshift;

    };

    T interpolateForwardImp(const T &time) const {
        return mySpline.getNode(time);
    }

    T interpolateForwardDerivativeImp(const T &time) const {
        return mySpline.getSlope(time);
    }



};

#endif //CODE_TERMSTRUCTUREFORWARDHERMITE_H
