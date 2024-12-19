//
// Created by Jonas Wolff on 31/10/2024.
//

#ifndef CODE_TERMSTRUCTUREHERMITE_H
#define CODE_TERMSTRUCTUREHERMITE_H


#include "termstructureStore.h"
#include "interpolation.h"
#include "optimizer.h"
#include "rates.h"

template <typename T>
class termstructureHermite : public termstructureStore<T, termstructureHermite<T>> {
private:

    hermiteMCspline<T> mySpline;


    vector<T> timepoints = {T(1.),T(5.),T(15.),T(30.)};
    //hermitespline<T> mySpline;

public:
    termstructureHermite(const string& name) : termstructureStore<T, termstructureHermite<T>>(name){};
    termstructureHermite(const string& name, vector<pair<bond, T>> &bondData, date priceDate);

    void calibrate(){

        // container to keep data in
        vector<tuple<T, T, T>> myData;

        // Find yields and forwards where forwards exists with start = 0
        for (auto [myStart, myEnd, myForwardRate] : termstructureStore<T, termstructureHermite<T>>::getAllForwards()){

            // Check if the forward is from today
            if (myStart != 0.)
                continue;

            // Get yield, this will throw error if doesnt exists
            T myYield = termstructureStore<T, termstructureHermite<T>>::getYield(myEnd);

            // Calculate slope
            T mySlope = (myForwardRate - myYield)/(myEnd);

            // Make tuple
            tuple<T, T, T> myPoint(myEnd, myYield, mySlope);

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

    // Override interpolateYield with Hermite spline interpolation logic
    T interpolateYieldImp(const T &time) const {

        return mySpline.getNode(time);

    };

    T interpolateForwardImp(const T &time) const {
        return time*mySpline.getSlope(time) +  mySpline.getNode(time);
    }

    T interpolateForwardDerivativeImp(const T &time) const {
        return 2.*mySpline.getSlope(time) + time*mySpline.getConvexity(time);
    }

};


#endif //CODE_TERMSTRUCTUREHERMITE_H
