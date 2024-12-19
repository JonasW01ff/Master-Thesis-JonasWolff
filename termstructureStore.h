//
// Created by Jonas Wolff on 31/10/2024.
//

#ifndef CODE_TERMSTRUCTURESTORE_H
#define CODE_TERMSTRUCTURESTORE_H

#include <vector>
#include <string>
#include <stdexcept>
#include "numbers.h"
#include "matrix.h"
#include "bond.h"

template <class T, class myTermStruct>
class termstructureStore {
    static_assert(is_floating_point_v<T> or std::is_same_v<T, dual<num>> or is_same_v<T, adjointIntegral>, "T must be either double or dual");

private:
    // Helper methods
    bool yieldExists(const T time) const{

        for (pair myYield : yields)
            if (myYield.first == time) return true;

        return false;

    };
    bool forwardExists(const T start,const T end) const{

        for (tuple myForward : forwards)
            if (get<1>(myForward) == start and get<2>(myForward) == end) return true;

        return false;

    };

    // Data Members
    std::string storeName;
    vector<pair<T, T>> yields; // Using Matrix to store yields
    vector<tuple<T, T, T>> forwards; // start, end, forwardRate
    vector<tuple<T, T, T>> forwardsDerivative; // start, end, forwardRate

public:
    // Constructors
    termstructureStore() = default;
    termstructureStore(const std::string& name){

        // Set model name
        storeName = name;

    };

    // Methods to add and retrieve yields
    void addYield(const T time, const T yield){

        // Creat time and yield pair
        pair<T, T> myYield(time, yield);

        // add myYield
        yields.push_back(myYield);

    };
    T getYield(const T time) const{

        for (pair myYield : yields)
            if (myYield.first == time) return myYield.second;
        throw runtime_error("time point is not found on the stored yields");

    };

    // Methods to add and retrieve forwards
    void addForward(T start, T end, T forwardRate){

        // Make tuple
        tuple<T, T, T> myForward(start, end, forwardRate);

        // Insert forward rate
        forwards.push_back(myForward);

    };
    T getForward(T start, T end) const{

        for (auto [myStart, myEnd, myForward] : forwards)
            if (myStart == start and myEnd == end) return myForward;

        throw runtime_error("time point is not found on the stored forwards");

    };

    // Methods to add and retrieve forwards
    void addForwardDerivative(T start, T end, T forwardRateD){

        // Make tuple
        tuple<T, T, T> myForwardD(start, end, forwardRateD);

        // Insert forward rate
        forwardsDerivative.push_back(myForwardD);

    };
    T getForwardDerivative(T start, T end) const{

        for (auto [myStart, myEnd, myForwardD] : forwardsDerivative)
            if (myStart == start and myEnd == end) return myForwardD;

        throw runtime_error("time point is not found on the stored forwards derivative");

    };

    // Methods to get all yields or forwards
    vector<pair<T, T>> getAllYields() const {
        return yields;
    };
    vector<tuple<T, T, T>> getAllForwards() const {
        return forwards;
    };

    vector<tuple<T, T, T>>& getAllForwardsRef() {
        return forwards;
    };

    vector<tuple<T, T, T>> getAllForwardsDerivative() const {
        return forwardsDerivative;
    };

    vector<tuple<T, T, T>>& getAllForwardsDerivativeRef() {
        return forwardsDerivative;
    };

    // Methods that clear all parameters
    void clearAllParameters() {

        // Reset parameters
        yields.clear();
        forwards.clear();
        forwardsDerivative.clear();
    };

    // Setters and getters for the name
    void setName(const string& name){
        storeName = name;
    };
    string getName() const{
        return storeName;
    };

    // Method for interpolation
    T interpolateYield(const T &time) const {

        return static_cast<const myTermStruct*>(this)->interpolateYieldImp(time);

    };

    // Method for interpolation
    T interpolateForward(const T &time) const {
        return static_cast<const myTermStruct*>(this)->interpolateForwardImp(time);
    };

    // Method for interpolation
    T interpolateForwardDerivative(const T &time) const {
        return static_cast<const myTermStruct*>(this)->interpolateForwardDerivativeImp(time);
    };

    // Method for interpolation
    T interpolateDiscount(const T &time) const{
        return exp(-interpolateYield(time)*time);
    }

    T interpolateDiscountDerivative(const T &time) const{
        return -interpolateDiscount(time)*interpolateForward(time);
    }

    T interpolateDiscountDerivative(const T &df, const T &forward) const{
        return -df*forward;
    }

    T interpolateDiscountSecondDerivative(const T &time) const{
        T df = interpolateDiscount(time);
        T forward = interpolateForward(time);
        T forwardT = interpolateForwardDerivative(time);
        return -df*forwardT +df*forward*forward;
    }

    T interpolateDiscountSecondDerivative(const T&df,const T &forward, const T &forwardT) const{
        return -df*forwardT +df*forward*forward;
    }

    // method for pricing a bond
    T getBondPrice(const bond myBond, const date priceDate) const {

        // Get specs
        date firstcouponDate = myBond.getFirstcouponDate();
        const date& intrestStartDate = myBond.getIntrestStartDate();
        const date& maturityDate = myBond.getMaturityDate();
        const months& couponFrequency = myBond.getCouponFrequency();
        const T& couponSize = T(myBond.getCouponSize());

        // Check first coupon date
        firstcouponDate = date(firstcouponDate.year(), maturityDate.month(), maturityDate.day());


        // Get accrued interest
        T ai;
        date lastcouponDate, nextcouponDate;
        if (priceDate >= firstcouponDate){

            // Get time of last coupon
            lastcouponDate = firstcouponDate;
            while (lastcouponDate <= priceDate)
                lastcouponDate += couponFrequency;
            lastcouponDate -= couponFrequency;


            // Get accrued interest (priceDate in also counted inclusive)
            nextcouponDate = date(lastcouponDate + couponFrequency);
            if (couponFrequency.count() == 0) throw runtime_error("couponFrequency cant be zero.");
            ai = couponSize*(priceDate- lastcouponDate).count()/(nextcouponDate - lastcouponDate).count();

        } else if (priceDate >= intrestStartDate) {

            // Get accrued interest
            date firstcouponDateStart(firstcouponDate - couponFrequency);
            date extracouponDateStart(firstcouponDate - 2*couponFrequency);
            if (couponFrequency.count() == 0) throw runtime_error("couponFrequency cant be zero.");
            if (priceDate>=firstcouponDateStart){
                ai = couponSize*(firstcouponDateStart-intrestStartDate).count()/(firstcouponDateStart-extracouponDateStart).count();
                ai += couponSize*(priceDate- firstcouponDateStart).count()/(firstcouponDate - firstcouponDateStart).count();}
            else
                ai = couponSize*(priceDate- firstcouponDateStart).count()/(firstcouponDate - firstcouponDateStart).count();

        } else {

            // No accrued interest
            ai = 0;
        }

        // Discount future payments
        T cleanPrice(0.);
        nextcouponDate = firstcouponDate;
        while (nextcouponDate <= priceDate)
            nextcouponDate += couponFrequency;

        T yearsstoCoupon;
        while (nextcouponDate <= maturityDate){
            yearsstoCoupon = priceDate.yearsuntil(nextcouponDate);
            cleanPrice += couponSize*exp(-interpolateYield(yearsstoCoupon)*yearsstoCoupon);
            nextcouponDate += couponFrequency;
        }

        // Check that coupon on maturity date has been given
        if (nextcouponDate - couponFrequency != maturityDate)
            throw runtime_error("Coupon schedule does not match bond maturity");

        yearsstoCoupon = priceDate.yearsuntil(maturityDate);
        cleanPrice += 100.0*exp(-interpolateYield(yearsstoCoupon)*yearsstoCoupon);

        return cleanPrice;// - ai;

    }


    T getYTM(const mat<T>& CF, const date& priceDate, const T &myPrice) {

        if (CF.nRows() != 2)
            throw runtime_error("This processes one bond at a time");

        // Set variables
        T myRateLow = 0.0001;
        T myRateHigh = 2.;

        // Dynamically adjust accuarcy
        T TSQRTEPS = sqrt(numeric_limits<T>::epsilon());
        auto calcPrice = [&CF](const T& myRate)->T {
                T tempPrice = 0;
                T* row0 = CF[0];
                T* row1 = CF[1];
                for (int i = 0; i < CF.nCols(); i++) {
                    tempPrice += row1[i] * exp(-myRate * row0[i]);
                }
                return tempPrice;
        };

        // Commence bisection
        T ytm = bisection<T>(calcPrice,myRateLow, myRateHigh);

        return ytm;
    }

    // Macaulay duration
    T getDuration(bond& myBond, const date& priceDate) {
        mat<T> CF = myBond.getCashFlow(priceDate);
        T myPrice = getBondPrice(myBond);
        T dur = T(0.);
        T* row0 = CF[0];
        T* row1 = CF[1];
        T ytm = getYTM(CF, priceDate, myPrice);
        for (int i = 0; i < CF.nCols(); i++) {
            dur += row0[i] * row1[i] * exp(-row0[i]*ytm);
        }
        dur /= myPrice;
        return dur;
    }

    // Modified Duration
    T getMD(bond& myBond, const date& priceDate){
        mat<T> CF = myBond.getCashFlow(priceDate);
        T myPrice = getBondPrice(myBond);
        T dur = T(0.);
        T* row0 = CF[0];
        T* row1 = CF[1];
        T ytm = getYTM(CF, priceDate, myPrice);
        for (int i = 0; i < CF.nCols(); i++) {
            dur += row0[i] * row1[i] * exp(-row0[i]*ytm);
        }
        dur /= myPrice;
        return dur/(1.+ytm/T(myBond.getCouponFrequency()*12));
    }



};

#endif // CODE_TERMSTRUCTURESTORE_H

