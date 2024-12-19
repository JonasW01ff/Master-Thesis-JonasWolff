#pragma once
#include "matrix.h"
#include "interpolation.h"
#include "date.h"

using namespace std;

// Overnight Index Swap
template<class T>
class ois
{
private:
	mat<T> df;

	hermiteMCspline<T> mySpline;

public:
	ois(const vector<pair<T,T>>& marketrates);

    ois(const date& curvedate);

    T getMarketRate(const T& tenor){
        return mySpline.getNode(tenor);
    };

    T getEffectiveRate(const date& fromDate, const date& toDate){
        // Includes fromDate, excludes toDate
        // Act/360 day counting
        T tenor = T((toDate -fromDate).count())/360.;
        T marketrate = mySpline.getNode(tenor);
        T df = 1./(1.+marketrate*tenor);

        return -log(df)/tenor;
    };

    T getDiscount(const date& fromDate, const date& toDate){
        // Includes fromDate, excludes toDate
        // Act/360 day counting
        T tenor = T((toDate -fromDate).count())/360.;
        T marketrate = mySpline.getNode(tenor);
        return 1./(1.+marketrate*tenor);

    };



};

