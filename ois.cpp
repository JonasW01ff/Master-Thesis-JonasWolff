#include "ois.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

// Data fetching
vector<double> getRatesForDate(const std::string& filePath, int year, int month, int day) {

}


template<class T>
ois<T>::ois(const vector<pair<T,T>>& marketrates) {
	// Row 1 is tenor
	// Row 2 is rates

	// Sanity check
	if (marketrates.size() <= 2)
		throw invalid_argument("At least 2 tenors need to be provided for OIS curve");

	// Allocate
	vector<tuple<T, T, T>> splineknots(marketrates.size());

	// Find first slope
	auto& [myX, myY] = marketrates[0];
	auto& [myX1, myY1] = marketrates[1];
	auto& [myX2, myY2] = marketrates[2];
    T mySlope =  (myY1-myY)/(myX1-myX); //polynimalslope2nd<T>(myX, myX1, myX2, myY, myY1, myY2, myX);
	//T mySlope = myY2 - myY + (myY1 - myY + myX * (myY2 - myY)) / (myX * myX - myX) * (T(2.) * myX - T(1.));
	splineknots[0] = make_tuple(myX, myY, mySlope);


	// Find intermidiary slopes
	for (int i = 1; i < marketrates.size()-1; i++) {
		auto& [myX, myY] = marketrates[i-1];
		auto& [myX1, myY1] = marketrates[i];
		auto& [myX2, myY2] = marketrates[i+1];
        T mySlope = (myY2-myY)/(myX2-myX); //polynimalslope2nd<T>(myX, myX1, myX2, myY, myY1, myY2, myX1);
		//mySlope = myY2 - myY + (myY1 - myY + myX1 * (myY2 - myY)) / (myX1 * myX1 - myX1) * (T(2.) * myX1 - T(1.));
		splineknots[i] = make_tuple(myX1, myY1, mySlope);
	}

	// Find Last slope
	auto& [myXend, myYend] = marketrates[marketrates.size() - 3];
	auto& [myX1end, myY1end] = marketrates[marketrates.size() - 2];
	auto& [myX2end, myY2end] = marketrates[marketrates.size() - 1];
    mySlope = (myY2end-myY1end)/(myX2end-myXend);//polynimalslope2nd<T>(myXend, myX1end, myX2end, myYend, myY1end, myY2end, myX2end);
	//mySlope = myY2end - myYend + (myY1end - myYend + myX2end * (myY2end - myYend)) / (myX2end * myX2end - myX2end) * (T(2.) * myX2end - T(1.));
	splineknots[marketrates.size() - 1] = make_tuple(myX2end, myY2end, mySlope);

	mySpline.calibrate(splineknots);


}

template<class T>
vector<pair<T,T>> ctorHelper(const date &curvedate) {
    int myYear = curvedate.year().operator int();
    int myMonth = curvedate.month().operator unsigned int();
    int myDay = curvedate.day().operator unsigned int();
    static const string estrlocation = "/at/your/location.data";
    vector<T> fileTenors = {0.,1./12., 3./12., 6./12., 9./12., 1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.};
    vector<double> fileRates = getRatesForDate(estrlocation, myYear, myMonth, myDay);
    vector<pair<T,T>> marketrates(fileTenors.size());
    for (int i=0; i<fileTenors.size();i++)
        marketrates[i] = make_pair(fileTenors[i], fileRates[i]);
    return marketrates;
}

template<class T>
ois<T>::ois(const date &curvedate): ois(ctorHelper<T>(curvedate)) {};

template class ois<float>;
template class ois<double>;