//
// Created by Jonas Wolff on 23/09/2024.
//

#include "interpolation.h"
#include "numbers.h"

template<class T1>
void hermitespline<T1>::calibrate(const vector<tuple<T1, T1, T1>> &data) {

    // Check that at least 3 knots are given
    if (data.size() < 3)
        throw runtime_error("Hermite spline needs at least 3 knots to calibrate");

    // Allocate memory for parameters
    params.resize(data.size()-1);

    // Calculate coefficients
    for (auto It = data.begin(); It<data.end()-1; It++){

        //tuple<T1, T1, T1> pointA;
        tuple<T1, T1, T1> pointB = *(It);
        tuple<T1, T1, T1> pointC = *(It+1);

        int i = distance(data.begin(), It);

        T1 h = get<0>(pointC)- get<0>(pointB);

        // current interval startpoint
        get<0>(params[i]) = get<0>(pointB); //*xVal;

        // a
        get<1>(params[i]) = get<1>(pointB); //*yVal;

        // b
        get<2>(params[i]) = get<2>(pointB); //*sVal;

        // c
        get<3>(params[i]) = (3 * (get<1>(pointC) - get<1>(pointB)) - h * (2 * get<2>(pointB) + get<2>(pointC))) / (h * h);

        // d
        get<4>(params[i]) = (2 * (get<1>(pointB) - get<1>(pointC)) + h * (get<2>(pointB) + get<2>(pointC))) / (h * h * h);

        //get<4>(params[i]) = get<2>(pointC)*h -2*(get<1>(pointC) - get<1>(pointB)) + get<2>(pointB)*h;
        // *sValNext - 2*(*yValNext - *yVal) + *sVal;

    }

    // Set intervals endpoint
    intervalsEnd = get<0>(*(data.end()-1));

    // Mark spline as calibrated
    this->setStatus(true);

}


template<class T1>
inline T1 hermitespline<T1>::getNode(const T1 &location) const {


    T1 t, t2, t3; // Unit distance inside interpolation interval
    const T1 &firstX = get<0>(params[0]);
    if (location < firstX){
        const T1 &firstA = get<1>(params[0]);
        const T1 &firstB = get<2>(params[0]);
        const T1 &firstC = get<3>(params[0]);
        const T1 &firstD = get<4>(params[0]);
        t = location - firstX;
        t2 = t*t;
        t3 = t2*t;
        return firstA + firstB*t + firstC*t2 + firstD*t3;
    }

    // Get correct interval
    int i = 0;
    for (auto& point : params) {
        if (get<0>(point) < location) i++;
        else break;
    }


    const T1 &priorX = get<0>(params[i-1]);
    const T1 &priorA = get<1>(params[i-1]);
    const T1 &priorB = get<2>(params[i-1]);
    const T1 &priorC = get<3>(params[i-1]);
    const T1 &priorD = get<4>(params[i-1]);
    const T1 &lastX = get<0>(*(params.end()-1));

    // Check if we are outside interval
    if (location >= intervalsEnd)
        return priorA + priorB + priorC + priorD; //params[1,i-1] + params[2,i-1] + params[3,i-1] + params[4,i-1];


    t = location- priorX;

    // return interpolated node
    t2 = t*t;
    t3 = t2*t;
    return priorA + priorB*t + priorC*t2 + priorD*t3;
    //return params[1,i-1] + params[2,i-1]*t + params[3,i-1]*pow(t,2) + params[4,i-1]*pow(t,3);
}

template<class T1>
inline T1 hermitespline<T1>::getSlope(const T1 &location) const {

    // Get correct interval
    int i = 0;
    for (auto& point : params) {
        if (get<0>(point) < location) i++;
        else break;
    }

    T1 t; // Unit distance inside interpolation interval
    T1 firstX = get<0>(params[0]);
    T1 firstA = get<1>(params[0]);
    T1 firstB = get<2>(params[0]);
    T1 firstC = get<3>(params[0]);
    T1 firstD = get<4>(params[0]);
    if (location <= firstX){
        t = location - firstX;
        return firstB + 2.*firstC*t + 3.*firstD*t*t;
    }


    T1 priorX = get<0>(params[i-1]);
    T1 priorA = get<1>(params[i-1]);
    T1 priorB = get<2>(params[i-1]);
    T1 priorC = get<3>(params[i-1]);
    T1 priorD = get<4>(params[i-1]);
    T1 lastX = get<0>(*(params.end()-1));

    // Check if we are outside interval
    if (location >= intervalsEnd)
        return priorB + 2.*priorC + 3.*priorD;//params[1,i-1] + params[2,i-1] + params[3,i-1] + params[4,i-1];
    else if (location >= lastX ) // params[0, params.nCols() - 1])
        t = location - priorX; //(location- params[0, i-1])/(intervalsEnd - params[0, i-1]);
    else
        t = location - priorX;

    // return interpolated node
    return priorB + 2.*priorC*t + 3.*priorD*t*t;
    //return params[1,i-1] + params[2,i-1]*t + params[3,i-1]*pow(t,2) + params[4,i-1]*pow(t,3);
}


template<class T1>
inline T1 hermitespline<T1>::getConvexity(const T1 &location) const {

    // Get correct interval
    int i = 0;
    for (auto point : params) {
        if (get<0>(point) < location) i++;
        else break;
    }

    T1 t; // Unit distance inside interpolation interval
    T1 firstX = get<0>(params[0]);
    T1 firstA = get<1>(params[0]);
    T1 firstB = get<2>(params[0]);
    T1 firstC = get<3>(params[0]);
    T1 firstD = get<4>(params[0]);
    if (location <= firstX){
        t = location - firstX;
        return 2.*firstC + 6.*firstD*t;
    }


    T1 priorX = get<0>(params[i-1]);
    T1 priorA = get<1>(params[i-1]);
    T1 priorB = get<2>(params[i-1]);
    T1 priorC = get<3>(params[i-1]);
    T1 priorD = get<4>(params[i-1]);
    T1 lastX = get<0>(*(params.end()-1));

    // Check if we are outside interval
    if (location >= intervalsEnd)
        return 2.*priorC + 6.*priorD ;//params[1,i-1] + params[2,i-1] + params[3,i-1] + params[4,i-1];
    else if (location >= lastX )// params[0, params.nCols() - 1])
        t = location - priorX; //(location- params[0, i-1])/(intervalsEnd - params[0, i-1]);
    else
        t = location- priorX;

    // return interpolated node
    return 2.*priorC + 6.*priorD*t;
    //return params[1,i-1] + params[2,i-1]*t + params[3,i-1]*pow(t,2) + params[4,i-1]*pow(t,3);
}

template<class T1>
inline T1 hermitespline<T1>::getIntegral(const T1 &from, const T1 &location) const {

    T1 myInt(0.);
    if (from != 0.)
        throw runtime_error("This is not implemented");

    // Integrate up to start
    const auto& [firstX, firstA, firstB, firstC, firstD] = params[0];
    T1 to = min(firstX, location);
    T1 delta = from-firstX;
    T1 delta2 = delta*delta,
            delta3 = delta2*delta,
            delta4 = delta3*delta;
    myInt -= firstA*delta
            + .5*firstB*delta2
            + 1./3.*firstC*delta3
            + 1./4.*firstD*delta4;
    delta = to-firstX;
    delta2 = delta*delta,
    delta3 = delta2*delta,
    delta4 = delta3*delta;
    myInt += firstA*delta
             + .5*firstB*delta2
             + 1./3.*firstC*delta3
             + 1./4.*firstD*delta4;

    int i = 0;
    T1 priorX;
    T1 priorA, priorB, priorC, priorD;
    T1 h;
    for (const auto& [myX, myA, myB, myC, myD] : params){

        if (i==0) {
            priorX = myX, priorA = myA, priorB = myB, priorC=myC, priorD=myD;
            i++;
            continue;
        }

        if (location <= priorX) break;

        to = min(location, myX);
        delta = to-priorX;
        delta2 = delta*delta,
        delta3 = delta2*delta,
        delta4 = delta3*delta;
        myInt += priorA*delta
                 + .5*priorB*delta2
                 + 1./3.*priorC*delta3
                 + 1./4.*priorD*delta4;

        if (i == params.size()-1 && location > myX) {
            to = min(location, intervalsEnd);
            delta = to-myX;
            delta2 = delta*delta,
            delta3 = delta2*delta,
            delta4 = delta3*delta;
            myInt += myA*delta
                     + .5*myB*delta2
                     + 1./3.*myC*delta3
                     + 1./4.*myD*delta4;
        }

        if (i == params.size()-1 && location > myX && location > intervalsEnd){
            h = myX - priorX;
            delta2 = h*h;
            delta3 = delta2*h;
            myInt += (location - myX)*(myA+
                                            myB*h+
                                            myC*delta2+
                                            myD*delta3);
        }

        priorX = myX, priorA = myA, priorB = myB, priorC=myC, priorD=myD;

        i++;
    }

    return myInt;

}

template class hermitespline<double>;
template class hermitespline<float>;
template class hermitespline<adjointIntegral>;
template class hermitespline<dual<double>>;
template class hermitespline<dual<float>>;