//
// Created by Jonas Wolff on 23/09/2024.
//

#ifndef CODE_INTERPOLATION_H
#define CODE_INTERPOLATION_H
#include "matrix.h"

template<class T1>
class interpolation {

    // Calibration status
    bool status = false;

public:

    inline virtual T1 getNode(const T1 &location) const = 0;

    inline virtual T1 getNode(const T1 &&location) const = 0;

    inline virtual T1 getSlope(const T1 &location) const = 0;

    inline virtual T1 getSlope(const T1 &&location) const = 0;

    inline virtual T1 getConvexity(const T1 &location) const = 0;

    inline virtual T1 getConvexity(const T1 &&location) const = 0;

    inline virtual T1 getIntegral(const T1 &from, const T1 &location) const = 0;

    inline virtual T1 getIntegral(const T1 &&from, const T1 &&location) const = 0;

    void setStatus(bool myStatus){ status = myStatus; };

    bool getStatus() {return status;};

};


template<class T1>
class hermitespline: public interpolation<T1>{
protected:

    // Rows are interval dividers, a, b,c, d
    //      for p(x) = a + b*t + c*t^2 + d*t^3, 0<=t<=1
    vector<array<T1,5>> params;
    T1 intervalsEnd;

public:

    hermitespline(){};

    ~hermitespline(){};

    template<class fpT,typename = enable_if<is_same_v<T1, adjointIntegral> && is_floating_point_v<fpT>>>
    hermitespline(const hermitespline<fpT> &rhsSpline){

        // Allocate parameters
        params.resize(rhsSpline.params.size());

        for (int i=0; i<params.size(); i++){
            get<0>(params[i]).setVal(get<0>(rhsSpline.params[i]));
            get<1>(params[i]).setVal(get<1>(rhsSpline.params[i]));
            get<2>(params[i]).setVal(get<2>(rhsSpline.params[i]));
            get<3>(params[i]).setVal(get<3>(rhsSpline.params[i]));
            get<4>(params[i]).setVal(get<4>(rhsSpline.params[i]));
        }

        intervalsEnd = rhsSpline.intervalsEnd;
    }

    decltype(params) getParams() const{
        return params;
    }

    decltype(params)& getParamsRef(){
        return params;
    }

    T1 getInervalsEnd() const{
        return intervalsEnd;
    }

    template<typename = enable_if<is_floating_point_v<T1>>>
    bool bumpParam(const T1& bump, size_t knotidx, size_t knotsubidx){
        if (params.size() <= knotidx){
            return false;
        } else {
            if (knotsubidx > 4)
                throw runtime_error("param does not exists");
            params[knotidx][knotsubidx] += bump;
            return true;
        }
    }

    void bumpIntervalsEnd(const T1& bump){
        intervalsEnd += bump;
    }

    inline T1 getNode(const T1 &location) const override;

    inline T1 getNode(const T1 &&location) const override{
        return getNode(location);
    };

    inline T1 getSlope(const T1 &location) const override;

    inline T1 getSlope(const T1 &&location) const override{
        return getSlope(location);
    };

    inline T1 getConvexity(const T1 &location) const override;

    inline T1 getConvexity(const T1 &&location) const override{
        return getConvexity(location);
    };

    inline T1 getIntegral(const T1 &from, const T1 &location) const override;

    inline T1 getIntegral(const T1 &&from,const T1 &&location) const override{
        return getIntegral(from, location);
    };

    void calibrate(const vector<tuple<T1, T1, T1>> &data);

};

// Monotonicity Constrained Hermite Spline
template<class T1>
class hermiteMCspline: public hermitespline<T1>{

public:

    hermiteMCspline()= default;

    ~hermiteMCspline()= default;

    template<class fpT,typename = enable_if<is_same_v<T1, adjointIntegral> && is_floating_point_v<fpT>>>
    hermiteMCspline(const hermiteMCspline<fpT> &rhsSpline){
        static_assert(is_same_v<T1, adjointIntegral>);
        static_assert(is_floating_point_v<fpT>);

        // Get params
        auto rhsParams = rhsSpline.getParams();

        // Allocate parameters
        hermitespline<T1>::params.resize(rhsParams.size());
        for (int i=0; i<hermitespline<T1>::params.size(); i++){
            get<0>(hermitespline<T1>::params[i]).setVal(get<0>(rhsParams[i]));
            get<1>(hermitespline<T1>::params[i]).setVal(get<1>(rhsParams[i]));
            get<2>(hermitespline<T1>::params[i]).setVal(get<2>(rhsParams[i]));
            get<3>(hermitespline<T1>::params[i]).setVal(get<3>(rhsParams[i]));
            get<4>(hermitespline<T1>::params[i]).setVal(get<4>(rhsParams[i]));
        }

        hermitespline<T1>::intervalsEnd = rhsSpline.getInervalsEnd();
    }

    // Override calibrate to alter slopes
    void calibrate(vector<tuple<T1, T1, T1>> &data) {

        // Push slopes into Boor-Schwartz monotonicity box
        for (auto It = data.begin()+1; It<data.end()-1; It++){
            T1 sValA = get<2>(*(It-1));
            T1 sValB = get<2>(*It);
            T1 sValC = get<2>(*(It+1));

            if (sValB >= 0) {
                get<2>(*It) = min(max(T1(0.), sValB), T1(3.) * min(abs(sValB - sValA), abs(sValB - sValC)));
            }
            else
                get<2>(*It) = max(min(T1(0.), sValB),T1(-3.)*min(abs(sValB - sValA), abs(sValB - sValC)));
        }

        // Invoke regular hermite spline
        hermitespline<T1>::calibrate(data);
    }

};

template<class T>
T polynimalslope2nd(const T&x1, const T&x2, const T&x3, const T&y1, const T&y2, const T&y3, const T&loc){
    T mySlope = T(2.)*(y1*(loc-x2)*(loc-x3)/(x1-x2)/(x1-x3)+y2*(loc-x1)*(loc-x3)/(x2-x1)/(x2-x3)+y3*(loc-x1)*(loc-x2)/(x3-x1)/(x3-x2));
    return mySlope;
}

#endif //CODE_INTERPOLATION_H
