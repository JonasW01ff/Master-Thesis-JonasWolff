//
// Created by Jonas Wolff on 22/11/2024.
//

#ifndef CODE_DISTRIBUTIONFUNCTIONS_H
#define CODE_DISTRIBUTIONFUNCTIONS_H
#include <math.h>

using namespace std;

// Naive
inline float NormCDF(const float &x){
    return erfc(-x/M_SQRT2)/2.;//.5*(1+ erf(x/M_SQRT2));
}

inline double NormCDF(const double &x){
    return erfc(-x/M_SQRT2)/2.;//.5*(1+ erf(x/M_SQRT2));
}

static constexpr array<num, 4> a = {2.50662823884,-18.61500062529,41.39119773534, -25.44106049637};
static constexpr array<num, 4> b = {-8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833};
static constexpr array<num, 9> c = {0.3374754822726147, 0.976169019091718, 0.1607979714918209,
                      0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
                      0.0000321767881768, 0.0000002888167364, 0.0000003960315187};

// Glasserman 2004
inline num NormCDFI(const num &p){
    num y = p - .5;
    num r, x;
    if (abs(y)<.42){
        r = y*y;
        x = y*(((a[3]*r+a[2])*r+a[1])*r+a[0])/((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1.);
    } else {
        r = p;
        if (y>0)
            r = 1-p;
        r=log(-log(r));
        x = c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]+r*(c[5]+r*(c[6]+r*(c[7]+r*c[8])))))));
        if (y<0)
            x = -x;
    }
    return x;
}

template<class T1>
auto NormPDF(const T1 &x){
    return M_SQRT1_2*M_2_SQRTPI*.5*exp(-x*x/2.);
}



#endif //CODE_DISTRIBUTIONFUNCTIONS_H
