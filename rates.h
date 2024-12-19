//
// Created by Jonas Wolff on 01/09/2024.
//

#ifndef CODE_RATES_H
#define CODE_RATES_H

#include "bond.h"
#include "string.h"
#include "array"
#include "matrix.h"
#include "optimizer.h"
#include "sobol.h"

using namespace std;
using namespace chrono;

template<class Tfp>
class rates{

    string yieldCurveModelType;
    // "NelsonModel", "Lagrange"

    mat<Tfp> paramsNelsonSiegel; // [beta, tau]

    // Discount factors: Row 1 are dates, rows 2 are discount factors
    mat<Tfp> paramsLagrange;

    sobol1D mySobol;

    bool paramsIsSet = false;
    bool modelIsCalibrated = false;

public:

    rates();

    // Calibrate yield curve with german government using lagrangian solvel Cd=p with least variance in d
    void setLinearSquareErrorYieldCurveFrgBonds(vector<bond> &bondsFRG, const vector<Tfp> &prices, const date &priceDate);

    // Calibrate yield curve with german government using lagrangian solvel Cd=p with least variance in d
    // only set discount factors at whole years.
    void setSparseLinearSquareErrorYieldCurveFrgBonds(vector<bond> &bondsFRG, const vector<Tfp> &prices, const date &priceDate);

    // Calibrate yield curve with german government bonds usign Nelson-SiegelSvensson model
    void setNSSYieldCurveFrgBonds(vector<bond> &bondsFRG, const vector<Tfp> &prices, const date &priceDate);

    mat<Tfp> getBondPriceYieldCurveGradient(bond &myBond, const date &priceDate, const Tfp &price) const;

    Tfp getYieldCurve(const Tfp &maturity) const;

    dual<num> getNSSYieldCurveDual(const dual<num> & maturity) const;

    Tfp getBondPrice(bond &myBond, const date priceDate) const;

    mat<Tfp> getYieldCurveGradient(const Tfp &maturity) const;

    const mat<Tfp>getNelsonSiegelParams() const {return paramsNelsonSiegel;};

    void setNelsonSiegelParams(const mat<Tfp> &myParams) {paramsNelsonSiegel = myParams;};

};

#endif //CODE_RATES_H
