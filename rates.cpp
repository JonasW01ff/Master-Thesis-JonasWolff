//
// Created by Jonas Wolff on 01/09/2024.
//

#include "rates.h"
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <numbers>

template<class Tfp>
rates<Tfp>::rates() {};

template<class Tfp>
Tfp rates<Tfp>::getBondPrice(bond &myBond, const date priceDate) const {

    // Get specs
    date firstcouponDate = myBond.getFirstcouponDate();
    const date& intrestStartDate = myBond.getIntrestStartDate();
    const date& maturityDate = myBond.getMaturityDate();
    const months& couponFrequency = myBond.getCouponFrequency();
    const Tfp couponSize(myBond.getCouponSize());

    // Check first coupon date
    firstcouponDate = date(firstcouponDate.year(), maturityDate.month(), maturityDate.day());


    // Get accrued interest
    Tfp ai;
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
    Tfp cleanPrice(0);
    nextcouponDate = firstcouponDate;
    while (nextcouponDate <= priceDate)
        nextcouponDate += couponFrequency;

    Tfp yearsstoCoupon;
    while (nextcouponDate <= maturityDate){
        yearsstoCoupon = priceDate.yearsuntil(nextcouponDate);
        cleanPrice += couponSize*exp(-getYieldCurve(yearsstoCoupon)*yearsstoCoupon);
        nextcouponDate += couponFrequency;
    }

    // Check that coupon on maturity date has been given
    if (nextcouponDate - couponFrequency != maturityDate)
        throw runtime_error("Coupon schedule does not match bond maturity");

    yearsstoCoupon = priceDate.yearsuntil(maturityDate);
    cleanPrice += 100.0*exp(-getYieldCurve(yearsstoCoupon)*yearsstoCoupon);

    return cleanPrice;// - ai;

}

template<class Tfp>
mat<Tfp> rates<Tfp>::getBondPriceYieldCurveGradient(bond &myBond, const date &priceDate, const Tfp &price) const {

    // Get specs
    const date& firstcouponDate = myBond.getFirstcouponDate();
    const date& intrestStartDate = myBond.getIntrestStartDate();
    const date& maturityDate = myBond.getMaturityDate();
    const months& couponFrequency = myBond.getCouponFrequency();
    const Tfp couponSize(myBond.getCouponSize());


    date lastcouponDate, nextcouponDate;

    // Discount future payments
    mat<Tfp> gradient(6,1);
    nextcouponDate = firstcouponDate;
    while (nextcouponDate <= priceDate)
        nextcouponDate += couponFrequency;

    Tfp yearsstoCoupon, yield_t;
    mat<Tfp> yieldGradient;
    while (nextcouponDate <= maturityDate){
        yearsstoCoupon = priceDate.yearsuntil(nextcouponDate);
        yieldGradient = getYieldCurveGradient(yearsstoCoupon);
        yield_t = getYieldCurve(yearsstoCoupon);

        // gradient
        gradient += yieldGradient*yearsstoCoupon*couponSize*exp(-yield_t*yearsstoCoupon);

        nextcouponDate += couponFrequency;

    }
    yearsstoCoupon = priceDate.yearsuntil(maturityDate);
    yieldGradient = getYieldCurveGradient(yearsstoCoupon);
    yield_t = getYieldCurve(yearsstoCoupon);

    // gradient at bullet
    gradient += yieldGradient*yearsstoCoupon*Tfp(100.)*exp(-yield_t*yearsstoCoupon);
    return gradient*Tfp(2.)*(price - getBondPrice(myBond, priceDate));

}

// data fetching
std::vector<double> getNSStartPramsAtDate(const date & mydate) {

}


// Set Nelson Siegel Svensson Yield Curve
template<class Tfp>
void rates<Tfp>::setNSSYieldCurveFrgBonds(vector<bond> &bondsFRG, const vector<Tfp> &prices, const date &priceDate) {

    // Check size
    if (bondsFRG.size() != prices.size())
        throw runtime_error("bonds vector did not match price vector");

    // Set model type
    this->yieldCurveModelType = "NelsonModel";

    // Set intial parameters
    /*mat<Tfp> initParams(6,1,Tfp(0.5)); //ECB
    initParams.setEntry(1,0,Tfp(3.15));
    initParams.setEntry(2,0,Tfp(-1.13));
    initParams.setEntry(3,0,Tfp(6.28));
    initParams.setEntry(4,0,Tfp(1.2));
    initParams.setEntry(5,0,Tfp(12.7));*/
    vector<double> initParamsECB;
    initParamsECB = getNSStartPramsAtDate(priceDate);
    mat<Tfp> initParams(6,1); // Best
    initParams[0,0] = initParamsECB[1];
    initParams[1,0] = initParamsECB[2];
    initParams[2,0] = initParamsECB[3];
    initParams[3,0] = initParamsECB[4];
    initParams[4,0] = initParamsECB[5];
    initParams[5,0] = initParamsECB[6];

    setNelsonSiegelParams(initParams);
    paramsIsSet = true;
    // minimize yieldcurvegradient
    auto fsqrErr = [&] (bond &myBond,const Tfp& myPrice) {
        return pow(myPrice - getBondPrice(myBond,priceDate) ,2);};
    Tfp sqrError = inner_product(bondsFRG.begin(), bondsFRG.end(),
                                    prices.begin(), Tfp(0.0),
                                    std::plus<>(), fsqrErr)/
                                            bondsFRG.size();

    // Minimize square error
    function<Tfp(const mat<Tfp>&)> minimizefunc = [&bondsFRG, &prices,&fsqrErr, this](const mat<Tfp> &myParmas)->Tfp{
        setNelsonSiegelParams(myParmas);
        return inner_product(bondsFRG.begin(), bondsFRG.end(),
                             prices.begin(), Tfp(0.0),
                             std::plus<>(), fsqrErr)/
               bondsFRG.size();
    };
    auto gradfunc = [&bondsFRG, &prices, this, &priceDate](const mat<Tfp> &params, const Tfp&,mat<Tfp> &grad){
        setNelsonSiegelParams(params);
        for (int i=0; i<grad.size(); i++)
            grad[i,0] = 0.;
        for (int i=0; i<bondsFRG.size(); i++){
            bond &myBond = bondsFRG[i];
            const Tfp &myPrice = prices[i];
            Tfp estPrice = getBondPrice(myBond, priceDate);
            grad -= getBondPriceYieldCurveGradient(myBond, priceDate, myPrice);
        }
        grad = grad / Tfp(bondsFRG.size());
    };

    bfgs<Tfp> myBFGS;
    int i= 0;
    Tfp olderr = minimizefunc(getNelsonSiegelParams());
    Tfp newerr;
    mat<Tfp> oldparams = initParams;
    while (olderr > .1) {
        if (i>=1)
            break;
        i++;
        //cout << "Using new NSS start params ERROR=" + to_string(olderr) << endl;
        //if (i>10000)
        //    return;//throw runtime_error("Nelson Siegel Svenson did not converge: ERROR=" + to_string(minimizefunc(getNelsonSiegelParams())));
        if (i>1) {
            for (int i = 0; i < initParams.size(); i++)
                initParams[i, 0] = oldparams[i, 0] += mySobol.getN();
        }
        /*initParams[0,0] = min(max(initParams[0,0], Tfp(-15.)),Tfp(15.));
        initParams[1,0] = min(max(initParams[1,0], Tfp(-15.)),Tfp(15.));
        initParams[2,0] = min(max(initParams[2,0], Tfp(-15.)),Tfp(15.));
        initParams[3,0] = min(max(initParams[3,0], Tfp(-15.)),Tfp(15.));*/
        //initParams[4,0] = min(abs(initParams[4,0]),Tfp(5.));
        //initParams[5,0] = min(abs(initParams[5,0]), Tfp(15.));
        //if (initParams[4,0]>initParams[5,0])
        //    swap(initParams[4,0], initParams[5,0]);
        int its;
        Tfp fret = minimizefunc(initParams);
        auto minimizefunccopy = minimizefunc;
        auto gradfunccopy = gradfunc;
        dfpmin<Tfp>(initParams,numEPS*2*2,its, fret,move(minimizefunccopy),move(gradfunccopy));
        //myBFGS.minimize(minimizefunc, initParams);
        newerr = fret;//minimizefunc(initParams);
        if (newerr<olderr) {
            cout << "Evolution " +to_string(i)+ ". Using new NSS start params ERROR=" + to_string(olderr) << endl;
            olderr = newerr;
            oldparams = initParams;
        }
    }
    if (olderr > 1.0)
        throw runtime_error("NSS err too great.");
    cout << "Ending with NSS with ERROR=" + to_string(olderr) << endl;
    return;

    int rounds = 1;
    int maxrounds = 1000;
    Tfp bestSqrError = sqrError;
    mat<Tfp> bestParams = initParams;
    int triesSinceBest = 0;
    while (sqrError > 0.1  and rounds<maxrounds){
        rounds++;


        // Find sum gradient for each bond
        mat<Tfp> gradient(6,1);
        for (int i=0; i<bondsFRG.size(); i++){
            bond &myBond = bondsFRG[i];
            const Tfp &myPrice = prices[i];
            gradient += getBondPriceYieldCurveGradient(myBond, priceDate, myPrice);
        }



        gradient = gradient*(Tfp(1.0)/bondsFRG.size());
        Tfp normconst = (mat<Tfp>(6,1)+Tfp(1.0)).T().matmul(gradient)[0,0]/gradient.size();
        if (abs(normconst) > 10)
            gradient.setEntries([&normconst](const Tfp &a){return a/abs(normconst);});


        // update params
        mat<Tfp> updatedParams;
        if (sqrError > 10)
            updatedParams = getNelsonSiegelParams()-gradient*Tfp(0.05);
        else if (sqrError > 3)
            updatedParams = getNelsonSiegelParams()-gradient*Tfp(0.01);
        else if (sqrError > 1)
            updatedParams = getNelsonSiegelParams()-gradient*Tfp(0.01);
        else if (sqrError > 0.6)
            updatedParams = getNelsonSiegelParams()-gradient*Tfp(0.001);
        else
            updatedParams = getNelsonSiegelParams()-gradient*Tfp(0.0001);

        if (updatedParams[4,0] < 0.3)
            updatedParams[4,0] = 0.3;
        if (updatedParams[5,0] < 0.3)
            updatedParams[5,0] = 0.3;

        setNelsonSiegelParams(updatedParams);



        sqrError = inner_product(bondsFRG.begin(), bondsFRG.end(),
                                 prices.begin(), Tfp(0.0),
                                 std::plus<>(), fsqrErr)/
                   bondsFRG.size();

        if (sqrError < bestSqrError && sqrError>0.01 && not isnan(to_num(sqrError))){
            triesSinceBest = 0;
            bestSqrError = sqrError;
            bestParams = getNelsonSiegelParams();
        } else
            triesSinceBest++;


    }
    setNelsonSiegelParams(bestParams);
    modelIsCalibrated = true;

    if (rounds >= maxrounds && sqrError > 1.0)
        modelIsCalibrated = false;
        throw runtime_error("Yield Curve could not calibrate. The average square error was:" + to_string(bestSqrError));

}

template<class Tfp>
dual<num> rates<Tfp>::getNSSYieldCurveDual(const dual<num> &maturity) const {
    // Return yield curve point
    const Tfp* beta = paramsNelsonSiegel[0];
    const Tfp* tau = paramsNelsonSiegel[4];
    dual<num> res(to_num(beta[0]),num(0));
    dual<num> exptemp1 = exp(-maturity/to_num(tau[0]));
    dual<num> exptemp2 = exp(-maturity/to_num(tau[1]));
    dual<num> temp = (1.-exptemp1)/maturity*to_num(tau[0]);
    res = res + temp*to_num(beta[1]);
    res = res + (temp- exptemp1)*to_num(beta[2]);
    res = res + ((1.0-exptemp2)/maturity*to_num(tau[1])- exptemp2 )*to_num(beta[3]);

    return res/100.;
}

template<class Tfp>
Tfp rates<Tfp>::getYieldCurve(const Tfp &maturity) const {

    // Check if model type is known
    if (yieldCurveModelType != "NelsonModel" and yieldCurveModelType != "Lagrange") {
        throw runtime_error("Yield curve model is not supported");
    }

    // Check if parameters is set
    if (paramsIsSet == false)
        throw runtime_error("No parameter is set for Yield Curve Model");

    if (yieldCurveModelType == "NelsonModel"){
        // Return yield curve point
        const Tfp* beta = paramsNelsonSiegel[0];
        const Tfp* tau = paramsNelsonSiegel[4];
        Tfp res = beta[0];
        Tfp exptemp1 = exp(-maturity/tau[0]);
        Tfp exptemp2 = exp(-maturity/tau[1]);
        Tfp temp = (1.-exptemp1)/maturity*tau[0];
        res += beta[1]*temp;
        res += beta[2]*(temp- exptemp1);
        res += beta[3]*((1.0-exptemp2)/maturity*tau[1]- exptemp2 );

        return res/100.;
    } else if (yieldCurveModelType == "Lagrange") {

        const auto maturityUpperBound = upper_bound(paramsLagrange.begin(),
                                                    paramsLagrange.begin()+ paramsLagrange.nCols(),
                                                    maturity);
        const auto maturityLowerBound = maturityUpperBound - 1;

        const auto discountUpperBound = maturityUpperBound + paramsLagrange.nCols();
        const auto discountLowerBound = discountUpperBound - 1;

        // Check that we not out of bound
        if (maturityLowerBound < paramsLagrange.begin())
            return paramsLagrange[1,0];
        else if (maturityUpperBound >= paramsLagrange.begin()+ paramsLagrange.nCols())
            return paramsLagrange[1,paramsLagrange.nCols()-1];

        return (*discountUpperBound - *discountLowerBound)*(maturity - *maturityLowerBound)/
                            (*maturityUpperBound - *maturityLowerBound) + *discountLowerBound;

    } else {
        throw runtime_error("Model is not supported");
    }

}

template<class Tfp>
mat<Tfp> rates<Tfp>::getYieldCurveGradient(const Tfp &maturity) const{

    // Check model type
    if (yieldCurveModelType != "NelsonModel")
        throw runtime_error("Model Not supported for Yield Curve Gradient.");

    Tfp beta0, beta1, beta2, beta3, tau1, tau2;

    const Tfp* beta = paramsNelsonSiegel[0];
    const Tfp* tau = paramsNelsonSiegel[4];

    Tfp frac1 = maturity/tau[0];
    Tfp frac2 = maturity/tau[1];

    Tfp exp1 = exp(-frac1);
    Tfp exp2 = exp(-frac2);

    Tfp bracket1 = (1.0 - exp1)/(frac1);
    Tfp bracket2 = (1.0 - exp1)/(frac2);

    Tfp delta1 = (1.0 - exp1*(frac1+1.0))/maturity;
    Tfp delta2 = (1.0 - exp2*(frac2+1.0))/maturity;

    beta0 = 1.;
    beta1 = bracket1;
    beta2 = bracket1 - exp1;
    beta3 = bracket2 - exp2;

    tau1 = beta[0]*delta1 + beta[1]*(delta1+exp1*frac1/tau[0]);
    tau2 = beta[2]*(delta2 + exp2*frac2/tau[1]);

    mat<Tfp> res(6,1);
    res.setEntry(0,0,beta0);
    res.setEntry(0,1,beta1);
    res.setEntry(0,2,beta2);
    res.setEntry(0,3,beta3);
    res.setEntry(0,4,tau1);
    res.setEntry(0,5, tau2);

    return res;
}

template<class Tfp>
void rates<Tfp>::setLinearSquareErrorYieldCurveFrgBonds(vector<bond> &bondsFRG, const vector<Tfp> &prices,
                                                   const date &priceDate) {

    // Get cash flow matrix with payment dates
    mat<Tfp> cashFlowMatrix(getCashFlowMatrix(bondsFRG, priceDate));

    // Separate payment dates
    vector<Tfp> paymentDates(cashFlowMatrix.nCols());
    swap_ranges(cashFlowMatrix.beginMutable(),
         cashFlowMatrix.getIterator(1,0),
         paymentDates.begin());
    cashFlowMatrix.removeRow(0);

    // Create price vector
    mat<Tfp> priceVector(cashFlowMatrix.nRows(),1);
    for (int i=0; i<cashFlowMatrix.nRows(); i++)
        priceVector[i,0] = prices[i];

    // Allocate discount vector
    mat<Tfp> d;

    // Check if squared cash flow matrix is invertible
    mat<Tfp> squaredCashFlow = (cashFlowMatrix.T()).matmul(cashFlowMatrix);
    if (squaredCashFlow.det() != 0){
        mat<Tfp> squaredInverted = squaredCashFlow;
        squaredInverted.inv();
        d = squaredInverted.matmul(cashFlowMatrix.T()).matmul(priceVector);
        return;
    }

    // Project the price vector onto our image
    mat<Tfp> projectedPrices = cashFlowMatrix.project(priceVector);

    // Clean up cashflow matrix
    mat<Tfp> cashFlowReduced = cashFlowMatrix.T();
    vector<int> rmvIdx = cashFlowReduced.makeColumnIdependent();
    cashFlowReduced.transpose();


    // Remove excess bond prices
    { // projectedPricesSize should be accessible later.
    int rmvCount = 0;
    int projectedPricesSize = projectedPrices.size();
    for (int i = 0; i < projectedPricesSize; i++) {
        const auto it = find(rmvIdx.begin(), rmvIdx.end(), i);
        if (it == rmvIdx.end())
            continue;
        projectedPrices.removeRow(i - rmvCount);
        rmvCount++;
    }
    }

    // Test rank
    mat<Tfp> myCopy = cashFlowReduced;
    mat<Tfp> myStats = myCopy.makeReducedEchelon();

    // Make tre tridiagonal moving average matrix
    mat<Tfp> matMA(cashFlowReduced.nCols(), cashFlowReduced.nCols(), Tfp(0.));
    matMA[0,0] = 1.1;
    matMA[0,1] = -1.;
    for (int i = 1; i< matMA.nRows(); i++) {
        Tfp normalizingConst(0.);
        //Tfp mat = priceDate.yearsuntil(paymentDates[i]);
        for (int j = 0 ; j< matMA.nCols(); j++)
            normalizingConst += (i==j) ? Tfp(0.) :   exp( -abs((paymentDates[i]- paymentDates[j])/365.));

        for (int j = 0 ; j< matMA.nCols(); j++)
            matMA[i,j] = (i==j) ? normalizingConst : -exp( -abs((paymentDates[i]- paymentDates[j])/365.));
    }

    //for (int j = 0 ; j< matMA.nCols(); j++)
    //    matMA[matMA.nCols()-1,j] = 0.1;

    //matMA[matMA.nCols()-1,matMA.nCols()-2] = -1.;
    //matMA[matMA.nCols()-1,matMA.nCols()-1] = 1.1;

    mat<Tfp> matMASqrdInv = matMA.T().matmul(matMA);
    matMASqrdInv.inv();
    mat<Tfp> bigInversion = cashFlowReduced.matmul(matMASqrdInv).matmul(cashFlowReduced.T());
    bigInversion.inv();


    mat<Tfp> dHat = matMASqrdInv.matmul(cashFlowReduced.T()).matmul(bigInversion).matmul(projectedPrices);


    for (int i =0; i< dHat.size(); i++){
        //cout << paymentDates[i] -Tfp(priceDate) << ":" << *dHat[i] << endl;
        cout << paymentDates[i] -Tfp(num(priceDate)) << ":" << round(*dHat[i]*10000) << endl;
    }


    mat<Tfp> priceDiff = cashFlowMatrix.matmul(dHat) - priceVector;

    // Set parameter
    paramsLagrange = mat<Tfp>(2, paymentDates.size());
    if constexpr (is_same_v<Tfp, adjointIntegral>){
        transform(paymentDates.begin(), paymentDates.begin()+paymentDates.size(),
                  paramsLagrange.beginMutable(),
                  [&priceDate](const Tfp &a){return priceDate.yearsuntil(a.getVal());});
    } else {
        transform(paymentDates.begin(), paymentDates.begin()+paymentDates.size(),
                  paramsLagrange.beginMutable(),
                  [&priceDate](const Tfp &a){return priceDate.yearsuntil(a);});
    }
    transform(dHat.begin(), dHat.end(),
              paramsLagrange.begin(),
              paramsLagrange.beginMutable()+paymentDates.size(),
              [](Tfp const &factor, const Tfp&mat){return -log(factor)/mat;});

    // Sort params
    auto datesBegin = paramsLagrange.beginMutable();
    auto datesEnd = datesBegin + paramsLagrange.nCols();
    int count = 0;
    while (count < 1000000){
        auto dateToSwap = adjacent_find(datesBegin, datesEnd, std::greater<>());
        if (dateToSwap == datesEnd) break;
        auto factorToSwap = dateToSwap + paramsLagrange.nCols();
        swap (*dateToSwap, *(dateToSwap+1));
        swap(*factorToSwap, *(factorToSwap+1));
        count++;
    }

    // Check if params are sorted
    if (count >= 1000000)
        throw runtime_error("Parameters for lagrange yield curve model did not get set due to sorting errors.");

    // Notify all of calibration
    paramsIsSet = true;
    modelIsCalibrated = true;
    yieldCurveModelType = "Lagrange";

}

template<class Tfp>
void rates<Tfp>::setSparseLinearSquareErrorYieldCurveFrgBonds(vector<bond> &bondsFRG, const vector<Tfp> &prices,
                                                   const date &priceDate) {

    // Get cash flow matrix with payment dates
    mat<Tfp> cashFlowMatrix(getCashFlowMatrix(bondsFRG, priceDate));

    // Separate payment dates
    vector<Tfp> paymentDates(cashFlowMatrix.nCols());
    swap_ranges(cashFlowMatrix.beginMutable(),
                cashFlowMatrix.getIterator(1,0),
                paymentDates.begin());
    cashFlowMatrix.removeRow(0);

    // Create projection of yearly discount factors
    vector<int> discountdates = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30};
    mat<Tfp> dP(cashFlowMatrix.nCols(), discountdates.size());
    for (int i=0; i<dP.nRows(); i++){
        for (int j=0; j<dP.nCols();j++){
            // Find two closest discount factor times and then find linear combination
            date myDate = date(to_num(paymentDates[i]));
            date myMaturityDate = priceDate + years(discountdates[j]);
            Tfp delta(myDate.yearsuntil(myMaturityDate));

            // Only proceed at first index greater than current date
            if (myDate >= myMaturityDate)
                continue;

            // Now write the linea combination
            dP[i,j] = abs(delta);

            // Repeat for the prior date
            myDate = date(to_num(paymentDates[i]));
            myMaturityDate = priceDate + years(discountdates[j-1]);
            delta = myDate.yearsuntil(myMaturityDate);
            dP[i,j-1] = abs(delta);

            // Normalize linear combination
            Tfp normalizingconst = dP[i,j]+ dP[i,j-1];
            dP[i,j] = (normalizingconst - dP[i,j])/normalizingconst;
            dP[i,j-1] = (normalizingconst - dP[i,j-1])/normalizingconst;
            break;

        }
    }

    // Apply sparse discount matrix
    cashFlowMatrix = cashFlowMatrix.matmul(dP);

    // Create price vector
    mat<Tfp> priceVector(cashFlowMatrix.nRows(),1);
    for (int i=0; i<cashFlowMatrix.nRows(); i++)
        priceVector[i,0] = prices[i];

    // Allocate discount vector
    mat<Tfp> d;

    // Check if squared cash flow matrix is invertible
    mat<Tfp> squaredCashFlow = (cashFlowMatrix.T()).matmul(cashFlowMatrix);

    mat<Tfp> squaredInverted = squaredCashFlow;
    squaredInverted.inv();
    d = squaredInverted.matmul(cashFlowMatrix.T()).matmul(priceVector);


    for (int i =0; i< d.size(); i++){
        //cout << paymentDates[i] -Tfp(priceDate) << ":" << *dHat[i] << endl;
        cout << paymentDates[i] -Tfp(num(priceDate)) << ":" << round(*d[i]*10000) << endl;
    }

    mat<Tfp> priceDiff = cashFlowMatrix.matmul(d) - priceVector;

    // Get all discount factors
    d = dP.matmul(d);

    /*for (int i =0; i< d.size(); i++){
        //cout << paymentDates[i] -Tfp(priceDate) << ":" << *dHat[i] << endl;
        cout << paymentDates[i] -Tfp(priceDate) << ":" << round(*d[i]*10000) << endl;
    }*/


    // Set parameter
    paramsLagrange = mat<Tfp>(2, paymentDates.size());
    transform(paymentDates.begin(), paymentDates.begin()+paymentDates.size(),
              paramsLagrange.beginMutable(),
              [&priceDate](const Tfp &a){return priceDate.yearsuntil(to_num(a));});
    transform(d.begin(), d.end(),
              paramsLagrange.begin(),
              paramsLagrange.beginMutable()+paymentDates.size(),
              [](Tfp const &factor, const Tfp&mat){return -log(factor)/mat;});

    // Sort params
    auto datesBegin = paramsLagrange.beginMutable();
    auto datesEnd = datesBegin + paramsLagrange.nCols();
    int count = 0;
    while (count < 1000000){
        auto dateToSwap = adjacent_find(datesBegin, datesEnd, std::greater<>());
        if (dateToSwap == datesEnd) break;
        auto factorToSwap = dateToSwap + paramsLagrange.nCols();
        swap (*dateToSwap, *(dateToSwap+1));
        swap(*factorToSwap, *(factorToSwap+1));
        count++;
    }

    // Check if params are sorted
    if (count >= 1000000)
        throw runtime_error("Parameters for lagrange yield curve model did not get set due to sorting errors.");

    // Notify all of calibration
    paramsIsSet = true;
    modelIsCalibrated = true;
    yieldCurveModelType = "Lagrange";

}

template class rates<float>;
template class rates<double>;
template class rates<adjointIntegral>;
