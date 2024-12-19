//
// Created by Jonas Wolff on 21/11/2024.
//

#ifndef CODE_SHORTRATEHWEV_H
#define CODE_SHORTRATEHWEV_H

#include "shortrate.h"
#include "termstructureStore.h"
#include "termstructureForwardHermite.h"
#include "quadrature.h"
#include "distributionFunctions.h"
#include "adjointNumeric.h"
#include <functional>
#include "optimizer.h"
#include "sobol.h"
#include <random>
#include "shortrateBachelier.h"
#include "ois.h"



// HULL-WHITE EXTENDED VASICEK
template<class T1, class termstructure>
class shortrateHWEV: public shortrate<shortrateHWEV<T1, termstructure>>{
private:

    //termstructureForwardHermite<T1> mytermstructure;
    termstructure mytermstructure;
    T1 sigma;
    T1 a;
    adjointIntegral sigmaAdj;
    adjointIntegral aAdj;
    sobol1D mySobol;
    std::random_device rd;
    std::mt19937 gen;
    uniform_real_distribution<> myMersenne;
    num dVoldBach = 0;
    num drstardVol = 0;

public:

    shortrateHWEV(termstructure myTerm): mytermstructure(myTerm), gen(std::mt19937(rd())), myMersenne(0,1), mySobol(){};

    // Construct adjoint from numeric
    template<typename = enable_if<is_same_v<T1, adjointIntegral> and is_same_v<termstructure,termstructureForwardHermite<adjointIntegral>>>>
    explicit shortrateHWEV(const shortrateHWEV<num, termstructureForwardHermite<num>> &rhsSRM) :
                mytermstructure(rhsSRM.getTermStructure()),
                gen(std::mt19937(rd())),
                myMersenne(0,1),
                mySobol(),
                sigma(rhsSRM.getVol(0.1)),
                a(rhsSRM.getMeanReversion(0.1))
                {
        static_assert(is_same_v<T1, adjointIntegral> and is_same_v<termstructure,termstructureForwardHermite<adjointIntegral>>);
    };

    template<typename = enable_if<is_same_v<T1, adjointIntegral> and is_same_v<termstructure,termstructureForwardHermite<adjointIntegral>>>>
    explicit shortrateHWEV(const shortrateHWEV<num, termstructureForwardHermite<num>> &rhsSRM, vector<pair<bond, num>> &bondData, date priceDate) :
            mytermstructure(rhsSRM.getTermStructure(), bondData, priceDate),
            gen(std::mt19937(rd())),
            myMersenne(0,1),
            mySobol(),
            sigma(rhsSRM.getVol(0.1)),
            a(rhsSRM.getMeanReversion(0.1))
    {
        static_assert(is_same_v<T1, adjointIntegral> and is_same_v<termstructure,termstructureForwardHermite<adjointIntegral>>);
    };

    ~shortrateHWEV() = default;

    shortrateHWEV<T1,termstructureForwardHermite<T1>>& operator=(const shortrateHWEV<T1,termstructureForwardHermite<T1>>&rhsSRM){
        sigma = rhsSRM.getVol(0.1);
        a = rhsSRM.getMeanReversion(0.1);
        mytermstructure = rhsSRM.getTermStructure();
    }

    shortrateHWEV<T1,termstructureForwardHermite<T1>>& operator=(shortrateHWEV<T1,termstructureForwardHermite<T1>>&&rhsSRM){
        sigma = move(rhsSRM.getVol(0.1));
        a = move(rhsSRM.getMeanReversion(0.1));
        mytermstructure = move(rhsSRM.getTermStructure());
    };

    num getdVoldBach() const {
        return dVoldBach;
    }

    num getdrstardVol() const {
        return drstardVol;
    }

    termstructure getTermStructure() const {
        return mytermstructure;
    }

    termstructure& getTermStructureRef() {
        return mytermstructure;
    }

    void setTermStructure(const termstructureForwardHermite<T1>& rhsTS){
        mytermstructure = rhsTS;
    }

    void setVol(const T1 &myVol){
        sigma = abs(myVol);
    }

    const T1 getVol(const num &time) const{
        return abs(sigma);
    }

    T1& getVolRef(const num &time) {
        return sigma;
    }

    void setMeanReversion(const T1 &mya){
        a=mya;
    }

    const T1& getMeanReversion(const num &time) const{
        return a;
    }

    T1& getMeanReversionRef(const num &time){
        return a;
    }

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    auto getmean(const timeFP &time){

        // Save as Expresion Template if needed
        auto myB = B(timeFP(0),time);
        auto myBT = BT(timeFP(0),time);
        auto sigma2 = sigma*sigma;
        auto g = sigma2/2.*myB*myB;
        auto gt = sigma2*myB*myBT;
        auto fT = mytermstructure.interpolateForwardDerivative(T1(time));
        auto f = mytermstructure.interpolateForward(T1(time));
        return fT+gt+a*(f+g);
    }

    num getmeannum(const num &time){
        // Save as Expresion Template if needed
        num myB = Bnum(0,time);
        num myBT = BTnum(0,time);
        num sigma2;
        num numa;
        if constexpr (is_same_v<T1, adjointIntegral>){
            sigma2 = sigma.getVal();
            sigma2 *= sigma2;
            numa = a.getVal();
        } else {
            sigma2 = sigma * sigma;
            numa = a;
        }
        num g = sigma2/2.*myB*myB;
        num gt = sigma2*myB*myBT;
        num fT, f;
        if constexpr (is_same_v<T1, adjointIntegral>){
            fT = mytermstructure.interpolateForwardDerivative(T1(time)).getVal();
            f = mytermstructure.interpolateForward(T1(time)).getVal();
        } else {
            fT = mytermstructure.interpolateForwardDerivative(time);
            f = mytermstructure.interpolateForward(time);
        }

        return fT+gt+numa*(f+g);
    }

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    auto getmeanAdj(const timeFP &time){
        // Save as Expresion Template if needed
        auto myB = BAdj(0,time);
        auto myBT = BTAdj(0,time);
        auto sigma2 = sigmaAdj*sigmaAdj;
        auto a2 = aAdj*aAdj;
        auto g = sigma2/2.*myB*myB;
        auto gt = sigma2*myB*myBT;
        auto fT = mytermstructure.interpolateForwardDerivative(time);
        auto f = mytermstructure.interpolateForward(time);
        return fT+gt+aAdj*(f+g);
    }

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    auto B(const timeFP &t, const timeFP &T){
        return 1./a*(1.-exp(-a*(T-t)));
    }

    num Bnum(const num &t, const num &T){
        if constexpr (is_same_v<T1, adjointIntegral>)
            return 1./a.getVal()*(1.-exp(-a.getVal()*(T-t)));
        else
            return 1./a*(1.-exp(-a*(T-t)));
    }

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    auto BT(const timeFP &t, const timeFP &T){
        return exp(-a*(T-t));
    }

    num BTnum(const num &t, const num &T){
        if constexpr (is_same_v<T1, adjointIntegral>)
            return exp(-a.getVal()*(T-t));
        else
            return exp(-a*(T-t));
    }

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    auto BAdj(const timeFP &t, const timeFP &T){
        return 1./aAdj*(1.-exp(-aAdj*(T-t)));
    }

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    auto BTAdj(const timeFP &t, const timeFP &T){
        return exp(-aAdj*(T-t));
    }

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    auto A(const timeFP&t, const timeFP &T){
        return adaptQuadrature([this,&T](const num &s){
            T1 myB = B(s,T);
            return .5*sigma*sigma*myB*myB-getmean(s)*myB; }, t, T, 100000000, numEPS*2*2);
    };

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    auto ZCB(const timeFP &t, const timeFP &T, const T1 &r){
        //return exp(A(t,T)-B(t,T)*r);

        auto ratio = mytermstructure.interpolateDiscount(T)/((t>0)?mytermstructure.interpolateDiscount(t):timeFP(1.));
        auto myB = B(t,T);
        auto arg = myB*mytermstructure.interpolateForward(t)-sigma*sigma/(4.*a)*myB*myB*(1.-exp(-2.*a*t))-myB*r;
        return ratio*exp(arg);
    }

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    auto ZCBAdj(const timeFP &t, const timeFP &T, const T1 &r){
        auto ratio = mytermstructure.interpolateDiscount(T)/((t>0)?mytermstructure.interpolateDiscount(t):timeFP(1.));
        auto myB = BAdj(t,T);
        auto arg = myB*mytermstructure.interpolateForward(t)-sigmaAdj*sigmaAdj/(4.*aAdj)*myB*myB*(1.-exp(-2.*aAdj*t))-myB*r;
        return ratio*exp(arg);
    }

    template<typename = enable_if<!is_same_v<T1, adjointIntegral>>>
    auto ZCBAdj(const num &t, const num &T, const adjointIntegral &r){
        auto ratio = mytermstructure.interpolateDiscount(T)/((t>0)?mytermstructure.interpolateDiscount(t):num(1.));
        auto myB = BAdj(t,T);
        auto arg = myB*mytermstructure.interpolateForward(t)-sigmaAdj*sigmaAdj/(4.*aAdj)*myB*myB*(1.-exp(-2.*aAdj*t))-myB*r;
        return ratio*exp(arg);
    }

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    auto call(const timeFP &optmat, const timeFP &bondmat, const T1 &K){
        auto temp = (1.-exp(-a*(bondmat-optmat)));
        auto SIGMA2 = sigma*sigma/(2.*a*a*a)*(1.-exp(-2.*a*optmat))*temp*temp;
        auto SIGMA = sqrt(SIGMA2);
        auto pT2 = mytermstructure.interpolateDiscount(bondmat);
        auto pT1 = mytermstructure.interpolateDiscount(optmat);
        auto d2 = log(pT2/pT1);
        d2 -= .5*SIGMA2;
        d2 /= SIGMA;
        auto d1 = d2 + SIGMA;
        return pT2*NormCDF(d1)-K*pT1*NormCDF(d2);
    }

    auto callAdj(const num &optmat, const num &bondmat, const adjointIntegral &K){
        auto temp = (1.-exp(-aAdj*(bondmat-optmat)));
        auto SIGMA2 = sigmaAdj*sigmaAdj/(2.*aAdj*aAdj*aAdj)*(1.-exp(-2.*aAdj*optmat))*temp*temp;
        auto SIGMA = sqrt(SIGMA2);
        auto pT2 = mytermstructure.interpolateDiscount(bondmat);
        auto pT1 = mytermstructure.interpolateDiscount(optmat);
        adjointIntegral d2(log(pT2/pT1));
        d2 -= .5*SIGMA2;
        d2 /= SIGMA;
        auto d1 = d2 + SIGMA;
        return pT2*NormCDF(d1)-K*pT1*NormCDF(d2);
    }

    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, dual<num>> or is_same_v<timeFP, T1>>>
    T1 swaption(const timeFP& optmat, const timeFP& swapmat, const int &floatfreq, const T1& strike){

        // Get timeline for cashflows
        size_t npaydays = int(to_num(swapmat-optmat))*floatfreq;
        vector<timeFP> paydays(npaydays, optmat);
        const timeFP matincrements = (swapmat - optmat)/npaydays;
        for (int i=0; i< npaydays; i++){
            paydays[i] += matincrements*(i+1);
        }


        // Find error function for r^*
        function<T1(const T1&)> errfunc = [this,&paydays,&strike, &floatfreq,&optmat, &swapmat](const T1& r)->T1{
            T1 res(0.);
            for (int i=0;i<paydays.size();i++)
                res += this->ZCB(optmat, paydays[i],r);
            res *= floatfreq;

            res *= strike;
            res += ZCB(optmat, swapmat,r) - 1.;
            return res;
        };

        // Solve r^*
        T1 rlow(-.9);
        T1 rhigh(.9);
        T1 rstar;
        try {
            rstar = bisection<T1>(move(errfunc), rlow, rhigh);
        } catch(exception e) {
            return NAN;
        }


        // Calculate swaption price (optmat=0 since we are discounting all the way back to t with tenor=0)
        // Strike is unchanged due to black scholes formula.
        T1 swaptionprice(0);
        for (int i=0;i<npaydays-1;i++){
            swaptionprice += call(optmat,paydays[i], ZCB(optmat,paydays[i],rstar));
        }
        T1 lastcall = call(optmat, paydays[npaydays-1],ZCB(optmat,paydays[npaydays-1],rstar));
        swaptionprice += lastcall;
        swaptionprice *= floatfreq;
        swaptionprice *= strike;
        swaptionprice += lastcall;

        T1 At = 0;
        for (int i=0; i< npaydays; i++) {
            At +=  mytermstructure.interpolateDiscount(paydays[i]); //ZCB(t,paydays[i], r0);
        }
        At *= 1./floatfreq;

        return swaptionprice;//*At;


    }

    template<typename = enable_if<is_same_v<T1, adjointIntegral>>>
    adjointIntegral swaptionAdj(const num& optmat, const num& swapmat, const int &floatfreq, const num& strike,
                                adjointIntegral& rstarAdj){

        // Get timeline for cashflows
        size_t npaydays = int(to_num(swapmat-optmat))*floatfreq;
        vector<num> paydays(npaydays, optmat);
        const num matincrements = (swapmat - optmat)/npaydays;
        for (int i=0; i< npaydays; i++){
            paydays[i] += matincrements*(i+1);
        }

        // Find error function for r^*
        function<num(const num&)> errfunc = [this,&paydays,&strike, &floatfreq,&optmat, &swapmat](const num& r)->num{
            num res(0.);
            for (int i=0;i<paydays.size();i++)
                res += this->ZCB(optmat, paydays[i],r);
            res *= floatfreq;
            res *= strike;
            res += ZCB(optmat, swapmat,r) - 1.;
            return res;
        };

        // Find error function for r^*
        function<adjointIntegral(adjointIntegral&)> errfuncAdj = [this,&paydays,&strike, &floatfreq,&optmat, &swapmat](adjointIntegral& r)->adjointIntegral{
            adjointIntegral res(0.);
            for (int i=0;i<paydays.size();i++)
                res += this->ZCBAdj(optmat, paydays[i],r);
            res *= floatfreq;
            res *= strike;
            res += ZCBAdj(optmat, swapmat,r) - 1.;
            return res;
        };

        // Solve r^*
        num rlow(-.9);
        num rhigh(.9);
        T1 rstar;
        try {
            rstar = bisection<T1>(move(errfunc), rlow, rhigh);
            num oldsigmaAdj = sigmaAdj.getAdjoint();
            num oldaAdj = aAdj.getAdjoint();
            sigmaAdj.setAdjoint(0.);
            aAdj.setAdjoint(0.);
            adjointIntegral myrstarAdj(rstar);
            myrstarAdj.setPause();
            adjointIntegral errAdj =  errfuncAdj(myrstarAdj);
            errAdj.setAdjoint(1.);
            errAdj.propagateToPause();
            drstardVol = -sigmaAdj.getAdjoint()/myrstarAdj.getAdjoint();
            cout << "drstar/dVol=" << drstardVol << endl;
            // Recreate parameters
            sigmaAdj.setAdjoint(oldsigmaAdj);
            aAdj.setAdjoint(oldaAdj);
        } catch(exception &e) {
            throw runtime_error("Bisection did not converge");
        }

        // Calculate swaption price (optmat=0 since we are discounting all the way back to t with tenor=0)
        // Strike is unchanged due to black scholes formula.
        adjointIntegral swaptionprice(0);
        rstarAdj = rstar;
        for (int i=0;i<npaydays-1;i++){
            swaptionprice += callAdj(optmat,paydays[i], ZCBAdj(optmat,paydays[i],rstarAdj));
        }
        adjointIntegral lastcall = callAdj(optmat, paydays[npaydays-1],ZCBAdj(optmat,paydays[npaydays-1],rstarAdj));
        swaptionprice += lastcall;
        swaptionprice *= floatfreq;
        swaptionprice *= strike;
        swaptionprice += lastcall;

        T1 At = 0;
        for (int i=0; i< npaydays; i++) {
            At +=  mytermstructure.interpolateDiscount(paydays[i]); //ZCB(t,paydays[i], r0);
        }
        At *= 1./floatfreq;

        return swaptionprice;//*At;


    }

    T1 calibrate(const vector<pair<bond, T1>> &data, const T1& swaptionbachvol, const T1& swaprate, const date &pricingDate,
                 const date &swaptionmat, const date&swapmat, const int &floatfreq, bool calcVolSinsitivity = false){

        // Extract bonds and prices
        vector<bond> myBonds;
        myBonds.reserve(data.size());
        mat<T1> myPrices(data.size(),1);
        for (int i=0;i<data.size();i++) {
            myBonds.push_back(data[i].first);
            //myPrices[i,0] = data[i].second;
            myPrices[i,0] = mytermstructure.getBondPrice(data[i].first, pricingDate);
        }

        // Generate cashflow matrix
        mat<T1> CF = getCashFlowMatrix(myBonds, pricingDate);

        // Extract dates
        vector<T1> paydays(CF.nCols());
        CF.getRow(0,paydays);
        CF.removeRow(0);

        // Get current short rate. No division by 0.
        T1 r0 = mytermstructure.interpolateYield(.00001);

        // Get tenor for swaption dates
        num T0 = pricingDate.yearsuntil(swaptionmat);
        num TN = pricingDate.yearsuntil(swapmat);
        num t = 0.;

        // Get annuity
        ois<T1> myOIS(pricingDate);
        size_t npaydays = swaptionmat.yearsuntil(swapmat)*floatfreq;
        vector<num> swappaydays(npaydays);
        num optmat = pricingDate.yearsuntil(swaptionmat);
        for (int i=0; i< npaydays; i++){
            swappaydays[i] = optmat + swaptionmat.yearsuntil(swapmat)/npaydays*(i+1);
        }
        T1 At = 0;
        for (int i=0; i< npaydays; i++) {
            date toDate = pricingDate + years(int(floor(swappaydays[i]))) + days(int(364*(swappaydays[i]-floor(swappaydays[i]))));
            At += 1./floatfreq * myOIS.getDiscount(pricingDate, toDate); //ZCB(t,paydays[i], r0);
        }

        // Get strike from bachelier
        T1 swaptionBachPrice = bachelier::getSwaptionPrice(t,T0, swaprate,abs(swaptionbachvol), swaprate, At);

        // Make error function
        auto calcErrVal = [this, &T0, &floatfreq, &TN, &swaptionBachPrice, &swaprate, &t]
                (const mat<T1> params)->T1{
            // Set parameters
            a = *params.begin();
            sigma = *(params.begin()+1);

            // swaption misspricing
            T1 err = swaptionBachPrice - this->swaption(T0, TN, floatfreq,swaprate);

            return err*err;
        };

        auto calcErrFDGrad = [this, &floatfreq, &T0, &TN, &swaprate, &swaptionBachPrice]
                (const mat<T1> &params, const T1& fval, mat<T1> &grad)->void {
            // Set parameters
            a = *params.begin() + numSQRTEPS;
            sigma = *(params.begin()+1);

            // swaption misspricing
            T1 err = swaptionBachPrice - this->swaption(T0, TN, floatfreq,swaprate);

            // Calculate initial error
            grad[0,0] = (err*err - fval)/numSQRTEPS;

            // Set parameters
            a = *params.begin();
            sigma = *(params.begin()+1) + numSQRTEPS;

            // swaption misspricing
            err = swaptionBachPrice - this->swaption(T0, TN, floatfreq,swaprate);

            // Calculate initial error
            grad[1,0] = (err*err - fval)/numSQRTEPS;
        };


        // Set initial parameters
        mat<T1> params(2,1);
        params[0,0] = 0.3;//(a)? a:0.3; // a
        params[1,0] = 0.06;//(sigma)? sigma: 0.06; // sigma

        // Minimize error
        int its;
        T1 error;
        dfpmin<T1>(params, numEPS*2*2, its, error, calcErrVal,calcErrFDGrad);

        // Check error
        /*cout << "Hull White FD Error: " << error << endl;
        if (false){
            swaptionBachPrice = bachelier::getSwaptionPrice(t,T0, swaprate,abs(swaptionbachvol)+numSQRTSQRTEPS, swaprate, At);;
            mat<T1> oldparams = params;
            T1 oldsigma = sigma;
            T1 newerr;
            dfpmin<T1>(params, numEPS*2*2, its, newerr, calcErrVal,calcErrFDGrad);
            dVoldBach = (sigma-oldsigma)/numSQRTSQRTEPS;
            cout << "dVol/dBach=" << dVoldBach << endl;
            a = oldparams[0,0];
            sigma = oldsigma;
        }
        return error;*/

        // Get strike from bachelier
        if (calcVolSinsitivity){
            adjointIntegral swaptionbachvolAdj(swaptionbachvol);
            adjointIntegral swaptionBachPriceAdj = bachelier::getSwaptionPrice(t,T0, swaprate,abs(swaptionbachvolAdj), swaprate, At);;
            aAdj = adjointIntegral(a);
            sigmaAdj = adjointIntegral(sigma);
            adjointIntegral rstarAdj(0.);

            auto calcErrAdj = [this, &T0, &floatfreq, &TN, &swaptionBachPriceAdj, &swaprate, &t, &rstarAdj]
                    ()->adjointIntegral{

                // swaption misspricing
                adjointIntegral err = swaptionBachPriceAdj - swaptionAdj( T0, TN, floatfreq,swaprate,rstarAdj);

                return err*err;
            };

            adjointIntegral myERRERR = calcErrAdj();
            myERRERR.setAdjoint(1.);
            myERRERR.propagateAll();
            dVoldBach = -(swaptionbachvolAdj.getAdjoint())/(sigmaAdj.getAdjoint()+rstarAdj.getAdjoint()*drstardVol);
            cout << "dVol/dBach=" << dVoldBach << endl;
            aAdj.setAdjoint(0.);
            sigmaAdj.setAdjoint(0.);
        }


        return error;
    }

    T1 calibrateAAD(const vector<pair<bond, T1>> &data, const date &pricingDate){

        // Extract bonds and prices
        vector<bond> myBonds;
        myBonds.reserve(data.size());
        mat<T1> myPrices(data.size(),1);
        for (int i=0;i<data.size();i++) {
            myBonds.push_back(data[i].first);
            //myPrices[i,0] = data[i].second;
            myPrices[i,0] = mytermstructure.getBondPrice(data[i].first, pricingDate);
        }

        // Generate cashflow matrix
        mat<T1> CF = getCashFlowMatrix(myBonds, pricingDate);

        // Extract dates
        vector<T1> paydays(CF.nCols());
        CF.getRow(0,paydays);
        CF.removeRow(0);

        // Get current short rate. No division by 0.
        T1 r0 = mytermstructure.interpolateYield(.00001);

        // Make error function
        auto calcErrVal = [&CF, &myPrices, &paydays, &r0, this]
                (const mat<T1> params)->T1{
            // Set parameters
            a = *params.begin();
            sigma = *(params.begin()+1);

            // Calculate deviation from yield
            mat<T1> err(29,1);
            for (int i=1; i<30;i++)
                err[i-1,0] = abs(mytermstructure.interpolateForwardDerivative(i)- getmean(i));

            // Calculate initial error
            return err.norm();
        };

        // Make adjoint error function
        auto calcErrAdj = [&CF, &myPrices, &paydays, &r0, this]
                (mat<adjointIntegral> &params)->adjointIntegral{
            // Set parameters
            aAdj = *params.begin();
            sigmaAdj = *(params.begin()+1);

            // Calculate deviation from yield
            mat<adjointIntegral> err(29,1);
            for (int i=1; i<30;i++)
                err[i-1,0] = abs(mytermstructure.interpolateYield(i)- getmeanAdj(i));

            // Calculate initial error
            return err.norm();
        };

        // Set initial parameters
        mat<T1> params(2,1);
        params[0,0] = .5; // a
        params[1,0] = .04; // sigma

        // Minimize error
        int its =0;
        T1 error;
        dfpminAAD(params,numEPS*2*2,its,error,calcErrAdj,calcErrVal);

        // Check error
        cout << "Hull White AAD Error: " << error << endl;
        cout << its << endl;

        return error;
    }

    void fillRatesPathSobol(vector<pair<num, T1>> &data) {
        // row 1 is dates
        // row 2 is the rates returned

        // Set short rate to short end of yield curve
        T1 lastrate = mytermstructure.interpolateYield(0.001);
        // Last date is now
        T1 lastdate = 0;
        for (int i = 0; i < data.nCols(); i++) {
            num todate = data[i].first;
            T1 myB = B(lastdate, todate);
            T1 mean = lastrate * BT(lastdate, todate);
            mean += adaptQuadrature([this, lastdate, todate](num& s) {
                return BT(s, todate) * getmean(s);
                },0,todate, 100000000);
            T1 sig = abs(sigma * myB);
            data[i].second = mean + mySobol.getN() * sig;
            lastdate = todate;
        }

        return;

    }

    // Calculate one simulation of ZCBs
    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    void getZCBsimulSobol(vector<pair<timeFP, T1>>& data, const timeFP simuldate) {
        // Row 1 is dates
        // Row 2 is ZCB

        // Simulate rate
        T1 r0 = mytermstructure.interpolateYield(T1(0.0001));
        num myN = mySobol.getN();
        T1 myB = B(timeFP(0), simuldate);
        T1 mean = r0 * BT(T1(0), T1(simuldate));
        if constexpr (is_same_v<T1, adjointIntegral>)
            mean += adaptQuadrature<adjointIntegral>([this, simuldate](const adjointIntegral &s) -> adjointIntegral {
                return BT(s, adjointIntegral(simuldate)) * getmean(s);
            },adjointIntegral(0),adjointIntegral(simuldate),1000000000, numSQRTEPS);
        else
            mean += adaptQuadrature<num>([this, simuldate](const num& s)->num {
                return BTnum(s, simuldate) * getmeannum(s);
                },num(0),num(simuldate),100000000);
        T1 sig = abs(sigma * myB);
        T1 rt = myN * sig + mean;

        // Calculate ZCBs with simulated rt
        T1 myForward;
        timeFP maturity;
        for (int i = 0; i < data.size(); i++) {
            // Check that date makes sense
            maturity = data[i].first;
            //if (maturity <= simuldate)
            //    throw invalid_argument("Bond is matured at simulation time");

            // Bjorks formular
            data[i].second = mytermstructure.interpolateDiscount(T1(maturity)) / mytermstructure.interpolateDiscount(T1(simuldate));
            myB = B(simuldate, maturity);
            myForward = mytermstructure.interpolateForward(T1(simuldate));
            data[i].second *= exp(myB * myForward - sigma * sigma / (4. * a) * myB * myB * (1 - exp(-2. * a * simuldate)) - myB * rt);
        }

    }

    // Calculate one simulation of ZCBs (const incorrectness, myRNG is changed)
    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    void getZCBsimulMersenne(vector<pair<num, T1>>& data, const timeFP simuldate) {
        // Row 1 is dates
        // Row 2 is ZCB

        // Simulate rate
        T1 r0 = mytermstructure.interpolateYield(0.0001);
        num myN = NormCDFI(myMersenne(gen));
        T1 myB = B(0, simuldate);
        T1 mean = r0 * BT(0, simuldate);
        mean += adaptQuadrature([this, simuldate](const num& s) {
            return BT(s, simuldate) * getmean(s);
            },0, simuldate, 100000000);
        T1 sig = abs(sigma * myB);
        T1 rt = myN * sig + mean;

        // Calculate ZCBs with simulated rt
        T1 myForward;
        num maturity;
        for (int i = 0; i < data.size(); i++) {
            // Check that date makes sense
            maturity = data[i].first;
            if (maturity <= simuldate)
                throw invalid_argument("Bond is matured at simulation time");

            // Bjorks formular
            data[i].second = mytermstructure.interpolateDiscount(maturity) / mytermstructure.interpolateDiscount(simuldate);
            myB = B(simuldate, maturity);
            myForward = mytermstructure.interpolateForward(simuldate);
            data[i].second *= exp(myB * myForward - sigma * sigma / (4. * a) * myB * myB * (1 - exp(-2. * a * simuldate)) - myB * rt);
        }

    }

    // calculate on simulation of the bond prices at simulation date (const incorectness myRNG changed)
    template<class timeFP, typename = enable_if<is_same_v<timeFP, num> or is_same_v<timeFP, T1>>>
    void getBONDsimulSobol(const mat<T1>& CF, const vector<timeFP>& tenors, const timeFP simuldate, vector<T1> &bondprices) {

        // Extract cashflow
        //mat<num> CF = getCashFlowMatrix(mybonds, simuldate);
        // Get tenors
        //vector<num> tenors(CF.nCols());
        //CF.getRow(0, tenors);
        //CF.removeRow(0);
        // Set up structure for discount factors
        vector<pair<timeFP, T1>> df(CF.nCols());
        for (int i = 0; i < CF.nCols(); i++)
            df[i].first = tenors[i];

        // generate one simulation of random ZCB
        getZCBsimulSobol(df, simuldate);

        // Convert df to matrix
        mat<T1> mydf(df.size(), 1);
        for (int i = 0; i < df.size(); i++)
            mydf[i, 0] = df[i].second;

        // Calculate bond prices
        mat<T1> tempprices = CF.matmul(mydf);
        for (int i=0; i<tempprices.size(); i++)
            bondprices[i] = *(tempprices.begin()+i);
    }

    // calculate on simulation of the bond prices at simulation date (const incorectness myRNG changed)
    void getBONDsimulMersenne(const mat<T1>& CF, const vector<num>& tenors, const num simuldate, vector<T1>& bondprices) {

        // Extract cashflow
        //mat<num> CF = getCashFlowMatrix(mybonds, simuldate);

        // Get tenors
        //vector<num> tenors(CF.nCols());
        //CF.getRow(0, tenors);
        //CF.removeRow(0);

        // Set up structure for discount factors
        vector<pair<num, T1>> df(CF.nCols());
        for (int i = 0; i < CF.nCols(); i++)
            df[i].first = tenors[i];

        // generate one simulation of random ZCB
        getZCBsimulMersenne(df, simuldate);

        // Convert df to matrix
        mat<T1> mydf(df.size(), 1);
        for (int i = 0; i < df.size(); i++)
            mydf[i, 0] = df[i].second;

        // Calculate bond prices
        // Calculate bond prices
        mat<T1> tempprices = CF.matmul(mydf);
        for (int i=0; i<tempprices.size(); i++)
            bondprices[i] = *(tempprices.begin()+i);
    }

};

#endif //CODE_SHORTRATEHWEV_H
