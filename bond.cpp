//
// Created by Jonas Wolff on 31/08/2024.
//

#include "bond.h"

mat<num> bond::getCashFlow(const date& priceDate) {

    // Check if we already got cashflow matrix
    if (cashflowMatrixDate == priceDate)
        return cashflowMatrix;

    // Get cashflow length
    years yearsToMaturity = maturityDate.year() - priceDate.year();

    // Check maturity
    if (maturityDate <= priceDate)
        throw runtime_error("Matured bonds has no cashflow");


    // Initiate cash flow
    mat cashflow(2,yearsToMaturity.count()+1, couponSize);

    // Get specs
    const date& firstcouponDate = getFirstcouponDate();
    const date& maturityDate = getMaturityDate();
    const months& couponFrequency = getCouponFrequency();
    const num& couponSize = getCouponSize();

    // Return only principle if coupon is zero
    if (couponSize == 0){
        cashflow.resize(2,1);
        cashflow[0,0] = priceDate.yearsuntil(maturityDate);
        cashflow[1,0] = 100.;
        return cashflow;
    }

    // Get next coupon date
    date nextcouponDate = firstcouponDate;


    // Loop forward to present
    while (nextcouponDate <= priceDate)
        nextcouponDate += couponFrequency;

    vector<num>::iterator cashflowIt = cashflow.beginMutable();
    while (nextcouponDate <= maturityDate){
        // Set payment
        (*cashflowIt) = priceDate.yearsuntil(nextcouponDate);
        nextcouponDate += couponFrequency;
        cashflowIt++;
    }


    // Resize vector
    cashflow.resize(2, cashflowIt- cashflow.begin());

    // Add principle
    cashflowIt = cashflow.endMutable();
    cashflowIt--;
    (*cashflowIt) += 100;


    cashflowMatrix = cashflow;
    cashflowMatrixDate = priceDate;


    return cashflow;
}

const num bond::getAcrruedIntrest(const date &priceDate) const {
    // Get specs
    const date& firstcouponDate = getFirstcouponDate();
    const date& intrestStartDate = getIntrestStartDate();
    const date& maturityDate = getMaturityDate();
    const months& couponFrequency = getCouponFrequency();
    const num& couponSize = getCouponSize();

    // Get accrued interest
    num ai;
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

    return ai;
}

// https://www.eurex.com/ex-en/data/clearing-files/notified-deliverable-bonds-conversion-factors
num bond::conversionfactor(const date &DD, const num futc) const{

    // Euro denominated bonds

    // Find first coupon
    date NCD = maturityDate;
    while (NCD> DD)
        NCD -= years(1);
    if (NCD<=DD)
        NCD += years(1);

    // Definitions
    const date NCD1y = NCD+years(-1);
    const date NCD2y = NCD+years(-2);
    const date LCD = (NCD1y<firstcouponDate)?intrestStartDate:NCD1y;
    const num delta_e = (NCD1y-DD).count();
    const num act1 = ((delta_e<num(0))? NCD - NCD1y: NCD1y-LCD).count();
    const num delta_i = (NCD1y - LCD).count();
    const num act2 = ((delta_i<num(0))?NCD-NCD1y:NCD1y-NCD2y).count();
    const num f = 1.+num(delta_e)/act1;
    const num c = couponSize;
    int n = NCD.yearsuntil(maturityDate);

    num futc100 = futc/num(100.);
    num c100 = c/num(100.);
    num futc100p1 = 1.+futc100;
    num futc100p1INV = 1./futc100p1;
    num deltaEact1 = delta_e/act1;
    num deltaIact2 = delta_i/act2;

    num futc100p1INV_N = pow(futc100p1INV,n);
    num term1 = futc100p1-futc100p1INV_N;
    num term2 = futc100*deltaIact2+c/futc*term1+futc100p1INV_N;
    term2 *= pow(futc100p1INV, f);

    num res = term2 - c100*(deltaIact2-deltaEact1);
    return res;


}


mat<num> getCashFlowMatrix(vector<bond> &bondsFRG, const date &priceDate){

    // Im assuming we won't have more than 8 payments a year for 50 years
    mat<num> cashFlowMatrix(bondsFRG.size()+1, 8*50);

    vector<num>::const_iterator newEnd;
    vector<mat<num>> myCFS(bondsFRG.size());
    for (int i = 0; i< bondsFRG.size(); i++){

        mat<num> bondCashFlow = bondsFRG[i].getCashFlow(priceDate);
        for (int j=0;  j< bondCashFlow.nCols(); j++)
            cashFlowMatrix[i,j] = bondCashFlow[0,j];

    }

    // Get dates in right order
    sort(cashFlowMatrix.beginMutable(), cashFlowMatrix.endMutable());

    // Find unique dates
    newEnd = unique(cashFlowMatrix.beginMutable(), cashFlowMatrix.endMutable());

    // 0 will be first since smallest, so throw to the back
    rotate(cashFlowMatrix.beginMutable(),cashFlowMatrix.beginMutable()+1,cashFlowMatrix.endMutable());

    // Overwrite old dates
    vector<num>::iterator fillStart = cashFlowMatrix.getIterator(0, newEnd-cashFlowMatrix.begin());
    std::fill(fillStart, cashFlowMatrix.endMutable(), 0.0);

    // Resize timeline to the actual number of unique paydays
    cashFlowMatrix.resize(cashFlowMatrix.nRows(),newEnd-cashFlowMatrix.begin()-1);

    // Add cash flows now
    for (int i=0; i< bondsFRG.size(); i++){
        mat bondCashFlow = bondsFRG[i].getCashFlow(priceDate);
        int ii = 0;
        for (int j=0; j < cashFlowMatrix.nCols(); j++)
            if (cashFlowMatrix[0, j] == bondCashFlow[0,ii])
            {cashFlowMatrix[i+1,j] = bondCashFlow[1,ii]; ii++;}

        if (ii != bondCashFlow.nCols())
            throw runtime_error("Bond did not insert all payments to the cashflow matrix");

    }

    return cashFlowMatrix;

};