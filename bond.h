//
// Created by Jonas Wolff on 31/08/2024.
//

#ifndef CODE_BOND_H
#define CODE_BOND_H

#include <vector>
#include <algorithm>
#include "date.h"
#include "chrono"
#include "matrix.h"

using namespace std;
using namespace chrono;

class bond {

    const date issuanceDate;
    const date maturityDate;
    const date intrestStartDate;
    const date firstcouponDate;

    const num couponSize;
    const chrono::months couponFrequency;

    // daycountSchem may be:
    // 1 Actual/Actual (ICMA)
    const int daycountScheme;

    mat<num> cashflowMatrix;
    date cashflowMatrixDate;


public:

    bond();

    bond(date pissuanceDate, date pmaturityDate, date pintrestStartDate, date pfirstcouponDate,
         num pcouponSize, chrono::months pcouponFrequency, int pdaycountScheme) :
            issuanceDate(pissuanceDate), maturityDate(pmaturityDate),
            firstcouponDate(pfirstcouponDate), intrestStartDate(pintrestStartDate),
            couponSize(pcouponSize), couponFrequency(pcouponFrequency),
            daycountScheme(pdaycountScheme) {};


    bond (const bond& rhs):
            issuanceDate(rhs.issuanceDate), maturityDate(rhs.maturityDate),
            firstcouponDate(rhs.firstcouponDate), intrestStartDate(rhs.intrestStartDate),
            couponSize(rhs.couponSize), couponFrequency(rhs.couponFrequency),
            daycountScheme(rhs.daycountScheme) {};

    const date& getIssuanceDate() const {return issuanceDate;};

    const date& getMaturityDate() const {return maturityDate;};

    const date& getIntrestStartDate() const {return intrestStartDate;};

    const date& getFirstcouponDate() const {return firstcouponDate;};

    const num& getCouponSize() const {return couponSize;};

    const months& getCouponFrequency() const {return couponFrequency;};

    const int& getDaycountScheme() const {return daycountScheme;};

    const num getAcrruedIntrest(const date &priceDate) const;

    mat<num> getCashFlow(const date& priceDate);

    num conversionfactor(const date &delivery, const num futC) const;

};

mat<num> getCashFlowMatrix(vector<bond> &bondsFRG, const date &priceDate);

#endif //CODE_BOND_H
