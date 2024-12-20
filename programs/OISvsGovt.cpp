//
// Created by Jonas Wolff on 06/12/2024.
//

using namespace matplot;

int main() {
    date priceDate(2024y / September / 2d);

    ois <num> estrcurve(priceDate);

    ThreadPool *mypool = &ThreadPool::getThreadPool();
    mypool->addMaxThreads();

    vector <bond> myBonds = getBondData();
    vector <num> myPrices = getBondPriceToday();
    cout << num(priceDate) << endl;

    // add accrued interest.
    for (int i = 0; i < myPrices.size(); i++) {
        myPrices[i] += myBonds[i].getAcrruedIntrest(priceDate);
    }

    // Discounted cashflows are dirty prices
    //myRatesLagrange.setLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, priceDate);
    //myRatesNelson.setNSSYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));
    //myRatesSparse.setSparseLinearSquareErrorYieldCurveFrgBonds(myBonds, myPrices, date(2024y,September, 2d));
    vector <pair<bond, num>> myBondData;
    for (int i = 0; i < myBonds.size(); i++)
        myBondData.push_back(make_pair(myBonds[i], myPrices[i]));

    // Start the timer
    auto start_time = std::chrono::high_resolution_clock::now();
    termstructureForwardHermite <num> myHermite("FRG", myBondData, priceDate);
    auto end_hermite_time = std::chrono::high_resolution_clock::now();


    vector <num> oisrates;
    vector <num> govtcurve;
    vector <num> oistenors;
    for (int i = 1; i < 365 * 11; i += 1) {
        date toDate = priceDate + days(i);
        num tenor = priceDate.yearsuntil(toDate);
        oistenors.push_back(tenor);
        oisrates.push_back(estrcurve.getEffectiveRate(priceDate, toDate));
        govtcurve.push_back(myHermite.interpolateYield(tenor));
        //cout << tenor << "," << estrcurve.getDiscount(priceDate, toDate) << endl;
    }

    plot(oistenors, oisrates, oistenors, govtcurve);
    title("OIS curve vs Govt Curve");
    ylabel("%");
    xlabel("tenor");
    text(9, 0.039, "Date: " + priceDate.tostring());
    ::matplot::legend({"OIS", "Govt"});
    save("OIScurve.png");
    return 1;
}