//
// Created by Jonas Wolff on 26/11/2024.
//

using namespace std;
using namespace matplot;

int main(){

    vector<double> soboldescrp;
    vector<double> stldescrp;
    vector<double> nsimul;

    for (int nSim=5;nSim<=1000; nSim+=1) {

        sobol1D rng = sobol1D();
        vector<vector<double>> sobolU(nSim, vector<double>(2));
        for (int i = 0; i < sobolU.size(); i++) {
            sobolU[i][0] = rng.getU<0>();
            sobolU[i][1] = rng.getU<1>();
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::random_device rd2;
        std::mt19937 gen2(rd2());
        uniform_real_distribution<> dis(0, 1);
        vector<vector<double>> stlU(nSim, vector<double>(2));
        for (int i = 0; i < stlU.size(); i++) {
            stlU[i][0] = dis(gen);
            stlU[i][1] = dis(gen2);
        }

        nsimul.push_back(nSim);
        soboldescrp.push_back(descrepency(sobolU));
        stldescrp.push_back(descrepency(stlU));

    }

    auto h = plot(nsimul, soboldescrp, nsimul, stldescrp, ":");
    ylim({0,0.1});
    xlabel("RNG sequence length");
    ylabel("Descrepency");
    ::matplot::legend({"Sobol", "STL: Mersennetwister"});
    save("DescrpencyComparison.png");
    return 1;
}