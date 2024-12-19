//
// Created by Jonas Wolff on 23/11/2024.
//

#ifndef CODE_PLOT_H
#define CODE_PLOT_H
#include "vector"
#include "algorithm"

using namespace std;

vector<pair<num, int>> easyhist(vector<num> &mydata, int nbins){

    // get max min
    num low = *min_element(mydata.begin(), mydata.end());
    num high = *max_element(mydata.begin(), mydata.end());
    num stepsize = (high-low)/nbins;

    // Compute buckets
    vector<num> edges(nbins);
    for (int i=0; i<= edges.size(); i++){
        edges[i] = num(i+1)*stepsize + low;
    }

    // Count
    vector<num> tempdata=mydata;
    sort(tempdata.begin(), tempdata.end());
    int binidx = 0;
    vector<int> freqs(nbins);
    for (auto y: tempdata){
        while (edges[binidx]<=y && binidx < edges.size()){
            binidx++;
        }
        freqs[binidx]++;
    }

    // Create data structure
    vector<pair<num, int>> output(nbins);
    output[0].first = (low + edges[0])/2.;
    output[0].second = freqs[0];
    for (int i=1; i<nbins; i++){
        output[i].first = (edges[i-1]+edges[i])/2.;
        output[i].second = freqs[i];
    }

    return output;

}

#endif //CODE_PLOT_H
