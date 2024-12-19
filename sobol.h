//
// Created by Jonas Wolff on 22/11/2024.
//

#ifndef CODE_SOBOL_H
#define CODE_SOBOL_H

#include <Accelerate/Accelerate.h>
#include <vector>
#include "distributionFunctions.h"

using namespace std;

// m = {1}, a = {}, s=1
constexpr auto gendirection1(){
    array<unsigned int, 32> m2;
    m2[0] = 1;
    for (unsigned int i=1; i<m2.size(); i++){
        m2[i] = m2[i-1] << 1; // Multiply by 2
        //m2[i] ^= m2[i-2] << 2; // Multiply by 2^2
        m2[i] ^= m2[i-1];
    }

    for (int i=0;i<m2.size();i++)
        m2[i] = m2[i] <<(32-i-1);

    return m2;
}

// last a is always 1, so they ormit its mention (a = 1000000) and we only have one entry.
// m = {1, 3}, a = {1}, s=2
constexpr auto gendirection2(){
    array<unsigned int, 32> m2;
    m2[0] = 1;
    m2[1] = 3;
    for (unsigned int i=2; i<m2.size(); i++){
        m2[i] = m2[i-1] << 1; // Multiply by 2
        m2[i] ^= m2[i-2] << 2; // Multiply by 2^2
        m2[i] ^= m2[i-2];

    }

    for (int i=0;i<m2.size();i++)
        m2[i] = m2[i] << (32-i-1);

    return m2;
}

// m = {1,3,1}, a = {1,0}, s=3
constexpr auto gendirection3(){
    constexpr unsigned int s = 3;
    array<bool, 2> a = {1,0};
    array<unsigned int, 32> m;
    m[0] = 1;
    m[1] = 3;
    m[2] = 1;
    for (unsigned int i=s; i<m.size(); i++){
        m[i] = m[i-s] << s;
        m[i] ^= m[i-s];
        for (unsigned int k=1; k<=a.size(); k++) {
            if (a[k - 1])
                m[i] ^= m[i - k] << k; // *2^k
        }

        /*
        m2[i] = m2[i-1] << 1; // Multiply by 2
        //m2[i] = 0* m2[i-2] * 2*2 * 3; // Multiply by 2^2
        m2[i] ^= m2[i-3] << 3; // Multiply by 23
        m2[i] ^= m2[i-3];*/
    }
    return m;
}

// Directions gotten from https://web.maths.unsw.edu.au/~fkuo/sobol/
constexpr auto direction1 = gendirection1();
constexpr auto direction2 = gendirection2();
//constexpr auto direction3 = gendirection3();
constexpr array<decltype(direction1)*, 2> directions = {&direction1, &direction2};


inline vector<num> SOBOLRNG1D(size_t nb){
    vector<num> res(nb);
    unsigned temp = 0;
    unsigned int idx = 0;
    unsigned int unit = 1;
    for (unsigned int i = 0;i<nb; i++){

        // Find first 1 bit in i from the right. idx is the index
        idx = 0;
        unit = 1;
        while (i&unit) {
            unit <<= 1;
            idx++;
        }

        // XOR the directional number
        temp ^= direction1[idx];

        // Normalise it to a unit interval
        res[i] = num(temp);
    }
    return res;

}

inline num descrepency(const vector<vector<num>> &seq) {

    // Set seed
    srand(0);

    int dim = seq[0].size();
    int samplesize = seq.size();

    num eps = 1e-10;
    num delta = 1;
    vector<num> d_rand(dim);
    int count = 0;
    int B_idx = 0;
    bool inB;
    num volB;
    num culmErr = 0;
    num currErr = 0;
    num worstErr = 0;
    while (delta > eps)
    {
        // Create box B
        for (int i=0; i< dim;i++)
            d_rand[i] = num(rand()) / num(RAND_MAX);

        // Calculate volume of B
        volB = 1.;
        for (int i=0; i<dim; i++){
            volB *= d_rand[i];
        }

        // Check how many is in box
        inB = true;
        count = 0;
        for (int i = 0; i < samplesize; i++) {
            inB = true;
            for (int j = 0; j < dim; j++)
                inB *= (seq[i][j] <= d_rand[j]);

            count += inB;
        }

        // Calculate error
        currErr = abs(num(count) / num(samplesize) - volB);
        if (culmErr > worstErr) {
            worstErr = currErr;
        }

        if (B_idx > 2) {
            // Calculate average error
            delta = abs(culmErr / (B_idx - 1) - (culmErr + currErr) / B_idx);
        }

        // accumulate error
        culmErr += currErr;

        // Increment idx
        B_idx++;
    }

    return worstErr;

}

static constexpr double M_2PI = M_PI*2;

static inline __attribute__((noinline)) num sinNOSIMD(num x) {
    return sin(x);
}

static inline __attribute__((noinline)) num cosNOSIMD(num x) {
    return cos(x);
}

class sobol1D {
private:
    array<unsigned int,directions.size()> lastunif = {}; // initializes to 0
    num spareN;
    bool isspareN = false;
    num lastN;
    num isAntithetic = false;
    unsigned int unit;
    unsigned int oneidx;
    array<unsigned int, directions.size()> i = {};
    const num norm = pow(2,-32);
public:

    // dim starts from 0
    template<size_t dim>
    num getU(){
        // Find first 1 bit in i from the right. idx is the index
        oneidx = 0;
        unit = 1;
        while (i[dim]&unit) {
            unit <<= 1;
            oneidx++;
        }

        decltype(direction1)& mydirection = *(directions[dim]);

        // XOR the directional number
        lastunif[dim] ^= mydirection[oneidx];

        // our index is incremented
        i[dim]++;

        // Normalise it to a unit interval
        return num(lastunif[dim])*norm;
    }

    num getNMoro(){
        num U = getU<0>();
        // Check for antithetic sampling
        if (isAntithetic) {
            isAntithetic = false;
            return -lastN;
        }
        isAntithetic = true;
        lastN = NormCDFI(U);
        return lastN;
    }

    num getNBoxMullerNOSIMD(){

        // Check for antithetic sampling
        if (isAntithetic) {
            isAntithetic = false;
            return -lastN;
        }
        isAntithetic = true;

        // Return last N if saved
        if (isspareN){
            isspareN = false;
            lastN = isspareN;
            return spareN;
        }

        isspareN = true;
        num U2 = getU<1>();
        num U1 = getU<0>();
        U2 *= M_2PI;
        num coef = sqrt(-log(U1)*2);
        spareN = coef*sinNOSIMD(U2);
        num res = coef * cosNOSIMD(U2);
        lastN = res;
        return res;
    }

    num getNBoxMuller(){

        // Check for antithetic sampling
        if (isAntithetic) {
            isAntithetic = false;
            return -lastN;
        }
        isAntithetic = true;

        // Return last N if saved
        if (isspareN){
            isspareN = false;
            lastN = spareN;
            return spareN;
        }

        isspareN = true;
        num U2 = getU<1>();
        num U1 = getU<0>();
        U2 *= M_2PI;
        num coef = sqrt(-log(U1)*2);
        spareN = coef*sin(U2);
        num res = coef * cos(U2);
        lastN = res;
        return res;
    }

    num getN(){
        return getNBoxMuller();
    }

};

#endif //CODE_SOBOL_H
