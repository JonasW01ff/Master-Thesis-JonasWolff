//
// Created by Jonas Wolff on 21/11/2024.
//

#ifndef CODE_QUADRATURE_H
#define CODE_QUADRATURE_H

#include <Accelerate/Accelerate.h>
#include <functional>
#include "settings.h"
#include "adjointNumeric.h"

using namespace std;

// Simpson's Adaptive Quadrature Algorithm
template<class Tfp = num, typename = enable_if<is_floating_point_v<Tfp> or is_same_v<Tfp, adjointIntegral>>>
Tfp adaptQuadrature(std::function<Tfp(const Tfp&)> f, const Tfp& low, const Tfp& high,
                  const size_t n=100000, const num eps=sqrt(numeric_limits<Tfp>::epsilon())) {

    // Check fat fingers
    if (Tfp(n)<Tfp(1.) or Tfp(eps) > Tfp(1.))
        throw runtime_error("adaptQuadrature: Are you sure that exponent has the right sign?");

    // a = low, b= high
    Tfp delta = high-low;
    Tfp sigma(0);
    Tfp h=delta*.5;
    Tfp c = (high+low)/2.;
    size_t k = 1;
    Tfp abar = f(low);
    Tfp bbar = f(high);
    Tfp cbar = f(c);
    Tfp S = (abar+4.*cbar+bbar)*h/3., Sp, zbar, Spp, vbar;
    array<Tfp, 6> vk = {low,h, abar,cbar,bbar, S};
    array<Tfp, 6> vkm1;
    Tfp ybar;
    int its = 0;
    while(k<=n){
        h = vk[1]/2.;
        ybar = f(vk[0]+h);
        Sp = (vk[2]+4*ybar+vk[3])*h/3.;
        zbar = f(vk[0]+3.*h);
        Spp = (vk[3]+4.*zbar+vk[4])*h/3;
        if (abs(Sp+Spp-vk[5])<30.*eps*h/delta){
            sigma += Sp + Spp + (Sp + Spp - vk[5])/15.;
            k -= 1;
            swap(vk,vkm1);
            if (k<=0) {
                return sigma;
            }
        } else {
            if (k>=n)
                throw runtime_error("Adaptive Quadrature did not converge");
            vbar = vk[4];
            vkm1 = {vk[0], h, vk[2], ybar, vk[3], Sp};
            k += 1;
            vk = {vkm1[0] +num(2.)*h, h, vkm1[4], zbar, vbar, Spp};
        }
    }
}


// Simpson's Adaptive Quadrature Algorithm With Adjoint Differentiaiton
num adaptQuadratureDiff(std::function<adjointIntegral(const adjointIntegral&)> fadj,
                    std::function<num(const num&)> f,
                    const num& low, const num& high,
                    const num& lowD, const num &highD,
                    const size_t n=1000000, const num eps=numSQRTEPS) {
    // Implement Leibniz Integral rule
    num deriv(0.);

    // Interval end dependence
    if(highD != 0.)
        deriv += f(high)*highD;

    // Interval start dependence
    if(lowD != 0.)
        deriv -= f(low)*lowD;

    // Setup numeric function of derivative integrand
    auto fD = [&fadj](const num &x){
        adjointIntegral xadj(x);
        adjointIntegral fDres = fadj(xadj);
        fDres.setAdjoint(1.);
        xadj.setPause();
        fDres.propagateToPause();
        return xadj.getAdjoint();
    };

    // intrgrant dependence
    num intgranddep = adaptQuadrature<num>(fD,low, high, n, eps);
    deriv += intgranddep;
    return deriv;
}
#endif //CODE_QUADRATURE_H
