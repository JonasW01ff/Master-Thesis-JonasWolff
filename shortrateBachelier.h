//
// Created by Jonas Wolff on 06/12/2024.
//

#ifndef CODE_SHORTRATEBACHELIER_H
#define CODE_SHORTRATEBACHELIER_H

#include "optimizer.h"
#include "distributionFunctions.h"

namespace bachelier{

    template<class T>
    T getSwaptionPrice(const T &t, const T &T0,const T &St, const T&sigma, const T& K, const T& At){

        T Narg = (K-St)/(sigma*sqrt(T0-t));
        return ((K-St)* NormCDF(Narg)+sigma*sqrt(T0-t)* NormPDF(Narg))*At;
    }

    template<class T>
    adjointIntegral getSwaptionPrice(const T &t, const T &T0,const T &St, const adjointIntegral&sigma, const T& K, const T& At){

        auto Narg = (K-St)/(sigma*sqrt(T0-t));
        return ((K-St)* NormCDF(Narg)+sigma*sqrt(T0-t)* NormPDF(Narg))*At;
    }


}

#endif //CODE_SHORTRATEBACHELIER_H
