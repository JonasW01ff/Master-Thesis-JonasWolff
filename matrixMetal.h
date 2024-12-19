//
// Created by Jonas Wolff on 14/11/2024.
//

#ifndef CODE_MATRIXMETAL_H
#define CODE_MATRIXMETAL_H

#include <Accelerate/Accelerate.h>
#include "settings.h"


class matrixMetal {
private:
    bool isLoaded = false;

public:

    matrixMetal();
    ~matrixMetal() = default;

    bool getIsLoaded() const{
        return isLoaded;
    };

    void madd(const num* A, const num *B, num *C, size_t n);



};

static matrixMetal myMatrixMetal = matrixMetal();

#endif //CODE_MATRIXMETAL_H
