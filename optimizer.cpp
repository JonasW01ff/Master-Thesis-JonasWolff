//
// Created by Jonas Wolff on 24/09/2024.
//

#include "optimizer.h"

// Finite differences
template<class T1, class T2>
void func_fd(
        function<void(const T1 &xBeginIt,const T1 &xEndIt, T2 &yBeginIt, T2 &yEndIt)> f,
        T1 &xBeginIt,
        T1 &xEndIt,
        const T2 &yBeginIt, // Pre-calculated y value at x
        const T2 &yEndIt,
        T2 &yDeltaBegin, // Row denominated transposed Jacobian matrix
        const num &sqrteps
        )
{
    // Find underlying types
    typedef remove_reference_t<decltype(*xBeginIt)> TX;
    typedef remove_reference_t<decltype(*yBeginIt)> TY;

    // Find number of variables
    //size_t n = distance(xBeginIt, xEndIt);

    // Fund number of outputs
    size_t m = distance(yBeginIt, yEndIt);

    T2 yDeltaIt = yDeltaBegin;
    for (const T1 xIt = xBeginIt; xIt < xEndIt; advance(xIt,1)){

        // Increment of x
        TX delta_x = sqrteps*abs(*xIt);

        // Check if we hit machine epsilon
        if (!delta_x) delta_x = sqrteps;

        // Save the old x
        TX myX = *xIt;

        // Increment x
        *xIt += delta_x;

        // Calculate increment
        f(xBeginIt, xEndIt, yDeltaIt, yDeltaIt+m);

        // Subtract old value
        transform(yDeltaIt, yDeltaIt+m,yBeginIt,yDeltaIt, [delta_x](TY y1, TY y0){return (y1-y0)/delta_x;});

        // Go to next row
        advance(yDeltaIt, m);

        // Reset x
        *xIt = myX;

    }

}

/*template<class T1>
T1 heat<T1>::minimize(auto startIt, auto endIt) {

    // Get dimensions
    int m = 20;
    int n = distance(endIt, startIt);

    // Check iterators location
    if (n <= 0 or n>100000)
        throw runtime_error("Iterators look incorrect");

    mat<T1> U(n,m);

    // Populate with random numbers
    U.setEntries(rand());


}*/