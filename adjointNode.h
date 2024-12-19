//
// Created by Jonas Wolff on 16/11/2024.
// Inspired (strongly) from Antoine Savines AAD library

#ifndef CODE_ADJOINTNODE_H
#define CODE_ADJOINTNODE_H


#include <cstdio>
#include <stdexcept>
#include "settings.h"

class adjointNode {
private:

    // We do not store the value here since it can be discarded early. This node is saved on a blocklist.

    //
    num adjoint = 0;

    // We save them on blocklist
    num* partials;

    // Where adjoins of arguments live. We need double pointer since adjoints will not be adjacent in memory.
    num** argadjoints;

    // Number of argument
    const size_t n;

    // is leaf?
    bool isleaf;

    // do we pause here?
    bool pause=false;

public:

    // Default constructor
    adjointNode(): n(0){};

    // Created in forwards pass. Only function evaluation is known so far.
    adjointNode(const size_t &myN) : n(myN){};

    // getter
    const size_t& getN( )const { return n; };

    // getter
    num& getAdjoint() {return adjoint;};

    // setter
    void setAdjoint(const num &adj){adjoint = adj;};

    // getter
    num* getPartials() const {return partials;};

    // setter
    void setPartials(num* myPar){partials = myPar;};

    // setter
    template<size_t idxPartial>
    void setPartial(const num& par){partials[idxPartial] = par;};

    // getter
    num** getArgAdjoint(){return argadjoints;};

    void setArgAdjoint(num** myArg){argadjoints = myArg;};

    // setter
    template<size_t idxAdjoint>
    void setArgAdjoint(num* myArg){argadjoints[idxAdjoint] = myArg;};

    // getter
    bool getIsLeaf(){return isleaf;};

    // setter
    void setIsLeaf(bool myisleaf){isleaf=myisleaf;};

    // mark as stopping point
    void setPause(bool mypause){pause=mypause;};

    // do we pause
    bool getPause(){return pause;};

    // push adjoint down to arguments
    void pushAdjoints()const{
        // max arguments of 255
        for (char i = 0; i<n; i++) {
            if (*(argadjoints + i) and *(partials + i))
                **(argadjoints + i) += *(partials + i) * adjoint; // first * are adjascent, second * points to adjoint
        }
    }


};


#endif //CODE_ADJOINTNODE_H
