//
// Created by Jonas Wolff on 16/11/2024.
// Inspired (strongly) from Antoine Savines AAD library

#ifndef CODE_ADJOINTTAPE_H
#define CODE_ADJOINTTAPE_H

#include "adjointNode.h"
#include "blockList.h"

constexpr size_t BLOCKSIZE = 1000000;
constexpr size_t LISTSIZE = 4;

class adjointTape{
private:

    blockList<num*, BLOCKSIZE, LISTSIZE> poolArgAdjoints;
    blockList<num, BLOCKSIZE, LISTSIZE> poolPartials;
    blockList<adjointNode, BLOCKSIZE, LISTSIZE> poolNodes;

    friend class adjointIntegral;

public:

    template<size_t n> // n needs to be known at compile time.
    adjointNode* record(){

        // Make node
        adjointNode* myNode = poolNodes.emplace(n);

        // Allocate derivatives
        num* myPartials = poolPartials.emplaceNCopy<n>(0.);

        // Map derivatives
        myNode->setPartials(myPartials);

        // Allocate adjacent pointers to adjoints of arguments
        num** myArgAdjoints = poolArgAdjoints.emplaceNCopy<n>(nullptr);

        // Map adjoints of arguments
        myNode->setArgAdjoint(myArgAdjoints);

        return myNode;

    }

    adjointNode* begin(){
        return poolNodes.begin();
    }

    adjointNode* end(){
        return poolNodes.end();
    }

    adjointNode* back(){
        return poolNodes.back();
    }


};

#endif //CODE_ADJOINTTAPE_H
