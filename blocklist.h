//
// Created by Jonas Wolff on 13/11/2024.
//

#ifndef CODE_BLOCKLIST_H
#define CODE_BLOCKLIST_H
#include <cstddef>
#include <stdexcept>
#include <list>
#include <mutex>
#include <utility>
#include <array>
#include <iostream>

using namespace std;

template<class T, size_t myBlocksize, size_t myMaxblocks>
class blockList {
private:

    const size_t blocksize = myBlocksize;
    const size_t maxblocks = myMaxblocks;
    list<array<T, myBlocksize>> blocks;
    mutex myMutex;
    T* currObj = 0;
    T* startBlock = 0;
    T* endBlock = 0;
    T* startBlockList = 0;
    T* endBlockList = 0;
    size_t blockidx = 0;

    // Allocate Block
    T* allocateBlock() {

        // Increment blockidx
        blockidx++;

        cout << "Allcation block with current size=" << blocks.size() << endl;

        // Check that we have not reached maxblocks
        if (blocks.size() == myMaxblocks)
            throw runtime_error("Too many blocks have been allocated");

        if (blockidx < 0 or blockidx > blocks.size())
            throw runtime_error("Block idx inavlid");

        /*
        // Create block
        array<T, myBlocksize> myBlock;

        // Push back
        blocks.push_back(move(myBlock));*/
        // Quicker
        //cout << "Size is: " << distance(blocks.begin(),blocks.end()) << endl;
        if (blockidx >= blocks.size()) {
            cout << "Emplace back commencing: size=" << blocks.size() << "; idx=" << blockidx << endl;
            blocks.emplace_back();
            //if (distance(blocks.end(), blocks.begin()))
            //    throw runtime_error("Emplace back on list in blocklist made corruption: Increase block size");
            cout << "Emplace back done: size=" << blocks.size() << "; idx=" << blockidx << endl;
        }

        // Get pointers
        auto blockIt = blocks.begin();
        //cout << "Advancing n_steps: " << blockidx << endl;
        //cout << "Size is: " << distance(blocks.begin(),blocks.end()) << endl;
        advance(blockIt,blockidx);
        //cout << "Block advanced" << endl;
        startBlock = blockIt->data(); //data(blocks.back());
        currObj = startBlock;
        endBlock = startBlock + myBlocksize;
        if (endBlockList<endBlock)
            endBlockList = endBlock;
        cout << "Block allocated: idx=" << blockidx << endl;
        return currObj;
    }

    // Deallocate Block
    T* deallocateBlock() {

        cout << "Deallocating block with current size " << blocks.size() << endl;

        // Decrement block index
        blockidx--;

        // Check if is empty
        if (blocks.size() <= 1)
            throw runtime_error("Cannot pop empty blocks or blocks of size 1");

        if (blockidx <0 or blockidx>=blocks.size())
            throw runtime_error("blockidx invalid");

        // deallocateBlock
        auto blockIt = blocks.begin();
        advance(blockIt, blockidx);
        endBlockList = blockIt->data() + myBlocksize; //data(blocks.back()) + myBlocksize
        //cout << "Blocksize: " << blockidx << endl;
        startBlock = blockIt->data();
        endBlock = endBlockList;
        currObj = endBlock;
        //cout << "Blocksize: " << blockidx << endl;
        return currObj;

    }

public:

    // Constructor
    blockList() {
        // Make first block
        blocks.emplace_back();
        startBlockList = blocks.front().data();
        endBlockList = blocks.back().data() + myBlocksize;
        startBlock = startBlockList;
        endBlock = endBlockList;
        currObj = startBlock;
        blockidx = 0;
    };

    // Destructor
    ~blockList() = default;


    // Allocate
    T* allocate(T& Obj) {

        {
            // RAII
            lock_guard<mutex> lg(myMutex);

            // Check if currObj is at end of array
            if (currObj == endBlock) {
                currObj = allocateBlock();
            }

            // Insert object
            *currObj = move(Obj);
        }

        return currObj++;
    }

    // Move currObj one Back
    T* deallocateOne() {
        // RAII
        lock_guard<mutex> lg(myMutex);

        if (currObj<startBlock or currObj >endBlock)
            throw runtime_error("WTF");

        // deallocateBlock returns back() + 1
        if (currObj > startBlock)
            return --currObj;
        while (currObj == startBlock)
            currObj = deallocateBlock();

        return --currObj;

    }

    // Get Prior
    T* getPrior() {

        // deallocateBlock returns back() + 1
        if (currObj > startBlock)
            return (currObj-1);

        if (blocks.size() <= 1)
            throw runtime_error("Cant move back from first object");

        auto blockIt = blocks.begin();
        advance(blockIt,blockidx-1);
        return &blockIt->back();

    }

    // Get prior
    T* deallocateToInc(T* myObj) {

        // RAII
        lock_guard<mutex> lg(myMutex);

        // Deallocate block until we in correct block
        while ((myObj < startBlock or myObj > endBlock) && blocks.size() != 1) {
            currObj = deallocateBlock();
        }
        // Check if block where found
        if (myObj > endBlock or myObj < startBlock)
            throw out_of_range("Object is not on the block list");

        // You cannot move back from first object
        //cout << "block size = " << distance(startBlock, currObj) << endl;
        if (myObj <= startBlockList){
            currObj = myObj;
            return nullptr;
        }


        currObj = myObj;
        return getPrior();

    }

    T* deallocateToStart() {

        // RAII
        lock_guard<mutex> lg(myMutex);

        // Deallocate block until we in correct block
        while ((startBlockList < startBlock or startBlockList > endBlock) && blocks.size() != 1) {
            currObj = deallocateBlock();
        }
        // Check if block where found
        if (startBlockList > endBlock or startBlockList < startBlock)
            throw out_of_range("Object is not on the block list");

        // You cannot move back from first object
        //cout << "block size = " << distance(startBlock, currObj) << endl;
        currObj = startBlockList;
        return currObj;

    }

    // Emplace
    template<class... Types>
    T* emplace(Types ...args) {

        {
            // RAII
            lock_guard<mutex> lg(myMutex);

            // Check if currObj is at end of array
            if (currObj == endBlock) {
                currObj = allocateBlock();
            }

            // Construct object in place
            new (currObj)T(std::forward<Types>(args)...);
        }

        T* res = currObj;
        currObj++;

        return res;

    }

    template<size_t N, class... Types>
    T* emplaceNCopy(Types ...args) {

        // Check if currObj is at end of array
        if (currObj == endBlock) {
            currObj = allocateBlock();
        }

        // Save position
        T* res = currObj;

        {
            // RAII
            lock_guard<mutex> lg(myMutex);

            for (int i = 0; i < N; i++) {
                // Check if currObj is at end of array
                if (currObj == endBlock) {
                    currObj = allocateBlock();
                }

                // Construct object in place
                new(currObj)T(std::forward<Types>(args)...);
                currObj++;
            }
        }


        return res;

    }

    template<size_t N, class... Types>
    T* emplaceNNone() {
        // RAII
        lock_guard<mutex> lg(myMutex);

        for (int i = 0; i < N; i++) {
            // Check if currObj is at end of array
            if (currObj == endBlock) {
                currObj = allocateBlock();
            }

            // Memory is already allocated
            currObj++;
        }

        return currObj;

    }

    // Deallocate
    T* deallocateTo(T* myMark) {

        {

            // RAII
            lock_guard<mutex> lg(myMutex);

            // Deallocate blocks
            while (myMark < startBlock or myMark > endBlock)
                deallocateBlock();

            // Set currObj to mark
            currObj = myMark;
        }

        return currObj;

    }

    // get start
    T* begin() {
        return startBlockList;
    }

    // get end
    T* end() {
        // currObj is always positioned at the next free piece of memory and thus is the end position
        return currObj;
    }

    // get back
    T* back() {
        return currObj - 1;
    }

};


#endif //CODE_BLOCKLIST_H