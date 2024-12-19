//
// Created by Jonas Wolff on 01/11/2024.
// but heavily inspired by Antoine Savine's ConcurrentQueue

#ifndef CODE_CONCURRENTQUEUE_H
#define CODE_CONCURRENTQUEUE_H

#include <queue>
#include <mutex>
#include <condition_variable>

using namespace std;

template<class T>
class concurrentQueue {
private:

    queue<T> myQueue;
    mutable mutex myMutex;
    condition_variable myCV;
    bool myStop = false;

public:

    // Constructor
    concurrentQueue() = default;

    // Destructor
    ~concurrentQueue(){
        stop();
    };

    void stop(){
        // make new scope
        {
            lock_guard<mutex> lg(myMutex);
            myStop = true;
        }
        myCV.notify_all();
    }

    bool empty() const {
        // Acquire lock so that others can finish inserting
        lock_guard<mutex> lg(myMutex);

        return myQueue.empty();
    }

    void push(T myObj){

        // New scope
        {
            lock_guard<mutex> lg(myMutex);

            // Insert object
            myQueue.push(move(myObj));

        }

        // Notify one that queue has been changed
        myCV.notify_one();
    }

    bool pop(T &myObj){

        // Acquire mutex so that que does not get empty after empty check but before pop
        unique_lock<mutex> ul(myMutex);

        // Check if que is empty
        myCV.wait(ul,[this](){return !myQueue.empty() or myStop;});

        // Check if queue is stopped
        if (myStop) {
            ul.unlock();
            return false;
        }

        // pop returns RValue and gets moved to myObj. myQueue is now shorter
        myObj =  move(myQueue.front());
        myQueue.pop();

        // unlock
        ul.unlock();

        // Return that queue got poped
        return true;

    }

    bool popNoWait(T &myObj){

        // Acquire mutex so that que does not get empty after empty check but before pop
        unique_lock<mutex> ul(myMutex);

        // Check if queue is stopped
        if (myQueue.empty()) {
            ul.unlock();
            return false;
        }

        // pop returns RValue and gets moved to myObj. myQueue is now shorter
        myObj =  move(myQueue.front());
        myQueue.pop();

        // unlock
        ul.unlock();

        // Return that queue got poped
        return true;

    }

    void clear(){

        // Get lock
        lock_guard<mutex> lg(myMutex);

        // Copy and move idiom (maybe stay safe when do concurrency)
        concurrentQueue<T> newQueue;
        myQueue = move(newQueue);

    }

};


#endif //CODE_CONCURRENTQUEUE_H
