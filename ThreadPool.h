//
// Created by Jonas Wolff on 01/11/2024.
// but heavily inspired by Antoine Savine's ConcurrentQueue

#ifndef CODE_THREADPOOL_H
#define CODE_THREADPOOL_H

#include <functional>
#include <future>
#include <thread>
#include <chrono>

#include "concurrentQueue.h"

using namespace std;

class ThreadPool {

    vector<thread> myThreads;

    concurrentQueue<packaged_task<void()>> myTaskQueue;

    thread_local static size_t myThreadNumber;


    // Race condyion may happen but 2 threads stopping a que will still stop it
    // Thou if ever changed we might make it atomic
    bool stopPool = false;

    // Make threadpool a singleton
    ThreadPool(){};

    void myWorkerFunction(int threadNumber){
        myThreadNumber = threadNumber;

        packaged_task<void()> myTask;

        while(!stopPool) {
            if(myTaskQueue.pop(myTask))
                myTask();
        }
    }

public:

    // make threadpool a singleton
    // CSTR
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool(ThreadPool &&) = delete;
    // ASGN
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool& operator=(ThreadPool &&) = delete;

    // Construct Threadpool (Meyers Singleton Pattern)
    static ThreadPool& getThreadPool(){
        static ThreadPool GlobalThreadPool;
        return GlobalThreadPool;

    };

    ~ThreadPool(){
        stopPool = true;
        myTaskQueue.stop();
        for (auto& t: myThreads)
            t.join();
    };

    void addThreads(int n){
        if (n+myThreads.size()+1 >= thread::hardware_concurrency())
            throw runtime_error("Too many threads");

        for (int i = 0; i < n ; i++){
            myThreads.emplace_back(&ThreadPool::myWorkerFunction, this, myThreads.size()+1);
        }

    };

    void addMaxThreads(){
        for (int i = myThreads.size(); i < thread::hardware_concurrency()-1 ; i++){
            myThreads.emplace_back(&ThreadPool::myWorkerFunction, this, myThreads.size()+1);
        }

    };

    void helpQueueWhile(future<void> &fut){
        packaged_task<void()> myTask;
        while (fut.wait_for(0s) != future_status::ready){
            if(!myTaskQueue.popNoWait(myTask))
                break;
            myTask();
        }
        fut.get();

    }


    future<void> spawnTask(function<void()> &&myCallable){

        packaged_task<void()> myTask(move(myCallable));

        future<void> myFuture = myTask.get_future();

        myTaskQueue.push(move(myTask));

        return move(myFuture);
    };

    size_t getThreadIdx() const{return myThreadNumber;};

};



#endif //CODE_THREADPOOL_H
