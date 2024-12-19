//
// Created by Jonas Wolff on 01/11/2024.
// but heavily inspired by Antoine Savine's ConcurrentQueue

#include "ThreadPool.h"

thread_local size_t ThreadPool::myThreadNumber = 0;
