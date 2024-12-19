//
// Created by Jonas Wolff on 16/11/2024.
// Inspired (strongly) from Antoine Savines AAD library

#include "adjointTape.h"
#include "adjointNumeric.h"

// Start Tape
adjointTape myTape;
thread_local adjointTape* adjointIntegral::localTape = &myTape;

// When new thread starts, copy old tape
//thread_local adjointTape* localTape = myTape;