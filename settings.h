#pragma once
#include <limits>
#include <cmath>
typedef float num;
constexpr num numEPS = std::numeric_limits<num>::epsilon();
static num numSQRTEPS = sqrt(numEPS); // Made a mistake here but I have to live with it
static num numSQRTSQRTEPS = sqrt(sqrt(numEPS)); // Made a mistake here but I have to live with it
/*
typedef float num;
constexpr num numEPS = std::numeric_limits<num>::epsilon();
constexpr num numSQRTEPS = 1e-4; // Made a mistake here but I have to live with it
constexpr num numSQRTSQRTEPS = 1e-2; // Made a mistake here but I have to live with it
 */