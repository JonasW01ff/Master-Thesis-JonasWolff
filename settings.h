#pragma once
#include <limits>
#include <cmath>
typedef float num;
constexpr num numEPS = std::numeric_limits<num>::epsilon();
static num numSQRTEPS = sqrt(numEPS); 
static num numSQRTSQRTEPS = sqrt(sqrt(numEPS)); 
