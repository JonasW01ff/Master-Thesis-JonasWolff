//
// Created by Jonas Wolff on 14/11/2024.
//

#include "matrixMetal.h"
#include <Foundation/Foundation.h>
#import <Metal/Metal.h>

thread_local id<MTLDevice> device;
thread_local id<MTLLibrary> myMetalLib;
thread_local id<MTLFunction> myMadd;
thread_local id<MTLComputePipelineState> myMaddCPSO;
thread_local NSError *err;
thread_local id<MTLCommandQueue> myComQueue;
thread_local id<MTLCommandBuffer> myComBuffer;
thread_local id<MTLComputeCommandEncoder> myComEncoder;

thread_local id<MTLBuffer> myBufferA;
thread_local id<MTLBuffer> myBufferB;
thread_local id<MTLBuffer> myBufferC;
thread_local size_t buffersize = 10e+6;

matrixMetal::matrixMetal() {

    // Find default device
    device = MTLCreateSystemDefaultDevice();

    //NSURL* myUrl = [NSURL fileURLWithPath:@"shader/products/default"];
    NSURL* myUrl = [NSURL fileURLWithPath:@"default.metallib"];
    myMetalLib = [device newLibraryWithURL:myUrl error:&err];

    // Get matrix adder
    myMadd = [myMetalLib newFunctionWithName:@"madd"];

    // Set adder ComputePipelineStateObject
    myMaddCPSO = [device newComputePipelineStateWithFunction:myMadd error:&err];

    // Set GPU command queue
    myComQueue = [device newCommandQueue];

    // Create buffers
    myBufferA = [device newBufferWithLength:buffersize options:MTLResourceStorageModeShared];
    myBufferB = [device newBufferWithLength:buffersize options:MTLResourceStorageModeShared];
    myBufferC = [device newBufferWithLength:buffersize options:MTLResourceStorageModeShared];

    if (myMetalLib == nil) {
        NSLog(@"Metal Library not found.");
        return;
    }

    // Check if we found madd
    if (myMadd == nil) {
        NSLog(@"Madd not found");
        return;
    }

    isLoaded = true;

}

void matrixMetal::madd(const double *A, const double *B, double *C, size_t n) {

    // Check if we can comply with maxbuffer size
    if (n>buffersize)
        NSLog(@"Input exceeded buffersize.");

    // Set command queue buffer
    myComBuffer = [myComQueue commandBuffer];

    // Set compute encoder for the buffer
    myComEncoder = [myComBuffer computeCommandEncoder];

    // Set encode
    [myComEncoder setComputePipelineState:myMaddCPSO];
    [myComEncoder setBuffer:myBufferA offset:0 atIndex:0];
    [myComEncoder setBuffer:myBufferB offset:0 atIndex:1];
    [myComEncoder setBuffer:myBufferC offset:0 atIndex:2];

    // Populate buffer
    float *bA = static_cast<float *>(myBufferA.contents);
    float *bB = static_cast<float *>(myBufferB.contents);

    vDSP_vdpsp(A,1,bA,1,n);
    vDSP_vdpsp(B,1,bB,1,n);

    //memcpy(myBufferA.contents, A, n);
    //memcpy(myBufferB.contents, B, n);

    // Set grid size
    NSUInteger taskpergrid = n;
    MTLSize gridSize = MTLSizeMake(n,1,1);

    // Set threads per thread group
    NSUInteger tGroupSize = myMaddCPSO.maxTotalThreadsPerThreadgroup;
    if (tGroupSize > gridSize.width)
        tGroupSize = gridSize.width;
    MTLSize threadGroupSize = MTLSizeMake(tGroupSize,1,1);

    // Dispatch
    [myComEncoder dispatchThreads:gridSize threadsPerThreadgroup:threadGroupSize];

    // End encoding
    [myComEncoder endEncoding];

    // Commit
    [myComBuffer commit];

    // Wait
    [myComBuffer waitUntilCompleted];

    // Read results
    float *bC = static_cast<float*>(myBufferC.contents);
    vDSP_vspdp(bC,1,C,1,n);

    //[myComEncoder release];
    //[myComBuffer release];

}
