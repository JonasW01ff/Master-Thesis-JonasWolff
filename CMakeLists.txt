cmake_minimum_required(VERSION 3.26)
project(Code)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g")
set(CMAKE_CXX_STANDARD 23)

find_package(Matplot++ REQUIRED)

# Pass the current PATH to the runtime environment
set(ENVIRONMENT
        PATH=${CURRENT_PATH}
)

add_executable(Code settings.h
        distributionFunctions.h
        sobol.h
        main.cpp
        date.h
        date.cpp
        bond.h
        bond.cpp
        rates.h
        rates.cpp
        matrix.h
        matrix.cpp
        interpolation.h
        interpolation.cpp
        optimizer.h
        optimizer.cpp
        numbers.h
        numbers.cpp
        termstructureStore.h
        termstructureHermite.cpp
        termstructureHermite.h
        termstructureForwardHermite.h
        blockList.h
        concurrentQueue.h
        ThreadPool.cpp
        ThreadPool.h
        adjointNode.h
        adjointTape.h
        adjointTape.cpp
        adjointNumeric.h
        #matrixMetal.mm
        #matrixMetal.h
        shortrate.h
        shortrateHWEV.h
        quadrature.h
        plot.h
        ois.h
        ois.cpp
        shortrateBachelier.h
)

set(METAL_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/shaders/products/release/default.metallib")

add_custom_command(
        TARGET Code POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${METAL_LIBRARY_PATH} $<TARGET_FILE_DIR:Code>
        COMMENT "Copying matrix.metallib to the build directory"
)

set(ENVIRONMENT
        PATH=${GNUPLOT_DIR}:$ENV{PATH}
)

# Link the Accelerate framework
target_link_libraries(Code "-framework Accelerate" "-framework Metal" "-framework MetalPerformanceShaders" "-framework Foundation" Matplot++::matplot)
