cmake_minimum_required(VERSION 3.17)
project(ME-Illumination)

set(CMAKE_CXX_STANDARD 17)

# my attempt at compiler optimization
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fopenmp -march=skylake-avx512")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
# end attempt

IF (WIN32)
    #nothing
ELSE()
    link_libraries(stdc++fs)
ENDIF()

add_executable(ME-Illumination main.cpp main.h illumination.cpp illumination.h lab.cpp lab.h test.cpp cvt.cpp cvt.h)