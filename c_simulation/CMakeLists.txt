cmake_minimum_required(VERSION 3.16)
project(MyProject CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable debug symbols and optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fdiagnostics-color=always")

# Include directory
include_directories(${CMAKE_SOURCE_DIR}/inc)

# Source files
add_executable(sim
    inc/common.cpp
    sim.cpp
)

# Destination directory
set_target_properties(sim PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)