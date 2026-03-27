cmake_minimum_required(VERSION 3.18)

project(MolPG
    VERSION 1.0
    DESCRIPTION "Get point group of a molecule."
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# if(NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type" FORCE)
# endif()

find_package(Eigen3 3.4 REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 found: Eigen3::Eigen")
else()
    message(FATAL_ERROR "Eigen3 not found. Please set EIGEN3_INCLUDE_DIR or CPLUS_INCLUDE_PATH")
endif()

find_package(fmt 8.0 REQUIRED CONFIG)
if(fmt_FOUND)
    message(STATUS "fmt found: fmt::fmt")
else()
    message(FATAL_ERROR "fmt not found. Please set fmt_DIR or CPLUS_INCLUDE_PATH and LIBRARY_PATH")
endif()

file(GLOB SOURCES
    src/*.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME} PRIVATE 
    Eigen3::Eigen
    fmt::fmt
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "")
message(STATUS "=== Configuration Summary ===")
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Eigen3: Eigen3::Eigen")
message(STATUS "fmt: fmt::fmt")
message(STATUS "Extra optimizations: ${ENABLE_EXTRA_OPTIMIZATIONS}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "=================================")
