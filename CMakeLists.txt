
project(VideoNLBayes)

message(STATUS "Configuring main project: " ${CMAKE_PROJECT_NAME})

set(PROJECT_VERSION_MAJOR  0)
set(PROJECT_VERSION_MINOR  1)
set(PROJECT_VERSION_PATCH  0)

set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
    CACHE INTERNAL "The project version.")

###############################################################################

set(CMAKE_BUILD_TYPE Debug)

###############################################################################

cmake_minimum_required(VERSION 2.8.7)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/")
#set(CMAKE_VERBOSE_MAKEFILE ON)

###############################################################################

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
	message("OPENMP FOUND")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

###############################################################################

add_subdirectory(src/NlBayes)
add_subdirectory(src)
add_subdirectory(src/unit_tests)
