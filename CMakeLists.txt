cmake_policy(SET CMP0048 NEW) # avoid annoing warnings

project(nldct VERSION 1.0 LANGUAGES C CXX)
message(STATUS "Configuring main project: " ${PROJECT_NAME})

# flags #######################################################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ggdb")

set(CMAKE_BUILD_TYPE Release)

# cmake global  ###############################################################

cmake_minimum_required(VERSION 2.8.2)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")
set(CMAKE_VERBOSE_MAKEFILE OFF)

# dependencies ################################################################

find_package(CBLAS REQUIRED)

find_package( OpenMP REQUIRED)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# add project directories #####################################################

add_subdirectory(lib/iio)
add_subdirectory(src/core)
add_subdirectory(src)

# configure cpack #############################################################

set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})

set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES
"\\\\.*\\\\.swp"
"\\\\*\\\\.o"
"^${PROJECT_SOURCE_DIR}/build/"
"^${PROJECT_SOURCE_DIR}/results/"
"^${PROJECT_SOURCE_DIR}/doc/"
"git")

set(CPACK_SOURCE_PACKAGE_FILE_NAME
	"${PROJECT_NAME}-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}-src")

include(CPack)
# make package <-- TODO : configure properly
# make package_source
