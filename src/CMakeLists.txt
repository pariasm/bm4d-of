set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# internal includes
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/utils)
include_directories(${CMAKE_SOURCE_DIR}/core)

# external excludes
#include_directories(${CBLAS_INCLUDE_DIR})
include_directories(${IIO_INCLUDE_DIR})

add_executable(nldct main_vnlb.cpp)

#link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(nldct vnldct ${CBLAS_LIBRARIES})

configure_file(scripts/nldct-mp4.sh ${EXECUTABLE_OUTPUT_PATH}/nldct-mp4.sh COPYONLY)

# These are not necessary
# ${IIO_LIBRARIES})

