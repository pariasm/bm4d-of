set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_executable(bm4d-of main_vnlb.cpp)
add_executable(awgn awgn.cpp)
add_executable(psnr psnr.cpp)

target_link_libraries(bm4d-of vnldct ${CBLAS_LIBRARIES})
target_link_libraries(awgn vnldct ${CBLAS_LIBRARIES})
target_link_libraries(psnr vnldct ${CBLAS_LIBRARIES})

# These are not necessary
# ${IIO_LIBRARIES})

### # internal includes
### include_directories(${CMAKE_SOURCE_DIR})
### include_directories(${CMAKE_SOURCE_DIR}/utils)
### include_directories(${CMAKE_SOURCE_DIR}/core)

# external excludes
#include_directories(${CBLAS_INCLUDE_DIR})
#include_directories(${IIO_INCLUDE_DIR})

