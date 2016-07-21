
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(../utils) # test LibVideo
include_directories(../core)   # test nldct

#message("SOURCE_DIR " ${CMAKE_SOURCE_DIR})

add_executable(test_lib_videoT test_lib_videot.cpp)
add_executable(patch_group     patch_group.cpp)
add_executable(test_dct        test_dct_handler.cpp)

include_directories(${LAPACK_INCLUDES})

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(test_lib_videoT videoNlBayes ${CBLAS_LIBRARIES})
target_link_libraries(patch_group     videoNlBayes ${CBLAS_LIBRARIES})
target_link_libraries(test_dct        videoNlBayes ${CBLAS_LIBRARIES})

