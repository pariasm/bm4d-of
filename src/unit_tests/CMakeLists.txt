
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(../Utilities) # test LibVideo
include_directories(../NlBayes)   # test VideoNLBayes

#message("SOURCE_DIR " ${CMAKE_SOURCE_DIR})

add_executable(test_libvideoT testLibVideoT.cpp)
add_executable(test_libmatrix testLibMatrix.cpp)
add_executable(test_videonlb  testVideoNLBayes.cpp)
add_executable(patch_group    patch_group.cpp)

include_directories(${LAPACK_INCLUDES})

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(test_libvideoT videoNlBayes png ${LAPACK_LIBRARIES})
target_link_libraries(test_libmatrix videoNlBayes png ${LAPACK_LIBRARIES})
target_link_libraries(test_videonlb  videoNlBayes png ${LAPACK_LIBRARIES})
target_link_libraries(patch_group    videoNlBayes png ${LAPACK_LIBRARIES})

