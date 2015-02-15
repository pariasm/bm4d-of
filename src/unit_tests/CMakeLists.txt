
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(../Utilities) # test LibVideo
include_directories(../NlBayes)   # test VideoNLBayes

message("SOURCE_DIR " ${CMAKE_SOURCE_DIR})

add_executable(testlibvideo testLibVideo.cpp)
add_executable(testlibvideot testLibVideoT.cpp)
add_executable(testvideonlb testVideoNLBayes.cpp)

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(testlibvideo videoNlBayes png)
target_link_libraries(testlibvideot videoNlBayes png)
target_link_libraries(testvideonlb videoNlBayes png)

