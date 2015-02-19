
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(../Utilities) # test LibVideo
include_directories(../NlBayes)   # test VideoNLBayes

message("SOURCE_DIR " ${CMAKE_SOURCE_DIR})

add_executable(test_libvideo testLibVideo.cpp)
add_executable(test_libvideoT testLibVideoT.cpp)
add_executable(test_videonlb testVideoNLBayes.cpp)
add_executable(test_videonlm testVideoNLMeans.cpp)

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(test_libvideo videoNlBayes png)
target_link_libraries(test_libvideoT videoNlBayes png)
target_link_libraries(test_videonlb videoNlBayes png)
target_link_libraries(test_videonlm videoNlBayes png)

