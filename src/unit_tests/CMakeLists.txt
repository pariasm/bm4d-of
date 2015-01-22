
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(../Utilities)

message("SOURCE_DIR " ${CMAKE_SOURCE_DIR})

add_executable(testlibvideo testLibVideo.cpp)

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(testlibvideo videoNlBayes png)

