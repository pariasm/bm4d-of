
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Utilities)
include_directories(${CMAKE_SOURCE_DIR}/NlBayes)

message("SOURCE_DIR " ${CMAKE_SOURCE_DIR})

add_executable(nlbayes main.cpp)

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(nlbayes videoNlBayes png)

