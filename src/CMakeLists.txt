
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Utilities)
include_directories(${CMAKE_SOURCE_DIR}/NlBayes)

add_executable(vnlbayes main_vnlb.cpp)

#message("-- LAPACK libs: " ${LAPACK_LIBRARIES})
#message("-- LAPACK incs: " ${LAPACK_INCLUDES})
include_directories(${LAPACK_INCLUDES})
include_directories(${CBLAS_INCLUDE_DIR})

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(vnlbayes videoNlBayes png
	${LAPACK_LIBRARIES} ${CBLAS_LIBRARIES} ${IDDIST_LIBRARIES})

