
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/src/utils)
include_directories(${IIO_INCLUDE_DIR})
include_directories(${FFTW_INCLUDES})

add_library(vnldct SHARED
	nldct.cpp
	matrix_funs.cpp
	../utils/lib_image.cpp
	../utils/lib_videot.cpp
	../utils/mt19937ar.c
	../utils/utilities.cpp
	)

target_link_libraries(vnldct iio ${FFTWF_LIBRARIES})

set_target_properties(vnldct PROPERTIES
	SOVERSION ${PROJECT_VERSION})

install(TARGETS vnldct
    LIBRARY DESTINATION lib
    COMPONENT libs)
