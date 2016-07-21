
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/src/utils)
include_directories(${IIO_INCLUDE_DIR})
include_directories(${FFTW_INCLUDES})

add_library(videoNlBayes SHARED
	nldct.cpp
	matrix_funs.cpp
	dct_handler.cpp
	../utils/lib_image.cpp
	../utils/LibVideoT.cpp
	../utils/mt19937ar.c
	../utils/utilities.cpp
	)

target_link_libraries(videoNlBayes iio ${FFTWF_LIBRARIES})

set_target_properties(videoNlBayes PROPERTIES
	SOVERSION ${PROJECT_VERSION})

install(TARGETS videoNlBayes
    LIBRARY DESTINATION lib
    COMPONENT libs)
