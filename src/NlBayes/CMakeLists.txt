
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/src/Utilities)

add_library(videoNlBayes SHARED
	NlBayes.cpp
	VideoNLBayes.cpp
	LibMatrix.cpp
	../Utilities/io_png.c
	../Utilities/LibImages.cpp
	../Utilities/LibVideo.cpp
	../Utilities/LibVideoT.hpp
	../Utilities/LibVideoUtils.cpp
	../Utilities/mt19937ar.c
	../Utilities/Utilities.cpp
	)

set_target_properties(videoNlBayes PROPERTIES
	SOVERSION ${PROJECT_VERSION})

install(TARGETS videoNlBayes
    LIBRARY DESTINATION lib
    COMPONENT libs)
