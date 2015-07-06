
project(id_dist Fortran)

message(STATUS "Configuring project: " ${CMAKE_PROJECT_NAME})

set(PROJECT_VERSION_MAJOR  0)
set(PROJECT_VERSION_MINOR  4)

set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
    CACHE INTERNAL "The project version.")

###############################################################################

# NOTE: -fPIC is added so that the static library can be included into
#       shared libraries.
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -w -fPIC")
set(CMAKE_BUILD_TYPE Release)

###############################################################################

cmake_minimum_required(VERSION 2.8.2)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")

###############################################################################

#find_package( CBLAS REQUIRED )
find_package( LAPACK REQUIRED )
message(STATUS "LAPACK_LIBRARIES " ${LAPACK_LIBRARIES})

###############################################################################

add_library(id STATIC
            src/prini.f src/dfft.f src/id_rand.f src/id_rtrans.f
            src/idd_snorm.f src/idd_sfft.f src/idd_frm.f src/idd_house.f
            src/idd_qrpiv.f src/idd_id.f src/idd_svd.f src/idd_id2svd.f
            src/iddr_rid.f src/iddr_rsvd.f src/iddp_rid.f src/iddp_rsvd.f
            src/iddr_aid.f src/iddr_asvd.f src/iddp_aid.f src/iddp_asvd.f
            src/idz_snorm.f src/idz_sfft.f src/idz_frm.f src/idz_house.f
            src/idz_qrpiv.f src/idz_id.f src/idz_svd.f src/idz_id2svd.f
            src/idzr_rid.f src/idzr_rsvd.f src/idzp_rid.f src/idzp_rsvd.f
            src/idzr_aid.f src/idzr_asvd.f src/idzp_aid.f src/idzp_asvd.f)

