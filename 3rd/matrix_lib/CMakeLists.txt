cmake_minimum_required(VERSION 3.0)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(MATRIX_LIB_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/jama_cholesky.h ${CMAKE_CURRENT_SOURCE_DIR}/jama_eig.h ${CMAKE_CURRENT_SOURCE_DIR}/jama_lu.h ${CMAKE_CURRENT_SOURCE_DIR}/jama_qr.h
${CMAKE_CURRENT_SOURCE_DIR}/jama_svd.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_array1d.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_array1d_utils.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_array2d.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_array2d_utils.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_array3d.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_array3d_utils.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_cmat.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_fortran_array1d.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_fortran_array1d_utils.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_fortran_array2d.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_fortran_array2d_utils.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_fortran_array3d.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_fortran_array3d_utils.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_i_refvec.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_math_utils.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_sparse_matrix_csr.h
${CMAKE_CURRENT_SOURCE_DIR}/tnt_stopwatch.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_subscript.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_vec.h ${CMAKE_CURRENT_SOURCE_DIR}/tnt_version.h PARENT_SCOPE)

set(INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE)

ADD_LIBRARY(matrix_lib INTERFACE)# ${MATRIX_LIB_HEADERS})

#ADD_CUSTOM_TARGET(matrix_lib SOURCES ${MATRIX_LIB_HEADERS})

INSTALL(TARGETS matrix_lib DESTINATION ${CMAKE_SOURCE_DIR}/includes/matrix_lib EXPORT matrix_lib_target)
INSTALL(FILES ${MATRIX_LIB_HEADERS} DESTINATION ${CMAKE_SOURCE_DIR}/includes/matrix_lib)
EXPORT(TARGETS matrix_lib FILE ${CMAKE_BINARY_DIR}/lib/matrix_lib/matrix_lib-targets.cmake)
#CONFIGURE_FILE(matrix_lib-targets.cmake.in lib/matrix_lib/matrix_lib-targets.cmake @ONLY)
