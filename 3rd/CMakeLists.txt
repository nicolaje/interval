
set(MATRIX_LIB_HEADERS PARENT_SCOPE)
add_subdirectory(matrix_lib)
foreach( ccc ${MATRIX_LIB_HEADERS})
  message("ccc ${ccc}")
endforeach()
set(MATRIX_LIB_HEADERS ${MATRIX_LIB_HEADERS} PARENT_SCOPE)
set(INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE)
