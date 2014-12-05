INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(INTERVAL_HEADERS interval.h box.h iboolean.h imatrix.h rmatrix.h)
SET(INTERVAL_SRC interval.cpp box.cpp iboolean.cpp imatrix.cpp rmatrix.cpp)

message("Test, trouve-t'on matrix lib?")
message(" test ${MATRIX_LIB_HEADERS}")
message("alors?")

if (TARGET matrix_lib)
  message("yeah")
  message(" c'est: ${matrix_lib}")
else (TARGET matrix_lib)
  message("booouh")
endif(TARGET matrix_lib)

message ("ZZZ? :(")

foreach( zzz ${MATRIX_LIB_HEADERS})
  message("zzz ${zzz}")
endforeach()


add_library(interval STATIC ${INTERVAL_HEADERS} ${INTERVAL_SRC} ${MATRIX_LIB_HEADERS})
#target_link_libraries(interval)
