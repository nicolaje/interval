INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(INTERVAL_HEADERS interval.h box.h iboolean.h imatrix.h rmatrix.h)
SET(INTERVAL_SRC interval.cpp box.cpp iboolean.cpp imatrix.cpp rmatrix.cpp)

add_library(interval STATIC ${INTERVAL_HEADERS} ${INTERVAL_SRC})
