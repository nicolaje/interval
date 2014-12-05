cmake_minimum_required(VERSION 2.8)

project(LIB_INTERVAL)
set(MATRIX_LIB_HEADERS 0)

add_subdirectory(3rd)

set(MATRIX_LIB_HEADERS PARENT_SCOPE)
set(INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE)

add_subdirectory(src)


include_directories("." ${INCLUDE_DIRS})
