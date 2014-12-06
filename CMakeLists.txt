cmake_minimum_required(VERSION 2.8)

project(LIB_INTERVAL)

include_directories("." ${CMAKE_CURRENT_LIST_DIR}/3rd/matrix_lib)
add_subdirectory(src)
