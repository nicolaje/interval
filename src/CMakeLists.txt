file(GLOB interval_SRC "*.h" "*.cpp")

add_library(interval STATIC ${matrix_lib_SRC} ${interval_SRC} )
