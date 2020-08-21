set(CMAKE_CXX_STANDARD 11)

set(EXECUTABLE_SRC_MAIN sampen.cpp)
set(EXECUTABLE_SRC_VAR sampen_var.cpp)
set(HEAD_LIST "kdtree.h\;random_sampler.h\;RangeTree2.h\;sampen_calculator.h\;tensor.h\;utils.h")
set(LIB_SRC_LIST random_sampler.cpp utils.cpp sampen_calculator.cpp kdtree.cpp)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")

message(STATUS "Headers: " ${HEAD_LIST})
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
include_directories($ENV{HOME}/local/include)
link_directories($ENV{HOME}/local/lib)

add_library(libsampen SHARED ${LIB_SRC_LIST})
set_target_properties(libsampen PROPERTIES OUTPUT_NAME "sampen")
set_target_properties(libsampen PROPERTIES VERSION 1.0 SUBVERSION 1)
set_target_properties(libsampen PROPERTIES PUBLIC_HEADER ${HEAD_LIST})
target_link_libraries(libsampen gsl gslcblas)
install(TARGETS libsampen 
        LIBRARY 
            DESTINATION lib
        PUBLIC_HEADER 
            DESTINATION include
        )

add_executable(sampen ${EXECUTABLE_SRC_MAIN})
target_link_libraries(sampen libsampen)
install(TARGETS sampen
        RUNTIME 
            DESTINATION bin 
        )
add_executable(sampen_var ${EXECUTABLE_SRC_VAR})
target_link_libraries(sampen_var libsampen)

add_executable(test_utils test_utils.cpp)
target_link_libraries(test_utils libsampen)