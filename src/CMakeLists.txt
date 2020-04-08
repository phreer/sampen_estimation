set(CMAKE_CXX_STANDARD 11)

set(EXECUTABLE_SRC_MAIN sampen.cpp)
set(EXECUTABLE_SRC_VAR sampen_var.cpp)
set(EXECUTABLE_SRC_DYNAMIC_ENTROPY dynamic_entropy.cpp)
set(HEAD_SRC_LIST entropy.h random_sampler.h RangeTree2.h sampen_calculator.h utils.h kdtree.h entropy_calculator.h)
set(LIB_SRC_LIST random_sampler.cpp utils.cpp sampen_calculator.cpp kdtree.cpp entropy_calculator.cpp)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")

message(STATUS ${CMAKE_BUILD_TYPE})

include_directories($ENV{HOME}/local/include)
link_directories($ENV{HOME}/local/lib)

add_library(libsampen SHARED ${LIB_SRC_LIST})
set_target_properties(libsampen PROPERTIES OUTPUT_NAME "sampen")
set_target_properties(libsampen PROPERTIES VERSION 1.0 SUBVERSION 1)
target_link_libraries(libsampen gsl gslcblas)

add_executable(sampen ${EXECUTABLE_SRC_MAIN})
target_link_libraries(sampen libsampen)

add_executable(sampen_var ${EXECUTABLE_SRC_VAR})
target_link_libraries(sampen_var libsampen)

add_executable(dynamic_entropy ${EXECUTABLE_SRC_DYNAMIC_ENTROPY})
target_link_libraries(dynamic_entropy libsampen)