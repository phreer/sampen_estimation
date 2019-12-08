set(CMAKE_CXX_STANDARD 11)

set(EXECUTABLE_SRC_LIST sampen.cpp)
set(HEAD_SRC_LIST entropy.h random_sampler.h RangeTree2.h sampen_calculator.h utils.h)
set(LIB_SRC_LIST random_sampler.cpp utils.cpp sampen_calculator.cpp)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3")

add_compile_options(-g -Wall)
include_directories($ENV{HOME}/local/include)
link_directories($ENV{HOME}/local/lib)

add_library(libsampen SHARED ${LIB_SRC_LIST})
set_target_properties(libsampen PROPERTIES OUTPUT_NAME "sampen")
set_target_properties(libsampen PROPERTIES VERSION 1.0 SUBVERSION 1)

add_executable(sampen ${EXECUTABLE_SRC_LIST})
target_link_libraries(sampen libsampen)
target_link_libraries(libsampen gsl gslcblas)