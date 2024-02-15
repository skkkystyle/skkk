cmake_minimum_required(VERSION 3.26)
project(LongNumskkk LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(LongNumskkk SHARED
	include/LongNumskkk.hpp
	src/LongNumskkk.cpp
)

add_executable(pi pi.cpp)
target_link_libraries(pi LongNumskkk)

add_executable(test testers.cpp)
target_link_libraries(test LongNumskkk)

add_executable(cat cat.cpp)

add_custom_target(runpi COMMAND pi)
add_custom_target(runtest COMMAND test)
add_custom_target(runcat COMMAND cat)