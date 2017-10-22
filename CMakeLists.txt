cmake_minimum_required(VERSION 2.8)

project(VMID CXX)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_FLAGS "-Wall")

include_directories(vmid)
include_directories(vmid/third_party)

add_executable(VMID test/test.cpp)
