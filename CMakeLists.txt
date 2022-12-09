
cmake_minimum_required(VERSION 2.8)
 
PROJECT(WriteSTL)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(ReadSTL MACOSX_BUNDLE writeSTL.cpp )
 
target_link_libraries(ReadSTL ${VTK_LIBRARIES})
