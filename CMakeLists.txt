cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE  "Release")
INCLUDE_DIRECTORIES(  "$ENV{HOME}/include")
LINK_DIRECTORIES(     "$ENV{HOME}/lib")
SET(CMAKE_CXX_FLAGS   "-Wunused-result -Wno-dev")

#Velocity extension library
project(Velexlib)
file(                  GLOB_RECURSE source_files sources/*)
file(                  GLOB_RECURSE header_files sources/*.h)
ADD_LIBRARY(           Velexlib SHARED ${source_files})
find_library(          Commons NAMES Commons HINTS "$ENV{HOME}/lib")
target_link_libraries( Velexlib Commons)
INSTALL(               TARGETS Velexlib LIBRARY  DESTINATION "$ENV{HOME}/lib")
INSTALL(               FILES ${header_files} DESTINATION "$ENV{HOME}/include")

#Executable
project(main)
add_executable(        velext sources/velext.c)
target_link_libraries( velext Commons Velexlib)
INSTALL(               TARGETS velext RUNTIME DESTINATION "$ENV{HOME}/bin")


 



