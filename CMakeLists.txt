# Set the minimum version of cmake (and default installation rpath if needed)
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
if (APPLE)
    set(CMAKE_MACOSX_RPATH TRUE)
endif ()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# One may want to impose the C-compiler (for example with a specific version)
#set(CMAKE_C_COMPILER "/usr/bin/clang-3.8")
#set(CMAKE_C_COMPILER_ID Clang)

# Print a warning if the os is not unix (linux and mac to be more precise)
if (NOT UNIX)
    message("Warning: the VelocityExtension library is currently intended for ")
    message("Unix systems since it has been tested merely on Linux and ")
    message("Mac OS. The commands 'cmake ..' and 'make' may not work.")
endif ()

# Define the project of the VelocityExtension library depending on cmake version
if (CMAKE_VERSION VERSION_GREATER 3.0)
    project(VelocityExtension LANGUAGES C)
else ()
    project(VelocityExtension C)
endif ()
set(CMAKE_BUILD_TYPE "Release")

# Print a warning if the compiler is not gcc or clang (apple clang on mac)
if (NOT ${CMAKE_C_COMPILER_ID} MATCHES GNU
                 AND NOT ${CMAKE_C_COMPILER_ID} MATCHES Clang
                              AND NOT ${CMAKE_C_COMPILER_ID} MATCHES AppleClang)
    message("Warning: the VelocityExtension library is currently intended to ")
    message("be compiled by gcc or clang since it has been tested merely on ")
    message("these two compilers. The commands 'cmake ..' and 'make' may not ")
    message("work.")
endif ()

# Add the source files to the Mshdistance library which also uses
#    * Headers from the C standard library:
#         - assert.h
#         - ctype.h
#         - float.h
#         - limits.h
#         - math.h (the associated static library needs to be linked)
#         - signal.h
#         - stdio.h
#         - stdlib.h
#         - string.h
#
#    * Headers from the external (dynamic) Commons library, which needs to be
#      properly linked:
#         - chrono.h
#         - libmesh5.h
#         - sparse.h
#
#    * Headers from the local source directory:
#        - velext.h
#        - vl_calls.h
add_library(VelocityExtension SHARED sources/boulep.c
                                     sources/hash.c
                                     sources/inout.c
                                     sources/packing.c
                                     sources/velex1_2d.c
                                     sources/velex1_3d.c
                                     #sources/velext.c
                                     sources/vl_calls.c)

# Find the Commons library (all the other libraries are its dependancies)
find_library(COMMONS_LIBRARY NAMES Commons HINTS "$ENV{HOME}/lib/iscd")

# Find the Math library (libm.a on Ubuntu 16.04.4, add -lm with gcc compiler)
find_library(MATH_LIBRARY m)

# Set all the targets and properties of the VelocityExtension library
target_include_directories(VelocityExtension
    PUBLIC
        $<BUILD_INTERFACE:$ENV{HOME}/include/iscd/Commons>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources>)
        # Remark: the user can add directly the path(s) for header file(s) in
        #         the above target_include_directories command (such as the ones
        #         associated with the Commons library)

if (CMAKE_VERSION VERSION_GREATER 3.8)
    target_compile_features(VelocityExtension PRIVATE c_std_99)
elseif (CMAKE_VERSION VERSION_GREATER 3.1)
    set_property(TARGET VelocityExtension PROPERTY C_STANDARD 99)
else ()
    target_compile_options(VelocityExtension PRIVATE -std=gnu99)
endif ()

#if (${CMAKE_C_COMPILER_ID} MATCHES GNU)
    #target_compile_options(VelocityExtension PRIVATE -Wno-char-subscripts)
    #target_compile_options(VelocityExtension PRIVATE -Wpedantic -Wall -Wextra)
#elseif (${CMAKE_C_COMPILER_ID} MATCHES Clang
#                                  OR ${CMAKE_C_COMPILER_ID} MATCHES AppleClang)
    #if (APPLE)
        # Important remark: on apple clang, the 'info' global variable of
        #                   velext.c causes problems when linking the
        #                   VelocityExtension library because it is not defined
        #                   inside the library but only in mshdist(.c) software
        #target_link_libraries(VelocityExtension PRIVATE -Wl,-U,_info)
    #endif ()
#else ()
    #target_compile_options(VelocityExtension PRIVATE -g)
    #target_compile_options(VelocityExtension PRIVATE -Weverything)
#endif ()

# Remark: the user can add directly here the path for linking the math and
#         and Commons libraries (if they have not been found with the above
#         find_library commands)
target_link_libraries(VelocityExtension PUBLIC ${MATH_LIBRARY}
                                               ${COMMONS_LIBRARY})

# Set the installation properties of the VelocityExtension library
# Remark: the user can change here the paths where the VelocityExtension library
#         and associated header files will be installed at the 'cmake install'
#         command
install(TARGETS VelocityExtension LIBRARY DESTINATION "$ENV{HOME}/lib/iscd")
install(FILES sources/velext.h
              sources/vl_calls.h
                        DESTINATION "$ENV{HOME}/include/iscd/VelocityExtension")

# Create and install the corresponding velext executable
if (CMAKE_VERSION VERSION_GREATER 3.0)
    project(velext LANGUAGES C)
else ()
    project(velext C)
endif ()
set(CMAKE_BUILD_TYPE "Release")

add_executable(velext sources/velext.c)

if (CMAKE_VERSION VERSION_GREATER 3.8)
    target_compile_features(velext PRIVATE c_std_99)
elseif (CMAKE_VERSION VERSION_GREATER 3.1)
    set_property(TARGET velext PROPERTY C_STANDARD 99)
else ()
    target_compile_options(velext PRIVATE -std=gnu99)
endif ()

if (${CMAKE_C_COMPILER_ID} MATCHES GNU)
    target_compile_options(velext PRIVATE -Wno-unused-result)
    #target_compile_options(velext PRIVATE -Wpedantic -Wall -Wextra)
#else ()
    #target_compile_options(velext PRIVATE -g)
    #target_compile_options(velext PRIVATE -Weverything)
endif ()

target_link_libraries(velext PUBLIC VelocityExtension)
install(TARGETS velext RUNTIME DESTINATION "$ENV{HOME}/bin")

