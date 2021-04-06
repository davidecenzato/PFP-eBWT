CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

# Set a default build type if none was specified
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# About this project
# ------------------------------------------------------------------------------
project(BigeBWT)
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "1")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Set environment
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install)


find_package(Git)
if(GIT_FOUND)
    message("git found: ${GIT_EXECUTABLE}")
else()
    message(WARNING "git not found. Cloning of submodules will not work.")
endif()



# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include") # This is an hack because include(GUIInstallDirs) doesn't work

add_subdirectory(thirdparty)


# Configure the compiler with the appropriate flags
# ------------------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  include(ConfigureCompilerClang)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  include(ConfigureCompilerGcc)
else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message(FATAL_ERROR "Only the compiler gcc and clang are supported")
endif()


# Targets
include_directories(${PROJECT_SOURCE_DIR})

add_executable(circpfp.x circpfp.cpp utils.c xerrors.c)
target_link_libraries(circpfp.x malloc_count z pthread)

add_executable(circpfpr.x circpfpmw.cpp utils.c xerrors.c)
target_link_libraries(circpfpr.x malloc_count z pthread)

add_executable(circpfpNT.x circpfp.cpp utils.c)
target_link_libraries(circpfpNT.x malloc_count z)
target_compile_options(circpfpNT.x PUBLIC "-DNOTHREADS")

add_executable(circpfprNT.x circpfpmw.cpp utils.c)
target_link_libraries(circpfprNT.x malloc_count z)
target_compile_options(circpfprNT.x PUBLIC "-DNOTHREADS")

add_executable(bebwtNT.x ebwt.cpp utils.c csais.cpp)
target_link_libraries(bebwtNT.x malloc_count dl z gsacak sdsl  divsufsort divsufsort64)

add_executable(bebwtNT64.x ebwt.cpp utils.c csais.cpp)
target_link_libraries(bebwtNT64.x malloc_count dl z gsacak64 sdsl  divsufsort divsufsort64)
target_compile_options(bebwtNT64.x PUBLIC "-DM64")
target_compile_options(bebwtNT64.x PUBLIC "-DP64")

configure_file(${PROJECT_SOURCE_DIR}/bebwt ${PROJECT_BINARY_DIR}/bebwt)