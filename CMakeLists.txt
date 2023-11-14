cmake_minimum_required(VERSION 3.10)

project(C4Hex VERSION 1.0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(macrOS)
include(GNUInstallDirs)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build this in a source directory (or any directory with a CMakeLists.txt file). "
                        "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

### Build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
  "Debug" "Release" "RelWithDebInfo")
endif()

### Standalone?
if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(C4Hex_STANDALONE TRUE)
    message("Configuring C4Hex as standalone project...")
else()
    set(C4Hex_STANDALONE FALSE)
    message("Configuring C4Hex inside another cmake project...")
endif()

### Options
option(C4Hex_ENABLE_LOGGING    "Enable logging for C4Hex" ${C4Hex_STANDALONE})
option(C4Hex_BUILD_CLI         "Build CLI app for C4Hex"  ${C4Hex_STANDALONE})
option(C4Hex_SUBMODULES_MANUAL "Skip automatic submodule download" OFF)
option(BUILD_SHARED_LIBS       "Build libraries as shared as opposed to static" ON)

set(MC3D_ENABLE_LOGGING ${C4Hex_ENABLE_LOGGING})
set(QGP3D_ENABLE_LOGGING ${C4Hex_ENABLE_LOGGING})

### Global cmake settings
set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_EXTENSIONS OFF)
if (C4Hex_STANDALONE)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/lib/)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/lib/)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/bin/)
endif()

### Load external libs
if(NOT ${C4Hex_SUBMODULES_MANUAL})
    add_subdirectory(extern)
endif()

### Add library
# lib target: C4Hex::C4Hex
add_subdirectory(src)

### Add CLI
if (C4Hex_BUILD_CLI)
    add_subdirectory(cli)
endif()

### Fake successful finder run if compiling as a dependent project.
if (NOT C4Hex_STANDALONE)
    set(C4Hex_FOUND true PARENT_SCOPE)
    set(C4Hex_LIBRARIES C4Hex PARENT_SCOPE)
    set(C4Hex_LIBRARY C4Hex PARENT_SCOPE)
    set(C4Hex_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include/" PARENT_SCOPE)
    set(C4Hex_LIBRARY_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" CACHE PATH "The directory where the C4Hex libraries can be found.")
endif()
