### Find dependencies
# qgp3d
if (NOT TARGET QGP3D::QGP3D)
    if(EXISTS "${PROJECT_SOURCE_DIR}/extern/QGP3D/CMakeLists.txt")
        add_subdirectory(${PROJECT_SOURCE_DIR}/extern/QGP3D extern/QGP3D EXCLUDE_FROM_ALL)
    else()
        find_package(QGP3D REQUIRED)
    endif()
endif()
list(APPEND C4Hex_LIB_LIST QGP3D::QGP3D)

# Clp
find_package(CLP REQUIRED)
list(APPEND C4Hex_LIB_LIST_PRV ${CLP_LIBRARIES})

# nlopt
find_package(NLOPT)
if (NLOPT_FOUND)
    list(APPEND C4Hex_COMPILE_DEFINITIONS_PRV "C4HEX_WITH_TLC")
    list(APPEND C4Hex_LIB_LIST_PRV "${NLOPT_LIBRARIES}")
endif()

#eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
list(APPEND C4Hex_LIB_LIST Eigen3::Eigen)

#lbfgspp
if (NOT TARGET lbfgspp)
    if(EXISTS "${PROJECT_SOURCE_DIR}/extern/LBFGSpp/CMakeLists.txt")
        add_subdirectory(${PROJECT_SOURCE_DIR}/extern/LBFGSpp extern/LBFGSpp EXCLUDE_FROM_ALL)
    else()
        find_package(lbfgspp REQUIRED)
    endif()
endif()
list(APPEND C4Hex_LIB_LIST lbfgspp)

### Add all source files
list(APPEND C4Hex_SOURCE_LIST
    "Algorithm/EmbeddingCollapser.cpp"
    "Algorithm/HexExtractor.cpp"
    "Algorithm/IGMInitializer.cpp"
    "Algorithm/IGMUntangler.cpp"
    "Algorithm/MCCollapser.cpp"
    "Algorithm/MCSmoother.cpp"
    "Algorithm/MCSplitter.cpp"
    "Algorithm/PathRouter.cpp"
    "Algorithm/SurfaceRouter.cpp"
    "Interface/HexRemesher.cpp"
    "Interface/IGMGenerator.cpp")

### Create target
add_library(C4Hex ${C4Hex_SOURCE_LIST})
add_library(C4Hex::C4Hex ALIAS C4Hex)

### Configure install paths
include(CMakePackageConfigHelpers)
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR "CMake")
else()
  set(DEF_INSTALL_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/C4Hex")
endif()
set(C4Hex_INSTALL_CMAKE_DIR "${DEF_INSTALL_CMAKE_DIR}" CACHE PATH "Installation directory for CMake files")
set(C4Hex_INSTALL_LIB_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "Installation directory for libraries")
set(C4Hex_INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Installation directory for header files")
# ... for the build tree
if (NLOPT_FOUND)
  set(CONF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../extern/tlc" "${CMAKE_CURRENT_SOURCE_DIR}/../include")
else()
  set(CONF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../include")
endif()
set(CONF_LIB_DIRS "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
configure_package_config_file(../C4HexConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/C4HexConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
# ... for the install tree
set(CONF_INCLUDE_DIRS "${C4Hex_INSTALL_INCLUDE_DIR}")
set(CONF_LIB_DIRS "${C4Hex_INSTALL_LIB_DIR}")
configure_package_config_file(../C4HexConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/C4HexConfig.cmake"
  INSTALL_DESTINATION "${C4Hex_INSTALL_CMAKE_DIR}")

### Include own headers with public access
target_include_directories(C4Hex
                           PUBLIC
                           "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
                           "$<INSTALL_INTERFACE:${C4Hex_INSTALL_INCLUDE_DIR}>")

target_include_directories(C4Hex SYSTEM
                           PRIVATE
                           "$<BUILD_INTERFACE:${CLP_INCLUDE_DIRS}>")

if (NLOPT_FOUND)
  target_include_directories(C4Hex
                             PUBLIC
                             "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/tlc/>"
                             "$<INSTALL_INTERFACE:${C4Hex_INSTALL_INCLUDE_DIR}>")
endif()


### Link with dependencies
target_link_libraries_system(C4Hex PUBLIC ${C4Hex_LIB_LIST})
target_link_libraries_system(C4Hex PRIVATE ${C4Hex_LIB_LIST_PRV})

### Set compile options
# properties
set_target_properties(C4Hex PROPERTIES
                      DEBUG_POSTFIX "d"
                      VERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}
                      SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED ON)
# compile options
list(APPEND C4Hex_COMPILE_OPTIONS_PRV "-Wall" "-Wextra" "-Wshadow" "-Wunused" "-Wunused-function")
target_compile_options(C4Hex PRIVATE ${C4Hex_COMPILE_OPTIONS_PRV})
# preprocessor defines
target_compile_definitions(C4Hex PRIVATE ${C4Hex_COMPILE_DEFINITIONS_PRV})

### Install and export
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/C4Hex"
        DESTINATION "${C4Hex_INSTALL_INCLUDE_DIR}")

install(FILES
        "${PROJECT_SOURCE_DIR}/extern/tlc/TLC.h"
        DESTINATION "${C4Hex_INSTALL_INCLUDE_DIR}")

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/C4HexConfig.cmake"
        DESTINATION "${C4Hex_INSTALL_CMAKE_DIR}")

install(TARGETS       C4Hex
        EXPORT        C4HexTargets
        LIBRARY       DESTINATION "${C4Hex_INSTALL_LIB_DIR}"
        ARCHIVE       DESTINATION "${C4Hex_INSTALL_LIB_DIR}")

install(EXPORT C4HexTargets
        FILE C4HexTargets.cmake
        NAMESPACE C4Hex::
        DESTINATION "${C4Hex_INSTALL_CMAKE_DIR}")

export(EXPORT C4HexTargets
       NAMESPACE C4Hex::)
