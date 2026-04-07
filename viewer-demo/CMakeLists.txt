set(VIEWER_NAME "c4hex_viewer_demo")

add_executable(${VIEWER_NAME} main.cpp)
target_link_libraries(${VIEWER_NAME} MC3D::MC3D)
set_target_properties(${VIEWER_NAME} PROPERTIES
                        CXX_STANDARD 17
                        CXX_STANDARD_REQUIRED ON
                        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
