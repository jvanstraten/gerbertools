cmake_minimum_required(VERSION 3.16.0)
project(gerber CXX)

add_executable(
    gerber
    ${CMAKE_CURRENT_SOURCE_DIR}/src/clipper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/coord.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/path.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/plot.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/aperture.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/aperture_macro.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gerber.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ncdrill.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pcb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
)
target_include_directories(
    gerber
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)
