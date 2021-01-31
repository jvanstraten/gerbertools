cmake_minimum_required(VERSION 3.4...3.18)
project(gerbertools CXX)
include(GNUInstallDirs)
set (CMAKE_CXX_STANDARD 11)

# Target with all the normal C++ files.
add_library(gerbertools_objlib OBJECT
    ${CMAKE_CURRENT_SOURCE_DIR}/src/clipper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/coord.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/path.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/plot.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/aperture.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/aperture_macro.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gerber.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ncdrill.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pcb.cpp
)
set_property(
    TARGET gerbertools_objlib
    PROPERTY POSITION_INDEPENDENT_CODE ON
)
# target_compile_options(
#     gerbertools_objlib
#     PUBLIC -Wall -Wextra -pedantic -Werror
# )
target_include_directories(
    gerbertools_objlib
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Library.
add_library(
    gerbertools SHARED
)
target_link_libraries(
    gerbertools
    PUBLIC gerbertools_objlib
)
install(
    TARGETS gerbertools
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install (
    DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.hpp*"
)

# Command line interface.
add_executable(
    gerbertools_cli
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
)
set_property(
    TARGET gerbertools_cli
    PROPERTY OUTPUT_NAME gerbertools
)
target_link_libraries(
    gerbertools_cli
    PUBLIC gerbertools_objlib
)
install(
    TARGETS gerbertools_cli
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

# Python extension.
add_subdirectory(python/pybind11)
pybind11_add_module(
    _gerbertools
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pymod.cpp
)
target_link_libraries(
    _gerbertools
    PUBLIC gerbertools_objlib
)
if (PYTHON_EXT)
    install(
        FILES "$<TARGET_FILE:_gerbertools>"
        DESTINATION "${PYTHON_DIR}"
        RENAME "${PYTHON_EXT}"
    )
endif()
