project(OnzaFDTD)
cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${OnzaFDTD_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})
include(MacroOutOfSourceBuild)
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build.")
add_subdirectory(src)


# If doxygen exists, use the doc/CMakeLists.txt to generate code documentation.
# find_package( Doxygen )
# if( DOXYGEN_FOUND )
#   # This is the directory that contains all the source files we want to have
#   # documented.
#   set( DOXYGEN_SOURCE_DIR ${OnzaFDTD_SOURCE_DIR} )
#   # This is the directory that holds the doxygen doxyfile template (doxy.conf.in)
#   set( DOXYGEN_INPUT_DIR  ${OnzaFDTD_SOURCE_DIR}/doc )
#   # This is the directory where the generated code documentation goes.
#   set( DOXYGEN_OUTPUT_DIR ${OnzaFDTD_BINARY_DIR}/doc )
#   # It is possible to put adddditional CMake commands in sub-directory
#   # CMakeLists.txt, but the add_subdirectory command must be used to specify
#   # which sub-directories should be scanned.
#   add_subdirectory( doc )
# else()
#   message( WARNING "Doxygen not found - code documentation won't be created" )
# endif()


# Directories passed to the include_directories() command will cause all the
# commands that require compilation, e.g. add_executable() or add_library(),
# which follow in the current CMakeLists.txt or sub-directory CMakeLists to have
# the given direction in the compilation command.
#
# For instance, it adds a -I<directory_name> to the gcc commands.
#
# Therefore, add directories with your own header files and header file
# directories of libraries that you use.
# include_directories( ${OnzaFDTD_SOURCE_DIR} )


# # Add any new source files to this list.
# set( OnzaFDTD_SRC
#   src/onza-fdtd.cc
#   )
# # Build the oop executable from the OnzaFDTD_SRC source files.
# add_executable( run-onza ${OnzaFDTD_SRC} onza-fdtd.cc )


# Here is where we specify what gets installed and where.  If we make the
# destinations relative paths, they will be relative to the CMAKE_INSTALL_PREFIX
# configuration variable.
# install( TARGETS OnzaFDTD
#   RUNTIME DESTINATION bin
#   )
