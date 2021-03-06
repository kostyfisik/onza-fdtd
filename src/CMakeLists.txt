# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# macro returns list of subdirs names of directory ${curdir}
macro(SUBDIRLIST result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
        set(dirlist ${dirlist} ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()
# use all subdirs in OnzaFDTD "src" direcory
subdirlist(SUBDIRS ${OnzaFDTD_SOURCE_DIR}/src)
foreach(subdir ${SUBDIRS})
    add_subdirectory(${subdir})
endforeach()

# Create a variable called OnzaFDTD_SOURCES containing all .cc files:
file(GLOB OnzaFDTD_SOURCES *.cc) 
# Create an executable file called run-onzf-fdtd from sources:
add_executable(run-onza-fdtd ${OnzaFDTD_SOURCES})
# subdirs names are synonyms for librarys names
message("Searching for MPI...")
find_package(MPI)
message("Searching for blitz...")

find_path(BLITZ_INCLUDE_DIR blitz/blitz.h)
find_library(BLITZ_LIBRARIES blitz)
set(BLITZ_FOUND 0)
if(BLITZ_INCLUDE_DIR)
  if(BLITZ_LIBRARIES)
    set(BLITZ_FOUND 1)
    message(STATUS "Found Blitz++")
  endif(BLITZ_LIBRARIES)
endif(BLITZ_INCLUDE_DIR)

if (MPI_FOUND AND BLITZ_FOUND)
  include_directories(${MPI_INCLUDE_DIRS} ${BLITZ_INCLUDE_DIR})
  target_link_libraries(run-onza-fdtd ${SUBDIRS} ${MPI_LIBRARIES} ${BLITZ_LIBRRIES} -lblitz)
  install( TARGETS run-onza-fdtd   DESTINATION ./    )
else()
  message( FATAL_ERROR "Onza FDTD needs MPI libs installed!" )
endif()
