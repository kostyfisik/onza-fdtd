# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Create a variable called OnzaFDTD_SOURCES containing all .cc files:
add_subdirectory(mpi-decomposition)
file(GLOB OnzaFDTD_SOURCES *.cc) 
#file(GLOB OnzaFDTD_H_FILES *.h) 
# Create an executable file called run-onzf-fdtd from sources:
#link_directories(${OnzaFDTD_SOURCES}/mpi-decomposition)
# add_executable(run-onza-fdtd ${OnzaFDTD_SOURCES} ${OnzaFDTD_H_FILES})
add_executable(run-onza-fdtd ${OnzaFDTD_SOURCES})
target_link_libraries(run-onza-fdtd halo-exchange-process)
install( TARGETS run-onza-fdtd
   DESTINATION ./
  )