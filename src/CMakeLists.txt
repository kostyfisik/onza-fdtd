# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Create a variable called helloworld_SOURCES containing all .cpp files:
file(GLOB OnzaFDTD_SOURCES *.cc) 
# Create an executable file called run-onzf-fdtd from sources:
add_executable(run-onza-fdtd ${OnzaFDTD_SOURCES})
install( TARGETS run-onza-fdtd
   DESTINATION ./
  )