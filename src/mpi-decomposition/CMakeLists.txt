# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
file(GLOB current_dir_src *.cc) 
file(GLOB sim_dir ../simulation-core/*.cc) 
file(GLOB profiling_dir ../profiling/*.cc) 
get_filename_component(lib_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
include_directories(../simulation-core)
include_directories(../profiling)
add_library(${lib_name} ${current_dir_src} ${sim_dir} ${profiling_dir})