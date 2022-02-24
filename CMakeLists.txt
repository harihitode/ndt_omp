cmake_minimum_required(VERSION 3.5)
project(ndt_ocl)

add_definitions(-std=c++14)
#set(CMAKE_CXX_FLAGS "-std=c++14 -MM -E")
#add_compile_options(-march=native)

# pcl 1.7 causes a segfault when it is built with debug mode
set(CMAKE_BUILD_TYPE "RELEASE")

if($ENV{ROS_VERSION} EQUAL 1)
# ROS1
find_package(catkin REQUIRED COMPONENTS
  roscpp
  pcl_ros
)

find_package(PCL 1.7 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

message(STATUS "PCL_INCLUDE_DIRS:" ${PCL_INCLUDE_DIRS})
message(STATUS "PCL_LIBRARY_DIRS:" ${PCL_LIBRARY_DIRS})
message(STATUS "PCL_DEFINITIONS:" ${PCL_DEFINITIONS})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

###################################
## catkin specific configuration ##
###################################
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ndt_ocl
)

###########
## Build ##
###########
include_directories(include)
include_directories(
  ${catkin_INCLUDE_DIRS}
)

add_library(ndt_ocl
  src/pclocl/voxel_grid_covariance_ocl.cpp
  src/pclocl/ndt_ocl.cpp
)

#############
## INSTALL ##
#############

install(
  TARGETS
    ndt_ocl
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
)

# install headers
install(DIRECTORY include/pclocl
        DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION})
else()
# ROS2
find_package(ament_cmake_auto REQUIRED)
find_package(OpenCL REQUIRED)
ament_auto_find_build_dependencies()

find_package(PCL 1.7 REQUIRED)
# include_directories(${PCL_INCLUDE_DIRS})
# link_directories(${PCL_LIBRARY_DIRS})
# add_definitions(${PCL_DEFINITIONS})

message(STATUS "PCL_INCLUDE_DIRS:" ${PCL_INCLUDE_DIRS})
message(STATUS "PCL_LIBRARY_DIRS:" ${PCL_LIBRARY_DIRS})
message(STATUS "PCL_DEFINITIONS:" ${PCL_DEFINITIONS})

find_package(ament_index_cpp REQUIRED)
find_package(OpenMP)

###########
## Build ##
###########

ament_auto_add_library(ndt_ocl SHARED
  src/pclocl/voxel_grid_covariance_ocl.cpp
  src/pclocl/ndt_ocl.cpp
  include/pclocl/voxel_grid_covariance_ocl.h
  include/pclocl/voxel_grid_covariance_ocl_impl.hpp
  include/pclocl/ndt_ocl.h
  include/pclocl/ndt_ocl_impl.hpp
)

target_link_libraries(ndt_ocl ${PCL_LIBRARIES})

if(OpenMP_CXX_FOUND)
  target_link_libraries(ndt_ocl OpenMP::OpenMP_CXX)
else()
  message(WARNING "OpenMP not found")
endif()

ament_auto_package()

install(DIRECTORY kernel DESTINATION share/${PROJECT_NAME})
endif()
