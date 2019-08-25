cmake_minimum_required(VERSION 2.8.3)
project(visibility-polygon)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++11)

## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)

find_package(Eigen3 REQUIRED)


# Avoid valgrind error due to overflow error, cf. https://bitbucket.org/ompl/ompl/issues/543
add_compile_definitions(BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS)
add_compile_definitions(BOOST_ALL_NO_LIB)

find_package(Boost REQUIRED)

find_path(decomp_INCLUDE_DIRS convex_decomposition.hpp thirdParty/decomp/source/decomp )
add_subdirectory( thirdParty/decomp )

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
    include
    ${EIGEN3_INCLUDE_DIRS}
    ${BOOST_INCLUDE_DIRS}
    ${decomp_INCLUDE_DIRS}
)

SET(visibility-polygon_LIBRARY_OUTPUT_DIRECTORY ${visibility-polygon_SOURCE_DIR}/lib CACHE PATH "Target for the libraries")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${visibility-polygon_LIBRARY_OUTPUT_DIRECTORY})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${visibility-polygon_LIBRARY_OUTPUT_DIRECTORY})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${visibility-polygon_RUNTIME_OUTPUT_DIRECTORY})

## Declare a C++ library
add_library(${PROJECT_NAME} SHARED
   src/VisibilityPolygon.cpp
    ${decomp_SOURCES}
)


## Specify libraries to link a library or executable target against
target_link_libraries(${PROJECT_NAME}
    ${EIGEN3_LIBRARIES}
    ${BOOST_LIBRARIES}
)

#############
## Install ##
#############

# all install targets should use catkin DESTINATION variables
# See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
# install(PROGRAMS
#   scripts/my_python_script
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark executables and/or libraries for installation
install(TARGETS ${PROJECT_NAME}
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        LIBRARY DESTINATION lib
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

## Mark cpp header files for installation
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
#   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
# )
INSTALL(DIRECTORY ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/include
        FILES_MATCHING PATTERN "*.h" )

## Mark other files for installation (e.g. launch and bag files, etc.)
# install(FILES
#   # myfile1
#   # myfile2
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )




#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_ob1.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
