cmake_minimum_required(VERSION 3.1)

SET(CMAKE_CXX_STANDARD 17)

project(CubeAnalysis VERSION 1.0.0)
SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-pthread")

find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  message("Root has been found")
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

add_subdirectory(CubeRecon/io)
add_subdirectory(CubeRecon/tools)

#add_executable(clark src/Clark.cxx)
#target_link_libraries(clark LINK_PUBLIC cuberecon_io cuberecon_tools ${ROOT_LIBRARIES})

#add_executable(lownu src/Lownu.cxx)
#target_link_libraries(lownu LINK_PUBLIC  cuberecon_io cuberecon_tools ${ROOT_LIBRARIES})

#add_executable(trackTest src/TrackTest.cxx)
#target_link_libraries(trackTest LINK_PUBLIC  cuberecon_io cuberecon_tools ${ROOT_LIBRARIES})

add_executable(EventSelection src/EventSelection.cxx)
target_link_libraries(EventSelection LINK_PUBLIC  cuberecon_io cuberecon_tools ${ROOT_LIBRARIES})
