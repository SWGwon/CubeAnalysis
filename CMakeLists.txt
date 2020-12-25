cmake_minimum_required(VERSION 3.1)

SET(CMAKE_CXX_STANDARD 17)

project(CubeAnalysis VERSION 1.0.0)

find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  message("Root has been found")
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

set(source
    src/main.cxx
    src/Event.cxx
    src/Primaries.cxx
    src/Object.cxx
    )

add_subdirectory(CubeRecon/io)
add_subdirectory(CubeRecon/tools)

add_executable(CubeAnalysis ${source})
target_link_libraries(CubeAnalysis LINK_PUBLIC cuberecon_io cuberecon_tools ${ROOT_LIBRARIES})

add_executable(clark src/Clark.cxx)
target_link_libraries(clark LINK_PUBLIC cuberecon_io cuberecon_tools ${ROOT_LIBRARIES})
