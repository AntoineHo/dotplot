cmake_minimum_required (VERSION 2.8.12)
PROJECT(dotplotter)

## Sets version number
set (MAJOR 0)
set (MINOR 1)

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
FILE(GLOB SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)

ADD_EXECUTABLE(dotplotter ${SOURCES})

TARGET_LINK_LIBRARIES(dotplotter ${Boost_LIBRARIES})
