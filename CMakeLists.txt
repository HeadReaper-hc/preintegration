cmake_minimum_required(VERSION 3.14)
project(preintegration)
set(CMAKE_CXX_STANDARD 14)

#eigen path
include_directories("/usr/include/eigen3")

#find_package Sophus
find_package(Sophus REQUIRED)
include_directories(${Sophus_INCLUDE_DIRS})

#add dir<include> to the header path
include_directories("include")

#find_package Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

#find_package Ceres
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

file(GLOB CPPFILE src/*.cpp)

add_executable(main main.cpp ${CPPFILE})
target_link_libraries(main ${Sophus_LIBRARIES} ${Boost_LIBRARIES} ${CERES_LIBRARIES})