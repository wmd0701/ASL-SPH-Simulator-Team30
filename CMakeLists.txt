cmake_minimum_required(VERSION 3.0.0)

project(SPH_sloshing)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(run)

target_sources(run PRIVATE main.c)


