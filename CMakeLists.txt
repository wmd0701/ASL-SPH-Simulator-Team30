cmake_minimum_required(VERSION 3.10.0)

project(SPH_sloshing)

add_executable(run)

target_include_directories(run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(run PRIVATE main.c)


