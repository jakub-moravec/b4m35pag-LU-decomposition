cmake_minimum_required(VERSION 3.8)
project(lu_hw1)

#Bring the headers, such as Student.h into the project
include_directories(include)

set(CMAKE_CXX_STANDARD 11)

#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "*.cpp")

add_executable(lu ${SOURCES})