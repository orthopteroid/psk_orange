cmake_minimum_required(VERSION 3.7)
project(psk_orange)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp PSKOrange.h)
add_executable(psk_orange ${SOURCE_FILES})