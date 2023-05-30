cmake_minimum_required(VERSION 3.7)
project(psk_orange)

set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES main.cpp PSKOrange.h)

add_executable(testencode ${SOURCE_FILES})
target_compile_definitions(testencode PUBLIC TARGET_TESTENCODE)

add_executable(dumpdecode ${SOURCE_FILES})
target_compile_definitions(dumpdecode PUBLIC TARGET_DUMPDECODE)

add_executable(debugdecode ${SOURCE_FILES})
target_compile_definitions(debugdecode PUBLIC TARGET_DEBUGDECODE)

add_executable(encode ${SOURCE_FILES})
target_compile_definitions(encode PUBLIC TARGET_ENCODE)

add_executable(decode ${SOURCE_FILES})
target_compile_definitions(decode PUBLIC TARGET_DECODE)
