cmake_minimum_required (VERSION 3.4)

project (cursor_demo VERSION 0.1.0)

find_package (SeqAn3 REQUIRED HINTS lib/seqan3/build_system)

add_executable ("${PROJECT_NAME}" src/main.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC seqan3::seqan3)

