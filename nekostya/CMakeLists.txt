cmake_minimum_required(VERSION 3.14)
project(SPH-STAR)

set(CMAKE_CXX_STANDART 17)

add_executable(${PROJECT_NAME} kotyk.cpp Particle.cpp)

find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)
target_link_libraries(${PROJECT_NAME} sfml-graphics sfml-window sfml-system)

