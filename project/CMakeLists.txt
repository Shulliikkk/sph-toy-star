cmake_minimum_required(VERSION 3.14)
project(SPH-STAR)

set(CMAKE_CXX_STANDART 17)
set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_BUILD_TYPE Debug)

set(SOURCES
	src/kotyk.cpp
	src/Particle.cpp
	src/Timer.cpp
	src/Processing.cpp
	src/Visualisation.cpp
	src/Menu.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)
target_link_libraries(${PROJECT_NAME} sfml-graphics sfml-window sfml-system)

