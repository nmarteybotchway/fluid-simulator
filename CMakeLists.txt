cmake_minimum_required(VERSION 3.10)
project(fluid_simulation)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Find SFML package
find_package(SFML 2.5 REQUIRED graphics window system)

# Add your source files
add_executable(fluid_simulation main.cpp
        Fluid.cpp
        Fluid.h
        Grid.h)

# Link SFML libraries to your project
target_link_libraries(fluid_simulation sfml-graphics sfml-window sfml-system)
