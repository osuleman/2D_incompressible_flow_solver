#set(VTK_DIR "~/Documents/Programming/VTK/VTK-6.1.0/VTK_Build")
set(CMAKE_C_FLAGS "-Wall -g -std=c++11")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
file( GLOB SOURCES *.h *.cpp)
add_executable(Main ${SOURCES})
target_link_libraries(Main ${VTK_LIBRARIES})


