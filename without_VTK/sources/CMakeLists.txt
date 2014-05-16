#set(VTK_DIR "~/Documents/Programming/VTK/VTK-6.1.0/VTK_Build")
set(CMAKE_C_FLAGS "-Wall -g -std=c++0x")
set(CMAKE_BUILD_TYPE Debug)
file( GLOB SOURCES *.h *.cpp)
add_executable(Main ${SOURCES})
add_executable(GridRefinementStudy ${SOURCES})


add_custom_target(clean-all
   COMMAND ${CMAKE_BUILD_TOOL} clean
   COMMAND ${CMAKE_COMMAND} -P clean-all.cmake
)
