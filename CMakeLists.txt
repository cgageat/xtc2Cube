cmake_minimum_required(VERSION 2.8)
project(MyProject)



# check existence of c++11 compiler
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} was n OC++11 support. Please use a different C++ compiler.")
endif()

#Inclusion de Boost

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 -O3 -g -lmkl_intel_thread -lmkl_core -liomp5 -lpthread") #-pg for debug
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Ofast -march=native -fopenmp") #-pg for debug

include(FindBoost)

find_package(

        Boost

        1.36.0

        REQUIRED program_options
        REQUIRED filesystem
        REQUIRED system

)

message( STATUS ${Boost_LIBRARIES} )



#inclusion de Openmp

find_package(OpenMP)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#Configuration de l'exécutable

file(

        GLOB_RECURSE

        source_files

        ${PROJECT_SOURCE_DIR}/src/modifier/*
        ${PROJECT_SOURCE_DIR}/src/constants/*
        ${PROJECT_SOURCE_DIR}/src/functions/*
        ${PROJECT_SOURCE_DIR}/main.cpp
        ${PROJECT_SOURCE_DIR}/src/system/*
        ${PROJECT_SOURCE_DIR}/src/parser/*
        ${PROJECT_SOURCE_DIR}/src/energy/*
        ${PROJECT_SOURCE_DIR}/src/inout/*
        ${PROJECT_SOURCE_DIR}/src/integrator/*
        ${PROJECT_SOURCE_DIR}/src/geometry/*
)

include(ExternalProject)

# xdrfile
ExternalProject_Add(
        xdrfile
        URL ${PROJECT_SOURCE_DIR}/xdrfile-1.1.4
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/xdrfile-cmake
        CONFIGURE_COMMAND <SOURCE_DIR>/configure
        BUILD_COMMAND make -j4
        INSTALL_COMMAND ""
)

SET(XDRFILE_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/xdrfile-cmake/src/xdrfile-build/src/.libs/libxdrfile.a)


#Déclaration de l'exécutable
include_directories(${XDRFILE_INCLUDES})

include_directories(${Boost_INCLUDE_DIRS})

add_executable(

        xtc2Cube

        ${source_files}

)

target_link_libraries(
        xtc2Cube
        ${Boost_LIBRARIES}
)

target_link_libraries(xtc2Cube ${XDRFILE_LIBRARIES})

add_dependencies(xtc2Cube xdrfile)
