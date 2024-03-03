cmake_minimum_required(VERSION 3.1)

PROJECT(dustAgglomerates)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

option (DO_PLOT "Plot particles?" ON)

set(CMAKE_VERBOSE_MAKEFILE ON) # verbose on



if(NOT DO_PLOT)
find_package(CUDA REQUIRED)
endif()

if(DO_PLOT)


if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/usr/lib/x86_64-linux-gnu/cmake/Qt5Widgets")
else()
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "D:/Programme/QT/5.9.3/msvc2015_64/lib/cmake/Qt5Widgets")
endif()

find_package(OpenGL REQUIRED)

# Find the QtWidgets library
find_package( Qt5Widgets REQUIRED )
find_package( Qt5Gui REQUIRED )
find_package( Qt5OpenGL REQUIRED )


set( CMAKE_AUTOMOC ON )  # to automoc remove the need to wrap by hand



if(CMAKE_SYSTEM_NAME STREQUAL "Linux")

    FIND_PATH( GLEW_INCLUDE_DIR GL/glew.h
    /usr/include
    /usr/local/include
    /sw/include
    /opt/local/include
    ${GLEW_ROOT_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
    DOC "The directory where GL/glew.h resides")
    set(GLEW_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/GL)


    FIND_LIBRARY( GLEW_LIBRARY
    NAMES GLEW glew
    PATHS
    /usr/lib64
    /usr/lib
    /usr/local/lib64
    /usr/local/lib
    /sw/lib
    /opt/local/lib
    ${GLEW_ROOT_DIR}/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/linLib)
    #DOC "The GLEW library")
    SET(GLEW_LIBRARIES ${GLEW_LIBRARY})
    find_package_handle_standard_args(GLEW DEFAULT_MSG GLEW_LIBRARY GLEW_INCLUDE_PATH)
    find_package(GLEW REQUIRED)

  # otherwise some bullet internal headers don't find friends..
  include_directories(${CUDA_INCLUDE}) # ${CMAKE_CURRENT_SOURCE_DIR}/linlib/include
else()
  #windows
  link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/winLib)
endif()

endif() # DO_PLOT

if(DO_PLOT)
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
else()
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/nvidiaHelperFuncs ${CMAKE_CURRENT_SOURCE_DIR}/cub)
endif()


## main file to build

if(DO_PLOT)
    set(MAIN "Main/Main.cpp")
else()
    #set(MAIN "Main/main_test_timestep.cpp")
    #set(MAIN "Main/main_determine_cor.cpp")
    #set(MAIN "Main/main_determine_v_s_r.cpp")
    #set(MAIN "Main/main_fit_parameters.cpp")
    #set(MAIN "Main/main_cuda_compactor.cpp")
    #set(MAIN "Main/main_cuda_compression_sw.cpp")
    #set(MAIN "Main/main_cuda_compression_sw_fixed_size.cpp")
    #set(MAIN "Main/main_cuda_copy_test.cpp")
    set(MAIN "Main/main_cuda_agg_collider.cpp")
    #set(MAIN "Main/main_cuda_agg_collider_generator.cpp")
    #set(MAIN "Main/main_coll_analyzer.cpp")
    #set(MAIN "Main/main_cuda_sphere_compaction.cpp")
    #set(MAIN "Main/main_cuda_load.cpp")
endif()

message(STATUS "PLOT?  " ${DO_PLOT} "  " ${MAIN})

if(NOT DO_PLOT)
FILE(GLOB CU_SOURCES "*.cu")
FILE(GLOB CUH_SOURCES "*.cuh")
endif()


if(DO_PLOT)
FILE(GLOB CPP_SOURCES "*.cpp")
FILE(GLOB H_SOURCES "*.h")
else()
FILE(GLOB CPP_SOURCES "*.cpp"  "cuda/*.cpp")
FILE(GLOB H_SOURCES "*.h" "nvidiaHelperFuncs/*.h"  "cuda/*.h")
endif()

set(CPP_SOURCES ${MAIN} ${CPP_SOURCES})
include_directories(${CMAKE_CURRENT_BINARY_DIR}/Main)

if(DO_PLOT)

    list(REMOVE_ITEM H_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/SimulationCuda.h)
    list(REMOVE_ITEM CPP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/SimulationCuda.cpp)
    list(REMOVE_ITEM H_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CudaWallInteraction.h)
    list(REMOVE_ITEM CPP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CudaWallInteraction.cpp)
    ## include Widget files
    FILE(GLOB CPP_VISUALIZATION_SOURCES "visualization/*.cpp" "visualization/subdir/*.cpp")
    FILE(GLOB H_VISUALIZATION_SOURCES "visualization/*.h" "visualization/subdir/*.h")
    set(CPP_SOURCES ${CPP_SOURCES} ${CPP_VISUALIZATION_SOURCES})
    set(H_SOURCES ${H_SOURCES} ${H_VISUALIZATION_SOURCES})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/visualization)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/visualization/subdir)

    set(QT_LIBS Qt5::Widgets Qt5::OpenGL Qt5::Gui)

else() # !DO_PLOT

    list(REMOVE_ITEM H_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/nvidiaHelperFuncs/helper_gl.h)

endif() # DO_PLOT

if(NOT DO_PLOT)
# Pass options to NVCC
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-O3 -gencode arch=compute_52,code=sm_52")
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_50,code=sm_50" "-x cu -O3 --ptxas-options=-v" "-lineinfo")
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_37,code=sm_37" "-x cu -O3 --ptxas-options=-v" "-lineinfo")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_50,code=sm_50" "-x cu -g --ptxas-options=-v" "-lineinfo" "-fmad=false")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   add_compile_options(-std=c++11)
   set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -std=c++11")
else()
    SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}" )
endif()


if(DO_PLOT)
ADD_EXECUTABLE(${PROJECT_NAME} ${CPP_SOURCES} ${H_SOURCES})
else()
#SET(CUDA_PROPAGATE_HOST_FLAGS OFF)
CUDA_ADD_EXECUTABLE(${PROJECT_NAME} ${CU_SOURCES} ${CPP_SOURCES} ${CUH_SOURCES} ${H_SOURCES})
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)



if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
	if(DO_PLOT)
		set(GLEW_LIBS GL glut GLU pthread)
	endif()
    set(LIBS ${LIBS} ${GLEW_LIBRARIES} ${GLEW_LIBS} ${CUDA_LIBRARIES}) #X11 Xxf86vm Xi GL glfw3 GLEW Xrandr pthread assimp BulletDynamics BulletCollision LinearMath fmodex64 freeimage gsl gslcblas
        target_link_libraries(${PROJECT_NAME} ${LIBS} ${QT_LIBS})
else()
	if(DO_PLOT)
		set(GLEW_LIBS glew64 opengl32 GLU32)
	endif()
	
    set(LIBS ${LIBS} ${GLEW_LIBS} ${CUDA_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${LIBS} ${QT_LIBS})
endif()
