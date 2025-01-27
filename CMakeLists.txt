#cmake3 -DOPENCL_LIBRARIES=/usr/lib64/libOpenCL.so -DOPENCL_INCLUDE_DIR=/usr/local/cuda/include/ .. 
cmake_minimum_required(VERSION 3.5)
project(OCLCHIP)

set(CMAKE_CXX_STANDARD 11)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -D_FORCE_INLINES")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wpedantic -Wall")

#set(GPU_ARCHITECTURE "sm_61" CACHE STRING "Target GPU architecture")

#set(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda-9.2")
#set(CUDA_LIBRARY_DIR /usr/local/cuda-9.2/lib64/)

#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${CUDA_NVCC_DEBUG_FLAGS} -arch=${GPU_ARCHITECTURE} -use_fast_math)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#find_package(CUDA 9.2 REQUIRED)
#find_package(Thrust)
set(BOOST_ROOT "/home/qianjiaqiang/boost_1_66_0")
set(BOOST_NO_SYSTEM_PATHS FALSE)
find_package(OpenCL REQUIRED)
find_package(HTSlib)
find_package(Threads)
find_package(ZLIB REQUIRED)
find_package(Boost REQUIRED COMPONENTS serialization)

set(CHIP_SOURCE_DIR ${CMAKE_SOURCE_DIR})
set(GFF_SOURCE_DIR ${CMAKE_SOURCE_DIR}/third_party/gff)
#set(COMPUTE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/third_party/compute)
file(GLOB GFF_SRCS ${GFF_SOURCE_DIR}/src/*.cpp)
list(APPEND Boost_LIBRARIES "/home/qianjiaqiang/boost_1_66_0/stage/lib/libboost_serialization.so")
include_directories(
	#        ${CUDA_INCLUDE_DIRS}
	#        ${THRUST_INCLUDE_DIR}
	${OPENCL_INCLUDE_DIR}
        ${HTSlib_INCLUDE_DIR}
        ${CHIP_SOURCE_DIR}/include
        ${ZLIB_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIR}
        ${GFF_SOURCE_DIR}/include
	#${COMPUTE_SOURCE_DIR}/include
)

#cuda_add_executable(paean ${CHIP_SOURCE_DIR}/src/gene.cu ${GFF_SRCS})
add_executable(paean ${CHIP_SOURCE_DIR}/src/gene.cpp ${GFF_SRCS} )
add_executable(tg ${CHIP_SOURCE_DIR}/src/tg.cpp ${GFF_SRCS} )
target_link_libraries(paean ${HTSlib_LIBRARY} ${ZLIB_LIBRARIES} ${OPENCL_LIBRARIES} ${Boost_LIBRARIES} Threads::Threads)
target_link_libraries(tg ${HTSlib_LIBRARY} ${ZLIB_LIBRARIES} ${OPENCL_LIBRARIES} ${Boost_LIBRARIES} Threads::Threads)
