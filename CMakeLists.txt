cmake_minimum_required(VERSION 3.20.5)
project(Subsystem_Sz LANGUAGES CXX)
add_subdirectory(./Subsystem_Sz)
find_package(IntelDPCPP REQUIRED)
find_package(Boost REQUIRED)
find_package(IntelSYCL REQUIRED)

set(EXECUTABLE_NAME $ENV{PLOT_MH_EXE_FILE})

add_executable(${EXECUTABLE_NAME} 
                main.cpp)

link_libraries("${MKLROOT}")
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${MKL_INCLUDE_DIRS})
target_compile_options(${EXECUTABLE_NAME} PUBLIC -DMKL_ILP64 -I"${MKLROOT}/include" -fiopenmp -O3 -xHost)
target_link_libraries(${EXECUTABLE_NAME}
                        Subsystem_Sz
                        ${MKL_LIBRARIES}
                        mkl_intel_ilp64
                        mkl_intel_thread
                        mkl_core
                        iomp5
                        m
                        pthread
                        iomp5
                        dl
                    )
                        
