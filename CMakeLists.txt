cmake_minimum_required(VERSION 3.20.5)
project(Subsystem_Sz LANGUAGES CXX)
add_subdirectory(./Subsystem_Sz)
find_package(IntelDPCPP REQUIRED)
find_package(Boost REQUIRED)
find_package(IntelSYCL REQUIRED)

add_executable(main_app 
                main.cpp)

link_libraries("${MKLROOT}")
target_include_directories(main_app PRIVATE ${MKL_INCLUDE_DIRS})
target_compile_options(main_app PUBLIC -DMKL_ILP64 -I"${MKLROOT}/include" -fiopenmp -O3 -xHost)
target_link_libraries(main_app
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
                        
