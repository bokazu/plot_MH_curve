set(FILES ./src/Int_Hamiltonian.cpp ./src/Iso_Hamiltonian.cpp ./src/Jset.cpp ./src/MEIGEN.cpp ./src/Subsystem_Sz.cpp ./src/Tot_Sz.cpp)

cmake_minimum_required(VERSION 3.20.5)
project(Subsystem_Sz
    LANGUAGES CXX)
add_library(Subsystem_Sz STATIC ${FILES})
target_compile_options(Subsystem_Sz PRIVATE -DMKL_ILP64 -I"${MKLROOT}/include" -fiopenmp -O3 -xHost)
target_include_directories(Subsystem_Sz INTERFACE ./include)