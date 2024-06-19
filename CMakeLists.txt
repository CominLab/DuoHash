cmake_minimum_required(VERSION 3.10)
project(DuoHash)

# Imposta il compilatore g++
set(CMAKE_CXX_COMPILER g++-13)

# Specifica la versione di C++ da utilizzare
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Aggiungi le opzioni di ottimizzazione
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Aggiungi i file header e source
include_directories(${PROJECT_SOURCE_DIR}/include)

# Definisci le directory dei sorgenti
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Crea la libreria statica
add_library(DuoHash STATIC ${SOURCES})

# Aggiungi la libreria agli include
target_include_directories(DuoHash PUBLIC ${PROJECT_SOURCE_DIR}/include)
