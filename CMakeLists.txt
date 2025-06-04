# Trieste-LF-Analyses/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(Trieste-LF-Analyses CXX)

set(INCLUDE_ROOT "AnalysisUtils")
set(SRC_ROOT "LFInvMassFitter")
set(ANALYSIS_STE "AnalysisSte")

# Add the include directories to the project
include_directories(${INCLUDE_ROOT})
#include_directories(${CMAKE_SOURCE_DIR})

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RooFit RooFitCore)
include(${ROOT_USE_FILE}) # Importa le variabili e le funzioni di ROOT

add_executable(PhiK0SAnalysis ${ANALYSIS_STE}/PhiK0SAnalysis.cpp ${SRC_ROOT}/Fitter.cpp ${SRC_ROOT}/SignalExtraction.cpp)

# Collega le librerie ROOT all'eseguibile
target_link_libraries(PhiK0SAnalysis ${ROOT_LIBRARIES})

# Target di installazione (opzionale)
#install(TARGETS mergeTree DESTINATION bin)