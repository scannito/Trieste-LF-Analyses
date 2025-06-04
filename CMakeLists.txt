# Trieste-LF-Analyses/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(Trieste-LF-Analyses CXX)

# Set standard C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RooFit RooFitCore)
include(${ROOT_USE_FILE})

# Add subdirectories
add_subdirectory(AnalysisUtils)
add_subdirectory(AnalysisCore)
add_subdirectory(AnalysisSte)

# Source/include path
set(INCLUDE_ROOT "${CMAKE_SOURCE_DIR}/AnalysisUtils")
set(SRC_ROOT "${CMAKE_SOURCE_DIR}/AnalysisCore/src")
set(HEADER_ROOT "${CMAKE_SOURCE_DIR}/AnalysisCore/include")
set(ANALYSIS_STE "${CMAKE_SOURCE_DIR}/AnalysisSte")

# Fitter headers & LinkDef
set(FITTER_HEADERS ${HEADER_ROOT}/Fitter.h)
set(LINKDEF ${HEADER_ROOT}/LinkDef.h)

# Generate ROOT dictionary
ROOT_GENERATE_DICTIONARY(G__Fitter ${FITTER_HEADERS} LINKDEF ${LINKDEF})

# Add dictionary and sources to the library
add_library(FitterDict ${SRC_ROOT}/Fitter.cpp ${CMAKE_BINARY_DIR}/G__Fitter.cxx)

# Add the include directories to the project
include_directories(${INCLUDE_ROOT} ${HEADER_ROOT} ${CMAKE_BINARY_DIR})

# Link ROOT to the FitterDict library
target_link_libraries(FitterDict ${ROOT_LIBRARIES})

# Define the executable
add_executable(PhiK0SAnalysis ${ANALYSIS_STE}/PhiK0SAnalysis.cpp)

# Link both your library and ROOT to the executable
target_link_libraries(PhiK0SAnalysis FitterDict ${ROOT_LIBRARIES})

# Target di installazione (opzionale)
#install(TARGETS mergeTree DESTINATION bin)