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

# Source/include path
set(ANALYSIS_UTILS "${CMAKE_SOURCE_DIR}/AnalysisUtils")
set(ANALYSIS_STE "${CMAKE_SOURCE_DIR}/AnalysisSte")

# Define the executable
add_executable(PhiK0SAnalysis ${ANALYSIS_STE}/PhiK0SAnalysis.cpp)

# Link both your library and ROOT to the executable
target_link_libraries(PhiK0SAnalysis PRIVATE AnalysisCore AnalysisUtils ${ROOT_LIBRARIES})