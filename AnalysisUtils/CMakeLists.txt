# Create interface library for headers
add_library(AnalysisUtils INTERFACE)

target_include_directories(AnalysisUtils INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})