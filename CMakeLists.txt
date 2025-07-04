# Trieste-LF-Analyses/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(Trieste-LF-Analyses CXX)

# Set standard C++
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RooFit RooFitCore)
include(${ROOT_USE_FILE})

# Add subdirectories
add_subdirectory(AnalysisUtils)
add_subdirectory(AnalysisCore)
add_subdirectory(AnalysisSte)

# Fetch rapidjson
#include(FetchContent)
#FetchContent_Declare(rapidjson GIT_REPOSITORY https://github.com/Tencent/rapidjson.git GIT_TAG master)
#FetchContent_MakeAvailable(rapidjson)

#target_link_libraries(PhiK0SAnalysis PRIVATE ${rapidjson_SOURCE_DIR}/include)

# Generate compile_commands.json
#set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Opzione configurabile da riga di comando
#set(ANALYSIS_DIR "" CACHE STRING "Name of analysis directory to include and copy")

# Check che sia stata definita
#if(NOT ANALYSIS_DIR)
    #message(FATAL_ERROR "You must define -DANALYSIS_DIR=AnalysisSte or similar")
#endif()

# Path assoluto
#set(ANALYSIS_PATH "${CMAKE_SOURCE_DIR}/${ANALYSIS_DIR}")

# Verifica che esista davvero
#if(NOT EXISTS "${ANALYSIS_PATH}")
    #message(FATAL_ERROR "Directory ${ANALYSIS_PATH} does not exist!")
#endif()

# Crea l'eseguibile con il file specificato
#file(GLOB ANALYSIS_SOURCES "${ANALYSIS_PATH}/*.cpp")
#add_executable(${ANALYSIS_DIR}_exe ${ANALYSIS_SOURCES})

# Includi le directory necessarie
#target_include_directories(${ANALYSIS_DIR}_exe PRIVATE
    #${CMAKE_SOURCE_DIR}/AnalysisUtils
    #${CMAKE_SOURCE_DIR}/Thirdparty/rapidjson/include
#)

# Linka librerie
#target_link_libraries(${ANALYSIS_DIR}_exe PRIVATE AnalysisCore AnalysisUtils ${ROOT_LIBRARIES})

# Copia la directory selezionata nel build dir
#add_custom_command(TARGET ${ANALYSIS_DIR}_exe POST_BUILD
    #COMMAND ${CMAKE_COMMAND} -E copy_directory
    #"${ANALYSIS_PATH}" "${CMAKE_BINARY_DIR}/${ANALYSIS_DIR}"
    #COMMENT "Copying selected analysis directory ${ANALYSIS_DIR} to build directory"
#)