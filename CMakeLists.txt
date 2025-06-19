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
#set(ANALYSIS_UTILS "${CMAKE_SOURCE_DIR}/AnalysisUtils")
set(ANALYSIS_CORE "${CMAKE_SOURCE_DIR}/AnalysisCore")
set(ANALYSIS_STE "${CMAKE_SOURCE_DIR}/AnalysisSte")
set(RAPIDJSON "${CMAKE_SOURCE_DIR}/Thirdparty/rapidjson/include")

# Copy config.json files from AnalysisSte to build directory
configure_file(${ANALYSIS_STE}/configK0S.json ${CMAKE_BINARY_DIR}/configK0S.json COPYONLY)
configure_file(${ANALYSIS_STE}/configPi.json ${CMAKE_BINARY_DIR}/configPi.json COPYONLY)

# Define the executable
add_executable(PhiK0SAnalysis1 ${ANALYSIS_STE}/PhiK0SAnalysis1.cpp ${ANALYSIS_CORE}/src/ObjectHandler.cpp)
add_executable(PhiK0SAnalysis2 ${ANALYSIS_STE}/PhiK0SAnalysis2.cpp ${ANALYSIS_CORE}/src/ObjectHandler.cpp ${ANALYSIS_CORE}/src/LFInvMassFitter.cpp)
add_executable(PhiK0SAnalysis3 ${ANALYSIS_STE}/PhiK0SAnalysis3.cpp ${ANALYSIS_CORE}/src/ObjectHandler.cpp ${ANALYSIS_CORE}/src/LFInvMassFitter.cpp)

# Include directories
target_include_directories(PhiK0SAnalysis1 PRIVATE ${ANALYSIS_CORE} ${RAPIDJSON})
target_include_directories(PhiK0SAnalysis2 PRIVATE ${ANALYSIS_CORE} ${RAPIDJSON})
target_include_directories(PhiK0SAnalysis3 PRIVATE ${ANALYSIS_CORE} ${RAPIDJSON})

# Link both your library and ROOT to the executable
target_link_libraries(PhiK0SAnalysis1 PRIVATE LFInvMassFitter ${ROOT_LIBRARIES})
target_link_libraries(PhiK0SAnalysis2 PRIVATE LFInvMassFitter ${ROOT_LIBRARIES})
target_link_libraries(PhiK0SAnalysis3 PRIVATE LFInvMassFitter ${ROOT_LIBRARIES})

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