# TriesteLFAnalyses/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(TriesteLFAnalyses CXX)

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

# Enable compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Define macro if it's an open source project
add_compile_definitions(PVS_STUDIO_FREE)

# PVS-Studio integration
find_program(PVS_STUDIO_ANALYZER pvs-studio-analyzer)
find_program(PLOG_CONVERTER plog-converter)

# Symlink automatico di compile_commands.json nella root (se non già esistente)
if(CMAKE_EXPORT_COMPILE_COMMANDS)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_BINARY_DIR}/compile_commands.json
      ${CMAKE_SOURCE_DIR}/compile_commands.json
    RESULT_VARIABLE _symlink_result
  )
endif()

if(PVS_STUDIO_ANALYZER AND PLOG_CONVERTER)
  add_custom_target(pvs-analysis
    COMMENT "Running PVS-Studio static analysis..."
    COMMAND ${PVS_STUDIO_ANALYZER} analyze
            -f ${CMAKE_BINARY_DIR}/compile_commands.json
            -o ${CMAKE_BINARY_DIR}/PVS-Studio.log
    COMMAND ${PLOG_CONVERTER}
            -a GA:1,2
            -t html
            -o ${CMAKE_BINARY_DIR}/PVS-Studio.html
            ${CMAKE_BINARY_DIR}/PVS-Studio.log
    COMMAND ${PLOG_CONVERTER}
            -a GA:1,2
            -t sarif
            -o ${CMAKE_BINARY_DIR}/PVS-Studio.sarif
            ${CMAKE_BINARY_DIR}/PVS-Studio.log
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
else()
  message(WARNING "PVS-Studio tools not found: pvs-analysis target will not be available.")
endif()


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