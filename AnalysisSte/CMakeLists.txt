# Copy config.json files from AnalysisSte to build directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/configK0S.json ${CMAKE_BINARY_DIR}/configK0S.json COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/configPi.json ${CMAKE_BINARY_DIR}/configPi.json COPYONLY)

set(ANALYSIS_CORE "${CMAKE_SOURCE_DIR}/AnalysisCore")
set(RAPIDJSON "${CMAKE_SOURCE_DIR}/Thirdparty/rapidjson/include")

set(SOURCES_COMMON ${ANALYSIS_CORE}/src/JSONReader.cpp ${ANALYSIS_CORE}/src/THnSparseProjector.cpp ${ANALYSIS_CORE}/src/LFInvMassFitter.cpp)

foreach(i RANGE 1 3)
    set(name PhiK0SAnalysis${i})
    add_executable(${name} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cpp ${SOURCES_COMMON})
    target_include_directories(${name} PRIVATE ${ANALYSIS_CORE} ${RAPIDJSON})
    target_link_libraries(${name} PRIVATE LFInvMassFitter ${ROOT_LIBRARIES})
endforeach()