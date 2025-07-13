# Cerca tutti i file JSON che iniziano per config nel source dir (AnalysisSte/JSONConfigs)
file(GLOB CONFIG_JSON_FILES "${CMAKE_CURRENT_SOURCE_DIR}/JSONConfigs/config*.json")

# Copy all found JSON files to build/AnalysisSte directory
set(CONFIG_DST_DIR "${CMAKE_BINARY_DIR}/AnalysisSte")
file(MAKE_DIRECTORY ${CONFIG_DST_DIR})

# Cicla su ciascun file trovato
foreach(json_file ${CONFIG_JSON_FILES})
    get_filename_component(json_name ${json_file} NAME)
    configure_file(${json_file} ${CONFIG_DST_DIR}/${json_name} COPYONLY)
endforeach()

set(ANALYSIS_CORE "${CMAKE_SOURCE_DIR}/AnalysisCore")
set(RAPIDJSON "${CMAKE_SOURCE_DIR}/Thirdparty/rapidjson/include")

set(SOURCES_COMMON ${ANALYSIS_CORE}/src/JSONReader.cpp ${ANALYSIS_CORE}/src/THnSparseProjector.cpp 
                   ${ANALYSIS_CORE}/src/LFInvMassFitter.cpp ${ANALYSIS_CORE}/src/EfficiencyHandler.cpp)

foreach(i RANGE 1 2)
    set(name PhiK0SAnalysis${i})
    add_executable(${name} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cpp ${SOURCES_COMMON})
    target_include_directories(${name} PRIVATE ${ANALYSIS_CORE} ${RAPIDJSON})
    target_link_libraries(${name} PRIVATE LFInvMassFitter ${ROOT_LIBRARIES})
endforeach()