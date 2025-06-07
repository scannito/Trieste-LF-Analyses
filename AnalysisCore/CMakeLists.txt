include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Fitter headers & LinkDef
set(FITTER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/LFInvMassFitter.h)
set(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)

# Directory include per rootcling
#set(ROOT_DICTIONARY_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/AnalysisUtils ${CMAKE_CURRENT_SOURCE_DIR}/include)
#set(ROOT_DICT_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/AnalysisUtils ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Generate ROOT dictionary
ROOT_GENERATE_DICTIONARY(G__LFInvMassFitter ${FITTER_HEADERS} LINKDEF ${LINKDEF})

# Add dictionary and sources to the library
add_library(LFInvMassFitter src/LFInvMassFitter.cpp ${CMAKE_CURRENT_BINARY_DIR}/G__LFInvMassFitter.cxx)

# Include dirs for this target
target_include_directories(LFInvMassFitter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/AnalysisUtils ${CMAKE_SOURCE_DIR})

# Link ROOT to the FitterDict library
target_link_libraries(LFInvMassFitter PUBLIC AnalysisUtils ${ROOT_LIBRARIES})

# Copy the .pcm file to the build directory
add_custom_command(TARGET LFInvMassFitter POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${CMAKE_CURRENT_BINARY_DIR}/libLFInvMassFitter_rdict.pcm ${CMAKE_BINARY_DIR}/libLFInvMassFitter_rdict.pcm)