# Fitter headers & LinkDef
set(FITTER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/LFInvMassFitter.h)
set(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)

# Generate ROOT dictionary
ROOT_GENERATE_DICTIONARY(G__LFInvMassFitter ${FITTER_HEADERS} LINKDEF ${LINKDEF})

# Add dictionary and sources to the library
add_library(LFInvMassFitter src/LFInvMassFitter.cpp ${CMAKE_CURRENT_BINARY_DIR}/G__LFInvMassFitter.cxx)

# Include dirs for this target
target_include_directories(LFInvMassFitter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR})

# Link ROOT to the FitterDict library
target_link_libraries(LFInvMassFitter PUBLIC AnalysisUtils ${ROOT_LIBRARIES})