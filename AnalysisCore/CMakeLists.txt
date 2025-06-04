# Fitter headers & LinkDef
set(FITTER_HEADERS include/LFInvMassFitter.h)
set(LINKDEF LinkDef.h)

# Generate ROOT dictionary
ROOT_GENERATE_DICTIONARY(G__LFInvMassFitter ${FITTER_HEADERS} LINKDEF ${LINKDEF})

# Add dictionary and sources to the library
add_library(AnalysisCore src/LFInvMassFitter.cpp ${CMAKE_BINARY_DIR}/G__LFInvMassFitter.cxx)

# Link ROOT to the FitterDict library
target_link_libraries(AnalysisCore PUBLIC ${ROOT_LIBRARIES})