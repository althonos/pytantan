set(TANTAN_SOURCES
    "cbrc_linalg.cc"
    "LambdaCalculator.cc"
    "mcf_alphabet.cc"
    "mcf_tantan_options.cc"
)

add_library(libtantan ${TANTAN_SOURCES})
target_include_directories(libtantan PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(libtantan PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
