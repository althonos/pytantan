file(COPY ../vendor/tantan DESTINATION .)
file(COPY ./tantan DESTINATION .)
file(CONFIGURE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/tantan/src/version.hh CONTENT "\"${CMAKE_PROJECT_VERSION}\"")

set(TANTAN_CC ${CMAKE_CURRENT_BINARY_DIR}/tantan/src/tantan.cc)

add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/tantan ${CMAKE_CURRENT_BINARY_DIR}/tantan/build)
add_subdirectory(pytantan)
