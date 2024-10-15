cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C CXX)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/src)

# --- Detect C++17 flags -------------------------------------------------------

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag(-funroll-loops HAVE_UNROLL_LOOPS)
check_cxx_compiler_flag(-faligned-new HAVE_ALIGNED_NEW)
if(HAVE_UNROLL_LOOPS)
    add_compile_options(-funroll-loops)
endif()
if(HAVE_ALIGNED_NEW)
    add_compile_options(-faligned-new)
endif()

# --- Detect SIMD --------------------------------------------------------------

include("src/scripts/cmake/FindNEON.cmake")
include("src/scripts/cmake/FindSSE4.cmake")
include("src/scripts/cmake/FindAVX2.cmake")
include("src/scripts/cmake/CythonExtension.cmake")

if(HAVE_SSE4)
    string(STRIP "${SSE4_C_FLAGS}" IMPL_FLAGS)
    if(IMPL_FLAGS)
        add_compile_options(${IMPL_FLAGS})
    endif()
endif()

if(HAVE_AVX2)
    string(STRIP "${AVX2_C_FLAGS}" IMPL_FLAGS)
    if(IMPL_FLAGS)
        add_compile_options(${IMPL_FLAGS})
    endif()
endif()

if(HAVE_NEON)
    string(STRIP "${NEON_C_FLAGS}" IMPL_FLAGS)
    if(IMPL_FLAGS)
        add_compile_options(${IMPL_FLAGS})
    endif()
endif()

# --- Setup include folders ----------------------------------------------------

include_directories(${CMAKE_CURRENT_LIST_DIR}/vendor/tantan/src)
set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/src/pytantan)

# --- Compile source code ------------------------------------------------------

add_subdirectory(src)
