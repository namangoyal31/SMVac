cmake_minimum_required(VERSION 3.10)
project(SMVac VERSION 2.0.0 LANGUAGES CXX
        DESCRIPTION "Standard Model electroweak vacuum stability: analytical approximations vs numerical reference")

# C++14: the reference toolchain (MinGW.org GCC 6.3) lacks C++17 library support.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
  # Use IEEE double arithmetic (SSE2) on 32-bit x86 so results do not depend on x87 extended precision.
  if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    add_compile_options(-msse2 -mfpmath=sse)
  endif()
endif()

add_library(smvac STATIC
  src/matching.cpp
  src/beta_functions.cpp
  src/rge.cpp
  src/potential.cpp
  src/analytic.cpp
  src/bounce.cpp
  src/decay.cpp
  src/stability.cpp)
target_include_directories(smvac PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(smvac_point apps/smvac_point.cpp)
target_link_libraries(smvac_point PRIVATE smvac)
add_executable(smvac_scan apps/smvac_scan.cpp)
target_link_libraries(smvac_scan PRIVATE smvac)

include(CTest)
if(BUILD_TESTING)
  foreach(t test_beta_functions test_rge test_potential test_analytic test_bounce test_decay test_pipeline)
    add_executable(${t} tests/${t}.cpp)
    target_link_libraries(${t} PRIVATE smvac)
    add_test(NAME ${t} COMMAND ${t})
  endforeach()
endif()

# Validation program (not built by default): cmake --build build --target precision_check
add_executable(precision_check EXCLUDE_FROM_ALL validation/precision_check.cpp)
target_link_libraries(precision_check PRIVATE smvac)
