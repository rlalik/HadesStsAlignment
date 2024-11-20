cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    HadesStsAlignment
    VERSION 0.1.0
    DESCRIPTION "Sts Detector Alignment"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES CXX
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
)

include(FetchContent)

FetchContent_Declare(ProMille
  GIT_REPOSITORY https://github.com/rlalik/promille
  GIT_TAG renaming
)
FetchContent_MakeAvailable(ProMille)

FetchContent_Declare(
  hades-tools
  GIT_REPOSITORY https://git.gsi.de/r.lalik/hades-tools
  GIT_TAG master
)
FetchContent_MakeAvailable(hades-tools)
list(APPEND CMAKE_MODULE_PATH ${hades-tools_SOURCE_DIR}/cmake)

find_package(HYDRA)
find_package(ROOT)

# ---- Declare library ----

add_library(
    HadesStsAlignment_lib OBJECT
    source/lib.cpp
)

target_include_directories(
    HadesStsAlignment_lib ${warning_guard}
    PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>"
)

target_compile_features(HadesStsAlignment_lib PUBLIC cxx_std_17)

target_link_libraries(HadesStsAlignment_lib
  PUBLIC HYDRA::all
         ROOT::Core ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Physics ROOT::GenVector
         promille::promille)

# ---- Declare executable ----

add_executable(HadesStsAlignment_exe source/main.cpp)
add_executable(HadesStsAlignment::exe ALIAS HadesStsAlignment_exe)

set_property(TARGET HadesStsAlignment_exe PROPERTY OUTPUT_NAME HadesStsAlignment)
target_compile_features(HadesStsAlignment_exe PRIVATE cxx_std_17)
target_link_libraries(HadesStsAlignment_exe PRIVATE HadesStsAlignment_lib)


add_executable(HadesParamUpdater_exe source/param_updater.cpp)
add_executable(HadesParamUpdater::exe ALIAS HadesStsAlignment_exe)

set_property(TARGET HadesParamUpdater_exe PROPERTY OUTPUT_NAME HadesParamUpdater)
target_compile_features(HadesParamUpdater_exe PRIVATE cxx_std_17)
target_link_libraries(HadesParamUpdater_exe PRIVATE HadesStsAlignment_lib)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Examples ----

if(PROJECT_IS_TOP_LEVEL)
  option(BUILD_EXAMPLES "Build examples tree." "${HadesStsAlignment_DEVELOPER_MODE}")
  if(BUILD_EXAMPLES)
    add_subdirectory(example)
  endif()
endif()

# ---- Hydra module ----

option(BUILD_HYDRA_MODULE "Build HYDRA moduke." ON)
if (BUILD_HYDRA_MODULE)
  if (NOT DEFINED ENV{MYHADDIR})
    message(FATAL_ERROR "MYHADDIR must be exported for building Hydra module.")
  endif()

  configure_file(${CMAKE_SOURCE_DIR}/source/hydra/Makefile.in ${CMAKE_BINARY_DIR}/hydra/Makefile @ONLY)

  add_custom_target(ForwardAligner
    ALL
    COMMAND make -f ${CMAKE_BINARY_DIR}/hydra/Makefile --no-silent
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/source
    COMMENT "Running custom Makefile when hydra module changes"
  )

  install(PROGRAMS ${CMAKE_BINARY_DIR}/hydra/build/lib/libForwardAligner.so
        DESTINATION $ENV{MYHADDIR}/lib
  )
  install(FILES ${CMAKE_SOURCE_DIR}/source/hydra/hforwardaligner.h
        DESTINATION $ENV{MYHADDIR}/include
  )
endif()

# ---- Developer mode ----

if(NOT HadesStsAlignment_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of HadesStsAlignment"
  )
endif()

include(cmake/dev-mode.cmake)
