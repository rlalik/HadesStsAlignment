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

find_package(HYDRA)
find_package(ROOT)

include(FetchContent)
FetchContent_Declare(MillePede
  GIT_REPOSITORY https://gitlab.desy.de/claus.kleinwort/millepede-ii.git
  GIT_TAG main
)
FetchContent_MakeAvailable(MillePede)

add_library(Mille STATIC
    ${millepede_SOURCE_DIR}/Mille.cc
)
target_include_directories(Mille PUBLIC ${millepede_SOURCE_DIR})


add_custom_target(pede
  COMMAND make
  WORKING_DIRECTORY ${millepede_SOURCE_DIR}
)

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
  PUBLIC ROOT::Core ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Physics ROOT::GenVector
         HYDRA::Hydra HYDRA::Forward HYDRA::Sts HYDRA::FRpc HYDRA::Dst HYDRA::Emc
         HYDRA::Tof HYDRA::Rpc HYDRA::iTof HYDRA::Start HYDRA::Rich HYDRA::Mdc
         HYDRA::MdcTrackD HYDRA::MdcTrackG HYDRA::PionTracker HYDRA::Wall
         HYDRA::Kalman HYDRA::Shower HYDRA::Particle
         Mille)

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
