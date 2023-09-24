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

# ---- Declare executable ----

add_executable(HadesStsAlignment_exe source/main.cpp)
add_executable(HadesStsAlignment::exe ALIAS HadesStsAlignment_exe)

set_property(TARGET HadesStsAlignment_exe PROPERTY OUTPUT_NAME HadesStsAlignment)

target_compile_features(HadesStsAlignment_exe PRIVATE cxx_std_17)

target_link_libraries(HadesStsAlignment_exe PRIVATE HadesStsAlignment_lib)

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
