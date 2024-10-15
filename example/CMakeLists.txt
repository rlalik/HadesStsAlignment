cmake_minimum_required(VERSION 3.14)

project(mille_builderExamples CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)

# find ROOT
if (DEFINED ENV{ROOTSYS})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
endif()

find_package(ROOT QUIET REQUIRED COMPONENTS GenVector)

if(PROJECT_IS_TOP_LEVEL)
  find_package(mille_builder REQUIRED)
endif()

add_custom_target(run-examples)

function(add_example NAME)
  add_executable("${NAME}" "${NAME}.cpp")
  target_link_libraries("${NAME}" PRIVATE promille::promille)
  target_compile_features("${NAME}" PRIVATE cxx_std_17)
  add_custom_target("run_${NAME}" COMMAND "${NAME}" VERBATIM)
  add_dependencies("run_${NAME}" "${NAME}")
  add_dependencies(run-examples "run_${NAME}")
endfunction()

add_example(hades_basic_example)

target_sources(hades_basic_example PRIVATE)
target_link_libraries(hades_basic_example PRIVATE MP2::Mille ROOT::GenVector)

add_folders(Example)
