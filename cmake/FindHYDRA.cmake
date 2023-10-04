# cmake-format: off
# * Find HYDRA instalation
#  - Original author: F.Uhlig@gsi.de (fairroot.gsi.de)
#  - Modifications: Rafal.Lalik@uj.edu.pl
# cmake-format: on

# Different set of modules for different Hydra generation
# cmake-format: off
set(HYDRA1_MODULES
    Alignment Dst HadesGo4 Hodo Hydra Hyp Kick MdcGarfield MdcPid Mdc MdcTrackD MdcTrackG MdcTrackS
    MdcUtil Pairs Particle PhyAna Pid PidUtil QA Revt Rich RichUtil Rpc Shower ShowerTofino
    ShowerUtil Simulation Start Tofino Tof TofUtil Tools Trigger TriggerUtil Wall Ora OraSim OraUtil)

set(HYDRA2_MODULES
    Alignment Dst Emc EventDisplay Forward FRpc Hydra iTof Kalman MdcGarfield Mdc MdcTrackD
    MdcTrackG MdcUtil Online Particle PionTracker QA Revt Rich Rpc Shower ShowerUtil Simulation
    Start Sts Tof Tools Wall Ora OraSim OraUtil)
# cmake-format: on

set(HYDRA_GENERATION
    AUTO
    CACHE STRING "Hydra software generation")
set_property(CACHE HYDRA_GENERATION PROPERTY STRINGS AUTO;HYDRA2;HYDRA)

message(STATUS "Looking for Hydra version ${HYDRA_GENERATION}...")

if(NOT HADDIR)
  set(INTHADDIR $ENV{HADDIR})
endif(NOT HADDIR)

if(NOT INTHADDIR)
  message(FATAL_ERROR "HADDIR is not set. Please set HADDIR first.")
endif(NOT INTHADDIR)

if(NOT MYHADDIR)
  set(MYHADDIR $ENV{MYHADDIR})
endif(NOT MYHADDIR)

set(HYDRA_PACKAGE_SEARCHPATH ${INTHADDIR}/lib)

set(HYDRA_DEFINITIONS "")

set(HYDRA_INSTALLED_VERSION_TOO_OLD FALSE)

set(HYDRA_MAIN_LIBRARY HYDRA_MAIN_LIBRARY-NOTFOUND)

find_library(
  HYDRA_MAIN_LIBRARY
  NAMES Hydra
  PATHS ${HYDRA_PACKAGE_SEARCHPATH}
  NO_DEFAULT_PATH)

if(${HYDRA_MAIN_LIBRARY} MATCHES "HYDRA_MAIN_LIBRARY-NOTFOUND")
  message(
    STATUS "Hydra library not found. Please check your Hydra installation.")
  set(HYDRA_FOUND FALSE)
else(${HYDRA_MAIN_LIBRARY} MATCHES "HYDRA_MAIN_LIBRARY-NOTFOUND")
  message(STATUS "Looking for Hydra... - found ${INTHADDIR}")
  message(STATUS "   MYHADDIR = ${MYHADDIR}")
  set(HYDRA_FOUND TRUE)
endif(${HYDRA_MAIN_LIBRARY} MATCHES "HYDRA_MAIN_LIBRARY-NOTFOUND")

if(HYDRA_FOUND)
  set(HYDRA_LIBRARY_DIR ${INTHADDIR}/lib)
  set(HYDRA_INCLUDE_DIR ${INTHADDIR}/include)

  # if option AUTO search for libPid library which is available
  # in Hydra1, if not found mark as Hydra2
  if(HYDRA_GENERATION STREQUAL "AUTO")
    find_library(
      HYDRA_HYDRA1_TEST_LIBRARY
      NAMES Pid
      PATHS ${HYDRA_PACKAGE_SEARCHPATH}
      NO_DEFAULT_PATH)

      if (HYDRA_HYDRA1_TEST_LIBRARY)
        set(HYDRA_GENERATION "HYDRA")
      else()
        set(HYDRA_GENERATION "HYDRA2")
      endif()
      unset(HYDRA_HYDRA1_TEST_LIBRARY CACHE)
  endif()

  if(HYDRA_GENERATION STREQUAL "HYDRA2")
    set(HYDRA_AVAILABLE_MODULES ${HYDRA2_MODULES})
  elseif(HYDRA_GENERATION STREQUAL "HYDRA")
    set(HYDRA_AVAILABLE_MODULES ${HYDRA1_MODULES})
  endif()

  set(HYDRA_LIBRARIES)
  foreach(_cpt ${HYDRA_AVAILABLE_MODULES})
    find_library(
      HYDRA_${_cpt}_LIBRARY
      NAMES ${_cpt} lib{_cpt}
      HINTS ${HYDRA_LIBRARY_DIR})
    if(HYDRA_${_cpt}_LIBRARY)
      mark_as_advanced(HYDRA_${_cpt}_LIBRARY)
      list(APPEND HYDRA_LIBRARIES ${HYDRA_${_cpt}_LIBRARY})
      add_library(${_cpt} SHARED IMPORTED)
      target_include_directories(${_cpt} INTERFACE ${HYDRA_INCLUDE_DIR})
      target_link_directories(${_cpt} INTERFACE ${HYDRA_LIBRARY_DIR})
      set_target_properties(
        ${_cpt} PROPERTIES IMPORTED_CONFIGURATIONS RELEASE
        IMPORTED_LOCATION_RELEASE "${HYDRA_${_cpt}_LIBRARY}"
                   IMPORTED_SONAME_RELEASE "lib${_cpt}.so")
      add_library(HYDRA::${_cpt} ALIAS ${_cpt})
    else()
      if(HYDRA_FIND_REQUIRED_${_cpt})
        message(FATAL_ERROR "Hydra component ${_cpt} not found")
      elseif(${_cpt} IN_LIST HYDRA_FIND_COMPONENTS AND NOT HYDRA_FIND_QUIETLY)
        message(WARNING " Hydra component ${_cpt} not found")
      endif()
    endif()
  endforeach()
  if(HYDRA_LIBRARIES)
    list(REMOVE_DUPLICATES HYDRA_LIBRARIES)
  endif()

  # Make variables changeble to the advanced user
  mark_as_advanced(HYDRA_LIBRARY_DIR HYDRA_INCLUDE_DIR HYDRA_DEFINITIONS)

  # Set HYDRA_INCLUDES
  set(HYDRA_INCLUDES ${HYDRA_INCLUDE_DIR})

  set(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${HYDRA_LIBRARY_DIR})
endif(HYDRA_FOUND)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HYDRA DEFAULT_MSG HYDRA_LIBRARY_DIR HYDRA_INCLUDE_DIR)
