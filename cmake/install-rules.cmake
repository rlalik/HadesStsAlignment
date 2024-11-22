if (BUILD_APPS)
  install(
    TARGETS HadesStsAlignment_exe HadesParamUpdater_exe
    RUNTIME COMPONENT HadesStsAlignment_Runtime
  )
endif()

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
