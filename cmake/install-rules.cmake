install(
    TARGETS HadesStsAlignment_exe
    RUNTIME COMPONENT HadesStsAlignment_Runtime
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
