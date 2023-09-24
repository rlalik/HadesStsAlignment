# Parent project does not export its library target, so this CML implicitly
# depends on being added from it, i.e. the testing is done only from the build
# tree and is not feasible from an install location

project(HadesStsAlignmentTests LANGUAGES CXX)

# ---- Tests ----

add_executable(HadesStsAlignment_test source/HadesStsAlignment_test.cpp)
target_link_libraries(HadesStsAlignment_test PRIVATE HadesStsAlignment_lib)
target_compile_features(HadesStsAlignment_test PRIVATE cxx_std_17)

add_test(NAME HadesStsAlignment_test COMMAND HadesStsAlignment_test)

# ---- End-of-file commands ----

add_folders(Test)
