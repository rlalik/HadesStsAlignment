# Parent project does not export its library target, so this CML implicitly
# depends on being added from it, i.e. the testing is done only from the build
# tree and is not feasible from an install location

project(HadesStsAlignmentTests LANGUAGES CXX)

find_package(GTest)

if(TARGET GTest::gtest_main)
  # new cmake provided
  set(GTEST_TRG GTest::gtest_main)
else()
  # old cmake and pc files
  set(GTEST_TRG gtest gtest_main)
endif()

# ---- Tests ----

add_executable(HadesStsAlignment_test
    source/test_straw_residuals_model.cpp
    source/test_derivatives.cpp
)
target_link_libraries(HadesStsAlignment_test PRIVATE HadesStsAlignment_lib ROOT::GenVector promille::promille ${GTEST_TRG} ${FMT_TARGET} tabulate::tabulate)
target_compile_features(HadesStsAlignment_test PRIVATE cxx_std_17)

add_test(NAME HadesStsAlignment_test COMMAND HadesStsAlignment_test)

# ---- End-of-file commands ----

add_folders(Test)
