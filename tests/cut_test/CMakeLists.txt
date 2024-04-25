set(OBJS_CUT_TEST
    # cmake-format: sortable
    cut_test_alex_39.cpp
    cut_test_alex_40.cpp
    cut_test_alex_41.cpp
    cut_test_alex_42.cpp
    cut_test_alex_43.cpp
    cut_test_alex_44.cpp
    cut_test_alex_45.cpp
    cut_test_alex_46.cpp
    cut_test_alex_47.cpp
    cut_test_alex_48.cpp
    cut_test_alex_49.cpp
    cut_test_alex_50.cpp
    cut_test_alex_51.cpp
    cut_test_alex_52.cpp
    cut_test_alex_53.cpp
    cut_test_alex_54.cpp
    cut_test_alex_55.cpp
    cut_test_alex_56.cpp
    cut_test_alex_57.cpp
    cut_test_alex_58.cpp
    cut_test_alex_59.cpp
    cut_test_alex_60.cpp
    cut_test_alex_61.cpp
    cut_test_alex_62.cpp
    cut_test_axel_10.cpp
    cut_test_axel_8.cpp
    cut_test_axel_9.cpp
    cut_test_benedikt_1.cpp
    cut_test_double.cpp
    cut_test_facet_failed_cln.cpp
    cut_test_fluidfluid.cpp
    cut_test_generated_1858.cpp
    cut_test_generated_1860.cpp
    cut_test_generated_1890.cpp
    cut_test_generated_1910.cpp
    cut_test_generated_1970.cpp
    cut_test_generated_197489.cpp
    cut_test_generated_2010.cpp
    cut_test_generated_227469.cpp
    cut_test_generated_238343.cpp
    cut_test_generated_238425.cpp
    cut_test_generated_26182.cpp
    cut_test_generated_369096.cpp
    cut_test_generated_41534.cpp
    cut_test_generated_43244.cpp
    cut_test_generated_463638.cpp
    cut_test_generated_622320.cpp
    cut_test_generated_622829.cpp
    cut_test_generated_627558.cpp
    cut_test_generated_6890.cpp
    cut_test_generated_6920.cpp
    cut_test_generated_6923.cpp
    cut_test_generated_7019.cpp
    cut_test_generated_7022.cpp
    cut_test_generated_79216.cpp
    cut_test_geometry.cpp
    cut_test_hex8_twintri.cpp
    cut_test_intersection.cpp
    cut_test_intersection.cpp
    cut_test_levelset.cpp
    cut_test_levelset_no_loc_coord.cpp
    cut_test_loader.cpp
    cut_test_main.cpp
    cut_test_newlineimplementation.cpp
    cut_test_selfcut.cpp
    cut_test_shadan_6.cpp
    cut_test_simple.cpp
    cut_test_triangulateFacet.cpp
    cut_test_units.cpp
    cut_test_utils.cpp
    cut_test_volume.cpp
    )

add_executable(cut_test EXCLUDE_FROM_ALL ${OBJS_CUT_TEST})

# ignore selected cut test
option(IGNORE_UNSTABLE_CUTTEST "Ignore unstable cut test" OFF)
if(IGNORE_UNSTABLE_CUTTEST)
  message(STATUS "Ignoring unstable cut test!")
  target_compile_definitions(cut_test PRIVATE IGNORE_UNSTABLE_CUTTEST)
endif(IGNORE_UNSTABLE_CUTTEST)
set_target_properties(cut_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

baci_set_up_executable(cut_test)
