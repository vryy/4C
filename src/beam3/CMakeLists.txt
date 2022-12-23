set(OBJS_ELE_BEAM3
    # cmake-format: sortable
    beam3.cpp
    beam3_base.cpp
    beam3_evaluate.cpp
    beam3_input.cpp
    beam3eb.cpp
    beam3eb_evaluate.cpp
    beam3eb_input.cpp
    beam3k.cpp
    beam3k_evaluate.cpp
    beam3k_input.cpp
    beam3r.cpp
    beam3r_evaluate.cpp
    beam3r_input.cpp
    beam_discretization_runtime_vtu_output_params.cpp
    beam_discretization_runtime_vtu_writer.cpp
    triad_interpolation.cpp
    triad_interpolation_local_rotation_vectors.cpp
    )

baci_add_library(ele_beam3 ${OBJS_ELE_BEAM3})

set(LIBS_ELE_BEAM3
    # cmake-format: sortable
    alg_structure_new
    fem_general
    lib
    lib_beaminteraction
    linalg
    )

target_link_libraries(ele_beam3 PUBLIC ${LIBS_ELE_BEAM3})
