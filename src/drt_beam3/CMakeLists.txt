set(OBJS_ELE_BEAM3
    beam_discretization_runtime_vtu_writer.cpp
    beam_discretization_runtime_vtu_output_params.cpp
    triad_interpolation.cpp
    triad_interpolation_local_rotation_vectors.cpp
    beam3_base.cpp
    beam3.cpp
    beam3_evaluate.cpp
    beam3_input.cpp
    beam3r.cpp
    beam3r_evaluate.cpp
    beam3r_input.cpp
    beam3eb.cpp
    beam3eb_evaluate.cpp
    beam3eb_input.cpp
    beam3k.cpp
    beam3k_evaluate.cpp
    beam3k_input.cpp
    )

baci_add_library(ele_beam3 ${OBJS_ELE_BEAM3})
target_link_libraries(ele_beam3 PUBLIC drt_lib linalg drt_fem_general)
