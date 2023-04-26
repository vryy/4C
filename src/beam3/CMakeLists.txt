set(OBJS_ELE_BEAM3
    # cmake-format: sortable
    beam3.cpp
    beam3_base.cpp
    beam3_discretization_runtime_vtu_output_params.cpp
    beam3_discretization_runtime_vtu_writer.cpp
    beam3_euler_bernoulli.cpp
    beam3_euler_bernoulli_evaluate.cpp
    beam3_euler_bernoulli_input.cpp
    beam3_evaluate.cpp
    beam3_input.cpp
    beam3_kirchhoff.cpp
    beam3_kirchhoff_evaluate.cpp
    beam3_kirchhoff_input.cpp
    beam3_reissner.cpp
    beam3_reissner_evaluate.cpp
    beam3_reissner_input.cpp
    beam3_triad_interpolation.cpp
    beam3_triad_interpolation_local_rotation_vectors.cpp
    )

baci_add_library(ele_beam3 ${OBJS_ELE_BEAM3})

set(LIBS_ELE_BEAM3
    # cmake-format: sortable
    alg_structure_new
    lib
    lib_beaminteraction
    linalg
    )

baci_add_dependency(ele_beam3 ${LIBS_ELE_BEAM3})
