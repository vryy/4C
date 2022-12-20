set(OBJS_FBI
    # cmake-format: sortable
    ad_fbi_constraintbridge.cpp
    ad_fbi_constraintbridge_penalty.cpp
    beam_to_fluid_assembly_manager_factory.cpp
    beam_to_fluid_meshtying_pair_base.cpp
    beam_to_fluid_meshtying_pair_factory.cpp
    beam_to_fluid_meshtying_pair_gauss_point.cpp
    beam_to_fluid_meshtying_pair_mortar.cpp
    beam_to_fluid_meshtying_params.cpp
    beam_to_fluid_meshtying_vtk_output_params.cpp
    beam_to_fluid_meshtying_vtk_output_writer.cpp
    beam_to_fluid_mortar_manager.cpp
    constraintenforcer_fbi.cpp
    constraintenforcer_fbi_factory.cpp
    constraintenforcer_fbi_penalty.cpp
    fbi_calc_utils.cpp
    fluid_assembly_strategy.cpp
    fluidblockmatrix_assembly_strategy.cpp
    immersed_geometry_coupler_fbi.cpp
    immersed_geometry_coupler_fbi_binning.cpp
    immersedcoupler_factory.cpp
    partitioned_penaltycoupling_assembly_manager.cpp
    partitioned_penaltycoupling_assembly_manager_direct.cpp
    partitioned_penaltycoupling_assembly_manager_indirect.cpp
    )

baci_add_library(alg_fbi ${OBJS_FBI})

set(LIBS_ALG_FBI
    # cmake-format: sortable
    lib_beaminteraction
    linalg
    )

target_link_libraries(alg_fbi PUBLIC ${LIBS_ALG_FBI})
