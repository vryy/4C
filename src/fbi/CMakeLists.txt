set(OBJS_FBI
    # cmake-format: sortable
    fbi_adapter_constraintbridge.cpp
    fbi_adapter_constraintbridge_penalty.cpp
    fbi_beam_to_fluid_assembly_manager_factory.cpp
    fbi_beam_to_fluid_meshtying_pair_base.cpp
    fbi_beam_to_fluid_meshtying_pair_factory.cpp
    fbi_beam_to_fluid_meshtying_pair_gauss_point.cpp
    fbi_beam_to_fluid_meshtying_pair_mortar.cpp
    fbi_beam_to_fluid_meshtying_params.cpp
    fbi_beam_to_fluid_meshtying_vtk_output_params.cpp
    fbi_beam_to_fluid_meshtying_vtk_output_writer.cpp
    fbi_beam_to_fluid_mortar_manager.cpp
    fbi_calc_utils.cpp
    fbi_constraintenforcer.cpp
    fbi_constraintenforcer_factory.cpp
    fbi_constraintenforcer_penalty.cpp
    fbi_fluid_assembly_strategy.cpp
    fbi_fluidblockmatrix_assembly_strategy.cpp
    fbi_immersed_geometry_coupler.cpp
    fbi_immersed_geometry_coupler_binning.cpp
    fbi_immersedcoupler_factory.cpp
    fbi_partitioned_penaltycoupling_assembly_manager.cpp
    fbi_partitioned_penaltycoupling_assembly_manager_direct.cpp
    fbi_partitioned_penaltycoupling_assembly_manager_indirect.cpp
    )

baci_add_library(alg_fbi ${OBJS_FBI})

set(LIBS_ALG_FBI
    # cmake-format: sortable
    lib_beaminteraction
    linalg
    )

baci_add_dependency(alg_fbi ${LIBS_ALG_FBI})
