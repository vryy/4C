set(OBJS_PARTICLE_INTERACTION
    # cmake-format: sortable
    particle_interaction_base.cpp
    particle_interaction_dem.cpp
    particle_interaction_dem_adhesion.cpp
    particle_interaction_dem_adhesion_law.cpp
    particle_interaction_dem_adhesion_surface_energy.cpp
    particle_interaction_dem_contact.cpp
    particle_interaction_dem_contact_normal.cpp
    particle_interaction_dem_contact_rolling.cpp
    particle_interaction_dem_contact_tangential.cpp
    particle_interaction_dem_history_pairs.cpp
    particle_interaction_dem_neighbor_pairs.cpp
    particle_interaction_material_handler.cpp
    particle_interaction_runtime_writer.cpp
    particle_interaction_sph.cpp
    particle_interaction_sph_artificialviscosity.cpp
    particle_interaction_sph_boundary_particle.cpp
    particle_interaction_sph_density.cpp
    particle_interaction_sph_density_correction.cpp
    particle_interaction_sph_equationofstate.cpp
    particle_interaction_sph_equationofstate_bundle.cpp
    particle_interaction_sph_heatloss_evaporation.cpp
    particle_interaction_sph_heatsource.cpp
    particle_interaction_sph_kernel.cpp
    particle_interaction_sph_momentum.cpp
    particle_interaction_sph_momentum_formulation.cpp
    particle_interaction_sph_neighbor_pairs.cpp
    particle_interaction_sph_open_boundary.cpp
    particle_interaction_sph_phase_change.cpp
    particle_interaction_sph_pressure.cpp
    particle_interaction_sph_rigid_particle_contact.cpp
    particle_interaction_sph_surface_tension.cpp
    particle_interaction_sph_surface_tension_barrier_force.cpp
    particle_interaction_sph_surface_tension_interface_viscosity.cpp
    particle_interaction_sph_surface_tension_recoilpressure_evaporation.cpp
    particle_interaction_sph_temperature.cpp
    particle_interaction_sph_virtual_wall_particle.cpp
    )

baci_add_library(alg_particle_interaction ${OBJS_PARTICLE_INTERACTION})

set(LIBS_ALG_PARTICLE_INTERACTION
    # cmake-format: sortable
    alg_particle_wall
    inpar
    particleengine
    )

target_link_libraries(alg_particle_interaction PUBLIC ${LIBS_ALG_PARTICLE_INTERACTION})
