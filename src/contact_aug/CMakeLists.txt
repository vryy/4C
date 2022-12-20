set(OBJS_CONTACT_AUG
    # cmake-format: sortable
    contact_aug_active_set.cpp
    contact_aug_combo_strategy.cpp
    contact_aug_combo_strategy_switching.cpp
    contact_aug_integrator_policy.cpp
    contact_aug_integrator_policy_debug.cpp
    contact_aug_interface_assembly.cpp
    contact_aug_lagrange_interface.cpp
    contact_aug_lagrange_multiplier_function.cpp
    contact_aug_lagrange_strategy.cpp
    contact_aug_parallel_distribution_controller.cpp
    contact_aug_penalty_update.cpp
    contact_aug_plot.cpp
    contact_aug_potential.cpp
    contact_aug_projector.cpp
    contact_aug_steepest_ascent_interface.cpp
    contact_aug_steepest_ascent_sp_strategy.cpp
    contact_aug_steepest_ascent_strategy.cpp
    contact_aug_timemonitor.cpp
    contact_augmented_integrator.cpp
    contact_augmented_integrator_evaluate.cpp
    contact_augmented_interface.cpp
    contact_augmented_interface_tools.cpp
    contact_augmented_strategy.cpp
    contact_augmented_strategy_tools.cpp
    contact_integrator_utils.cpp
    nox_nln_contact_linearsystem.cpp
    )

baci_add_library(alg_contact_aug ${OBJS_CONTACT_AUG})

set(LIBS_ALG_CONTACT_AUG
    # cmake-format: sortable
    alg_contact
    )

target_link_libraries(alg_contact_aug PUBLIC ${LIBS_ALG_CONTACT_AUG})
