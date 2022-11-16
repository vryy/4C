set(OBJS_DRT_CONTACT
    contact_abstract_strategy.cpp
    contact_noxinterface.cpp
    contact_coupling2d.cpp
    contact_coupling3d.cpp
    contact_line_coupling.cpp
    contact_integrator.cpp
    contact_nitsche_integrator.cpp
    contact_nitsche_integrator_poro.cpp
    contact_nitsche_integrator_fsi.cpp
    contact_nitsche_integrator_fpi.cpp
    contact_nitsche_integrator_ssi.cpp
    contact_nitsche_integrator_ssi_elch.cpp
    contact_nitsche_integrator_tsi.cpp
    contact_ehl_integrator.cpp
    contact_integrator_factory.cpp
    contact_integrator_assemble.cpp
    contact_interface.cpp
    contact_interface_roundrobin.cpp
    contact_interface_tools.cpp
    contact_interface_assemble.cpp
    contact_interpolator.cpp
    contact_lagrange_strategy.cpp
    contact_monocoupled_lagrange_strategy.cpp
    contact_poro_lagrange_strategy.cpp
    contact_tsi_lagrange_strategy.cpp
    contact_tsi_interface.cpp
    contact_manager.cpp
    contact_strategy_factory.cpp
    contact_penalty_strategy.cpp
    contact_nitsche_strategy.cpp
    contact_nitsche_strategy_poro.cpp
    contact_nitsche_strategy_fsi.cpp
    contact_nitsche_strategy_fpi.cpp
    contact_nitsche_strategy_ssi.cpp
    contact_nitsche_strategy_ssi_elch.cpp
    contact_nitsche_strategy_tsi.cpp
    contact_utils.cpp
    contact_utils_parallel.cpp
    contact_wear_interface.cpp
    contact_wear_interface_tools.cpp
    contact_wear_lagrange_strategy.cpp
    meshtying_strategy_factory.cpp
    meshtying_abstract_strategy.cpp
    meshtying_noxinterface.cpp
    meshtying_contact_bridge.cpp
    meshtying_lagrange_strategy.cpp
    meshtying_poro_lagrange_strategy.cpp
    meshtying_manager.cpp
    meshtying_penalty_strategy.cpp
    selfcontact_binarytree.cpp
    selfcontact_binarytree_unbiased.cpp
    nox_nln_meshtying_linearsystem.cpp
    )

baci_add_library(alg_contact ${OBJS_DRT_CONTACT})

target_link_libraries(alg_contact PUBLIC drt_lib_mortar)

set(OBJS_DRT_CONTACT_LIB
    contact_analytical.cpp contact_element.cpp contact_node.cpp friction_node.cpp
    )

baci_add_library(drt_lib_contact ${OBJS_DRT_CONTACT_LIB})

target_link_libraries(drt_lib_contact PUBLIC drt_lib_mortar linalg)
