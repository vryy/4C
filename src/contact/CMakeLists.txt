set(OBJS_CONTACT
    # cmake-format: sortable
    contact_abstract_strategy.cpp
    contact_coupling2d.cpp
    contact_coupling3d.cpp
    contact_ehl_integrator.cpp
    contact_integrator.cpp
    contact_integrator_assemble.cpp
    contact_integrator_factory.cpp
    contact_interface.cpp
    contact_interface_assemble.cpp
    contact_interface_roundrobin.cpp
    contact_interface_tools.cpp
    contact_interpolator.cpp
    contact_lagrange_strategy.cpp
    contact_line_coupling.cpp
    contact_manager.cpp
    contact_meshtying_abstract_strategy.cpp
    contact_meshtying_contact_bridge.cpp
    contact_meshtying_lagrange_strategy.cpp
    contact_meshtying_manager.cpp
    contact_meshtying_noxinterface.cpp
    contact_meshtying_penalty_strategy.cpp
    contact_meshtying_poro_lagrange_strategy.cpp
    contact_meshtying_strategy_factory.cpp
    contact_monocoupled_lagrange_strategy.cpp
    contact_nitsche_integrator.cpp
    contact_nitsche_integrator_fpi.cpp
    contact_nitsche_integrator_fsi.cpp
    contact_nitsche_integrator_poro.cpp
    contact_nitsche_integrator_ssi.cpp
    contact_nitsche_integrator_ssi_elch.cpp
    contact_nitsche_integrator_tsi.cpp
    contact_nitsche_strategy.cpp
    contact_nitsche_strategy_fpi.cpp
    contact_nitsche_strategy_fsi.cpp
    contact_nitsche_strategy_poro.cpp
    contact_nitsche_strategy_ssi.cpp
    contact_nitsche_strategy_ssi_elch.cpp
    contact_nitsche_strategy_tsi.cpp
    contact_nitsche_utils.cpp
    contact_nox_nln_meshtying_linearsystem.cpp
    contact_noxinterface.cpp
    contact_penalty_strategy.cpp
    contact_poro_lagrange_strategy.cpp
    contact_selfcontact_binarytree.cpp
    contact_selfcontact_binarytree_unbiased.cpp
    contact_strategy_factory.cpp
    contact_tsi_interface.cpp
    contact_tsi_lagrange_strategy.cpp
    contact_utils.cpp
    contact_utils_parallel.cpp
    contact_wear_interface.cpp
    contact_wear_interface_tools.cpp
    contact_wear_lagrange_strategy.cpp
    )

baci_add_library(alg_contact ${OBJS_CONTACT})

set(LIBS_ALG_CONTACT
    # cmake-format: sortable
    alg_contact_aug
    alg_xstructure
    ele_scatra
    ele_so
    ele_so_plast
    lib_mortar
    )

baci_add_dependency(alg_contact ${LIBS_ALG_CONTACT})

set(OBJS_CONTACT_LIB
    # cmake-format: sortable
    contact_analytical.cpp
    contact_element.cpp
    contact_friction_node.cpp
    contact_friction_node.cpp
    contact_node.cpp
    )

baci_add_library(lib_contact ${OBJS_CONTACT_LIB})

set(LIBS_LIB_CONTACT
    # cmake-format: sortable
    alg_contact_aug
    lib_mortar
    linalg
    rebalance
    )

baci_add_dependency(lib_contact ${LIBS_LIB_CONTACT})
