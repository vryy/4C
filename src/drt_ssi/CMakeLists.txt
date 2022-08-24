set(OBJS_DRT_SSI
    ssi_dyn.cpp
    ssi_base.cpp
    ssi_clonestrategy.cpp
    ssi_coupling.cpp
    ssi_manifold_utils.cpp
    ssi_monolithic.cpp
    ssi_monolithic_assemble_strategy.cpp
    ssi_monolithic_contact_strategy.cpp
    ssi_monolithic_convcheck_strategies.cpp
    ssi_monolithic_dbc_handler.cpp
    ssi_monolithic_evaluate_OffDiag.cpp
    ssi_monolithic_meshtying_strategy.cpp
    ssi_partitioned.cpp
    ssi_partitioned_1wc.cpp
    ssi_partitioned_2wc.cpp
    ssi_resulttest.cpp
    ssi_utils.cpp
    ssi_str_model_evaluator_base.cpp
    ssi_str_model_evaluator_partitioned.cpp
    )

baci_add_library(alg_ssi ${OBJS_DRT_SSI})
target_link_libraries(alg_ssi PUBLIC drt_inpar drt_lib)
