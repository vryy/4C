set(OBJS_SSI
    # cmake-format: sortable
    ssi_base.cpp
    ssi_clonestrategy.cpp
    ssi_coupling.cpp
    ssi_dyn.cpp
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
    ssi_str_model_evaluator_base.cpp
    ssi_str_model_evaluator_partitioned.cpp
    ssi_utils.cpp
    )

baci_add_library(alg_ssi ${OBJS_SSI})

set(LIBS_ALG_SSI
    # cmake-format: sortable
    ele_membrane
    fem_general
    inpar
    lib
    )

baci_add_dependency(alg_ssi ${LIBS_ALG_SSI})
