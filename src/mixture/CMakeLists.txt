set(OBJS_MIXTURE
    # cmake-format: sortable
    growth_evolution.H
    growth_evolution_linear_cauchy.H
    mixture_constituent.cpp
    mixture_constituent_elasthyper.cpp
    mixture_constituent_elasthyper_damage.cpp
    mixture_constituent_elasthyper_elastin_membrane.cpp
    mixture_constituent_elasthyperbase.cpp
    mixture_constituent_remodelfiber_expl.cpp
    mixture_constituent_remodelfiber_impl.cpp
    mixture_constituent_remodelfiber_lib.cpp
    mixture_constituent_remodelfiber_material.cpp
    mixture_constituent_remodelfiber_material_exponential.cpp
    mixture_constituent_remodelfiber_material_exponential_active.cpp
    mixture_constituent_solidmaterial.cpp
    mixture_growth_strategy.cpp
    mixture_growth_strategy_anisotropic.cpp
    mixture_growth_strategy_isotropic.cpp
    mixture_growth_strategy_stiffness.cpp
    mixture_prestress_strategy.cpp
    mixture_prestress_strategy_constant.cpp
    mixture_prestress_strategy_isocyl.cpp
    mixture_prestress_strategy_iterative.cpp
    mixture_rule.cpp
    mixture_rule_growthremodel.cpp
    mixture_rule_simple.cpp
    remodelfiber.cpp
    )

baci_add_library(mixture ${OBJS_MIXTURE})

set(LIBS_MIXTURE
    # cmake-format: sortable
    mat
    )

target_link_libraries(mixture PUBLIC ${LIBS_MIXTURE})
