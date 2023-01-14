set(OBJS_POROFLUIDMULTIPHASE_LIB
    # cmake-format: sortable
    porofluidmultiphase_ele.cpp
    porofluidmultiphase_ele_boundary_calc.cpp
    porofluidmultiphase_ele_boundary_evaluate.cpp
    porofluidmultiphase_ele_boundary_factory.cpp
    porofluidmultiphase_ele_calc.cpp
    porofluidmultiphase_ele_calc_utils.cpp
    porofluidmultiphase_ele_evaluate.cpp
    porofluidmultiphase_ele_factory.cpp
    porofluidmultiphase_ele_parameter.cpp
    porofluidmultiphase_ele_porofluid_evaluator.cpp
    porofluidmultiphase_ele_porofluid_phasemanager.cpp
    porofluidmultiphase_ele_porofluid_variablemanager.cpp
    )

baci_add_library(ele_porofluidmultiphase ${OBJS_POROFLUIDMULTIPHASE_LIB})

set(LIBS_ALG_POROFLUIDMULTIPHASE
    # cmake-format: sortable
    lib
    )

target_link_libraries(ele_porofluidmultiphase PUBLIC ${LIBS_ALG_POROFLUIDMULTIPHASE})
