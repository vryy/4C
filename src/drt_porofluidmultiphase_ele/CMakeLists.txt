set(OBJS_DRT_POROFLUIDMULTIPHASE_LIB
    # cmake-format: sortable
    porofluid_evaluator.cpp
    porofluid_phasemanager.cpp
    porofluid_variablemanager.cpp
    porofluidmultiphase_ele.cpp
    porofluidmultiphase_ele_boundary_calc.cpp
    porofluidmultiphase_ele_boundary_evaluate.cpp
    porofluidmultiphase_ele_boundary_factory.cpp
    porofluidmultiphase_ele_calc.cpp
    porofluidmultiphase_ele_calc_utils.cpp
    porofluidmultiphase_ele_evaluate.cpp
    porofluidmultiphase_ele_factory.cpp
    porofluidmultiphase_ele_parameter.cpp
    )

baci_add_library(ele_porofluidmultiphase ${OBJS_DRT_POROFLUIDMULTIPHASE_LIB})

set(LIBS_ALG_POROFLUIDMULTIPHASE
    # cmake-format: sortable
    drt_lib
    )

target_link_libraries(ele_porofluidmultiphase PUBLIC ${LIBS_ALG_POROFLUIDMULTIPHASE})
