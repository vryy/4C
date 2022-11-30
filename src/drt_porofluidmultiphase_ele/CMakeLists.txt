set(OBJS_DRT_POROFLUIDMULTIPHASE_LIB
    porofluidmultiphase_ele.cpp
    porofluidmultiphase_ele_calc.cpp
    porofluidmultiphase_ele_calc_utils.cpp
    porofluidmultiphase_ele_factory.cpp
    porofluidmultiphase_ele_parameter.cpp
    porofluidmultiphase_ele_evaluate.cpp
    porofluid_phasemanager.cpp
    porofluid_variablemanager.cpp
    porofluid_evaluator.cpp
    porofluidmultiphase_ele_boundary_evaluate.cpp
    porofluidmultiphase_ele_boundary_factory.cpp
    porofluidmultiphase_ele_boundary_calc.cpp
    )

baci_add_library(ele_porofluidmultiphase ${OBJS_DRT_POROFLUIDMULTIPHASE_LIB})
target_link_libraries(ele_porofluidmultiphase PUBLIC drt_lib)
