set(OBJS_CARDIOVASCULAR0D
    # cmake-format: sortable
    cardiovascular0d.cpp
    cardiovascular0d_4elementwindkessel.cpp
    cardiovascular0d_arterialproxdist.cpp
    cardiovascular0d_dofset.cpp
    cardiovascular0d_manager.cpp
    cardiovascular0d_resulttest.cpp
    cardiovascular0d_syspulcirculation.cpp
    cardiovascularrespiratory0d_syspulperiphcirculation.cpp
    nox_nln_cardiovascular0d_linearsystem.cpp
    )

baci_add_library(cardiovascular0d ${OBJS_CARDIOVASCULAR0D})

set(LIBS_CARDIOVASCULAR0D
    # cmake-format: sortable
    ele_so
    lib
    mor
    solver
    )

target_link_libraries(cardiovascular0d PUBLIC ${LIBS_CARDIOVASCULAR0D})
