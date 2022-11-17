set(OBJS_DRT_CARDIOVASCULAR0D
    cardiovascular0d.cpp
    cardiovascular0d_manager.cpp
    cardiovascular0d_dofset.cpp
    cardiovascular0d_resulttest.cpp
    cardiovascular0d_4elementwindkessel.cpp
    cardiovascular0d_arterialproxdist.cpp
    cardiovascular0d_syspulcirculation.cpp
    cardiovascularrespiratory0d_syspulperiphcirculation.cpp
    nox_nln_cardiovascular0d_linearsystem.cpp
    )

baci_add_library(drt_cardiovascular0d ${OBJS_DRT_CARDIOVASCULAR0D})
target_link_libraries(drt_cardiovascular0d PUBLIC drt_lib)
