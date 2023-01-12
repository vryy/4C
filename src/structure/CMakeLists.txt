set(OBJS_STRUCTURE
    # cmake-format: sortable
    structure_aux.cpp
    structure_dyn_nln_drt.cpp
    structure_resulttest.cpp
    structure_timada.cpp
    structure_timada_create.cpp
    structure_timada_zienxie.cpp
    structure_timint.cpp
    structure_timint_ab2.cpp
    structure_timint_centrdiff.cpp
    structure_timint_create.cpp
    structure_timint_expl.cpp
    structure_timint_expleuler.cpp
    structure_timint_gemm.cpp
    structure_timint_genalpha.cpp
    structure_timint_impl.cpp
    structure_timint_impl_nox.cpp
    structure_timint_noxgroup.cpp
    structure_timint_noxlinsys.cpp
    structure_timint_ost.cpp
    structure_timint_prestress.cpp
    structure_timint_statics.cpp
    )

baci_add_library(alg_structure ${OBJS_STRUCTURE})

set(LIBS_ALG_STRUCTURE
    # cmake-format: sortable
    adapter
    cardiovascular0d
    inpar
    levelset
    lib
    timestepping
    )

target_link_libraries(alg_structure PUBLIC ${LIBS_ALG_STRUCTURE})
