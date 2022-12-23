set(OBJS_STRUCTURE
    # cmake-format: sortable
    strtimada.cpp
    strtimada_create.cpp
    strtimada_zienxie.cpp
    strtimint.cpp
    strtimint_ab2.cpp
    strtimint_centrdiff.cpp
    strtimint_create.cpp
    strtimint_expl.cpp
    strtimint_expleuler.cpp
    strtimint_gemm.cpp
    strtimint_genalpha.cpp
    strtimint_impl.cpp
    strtimint_impl_nox.cpp
    strtimint_noxgroup.cpp
    strtimint_noxlinsys.cpp
    strtimint_ost.cpp
    strtimint_prestress.cpp
    strtimint_statics.cpp
    stru_aux.cpp
    stru_dyn_nln_drt.cpp
    stru_resulttest.cpp
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
