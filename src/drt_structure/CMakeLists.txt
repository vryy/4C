set(OBJS_DRT_STRUCTURE
    strtimada.cpp
    strtimada_create.cpp
    strtimada_zienxie.cpp
    strtimint.cpp
    strtimint_create.cpp
    strtimint_ab2.cpp
    strtimint_centrdiff.cpp
    strtimint_expl.cpp
    strtimint_expleuler.cpp
    strtimint_gemm.cpp
    strtimint_genalpha.cpp
    strtimint_impl.cpp
    strtimint_impl_nox.cpp
    strtimint_noxgroup.cpp
    strtimint_noxlinsys.cpp
    strtimint_ost.cpp
    strtimint_statics.cpp
    strtimint_prestress.cpp
    stru_aux.cpp
    stru_resulttest.cpp
    stru_dyn_nln_drt.cpp
    )

baci_add_library(alg_structure ${OBJS_DRT_STRUCTURE})
target_link_libraries(alg_structure PUBLIC drt_lib drt_inpar drt_adapter)
