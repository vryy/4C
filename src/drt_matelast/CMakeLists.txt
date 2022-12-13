set(OBJS_DRT_MATELAST
    # cmake-format: sortable
    elast_aniso_structuraltensor_strategy.cpp
    elast_anisoactivestress_evolution.cpp
    elast_coup13apow.cpp
    elast_coup1pow.cpp
    elast_coup2pow.cpp
    elast_coup3pow.cpp
    elast_coupanisoexpo.cpp
    elast_coupanisoexpoactive.cpp
    elast_coupanisoexpobase.cpp
    elast_coupanisoexposhear.cpp
    elast_coupanisoexpotwocoup.cpp
    elast_coupanisoneohooke.cpp
    elast_coupanisoneohooke_VarProp.cpp
    elast_coupanisopow.cpp
    elast_coupblatzko.cpp
    elast_coupexppol.cpp
    elast_couplogmixneohooke.cpp
    elast_couplogneohooke.cpp
    elast_coupmooneyrivlin.cpp
    elast_coupneohooke.cpp
    elast_coupSaintVenantKirchhoff.cpp
    elast_coupsimopister.cpp
    elast_couptransverselyisotropic.cpp
    elast_coupvarga.cpp
    elast_iso1pow.cpp
    elast_iso2pow.cpp
    elast_isoanisoexpo.cpp
    elast_isoexpopow.cpp
    elast_isomooneyrivlin.cpp
    elast_isomuscle_blemker.cpp
    elast_isoneohooke.cpp
    elast_isoogden.cpp
    elast_isotestmaterial.cpp
    elast_isovarga.cpp
    elast_isovolaaagasser.cpp
    elast_isovolHUdependentneohooke.cpp
    elast_isoyeoh.cpp
    elast_remodelfiber.cpp
    elast_summand.cpp
    elast_vologden.cpp
    elast_volpenalty.cpp
    elast_volpow.cpp
    elast_volsussmanbathe.cpp
    visco_coupmyocard.cpp
    visco_fract.cpp
    visco_generalizedgenmax.cpp
    visco_genmax.cpp
    visco_isoratedep.cpp
    )

baci_add_library(drt_matelast ${OBJS_DRT_MATELAST})

set(LIBS_DRT_MATELAST
    # cmake-format: sortable
    drt_inpar
    drt_mat
    )

target_link_libraries(drt_matelast PUBLIC ${LIBS_DRT_MATELAST})
