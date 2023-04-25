set(OBJS_MATELAST
    # cmake-format: sortable
    matelast_aniso_structuraltensor_strategy.cpp
    matelast_anisoactivestress_evolution.cpp
    matelast_coup13apow.cpp
    matelast_coup1pow.cpp
    matelast_coup2pow.cpp
    matelast_coup3pow.cpp
    matelast_coupanisoexpo.cpp
    matelast_coupanisoexpoactive.cpp
    matelast_coupanisoexpobase.cpp
    matelast_coupanisoexposhear.cpp
    matelast_coupanisoexpotwocoup.cpp
    matelast_coupanisoneohooke.cpp
    matelast_coupanisoneohooke_VarProp.cpp
    matelast_coupanisopow.cpp
    matelast_coupblatzko.cpp
    matelast_coupexppol.cpp
    matelast_couplogmixneohooke.cpp
    matelast_couplogneohooke.cpp
    matelast_coupmooneyrivlin.cpp
    matelast_coupneohooke.cpp
    matelast_coupSaintVenantKirchhoff.cpp
    matelast_coupsimopister.cpp
    matelast_couptransverselyisotropic.cpp
    matelast_coupvarga.cpp
    matelast_iso1pow.cpp
    matelast_iso2pow.cpp
    matelast_isoanisoexpo.cpp
    matelast_isoexpopow.cpp
    matelast_isomooneyrivlin.cpp
    matelast_isomuscle_blemker.cpp
    matelast_isoneohooke.cpp
    matelast_isoogden.cpp
    matelast_isotestmaterial.cpp
    matelast_isovarga.cpp
    matelast_isovolaaagasser.cpp
    matelast_isovolHUdependentneohooke.cpp
    matelast_isoyeoh.cpp
    matelast_remodelfiber.cpp
    matelast_summand.cpp
    matelast_visco_coupmyocard.cpp
    matelast_visco_fract.cpp
    matelast_visco_generalizedgenmax.cpp
    matelast_visco_genmax.cpp
    matelast_visco_isoratedep.cpp
    matelast_vologden.cpp
    matelast_volpenalty.cpp
    matelast_volpow.cpp
    matelast_volsussmanbathe.cpp
    )

baci_add_library(matelast ${OBJS_MATELAST})

set(LIBS_MATELAST
    # cmake-format: sortable
    inpar
    mat
    )

baci_add_dependency(matelast ${LIBS_MATELAST})
