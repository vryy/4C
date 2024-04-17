/*----------------------------------------------------------------------*/
/*! \file
\brief auxiliary material for macro-scale elements in multi-scale simulations of scalar transport
problems. This material handles the communication between micro and macro materials

\level 2

*/
/*----------------------------------------------------------------------*/
#include "baci_mat_scatra_micro_macro_coupling.hpp"

#include "baci_mat_par_material.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::PAR::ScatraMicroMacroCoupling::ScatraMicroMacroCoupling(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : microfile_(*(matdata->Get<std::string>("MICROFILE"))),
      microdisnum_(*matdata->Get<int>("MICRODIS_NUM")),
      A_s_(*matdata->Get<double>("A_s"))
{
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::ScatraMicroMacroCoupling::ScatraMicroMacroCoupling() : matgp_() {}

FOUR_C_NAMESPACE_CLOSE
