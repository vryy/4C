/*----------------------------------------------------------------------*/
/*! \file
\brief auxiliary material for macro-scale elements in multi-scale simulations of scalar transport
problems

\level 2

*/
/*----------------------------------------------------------------------*/
#include "baci_mat_scatra_multiscale.hpp"

#include "baci_mat_par_material.hpp"

BACI_NAMESPACE_OPEN

/*--------------------------------------------------------------------*
 | constructor                                             fang 07/17 |
 *--------------------------------------------------------------------*/
MAT::PAR::ScatraMultiScale::ScatraMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata)
    : microfile_(*(matdata->Get<std::string>("MICROFILE"))),
      microdisnum_(matdata->GetInt("MICRODIS_NUM")),
      A_s_(matdata->GetDouble("A_s"))
{
  return;
}


/*--------------------------------------------------------------------*
 | construct empty material                                fang 07/17 |
 *--------------------------------------------------------------------*/
MAT::ScatraMultiScale::ScatraMultiScale() : matgp_() { return; }

BACI_NAMESPACE_CLOSE
