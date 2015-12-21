/*----------------------------------------------------------------------*/
/*!
\file elast_volgrowthpenalty.cpp
\brief


the input line should read
  MAT 1 ELAST_VolGrowthPenalty EPSILON 1. EXPONENT 2.

<pre>
Maintainer: Fabian Bräu
            braeu@lnm.mw.tum.de
            089/289 15236
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_volgrowthpenalty.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::VolGrowthPenalty::VolGrowthPenalty(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  eps_(matdata->GetDouble("EPSILON")),
  n_(matdata->GetDouble("EXPONENT"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VolGrowthPenalty::VolGrowthPenalty(MAT::ELASTIC::PAR::VolGrowthPenalty* params)
  : params_(params)
{
  // Initialize some variables
  eps_ = params_->eps_;
  n_ = params_->n_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolGrowthPenalty::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data,eps_);
  AddtoPack(data,n_);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolGrowthPenalty::UnpackSummand(const std::vector<char>& data,std::vector<char>::size_type& position)
{
  ExtractfromPack(position,data,eps_);
  ExtractfromPack(position,data,n_);

  return;
}


/*----------------------------------------------------------------------
 *                                                       braeu 10/2015  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolGrowthPenalty::AddDerivativesModified(
    LINALG::Matrix<3,1>& dPmodI,
    LINALG::Matrix<6,1>& ddPmodII,
    const LINALG::Matrix<3,1>& modinv,
    const int eleGID
)
{
  dserror("Not implemented so far!");

  return;
}
