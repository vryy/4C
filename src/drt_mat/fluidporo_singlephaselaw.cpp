/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_singlephaselaw.cpp

 \brief a material defining pressure-saturation relationship
        for a single phase within a multiphase porous fluid

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/


#include "fluidporo_singlephaselaw.H"

#include "matpar_bundle.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseLaw::FluidPoroPhaseLaw(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  numdof_(matdata->GetInt("NUMDOF")),
  presids_(matdata->Get<std::vector<int> >("PRESCOEFF"))
{
  // check if sizes fit
  if (numdof_ != (int)presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", numdof_, presids_->size());
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseLawLinear::FluidPoroPhaseLawLinear(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: FluidPoroPhaseLaw(matdata),
  reltensions_(matdata->GetDouble("RELTENSION")),
  sat0_(matdata->GetDouble("SATURATION_0"))
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawLinear::EvaluateSaturation(
    const std::vector<double>& pressure) const
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(), presids_->size());

  double presval = std::inner_product(presids_->begin(),presids_->end(),pressure.begin(),0.0);

  double saturation = sat0_ + reltensions_*presval;

  return saturation;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawLinear::EvaluateDerivOfSaturationWrtPressure(
    int doftoderive,
    const std::vector<double>& state) const
{
  // check if sizes fit
  if (state.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", state.size(), presids_->size());

  if((*presids_)[doftoderive]==0)
    return 0.0;

  double deriv = reltensions_;

  return deriv*(*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawLinear::EvaluateDerivOfPressureWrtSaturation(
    int doftoderive,
    double saturation) const
{

  if((*presids_)[doftoderive]==0)
    return 0.0;

  double deriv = 1.0/reltensions_;

  return deriv*(*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawLinear::EvaluateGenPressure(double saturation) const
{

  double presval = 1.0/reltensions_ *(saturation-sat0_);

  return presval;
}

/************************************************************************/
/************************************************************************/

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseLawTangent::FluidPoroPhaseLawTangent(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: FluidPoroPhaseLaw(matdata),
  reltensions_(matdata->GetDouble("RELTENSION")),
  exp_(matdata->GetDouble("EXP")),
  sat0_(matdata->GetDouble("SATURATION_0"))
{

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawTangent::EvaluateSaturation(
    const std::vector<double>& pressure) const
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(), presids_->size());

  double presval = std::inner_product(presids_->begin(),presids_->end(),pressure.begin(),0.0);

  double saturation = sat0_ - std::pow( 2/M_PI * std::atan(reltensions_*presval),exp_);

  return saturation;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawTangent::EvaluateDerivOfSaturationWrtPressure(
    int doftoderive,
    const std::vector<double>& state) const
{
  // check if sizes fit
  if (state.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", state.size(), presids_->size());

  if((*presids_)[doftoderive]==0)
    return 0.0;

  double presval = std::inner_product(presids_->begin(),presids_->end(),state.begin(),0.0);

  double deriv = - exp_*std::pow( 2/M_PI * std::atan(reltensions_*presval),exp_-1.0)
                   *2.0*reltensions_/(M_PI * (1.0+(reltensions_*presval)*(reltensions_*presval)));

  return deriv*(*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawTangent::EvaluateDerivOfPressureWrtSaturation(
    int doftoderive,
    double saturation) const
{

  if((*presids_)[doftoderive]==0)
    return 0.0;

  double deriv = -0.5 * M_PI/(reltensions_*exp_)*std::pow(sat0_-saturation,1.0/exp_-1.0)*
                  (1.0+std::pow(std::tan(0.5*M_PI*std::pow(sat0_-saturation,1.0/exp_)),2));

  return deriv*(*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawTangent::EvaluateGenPressure(double saturation) const
{

  double presval = 1.0/reltensions_ *std::tan(0.5*M_PI*std::pow(sat0_-saturation,1.0/exp_));

  return presval;
}
