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

#include "../drt_lib/drt_globalproblem.H"

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

/************************************************************************/
/************************************************************************/

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseLawByFunction::FluidPoroPhaseLawByFunction(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: FluidPoroPhaseLaw(matdata),
  functionID_saturation_(matdata->GetInt("FUNCTSAT")),
  functionID_pressure_(matdata->GetInt("FUNCTPRES"))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroPhaseLawByFunction::Initialize()
{
  if(Function(functionID_saturation_-1).NumberComponents()!=1)
    dserror("expected only one component for the saturation evaluation");
  if(Function(functionID_pressure_-1).NumberComponents()!=1)
    dserror("expected only one component for the pressure evaluation");

  if(not Function(functionID_pressure_-1).IsVariable(0,"S"))
    Function(functionID_pressure_-1).AddVariable(0,"S",0.0);
  if(not Function(functionID_saturation_-1).IsVariable(0,"dp"))
    Function(functionID_saturation_-1).AddVariable(0,"dp",0.0);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline DRT::UTILS::VariableExprFunction& MAT::PAR::FluidPoroPhaseLawByFunction::Function(int functnum) const
{
  try
  {
    DRT::UTILS::VariableExprFunction& funct =
        dynamic_cast<DRT::UTILS::VariableExprFunction&>(DRT::Problem::Instance()->Funct(functnum));

    return funct;
  }
  catch(std::bad_cast & exp)
  {
    dserror("Cast to VarExp Function failed! For phase law definition only 'VAREXPR' functions are allowed!\n"
        "Check your input file!");
    return dynamic_cast<DRT::UTILS::VariableExprFunction&>(DRT::Problem::Instance()->Funct(functnum));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateSaturation(
    const std::vector<double>& pressure) const
{
  // check if sizes fit
  if (pressure.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", pressure.size(), presids_->size());

  double presval = std::inner_product(presids_->begin(),presids_->end(),pressure.begin(),0.0);

  std::vector<std::pair<std::string,double> > variables(1);
  variables[0]=std::pair<std::string,double>("dp",presval);

  return Function(functionID_saturation_-1).Evaluate(0,variables);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateDerivOfSaturationWrtPressure(
    int doftoderive,
    const std::vector<double>& state) const
{
  // check if sizes fit
  if (state.size() != presids_->size())
    dserror("number of dofs %d does not fit to size of dof vector %d", state.size(), presids_->size());

  if((*presids_)[doftoderive]==0)
    return 0.0;

  double presval = std::inner_product(presids_->begin(),presids_->end(),state.begin(),0.0);
  std::vector<std::pair<std::string,double> > variables(1);
  variables[0]=std::pair<std::string,double>("dp",presval);

  std::vector<std::vector<double> > deriv = Function(functionID_saturation_-1).FctDer(0,variables);

  return deriv[0][0]*(*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateDerivOfPressureWrtSaturation(
    int doftoderive,
    double saturation) const
{

  if((*presids_)[doftoderive]==0)
    return 0.0;

  std::vector<std::pair<std::string,double> > variables(1);
  variables[0]=std::pair<std::string,double>("S",saturation);

  std::vector<std::vector<double> > deriv = Function(functionID_pressure_-1).FctDer(0,variables);

  return deriv[0][0]*(*presids_)[doftoderive];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseLawByFunction::EvaluateGenPressure(double saturation) const
{
  std::vector<std::pair<std::string,double> > variables(1);
  variables[0]=std::pair<std::string,double>("S",saturation);

  return Function(functionID_pressure_-1).Evaluate(0,variables);
}
