/*----------------------------------------------------------------------*/
/*! \file
 \brief a fluid material for porous multiphase flow defining one reaction (mass sources and sinks)

   \level 3

 *----------------------------------------------------------------------*/

#include "fluidporo_multiphase_singlereaction.H"

#include <vector>
#include "drt_globalproblem.H"
#include "standardtypes_cpp.H"
#include "matpar_bundle.H"

/*----------------------------------------------------------------------*
 *  constructor (public)                               vuong 08/16      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroSingleReaction::FluidPoroSingleReaction(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      numscal_(matdata->GetInt("NUMSCAL")),
      numvolfrac_(matdata->GetInt("NUMVOLFRAC")),
      totalnummultiphasedof_(matdata->GetInt("TOTALNUMDOF")),
      numfluidphases_(totalnummultiphasedof_ - 2 * numvolfrac_),
      scale_(matdata->Get<std::vector<int>>("SCALE")),
      coupling_(SetCouplingType(matdata)),
      functID_(matdata->GetInt("FUNCTID")),
      isinit_(false),
      scalarnames_(numscal_),
      pressurenames_(numfluidphases_),
      saturationnames_(numfluidphases_),
      porosityname_("porosity"),
      volfracnames_(numvolfrac_),
      volfracpressurenames_(numvolfrac_)

{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroSingleReaction::Initialize()
{
  switch (DRT::Problem::Instance()->NDim())
  {
    case 1:
      return InitializeInternal<1>();
    case 2:
      return InitializeInternal<2>();
    case 3:
      return InitializeInternal<3>();
    default:
      dserror("Unsupported dimension %d.", DRT::Problem::Instance()->NDim());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void MAT::PAR::FluidPoroSingleReaction::InitializeInternal()
{
  // safety check
  if (Function<dim>(functID_ - 1).NumberComponents() != 1)
    dserror("expected only one component for single phase reaction!");

  for (int k = 0; k < numscal_; k++)
  {
    // add scalar names
    {
      std::ostringstream temp;
      temp << k + 1;
      scalarnames_[k] = "phi" + temp.str();
    }
  }

  for (int k = 0; k < numfluidphases_; k++)
  {
    // add pressure names
    {
      std::ostringstream temp;
      temp << k + 1;
      pressurenames_[k] = "p" + temp.str();
    }

    // add saturation names
    {
      std::ostringstream temp;
      temp << k + 1;
      saturationnames_[k] = "S" + temp.str();
    }
  }


  // add additional volume fractions
  for (int k = 0; k < numvolfrac_; k++)
  {
    // add volume fraction names
    {
      std::ostringstream temp;
      temp << k + 1;
      volfracnames_[k] = "VF" + temp.str();
    }
    // add volume fraction pressure names
    {
      std::ostringstream temp;
      temp << k + 1;
      volfracpressurenames_[k] = "VFP" + temp.str();
    }
  }

  isinit_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroSingleReaction::EvaluateFunction(std::vector<double>& reacval,
    std::vector<std::vector<double>>& reacderivspressure,
    std::vector<std::vector<double>>& reacderivssaturation, std::vector<double>& reacderivsporosity,
    std::vector<std::vector<double>>& reacderivsvolfrac,
    std::vector<std::vector<double>>& reacderivsvolfracpressure,
    std::vector<std::vector<double>>& reacderivsscalar, const std::vector<double>& pressure,
    const std::vector<double>& saturation, const double& porosity,
    const std::vector<double>& volfracs, const std::vector<double>& volfracpressures,
    const std::vector<double>& scalar)
{
  switch (DRT::Problem::Instance()->NDim())
  {
    case 1:
      return EvaluateFunctionInternal<1>(reacval, reacderivspressure, reacderivssaturation,
          reacderivsporosity, reacderivsvolfrac, reacderivsvolfracpressure, reacderivsscalar,
          pressure, saturation, porosity, volfracs, volfracpressures, scalar);
    case 2:
      return EvaluateFunctionInternal<2>(reacval, reacderivspressure, reacderivssaturation,
          reacderivsporosity, reacderivsvolfrac, reacderivsvolfracpressure, reacderivsscalar,
          pressure, saturation, porosity, volfracs, volfracpressures, scalar);
    case 3:
      return EvaluateFunctionInternal<3>(reacval, reacderivspressure, reacderivssaturation,
          reacderivsporosity, reacderivsvolfrac, reacderivsvolfracpressure, reacderivsscalar,
          pressure, saturation, porosity, volfracs, volfracpressures, scalar);
    default:
      dserror("Unsupported dimension %d.", DRT::Problem::Instance()->NDim());
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void MAT::PAR::FluidPoroSingleReaction::EvaluateFunctionInternal(std::vector<double>& reacval,
    std::vector<std::vector<double>>& reacderivspressure,
    std::vector<std::vector<double>>& reacderivssaturation, std::vector<double>& reacderivsporosity,
    std::vector<std::vector<double>>& reacderivsvolfrac,
    std::vector<std::vector<double>>& reacderivsvolfracpressure,
    std::vector<std::vector<double>>& reacderivsscalar, const std::vector<double>& pressure,
    const std::vector<double>& saturation, const double& porosity,
    const std::vector<double>& volfracs, const std::vector<double>& volfracpressures,
    const std::vector<double>& scalar)
{
  // safety check if sizes fit
  CheckSizes(reacval, reacderivspressure, reacderivssaturation, reacderivsporosity,
      reacderivsvolfrac, reacderivsvolfracpressure, reacderivsscalar, pressure, saturation,
      porosity, volfracs, volfracpressures, scalar);

  std::vector<std::pair<std::string, double>> variables;
  variables.reserve(numfluidphases_ + numfluidphases_ + 1 + numscal_ + numvolfrac_ + numvolfrac_);

  std::vector<std::pair<std::string, double>> constants;
  constants.reserve(0);

  // set pressure values as variable
  for (int k = 0; k < numfluidphases_; k++)
    variables.push_back(std::pair<std::string, double>(pressurenames_[k], pressure[k]));

  // set saturation values as variable
  for (int k = 0; k < numfluidphases_; k++)
    variables.push_back(std::pair<std::string, double>(saturationnames_[k], saturation[k]));

  // set porosity value as variable
  variables.push_back(std::pair<std::string, double>(porosityname_, porosity));

  // set scalar values as variables
  for (int k = 0; k < numscal_; k++)
    variables.push_back(std::pair<std::string, double>(scalarnames_[k], scalar[k]));

  // set volfrac values as variables
  for (int k = 0; k < numvolfrac_; k++)
    variables.push_back(std::pair<std::string, double>(volfracnames_[k], volfracs[k]));

  // set volfrac pressure values as variables
  for (int k = 0; k < numvolfrac_; k++)
    variables.push_back(
        std::pair<std::string, double>(volfracpressurenames_[k], volfracpressures[k]));

  // evaluate the reaction term
  double curval = Function<dim>(functID_ - 1).Evaluate(0, variables, constants);
  // evaluate derivatives
  std::vector<double> curderivs(
      Function<dim>(functID_ - 1).EvaluateDerivative(0, variables, constants));

  // fill the output vector
  for (int k = 0; k < totalnummultiphasedof_; k++)
  {
    const int scale = (*scale_)[k];
    if (scale != 0)
    {
      // add the values from the reaction terms
      reacval[k] += scale * curval;

      // derivatives
      std::vector<double>& presk = reacderivspressure[k];
      std::vector<double>& satk = reacderivssaturation[k];
      for (int j = 0; j < numfluidphases_; j++)
      {
        presk[j] += scale * curderivs[j];
        satk[j] += scale * curderivs[numfluidphases_ + j];
      }
      reacderivsporosity[k] += scale * curderivs[2 * numfluidphases_];
      std::vector<double>& scalark = reacderivsscalar[k];
      for (int j = 0; j < numscal_; j++)
      {
        scalark[j] += scale * curderivs[2 * numfluidphases_ + 1 + j];
      }
      std::vector<double>& volfrack = reacderivsvolfrac[k];
      std::vector<double>& volfracpressurek = reacderivsvolfracpressure[k];
      for (int j = 0; j < numvolfrac_; j++)
      {
        volfrack[j] += scale * curderivs[2 * numfluidphases_ + 1 + numscal_ + j];
        volfracpressurek[j] +=
            scale * curderivs[2 * numfluidphases_ + 1 + numscal_ + numvolfrac_ + j];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *  check sizes of vectors                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroSingleReaction::CheckSizes(std::vector<double>& reacval,
    std::vector<std::vector<double>>& reacderivspressure,
    std::vector<std::vector<double>>& reacderivssaturation, std::vector<double>& reacderivsporosity,
    std::vector<std::vector<double>>& reacderivsvolfrac,
    std::vector<std::vector<double>>& reacderivsvolfracpressure,
    std::vector<std::vector<double>>& reacderivsscalar, const std::vector<double>& pressure,
    const std::vector<double>& saturation, const double& porosity,
    const std::vector<double>& volfracs, const std::vector<double>& volfracpressures,
    const std::vector<double>& scalar)
{
  //  dsassert(pressurenames_.size()==pressure.size(),"Invalid number of pressure values for this
  //  the fluid poro reaction material!");
  //  dsassert(pressurenames_.size()==saturation.size(),"Invalid number of pressure values for this
  //  the fluid poro reaction material!"); dsassert(scalarnames_.size()==scalar.size(),"Invalid
  //  number of pressure values for this the fluid poro reaction material!");

  if (numfluidphases_ != (int)pressure.size())
    dserror("Invalid number of pressure values for this fluid poro reaction material!");

  if (numfluidphases_ != (int)saturation.size())
    dserror("Invalid number of saturation values for this fluid poro reaction material!");
  if (numscal_ != (int)scalar.size())
    dserror("Invalid number of scalar values for this fluid poro reaction material!");
  if (numvolfrac_ != (int)volfracs.size())
    dserror("Invalid number of volfrac values for this fluid poro reaction material!");
  if (numvolfrac_ != (int)volfracpressures.size())
    dserror("Invalid number of volfrac values for this fluid poro reaction material!");

  if (totalnummultiphasedof_ != (int)reacderivsporosity.size())
    dserror(
        "Invalid length of vector for porosity derivatives for this fluid poro reaction material!");
  if (totalnummultiphasedof_ != (int)reacderivspressure.size())
    dserror(
        "Invalid length of vector for pressure derivatives for this fluid poro reaction material!");
  for (int k = 0; k < totalnummultiphasedof_; k++)
  {
    if (numfluidphases_ != (int)reacderivspressure[k].size())
      dserror(
          "Invalid length of vector for pressure derivatives for this fluid poro reaction "
          "material!");
  }
  if (totalnummultiphasedof_ != (int)reacderivssaturation.size())
    dserror(
        "Invalid length of vector for pressure derivatives for this fluid poro reaction material!");
  for (int k = 0; k < totalnummultiphasedof_; k++)
  {
    if (numfluidphases_ != (int)reacderivssaturation[k].size())
      dserror(
          "Invalid length of vector for pressure derivatives for this fluid poro reaction "
          "material!");
  }
  if (totalnummultiphasedof_ != (int)reacderivsscalar.size())
    dserror(
        "Invalid length of vector for scalar derivatives for this fluid poro reaction material!");
  for (int k = 0; k < totalnummultiphasedof_; k++)
  {
    if (numscal_ != (int)reacderivsscalar[k].size())
      dserror(
          "Invalid length of vector for scalar derivatives for this fluid poro reaction material!");
  }
  if (totalnummultiphasedof_ != (int)reacderivsvolfrac.size())
    dserror(
        "Invalid length of vector for vol frac derivatives for this fluid poro reaction material!");
  for (int k = 0; k < totalnummultiphasedof_; k++)
  {
    if (numvolfrac_ != (int)reacderivsvolfrac[k].size())
      dserror(
          "Invalid length of vector for scalar derivatives for this fluid poro reaction material!");
  }
  if (totalnummultiphasedof_ != (int)reacderivsvolfracpressure.size())
    dserror(
        "Invalid length of vector for vol frac derivatives for this fluid poro reaction material!");
  for (int k = 0; k < totalnummultiphasedof_; k++)
  {
    if (numvolfrac_ != (int)reacderivsvolfracpressure[k].size())
      dserror(
          "Invalid length of vector for scalar derivatives for this fluid poro reaction material!");
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
inline DRT::UTILS::VariableExprFunction<dim>& MAT::PAR::FluidPoroSingleReaction::Function(
    int functnum) const
{
  try
  {
    DRT::UTILS::VariableExprFunction<dim>& funct =
        dynamic_cast<DRT::UTILS::VariableExprFunction<dim>&>(
            DRT::Problem::Instance()->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(functnum));
    return funct;
  }
  catch (std::bad_cast& exp)
  {
    dserror(
        "Cast to VarExp Function failed! For phase law definition only 'VARFUNCTION' functions are "
        "allowed!\n"
        "Check your input file!");
    return dynamic_cast<DRT::UTILS::VariableExprFunction<dim>&>(
        DRT::Problem::Instance()->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(functnum));
  }
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                             vuong 08/16      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::FluidPoroSingleReaction::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FluidPoroSingleReaction(this));
}

/*----------------------------------------------------------------------*
 *  translate coupling type                             vuong 08/16      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroSingleReaction::porofluid_reaction_coupling
MAT::PAR::FluidPoroSingleReaction::SetCouplingType(Teuchos::RCP<MAT::PAR::Material> matdata)
{
  if (*(matdata->Get<std::string>("COUPLING")) == "scalar_by_function")
  {
    return porofluid_reac_coup_scalarsbyfunction;
  }
  else if (*(matdata->Get<std::string>("COUPLING")) == "no_coupling")
  {
    return porofluid_reac_coup_none;
  }
  else
  {
    return porofluid_reac_coup_none;
  }
}

/*----------------------------------------------------------------------*
  global instance of parameter class                   vuong 08/16     |
*----------------------------------------------------------------------*/
MAT::FluidPoroSingleReactionType MAT::FluidPoroSingleReactionType::instance_;

/*----------------------------------------------------------------------*
 *  Create material from given data                          vuong 08/16 |
 *----------------------------------------------------------------------*/

DRT::ParObject* MAT::FluidPoroSingleReactionType::Create(const std::vector<char>& data)
{
  MAT::FluidPoroSingleReaction* fluid_poro = new MAT::FluidPoroSingleReaction();
  fluid_poro->Unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroSingleReaction::FluidPoroSingleReaction() : params_(NULL) {}

/*----------------------------------------------------------------------*
 *   Create material with parameters                         vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroSingleReaction::FluidPoroSingleReaction(MAT::PAR::FluidPoroSingleReaction* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                           vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSingleReaction::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
 * unpack material                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSingleReaction::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::FluidPoroSingleReaction*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *  initialize                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSingleReaction::Initialize()
{
  params_->Initialize();
  return;
}

/*----------------------------------------------------------------------*
 *  set values in function                                 vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSingleReaction::EvaluateReaction(std::vector<double>& reacval,
    std::vector<std::vector<double>>& reacderivspressure,
    std::vector<std::vector<double>>& reacderivssaturation, std::vector<double>& reacderivsporosity,
    std::vector<std::vector<double>>& reacderivsvolfrac,
    std::vector<std::vector<double>>& reacderivsvolfracpressure,
    std::vector<std::vector<double>>& reacderivsscalar, const std::vector<double>& pressure,
    const std::vector<double>& saturation, const double& porosity,
    const std::vector<double>& volfracs, const std::vector<double>& volfracpressures,
    const std::vector<double>& scalar)
{
  params_->EvaluateFunction(reacval, reacderivspressure, reacderivssaturation, reacderivsporosity,
      reacderivsvolfrac, reacderivsvolfracpressure, reacderivsscalar, pressure, saturation,
      porosity, volfracs, volfracpressures, scalar);

  return;
}
