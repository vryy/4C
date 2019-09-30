/*----------------------------------------------------------------------*/
/*! \file
\brief  fluid material for poroelasticity problems

\maintainer Johannes Kremheller

\level 2
 *-----------------------------------------------------------------------*/

#include <vector>
#include "fluidporo.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*
 *  constructor (public)                               vuong 06/11      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoro::FluidPoro(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(matdata->GetDouble("DYNVISCOSITY")),
      density_(matdata->GetDouble("DENSITY")),
      permeability_(matdata->GetDouble("PERMEABILITY")),
      type_(undefined),
      varyingpermeability_(false),
      permeabilityfunc_(MAT::PAR::pf_undefined),
      permeabilitycorrectionfactor_(1.0),
      initialporosity_(1.0)
{
  const std::string* typestring = matdata->Get<std::string>("TYPE");

  if (*typestring == "Darcy")
    type_ = darcy;
  else if (*typestring == "Darcy-Brinkman")
    type_ = darcy_brinkman;

  const std::string* pfuncstring = matdata->Get<std::string>("PERMEABILITYFUNCTION");

  if (*pfuncstring == "Const")
    permeabilityfunc_ = MAT::PAR::const_;
  else if (*pfuncstring == "Kozeny_Carman")
    permeabilityfunc_ = MAT::PAR::kozeny_karman;
  else
    dserror("Unknown Permeabilityfunction: %s", pfuncstring->c_str());
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                             vuong 06/11      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::FluidPoro::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FluidPoro(this));
}

/*----------------------------------------------------------------------*
  Set Initial Porosity (public)                          vuong 06/11     |
*----------------------------------------------------------------------*/
void MAT::PAR::FluidPoro::SetInitialPorosity(double initialporosity)
{
  initialporosity_ = initialporosity;

  if (permeabilityfunc_ == MAT::PAR::const_)
  {
    permeabilitycorrectionfactor_ = 1.0;
  }
  else if (permeabilityfunc_ == MAT::PAR::kozeny_karman)
  {
    // c = (phi0^3 / (1 - phi0^2))
    permeabilitycorrectionfactor_ = initialporosity_ * initialporosity_ * initialporosity_ /
                                    (1 - initialporosity_ * initialporosity_);
  }
  return;
}

/*----------------------------------------------------------------------*
                                                          vuong 06/11     |
*----------------------------------------------------------------------*/
MAT::FluidPoroType MAT::FluidPoroType::instance_;

/*----------------------------------------------------------------------*
 *                                                           vuong 06/11 |
 *----------------------------------------------------------------------*/

DRT::ParObject* MAT::FluidPoroType::Create(const std::vector<char>& data)
{
  MAT::FluidPoro* fluid_poro = new MAT::FluidPoro();
  fluid_poro->Unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *                                                          vuong 06/11 |
 *----------------------------------------------------------------------*/
MAT::FluidPoro::FluidPoro() : params_(NULL) {}

/*----------------------------------------------------------------------*
                                                         vuong 06/11 |
*----------------------------------------------------------------------*/
MAT::FluidPoro::FluidPoro(MAT::PAR::FluidPoro* params) : params_(params) {}

/*----------------------------------------------------------------------*
 *                                                          vuong 06/11 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoro::Pack(DRT::PackBuffer& data) const
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
 *                                                          vuong 06/11 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoro::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::FluidPoro*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *                                                          vuong 06/11 |
 *----------------------------------------------------------------------*/
double MAT::FluidPoro::ComputeReactionCoeff() const
{
  // check for zero or negative viscosity
  if (Viscosity() < EPS15) dserror("zero or negative viscosity");

  // check for zero or negative permeability
  if (Permeability() < EPS15) dserror("zero or negative permeability");

  // viscosity divided by permeability
  double reacoeff = Viscosity() / Permeability();

  return reacoeff;
}

/*----------------------------------------------------------------------*
 *                                                          vuong 06/11 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoro::ComputeReactionTensor(
    LINALG::Matrix<2, 2>& reactiontensor, double J, double porosity) const
{
  // viscosity divided by permeability
  double reacoeff = ComputeReactionCoeff();  // ReactionCoeff for const PermeabilityFunction

  reactiontensor.Clear();

  if (PermeabilityFunction() == MAT::PAR::kozeny_karman)
  {
    reacoeff *= (1 - porosity * porosity * J * J) / (porosity * porosity * porosity * J * J * J) *
                PermeabilityCorrectionFactor();
  }

  for (int i = 0; i < 2; i++) reactiontensor(i, i) = reacoeff;
  return;
}

/*----------------------------------------------------------------------*
 *                                                          vuong 06/11 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoro::ComputeReactionTensor(
    LINALG::Matrix<3, 3>& reactiontensor, double J, double porosity) const
{
  // viscosity divided by permeability
  double reacoeff = ComputeReactionCoeff();  // ReactionCoeff for const PermeabilityFunction

  reactiontensor.Clear();

  if (PermeabilityFunction() == MAT::PAR::kozeny_karman)
  {
    reacoeff *= (1 - porosity * porosity * J * J) / (porosity * porosity * porosity * J * J * J) *
                PermeabilityCorrectionFactor();
  }

  for (int i = 0; i < 3; i++) reactiontensor(i, i) = reacoeff;
  return;
}

/*----------------------------------------------------------------------*
 *                                                          vuong 06/11 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoro::ComputeLinMatReactionTensor(LINALG::Matrix<2, 2>& linreac_dphi,
    LINALG::Matrix<2, 2>& linreac_dJ, double J, double porosity) const
{
  // viscosity divided by permeability
  double reacoeff = ComputeReactionCoeff();

  linreac_dphi.Clear();
  linreac_dJ.Clear();

  if (PermeabilityFunction() == MAT::PAR::const_)
    return;  // Permeability is not a function of porosity or J
  else if (PermeabilityFunction() == MAT::PAR::kozeny_karman)
  {
    // d(isotropic_mat_reactiontensor)/d(phi) = reacoeff * [(J * phi)^2 - 3] / ( J^3 * phi^4 )
    // d(isotropic_mat_reactiontensor)/d(J) = reacoeff * [(J * phi)^2 - 3] / ( J^4 * phi^3 )

    double linreac_tmp = reacoeff * ((J * J * porosity * porosity) - 3.0) /
                         (J * J * J * porosity * porosity * porosity) *
                         PermeabilityCorrectionFactor();
    linreac_dphi(0, 0) = linreac_tmp / porosity;
    linreac_dJ(0, 0) = linreac_tmp / J;

    for (int i = 1; i < 2; i++)
    {
      linreac_dphi(i, i) = linreac_dphi(0, 0);
      linreac_dJ(i, i) = linreac_dJ(0, 0);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 *                                                          vuong 06/11 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoro::ComputeLinMatReactionTensor(LINALG::Matrix<3, 3>& linreac_dphi,
    LINALG::Matrix<3, 3>& linreac_dJ, double J, double porosity) const
{
  // viscosity divided by permeability
  double reacoeff = ComputeReactionCoeff();

  linreac_dphi.Clear();
  linreac_dJ.Clear();

  if (PermeabilityFunction() == MAT::PAR::const_)
    return;  // Permeability is not a function of porosity or J
  else if (PermeabilityFunction() == MAT::PAR::kozeny_karman)
  {
    // d(isotropic_mat_reactiontensor)/d(phi) = reacoeff * [(J * phi)^2 - 3] / ( J^3 * phi^4 )
    // d(isotropic_mat_reactiontensor)/d(J) = reacoeff * [(J * phi)^2 - 3] / ( J^4 * phi^3 )

    double linreac_tmp = reacoeff * ((J * J * porosity * porosity) - 3.0) /
                         (J * J * J * porosity * porosity * porosity) *
                         PermeabilityCorrectionFactor();
    linreac_dphi(0, 0) = linreac_tmp / porosity;
    linreac_dJ(0, 0) = linreac_tmp / J;

    for (int i = 1; i < 3; i++)
    {
      linreac_dphi(i, i) = linreac_dphi(0, 0);
      linreac_dJ(i, i) = linreac_dJ(0, 0);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 *                                                           vuong 06/11 |
 *----------------------------------------------------------------------*/
double MAT::FluidPoro::EffectiveViscosity() const
{
  // set zero viscosity and only modify it for Darcy-Stokes problems
  double viscosity = -1.0;
  if (Type() == PAR::darcy)
    viscosity = 0.0;
  else if (Type() == PAR::darcy_brinkman)
    viscosity = Viscosity();
  else
    dserror("Unknown flow type for porous flow");

  return viscosity;
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)     vuong 05/12|
 *----------------------------------------------------------------------*/
void MAT::FluidPoro::EvaluateViscStress(LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
    const LINALG::Matrix<6, 1>* glstrain, const int gp, Teuchos::ParameterList& params)

{
  dserror("macroscopic viscous stress not yet implemented for poroelasticity");
  return;
}
