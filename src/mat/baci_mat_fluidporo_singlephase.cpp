/*----------------------------------------------------------------------*/
/*! \file
 \brief a single phase material of a multiphase porous fluid

   \level 3

 *----------------------------------------------------------------------*/



#include "baci_mat_fluidporo_singlephase.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_fluidporo_singlephaseDof.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_mat_poro_density_law.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *  constructor (public)                               vuong 08/16      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroSinglePhase::FluidPoroSinglePhase(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), density_(*matdata->Get<double>("DENSITY")), isinit_(false)
{
  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // create density law
  densitylaw_ = MAT::PAR::PoroDensityLaw::CreateDensityLaw(*matdata->Get<int>("DENSITYLAWID"));

  // create permeability law
  relpermeabilitylaw_ = MAT::PAR::FluidPoroRelPermeabilityLaw::CreateRelPermeabilityLaw(
      *matdata->Get<int>("RELPERMEABILITYLAWID"));

  // create viscosity law
  viscositylaw_ =
      MAT::PAR::FluidPoroViscosityLaw::CreateViscosityLaw(*matdata->Get<int>("VISCOSITYLAWID"));

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      GLOBAL::Problem::Instance(probinst)->Materials()->ById(*matdata->Get<int>("DOFTYPEID"));

  switch (curmat->Type())
  {
    case INPAR::MAT::m_fluidporo_phasedof_diffpressure:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofDiffPressure(curmat));
      phasedof_ = static_cast<MAT::PAR::FluidPoroPhaseDofDiffPressure*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_fluidporo_phasedof_pressure:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofPressure(curmat));
      phasedof_ = static_cast<MAT::PAR::FluidPoroPhaseDofPressure*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_fluidporo_phasedof_saturation:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofSaturation(curmat));
      phasedof_ = static_cast<MAT::PAR::FluidPoroPhaseDofSaturation*>(curmat->Parameter());
      break;
    }
    default:
      FOUR_C_THROW("invalid pressure-saturation law for material %d", curmat->Type());
      break;
  }
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                             vuong 08/16      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::FluidPoroSinglePhase::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FluidPoroSinglePhase(this));
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                             vuong 08/16      |
 *----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroSinglePhase::Initialize()
{
  if (not isinit_)
  {
    phasedof_->Initialize();
    isinit_ = true;
  }
  return;
}

/*----------------------------------------------------------------------*
  global instance of parameter class                   vuong 08/16     |
*----------------------------------------------------------------------*/
MAT::FluidPoroSinglePhaseType MAT::FluidPoroSinglePhaseType::instance_;

/*----------------------------------------------------------------------*
 *  Create material from given data                          vuong 08/16 |
 *----------------------------------------------------------------------*/

CORE::COMM::ParObject* MAT::FluidPoroSinglePhaseType::Create(const std::vector<char>& data)
{
  MAT::FluidPoroSinglePhase* fluid_poro = new MAT::FluidPoroSinglePhase();
  fluid_poro->Unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroSinglePhase::FluidPoroSinglePhase() : params_(nullptr) {}

/*----------------------------------------------------------------------*
 *   Create material with parameters                         vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroSinglePhase::FluidPoroSinglePhase(MAT::PAR::FluidPoroSinglePhase* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                           vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSinglePhase::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
 * unpack material                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSinglePhase::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::FluidPoroSinglePhase*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *  initialize                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSinglePhase::Initialize()
{
  params_->Initialize();
  return;
}

/*----------------------------------------------------------------------*
 *  return dof type                                         vuong 08/16 |
 *----------------------------------------------------------------------*/
INPAR::MAT::MaterialType MAT::FluidPoroSinglePhase::PoroDofType() const
{
  return params_->phasedof_->Type();
}

/*----------------------------------------------------------------------*
 *  return phase law type                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
INPAR::MAT::MaterialType MAT::FluidPoroSinglePhase::PoroPhaseLawType() const
{
  return params_->phasedof_->PoroPhaseLawType();
}

/*----------------------------------------------------------------------*
 *  fill the dof matrix with the phase dofs                 vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSinglePhase::FillDoFMatrix(
    CORE::LINALG::SerialDenseMatrix& dofmat, int numphase) const
{
  params_->phasedof_->FillDoFMatrix(dofmat, numphase);
  return;
}

/*----------------------------------------------------------------------*
 *  Evaluate generalized pressure of a phase                vuong 08/16 |
 *----------------------------------------------------------------------*/
double MAT::FluidPoroSinglePhase::EvaluateGenPressure(
    int phasenum, const std::vector<double>& state) const
{
  return params_->phasedof_->EvaluateGenPressure(phasenum, state);
}

/*----------------------------------------------------------------------*
 *   Evaluate saturation of the phase                       vuong 08/16 |
 *----------------------------------------------------------------------*/
double MAT::FluidPoroSinglePhase::EvaluateSaturation(
    int phasenum, const std::vector<double>& state, const std::vector<double>& pressure) const
{
  return params_->phasedof_->EvaluateSaturation(phasenum, state, pressure);
}

/*--------------------------------------------------------------------------*
 *  Evaluate derivative of saturation w.r.t. pressure           vuong 08/16 |
 *---------------------------------------------------------------------------*/
double MAT::FluidPoroSinglePhase::EvaluateDerivOfSaturationWrtPressure(
    int phasenum, int doftoderive, const std::vector<double>& pressure) const
{
  return params_->phasedof_->EvaluateDerivOfSaturationWrtPressure(phasenum, doftoderive, pressure);
}

/*--------------------------------------------------------------------------*
 *  Evaluate 2nd derivative of saturation w.r.t. pressure  kremheller 05/17 |
 *---------------------------------------------------------------------------*/
double MAT::FluidPoroSinglePhase::EvaluateSecondDerivOfSaturationWrtPressure(int phasenum,
    int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure) const
{
  return params_->phasedof_->EvaluateSecondDerivOfSaturationWrtPressure(
      phasenum, firstdoftoderive, seconddoftoderive, pressure);
}

/*----------------------------------------------------------------------------------------*
 * Evaluate derivative of degree of freedom with respect to pressure          vuong 08/16 |
 *----------------------------------------------------------------------------------------*/
double MAT::FluidPoroSinglePhase::EvaluateDerivOfDofWrtPressure(
    int phasenum, int doftoderive, const std::vector<double>& state) const
{
  return params_->phasedof_->EvaluateDerivOfDofWrtPressure(phasenum, doftoderive, state);
}

/*----------------------------------------------------------------------*
 *  constructor (public)                               kremheller 10/17 |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroSingleVolFrac::FluidPoroSingleVolFrac(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      density_(*matdata->Get<double>("DENSITY")),
      diffusivity_(*matdata->Get<double>("DIFFUSIVITY")),
      scalardependentflux_(*matdata->Get<bool>("AddScalarDependentFlux")),
      numscal_(*matdata->Get<int>("NUMSCAL")),
      scalardiffs_(*(matdata->Get<std::vector<double>>("SCALARDIFFS"))),
      omega_half_(*(matdata->Get<std::vector<double>>("OMEGA_HALF"))),
      isinit_(false)
{
  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // safety checks
  if (numscal_ < 0) FOUR_C_THROW("NUMSCAL smaller than zero not possible");

  if (scalardependentflux_ && numscal_ == 0)
    FOUR_C_THROW("AddScalarDependentFlux has been set to YES, but NUMSCAL is equal to zero");

  if (scalardependentflux_ && scalardiffs_.size() == 0)
    FOUR_C_THROW(
        "AddScalarDependentFlux has been set to YES, but length of SCALARDIFFS is equal to zero");

  if (scalardependentflux_ && omega_half_.size() == 0)
    FOUR_C_THROW(
        "AddScalarDependentFlux has been set to YES, but length of OMEGA_HALF is equal to zero");

  if (!scalardependentflux_ && numscal_ > 0)
    FOUR_C_THROW("AddScalarDependentFlux has been set to NO, but NUMSCAL is greater than zero");

  if (!scalardependentflux_ && scalardiffs_.size() > 0)
    FOUR_C_THROW(
        "AddScalarDependentFlux has been set to NO, but length of SCALARDIFFS is greater than "
        "zero");

  if (!scalardependentflux_ && omega_half_.size() > 0)
    FOUR_C_THROW(
        "AddScalarDependentFlux has been set to NO, but length of OMEGA_HALF is greater than zero");
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                           kremheller 10/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::FluidPoroSingleVolFrac::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FluidPoroSingleVolFrac(this));
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                           kremheller 10/17 |
 *----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroSingleVolFrac::Initialize()
{
  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------*
  global instance of parameter class                   kremheller 10/17 |
*----------------------------------------------------------------------*/
MAT::FluidPoroSingleVolFracType MAT::FluidPoroSingleVolFracType::instance_;

/*----------------------------------------------------------------------*
 *  Create material from given data                    kremheller 10/17 |
 *----------------------------------------------------------------------*/

CORE::COMM::ParObject* MAT::FluidPoroSingleVolFracType::Create(const std::vector<char>& data)
{
  MAT::FluidPoroSingleVolFrac* fluid_poro = new MAT::FluidPoroSingleVolFrac();
  fluid_poro->Unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                             kremheller 10/17 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroSingleVolFrac::FluidPoroSingleVolFrac() : params_(nullptr) {}

/*----------------------------------------------------------------------*
 *   Create material with parameters                    kremheller 10/17 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroSingleVolFrac::FluidPoroSingleVolFrac(MAT::PAR::FluidPoroSingleVolFrac* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                      kremheller 10/17 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSingleVolFrac::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
 * unpack material                                      kremheller 10/17 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSingleVolFrac::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::FluidPoroSingleVolFrac*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *  initialize                                         kremheller 10/17 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroSingleVolFrac::Initialize()
{
  params_->Initialize();
  return;
}

/*----------------------------------------------------------------------*
 *  constructor (public)                               kremheller 02/18 |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroVolFracPressure::FluidPoroVolFracPressure(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      permeability_(*matdata->Get<double>("PERMEABILITY")),
      min_volfrac_(*matdata->Get<double>("MIN_VOLFRAC")),
      isinit_(false)
{
  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // create viscosity law
  viscositylaw_ =
      MAT::PAR::FluidPoroViscosityLaw::CreateViscosityLaw(*matdata->Get<int>("VISCOSITYLAWID"));
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                           kremheller 02/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::FluidPoroVolFracPressure::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FluidPoroVolFracPressure(this));
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                           kremheller 02/18 |
 *----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroVolFracPressure::Initialize()
{
  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------*
  global instance of parameter class                   kremheller 02/18 |
*----------------------------------------------------------------------*/
MAT::FluidPoroVolFracPressureType MAT::FluidPoroVolFracPressureType::instance_;

/*----------------------------------------------------------------------*
 *  Create material from given data                    kremheller 02/18 |
 *----------------------------------------------------------------------*/

CORE::COMM::ParObject* MAT::FluidPoroVolFracPressureType::Create(const std::vector<char>& data)
{
  MAT::FluidPoroVolFracPressure* fluid_poro = new MAT::FluidPoroVolFracPressure();
  fluid_poro->Unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                             kremheller 02/18 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroVolFracPressure::FluidPoroVolFracPressure() : params_(nullptr) {}

/*----------------------------------------------------------------------*
 *   Create material with parameters                    kremheller 02/18 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroVolFracPressure::FluidPoroVolFracPressure(MAT::PAR::FluidPoroVolFracPressure* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                      kremheller 02/18 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroVolFracPressure::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
 * unpack material                                      kremheller 02/18 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroVolFracPressure::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::FluidPoroVolFracPressure*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *  initialize                                         kremheller 02/18 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroVolFracPressure::Initialize()
{
  params_->Initialize();
  return;
}

FOUR_C_NAMESPACE_CLOSE
