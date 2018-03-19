/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_singlephase.cpp

 \brief a single phase material of a multiphase porous fluid

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/



#include <vector>
#include "fluidporo_singlephase.H"

#include "fluidporo_singlephaseDof.H"
#include "poro_density_law.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*
 *  constructor (public)                               vuong 08/16      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroSinglePhase::FluidPoroSinglePhase(Teuchos::RCP<MAT::PAR::Material> matdata) :
  Parameter(matdata),
  density_(matdata->GetDouble("DENSITY")),
  isinit_(false)
{
  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("Sorry dude, cannot work out problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("Sorry dude, no materials defined.");

  // create density law
  densitylaw_ = MAT::PAR::PoroDensityLaw::CreateDensityLaw(matdata->GetInt("DENSITYLAWID"));

  // create permeability law
  relpermeabilitylaw_ = MAT::PAR::FluidPoroRelPermeabilityLaw::CreateRelPermeabilityLaw(matdata->GetInt("RELPERMEABILITYLAWID"));

  // create viscosity law
  viscositylaw_ = MAT::PAR::FluidPoroViscosityLaw::CreateViscosityLaw(matdata->GetInt("VISCOSITYLAWID"));

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matdata->GetInt("DOFTYPEID"));

  switch (curmat->Type())
  {
  case INPAR::MAT::m_fluidporo_phasedof_diffpressure:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofDiffPressure(curmat));
    phasedof_ = static_cast<MAT::PAR::FluidPoroPhaseDofDiffPressure*>(curmat->Parameter());
    break;
  }
  case INPAR::MAT::m_fluidporo_phasedof_pressure:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofPressure(curmat));
    phasedof_ = static_cast<MAT::PAR::FluidPoroPhaseDofPressure*>(curmat->Parameter());
    break;
  }
  case INPAR::MAT::m_fluidporo_phasedof_saturation:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofSaturation(curmat));
    phasedof_ = static_cast<MAT::PAR::FluidPoroPhaseDofSaturation*>(curmat->Parameter());
    break;
  }
  default:
    dserror("invalid pressure-saturation law for material %d", curmat->Type());
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
  if(not isinit_)
  {
    phasedof_->Initialize();
    isinit_=true;
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

DRT::ParObject* MAT::FluidPoroSinglePhaseType::Create(const std::vector<char> & data)
{
  MAT::FluidPoroSinglePhase* fluid_poro = new MAT::FluidPoroSinglePhase();
  fluid_poro->Unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroSinglePhase::FluidPoroSinglePhase() :
  params_(NULL)
{
}

/*----------------------------------------------------------------------*
*   Create material with parameters                         vuong 08/16 |
*----------------------------------------------------------------------*/
MAT::FluidPoroSinglePhase::FluidPoroSinglePhase(MAT::PAR::FluidPoroSinglePhase* params) :
  params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                           vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroSinglePhase::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL)
    matid = params_->Id(); // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
* unpack material                                           vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroSinglePhase::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  if (DRT::Problem::Instance()->Materials()->Num() != 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::FluidPoroSinglePhase*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }

  if (position != data.size())
  dserror("Mismatch in size of data %d <-> %d",data.size(),position);
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
    Epetra_SerialDenseMatrix& dofmat,
    int numphase) const
{
  params_->phasedof_->FillDoFMatrix(dofmat,numphase);
  return;
}

/*----------------------------------------------------------------------*
 *  Evaluate generalized pressure of a phase                vuong 08/16 |
*----------------------------------------------------------------------*/
double MAT::FluidPoroSinglePhase::EvaluateGenPressure(
    int phasenum,
    const std::vector<double>& state) const
{
  return params_->phasedof_->EvaluateGenPressure(phasenum,state);
}

/*----------------------------------------------------------------------*
 *   Evaluate saturation of the phase                       vuong 08/16 |
*----------------------------------------------------------------------*/
double MAT::FluidPoroSinglePhase::EvaluateSaturation(
    int phasenum,
    const std::vector<double>& state,
    const std::vector<double>& pressure) const
{
  return params_->phasedof_->EvaluateSaturation(phasenum,state,pressure);
}

/*--------------------------------------------------------------------------*
 *  Evaluate derivative of saturation w.r.t. pressure           vuong 08/16 |
*---------------------------------------------------------------------------*/
double MAT::FluidPoroSinglePhase::EvaluateDerivOfSaturationWrtPressure(
    int phasenum,
    int doftoderive,
    const std::vector<double>& pressure) const
{
  return params_->phasedof_->EvaluateDerivOfSaturationWrtPressure(phasenum,doftoderive,pressure);
}

/*--------------------------------------------------------------------------*
 *  Evaluate 2nd derivative of saturation w.r.t. pressure  kremheller 05/17 |
*---------------------------------------------------------------------------*/
double MAT::FluidPoroSinglePhase::EvaluateSecondDerivOfSaturationWrtPressure(
    int phasenum,
    int firstdoftoderive,
    int seconddoftoderive,
    const std::vector<double>& pressure) const
{
  return params_->phasedof_->EvaluateSecondDerivOfSaturationWrtPressure(phasenum,firstdoftoderive,seconddoftoderive,pressure);
}

/*----------------------------------------------------------------------------------------*
 * Evaluate derivative of degree of freedom with respect to pressure          vuong 08/16 |
*----------------------------------------------------------------------------------------*/
double MAT::FluidPoroSinglePhase::EvaluateDerivOfDofWrtPressure(
    int phasenum,
    int doftoderive,
    const std::vector<double>& state) const
{
  return params_->phasedof_->EvaluateDerivOfDofWrtPressure(phasenum,doftoderive,state);
}

/*----------------------------------------------------------------------*
 *  constructor (public)                               kremheller 10/17 |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroSingleVolFrac::FluidPoroSingleVolFrac(Teuchos::RCP<MAT::PAR::Material> matdata) :
  Parameter(matdata),
  density_(matdata->GetDouble("DENSITY")),
  diffusivity_(matdata->GetDouble("DIFFUSIVITY")),
  scalardependentflux_(matdata->GetInt("AddScalarDependentFlux")),
  numscal_(matdata->GetInt("NUMSCAL")),
  scalardiffs_(*(matdata->Get<std::vector<double> >("SCALARDIFFS"))),
  omega_half_(*(matdata->Get<std::vector<double> >("OMEGA_HALF"))),
  isinit_(false)
{
  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("Sorry dude, cannot work out problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("Sorry dude, no materials defined.");

  // safety checks
  if(numscal_ < 0)
    dserror("NUMSCAL smaller than zero not possible");

  if(scalardependentflux_ && numscal_ == 0)
    dserror("AddScalarDependentFlux has been set to YES, but NUMSCAL is equal to zero");

  if(scalardependentflux_ && scalardiffs_.size() == 0)
    dserror("AddScalarDependentFlux has been set to YES, but length of SCALARDIFFS is equal to zero");

  if(scalardependentflux_ && omega_half_.size() == 0)
    dserror("AddScalarDependentFlux has been set to YES, but length of OMEGA_HALF is equal to zero");

  if(!scalardependentflux_ && numscal_ > 0)
    dserror("AddScalarDependentFlux has been set to NO, but NUMSCAL is greater than zero");

  if(!scalardependentflux_ && scalardiffs_.size() > 0)
    dserror("AddScalarDependentFlux has been set to NO, but length of SCALARDIFFS is greater than zero");

  if(!scalardependentflux_ && omega_half_.size() > 0)
    dserror("AddScalarDependentFlux has been set to NO, but length of OMEGA_HALF is greater than zero");

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

DRT::ParObject* MAT::FluidPoroSingleVolFracType::Create(const std::vector<char> & data)
{
  MAT::FluidPoroSingleVolFrac* fluid_poro = new MAT::FluidPoroSingleVolFrac();
  fluid_poro->Unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                             kremheller 10/17 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroSingleVolFrac::FluidPoroSingleVolFrac() :
  params_(NULL)
{
}

/*----------------------------------------------------------------------*
*   Create material with parameters                    kremheller 10/17 |
*----------------------------------------------------------------------*/
MAT::FluidPoroSingleVolFrac::FluidPoroSingleVolFrac(MAT::PAR::FluidPoroSingleVolFrac* params) :
  params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                      kremheller 10/17 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroSingleVolFrac::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL)
    matid = params_->Id(); // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
* unpack material                                      kremheller 10/17 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroSingleVolFrac::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  if (DRT::Problem::Instance()->Materials()->Num() != 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::FluidPoroSingleVolFrac*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }

  if (position != data.size())
  dserror("Mismatch in size of data %d <-> %d",data.size(),position);
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
MAT::PAR::FluidPoroVolFracPressure::FluidPoroVolFracPressure(Teuchos::RCP<MAT::PAR::Material> matdata) :
  Parameter(matdata),
  permeability_(matdata->GetDouble("PERMEABILITY")),
  min_volfrac_(matdata->GetDouble("MIN_VOLFRAC")),
  isinit_(false)
{
  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("Sorry dude, cannot work out problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("Sorry dude, no materials defined.");

  // create viscosity law
  viscositylaw_ = MAT::PAR::FluidPoroViscosityLaw::CreateViscosityLaw(matdata->GetInt("VISCOSITYLAWID"));

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

DRT::ParObject* MAT::FluidPoroVolFracPressureType::Create(const std::vector<char> & data)
{
  MAT::FluidPoroVolFracPressure* fluid_poro = new MAT::FluidPoroVolFracPressure();
  fluid_poro->Unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                             kremheller 02/18 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroVolFracPressure::FluidPoroVolFracPressure() :
  params_(NULL)
{
}

/*----------------------------------------------------------------------*
*   Create material with parameters                    kremheller 02/18 |
*----------------------------------------------------------------------*/
MAT::FluidPoroVolFracPressure::FluidPoroVolFracPressure(MAT::PAR::FluidPoroVolFracPressure* params) :
  params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                      kremheller 02/18 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroVolFracPressure::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL)
    matid = params_->Id(); // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
* unpack material                                      kremheller 02/18 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroVolFracPressure::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  if (DRT::Problem::Instance()->Materials()->Num() != 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::FluidPoroVolFracPressure*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }

  if (position != data.size())
  dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

/*----------------------------------------------------------------------*
 *  initialize                                         kremheller 02/18 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroVolFracPressure::Initialize()
{
  params_->Initialize();
  return;
}
