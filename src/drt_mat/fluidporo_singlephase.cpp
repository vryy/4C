/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_singlephase.cpp

 \brief a single phase material of a multiphase porous fluid

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
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
  viscosity_(matdata->GetDouble("DYNVISCOSITY")),
  density_(matdata->GetDouble("DENSITY")),
  permeability_(matdata->GetDouble("PERMEABILITY")),
  bulkmodulus_(matdata->GetDouble("BULKMODULUS")),
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
  case INPAR::MAT::m_fluidporo_phasedof_pressuresum:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofPressureSum(curmat));
    phasedof_ = static_cast<MAT::PAR::FluidPoroPhaseDofPressureSum*>(curmat->Parameter());
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
 *  fill the dof matrix with the phase dofs                 vuong 08/16 |
*----------------------------------------------------------------------*/
INPAR::MAT::MaterialType MAT::FluidPoroSinglePhase::PoroDofType() const
{
  return params_->phasedof_->Type();
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
    const std::vector<double>& state) const
{
  return params_->phasedof_->EvaluateDerivOfSaturationWrtPressure(phasenum,doftoderive,state);
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
