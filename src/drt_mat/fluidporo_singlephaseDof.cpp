/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_singlephaseDof.cpp

 \brief a material defining the degree of freedom of a single phase of
        a multiphase porous fluid

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "fluidporo_singlephaseDof.H"

#include "fluidporo_singlephaselaw.H"

#include "matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 *  constructor (public)                               vuong 08/16      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseDof::FluidPoroPhaseDof(Teuchos::RCP<MAT::PAR::Material> matdata) :
  Parameter(matdata)
{
}

/*----------------------------------------------------------------------*
 *  factory method for phase dof                       vuong 08/16      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseDof* MAT::PAR::FluidPoroPhaseDof::CreatePhaseDof(
    int phasedofId)
{
  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("Sorry dude, cannot work out problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("Sorry dude, no materials defined.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat = DRT::Problem::Instance(probinst)->Materials()->ById(phasedofId);

  // phase law
  MAT::PAR::FluidPoroPhaseDof* phasedof = NULL;

   switch (curmat->Type())
   {
   case INPAR::MAT::m_fluidporo_phasedof_diffpressure:
   {
     if (curmat->Parameter() == NULL)
       curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofDiffPressure(curmat));
     phasedof = static_cast<MAT::PAR::FluidPoroPhaseDofDiffPressure*>(curmat->Parameter());
     break;
   }
   case INPAR::MAT::m_fluidporo_phasedof_pressure:
   {
     if (curmat->Parameter() == NULL)
       curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofPressure(curmat));
     phasedof = static_cast<MAT::PAR::FluidPoroPhaseDofPressure*>(curmat->Parameter());
     break;
   }
   case INPAR::MAT::m_fluidporo_phasedof_saturation:
   {
     if (curmat->Parameter() == NULL)
       curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofSaturation(curmat));
     phasedof = static_cast<MAT::PAR::FluidPoroPhaseDofSaturation*>(curmat->Parameter());
     break;
   }
   default:
     dserror("invalid pressure-saturation law for material %d", curmat->Type());
     break;
   }

  return phasedof;
}

/************************************************************************/
/************************************************************************/

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseDofDiffPressure::FluidPoroPhaseDofDiffPressure(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: FluidPoroPhaseDof(matdata),
  diffpresCoeffs_(matdata->Get<std::vector<int> >("PRESCOEFF")),
  phaselawId_(matdata->GetInt("PHASELAWID"))
{
  phaselaw_ = MAT::PAR::FluidPoroPhaseLaw::CreatePhaseLaw(phaselawId_);
  return;
}

/*----------------------------------------------------------------------*
 *  Initialize                                               vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroPhaseDofDiffPressure::Initialize()
{
  phaselaw_->Initialize();
  return;
}

/*----------------------------------------------------------------------*
 *  return phase law type                                  vuong 08/16 |
*----------------------------------------------------------------------*/
INPAR::MAT::MaterialType MAT::PAR::FluidPoroPhaseDofDiffPressure::PoroPhaseLawType() const
{
  return phaselaw_->Type();
}

/*----------------------------------------------------------------------*
 *  fill the dof matrix with the phase dofs                 vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroPhaseDofDiffPressure::FillDoFMatrix(
    Epetra_SerialDenseMatrix& dofmat,
    int numphase) const
{
  // safety check
  if((int)diffpresCoeffs_->size() != dofmat.N())
    dserror("Number of phases given by the poro singlephase material %i "
        "does not match number of DOFs (%i phases and %i DOFs)!",
        phaselaw_->Id(), diffpresCoeffs_->size(), dofmat.N());

  // fill pressure coefficients into matrix
  for(size_t i=0; i<diffpresCoeffs_->size();i++)
  {
    const int val = (*diffpresCoeffs_)[i];
    if(val!=0)
      dofmat(numphase,i) = val;
  }
}

/*----------------------------------------------------------------------*
 *  Evaluate generalized pressure of a phase                vuong 08/16 |
*----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofDiffPressure::EvaluateGenPressure(
    int phasenum,
    const std::vector<double>& state) const
{
  // return the corresponding dof value
  return state[phasenum];
}


/*----------------------------------------------------------------------*
 *   Evaluate saturation of the phase                       vuong 08/16 |
*----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofDiffPressure::EvaluateSaturation(
    int phasenum,
    const std::vector<double>& state,
    const std::vector<double>& pressure) const
{
  // call the phase law
  return  phaselaw_->EvaluateSaturation(pressure);
}


/*--------------------------------------------------------------------------*
 *  Evaluate derivative of saturation w.r.t. pressure           vuong 08/16 |
*---------------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofDiffPressure::EvaluateDerivOfSaturationWrtPressure(
    int phasenum,
    int doftoderive,
    const std::vector<double>& state) const
{
  // call the phase law
  return phaselaw_->EvaluateDerivOfSaturationWrtPressure(doftoderive,state);
}

/*--------------------------------------------------------------------------*
 *  Evaluate 2nd derivative of saturation w.r.t. pressure  kremheller 05/17 |
*---------------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofDiffPressure::EvaluateSecondDerivOfSaturationWrtPressure(
    int phasenum,
    int firstdoftoderive,
    int seconddoftoderive,
    const std::vector<double>& pressure) const
{
  // call the phase law
  return phaselaw_->EvaluateSecondDerivOfSaturationWrtPressure(firstdoftoderive,seconddoftoderive,pressure);
}

/*----------------------------------------------------------------------------------------*
 * Evaluate derivative of degree of freedom with respect to pressure          vuong 08/16 |
*----------------------------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofDiffPressure::EvaluateDerivOfDofWrtPressure(
    int phasenum,
    int doftoderive,
    const std::vector<double>& state) const
{
  // derivative is the corresponding coefficient
  return  (*diffpresCoeffs_)[doftoderive];
}

/************************************************************************/
/************************************************************************/

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseDofPressure::FluidPoroPhaseDofPressure(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: FluidPoroPhaseDof(matdata),
  phaselawId_(matdata->GetInt("PHASELAWID"))
{
  phaselaw_ = MAT::PAR::FluidPoroPhaseLaw::CreatePhaseLaw(phaselawId_);
  return;
}

/*----------------------------------------------------------------------*
 *  Initialize                                               vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroPhaseDofPressure::Initialize()
{
  phaselaw_->Initialize();
  return;
}

/*----------------------------------------------------------------------*
 *  return phase law type                                  vuong 08/16 |
*----------------------------------------------------------------------*/
INPAR::MAT::MaterialType MAT::PAR::FluidPoroPhaseDofPressure::PoroPhaseLawType() const
{
  return phaselaw_->Type();
}

/*----------------------------------------------------------------------*
 *  fill the dof matrix with the phase dofs                 vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroPhaseDofPressure::FillDoFMatrix(
    Epetra_SerialDenseMatrix& dofmat,
    int numphase) const
{
  // just mark the corresponding entry in the matrix
  dofmat(numphase,numphase) = 1.0;
}

/*----------------------------------------------------------------------*
 *  Evaluate generalized pressure of a phase                vuong 08/16 |
*----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofPressure::EvaluateGenPressure(
    int phasenum,
    const std::vector<double>& state) const
{
  // return the corresponding dof value
  return state[phasenum];
}


/*----------------------------------------------------------------------*
 *   Evaluate saturation of the phase                       vuong 08/16 |
*----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofPressure::EvaluateSaturation(
    int phasenum,
    const std::vector<double>& state,
    const std::vector<double>& pressure) const
{
  // call the phase law
  return phaselaw_->EvaluateSaturation(pressure);
}


/*--------------------------------------------------------------------------*
 *  Evaluate derivative of saturation w.r.t. pressure           vuong 08/16 |
*---------------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofPressure::EvaluateDerivOfSaturationWrtPressure(
    int phasenum,
    int doftoderive,
    const std::vector<double>& state) const
{
  // call the phase law
  return phaselaw_->EvaluateDerivOfSaturationWrtPressure(doftoderive,state);
}

/*--------------------------------------------------------------------------*
 *  Evaluate 2nd derivative of saturation w.r.t. pressure  kremheller 05/17 |
*---------------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofPressure::EvaluateSecondDerivOfSaturationWrtPressure(
    int phasenum,
    int firstdoftoderive,
    int seconddoftoderive,
    const std::vector<double>& pressure) const
{
  // call the phase law
  return phaselaw_->EvaluateSecondDerivOfSaturationWrtPressure(firstdoftoderive,seconddoftoderive,pressure);
}

/*----------------------------------------------------------------------------------------*
 * Evaluate derivative of degree of freedom with respect to pressure          vuong 08/16 |
*----------------------------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofPressure::EvaluateDerivOfDofWrtPressure(
    int phasenum,
    int doftoderive,
    const std::vector<double>& state) const
{
  double presurederiv = 0.0;

  // respective derivative of w.r.t. is either 0 or 1
  if(phasenum==doftoderive)
    presurederiv = 1.0;

  return presurederiv;
}

/************************************************************************/
/************************************************************************/

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroPhaseDofSaturation::FluidPoroPhaseDofSaturation(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: FluidPoroPhaseDof(matdata),
  phaselawId_(matdata->GetInt("PHASELAWID"))
{

  phaselaw_ = MAT::PAR::FluidPoroPhaseLaw::CreatePhaseLaw(phaselawId_);
  return;
}

/*----------------------------------------------------------------------*
 *  Initialize                                               vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroPhaseDofSaturation::Initialize()
{
  phaselaw_->Initialize();
  return;
}

/*----------------------------------------------------------------------*
 *  return phase law type                                  vuong 08/16 |
*----------------------------------------------------------------------*/
INPAR::MAT::MaterialType MAT::PAR::FluidPoroPhaseDofSaturation::PoroPhaseLawType() const
{
  return phaselaw_->Type();
}

/*----------------------------------------------------------------------*
 *  fill the dof matrix with the phase dofs                 vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroPhaseDofSaturation::FillDoFMatrix(
    Epetra_SerialDenseMatrix& dofmat,
    int numphase) const
{
  // get pressure coefficients of phase law
  const std::vector<int>* presIDs = phaselaw_->PresIds();

  // safety check
  if((int)presIDs->size() != dofmat.N())
    dserror("Number of phases given by the poro phase law material %i "
        "does not match number of DOFs (%i phases and %i DOFs)!",
        phaselaw_->Id(), presIDs->size(), dofmat.N());

  // fill pressure coefficients of phase law into matrix
  for(size_t i=0; i<presIDs->size();i++)
  {
    const int val = (*presIDs)[i];
    if(val!=0)
      dofmat(numphase,i) = val;
  }
}

/*----------------------------------------------------------------------*
 *  Evaluate generalized pressure of a phase                vuong 08/16 |
*----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofSaturation::EvaluateGenPressure(
    int phasenum,
    const std::vector<double>& state) const
{
  // evaluate the phase law for the generalized (i.e. some differential pressure)
  // the phase law depends on
  return phaselaw_->EvaluateGenPressure(state[phasenum]);
}


/*----------------------------------------------------------------------*
 *   Evaluate saturation of the phase                       vuong 08/16 |
*----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofSaturation::EvaluateSaturation(
    int phasenum,
    const std::vector<double>& state,
    const std::vector<double>& pressure) const
{
  // get the corresponding dof value
  return state[phasenum];
}


/*--------------------------------------------------------------------------*
 *  Evaluate derivative of saturation w.r.t. pressure           vuong 08/16 |
*---------------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofSaturation::EvaluateDerivOfSaturationWrtPressure(
    int phasenum,
    int doftoderive,
    const std::vector<double>& state) const
{
  // call the phase law
  return phaselaw_->EvaluateDerivOfSaturationWrtPressure(doftoderive,state);
}

/*--------------------------------------------------------------------------*
 *  Evaluate 2nd derivative of saturation w.r.t. pressure  kremheller 05/17 |
*---------------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofSaturation::EvaluateSecondDerivOfSaturationWrtPressure(
    int phasenum,
    int firstdoftoderive,
    int seconddoftoderive,
    const std::vector<double>& pressure) const
{
  // call the phase law
  return phaselaw_->EvaluateSecondDerivOfSaturationWrtPressure(firstdoftoderive,seconddoftoderive,pressure);
}

/*----------------------------------------------------------------------------------------*
 * Evaluate derivative of degree of freedom with respect to pressure          vuong 08/16 |
*----------------------------------------------------------------------------------------*/
double MAT::PAR::FluidPoroPhaseDofSaturation::EvaluateDerivOfDofWrtPressure(
    int phasenum,
    int doftoderive,
    const std::vector<double>& state) const
{
  // call the phase law for the derivative
  return phaselaw_->EvaluateDerivOfSaturationWrtPressure(doftoderive,state);
}

/************************************************************************/
/************************************************************************/
