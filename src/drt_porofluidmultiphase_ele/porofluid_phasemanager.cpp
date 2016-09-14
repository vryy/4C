/*----------------------------------------------------------------------*/
/*!
 \file porofluid_phasemanager.cpp

 \brief manager class for handling the phases and their dofs on element level

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#include "porofluid_phasemanager.H"

#include "porofluid_reactionevaluator.H"
#include "porofluidmultiphase_ele_calc_utils.H"

#include "../drt_mat/fluidporo_multiphase.H"
#include "../drt_mat/fluidporo_singlephase.H"

#include <Epetra_SerialDenseSolver.h>


/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidPhaseManager::PoroFluidPhaseManager(int numphases):
    dof2pres_(Teuchos::null),
    numphases_(numphases),
    genpressure_(numphases,0.0),
    pressure_(numphases,0.0),
    saturation_(numphases,0.0),
    solidpressure_(0.0),
    pressurederiv_(Teuchos::null),
    saturationderiv_(Teuchos::null),
    solidpressurederiv_(Teuchos::null),
    solidpressurederivderiv_(Teuchos::null),
    statesset_(false),
    reactionevaluator_(Teuchos::null)
{
  // initialize matrixes and vectors
  pressurederiv_= Teuchos::rcp(new Epetra_SerialDenseMatrix(numphases,numphases));
  saturationderiv_= Teuchos::rcp(new Epetra_SerialDenseMatrix(numphases,numphases));
  solidpressurederiv_= Teuchos::rcp(new Epetra_SerialDenseVector(numphases));
  solidpressurederivderiv_= Teuchos::rcp(new Epetra_SerialDenseMatrix(numphases,numphases));

  //  matrix holding the conversion from pressures and dofs
  dof2pres_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(numphases_,numphases_));

  return;
}

/*----------------------------------------------------------------------*
 | setup                                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidPhaseManager::Setup(
    const MAT::Material& material)
{
  // check the material
  if(material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase and
     material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("only poro multiphase and poro multiphase reactions material valid");

  // cast
  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  //safety check
  if(numphases_ != multiphasemat.NumMat())
    dserror("Number of phases given by the poro multiphase material does not match number "
        "of DOFs (%i phases and %i DOFs)!"
        , numphases_, multiphasemat.NumMat());

  //  matrix holding the conversion from pressures and dofs
  // reset
  dof2pres_->Scale(0.0);

  for(int iphase=0; iphase<numphases_; iphase++)
  {
    // get the single phase material
    const MAT::FluidPoroSinglePhase& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMultiMaterial(multiphasemat,iphase);

    // consistency checks
    if( singlephasemat.PoroDofType() == INPAR::MAT::m_fluidporo_phasedof_pressuresum
        and iphase!= numphases_-1)
      dserror("Only the last material in the list of poro multiphase materials needs to be of type 'PressureSum'!");
    if(singlephasemat.PoroDofType() !=  INPAR::MAT::m_fluidporo_phasedof_pressuresum
        and iphase== numphases_-1
        and numphases_!=1)
      dserror("The last material in the list of poro multiphase materials needs to be of type 'PressureSum'!");

    // fill the coefficients into matrix
    singlephasemat.FillDoFMatrix(*dof2pres_,iphase);
  }

  // invert dof2pres_ to get conversion from dofs to pressures
  {
    Epetra_SerialDenseSolver inverse;
    inverse.SetMatrix(*dof2pres_);
    int err = inverse.Invert();
    if (err != 0)
      dserror("Inversion of matrix for DOF transform failed with errorcode %d. Is your system of DOFs linear independent?",err);
  }

  if(material.MaterialType() == INPAR::MAT::m_fluidporo_multiphase_reactions)
  {
    reactionevaluator_ = Teuchos::rcp(new DRT::ELEMENTS::POROFLUIDEVALUATOR::ReactionEvaluator());
    reactionevaluator_->Setup(*this,material);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate pressures, saturations and derivatives at GP     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidPhaseManager::EvaluateGPState(
    const MAT::Material& material,
    const std::vector<double>& phinp)
{
  if(statesset_ == true)
    dserror("state has already been set!");
  if(numphases_ != (int)phinp.size())
    dserror("Length of phinp vector is not equal to the number of phases");

  // cast the material to multiphase material
  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  //resize all vectors
  genpressure_.resize(numphases_,0.0);
  pressure_.resize(numphases_,0.0);
  saturation_.resize(numphases_,0.0);

  // evaluate the pressures
  for(int iphase=0; iphase<numphases_; iphase++)
  {
    //get the single phase material
    const MAT::FluidPoroSinglePhase& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMultiMaterial(multiphasemat,iphase);

    // evaluate generalized pressure (i.e. some kind of linear combination of the true pressures)
    genpressure_[iphase] = singlephasemat.EvaluateGenPressure(iphase,phinp);
  }

  //! transform generalized pressures to true pressure values
  TransformGenPresToTruePres(genpressure_,pressure_);

  // explicit evaluation of saturation
  saturation_[numphases_-1] = 1.0;
  for(int iphase=0; iphase<numphases_-1; iphase++)
  {
    // get the single phase material
    const MAT::FluidPoroSinglePhase& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMultiMaterial(multiphasemat,iphase);

    saturation_[iphase] = singlephasemat.EvaluateSaturation(iphase,phinp,pressure_);
    // the saturation of the last phase is 1.0- (sum of all saturations)
    saturation_[numphases_-1] -= saturation_[iphase];
  }

  // solid pressure = sum (S_i*p_i)
  solidpressure_ = std::inner_product(saturation_.begin(),saturation_.end(),pressure_.begin(),0.0);

  // calculate the derivative of the pressure (actually first its inverse)
  for(int iphase=0; iphase<numphases_; iphase++)
    for(int jphase=0; jphase<numphases_; jphase++)
      (*pressurederiv_)(iphase,jphase)   = EvaluateDerivOfDofWrtPressure(material,iphase,jphase,phinp);

  // now invert the derivatives of the dofs w.r.t. pressure to get the derivatives
  // of the pressure w.r.t. the dofs
  {
    Epetra_SerialDenseSolver inverse;
    inverse.SetMatrix(*pressurederiv_);
    int err = inverse.Invert();
    if (err != 0)
      dserror("Inversion of matrix for pressure derivative failed with error code %d.",err);
  }

  // calculate derivatives of saturation w.r.t. pressure
  Epetra_SerialDenseMatrix deriv(numphases_,numphases_);
  for(int iphase=0; iphase<numphases_-1; iphase++)
  {
    for(int jphase=0; jphase<numphases_; jphase++)
    {
      const double saturationderiv = EvaluateDerivOfSaturationWrtPressure(material,iphase,jphase,phinp);
      deriv(iphase,jphase) = saturationderiv;
      // the saturation of the last phase is 1.0- (sum of all saturations)
      // -> the derivative of this saturation = -1.0 (sum of all saturation derivatives)
      deriv(numphases_-1,jphase) += -1.0*saturationderiv;
    }
  }

  // chain rule: the derivative of saturation w.r.t. dof =
  // (derivative of saturation w.r.t. pressure) * (derivative of pressure w.r.t. dof)
  saturationderiv_->Multiply('N','N',1.0,deriv,*pressurederiv_,0.0);

  // compute derivative of solid pressure w.r.t. dofs with product rule
  solidpressurederiv_->Scale(0.0);
  for(int iphase=0; iphase<numphases_; iphase++)
    for(int jphase=0; jphase<numphases_; jphase++)
      (*solidpressurederiv_)(iphase) +=   (*pressurederiv_)(jphase,iphase)*saturation_[jphase]
                                        + (*saturationderiv_)(jphase,iphase)*pressure_[jphase];

  // compute second derivative of solid pressure w.r.t. dofs with product rule
  //TODO also include second derivs of pressure and saturation
  solidpressurederivderiv_->Scale(0.0);
  for(int iphase=0; iphase<numphases_; iphase++)
    for(int jphase=0; jphase<numphases_; jphase++)
      for(int kphase=0; kphase<numphases_; kphase++)
        (*solidpressurederivderiv_)(jphase,kphase) +=
              (*pressurederiv_)(iphase,kphase)*(*saturationderiv_)(iphase,jphase)
            + (*saturationderiv_)(iphase,kphase)*(*pressurederiv_)(iphase,jphase);

  // done
  statesset_ = true;

  return;
}

/*----------------------------------------------------------------------*
 | evaluate pressures, saturations and derivatives at GP     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidPhaseManager::EvaluateGPReactions(
    const MAT::Material& material,
    double porosity,
    const std::vector<double>& scalars)
{
  if(not statesset_)
    dserror("Gauss point states have not been set!");

  if(reactionevaluator_ != Teuchos::null)
    reactionevaluator_->EvaluateGPState(
        *this,
        material,
         porosity,
        scalars);

  return;
}

/*----------------------------------------------------------------------*
 | zero all values at GP                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidPhaseManager::ClearGPState()
{
  // this is just a safety call. All quantities should be recomputed
  // in the next EvaluateGPState call anyway, but you never know ...

  //zero everything
  std::fill(genpressure_.begin(), genpressure_.end(), 0.0);
  std::fill(pressure_.begin(), pressure_.end(), 0.0);
  std::fill(saturation_.begin(), saturation_.end(), 0.0);
  solidpressure_=0.0;
  std::vector<double> zero(pressure_.size(),0.0);
  pressurederiv_->Scale(0.0);
  saturationderiv_->Scale(0.0);
  solidpressurederiv_->Scale(0.0);
  solidpressurederivderiv_->Scale(0.0);

  if(reactionevaluator_ != Teuchos::null)
    reactionevaluator_->ClearGPState();

  // states are reseted
  statesset_ = false;

  return;
}

/*----------------------------------------------------------------------*
 | reaction query                                            vuong 08/16 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::PoroFluidPhaseManager::IsReactive(int phasenum)
{
  if(reactionevaluator_==Teuchos::null)
    return false;
  else
    return reactionevaluator_->IsReactive(phasenum);
}

/*----------------------------------------------------------------------*
 | transform generalized pressures to true pressures        vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidPhaseManager::TransformGenPresToTruePres(
    const std::vector<double>& phinp,
    std::vector<double>& phi_transformed)
{
  //simple matrix vector product
  phi_transformed.resize(phinp.size());
  for(int i=0;i<numphases_;i++)
    for(int j=0;j<numphases_;j++)
      phi_transformed[i] += (*dof2pres_)(i,j)*phinp[j];
  return;
}

/*---------------------------------------------------------------------------*
 *  evaluate derivative of saturation with respect to pressure   vuong 08/16 |
*---------------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::EvaluateDerivOfSaturationWrtPressure(
    const MAT::Material& material,
    int phasenum,
    int doftoderive,
    const std::vector<double>& state) const
{
  //get the single phase material
  const MAT::FluidPoroSinglePhase& singlephasemat =
      POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,phasenum);

  return singlephasemat.EvaluateDerivOfSaturationWrtPressure(phasenum,doftoderive,state);
}

/*-----------------------------------------------------------------------------------*
 * evaluate derivative of degree of freedom with respect to pressure     vuong 08/16 |
*------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::EvaluateDerivOfDofWrtPressure(
    const MAT::Material& material,
    int phasenum,
    int doftoderive,
    const std::vector<double>& state) const
{
  //get the single phase material
  const MAT::FluidPoroSinglePhase& singlephasemat =
      POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,phasenum);

  return singlephasemat.EvaluateDerivOfDofWrtPressure(phasenum,doftoderive,state);
}


/*----------------------------------------------------------------------*
 * get solid pressure                                        vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::SolidPressure() const
{
  if(not statesset_)
    dserror("Gauss point states have not been set!");

  return solidpressure_;
}

/*----------------------------------------------------------------------*
 *    get saturation of phase                               vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::Saturation(int phasenum) const
{
  if(not statesset_)
    dserror("Gauss point states have not been set!");

  return saturation_[phasenum];
}

/*----------------------------------------------------------------------*
 *  get pressure of phase                                  vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::Pressure(int phasenum) const
{
  if(not statesset_)
    dserror("Gauss point states have not been set!");

  return pressure_[phasenum];
}

/*----------------------------------------------------------------------*
 *  get derivative of saturation of phase                   vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::SaturationDeriv(
    int phasenum,
    int doftoderive) const
{
  if(not statesset_)
    dserror("Gauss point states have not been set!");

  return (*saturationderiv_)(phasenum,doftoderive);
}

/*----------------------------------------------------------------------*
 *  get derivative of pressure of phase                      vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::PressureDeriv(
    int phasenum,
    int doftoderive) const
{
  if(not statesset_)
    dserror("Gauss point states have not been set!");

  return (*pressurederiv_)(phasenum,doftoderive);
}

/*----------------------------------------------------------------------*
 * get derivative of solid pressure                          vuong 08/16 |
*-----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::SolidPressureDeriv(
    int doftoderive) const
{
  if(not statesset_)
    dserror("Gauss point states have not been set!");

  return (*solidpressurederiv_)(doftoderive);
}

/*---------------------------------------------------------------------------*
 * get second derivative of solid pressure                       vuong 08/16 |
*---------------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::SolidPressureDerivDeriv(
    int doftoderive,
    int doftoderive2) const
{
  if(not statesset_)
    dserror("Gauss point states have not been set!");

  return (*solidpressurederivderiv_)(doftoderive,doftoderive2);
}


/*----------------------------------------------------------------------*
 *   get bulk modulus of phase 'phasenum'                    vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::Bulkmodulus(
    const MAT::Material& material,
    int phasenum) const
{

  //get the single phase material
  const MAT::FluidPoroSinglePhase& singlephasemat =
      POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,phasenum);

  return singlephasemat.Bulkmodulus();
}

/*----------------------------------------------------------------------*
 *   get density of phase 'phasenum'                    vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::Density(
    const MAT::Material& material,
    int phasenum) const
{

  //get the single phase material
  const MAT::FluidPoroSinglePhase& singlephasemat =
      POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,phasenum);

  return singlephasemat.Density();
}

/*----------------------------------------------------------------------*
 *   get reaction term of phase 'phasenum'                   vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::ReacTerm(
    int phasenum) const
{
  if(reactionevaluator_ == Teuchos::null)
    dserror("Reaction evaluator is null pointer!");

  return reactionevaluator_->GetReacTerm(phasenum);
}

/*----------------------------------------------------------------------*
 *   get reaction term of phase 'phasenum'                   vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::ReacDeriv(
    int phasenum,
    int doftoderive) const
{
  if(reactionevaluator_ == Teuchos::null)
    dserror("Reaction evaluator is null pointer!");

  return reactionevaluator_->GetReacDeriv(phasenum,doftoderive);
}


/*----------------------------------------------------------------------*
 * evaluate relative diffusivity of phase 'phasenum'         vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidPhaseManager::EvaluateRelDiffusivity(
    const MAT::Material& material,
    int phasenum) const
{
  if(not statesset_)
    dserror("Gauss point states have not been set!");

  //get the single phase material
  const MAT::FluidPoroSinglePhase& singlephasemat =
      POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,phasenum);

  return singlephasemat.Permeability()/singlephasemat.Viscosity();
}
