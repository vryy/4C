/*----------------------------------------------------------------------*/
/*!
 \file porofluid_phasemanager.cpp

 \brief manager class for handling the phases and their dofs on element level

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "porofluid_phasemanager.H"

#include "porofluid_variablemanager.H"
#include "porofluidmultiphase_ele_calc_utils.H"

#include "../drt_mat/fluidporo_multiphase.H"
#include "../drt_mat/fluidporo_singlephase.H"
#include "../drt_mat/structporo.H"
#include "../drt_mat/fluidporo_multiphase_reactions.H"
#include "../drt_mat/fluidporo_multiphase_singlereaction.H"

#include <Epetra_SerialDenseSolver.h>



/*----------------------------------------------------------------------*
 | factory method                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface >
DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface::CreatePhaseManager(
    const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
    int nsd,
    INPAR::MAT::MaterialType mattype,
    const POROFLUIDMULTIPHASE::Action&     action,
    int numphases)
{
  Teuchos::RCP< PhaseManagerInterface > phasemanager = Teuchos::null;

  // build the standard phase manager
  phasemanager = Teuchos::rcp(new PhaseManagerCore(numphases));

  return WrapPhaseManager(para,nsd,mattype,action,phasemanager);
}

/*----------------------------------------------------------------------*
 | factory method                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface >
DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface::WrapPhaseManager(
    const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
    int nsd,
    INPAR::MAT::MaterialType mattype,
    const POROFLUIDMULTIPHASE::Action&     action,
    Teuchos::RCP< PhaseManagerInterface > corephasemanager
    )
{
  Teuchos::RCP< PhaseManagerInterface > phasemanager = Teuchos::null;

  // determine action
  switch(action)
  {
  // calculate true pressures and saturation
  case POROFLUIDMULTIPHASE::calc_pres_and_sat:
  // calculate solid pressure
  case POROFLUIDMULTIPHASE::calc_solidpressure:
  {
    // no extensions needed
    phasemanager = corephasemanager;
    break;
  }
  // reconstruct velocities
  case POROFLUIDMULTIPHASE::recon_flux_at_nodes:
  {
    // derivatives needed
    phasemanager = Teuchos::rcp(new PhaseManagerDeriv(corephasemanager));
    // enhance by diffusion tensor
    switch(nsd)
    {
    case 1:
      phasemanager = Teuchos::rcp(new PhaseManagerDiffusion<1>(phasemanager));
      break;
    case 2:
      phasemanager = Teuchos::rcp(new PhaseManagerDiffusion<2>(phasemanager));
      break;
    case 3:
      phasemanager = Teuchos::rcp(new PhaseManagerDiffusion<3>(phasemanager));
      break;
    default:
      dserror("invalid dimension for creating phase manager!");
    }
    break;
  }
  // standard evaluate call
  case POROFLUIDMULTIPHASE::calc_mat_and_rhs:
  case POROFLUIDMULTIPHASE::calc_initial_time_deriv:
  case POROFLUIDMULTIPHASE::calc_fluid_coupl_mat:
  {
    // porosity (includes derivatves) needed
    phasemanager = Teuchos::rcp(new PhaseManagerDerivAndPorosity(corephasemanager));

    // enhance by diffusion tensor
    switch(nsd)
    {
    case 1:
      phasemanager = Teuchos::rcp(new PhaseManagerDiffusion<1>(phasemanager));
      break;
    case 2:
      phasemanager = Teuchos::rcp(new PhaseManagerDiffusion<2>(phasemanager));
      break;
    case 3:
      phasemanager = Teuchos::rcp(new PhaseManagerDiffusion<3>(phasemanager));
      break;
    default:
      dserror("invalid dimension for creating phase manager!");
    }

    if(mattype==INPAR::MAT::m_fluidporo_multiphase_reactions)
      // enhance by scalar handling capability
      phasemanager = Teuchos::rcp(new PhaseManagerReaction(phasemanager));
    break;
  }
  case POROFLUIDMULTIPHASE::calc_porosity:
  {
    // porosity (includes derivatves) needed
    phasemanager = Teuchos::rcp(new PhaseManagerDerivAndPorosity(corephasemanager));
    break;
  }
  default:
    dserror("unknown action %i for creating the phase manager!",action);
    break;
  } // switch(action)


  return phasemanager;
}

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::PhaseManagerCore(int numphases):
    numphases_(numphases),
    genpressure_(numphases,0.0),
    pressure_(numphases,0.0),
    saturation_(numphases,0.0),
    solidpressure_(0.0),
    invbulkmodulifluid_(numphases,0.0),
    invbulkmodulussolid_(0.0),
    ele_(NULL),
    isevaluated_(false),
    issetup_(false)
{
  return;
}

/*----------------------------------------------------------------------*
 | copy constructor                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::PhaseManagerCore(const PhaseManagerCore& old):
    numphases_(old.numphases_),
    genpressure_(old.genpressure_),
    pressure_(old.pressure_),
    saturation_(old.saturation_),
    solidpressure_(old.solidpressure_),
    invbulkmodulifluid_(old.invbulkmodulifluid_),
    invbulkmodulussolid_(old.invbulkmodulussolid_),
    ele_(NULL),
    isevaluated_(old.isevaluated_),
    issetup_(old.issetup_)
{
  return;
}

/*----------------------------------------------------------------------*
 | setup                                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::Setup(
    const DRT::Element* ele)
{

  dsassert(ele!=NULL,"Element is null pointer for setup of phase manager!");
  // save current element
  ele_ = ele;
  // get material
  const MAT::Material& material = *(ele_->Material());

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

  issetup_=true;
  return;
}

/*----------------------------------------------------------------------*
 | evaluate pressures, saturations and derivatives at GP     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::EvaluateGPState(
    double                 J,
    const VariableManagerMinAccess& varmanager)
{
  CheckIsSetup();

  if(isevaluated_ == true)
    dserror("state has already been set!");

  // get material
  dsassert(ele_->Material()!=Teuchos::null,"Material of element is null pointer!");
  const MAT::Material& material = *(ele_->Material());

  // access state vector
  const std::vector<double>& phinp = *varmanager.Phinp();

  if(numphases_ != (int)phinp.size())
    dserror("Length of phinp vector is not equal to the number of phases");

  // cast the material to multiphase material
  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  //resize all vectors
  genpressure_.resize(numphases_,0.0);
  pressure_.resize(numphases_,0.0);
  saturation_.resize(numphases_,0.0);
  invbulkmodulifluid_.resize(numphases_,0.0);

  // evaluate the pressures
  multiphasemat.EvaluateGenPressure(genpressure_,phinp);

  //! transform generalized pressures to true pressure values
  multiphasemat.TransformGenPresToTruePres(genpressure_,pressure_);

  // explicit evaluation of saturation
  multiphasemat.EvaluateSaturation(saturation_,phinp,pressure_);

  // solid pressure = sum (S_i*p_i)
  solidpressure_ = std::inner_product(saturation_.begin(),saturation_.end(),pressure_.begin(),0.0);

  //access second material in structure element
  dsassert(ele_->NumMaterial() > 1,"no second material defined for element ");
  dsassert(ele_->Material(1)->MaterialType() == INPAR::MAT::m_structporo,"invalid structure material for poroelasticity");

  // cast second material to poro material
  Teuchos::RCP<MAT::StructPoro> structmat = Teuchos::rcp_static_cast<MAT::StructPoro>(ele_->Material(1));

  invbulkmodulussolid_ = structmat->InvBulkmodulus();

  for(int iphase=0; iphase<numphases_; iphase++)
  {
    //get the single phase material
    const MAT::FluidPoroSinglePhase& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,iphase);
    invbulkmodulifluid_[iphase]=singlephasemat.InvBulkmodulus();
  }

  // done
  isevaluated_ = true;

  return;
}

/*----------------------------------------------------------------------*
 | zero all values at GP                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::ClearGPState()
{
  // this is just a safety call. All quantities should be recomputed
  // in the next EvaluateGPState call anyway, but you never know ...

  //zero everything
  std::fill(genpressure_.begin(), genpressure_.end(), 0.0);
  std::fill(pressure_.begin(), pressure_.end(), 0.0);
  std::fill(saturation_.begin(), saturation_.end(), 0.0);
  solidpressure_=0.0;
  std::vector<double> zero(pressure_.size(),0.0);
  invbulkmodulussolid_ = 0.0;
  std::fill(invbulkmodulifluid_.begin(), invbulkmodulifluid_.end(), 0.0);

  // states are reset
  isevaluated_ = false;

  return;
}

/*----------------------------------------------------------------------*
 * get solid pressure                                        vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::SolidPressure() const
{
  CheckIsEvaluated();

  return solidpressure_;
}

/*----------------------------------------------------------------------*
 *    get saturation of phase                               vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::Saturation(int phasenum) const
{
  CheckIsEvaluated();

  return saturation_[phasenum];
}

/*----------------------------------------------------------------------*
 *  get pressure of phase                                  vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::Pressure(int phasenum) const
{
  CheckIsEvaluated();

  return pressure_[phasenum];
}

/*----------------------------------------------------------------------*
 *    get saturation                                        vuong 08/16 |
*----------------------------------------------------------------------*/
const std::vector<double>& DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::Saturation() const
{
  CheckIsEvaluated();

  return saturation_;
}

/*----------------------------------------------------------------------*
 *  get pressure                                            vuong 08/16 |
*----------------------------------------------------------------------*/
const std::vector<double>& DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::Pressure() const
{
  CheckIsEvaluated();

  return pressure_;
}

/*----------------------------------------------------------------------*
 *   get bulk modulus of phase 'phasenum'                    vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::InvBulkmodulus(
    int phasenum) const
{
  CheckIsEvaluated();

  return invbulkmodulifluid_[phasenum];
}

/*----------------------------------------------------------------------*
 *  check if fluid phase 'phasenum' is incompressible  kremheller 06/17 |
*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::IncompressibleFluidPhase(
    int phasenum) const
{
  CheckIsEvaluated();

  return invbulkmodulifluid_[phasenum] < 1.0e-14;
}

/*----------------------------------------------------------------------*
 *   get inverse bulk modulus of solid phase                 vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::InvBulkmodulusSolid() const
{
  CheckIsEvaluated();

  return invbulkmodulussolid_;
}

/*----------------------------------------------------------------------*
 *  check if solid is incompressible                   kremheller 06/17 |
*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::IncompressibleSolid() const
{
  CheckIsEvaluated();

  return invbulkmodulussolid_ < 1.0e-14;
}

/*----------------------------------------------------------------------*
 *   get density of phase 'phasenum'                    vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::Density(
    const MAT::Material& material,
    int phasenum) const
{

  //get the single phase material
  const MAT::FluidPoroSinglePhase& singlephasemat =
      POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,phasenum);

  return singlephasemat.Density();
}

/*----------------------------------------------------------------------*
 *   get density of phase 'phasenum'                    vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::Density(
    int phasenum) const
{
  return Density(*Element()->Material(),phasenum);
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDeriv::PhaseManagerDeriv(
    Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager)
:   PhaseManagerDecorator(phasemanager),
    pressurederiv_(Teuchos::null),
    saturationderiv_(Teuchos::null),
    saturationderivderiv_(Teuchos::null),
    solidpressurederiv_(Teuchos::null),
    solidpressurederivderiv_(Teuchos::null)
{
  const int numphases = phasemanager_->NumPhases();
  // initialize matrixes and vectors
  pressurederiv_= Teuchos::rcp(new Epetra_SerialDenseMatrix(numphases,numphases));
  saturationderiv_= Teuchos::rcp(new Epetra_SerialDenseMatrix(numphases,numphases));
  saturationderivderiv_ = Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix>(numphases,Epetra_SerialDenseMatrix(numphases,numphases)));
  solidpressurederiv_= Teuchos::rcp(new Epetra_SerialDenseVector(numphases));
  solidpressurederivderiv_= Teuchos::rcp(new Epetra_SerialDenseMatrix(numphases,numphases));

  return;
}

/*----------------------------------------------------------------------*
 | evaluate pressures, saturations and derivatives at GP     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDeriv::EvaluateGPState(
    double                 J,
    const VariableManagerMinAccess& varmanager)
{
  // evaluate wrapped phase manager
  phasemanager_->EvaluateGPState(J,varmanager);

  // get material
  dsassert(phasemanager_->Element()->Material()!=Teuchos::null,"Material of element is null pointer!");
  const MAT::Material& material = *(phasemanager_->Element()->Material());

  // get number of phases
  const int numphases = phasemanager_->NumPhases();
  // get pressure
  const std::vector<double>& pressure = phasemanager_->Pressure();
  // get saturation
   const std::vector<double>& saturation = phasemanager_->Saturation();

  // access state vector
  const std::vector<double>& phinp = *varmanager.Phinp();

  //cast
  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  // calculate the derivative of the pressure (actually first its inverse)
  multiphasemat.EvaluateDerivOfDofWrtPressure(*pressurederiv_,phinp);

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
  Epetra_SerialDenseMatrix deriv(numphases,numphases);
  multiphasemat.EvaluateDerivOfSaturationWrtPressure(deriv,pressure);

  // chain rule: the derivative of saturation w.r.t. dof =
  // (derivative of saturation w.r.t. pressure) * (derivative of pressure w.r.t. dof)
  saturationderiv_->Multiply('N','N',1.0,deriv,*pressurederiv_,0.0);

  // calculate 2nd derivatives of saturation w.r.t. pressure
  // TODO: this should work for pressure und diffpressure DOFs, however not for
  //       saturation DOFs
  Teuchos::RCP<std::vector<Epetra_SerialDenseMatrix> > dummyderiv =
      Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix>(numphases,Epetra_SerialDenseMatrix(numphases,numphases)));
  multiphasemat.EvaluateSecondDerivOfSaturationWrtPressure(*dummyderiv,pressure);
  for (int i = 0; i < numphases; i++)
  {
    deriv.Multiply('T','N',1.0,*pressurederiv_,(*dummyderiv)[i],0.0);
    (*saturationderivderiv_)[i].Multiply('N','N',1.0,deriv,*pressurederiv_,0.0);
  }

  // compute derivative of solid pressure w.r.t. dofs with product rule
  solidpressurederiv_->Scale(0.0);
  for(int iphase=0; iphase<numphases; iphase++)
    for(int jphase=0; jphase<numphases; jphase++)
      (*solidpressurederiv_)(iphase) +=   (*pressurederiv_)(jphase,iphase)*saturation[jphase]
                                        + (*saturationderiv_)(jphase,iphase)*pressure[jphase];

  // compute second derivative of solid pressure w.r.t. dofs with product rule
  //TODO also include second derivs of pressure and saturation
  solidpressurederivderiv_->Scale(0.0);
  for(int iphase=0; iphase<numphases; iphase++)
    for(int jphase=0; jphase<numphases; jphase++)
      for(int kphase=0; kphase<numphases; kphase++)
        (*solidpressurederivderiv_)(jphase,kphase) +=
              (*pressurederiv_)(iphase,kphase)*(*saturationderiv_)(iphase,jphase)
            + (*saturationderiv_)(iphase,kphase)*(*pressurederiv_)(iphase,jphase);

  return;
}

/*----------------------------------------------------------------------*
 | zero all values at GP                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDeriv::ClearGPState()
{
  // this is just a safety call. All quantities should be recomputed
  // in the next EvaluateGPState call anyway, but you never know ...

  const int numphases = phasemanager_->NumPhases();
  phasemanager_->ClearGPState();

  //zero everything
  pressurederiv_->Scale(0.0);
  saturationderiv_->Scale(0.0);
  for (int iphase = 0; iphase < numphases; iphase++)
    (*saturationderivderiv_)[iphase].Scale(0.0);
  solidpressurederiv_->Scale(0.0);
  solidpressurederivderiv_->Scale(0.0);


  return;
}

/*----------------------------------------------------------------------*
 *  get derivative of saturation of phase                   vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDeriv::SaturationDeriv(
    int phasenum,
    int doftoderive) const
{
  phasemanager_->CheckIsEvaluated();

  return (*saturationderiv_)(phasenum,doftoderive);
}

/*----------------------------------------------------------------------*
 *  get 2nd derivative of saturation of phase          kremheller 05/17 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDeriv::SaturationDerivDeriv(
    int phasenum,
    int firstdoftoderive,
    int seconddoftoderive) const
{
  phasemanager_->CheckIsEvaluated();

  return (*saturationderivderiv_)[phasenum](firstdoftoderive,seconddoftoderive);
}

/*----------------------------------------------------------------------*
 *  get derivative of pressure of phase                      vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDeriv::PressureDeriv(
    int phasenum,
    int doftoderive) const
{
  phasemanager_->CheckIsEvaluated();

  return (*pressurederiv_)(phasenum,doftoderive);
}

/*----------------------------------------------------------------------*
 * get derivative of solid pressure                          vuong 08/16 |
*-----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDeriv::SolidPressureDeriv(
    int doftoderive) const
{
  phasemanager_->CheckIsEvaluated();

  return (*solidpressurederiv_)(doftoderive);
}

/*---------------------------------------------------------------------------*
 * get second derivative of solid pressure                       vuong 08/16 |
*---------------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDeriv::SolidPressureDerivDeriv(
    int doftoderive,
    int doftoderive2) const
{
  phasemanager_->CheckIsEvaluated();

  return (*solidpressurederivderiv_)(doftoderive,doftoderive2);
}


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::PhaseManagerDerivAndPorosity(
    Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager)
:   PhaseManagerDeriv(phasemanager),
    porosity_(0.0),
    J_(0.0),
    dporosity_dJ_(0.0),
    porosityderiv_(Teuchos::null)
{
  const int numphases = phasemanager_->NumPhases();
  // initialize matrixes and vectors
  porosityderiv_= Teuchos::rcp(new Epetra_SerialDenseVector(numphases));

  return;
}

/*----------------------------------------------------------------------*
 | evaluate pressures, saturations and derivatives at GP     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::EvaluateGPState(
    double                 J,
    const VariableManagerMinAccess& varmanager)
{
  // evaluate base class
  PhaseManagerDeriv::EvaluateGPState(J,varmanager);

  //access second material in structure element
  dsassert(phasemanager_->Element()->NumMaterial() > 1,"no second material defined for element ");
  dsassert(phasemanager_->Element()->Material(1)->MaterialType() == INPAR::MAT::m_structporo,"invalid structure material for poroelasticity");

  // cast second material to poro material
  Teuchos::RCP<MAT::StructPoro> structmat
  = Teuchos::rcp_static_cast<MAT::StructPoro>(phasemanager_->Element()->Material(1));

  // derivative of porosity w.r.t. solid pressure
  double dporosity_dp=0.0;

  //empty parameter list
  Teuchos::ParameterList params;

  J_ = J;

  //use structure material to evaluate porosity
  structmat->ComputePorosity( params,
                              phasemanager_->SolidPressure(),
                              J,
                              -1,
                              porosity_,
                              &dporosity_dp,
                              &dporosity_dJ_,
                              NULL,
                              NULL,
                              NULL,
                              false);

  // calculate the derivative of the porosity
  for(int iphase=0; iphase<phasemanager_->NumPhases(); iphase++)
      (*porosityderiv_)(iphase)   = dporosity_dp*SolidPressureDeriv(iphase);

  return;
}

/*----------------------------------------------------------------------*
 | zero all values at GP                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::ClearGPState()
{
  // this is just a safety call. All quantities should be recomputed
  // in the next EvaluateGPState call anyway, but you never know ...

  PhaseManagerDeriv::ClearGPState();

  //zero everything
  porosity_=0.0;
  J_ = 0.0;
  dporosity_dJ_ = 0.0;
  porosityderiv_->Scale(0.0);

  return;
}

/*----------------------------------------------------------------------*
 *  get porosity                                            vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::Porosity() const
{
  CheckIsEvaluated();

  return  porosity_;
}

/*----------------------------------------------------------------------*
 *  get Jacobian                                       kremheller 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::Jacobian() const
{
  CheckIsEvaluated();

  return J_;
}

/*----------------------------------------------------------------------*
 *  get derivative of porosity wrt jacobian            kremheller 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::PorosityDerivWrtJacobian() const
{
  CheckIsEvaluated();

  return dporosity_dJ_;
}

/*----------------------------------------------------------------------*
 *  get derivative of saturation of phase                   vuong 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::PorosityDeriv(
    int doftoderive) const
{
  phasemanager_->CheckIsEvaluated();

  return (*porosityderiv_)(doftoderive);
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::PhaseManagerReaction(
    Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager
    )
: PhaseManagerDecorator(phasemanager)
{
  return;
}

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::Setup(
    const DRT::Element* ele
    )
{
  // setup the wrapped class
  phasemanager_->Setup(ele);

  // get material
  const MAT::Material& material = *(phasemanager_->Element()->Material());

  // get number of phases
  const int numphases = phasemanager_->NumPhases();

  isreactive_.resize(numphases,false);

  // cast the material to multiphase material
  const MAT::FluidPoroMultiPhaseReactions& multiphasemat =
      dynamic_cast<const MAT::FluidPoroMultiPhaseReactions&>(material);

  // get number of reactions
  const int numreactions = multiphasemat.NumReac();

  for(int ireac=0; ireac<numreactions; ireac++)
  {
    // get the single phase material
    MAT::FluidPoroSingleReaction& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSingleReactionMatFromMultiReactionsMaterial(multiphasemat,ireac);

    for(int iphase=0; iphase<numphases; iphase++)
    {
      isreactive_[iphase] = isreactive_[iphase] or singlephasemat.IsReactive(iphase);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::EvaluateGPState(
    double                 J,
    const VariableManagerMinAccess& varmanager
    )
{
  phasemanager_->EvaluateGPState(J,varmanager);

  // get material
  dsassert(phasemanager_->Element()->Material()!=Teuchos::null,"Material of element is null pointer!");
  const MAT::Material& material = *(phasemanager_->Element()->Material());

  if(material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("Invalid material! Only MAT_FluidPoroMultiPhaseReactions material valid for reaction evaluation!");

  // cast the material to multiphase material
  const MAT::FluidPoroMultiPhaseReactions& multiphasemat =
      dynamic_cast<const MAT::FluidPoroMultiPhaseReactions&>(material);

  // get number of reactions
  const int numreactions = multiphasemat.NumReac();
  // get number of phases
  const int numphases = phasemanager_->NumPhases();

  reacterms_.clear();
  reactermsderivs_.clear();
  reactermsderivspressure_.clear();
  reactermsderivssaturation_.clear();
  reactermsderivsporosity_.clear();
  reacterms_.resize(numphases,0.0);
  reactermsderivs_.resize(numphases,std::vector<double>(numphases,0.0));
  reactermsderivspressure_.resize(numphases,std::vector<double>(numphases,0.0));
  reactermsderivssaturation_.resize(numphases,std::vector<double>(numphases,0.0));
  reactermsderivsporosity_.resize(numphases,0.0);

  for(int ireac=0; ireac<numreactions; ireac++)
  {
    // get the single phase material
    MAT::FluidPoroSingleReaction& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSingleReactionMatFromMultiReactionsMaterial(multiphasemat,ireac);

    // evaluate the reaction
    singlephasemat.EvaluateReaction(
        reacterms_,
        reactermsderivspressure_,
        reactermsderivssaturation_,
        reactermsderivsporosity_,
        phasemanager_->Pressure(),
        phasemanager_->Saturation(),
        phasemanager_->Porosity(),
        *varmanager.Scalarnp()
        );
  }

  for(int iphase=0; iphase<numphases; iphase++)
  {
    std::vector<double>& myphasederiv = reactermsderivs_[iphase];
    const std::vector<double>& myderivspressure = reactermsderivspressure_[iphase];
    const std::vector<double>& myderivssaturation = reactermsderivssaturation_[iphase];

    for(int doftoderive=0; doftoderive<numphases; doftoderive++)
    {
      for(int idof=0; idof<numphases; idof++)
        myphasederiv[doftoderive] +=   myderivspressure[idof]*phasemanager_->PressureDeriv(idof,doftoderive)
                                     + myderivssaturation[idof]*phasemanager_->SaturationDeriv(idof,doftoderive);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | zero all values at GP                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::ClearGPState()
{
  // this is just a safety call. All quantities should be recomputed
  // in the next EvaluateGPState call anyway, but you never know ...
  phasemanager_->ClearGPState();

  reacterms_.clear();
  reactermsderivs_.clear();
  reactermsderivspressure_.clear();
  reactermsderivssaturation_.clear();
  reactermsderivsporosity_.clear();

  return;
}

/*----------------------------------------------------------------------*
 | get the reaction term                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::ReacTerm(int phasenum) const
{
  phasemanager_->CheckIsEvaluated();
  return reacterms_[phasenum];
}

/*----------------------------------------------------------------------*
 | get the derivative of the reaction term                   vuong 08/16 |
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::ReacDeriv(int phasenum, int doftoderive) const
{
  phasemanager_->CheckIsEvaluated();
  return reactermsderivs_[phasenum][doftoderive];
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
template<int nsd>
DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::PhaseManagerDiffusion(
    Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager
    )
: PhaseManagerDecorator(phasemanager)
{
  return;
}

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
template<int nsd>
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::Setup(
    const DRT::Element* ele
    )
{
  // setup the wrapped class
  phasemanager_->Setup(ele);

  // get number of phases
  const int numphases = phasemanager_->NumPhases();

  // resize vectors to numphase
  constrelpermeability_.resize(numphases,false);
  constdynviscosity_.resize(numphases,false);
  relpermeabilities_.resize(numphases,false);
  derrelpermeabilities_.resize(numphases,false);

  const MAT::Material& material = *( phasemanager_->Element()->Material());

  // check material type
  if (material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase
      and material.MaterialType()
          != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("only poro multiphase material valid");

  for(int iphase=0; iphase<numphases; iphase++)
  {
    //get the single phase material
    const MAT::FluidPoroSinglePhase& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,iphase);
    constrelpermeability_[iphase]=singlephasemat.HasConstantRelPermeability();
    constdynviscosity_[iphase]=singlephasemat.HasConstantViscosity();
  }

  return;
}
/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
template<int nsd>
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::EvaluateGPState(
    double                 J,
    const VariableManagerMinAccess& varmanager
    )
{
  // evaluate wrapped manager
  phasemanager_->EvaluateGPState( J, varmanager);

  // get material
  dsassert(phasemanager_->Element()->Material()!=Teuchos::null,"Material of element is null pointer!");
  const MAT::Material& material = *( phasemanager_->Element()->Material());

  // get number of phases
  const int numphases = phasemanager_->NumPhases();

  // resize vectors
  permeabilitytensors_.resize(numphases);
  relpermeabilities_.resize(numphases);
  derrelpermeabilities_.resize(numphases);

  // check material type
  if (material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase
      and material.MaterialType()
          != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("only poro multiphase material valid");

    // cast to multiphase material
  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  for (int iphase = 0; iphase < numphases; iphase++)
  {
    // TODO only isotropic, constant permeability for now
    permeabilitytensors_[iphase].Clear();
    const double permeability = multiphasemat.Permeability();
    for (int i = 0; i < nsd; i++)
      (permeabilitytensors_[iphase])(i, i) = permeability;

    //get the single phase material
    const MAT::FluidPoroSinglePhase& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(multiphasemat,iphase);

    // evaluate relative permeabilities
    relpermeabilities_[iphase] =
        singlephasemat.RelPermeability(phasemanager_->Saturation(iphase));

    // evaluate derivative of relative permeability w.r.t. saturation
    derrelpermeabilities_[iphase] =
        singlephasemat.EvaluateDerivOfRelPermeabilityWrtSaturation(phasemanager_->Saturation(iphase));
  }

  return;
 }

/*----------------------------------------------------------------------*
 | zero all values at GP                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
template<int nsd>
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::ClearGPState()
{
  // this is just a safety call. All quantities should be recomputed
  // in the next EvaluateGPState call anyway, but you never know ...
  phasemanager_->ClearGPState();

  permeabilitytensors_.clear();
  relpermeabilities_.clear();
  derrelpermeabilities_.clear();

  return;
}

/*---------------------------------------------------------------------------*
 * get diffusion tensor                                         vuong 08/16 |
*---------------------------------------------------------------------------*/
template<int nsd>
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::PermeabilityTensor(
    int phasenum, LINALG::Matrix<nsd,nsd>& permeabilitytensor) const
{
  phasemanager_->CheckIsEvaluated();
  // make a hard copy for now
  permeabilitytensor = permeabilitytensors_[phasenum];
}

/*---------------------------------------------------------------------------*
 * check for constant rel permeability                      kremheller 02/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
bool DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::HasConstantRelPermeability(
    int phasenum) const
{
  CheckIsSetup();

  return constrelpermeability_[phasenum];
}
/*---------------------------------------------------------------------------*
 * get relative diffusivity of phase 'phasenum'             kremheller 02/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::RelPermeability(
    int phasenum) const
{
  phasemanager_->CheckIsEvaluated();

  return relpermeabilities_[phasenum];
}
/*---------------------------------------------------------------------------*
 * get relative diffusivity of phase 'phasenum'             kremheller 02/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::RelPermeabilityDeriv(
    int phasenum) const
{
  phasemanager_->CheckIsEvaluated();

  return derrelpermeabilities_[phasenum];
}
/*---------------------------------------------------------------------------*
 * check for constant dynamic viscosity                     kremheller 06/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
bool DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::HasConstantDynViscosity(
    int phasenum) const
{
  CheckIsSetup();

  return constdynviscosity_[phasenum];
}
/*---------------------------------------------------------------------------*
 * get dynamic viscosity of phase 'phasenum'                kremheller 02/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::DynViscosity(
    int phasenum,
    double abspressgrad) const
{
  phasemanager_->CheckIsEvaluated();

  return DynViscosity(*Element()->Material(),phasenum,abspressgrad);
}

/*----------------------------------------------------------------------*
 *   get dynamic viscosity of phase 'phasenum'         kremheller 02/17 |
*----------------------------------------------------------------------*/
template<int nsd>
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::DynViscosity(
    const MAT::Material& material,
    int phasenum,
    double abspressgrad) const
{

  //get the single phase material
  const MAT::FluidPoroSinglePhase& singlephasemat =
      POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,phasenum);

  return singlephasemat.Viscosity(abspressgrad);
}

/*---------------------------------------------------------------------------*
 * get derivative of dynamic viscosity of phase 'phasenum'  kremheller 06/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::DynViscosityDeriv(
    int phasenum,
    double abspressgrad) const
{
  phasemanager_->CheckIsEvaluated();

  return DynViscosityDeriv(*Element()->Material(),phasenum,abspressgrad);
}

/*---------------------------------------------------------------------------*
 * get derivative of dynamic viscosity of phase 'phasenum'  kremheller 06/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<nsd>::DynViscosityDeriv(
    const MAT::Material& material,
    int phasenum,
    double abspressgrad) const
{

  //get the single phase material
  const MAT::FluidPoroSinglePhase& singlephasemat =
      POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,phasenum);

  return singlephasemat.ViscosityDeriv(abspressgrad);
}

///*----------------------------------------------------------------------*
// *----------------------------------------------------------------------*/
//// template classes

template class DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<1>;
template class DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<2>;
template class DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<3>;

