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
#include "../drt_mat/scatra_mat_multiporo.H"
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
    int totalnumdofpernode,
    int numfluidphases)
{
  Teuchos::RCP< PhaseManagerInterface > phasemanager = Teuchos::null;

  // build the standard phase manager
  phasemanager = Teuchos::rcp(new PhaseManagerCore(totalnumdofpernode,numfluidphases));

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
  {
    // no extensions needed
    phasemanager = corephasemanager;
    break;
  }
  // calculate solid pressure
  case POROFLUIDMULTIPHASE::calc_solidpressure:
  {
    // we have volume fractions --> we need PhaseManagerDerivAndPorosity because solid pressure is calculated
    // differently
    if(corephasemanager->TotalNumDof() > corephasemanager->NumFluidPhases())
    {
      // porosity (includes derivatves) needed
      phasemanager = Teuchos::rcp(new PhaseManagerDerivAndPorosity(corephasemanager));

      if(mattype==INPAR::MAT::m_fluidporo_multiphase_reactions)
        // enhance by scalar handling capability
        phasemanager = Teuchos::rcp(new PhaseManagerReaction(phasemanager));

      switch(nsd)
      {
      case 1:
        phasemanager = Teuchos::rcp(new PhaseManagerVolFrac<1>(phasemanager));
        break;
      case 2:
        phasemanager = Teuchos::rcp(new PhaseManagerVolFrac<2>(phasemanager));
        break;
      case 3:
        phasemanager = Teuchos::rcp(new PhaseManagerVolFrac<3>(phasemanager));
        break;
      default:
        dserror("invalid dimension for creating phase manager!");
      }
    }
    else
    {
      // no extensions needed
      phasemanager = corephasemanager;
    }
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
  case POROFLUIDMULTIPHASE::calc_fluid_struct_coupl_mat:
  case POROFLUIDMULTIPHASE::calc_fluid_scatra_coupl_mat:
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

    if(corephasemanager->TotalNumDof() > corephasemanager->NumFluidPhases())
    {
      switch(nsd)
      {
      case 1:
        phasemanager = Teuchos::rcp(new PhaseManagerVolFrac<1>(phasemanager));
        break;
      case 2:
        phasemanager = Teuchos::rcp(new PhaseManagerVolFrac<2>(phasemanager));
        break;
      case 3:
        phasemanager = Teuchos::rcp(new PhaseManagerVolFrac<3>(phasemanager));
        break;
      default:
        dserror("invalid dimension for creating phase manager!");
      }
    }

    break;
  }
  case POROFLUIDMULTIPHASE::get_access_from_scatra:
  {
    // porosity (includes derivatives) needed
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
DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::PhaseManagerCore(int totalnumdofpernode, int numfluidphases):
    totalnumdofpernode_(totalnumdofpernode),
    numfluidphases_(numfluidphases),
    genpressure_(numfluidphases,0.0),
    volfrac_(totalnumdofpernode-numfluidphases,0.0),
    sumaddvolfrac_(0.0),
    pressure_(numfluidphases,0.0),
    saturation_(numfluidphases,0.0),
    density_(numfluidphases,0.0),
    solidpressure_(0.0),
    invbulkmodulifluid_(numfluidphases,0.0),
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
    totalnumdofpernode_(old.totalnumdofpernode_),
    numfluidphases_(old.numfluidphases_),
    genpressure_(old.genpressure_),
    volfrac_(old.volfrac_),
    sumaddvolfrac_(old.sumaddvolfrac_),
    pressure_(old.pressure_),
    saturation_(old.saturation_),
    density_(old.density_),
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
    const DRT::Element* ele,
    const int matnum)
{

  dsassert(ele!=NULL,"Element is null pointer for setup of phase manager!");
  // save current element
  ele_ = ele;
  // get material
  const MAT::Material& material = *(ele_->Material(matnum));

  // check the material
  if(material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase and
     material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("only poro multiphase and poro multiphase reactions material valid");

  // cast
  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  //safety check
  if(numfluidphases_ != multiphasemat.NumFluidPhases() ||
      totalnumdofpernode_ != multiphasemat.NumMat())
    dserror("Number of phases given by the poro multiphase material does not match number "
        "of DOFs (%i phases and %i Fluid DOFs, %i phases and %i Total DOFs)!"
        , numfluidphases_, multiphasemat.NumFluidPhases(), totalnumdofpernode_, multiphasemat.NumMat());

  issetup_=true;
  return;
}

/*----------------------------------------------------------------------*
 | evaluate pressures, saturations and derivatives at GP     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::EvaluateGPState(
    double                 J,
    const VariableManagerMinAccess& varmanager,
    const int matnum)
{
  CheckIsSetup();

  if(isevaluated_ == true)
    dserror("state has already been set!");

  // get material
  dsassert(ele_->Material(matnum)!=Teuchos::null,"Material of element is null pointer!");
  const MAT::Material& material = *(ele_->Material(matnum));

  // access state vector
  const std::vector<double>& phinp = *varmanager.Phinp();

  if(totalnumdofpernode_ != (int)phinp.size())
    dserror("Length of phinp vector is not equal to the number of dofs");

  // fluid primary variables
  const std::vector<double> fluid_phinp(&phinp[0],&phinp[numfluidphases_]);
  // volume fractions
  const std::vector<double> volfrac(&(phinp[numfluidphases_]),&(phinp[totalnumdofpernode_]));

  if(numfluidphases_ != (int)fluid_phinp.size())
    dserror("Length of fluid-phinp vector is not equal to the number of phases");

  // cast the material to multiphase material
  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  //resize all vectors
  genpressure_.resize(numfluidphases_,0.0);
  volfrac_.resize(totalnumdofpernode_-numfluidphases_,0.0);
  pressure_.resize(numfluidphases_,0.0);
  saturation_.resize(numfluidphases_,0.0);
  density_.resize(numfluidphases_,0.0);
  invbulkmodulifluid_.resize(numfluidphases_,0.0);

  // evaluate the volume fractions
  volfrac_ = volfrac;
  if(totalnumdofpernode_-numfluidphases_ != (int)volfrac.size())
    dserror("Length of volfrac vector is not equal to the number of volume fractions");
  sumaddvolfrac_ = 0.0;
  for (int ivolfrac = 0; ivolfrac < totalnumdofpernode_-numfluidphases_; ivolfrac++)
    sumaddvolfrac_ += volfrac_[ivolfrac];

  // evaluate the pressures
  multiphasemat.EvaluateGenPressure(genpressure_,fluid_phinp);

  //! transform generalized pressures to true pressure values
  multiphasemat.TransformGenPresToTruePres(genpressure_,pressure_);

  // explicit evaluation of saturation
  multiphasemat.EvaluateSaturation(saturation_,fluid_phinp,pressure_);

  // solid pressure = sum (S_i*p_i)
  solidpressure_ = std::inner_product(saturation_.begin(),saturation_.end(),pressure_.begin(),0.0);

  //access second material in structure element
  dsassert(ele_->NumMaterial() > 1,"no second material defined for element ");
  dsassert(ele_->Material(1)->MaterialType() == INPAR::MAT::m_structporo,"invalid structure material for poroelasticity");

  // cast second material to poro material
  Teuchos::RCP<MAT::StructPoro> structmat = Teuchos::rcp_static_cast<MAT::StructPoro>(ele_->Material(1));

  invbulkmodulussolid_ = structmat->InvBulkmodulus();

  for(int iphase=0; iphase<numfluidphases_; iphase++)
  {
    //get the single phase material
    const MAT::FluidPoroSinglePhase& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(material,iphase);
    invbulkmodulifluid_[iphase]=singlephasemat.InvBulkmodulus();

    //get density
    density_[iphase] = singlephasemat.Density();
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
  std::fill(volfrac_.begin(), volfrac_.end(), 0.0);
  std::fill(pressure_.begin(), pressure_.end(), 0.0);
  std::fill(saturation_.begin(), saturation_.end(), 0.0);
  std::fill(density_.begin(), density_.end(), 0.0);
  solidpressure_=0.0;
  sumaddvolfrac_=0.0;
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
 * recalculate solid pressure if volfracs are present  kremheller 09/17 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::RecalculateSolidPressure(
    const double porosity,
    const std::vector<double>& volfracpressure
    )
{
  CheckIsEvaluated();

  // p_s = (porosity)/(porosity+sumaddvolfrac) * sum_i^numfluidphases (S_i*p_i)
  //     + 1.0 / (porosity+sumaddvolfrac) sum_i=1^numvolfrac (volfrac_i*pressure_i)
  // with porosity = porosity' - sumaddvolfrac
  // porosity' := porosity if no volume fractions were present
  solidpressure_*= porosity/(sumaddvolfrac_+porosity);
  for(int ivolfrac = 0; ivolfrac < totalnumdofpernode_-numfluidphases_; ivolfrac++)
    solidpressure_ += volfrac_[ivolfrac]/(porosity+sumaddvolfrac_)*volfracpressure[ivolfrac];

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
 *    get volfracs                                     kremheller 08/17 |
*----------------------------------------------------------------------*/
const std::vector<double>& DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::VolFrac() const
{
  CheckIsEvaluated();

  return volfrac_;
}

/*----------------------------------------------------------------------*
 *    get volfrac of volfrac 'volfracnum'              kremheller 08/17 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::VolFrac(int volfracnum) const
{
  CheckIsEvaluated();

  return volfrac_[volfracnum];
}

/*----------------------------------------------------------------------*
 * get the sum of the additional volume fractions      kremheller 09/17 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerCore::SumAddVolFrac() const
{
  CheckIsEvaluated();

  return sumaddvolfrac_;
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
    int phasenum) const
{
  CheckIsEvaluated();

  return density_[phasenum];
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
  const int numfluidphases = phasemanager_->NumFluidPhases();
  // initialize matrixes and vectors
  pressurederiv_= Teuchos::rcp(new Epetra_SerialDenseMatrix(numfluidphases,numfluidphases));
  saturationderiv_= Teuchos::rcp(new Epetra_SerialDenseMatrix(numfluidphases,numfluidphases));
  saturationderivderiv_ = Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix>(numfluidphases,Epetra_SerialDenseMatrix(numfluidphases,numfluidphases)));
  solidpressurederiv_= Teuchos::rcp(new Epetra_SerialDenseVector(numfluidphases));
  solidpressurederivderiv_= Teuchos::rcp(new Epetra_SerialDenseMatrix(numfluidphases,numfluidphases));

  return;
}

/*----------------------------------------------------------------------*
 | evaluate pressures, saturations and derivatives at GP     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDeriv::EvaluateGPState(
    double                 J,
    const VariableManagerMinAccess& varmanager,
    const int matnum)
{
  // evaluate wrapped phase manager
  phasemanager_->EvaluateGPState(J,varmanager,matnum);

  // get material
  dsassert(phasemanager_->Element()->Material(matnum)!=Teuchos::null,"Material of element is null pointer!");
  const MAT::Material& material = *(phasemanager_->Element()->Material(matnum));

  // get number of phases
  const int numfluidphases = phasemanager_->NumFluidPhases();
  // get pressure
  const std::vector<double>& pressure = phasemanager_->Pressure();
  // get saturation
  const std::vector<double>& saturation = phasemanager_->Saturation();

  // access state vector
  const std::vector<double>& phinp = *varmanager.Phinp();

  if(phasemanager_->TotalNumDof()!= (int)phinp.size())
    dserror("Length of phinp vector is not equal to the number of dofs");

  // get fluid primary variable
  const std::vector<double> fluid_phinp(&phinp[0],&phinp[numfluidphases]);

  if(numfluidphases != (int)fluid_phinp.size())
    dserror("Length of phinp vector is not equal to the number of phases");

  //cast
  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  // calculate the derivative of the pressure (actually first its inverse)
  multiphasemat.EvaluateDerivOfDofWrtPressure(*pressurederiv_,fluid_phinp);

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
  Epetra_SerialDenseMatrix deriv(numfluidphases,numfluidphases);
  multiphasemat.EvaluateDerivOfSaturationWrtPressure(deriv,pressure);

  // chain rule: the derivative of saturation w.r.t. dof =
  // (derivative of saturation w.r.t. pressure) * (derivative of pressure w.r.t. dof)
  saturationderiv_->Multiply('N','N',1.0,deriv,*pressurederiv_,0.0);

  // calculate 2nd derivatives of saturation w.r.t. pressure
  // TODO: this should work for pressure und diffpressure DOFs, however not for
  //       saturation DOFs
  Teuchos::RCP<std::vector<Epetra_SerialDenseMatrix> > dummyderiv =
      Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix>(numfluidphases,Epetra_SerialDenseMatrix(numfluidphases,numfluidphases)));
  multiphasemat.EvaluateSecondDerivOfSaturationWrtPressure(*dummyderiv,pressure);
  for (int i = 0; i < numfluidphases; i++)
  {
    deriv.Multiply('T','N',1.0,*pressurederiv_,(*dummyderiv)[i],0.0);
    (*saturationderivderiv_)[i].Multiply('N','N',1.0,deriv,*pressurederiv_,0.0);
  }

  // compute derivative of solid pressure w.r.t. dofs with product rule
  solidpressurederiv_->Scale(0.0);
  for(int iphase=0; iphase<numfluidphases; iphase++)
    for(int jphase=0; jphase<numfluidphases; jphase++)
      (*solidpressurederiv_)(iphase) +=   (*pressurederiv_)(jphase,iphase)*saturation[jphase]
                                        + (*saturationderiv_)(jphase,iphase)*pressure[jphase];

  // compute second derivative of solid pressure w.r.t. dofs with product rule
  //TODO also include second derivs of pressure and saturation
  solidpressurederivderiv_->Scale(0.0);
  for(int iphase=0; iphase<numfluidphases; iphase++)
    for(int jphase=0; jphase<numfluidphases; jphase++)
      for(int kphase=0; kphase<numfluidphases; kphase++)
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

  const int numfluidphases = phasemanager_->NumFluidPhases();
  phasemanager_->ClearGPState();

  //zero everything
  pressurederiv_->Scale(0.0);
  saturationderiv_->Scale(0.0);
  for (int iphase = 0; iphase < numfluidphases; iphase++)
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
    dporosity_dp_(0.0),
    porosityderiv_(Teuchos::null)
{
  const int totalnumdof = phasemanager_->TotalNumDof();
  // initialize matrixes and vectors
  porosityderiv_= Teuchos::rcp(new Epetra_SerialDenseVector(totalnumdof));

  return;
}

/*----------------------------------------------------------------------*
 | evaluate pressures, saturations and derivatives at GP     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::EvaluateGPState(
    double                 J,
    const VariableManagerMinAccess& varmanager,
    const int matnum)
{
  // evaluate base class
  PhaseManagerDeriv::EvaluateGPState(J,varmanager,matnum);

  //access second material in structure element
  dsassert(phasemanager_->Element()->NumMaterial() > 1,"no second material defined for element ");
  dsassert(phasemanager_->Element()->Material(1)->MaterialType() == INPAR::MAT::m_structporo,"invalid structure material for poroelasticity");

  // cast second material to poro material
  Teuchos::RCP<MAT::StructPoro> structmat
  = Teuchos::rcp_static_cast<MAT::StructPoro>(phasemanager_->Element()->Material(1));

  //empty parameter list
  Teuchos::ParameterList params;

  J_ = J;

  //use structure material to evaluate porosity
  structmat->ComputePorosity( params,
                              phasemanager_->SolidPressure(),
                              J,
                              -1,
                              porosity_,
                              &dporosity_dp_,
                              &dporosity_dJ_,
                              NULL,
                              NULL,
                              NULL,
                              false);

  // Note:
  // for phase law density dependent, incompressible:
  // porosity = 1 - sumaddvolfrac - (1 - porosity_0 - sumaddvolfrac_0)/J
  // with porosity_0 = porosity_0' - sumaddvolfrac_0   (porosity_0' = element initial porosity)
  // porosity = 1 - sumaddvolfrac - (1 - porosity_0')/J

  // subtract additional volume fractions
  // porosity = porosity' - sumaddvolfracs
  porosity_ -= phasemanager_->SumAddVolFrac();

  // calculate the derivative of the porosity w.r.t. fluid phases
  for(int iphase=0; iphase<phasemanager_->NumFluidPhases(); iphase++)
      (*porosityderiv_)(iphase)   = dporosity_dp_*SolidPressureDeriv(iphase);

  // additional derivatives w.r.t. volume fractions
  for(int ivolfrac=0; ivolfrac<phasemanager_->TotalNumDof()-phasemanager_->NumFluidPhases(); ivolfrac++)
    (*porosityderiv_)(ivolfrac+phasemanager_->NumFluidPhases()) = - 1.0;

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
  dporosity_dp_ = 0.0;
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
 *  get Jacobian of deformation gradient               kremheller 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::JacobianDefGrad() const
{
  CheckIsEvaluated();

  return J_;
}

/*----------------------------------------------------------------------*
 *  get derivative of porosity wrt jacobian of defgrad kremheller 08/16 |
*----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::PorosityDerivWrtJacobianDefGrad() const
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
 * check if porosity depends on pressure               kremheller 07/17 |
*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::PorosityDependsOnFluid() const
{
  phasemanager_->CheckIsEvaluated();

  // this check is needed for speeding up calculations of element stiffness matrices
  return (fabs(dporosity_dp_) > 1.0e-10 || phasemanager_->TotalNumDof()-phasemanager_->NumFluidPhases() > 0);
}

/*----------------------------------------------------------------------*
 * check if porosity depends on JacobianDefGradient    kremheller 07/17 |
*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDerivAndPorosity::PorosityDependsOnStruct() const
{
  phasemanager_->CheckIsEvaluated();

  // this check is needed for speeding up calculations of element stiffness matrices
  return (fabs(dporosity_dJ_) > 1.0e-10);
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
: PhaseManagerDecorator(phasemanager),
  numscal_(0)
{
  return;
}

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::Setup(
    const DRT::Element* ele,
    const int matnum
    )
{
  // setup the wrapped class
  phasemanager_->Setup(ele,matnum);

  // get material
  const MAT::Material& material = *(phasemanager_->Element()->Material(matnum));

  // get total number of dofs
  const int totalnumdof = phasemanager_->TotalNumDof();

  isreactive_.resize(totalnumdof,false);

  // cast the material to multiphase material
  const MAT::FluidPoroMultiPhaseReactions& multiphasemat =
      dynamic_cast<const MAT::FluidPoroMultiPhaseReactions&>(material);

  // get number of reactions
  const int numreactions = multiphasemat.NumReac();

  // cast third material to scatra material
  // TODO: does this always work? probably yes if called with porofluidmultiphase element, but not when
  //       access from outside is requested with matnum =! 0, however, in that case we will not need the
  //       PhaseManagerReaction
  Teuchos::RCP<MAT::MatList> scatramat
  = Teuchos::rcp_static_cast<MAT::MatList>(phasemanager_->Element()->Material(2));

  // get number of scalars
  numscal_ = scatramat->NumMat();
  scalartophaseid_.resize(numscal_);

  // fill scalar to phase id vector
  if (scatramat->MaterialType() == INPAR::MAT::m_matlist or
      scatramat->MaterialType() == INPAR::MAT::m_matlist_reactions)
  {
    for (int k = 0;k<numscal_;++k)
    {
      int matid = scatramat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = scatramat->MaterialById(matid);
      if(singlemat->MaterialType() == INPAR::MAT::m_scatra_multiporo)
      {
        const Teuchos::RCP<const MAT::ScatraMatMultiPoro>& poromat
          = Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoro>(singlemat);
        scalartophaseid_[k] = poromat->PhaseID();
      }
      else
        dserror("only MAT_scatra_multiporo valid here");
    }
  }
  else if(scatramat->MaterialType() == INPAR::MAT::m_scatra_multiporo)
  {
    const Teuchos::RCP<const MAT::ScatraMatMultiPoro>& poromat
      = Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoro>(scatramat);
    scalartophaseid_[0] = poromat->PhaseID();
  }
  else
    dserror("only MAT_matlist_reactions, MAT_matlist or MAT_scatra_multiporo valid here");

  for(int ireac=0; ireac<numreactions; ireac++)
  {
    // get the single phase material
    MAT::FluidPoroSingleReaction& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSingleReactionMatFromMultiReactionsMaterial(multiphasemat,ireac);

    for(int iphase=0; iphase<totalnumdof; iphase++)
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
    const VariableManagerMinAccess& varmanager,
    const int matnum
    )
{
  phasemanager_->EvaluateGPState(J,varmanager,matnum);

  // get material
  dsassert(phasemanager_->Element()->Material(matnum)!=Teuchos::null,"Material of element is null pointer!");
  const MAT::Material& material = *(phasemanager_->Element()->Material(matnum));

  if(material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("Invalid material! Only MAT_FluidPoroMultiPhaseReactions material valid for reaction evaluation!");

  // cast the material to multiphase material
  const MAT::FluidPoroMultiPhaseReactions& multiphasemat =
      dynamic_cast<const MAT::FluidPoroMultiPhaseReactions&>(material);

  // get number of reactions
  const int numreactions = multiphasemat.NumReac();
  // get number of phases
  const int numfluidphases = phasemanager_->NumFluidPhases();
  // get total number of dofs
  const int totalnumdof = phasemanager_->TotalNumDof();

  // get the volume fraction vector
  const std::vector<double> volfrac = phasemanager_->VolFrac();

  reacterms_.clear();
  reactermsderivs_.clear();
  reactermsderivspressure_.clear();
  reactermsderivssaturation_.clear();
  reactermsderivsporosity_.clear();
  reactermsderivsscalar_.clear();
  reacterms_.resize(totalnumdof,0.0);
  reactermsderivs_.resize(totalnumdof,std::vector<double>(totalnumdof,0.0));
  reactermsderivspressure_.resize(totalnumdof,std::vector<double>(numfluidphases,0.0));
  reactermsderivssaturation_.resize(totalnumdof,std::vector<double>(numfluidphases,0.0));
  reactermsderivsporosity_.resize(totalnumdof,0.0);
  reactermsderivsscalar_.resize(totalnumdof,std::vector<double>(numscal_,0.0));
  reactermsderivsvolfrac_.resize(totalnumdof,std::vector<double>(totalnumdof-numfluidphases,0.0));

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
        reactermsderivsvolfrac_,
        reactermsderivsscalar_,
        phasemanager_->Pressure(),
        phasemanager_->Saturation(),
        phasemanager_->Porosity(),
        volfrac,
        *varmanager.Scalarnp()
        );
  }

  for(int jdof=0; jdof<totalnumdof; jdof++)
  {
    std::vector<double>& myphasederiv = reactermsderivs_[jdof];
    const std::vector<double>& myderivspressure = reactermsderivspressure_[jdof];
    const std::vector<double>& myderivssaturation = reactermsderivssaturation_[jdof];
    const std::vector<double>& myderivsvolfrac = reactermsderivsvolfrac_[jdof];
    const double myderivsporosity = reactermsderivsporosity_[jdof];

    // reaction derivs w.r.t. to fluid phases
    for(int doftoderive=0; doftoderive<numfluidphases; doftoderive++)
    {
      for(int idof=0; idof<numfluidphases; idof++)
        myphasederiv[doftoderive] +=   myderivspressure[idof]*phasemanager_->PressureDeriv(idof,doftoderive)
                                     + myderivssaturation[idof]*phasemanager_->SaturationDeriv(idof,doftoderive);
      if(phasemanager_->PorosityDependsOnFluid())
      {
        myphasederiv[doftoderive] += myderivsporosity*phasemanager_->PorosityDeriv(doftoderive);
      }
    }
    // reaction derivs w.r.t. to volume fraction phases
    for(int doftoderive=numfluidphases; doftoderive<totalnumdof; doftoderive++)
    {
      myphasederiv[doftoderive] += myderivsvolfrac[doftoderive-numfluidphases]
                                 + myderivsporosity*phasemanager_->PorosityDeriv(doftoderive);
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
  reactermsderivsscalar_.clear();
  reactermsderivsvolfrac_.clear();

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
 | get total number of scalars                         kremheller 07/17 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::NumScal() const
{
  phasemanager_->CheckIsSetup();
  return numscal_;
}

/*------------------------------------------------------------------------*
 | get the derivative of the reaction term wrt. porosity kremheller 07/17 |
 *------------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::ReacDerivPorosity(int phasenum) const
{
  phasemanager_->CheckIsEvaluated();
  return reactermsderivsporosity_[phasenum];
}

/*----------------------------------------------------------------------*
 | get the derivative of the reaction term wrt. scalar kremheller 07/17 |
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::ReacDerivScalar(int phasenum, int scaltoderive) const
{
  phasemanager_->CheckIsEvaluated();
  return reactermsderivsscalar_[phasenum][scaltoderive];
}

/*----------------------------------------------------------------------*
 | get the scalar to phase id                          kremheller 08/17 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerReaction::ScalarToPhaseID(int iscal) const
{
  phasemanager_->CheckIsEvaluated();
  return scalartophaseid_[iscal];
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
    const DRT::Element* ele,
    const int matnum
    )
{
  // setup the wrapped class
  phasemanager_->Setup(ele, matnum);

  // get number of phases
  const int numfluidphases = phasemanager_->NumFluidPhases();

  // resize vectors to numfluidphases
  constrelpermeability_.resize(numfluidphases,false);
  constdynviscosity_.resize(numfluidphases,false);
  relpermeabilities_.resize(numfluidphases,false);
  derrelpermeabilities_.resize(numfluidphases,false);

  const MAT::Material& material = *( phasemanager_->Element()->Material(matnum));

  // check material type
  if (material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase
      and material.MaterialType()
          != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("only poro multiphase material valid");

  for(int iphase=0; iphase<numfluidphases; iphase++)
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
    const VariableManagerMinAccess& varmanager,
    const int matnum
    )
{
  // evaluate wrapped manager
  phasemanager_->EvaluateGPState( J, varmanager, matnum);

  // get material
  dsassert(phasemanager_->Element()->Material(matnum)!=Teuchos::null,"Material of element is null pointer!");
  const MAT::Material& material = *( phasemanager_->Element()->Material(matnum));

  // get number of phases
  const int numfluidphases = phasemanager_->NumFluidPhases();

  // resize vectors
  permeabilitytensors_.resize(numfluidphases);
  relpermeabilities_.resize(numfluidphases);
  derrelpermeabilities_.resize(numfluidphases);

  // check material type
  if (material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase
      and material.MaterialType()
          != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("only poro multiphase material valid");

    // cast to multiphase material
  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  for (int iphase = 0; iphase < numfluidphases; iphase++)
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
    double abspressgrad,
    int matnum) const
{
  phasemanager_->CheckIsEvaluated();

  return DynViscosity(*Element()->Material(matnum),phasenum,abspressgrad);
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

/*----------------------------------------------------------------------*
 | constructor                                         kremheller 08/17 |
 *----------------------------------------------------------------------*/
template<int nsd>
DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<nsd>::PhaseManagerVolFrac(
    Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager
    )
: PhaseManagerDecorator(phasemanager)
{

  return;
}

/*----------------------------------------------------------------------*
 | constructor                                         kremheller 08/17 |
 *----------------------------------------------------------------------*/
template<int nsd>
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<nsd>::Setup(
    const DRT::Element* ele,
    const int matnum
    )
{
  // setup the wrapped class
  phasemanager_->Setup(ele, matnum);

  const int totalnumdof = phasemanager_->TotalNumDof();
  const int numfluidphases = phasemanager_->NumFluidPhases();
  const int numvolfrac = totalnumdof - numfluidphases;

  if(numfluidphases >= totalnumdof)
    dserror("We should not be here, total numdof is %d, numfluidphases is %d",
        phasemanager_->TotalNumDof(), phasemanager_->NumFluidPhases());

  difftensorsvolfrac_.resize(numvolfrac);
  volfracdensity_.resize(numvolfrac);
  volfracpressure_.resize(numvolfrac);
  scalardiffs_.resize(numvolfrac);

  hasaddscalardpendentflux_.resize(numvolfrac,false);

  // get material
  const MAT::Material& material = *(ele->Material(matnum));

  // check the material
  if(material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase and
     material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("only poro multiphase and poro multiphase reactions material valid");

  for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++)
  {
    //get the single phase material
    const MAT::FluidPoroSingleVolFrac& singlevolfracmat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSingleVolFracMatFromMaterial(material,ivolfrac+numfluidphases);

    // clear
    difftensorsvolfrac_[ivolfrac].Clear();

    // TODO only isotropic, constant diffusivity for now
    const double diffusivity = singlevolfracmat.Diffusivity();
    for (int i = 0; i < nsd; i++)
      (difftensorsvolfrac_[ivolfrac])(i, i) = diffusivity;

    volfracdensity_[ivolfrac] = singlevolfracmat.Density();
    volfracpressure_[ivolfrac] = singlevolfracmat.Pressure();

    if(singlevolfracmat.HasAddScalarDependentFlux())
    {
      if(phasemanager_->NumScal() != singlevolfracmat.NumScal())
        dserror("Wrong number of scalars for additional scalar dependent flux");
      hasaddscalardpendentflux_[ivolfrac] = true;
      scalardiffs_[ivolfrac].resize(singlevolfracmat.NumScal());
      scalardiffs_[ivolfrac] = singlevolfracmat.ScalarDiffs();
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 | constructor                                         kremheller 08/17 |
 *----------------------------------------------------------------------*/
template<int nsd>
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<nsd>::EvaluateGPState(
    double                 J,
    const VariableManagerMinAccess& varmanager,
    const int matnum
    )
{
  // evaluate wrapped manager
  phasemanager_->EvaluateGPState( J, varmanager, matnum);

  // safety check
  if(phasemanager_->InvBulkmodulusSolid() > 1.0e-14)
    dserror("So far volume fractions are only possible for an incompressible solid");

  // recalculate the solid pressure
  phasemanager_->RecalculateSolidPressure(phasemanager_->Porosity(),volfracpressure_);

  return;
 }

/*---------------------------------------------------------------------------*
 * get diffusion tensor for volume fractions                kremheller 08/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<nsd>::DiffTensorVolFrac(
    int volfracnum, LINALG::Matrix<nsd,nsd>& difftensorvolfrac) const
{
  phasemanager_->CheckIsEvaluated();
  // make a hard copy for now
  difftensorvolfrac = difftensorsvolfrac_[volfracnum];
}

/*---------------------------------------------------------------------------*
 * get densities for volume fractions                       kremheller 08/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<nsd>::VolFracDensity(
    int volfracnum) const
{
  phasemanager_->CheckIsEvaluated();

  return volfracdensity_[volfracnum];
}

/*---------------------------------------------------------------------------*
 * get pressures for volume fractions                       kremheller 08/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<nsd>::VolFracPressure(
    int volfracnum) const
{
  phasemanager_->CheckIsEvaluated();

  return volfracpressure_[volfracnum];
}

/*---------------------------------------------------------------------------*
 * get densities for volume fractions                       kremheller 08/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
bool DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<nsd>::HasAddScalarDependentFlux(
    int volfracnum) const
{
  phasemanager_->CheckIsEvaluated();

  return hasaddscalardpendentflux_[volfracnum];
}

/*---------------------------------------------------------------------------*
 * get densities for volume fractions                       kremheller 08/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
bool DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<nsd>::HasAddScalarDependentFlux(
    int volfracnum,
    int iscal) const
{
  phasemanager_->CheckIsEvaluated();

  return fabs((scalardiffs_[volfracnum])[iscal]) > 1.0e-12;
}

/*---------------------------------------------------------------------------*
 * get densities for volume fractions                       kremheller 08/17 |
*---------------------------------------------------------------------------*/
template<int nsd>
double DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<nsd>::ScalarDiff(
    int volfracnum,
    int iscal) const
{
  phasemanager_->CheckIsEvaluated();

  return (scalardiffs_[volfracnum])[iscal];
}

/*----------------------------------------------------------------------*
 | zero all values at GP                               kremheller 08/17 |
 *----------------------------------------------------------------------*/
template<int nsd>
void DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<nsd>::ClearGPState()
{
  // this is just a safety call. All quantities should be recomputed
  // in the next EvaluateGPState call anyway, but you never know ...
  phasemanager_->ClearGPState();

  return;
}

///*----------------------------------------------------------------------*
// *----------------------------------------------------------------------*/
//// template classes

template class DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<1>;
template class DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<2>;
template class DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerDiffusion<3>;
template class DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<1>;
template class DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<2>;
template class DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerVolFrac<3>;

