/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_artery_coupling_pair.cpp

 \brief one pair consisting of exactly one artery element and one poro-
        multiphase-scatra element which might be tied to each other

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_artery_coupling_pair.H"
#include "poromultiphase_scatra_artery_coupling_defines.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mat/fluidporo_multiphase.H"
#include "../drt_mat/scatra_mat_multiporo.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_parameter.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/cnst_1d_art.H"
#include "../headers/FAD_utils.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_parameter.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"

/*----------------------------------------------------------------------*
 | constructor                                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::PoroMultiPhaseScatraArteryCouplingPair() :
PoroMultiPhaseScatraArteryCouplingPairBase(),
coupltype_(CouplingType::type_undefined),
couplmethod_(INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::none),
isinit_(false),
ispreevaluated_(false),
isactive_(false),
funct_coupl_active_(false),
evaluate_in_ref_config_(true),
element1_(NULL),
element2_(NULL),
arterydiam_(0.0),
numdof_cont_(0),
numdof_art_(0),
numcoupleddofs_(0),
numfluidphases_(0),
numvolfrac_(0),
numscalcont_(0),
numscalart_(0),
nds_porofluid_(-1),
n_gp_(0),
arteryelelengthref_(0.0),
jacobi_(0.0),
pp_(0.0),
eta_a_(0.0),
eta_b_(0.0),
curr_segment_length_(0.0),
constant_part_evaluated_(false),
porosityname_("porosity"),
artpressname_("p_art"),
segmentid_(-1),
timefacrhs_art_(0.0),
timefacrhs_cont_(0.0)
{

  return;
}

/*----------------------------------------------------------------------*
 | init                                                kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::Init(
    std::vector< DRT::Element const*> elements,
    const Teuchos::ParameterList&      meshtyingparams,
    const std::vector<int>& coupleddofs_cont,
    const std::vector<int>& coupleddofs_art,
    const std::vector<std::vector<int>>& scale_vec,
    const std::vector<std::vector<int>>& funct_vec
    )
{

  // init stuff
  couplmethod_ =
      DRT::INPUT::IntegralValue<INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod>(meshtyingparams,"ARTERY_COUPLING_METHOD");

  evaluate_in_ref_config_ =
      DRT::INPUT::IntegralValue<int>(meshtyingparams,"EVALUATE_IN_REF_CONFIG");

  element1_ = elements[0];
  element2_ =  elements[1];

  // set coupling type
  if(element1_->ElementType().Name() == "ArteryType" && element2_->ElementType().Name() == "PoroFluidMultiPhaseType")
  {
    coupltype_ = type_porofluid;
    nds_porofluid_ = 0;
  }
  else if(element1_->ElementType().Name() == "TransportType" && element2_->ElementType().Name() == "TransportType")
  {
    coupltype_ = type_scatra;
    nds_porofluid_ = 2;
  }
  else
    dserror("Your selected coupling is not possible, type of element1: "+element1_->ElementType().Name() + ", type of element2: " + element2_->ElementType().Name());

  // get number of DOFs of artery or artery-scatra
  const DRT::Node* const* artnodes;
  artnodes = element1_->Nodes();

  numdof_art_ = element1_->NumDofPerNode(*artnodes[0]);
  for (int i = 1; i < element1_->NumNode(); i++)
    if(numdof_art_ != element1_->NumDofPerNode(*artnodes[i]))
      dserror("It is not possible to have different number of Dofs in artery discretization");

  // get number of DOFs of continuous ele (scatra or porofluid)
  const DRT::Node* const* contnodes;
  contnodes = element2_->Nodes();

  numdof_cont_ = element2_->NumDofPerNode(*contnodes[0]);
  for (int i = 1; i < element2_->NumNode(); i++)
    if(numdof_cont_ != element2_->NumDofPerNode(*contnodes[i]))
      dserror("It is not possible to have different number of Dofs in continuos discretization");

  // safety check
  if(numdof_art_ != (int)(scale_vec[0].size()))
    dserror("Wrong size of scale-vector (artery)");
  if(numdof_art_ != (int)(funct_vec[0].size()))
    dserror("Wrong size of function-vector (artery)");
  if(numdof_cont_ != (int)(scale_vec[1].size()))
    dserror("Wrong size of scale-vector (continuous discretization)");
  if(numdof_cont_ != (int)(funct_vec[1].size()))
    dserror("Wrong size of function-vector (continuous discretization)");

  // fill scale vector
  scale_vec_ = scale_vec;

  // fill function vector
  funct_vec_.resize(2);
  funct_vec_[0].resize(numdof_art_);
  funct_vec_[1].resize(numdof_cont_);
  FillFunctionVector(&funct_vec_[0], funct_vec[0], scale_vec_[0]);
  FillFunctionVector(&funct_vec_[1], funct_vec[1], scale_vec_[1]);

  // get the actually coupled dofs
  coupleddofs_cont_ = coupleddofs_cont;
  // Note: this will be overwritten in case of arteryscatra-scatra coupling
  volfracpressid_ = coupleddofs_cont;
  coupleddofs_art_  = coupleddofs_art;
  numcoupleddofs_ = coupleddofs_cont.size();

  // safety check
  for(int icont = 0; icont < numcoupleddofs_; icont++)
    if(coupleddofs_cont_[icont] >= numdof_cont_)
      dserror("You try to couple DOF %d, which is larger than the number of dofs of the continuous discretization", coupleddofs_cont_[icont]+1);
  for(int iart = 0; iart < numcoupleddofs_; iart++)
    if(coupleddofs_art_[iart] >= numdof_art_)
      dserror("You try to couple DOF %d, which is larger than the number of dofs of the artery discretization", coupleddofs_art_[iart]+1);

  // Set reference nodal positions for artery element
  for (unsigned int n=0;n<numnodesart_;++n)
  {
    const DRT::Node* node = element1_->Nodes()[n];
    for (unsigned int d=0;d<numdim_;++d)
      ele1posref_(numdim_*n+d) = node->X()[d];
  }

  // get length of 1D element
  static LINALG::Matrix<numdim_,1> arterypos0;
  for (unsigned int d=0;d<numdim_;++d)
    arterypos0(d) = ele1posref_(d);
  static LINALG::Matrix<numdim_,1> arterypos1;
  for (unsigned int d=0;d<numdim_;++d)
    arterypos1(d) = ele1posref_(numdim_+d);

  static LINALG::Matrix<numdim_,1> dist;
  dist.Update(-1.0, arterypos0, 1.0, arterypos1, 0.0);
  arteryelelengthref_ = dist.Norm2();

  // get initial direction of artery elemetn
  t0_.Update(1.0/arteryelelengthref_, dist, 0.0);

  // Set reference nodal positions for continuous discretization element
  for (unsigned int inode=0; inode<numnodescont_; ++inode)
  {
    const DRT::Node* node = element2_->Nodes()[inode];
    for (unsigned int idim=0; idim<numdim_; ++idim)
      ele2posref_(idim, inode) = node->X()[idim];
  }

  // set current nodal positions to reference nodal positions for continuous discretization element
  ele2pos_.Update(1.0,ele2posref_,0.0);

  // get penalty parameter
  pp_ = meshtyingparams.get<double>("PENALTY");

  // get out of here
  isinit_ = true;

  return;
}

/*----------------------------------------------------------------------*
 | setup the fluid managers and materials              kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::SetupFluidManagersAndMaterials(
    const std::string disname
    )
{

  // dummy parameter list
  DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter* para =
      DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance(disname);

  Teuchos::RCP<MAT::FluidPoroMultiPhase> multiphasemat = Teuchos::null;
  switch(coupltype_)
  {
  case type_porofluid:
  {
    multiphasemat =
        Teuchos::rcp_dynamic_cast<MAT::FluidPoroMultiPhase>(element2_->Material(nds_porofluid_));
    if(multiphasemat==Teuchos::null)
      dserror("cast to MAT::FluidPoroMultiPhase failed for artery-porofluid coupling!");
    for(int idof = 0; idof < numcoupleddofs_; idof++)
    {
      const int matid = multiphasemat->MatID(coupleddofs_cont_[idof]);
      Teuchos::RCP< MAT::Material> singlemat = multiphasemat->MaterialById(matid);

      // safety check
      if(singlemat->MaterialType() != INPAR::MAT::m_fluidporo_volfracpressure)
        dserror("You can only couple volume fraction pressures, your material is of type %d", singlemat->MaterialType());
    }
    // we have a coupling with scatra -> the scatra-material is the third material in the 2D/3D element
    if(element2_->NumMaterial() == 3)
    {
      if (element2_->Material(2)->MaterialType() == INPAR::MAT::m_matlist or
          element2_->Material(2)->MaterialType() == INPAR::MAT::m_matlist_reactions)
      {
        Teuchos::RCP<MAT::MatList> scatramat = Teuchos::rcp_static_cast<MAT::MatList>(element2_->Material(2));
        numscalcont_ = scatramat->NumMat();
      }
      else
        dserror("Only MAT_matlist is valid for poromultiphase-scatra material");
    }
    // we have a coupling with artery-scatra -> the artery-scatra-material is the second material in the artery element
    if(element1_->NumMaterial() == 2)
    {
      if (element1_->Material(1)->MaterialType() == INPAR::MAT::m_matlist)
      {
        Teuchos::RCP<MAT::MatList> artscatramat = Teuchos::rcp_static_cast<MAT::MatList>(element1_->Material(1));
        numscalart_ = artscatramat->NumMat();
      }
      else if(element1_->Material(1)->MaterialType() == INPAR::MAT::m_scatra)
        numscalart_ = 1;
      else
        dserror("Only MAT_matlist and MAT_scatra are valid for artery-scatra material");
    }

    Teuchos::RCP<MAT::Cnst_1d_art> arterymat = Teuchos::rcp_static_cast<MAT::Cnst_1d_art>(element1_->Material(0));
    if(arterymat==Teuchos::null)
      dserror("cast to artery material failed for porofluid-artery coupling!");
    arterydiam_ = arterymat->Diam();

    break;
  }
  case type_scatra:
  {
    // check if we actually have three materials
    if(element2_->NumMaterial()<3)
      dserror("no third material available");

    multiphasemat =
        Teuchos::rcp_dynamic_cast<MAT::FluidPoroMultiPhase>(element2_->Material(nds_porofluid_));

    if(multiphasemat==Teuchos::null)
      dserror("cast to MAT::FluidPoroMultiPhase failed for arteryscatra-scatra coupling!");

    Teuchos::RCP<MAT::MatList> scatramat = Teuchos::rcp_static_cast<MAT::MatList>(element2_->Material(0));
    if(scatramat==Teuchos::null)
      dserror("cast to ScatraMat failed for arteryscatra-scatra coupling!");

    for(int idof = 0; idof < numcoupleddofs_; idof++)
    {
      const int matid = scatramat->MatID(coupleddofs_cont_[idof]);
      Teuchos::RCP< MAT::Material> singlemat = scatramat->MaterialById(matid);
      if(singlemat->MaterialType() == INPAR::MAT::m_scatra_multiporo_volfrac)
      {
        const Teuchos::RCP<const MAT::ScatraMatMultiPoroVolFrac>& poromat
             = Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroVolFrac>(singlemat);
        volfracpressid_[idof] = poromat->PhaseID()+multiphasemat->NumVolFrac();
      }
      else
        dserror("You can only couple MAT::ScatraMatMultiPoroVolFrac, your material is of type %d", singlemat->MaterialType());
    }
    // get the artery scatra-material
    if (element1_->Material(0)->MaterialType() == INPAR::MAT::m_matlist)
    {
      Teuchos::RCP<MAT::MatList> artscatramat = Teuchos::rcp_static_cast<MAT::MatList>(element1_->Material(0));
      numscalart_ = artscatramat->NumMat();
    }
    else if(element1_->Material(0)->MaterialType() == INPAR::MAT::m_scatra)
      numscalart_ = 1;
    else
      dserror("Only MAT_matlist and MAT_scatra are valid for artery-scatra material");
    numscalcont_ = numdof_cont_;
    Teuchos::RCP<MAT::Cnst_1d_art> arterymat = Teuchos::rcp_static_cast<MAT::Cnst_1d_art>(element1_->Material(1));
    if(arterymat==Teuchos::null)
      dserror("cast to artery material failed for arteryscatra-scatra coupling!");
    arterydiam_ = arterymat->Diam();
    break;
  }
  default:
    dserror("Unknown coupling type");
    break;
  }

  numfluidphases_ = multiphasemat->NumFluidPhases();
  numvolfrac_ = multiphasemat->NumVolFrac();

  // create phase-manager
  phasemanager_ =
      DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface::CreatePhaseManager(
          *para,
          numdim_,
          multiphasemat->MaterialType(),
          POROFLUIDMULTIPHASE::Action::get_access_from_artcoupling,
          multiphasemat->NumMat(),
          multiphasemat->NumFluidPhases()
          );

  // setup phasemanager
  phasemanager_->Setup(element2_,nds_porofluid_);

  // create variablemanager
  variablemanager_ =
      DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<numdim_,numnodescont_>::CreateVariableManager(
          *para,
          POROFLUIDMULTIPHASE::Action::get_access_from_artcoupling,
          multiphasemat,
          multiphasemat->NumMat(),
          multiphasemat->NumFluidPhases()
          );

  // initialize the names used in functions
  InitializeFunctionNames();

  // initialize the functions
  for(int i = 0; i < 2; i++)
    for(unsigned int idof = 0; idof < funct_vec_[i].size(); idof++)
      if(funct_vec_[i][idof] != 0)
        InitializeFunction(funct_vec_[i][idof]);

  // set time fac for right hand side evaluation of coupling
  SetTimeFacRhs();

  return;
}

/*----------------------------------------------------------------------*
 | pre-evaluate                                        kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::PreEvaluate()
{

  if(!isinit_)
    dserror("MeshTying Pair has not yet been initialized");

  // number of checks
  // we check at NUMPROJCHECKS (default = 11) positions from [-0.99, ..., 0, ..., 0.99]
  // if these points on the artery element can be projected
  std::vector<double> check_proj(NUMPROJCHECKS);
  check_proj[0] = -0.99;
  const double dist = 1.98/((double)(NUMPROJCHECKS) - 1.0);
  for(int i = 1; i < NUMPROJCHECKS; i++)
    check_proj[i] = -0.99 + (double)(i)*dist;
  std::vector<bool> validprojections(NUMPROJCHECKS, false);

  std::vector<double> xi(numdim_);

  for (int i = 0; i < NUMPROJCHECKS; i++)
  {
    bool projection_valid = false;
    if(PROJOUTPUT)
      std::cout << "projection for " << check_proj[i] << " ========================================" << std::endl;
    Projection<double>(check_proj[i], xi, projection_valid);
    if(projection_valid)
    {
      isactive_ = true;
      validprojections[i] = true;
    }
  }

  // no projection found
  if(!isactive_)
  {
    return;
  }

  // choice of optimal Gauss-rule:
  // basically the N^(2)*N^(2) term is crucial
  // for (bi-, tri-)linear elements (only considered so far):
  // in 2D the highest possible polynomial order for this term is 4 since N^(2) can be quadratic for arbitrary integration in element
  // --> we need 3 gp for exact integration
  // in 3D the highest possible polynomial order for this term is 6 since N^(2) can be cubic for arbitrary integration in element
  // --> we need 4 gp for exact integration

  DRT::UTILS::IntegrationPoints1D gaussPoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GaussRule1D::intrule_line_3point);
  if(numdim_ == 3)
    gaussPoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GaussRule1D::intrule_line_4point);

  n_gp_ = gaussPoints.nquad;
  // define Gauss points and n_gp-sized quantities
  eta_.resize(n_gp_);
  eta_s_.resize(n_gp_);
  wgp_.resize(n_gp_);
  invJ_.resize(n_gp_);
  xi_.resize(n_gp_);
  for(int i_gp = 0; i_gp < n_gp_; i_gp++)
    xi_[i_gp].resize(numdim_);

  // integration segment [eta_a, eta_b] is created
  CreateIntegrationSegment(validprojections);

  // get jacobian determinant
  const double determinant = (eta_b_ - eta_a_)/2.0;
  jacobi_ = determinant*arteryelelengthref_/2.0;

  static LINALG::Matrix<1, numnodescont_> N2(true);                     // = N2
  static LINALG::Matrix<numdim_, numnodescont_> N2_xi(true);             // = N2,xi1

  // project the Gauss points --> those have to able to be projected
  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    // compute the coordinate trafo, weights and project Gauss points
    eta_[i_gp] = (eta_a_ + eta_b_)/2.0 + gaussPoints.qxg[i_gp][0] * determinant;
    eta_s_[i_gp] = eta_[i_gp];
    wgp_[i_gp] = gaussPoints.qwgt[i_gp];
    bool projection_valid = false;
    Projection<double>(eta_[i_gp], xi_[i_gp], projection_valid);
    if(!projection_valid)
      dserror("Gauss point could not be projected");

    // compute (dX/dxi)^-1
    Get2D3DShapeFunctions<double>(N2, N2_xi, xi_[i_gp]);
    invJ_[i_gp].MultiplyNT(N2_xi,ele2pos_);
    invJ_[i_gp].Invert();
  }

  ispreevaluated_ = true;

  return;
}

/*----------------------------------------------------------------------*
 | reset state                                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::ResetState(
    const std::vector<double>& contelephinp,
    const std::vector<double>& artelephinp,
    Teuchos::RCP<DRT::Discretization> contdis,
    Teuchos::RCP<DRT::Discretization> artdis
    )
{

  if(!ispreevaluated_)
    dserror("MeshTying Pair has not yet been pre-evaluated");

  if(contelephinp.size() != numnodescont_*numdof_cont_)
    dserror("Mismatch in size for continuous element, expected %d Dofs, but got %d", numnodescont_*numdof_cont_, contelephinp.size());

  if(artelephinp.size() != numnodesart_*numdof_art_)
    dserror("Mismatch in size for artery element, expected %d Dofs, but got %d", numnodesart_*numdof_art_, artelephinp.size());

  // reset
  contelephinp_ = contelephinp;
  artelephinp_  = artelephinp;

  DRT::Element::LocationArray la(contdis->NumDofSets());
  element2_->LocationVector(*contdis,la,false);


  switch(coupltype_)
  {
  case type_porofluid:
  {
    // extract element and node values of fluid
    variablemanager_->ExtractElementAndNodeValues(*element2_,*contdis,la,ele2pos_,0);

    // extract values of artery-scatra discretization
    if(numscalart_ > 0)
    {
      Teuchos::RCP<const Epetra_Vector> artscalarnp = artdis->GetState(ndsartery_scatra_, "one_d_artery_phinp");
      if (artscalarnp!=Teuchos::null)
      {
        DRT::Element::LocationArray la(artdis->NumDofSets());
        element1_->LocationVector(*artdis,la,false);
        // rebuild scalar vector
        eartscalarnp_.clear();
        eartscalarnp_.resize(numscalart_,LINALG::Matrix<numnodesart_,1>(true));
        // extract local values of artery-scatra field from global state vector
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<numnodesart_,1> >(*artscalarnp,eartscalarnp_,la[ndsartery_scatra_].lm_);
      }
      else
        dserror("Cannot get artery-scatra from artery discretization");
    }
    // extract values of continuous scatra discretization
    if(numscalcont_ > 0)
    {
      // get state vector from discretization
      Teuchos::RCP<const Epetra_Vector> contscalarnp = contdis->GetState(3, "scalars");
      if (contscalarnp!=Teuchos::null)
      {
        DRT::Element::LocationArray la(contdis->NumDofSets());
        element2_->LocationVector(*contdis,la,false);
        // rebuild scalar vector
        econtscalarnp_.clear();
        econtscalarnp_.resize(numscalcont_,LINALG::Matrix<numnodescont_,1>(true));
        // extract local values of continuous-scatra field from global state vector
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<numnodescont_,1> >(*contscalarnp,econtscalarnp_,la[3].lm_);
      }
      else
        dserror("Cannot get state vector 'scalars'");
    }
    break;
  }
  case type_scatra:
  {
    // extract element and node values of fluid
    variablemanager_->ExtractElementAndNodeValues(*element2_,*contdis,la,ele2pos_,2);
    //extract artery pressure
    Teuchos::RCP<const Epetra_Vector> artpressnp = artdis->GetState(ndsscatra_artery_, "one_d_artery_pressure");
    if (artpressnp!=Teuchos::null)
    {
      DRT::Element::LocationArray la(artdis->NumDofSets());
      element1_->LocationVector(*artdis,la,false);
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<numnodesart_,1> >(*artpressnp,earterypressurenp_,la[ndsscatra_artery_].lm_);
    }
    else
      dserror("Cannot get arterypressure from artery-scatra discretization");
    break;
  }
  default:
    dserror("Unknown coupling type");
    break;
  }

}

/*----------------------------------------------------------------------*
 | evaluate                                            kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::Evaluate(
    LINALG::SerialDenseVector* forcevec1,
    LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11,
    LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21,
    LINALG::SerialDenseMatrix* stiffmat22,
    LINALG::SerialDenseMatrix* D_ele,
    LINALG::SerialDenseMatrix* M_ele,
    const std::vector<double>& segmentlengths
    )
{

  if(!ispreevaluated_)
    dserror("MeshTying Pair has not yet been pre-evaluated");

  const int dim1 = artelephinp_.size();
  const int dim2 = contelephinp_.size();

  // resize and initialize variables to zero
  if (forcevec1 != NULL)
    forcevec1->Size(dim1);
  if (forcevec2 != NULL)
    forcevec2->Size(dim2);

  if (stiffmat11 != NULL)
    stiffmat11->Shape(dim1,dim1);
  if (stiffmat12 != NULL)
    stiffmat12->Shape(dim1,dim2);
  if (stiffmat21 != NULL)
    stiffmat21->Shape(dim2,dim1);
  if (stiffmat22 != NULL)
    stiffmat22->Shape(dim2,dim2);

  switch (couplmethod_)
  {
  case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::gpts:
  {
    EvaluateGPTS(dim1, dim2, forcevec1, forcevec2, stiffmat11, stiffmat12, stiffmat21, stiffmat22);
    break;
  }
  case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp:
  {
    EvaluateDM(dim1, dim2, D_ele, M_ele);
    break;
  }
  default:
    dserror("Unknown coupling type for line-based coupling");
    break;
  }

  if(funct_coupl_active_)
    EvaluateFunctionCoupling(segmentlengths, forcevec1, forcevec2, stiffmat11, stiffmat12, stiffmat21, stiffmat22);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate kappa                                      kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateKappa(
    LINALG::SerialDenseVector* kappa)
{

  if(!ispreevaluated_)
    dserror("MeshTying Pair has not yet been pre-evaluated");

  const int dim1 = numdof_art_*numnodesart_;

  // resize and initialize variables to zero
  if (kappa != NULL)
    kappa->Size(dim1);

  // Vectors for shape functions and their derivatives
  static LINALG::Matrix<1, numnodesart_> N1(true);         // = N1
  static LINALG::Matrix<1, numnodesart_> N1_eta(true);     // = N1,eta

  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    // Get constant values from projection
    const double w_gp = wgp_[i_gp];
    const double eta = eta_[i_gp];
    const double jac = jacobi_;

    // Update 1D shape functions
    Get1DShapeFunctions<double>(N1, N1_eta, eta);

    for(unsigned int inode = 0; inode < numnodesart_; inode++)
    {
      const double mykappa = w_gp*jac*N1(inode);
      for (int dof = 0; dof < numdof_art_; dof++)
        (*kappa)(inode*numdof_art_+dof) += mykappa;
    }
  }

  return;
}

/*------------------------------------------------------------------------*
 | element-id of artery element                          kremheller 05/18 |
 *------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
int POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::Ele1GID() const
{

  return element1_->Id();
}

/*------------------------------------------------------------------------*
 | element-id of 2D/3D-element                           kremheller 05/18 |
 *------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
int POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::Ele2GID() const
{

  return element2_->Id();
}

/*------------------------------------------------------------------------*
 | get the segment id                                    kremheller 05/18 |
 *------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
int POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::GetSegmentID() const
{

  return segmentid_;
}

/*------------------------------------------------------------------------*
 | set the segment id                                    kremheller 05/18 |
 *------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::SetSegmentID(
    const int& segmentid
    )
{

  segmentid_ = segmentid;
}

/*----------------------------------------------------------------------*
 | apply mesh movement (on artery)                     kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
double POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::ApplyMeshMovement(
    Teuchos::RCP<const Epetra_Vector> disp,
    Teuchos::RCP<DRT::Discretization> contdis
    )
{

  // nodal displacement values for ALE
  LINALG::Matrix<numdim_,numnodescont_> edispnp(true);

  if (disp!=Teuchos::null)
  {
    DRT::Element::LocationArray la(contdis->NumDofSets());
    element2_->LocationVector(*contdis,la,false);

    // construct location vector for displacement related dofs
    std::vector<int> lmdisp(numdim_*numnodescont_,-1);
    for (unsigned int inode=0; inode<numnodescont_; ++inode)
      for (unsigned int idim=0; idim<numdim_; ++idim)
        lmdisp[inode*numdim_+idim] = la[1].lm_[inode*numdim_+idim];

    // extract local values of displacement field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<numdim_,numnodescont_> >(*disp,edispnp,lmdisp);
  }
  else
    return (eta_b_ - eta_a_)/2.0*arteryelelengthref_;

  // update current configuration
  ele2pos_.Update(1.0,ele2posref_,1.0,edispnp,0.0);

  static LINALG::Matrix<1, numnodescont_> N2(true);                     // = N2
  static LINALG::Matrix<numdim_, numnodescont_> N2_xi(true);            // = N2,xi1
  static LINALG::Matrix<numdim_,numnodescont_> N2_XYZ(true);            // = N2,X
  static LINALG::Matrix<numdim_,numdim_> defgrad(true);                 // = dx/dX = F
  static LINALG::Matrix<numdim_,1> Ft0;                                 // = F*t0

  curr_segment_length_ = 0.0;
  // current segment length = \int_{\eta_a}^{eta_b} || F*t0 ||_2 ds
  // all under the assumption that artery element completely follows deformation of
  // underlying 2D/3D problem
  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    // get shape functions of continuous element
    const std::vector<double> xi = xi_[i_gp];
    Get2D3DShapeFunctions<double>(N2, N2_xi, xi);
    // dN/dX = dN/dxi * dxi/dX = dN/dxi * (dX/dxi)^-1
    N2_XYZ.Multiply(invJ_[i_gp],N2_xi);
    // dx/dX = x * N_XYZ^T
    defgrad.MultiplyNT(ele2pos_,N2_XYZ);
    Ft0.Multiply(defgrad,t0_);
    curr_segment_length_ += Ft0.Norm2()*wgp_[i_gp]*jacobi_;
  }

  return curr_segment_length_;
}

/*----------------------------------------------------------------------*
 | recompute eta and xi                                kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::RecomputeEtaAndXiInDeformedConfiguration(
    const std::vector<double>&              segmentlengths,
    std::vector<double>&                    myEta,
    std::vector<std::vector<double>>&       myXi)
{

  // NOTE: we assume that the 1D artery element completely follows the deformation of the underlying
  //       porous medium, so its length might change. Interaction between artery element and porous
  //       medium has to be evaluated in current/deformed configuration. However, Gauss points of the
  //       original projection (in reference configuration) cannot be used then anymore but we must
  //       define a new parameter space [-1, 1] which maps into the current configuration.
  //       First, we determine the new etaA and etaB of this segment as sum of segment lengths
  //       etaA_new = -1.0 + 2.0* ( \sum_{i=0}_{this_seg-1} l_i / total_ele_length )
  //       etaB_new = -1.0 + 2.0* ( \sum_{i=0}_{this_seg} l_i / total_ele_length )
  //       then GPs are distributed in the interval [etaA_new, etaB_new]
  //       The last step is to get the new projected xi_i in the 2D/3D parameter space of the second element
  //       For each new GP eta_new, this can be done by finding the point in reference configuration which
  //       deforms to the point in current configuration where the GP now lies as
  //       \int_{\eta_a}^{eta_s} || F*t0 ||_2 ds where eta_s is the unknown.
  //       Linearization of this nonlinear equation within the Newton loop is done with FAD

  // not necessary if we do not take into account mesh movement
  if(evaluate_in_ref_config_)
  {
    myEta = eta_;
    myXi  = xi_;
  }
  else
  {
    // current length of artery
    double curr_ele_length = 0.0;
    for(unsigned int iseg = 0; iseg < segmentlengths.size(); iseg++)
      curr_ele_length += segmentlengths[iseg];

    // length of segments [0, 1, ..., this_seg-1]
    double length_so_far = 0.0;
    for(int iseg = 0; iseg < segmentid_; iseg++)
      length_so_far += segmentlengths[iseg];

    // length of this segment
    const double curr_seg_length = segmentlengths[segmentid_];

    // get new etaA and etaB
    const double etaA = -1.0 + 2.0*(length_so_far/curr_ele_length);
    const double etaB = -1.0 + 2.0*((length_so_far+curr_seg_length)/curr_ele_length);

    DRT::UTILS::IntegrationPoints1D gaussPoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GaussRule1D::intrule_line_3point);
    if(numdim_ == 3)
      gaussPoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GaussRule1D::intrule_line_4point);

    // distribute new Gauss points
    const double determinant = (etaB - etaA)/2.0;
    for(int i_gp = 0; i_gp < n_gp_; i_gp++)
      myEta[i_gp] = (etaA + etaB)/2.0 + gaussPoints.qxg[i_gp][0] * determinant;

    // solution variable for Newton loop
    FAD eta_s = 0.0;
    eta_s.diff(0,1);  // independent variable 0 out of a total of 1

    // the GP loop
    bool converged = false;
    for (int i_gp = 0; i_gp < n_gp_; i_gp++)
    {
      // start value for Newton: last converged value of previous evaluation (should be pretty close)
      eta_s.val() = eta_s_[i_gp];
      const double desired_length = curr_seg_length*(myEta[i_gp]-etaA)/(etaB-etaA);
      double val = -1.0;
      // Netwon loop
      for(int istep = 0; istep < MESHMOVEMENTMAXITER; istep++)
      {
        // integrate \int_{\eta_a}^{eta_s} || F*t0 ||_2 ds
        const FAD curr_length = IntegrateLengthToEtaS(eta_s);

        val = curr_length.val()-desired_length;
        if(fabs(val) < CONVTOLNEWTONMESHMOVEMENT)
        {
          converged = true;
          break;
        }
        const double deriv = curr_length.fastAccessDx(0);
        // Newton update
        eta_s.val() -= val/deriv;
      }

      if(!converged)
        std::cout << "WARNING: could not find Gauss point position in reference configuration";
      // finally find new xi_i by projection eta_s in reference configuration
      bool projection_valid = false;
      Projection<double>(eta_s.val(), myXi[i_gp], projection_valid);
      if(!projection_valid)
        dserror("Gauss point could not be projected");
      // save the converged value
      eta_s_[i_gp] = eta_s.val();
    } // GP loop
  } // !evaluate_in_ref_config_

  return;
}

/*----------------------------------------------------------------------*
 | evaluate for gauss-point-to-segment -> directly assembled into       |
 | element stiffness matrix                            kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateGPTS(
    const int&                 dim1,
    const int&                 dim2,
    LINALG::SerialDenseVector* forcevec1,
    LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11,
    LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21,
    LINALG::SerialDenseMatrix* stiffmat22)
{

  if(numcoupleddofs_ > 0)
  {
    // we only have to this once since evaluated in reference configuration
    if(!constant_part_evaluated_)
    {
      // Vectors for shape functions and their derivatives
      static LINALG::Matrix<1, numnodesart_> N1(true);         // = N1
      static LINALG::Matrix<1, numnodesart_> N1_eta(true);     // = N1,eta

      static LINALG::Matrix<1, numnodescont_> N2(true);                     // = N2
      static LINALG::Matrix<numdim_, numnodescont_> N2_xi(true);             // = N2,xi1

      GPTS_stiffmat11_ = new LINALG::SerialDenseMatrix();
      GPTS_stiffmat12_ = new LINALG::SerialDenseMatrix();
      GPTS_stiffmat21_ = new LINALG::SerialDenseMatrix();
      GPTS_stiffmat22_ = new LINALG::SerialDenseMatrix();

      GPTS_stiffmat11_.Shape(dim1,dim1);
      GPTS_stiffmat12_.Shape(dim1,dim2);
      GPTS_stiffmat21_.Shape(dim2,dim1);
      GPTS_stiffmat22_.Shape(dim2,dim2);

      for (int i_gp = 0; i_gp < n_gp_; i_gp++)
      {
        // Get constant values from projection
        const double w_gp = wgp_[i_gp];
        const double eta = eta_[i_gp];
        const double jac = jacobi_;
        const std::vector<double> xi = xi_[i_gp];

        // Update shape functions and their derivatives for 1D and 2D/3D element
        Get1DShapeFunctions<double>(N1, N1_eta, eta);
        Get2D3DShapeFunctions<double>(N2, N2_xi, xi);

        // evaluate
        EvaluateGPTSStiff(w_gp, N1, N2, jac, pp_);
      }
    } //!constant_part_evaluated_

    UpdateGPTSStiff(*stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);
    EvaluateGPTSForce(*forcevec1, *forcevec2,*stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate D and M for mortar penalty case            kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateDM(
    const int&                 dim1,
    const int&                 dim2,
    LINALG::SerialDenseMatrix* D_ele,
    LINALG::SerialDenseMatrix* M_ele)
{

  if (D_ele != NULL)
    D_ele->Shape(dim1,dim1);
  if (M_ele != NULL)
    M_ele->Shape(dim1,dim2);

  if(numcoupleddofs_ > 0)
  {
    // we only have to this once since evaluated in reference configuration
    if(!constant_part_evaluated_)
    {
      // Vectors for shape functions and their derivatives
      static LINALG::Matrix<1, numnodesart_> N1(true);         // = N1
      static LINALG::Matrix<1, numnodesart_> N1_eta(true);     // = N1,eta

      static LINALG::Matrix<1, numnodescont_> N2(true);                     // = N2
      static LINALG::Matrix<numdim_, numnodescont_> N2_xi(true);             // = N2,xi1

      D_ = new LINALG::SerialDenseMatrix();
      M_ = new LINALG::SerialDenseMatrix();

      D_.Shape(dim1,dim1);
      M_.Shape(dim1,dim2);

      for (int i_gp = 0; i_gp < n_gp_; i_gp++)
      {
        // Get constant values from projection
        const double w_gp = wgp_[i_gp];
        const double eta = eta_[i_gp];
        const double jac = jacobi_;
        const std::vector<double> xi = xi_[i_gp];

        // Update shape functions and their derivatives for 1D and 2D/3D element
        Get1DShapeFunctions<double>(N1, N1_eta, eta);
        Get2D3DShapeFunctions<double>(N2, N2_xi, xi);

        EvaluateDM(w_gp, N1, N2, jac, pp_);
      }
    } //!constant_part_evaluated_

    UpdateDM(*D_ele, *M_ele);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate function coupling                          kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateFunctionCoupling(
    const std::vector<double>& segmentlengths,
    LINALG::SerialDenseVector* forcevec1,
    LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11,
    LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21,
    LINALG::SerialDenseMatrix* stiffmat22)
{

  std::vector<double> myEta(n_gp_);
  std::vector< std::vector<double> > myXi(n_gp_, std::vector<double>(numdim_, 0.0));

  // recompute eta and xi --> see note in this function
  RecomputeEtaAndXiInDeformedConfiguration(segmentlengths, myEta, myXi);

  // Vectors for shape functions and their derivatives
  static LINALG::Matrix<1, numnodesart_> N1(true);         // = N1
  static LINALG::Matrix<1, numnodesart_> N1_eta(true);     // = N1,eta

  static LINALG::Matrix<1, numnodescont_> N2(true);                     // = N2
  static LINALG::Matrix<numnodescont_, 1> N2_transpose(true);           // = N2^T

  static LINALG::Matrix<numdim_, numnodescont_> N2_xi(true);             // = N2,xi1
  static LINALG::Matrix<numdim_, numnodescont_> derxy(true);             // = N2,xi1

  static LINALG::Matrix<numdim_,numdim_> xjm;
  static LINALG::Matrix<numdim_,numdim_> xjm0;
  static LINALG::Matrix<numdim_,numdim_> xji;

  const double curr_seg_length = segmentlengths[segmentid_];

  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    // Get constant values from projection
    const double w_gp = wgp_[i_gp];
    const double eta = myEta[i_gp];
    const std::vector<double> xi = myXi[i_gp];

    const double jac = curr_seg_length/2.0;

    // Update shape functions and their derivatives for 1D and 2D/3D element
    Get1DShapeFunctions<double>(N1, N1_eta, eta);
    Get2D3DShapeFunctions<double>(N2, N2_xi, xi);
    N2_transpose.UpdateT(N2);

    xjm.MultiplyNT(N2_xi,ele2pos_);
    xjm0.MultiplyNT(N2_xi,ele2posref_);

    const double det = xji.Invert(xjm);
    // inverse of transposed jacobian "ds/dX"
    const double det0 = xjm0.Determinant();

    derxy.Multiply(xji,N2_xi);

    // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
    const double JacobianDefGradient = det/det0;

    // clear current gauss point data for safety
    phasemanager_->ClearGPState();
    variablemanager_->EvaluateGPVariables(N2_transpose,derxy);

    phasemanager_->EvaluateGPState(JacobianDefGradient,*variablemanager_,nds_porofluid_);

    EvaluateFunctionCoupling(w_gp, N1, N2, jac, *forcevec1, *forcevec2,*stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate element stiffness matrix for GPTS          kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateGPTSStiff(
    const double& w_gp,
    const LINALG::Matrix<1, numnodesart_>& N1,
    const LINALG::Matrix<1, numnodescont_>& N2,
    const double& jacobi,
    const double& pp
    )
{

  // Evaluate meshtying stiffness for artery element N_1^T * N_1
  for (unsigned int i = 0; i < numnodesart_; i++)
  {
    for (unsigned int j = 0; j < numnodesart_; j++)
    {
      const double stiff = pp * jacobi * w_gp * N1(i) * N1(j);
      for (int dof = 0; dof < numcoupleddofs_; dof++)
        GPTS_stiffmat11_(i*numdof_art_+coupleddofs_art_[dof], j*numdof_art_+coupleddofs_art_[dof]) += stiff;
    }
  }

  // Evaluate meshtying stiffness for artery element "mixed" N_1^T * (-N_2)
  for (unsigned int i = 0; i < numnodesart_; i++)
  {
    for (unsigned int j = 0; j < numnodescont_; j++)
    {
      const double stiff = pp * jacobi * w_gp * N1(i) * ( - N2(j) );
      for (int dof = 0; dof < numcoupleddofs_; dof++)
        GPTS_stiffmat12_(i*numdof_art_+coupleddofs_art_[dof], j*numdof_cont_+coupleddofs_cont_[dof]) += stiff;
    }
  }

  // Evaluate meshtying stiffness for continuous element "mixed" N_2^T * (-N_1)
  for (unsigned int i = 0; i < numnodescont_; i++)
  {
    for (unsigned int j = 0; j < numnodesart_; j++)
    {
      const double stiff = pp * jacobi * w_gp * N2(i) * ( - N1(j) );
      for (int dof = 0; dof < numcoupleddofs_; dof++)
        GPTS_stiffmat21_(i*numdof_cont_+coupleddofs_cont_[dof], j*numdof_art_+coupleddofs_art_[dof]) += stiff;
    }
  }

  // Evaluate meshtying stiffness for continuous element N_2^T * N_2
  for (unsigned int i = 0; i < numnodescont_; i++)
  {
    for (unsigned int j = 0; j < numnodescont_; j++)
    {
      const double stiff = pp * jacobi * w_gp * N2(i) * N2(j);
      for (int dof = 0; dof < numcoupleddofs_; dof++)
        GPTS_stiffmat22_(i*numdof_cont_+coupleddofs_cont_[dof], j*numdof_cont_+coupleddofs_cont_[dof]) += stiff;
    }
  }

  constant_part_evaluated_ = true;

  return;
}

/*----------------------------------------------------------------------*
 | evaluate mortar coupling matrices D and M for MP    kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateDM(
    const double& w_gp,
    const LINALG::Matrix<1, numnodesart_>& N1,
    const LINALG::Matrix<1, numnodescont_>& N2,
    const double& jacobi,
    const double& pp
    )
{

  // Evaluate element mortar coupling operator D = N_1^T * N_1
  for (unsigned int i = 0; i < numnodesart_; i++)
  {
    for (unsigned int j = 0; j < numnodesart_; j++)
    {
      const double D = jacobi * w_gp * N1(i) * N1(j);
      for (int dof = 0; dof < numcoupleddofs_; dof++)
        D_(i*numdof_art_+coupleddofs_art_[dof], j*numdof_art_+coupleddofs_art_[dof]) += D;
    }
  }

  // Evaluate element mortar coupling operator M = N_1^T * N_2
  for (unsigned int i = 0; i < numnodesart_; i++)
  {
    for (unsigned int j = 0; j < numnodescont_; j++)
    {
      const double M = jacobi * w_gp * N1(i) * N2(j);
      for (int dof = 0; dof < numcoupleddofs_; dof++)
        M_(i*numdof_art_+coupleddofs_art_[dof], j*numdof_cont_+coupleddofs_cont_[dof]) += M;
    }
  }

  constant_part_evaluated_ = true;

  return;
}

/*----------------------------------------------------------------------*
 | evaluate element force for GPTS                     kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateGPTSForce(
    LINALG::SerialDenseVector& forcevec1,
    LINALG::SerialDenseVector& forcevec2,
    const LINALG::SerialDenseMatrix& stiffmat11,
    const LINALG::SerialDenseMatrix& stiffmat12,
    const LINALG::SerialDenseMatrix& stiffmat21,
    const LINALG::SerialDenseMatrix& stiffmat22
    )
{

  const int dim1 = artelephinp_.size();
  const int dim2 = contelephinp_.size();

  // Evaluate meshtying forces for artery element
  for (int i = 0; i < dim1; i++)
    for (int j = 0; j < dim1; j++)
      forcevec1(i) -= stiffmat11(i,j) * artelephinp_[j];

  for (int i = 0; i < dim1; i++)
    for (int j = 0; j < dim2; j++)
      forcevec1(i) -= stiffmat12(i,j) * contelephinp_[j];

  // Evaluate meshtying forces for continuous-dis element
  for (int i = 0; i < dim2; i++)
    for (int j = 0; j < dim1; j++)
      forcevec2(i) -= stiffmat21(i,j) * artelephinp_[j];

  for (int i = 0; i < dim2; i++)
    for (int j = 0; j < dim2; j++)
      forcevec2(i) -= stiffmat22(i,j) * contelephinp_[j];

  return;
}

/*----------------------------------------------------------------------*
 | update element stiffness matrix for GPTS --> without valid volume    |
 | fraction pressure no coupling                       kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::UpdateGPTSStiff(
    LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22)
{

  stiffmat11.Update(1.0,GPTS_stiffmat11_,0.0);
  stiffmat12.Update(1.0,GPTS_stiffmat12_,0.0);
  stiffmat21.Update(1.0,GPTS_stiffmat21_,0.0);
  stiffmat22.Update(1.0,GPTS_stiffmat22_,0.0);

  for(int idof = 0; idof < numcoupleddofs_; idof++)
  {
    // this coupling is only possible if we also have an element with a valid volume fraction pressure,
    // i.e., if we also have a smeared representation of the neovasculature at this point
    // if not ---> corresponding matrices are set to zero
    if(!variablemanager_->ElementHasValidVolFracPressure(volfracpressid_[idof]-numfluidphases_-numvolfrac_))
    {
      // reset to zero for this dof
      for (unsigned int i = 0; i < numnodesart_; i++)
        for (unsigned int j = 0; j < numnodesart_; j++)
          stiffmat11(i*numdof_art_+coupleddofs_art_[idof], j*numdof_art_+coupleddofs_art_[idof]) = 0.0;

      for (unsigned int i = 0; i < numnodesart_; i++)
        for (unsigned int j = 0; j < numnodescont_; j++)
          stiffmat12(i*numdof_art_+coupleddofs_art_[idof], j*numdof_cont_+coupleddofs_cont_[idof]) = 0.0;

      for (unsigned int i = 0; i < numnodescont_; i++)
        for (unsigned int j = 0; j < numnodesart_; j++)
          stiffmat21(i*numdof_cont_+coupleddofs_cont_[idof], j*numdof_art_+coupleddofs_art_[idof]) = 0.0;

      for (unsigned int i = 0; i < numnodescont_; i++)
        for (unsigned int j = 0; j < numnodescont_; j++)
          stiffmat22(i*numdof_cont_+coupleddofs_cont_[idof], j*numdof_cont_+coupleddofs_cont_[idof]) = 0.0;

    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | update D and M for MP --> without valid volume fraction pressure     |
 | no coupling                                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::UpdateDM(
    LINALG::SerialDenseMatrix& D_ele,
    LINALG::SerialDenseMatrix& M_ele)
{

  D_ele.Update(1.0,D_,0.0);
  M_ele.Update(1.0,M_,0.0);

  for(int idof = 0; idof < numcoupleddofs_; idof++)
  {
    // this coupling is only possible if we also have an element with a valid volume fraction pressure,
    // i.e., if we also have a smeared representation of the neovasculature at this point
    // if not ---> corresponding matrices are set to zero
    if(!variablemanager_->ElementHasValidVolFracPressure(volfracpressid_[idof]-numfluidphases_-numvolfrac_))
    {
      // reset to zero for this dof
      for (unsigned int i = 0; i < numnodesart_; i++)
        for (unsigned int j = 0; j < numnodesart_; j++)
          D_ele(i*numdof_art_+coupleddofs_art_[idof], j*numdof_art_+coupleddofs_art_[idof]) = 0.0;

      for (unsigned int i = 0; i < numnodesart_; i++)
        for (unsigned int j = 0; j < numnodescont_; j++)
          M_ele(i*numdof_art_+coupleddofs_art_[idof], j*numdof_cont_+coupleddofs_cont_[idof]) = 0.0;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate function coupling at GP                    kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateFunctionCoupling(
    const double& w_gp,
    const LINALG::Matrix<1, numnodesart_>& N1,
    const LINALG::Matrix<1, numnodescont_>& N2,
    const double& jacobi,
    LINALG::SerialDenseVector& forcevec1,
    LINALG::SerialDenseVector& forcevec2,
    LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22)
{
  // resize
  std::vector<double> artscalarnpAtGP(numscalart_,0.0);
  std::vector<double> contscalarnpAtGP(numscalcont_,0.0);
  double artpressAtGP = 0.0;

  // get artery values at GP
  GetArteryValuesAtGP(N1, artpressAtGP, artscalarnpAtGP);
  // get scatra values at GP
  GetContScalarValuesAtGP(N2, contscalarnpAtGP);
  // NOTE: values of fluid held by managers

  // artery functions
  for(int i_art = 0; i_art < numdof_art_; i_art++)
  {
    if(funct_vec_[0][i_art] != 0)
    {
      // resize
      std::vector<double> artderivs(numdof_art_, 0.0);
      std::vector<double> contderivs(numdof_cont_,0.0);
      double functval = 0.0;
      // evaluate and assemble
      EvaluateFunctionAndDeriv(funct_vec_[0][i_art], artpressAtGP, artscalarnpAtGP, contscalarnpAtGP, functval, artderivs, contderivs);
      AssembleFunctionCouplingIntoForceStiffArt(i_art, w_gp, N1, N2, jacobi, scale_vec_[0][i_art], functval, artderivs, contderivs, forcevec1, stiffmat11, stiffmat12);
    }
  }
  // continuous discretization functions
  for(int i_cont = 0; i_cont < numdof_cont_; i_cont++)
  {
    if(funct_vec_[1][i_cont] != 0)
    {
      // resize
      std::vector<double> artderivs(numdof_art_, 0.0);
      std::vector<double> contderivs(numdof_cont_,0.0);
      double functval = 0.0;
      // evaluate and assemble
      EvaluateFunctionAndDeriv(funct_vec_[1][i_cont], artpressAtGP, artscalarnpAtGP, contscalarnpAtGP, functval, artderivs, contderivs);
      AssembleFunctionCouplingIntoForceStiffCont(i_cont, w_gp, N1, N2, jacobi, scale_vec_[1][i_cont], functval, artderivs, contderivs, forcevec2, stiffmat21, stiffmat22);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | set time fac for rhs terms                          kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::SetTimeFacRhs()
{

  switch(coupltype_)
  {
  case type_porofluid:
  {
    DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter* eleparams = DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance("porofluid");
    timefacrhs_art_  = 1.0;
    timefacrhs_cont_ = eleparams->TimeFacRhs();
    break;
  }
  case type_scatra:
  {
    DRT::ELEMENTS::ScaTraEleParameterTimInt* eleparams = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra");
    timefacrhs_art_  = eleparams->TimeFacRhs();
    timefacrhs_cont_ = eleparams->TimeFacRhs();
    break;
  }
  default:
    dserror("Unknown coupling type");
    break;
  }

  return;
}

/*----------------------------------------------------------------------*
 | integrate \int_{\eta_a}^{eta_s} || F*t0 ||_2 ds     kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
FAD POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::IntegrateLengthToEtaS(
    const FAD& eta_s)
{


  FAD length = 0.0;

  // define GPs
  DRT::UTILS::IntegrationPoints1D gaussPoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GaussRule1D::intrule_line_3point);
  if(numdim_ == 3)
    gaussPoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GaussRule1D::intrule_line_4point);

  static LINALG::TMatrix<FAD, 1, numnodescont_> N2(true);                     // = N2
  static LINALG::TMatrix<FAD, numdim_, numnodescont_> N2_xi(true);            // = N2,xi1

  static LINALG::TMatrix<FAD, numdim_, numdim_> InvJ(true);                   // (dX/dxi)^-1
  static LINALG::TMatrix<FAD,numdim_, numdim_> defGrad(true);                 // (dX/dx) = F
  static LINALG::TMatrix<FAD,numdim_,numnodescont_> N2_XYZ(true);             // = N2,X
  static LINALG::TMatrix<FAD,numdim_,1> Ft0(true);                            // = F*t0

  // t0
  static LINALG::TMatrix<FAD,numdim_,1> t0;
  for(unsigned int i = 0; i < numdim_; i++)
    t0(i).val() = t0_(i);
  // ele2posref
  static LINALG::TMatrix<FAD,numdim_,numnodescont_> ele2posref;
  for(unsigned int i = 0; i < numdim_; i++)
    for(unsigned int j = 0; j < numnodescont_; j++)
      ele2posref(i,j).val() = ele2posref_(i,j);
  // ele2pos
  static LINALG::TMatrix<FAD,numdim_,numnodescont_> ele2pos;
  for(unsigned int i = 0; i < numdim_; i++)
    for(unsigned int j = 0; j < numnodescont_; j++)
      ele2pos(i,j).val() = ele2pos_(i,j);

  const FAD determinant = (eta_s - eta_a_)/2.0;
  const FAD jacobi      = determinant*arteryelelengthref_/2.0;

  // integrate from etaA to eta_s
  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    const double w_gp = wgp_[i_gp];
    const FAD eta = (eta_s + eta_a_)/2.0 + gaussPoints.qxg[i_gp][0] * determinant;

    // project
    bool projection_valid = false;
    std::vector<FAD> xi(numdim_, 0.0);
    Projection<FAD>(eta, xi, projection_valid);
    if(!projection_valid)
      dserror("Gauss point could not be projected");

    Get2D3DShapeFunctions<FAD>(N2, N2_xi, xi);

    InvJ.MultiplyNT(N2_xi,ele2posref);
    InvJ.Invert();

    // dN/dX = dN/dxi * dxi/dX = dN/dxi * (dX/dxi)^-1
    N2_XYZ.Multiply(InvJ,N2_xi);
    // dx/dX = x * N_XYZ^T
    defGrad.MultiplyNT(ele2pos,N2_XYZ);
    Ft0.Multiply(defGrad,t0);
    const FAD Ft0Norm = FADUTILS::VectorNorm<numdim_>(Ft0);
    // finally get the length
    length += Ft0Norm*w_gp*jacobi;
  }

  return length;
}

/*----------------------------------------------------------------------*
 | get artery values (pressure and scalars) at GP      kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::GetArteryValuesAtGP(
    const LINALG::Matrix<1, numnodesart_>& N1,
    double&                               artpress,
    std::vector<double>&                  artscalar
    )
{

  switch(coupltype_)
  {
  case type_porofluid:
  {
    for(unsigned int i = 0; i < numnodesart_; i++)
    {
      artpress += N1(i)*artelephinp_[i];
      for(int i_scal = 0; i_scal < numscalart_; i_scal++)
        artscalar[i_scal] += N1(i)*eartscalarnp_[i_scal](i);
    }
    break;
  }
  case type_scatra:
  {
    for(unsigned int i = 0; i < numnodesart_; i++)
    {
      artpress += N1(i)*earterypressurenp_(i);
      for(int i_art = 0; i_art < numdof_art_; i_art++)
        artscalar[i_art] += N1(i)*artelephinp_[i*numdof_art_+i_art];
    }
    break;
  }
  default:
    dserror("Unknown coupling type");
    break;
  }

  return;
}

/*----------------------------------------------------------------------*
 | get scalar values (continuous dis) at GP            kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::GetContScalarValuesAtGP(
    const LINALG::Matrix<1, numnodescont_>& N2,
    std::vector<double>&                   contscalarnp
    )
{

  switch(coupltype_)
  {
  case type_porofluid:
  {
    for(unsigned int i = 0; i < numnodescont_; i++)
    {
      for(int i_cont = 0; i_cont < numscalcont_; i_cont++)
        contscalarnp[i_cont] += N2(i)*econtscalarnp_[i_cont](i);
    }
    break;
  }
  case type_scatra:
  {
    for(unsigned int i = 0; i < numnodescont_; i++)
    {
      for(int i_cont = 0; i_cont < numdof_cont_; i_cont++)
        contscalarnp[i_cont] += N2(i)*contelephinp_[i*numdof_cont_+i_cont];
    }
    break;
  }
  default:
    dserror("Unknown coupling type");
    break;
  }

  return;
}

/*----------------------------------------------------------------------*
 | assemble into element force and stiffness matrix (functions of       |
 | artery discretization)                              kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::AssembleFunctionCouplingIntoForceStiffArt(
    const int& i_art,
    const double& w_gp,
    const LINALG::Matrix<1, numnodesart_>& N1,
    const LINALG::Matrix<1, numnodescont_>& N2,
    const double& jacobi,
    const int& scale,
    const double& functval,
    const std::vector<double>& artderivs,
    const std::vector<double>& contderivs,
    LINALG::SerialDenseVector& forcevec1,
    LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12
    )
{
  const double myscale = (double)(scale);

  // rhs ---> +
  for (unsigned int i = 0; i < numnodesart_; i++)
    forcevec1(i*numdof_art_+i_art) += N1(i)*myscale*w_gp*jacobi*functval*timefacrhs_art_;

  // matrix --> -
  for (unsigned int i = 0; i < numnodesart_; i++)
    for (unsigned int j = 0; j < numnodesart_; j++)
      for(int j_art = 0; j_art < numdof_art_; j_art++)
        stiffmat11(i*numdof_art_+i_art, j*numdof_art_+j_art) -= N1(i)*N1(j)*myscale*w_gp*jacobi*artderivs[j_art]*timefacrhs_art_;

  // matrix --> -
  for (unsigned int i = 0; i < numnodesart_; i++)
    for (unsigned int j = 0; j < numnodescont_; j++)
      for(int j_cont = 0; j_cont < numdof_cont_; j_cont++)
        stiffmat12(i*numdof_art_+i_art, j*numdof_cont_+j_cont) -= N1(i)*N2(j)*myscale*w_gp*jacobi*contderivs[j_cont]*timefacrhs_art_;

  return;
}

/*----------------------------------------------------------------------*
 | assemble into element force and stiffness matrix (functions of       |
 | continuous discretization)                          kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::AssembleFunctionCouplingIntoForceStiffCont(
    const int& i_cont,
    const double& w_gp,
    const LINALG::Matrix<1, numnodesart_>& N1,
    const LINALG::Matrix<1, numnodescont_>& N2,
    const double& jacobi,
    const int& scale,
    const double& functval,
    const std::vector<double>& artderivs,
    const std::vector<double>& contderivs,
    LINALG::SerialDenseVector& forcevec2,
    LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22
    )
{
  const double myscale = (double)(scale);

  // rhs ---> +
  for (unsigned int i = 0; i < numnodescont_; i++)
    forcevec2(i*numdof_cont_+i_cont) += N2(i)*myscale*w_gp*jacobi*functval*timefacrhs_cont_;

  // matrix --> -
  for (unsigned int i = 0; i < numnodescont_; i++)
    for (unsigned int j = 0; j < numnodesart_; j++)
      for(int j_art = 0; j_art < numdof_art_; j_art++)
        stiffmat21(i*numdof_cont_+i_cont, j*numdof_art_+j_art) -= N2(i)*N1(j)*myscale*w_gp*jacobi*artderivs[j_art]*timefacrhs_cont_;

  // matrix --> -
  for (unsigned int i = 0; i < numnodescont_; i++)
    for (unsigned int j = 0; j < numnodescont_; j++)
      for(int j_cont = 0; j_cont < numdof_cont_; j_cont++)
        stiffmat22(i*numdof_cont_+i_cont, j*numdof_cont_+j_cont) -= N2(i)*N2(j)*myscale*w_gp*jacobi*contderivs[j_cont]*timefacrhs_cont_;

  return;
}

/*----------------------------------------------------------------------*
 | evaluate function values and derivatives            kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateFunctionAndDeriv(
    DRT::UTILS::VariableExprFunction* funct,
    const double& artpressnpAtGP,
    const std::vector<double>& artscalarnpAtGP,
    const std::vector<double>& scalarnpAtGP,
    double& functval,
    std::vector<double>& artderivs,
    std::vector<double>& contderivs)
{

  switch(coupltype_)
  {
  case type_porofluid:
  {
    // we have to derive w.r.t. fluid variables
    std::vector<std::pair<std::string,double> > variables;
    variables.reserve(numfluidphases_+numfluidphases_+1+numvolfrac_+numvolfrac_+1);

    // scalar variables are constants
    std::vector<std::pair<std::string,double> > constants;
    constants.reserve(numscalcont_+numscalart_+1);

    SetScalarValuesAsConstants(constants, artscalarnpAtGP, scalarnpAtGP);

    SetFluidValuesAsVariables(variables, artpressnpAtGP);

    // evaluate the reaction term
    functval = funct->Evaluate(0,variables,constants);
    // evaluate derivatives
    std::vector<double> curderivs(funct->EvaluateDerivative(0,variables,constants));

    EvaluateFluidDerivs(artderivs, contderivs, curderivs);

    break;
  }
  case type_scatra:
  {
    // scalars (both cont and art) are variables
    std::vector<std::pair<std::string,double> > variables;
    variables.reserve(numscalcont_+numscalart_);

    // fluid variables are constants
    std::vector<std::pair<std::string,double> > constants;
    constants.reserve(numfluidphases_+numfluidphases_+1+numvolfrac_+numvolfrac_+1+1);

    SetScalarValuesAsVariables(variables, artscalarnpAtGP, scalarnpAtGP);

    SetFluidValuesAsConstants(constants, artpressnpAtGP);

    // evaluate the reaction term
    functval = funct->Evaluate(0,variables,constants);
    // evaluate derivatives
    std::vector<double> curderivs(funct->EvaluateDerivative(0,variables,constants));

    EvaluateScalarDerivs(artderivs, contderivs, curderivs);

    break;
  }
  default:
    dserror("Unknown coupling type");
    break;
  }

  return;
}

/*----------------------------------------------------------------------*
 | set scalar values as constants for function         kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::SetScalarValuesAsConstants(
    std::vector<std::pair<std::string,double> >& constants,
    const std::vector<double>& artscalarnpAtGP,
    const std::vector<double>& scalarnpAtGP)
{

  // set scalar values as constant
  for (int k=0;k<numscalcont_;k++)
    constants.push_back(std::pair<std::string,double>(scalarnames_[k],scalarnpAtGP[k]));

  // set artery-scalar values as constant
  for (int k=0;k<numscalart_;k++)
    constants.push_back(std::pair<std::string,double>(artscalarnames_[k],artscalarnpAtGP[k]));

  // set artery diameter as constant
  constants.push_back(std::pair<std::string,double>("D",arterydiam_));

  return;
}

/*----------------------------------------------------------------------*
 | set fluid values as variables for function          kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::SetFluidValuesAsVariables(
    std::vector<std::pair<std::string,double> >& variables,
    const double& artpressnpAtGP)
{

  // set pressure values as variable
  for (int k=0;k<numfluidphases_;k++)
    variables.push_back(std::pair<std::string,double>(pressurenames_[k],phasemanager_->Pressure(k)));

  // set saturation values as variable
  for (int k=0;k<numfluidphases_;k++)
    variables.push_back(std::pair<std::string,double>(saturationnames_[k],phasemanager_->Saturation(k)));

  // set porosity value as variable
  variables.push_back(std::pair<std::string,double>(porosityname_,phasemanager_->Porosity()));

  // set volfrac values as variables
  for (int k=0;k<numvolfrac_;k++)
    variables.push_back(std::pair<std::string,double>(volfracnames_[k],phasemanager_->VolFrac(k)));

  // set volfrac pressure values as variables
  for (int k=0;k<numvolfrac_;k++)
    variables.push_back(std::pair<std::string,double>(volfracpressurenames_[k],phasemanager_->VolFracPressure(k)));

  // set artery pressure as variable
  variables.push_back(std::pair<std::string,double>(artpressname_,artpressnpAtGP));

  return;
}

/*----------------------------------------------------------------------*
 | set fluid values as variables for function          kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::SetFluidValuesAsConstants(
    std::vector<std::pair<std::string,double> >& constants,
    const double& artpressnpAtGP)
{

  // set pressure values as constants
  for (int k=0;k<numfluidphases_;k++)
    constants.push_back(std::pair<std::string,double>(pressurenames_[k],phasemanager_->Pressure(k)));

  // set saturation values as constants
  for (int k=0;k<numfluidphases_;k++)
    constants.push_back(std::pair<std::string,double>(saturationnames_[k],phasemanager_->Saturation(k)));

  // set porosity value as constants
  constants.push_back(std::pair<std::string,double>(porosityname_,phasemanager_->Porosity()));

  // set volfrac values as constants
  for (int k=0;k<numvolfrac_;k++)
    constants.push_back(std::pair<std::string,double>(volfracnames_[k],phasemanager_->VolFrac(k)));

  // set volfrac pressure values as constants
  for (int k=0;k<numvolfrac_;k++)
    constants.push_back(std::pair<std::string,double>(volfracpressurenames_[k],phasemanager_->VolFracPressure(k)));

  // set artery pressure as constant
  constants.push_back(std::pair<std::string,double>(artpressname_,artpressnpAtGP));

  // set artery diameter as constant
  constants.push_back(std::pair<std::string,double>("D",arterydiam_));

  return;
}

/*----------------------------------------------------------------------*
 | set scalar values as variables for function         kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::SetScalarValuesAsVariables(
    std::vector<std::pair<std::string,double> >& variables,
    const std::vector<double>& artscalarnpAtGP,
    const std::vector<double>& scalarnpAtGP)
{

  // set scalar values as variables
  for (int k=0;k<numscalcont_;k++)
    variables.push_back(std::pair<std::string,double>(scalarnames_[k],scalarnpAtGP[k]));

  // set artery-scalar values as variables
  for (int k=0;k<numscalart_;k++)
    variables.push_back(std::pair<std::string,double>(artscalarnames_[k],artscalarnpAtGP[k]));

  return;
}

/*----------------------------------------------------------------------*
 | build the fluid derivatives                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateFluidDerivs(
    std::vector<double>&       artderivs,
    std::vector<double>&       contderivs,
    const std::vector<double>& functderivs)
{

  // function derivs w.r.t. to fluid phases
  for(int doftoderive=0; doftoderive<numfluidphases_; doftoderive++)
  {
    for(int idof=0; idof<numfluidphases_; idof++)
      contderivs[doftoderive] +=   functderivs[idof]*phasemanager_->PressureDeriv(idof,doftoderive)
                                 + functderivs[idof+numfluidphases_]*phasemanager_->SaturationDeriv(idof,doftoderive);
    if(phasemanager_->PorosityDependsOnFluid())
      contderivs[doftoderive] += functderivs[2*numfluidphases_]*phasemanager_->PorosityDeriv(doftoderive);
  }
  // function derivs w.r.t. to volume fraction phases
  for(int doftoderive=numfluidphases_; doftoderive<numdof_cont_-numvolfrac_; doftoderive++)
  {
    // derivatives w.r.t. volume fractions directly appearing
    //                and porosity (since it depends on volfrac)
    contderivs[doftoderive] += functderivs[doftoderive+numfluidphases_+1]
                             + functderivs[2*numfluidphases_]*phasemanager_->PorosityDeriv(doftoderive);
  }
  // function derivs w.r.t. to volume fraction pressures
  for(int doftoderive=numfluidphases_+numvolfrac_; doftoderive<numdof_cont_; doftoderive++)
    contderivs[doftoderive] += functderivs[doftoderive+numfluidphases_+1];

  // function derivs w.r.t. to artery pressure
  artderivs[0] += functderivs[numfluidphases_+numfluidphases_+1+numvolfrac_+numvolfrac_];

  return;
}

/*----------------------------------------------------------------------*
 | build the scalar derivatives                        kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::EvaluateScalarDerivs(
    std::vector<double>&       artderivs,
    std::vector<double>&       contderivs,
    const std::vector<double>& functderivs)
{

  // derivatives after continuous scalars
  for(int doftoderive=0; doftoderive<numdof_cont_; doftoderive++)
    contderivs[doftoderive] += functderivs[doftoderive];

  // derivatives after artery scalars
  for(int doftoderive=0; doftoderive<numdof_art_; doftoderive++)
    artderivs[doftoderive] += functderivs[doftoderive+numdof_cont_];

  return;
}
/*----------------------------------------------------------------------*
 | create integration segment [eta_a, eta_b]           kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::CreateIntegrationSegment(
    const std::vector<bool>& validprojections)
{

  if(PROJOUTPUT)
    std::cout << "\nStarting with creation of integration segments ===========================================\n" << std::endl;
  bool allprojectionsvalid = true;
  const int numproj = validprojections.size();
  for(int i = 0; i < numproj; i++)
  {
    if(!validprojections[i])
    {
      allprojectionsvalid = false;
      break;
    }
  }

  // all points could be projected, the 1D element lies completely inside
  if(allprojectionsvalid)
  {
    eta_a_ = -1.0;
    eta_b_ =  1.0;
    if(PROJOUTPUT)
      std::cout << "all could be projected" << std::endl;
  }
  else
  {
    // the first few could be projected, segment will look like [-1.0, etaB]
    if(validprojections[0] && !validprojections[numproj-1])
    {
      if(PROJOUTPUT)
        std::cout << "first few could be projected" << std::endl;
      eta_a_ = -1.0;
      std::vector<double> intersections = GetAllInterSections();
      if(intersections.size() < 1)
        dserror("found not enough intersections");
      eta_b_ = *std::max_element(intersections.begin(), intersections.end());
    }
    // last few could be projected, segment will look like [etaA, 1.0]
    else if(validprojections[numproj-1] && !validprojections[0])
    {
      if(PROJOUTPUT)
        std::cout << "last few could be projected" << std::endl;
      eta_b_ = 1.0;
      std::vector<double> intersections = GetAllInterSections();
      if(intersections.size() < 1)
        dserror("found not enough intersections");
      eta_a_ = *std::min_element(intersections.begin(), intersections.end());
    }
    // middle few could be projected, segment will look like [etaA, etaB]
    else if(!validprojections[0] && !validprojections[numproj-1])
    {
      if(PROJOUTPUT)
        std::cout << "middle few could be projected" << std::endl;
      std::vector<double> intersections = GetAllInterSections();
      if(intersections.size() != 2)
        dserror("found too many or less intersections");
      eta_a_ = intersections[0];
      eta_b_ = intersections[1];
    }
    else
      dserror("this should not happen");
  }

  if(PROJOUTPUT)
  {
    std::cout << "Finished with creation of integration segments ===========================================" << std::endl;
    std::cout << "eta_a: " << eta_a_ << ", " << "eta_b: " << eta_b_ << std::endl;
  }

  // safety checks
  if(eta_a_ > eta_b_)
    dserror("something went terribly wrong, eta_a is bigger than eta_b");
  if(fabs(eta_a_-eta_b_) < SMALLESTSEGMENT)
    dserror("something went terribly wrong, found extremely small integration segment");

  return;
}

/*-----------------------------------------------------------------------------*
 | get all intersections of artery-element with 2D/3D-element kremheller 05/18 |
 *-----------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
std::vector<double> POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::GetAllInterSections()
{

  std::vector<double> intersections(0);

  std::vector<double> xi(numdim_,0.0);
  double eta = 0.0;
  bool projection_valid = true;

  if(numdim_ == 2)
  {
    if (numnodescont_ == 4)
    {
      for(unsigned int j = 0; j < numdim_; j++)
      {
        // project edge xi1 or xi2 = 1.0
        InterSectWith2D3D(xi, eta, j, 1.0, projection_valid);
        if (projection_valid && ProjectionNotYetFound(intersections, eta))
          intersections.push_back(eta);

        // project edge xi1 or xi2 = -1.0
        InterSectWith2D3D(xi, eta, j, -1.0, projection_valid);
        if (projection_valid && ProjectionNotYetFound(intersections, eta))
          intersections.push_back(eta);
      }
    }
    else
      dserror("only quad4 ele valid in 2D so far");
  }
  else if (numdim_ == 3)
  {
    if (numnodescont_ == 8)
    {
      for(unsigned int j = 0; j < numdim_; j++)
      {
        // project edge xi1 or xi2 or xi3 = 1.0
        InterSectWith2D3D(xi, eta, j, 1.0, projection_valid);
        if (projection_valid && ProjectionNotYetFound(intersections, eta))
          intersections.push_back(eta);

        // project edge xi1 or xi2 or xi3 = -1.0
        InterSectWith2D3D(xi, eta, j, -1.0, projection_valid);
        if (projection_valid && ProjectionNotYetFound(intersections, eta))
          intersections.push_back(eta);
      }
    }
    else
      dserror("only hex8 ele valid in 3D so far");
  }
  else
    dserror("Only numdim_ = 2, 3 is valid");

  std::sort(intersections.begin(), intersections.end());

  return intersections;
}

/*-----------------------------------------------------------------------------*
 | check for duplicate projections                            kremheller 05/18 |
 *-----------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
bool POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::ProjectionNotYetFound(
    const std::vector<double>& intersections,
    const double&              eta)
{

  for(unsigned int i = 0; i < intersections.size(); i++)
  {
    if(fabs(intersections[i]-eta) < ETAABTOL)
    {
      if(PROJOUTPUT)
        std::cout << "duplicate intersection found" << std::endl;
      return false;
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 | intersect with boundary of 2D/3D-element            kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::InterSectWith2D3D(
    std::vector<double>& xi,
    double&              eta,
    const int&           fixedPar,
    const double&        fixedAt,
    bool&                projection_valid
    )
{

  projection_valid = true;
  bool parallel = false;

  // Initialize limit for parameter values (interval [-limit, limit])
  const double limit = 1.0 + VALIDPROJTOL;

  // reset iteration variables
  eta = 0.0;
  if(numdim_ == 2)
  {
    if (numnodescont_ == 4)
    {
      if (fixedPar == 0) // xi1 fixed
      {
        xi[0] = fixedAt;
        xi[1] = 0.0;
      }
      else if (fixedPar == 1) // xi2 fixed
      {
        xi[0] = 0.0;
        xi[1] = fixedAt;
      }
      else
        dserror("wrong input for fixedPar");
    }
    else
      dserror("only quad4 ele valid in 2D so far");
  }
  else if (numdim_ == 3)
  {
    if (numnodescont_ == 8)
    {
      if (fixedPar == 0) // xi1 fixed
      {
        xi[0] = fixedAt;
        xi[1] = 0.0;
        xi[2] = 0.0;
      }
      else if (fixedPar == 1) // xi2 fixed
      {
        xi[0] = 0.0;
        xi[1] = fixedAt;
        xi[2] = 0.0;
      }
      else if (fixedPar == 2) // xi3 fixed
      {
        xi[0] = 0.0;
        xi[1] = 0.0;
        xi[2] = fixedAt;
      }
      else
        dserror("wrong input for fixedPar");
    }
    else
      dserror("only hex8 ele valid in 3D so far");
  }
  else
    dserror("Only numdim_ = 2, 3 is valid");

  if (PROJOUTPUT)
  {
    std::cout << "Projection output:" << std::endl;
    std::cout << "Start parameters eta: " << eta;
    if(numdim_ == 2)
      std::cout << ", xi1: " << xi[0] << ", xi2: " << xi[1] << std::endl;
    else if(numdim_ == 3)
      std::cout << ", xi1: " << xi[0] << ", xi2: " << xi[1] << ", xi3: "<< xi[2] << std::endl;

  }

  // Initialize function f and Jacobian J for Newton iteration
  LINALG::Matrix<numdim_,1> f(true);
  LINALG::Matrix<numdim_,numdim_> J(true);
  LINALG::Matrix<numdim_,numdim_> Jinv(true);

  // Vectors for shape functions and their derivatives
  static LINALG::Matrix<1, numnodesart_> N1(true);         // = N1
  static LINALG::Matrix<1, numnodesart_> N1_eta(true);     // = N1,eta

  static LINALG::Matrix<1, numnodescont_> N2(true);                     // = N2
  static LINALG::Matrix<numdim_, numnodescont_> N2_xi(true);             // = N2,xi1

  // Coords and derivatives of for 1D and 2D/3D element
  static LINALG::Matrix<numdim_, 1> r1(true);                                 // = r1
  static LINALG::Matrix<numdim_, 1> r1_eta(true);                             // = r1,eta

  static LINALG::Matrix<numdim_, 1> x2(true);                                 // = x2
  static LINALG::Matrix<numdim_, numdim_> x2_xi(true);                   // = x2,xi

  // Initial scalar residual (L2-norm of f)
  double residual;

  // Local newton iteration
  // -----------------------------------------------------------------
  int iter;

  for (iter = 0; iter < PROJMAXITER; iter++)
  {
    // Update shape functions and their derivatives for 1D and 2D/3D element
    Get1DShapeFunctions<double>(N1, N1_eta, eta);
    Get2D3DShapeFunctions<double>(N2, N2_xi, xi);

    // Update coordinates and derivatives for 1D and 2D/3D element
    ComputeArteryCoordsAndDerivsRef<double>(r1, r1_eta, N1, N1_eta);
    Compute2D3DCoordsAndDerivsRef<double>(x2, x2_xi, N2, N2_xi);

    // Evaluate f at current xi1, xi2, alpha
    f.Clear();
    for (unsigned int i = 0; i < numdim_; i++)
      f(i) = x2(i) - r1(i);

    // Compute scalar residuum
    residual = 0.0;
    for (unsigned int i = 0; i < numdim_; i++)
      residual += f(i)*f(i);
    residual = sqrt(residual);

    J.Clear();

    if (fixedPar == 0) // xi1 fixed --> we need x_{,xi2} (and x_{,xi3} in case of 3D)
    {
      for(unsigned int idim = 0; idim < numdim_; idim++)
        for(unsigned int jdim = 1; jdim < numdim_; jdim++)
          J(idim,jdim-1) = x2_xi(idim, jdim);
    }
    else if (fixedPar == 1) // xi2 fixed --> we need x_{,xi1} (and x_{,xi3} in case of 3D)
    {
      if(numdim_ == 2)
      {
        for(unsigned int jdim = 0; jdim < numdim_; jdim++)
          J(jdim,0) = x2_xi(jdim, 0);
      }
      else if(numdim_ == 3)
      {
        for(unsigned int jdim = 0; jdim < numdim_; jdim++)
        {
          J(jdim,0) = x2_xi(jdim, 0);
          J(jdim,1) = x2_xi(jdim, 2);
        }
      }
    }
    else if (fixedPar == 2) // xi3 fixed  --> we need x_{,xi1} (and x_{,xi2} in case of 3D)
    {
      for(unsigned int idim = 0; idim < numdim_; idim++)
        for(unsigned int jdim = 0; jdim < numdim_-1; jdim++)
          J(idim,jdim) = x2_xi(idim, jdim);
    }

    // fill dr_deta into Jacobian
    for(unsigned int idim = 0; idim < numdim_; idim++)
      J(idim, numdim_-1) = -r1_eta(idim);

    double jacdet = J.Determinant();

    // If det_J = 0 we assume, that the artery and the surface edge are parallel.
    // These projection is not needed due the fact that the contact interval can also be
    // identified by other projections
    parallel = fabs(jacdet) < COLINEARTOL;
    if (!parallel)
      jacdet = J.Invert();

    // Check if the local Newton iteration has converged
    if (residual < CONVTOLNEWTONPROJ && !parallel)
    {
      if (PROJOUTPUT)
      {
        std::cout << "Local Newton iteration converged after " << iter << " iterations" << std::endl;
        std::cout << "Found point at xi1: " << xi[0] << ", xi2: " << xi[1];
        if(numdim_ == 3)
          std::cout << ", xi3: " << xi[2];
        std::cout << ", eta: " << eta << " with residual: " << residual << std::endl;
        std::cout << "r1:\n" << r1 << ", x2:\n" << x2 << std::endl;
      }
      // Local Newton iteration converged
      break;
    }
    else if (PROJOUTPUT && iter > 0)
    {
      std::cout << "New point at xi1: " << xi[0] << ", xi2: " << xi[0];
      if(numdim_ == 3)
        std::cout << ", xi3: " << xi[2];
      std::cout << ", eta: " << eta << " with residual: " << residual <<  std::endl;
    }

    // Singular J
    if (parallel)
    {
      // Sort out
      if (PROJOUTPUT)
      {
        std::cout << "elementscolinear: det_J = " << jacdet << std::endl;
      }
      break;
    }
    // Regular J (inversion possible)

    if (fixedPar == 0) // xi1 fixed --> we have to update xi2 (and xi3 in case of 3D)
    {
      for(unsigned int idim = 1; idim < numdim_; idim++)
        for(unsigned int jdim = 0; jdim < numdim_; jdim++)
          xi[idim] += -J(idim-1, jdim) * f(jdim);
    }
    else if (fixedPar == 1) // xi2 fixed --> we have to update xi1 (and xi3 in case of 3D)
    {
      if(numdim_ == 2)
      {
        for(unsigned int jdim = 0; jdim < numdim_; jdim++)
          xi[0] += -J(0, jdim) * f(jdim);
      }
      else if(numdim_ == 3)
      {
        for(unsigned int jdim = 0; jdim < numdim_; jdim++)
        {
          xi[0] += -J(0, jdim) * f(jdim);
          xi[2] += -J(1, jdim) * f(jdim);
        }
      }
    }
    else if (fixedPar == 2)  // xi3 fixed --> we have to update xi1 (and xi2 in case of 3D)
    {
      for(unsigned int idim = 0; idim < numdim_-1; idim++)
        for(unsigned int jdim = 0; jdim < numdim_; jdim++)
          xi[idim] += -J(idim, jdim)*f(jdim);
    }
    //else if (fixedPar == 3)
    //{
    //  xi1 += -J(0, 0) * f(0) - J(0, 1) * f(1) - J(0, 2) * f(2);
    //  xi2 += -J(1, 0) * f(0) - J(1, 1) * f(1) - J(1, 2) * f(2);
    //  xi3 = 1.0 - xi1 - xi2;
    //}

    // update also eta
    for(unsigned int jdim = 0; jdim < numdim_; jdim++)
      eta += -J(numdim_-1, jdim) * f(jdim);

  }

  // Local Newton iteration unconverged after PROJMAXITER
  if (residual > CONVTOLNEWTONPROJ || parallel)
  {

    for(unsigned int idim = 0; idim < numdim_; idim++)
      xi[idim] = 1e+12;
    eta = 1e+12;

    if (PROJOUTPUT)
      std::cout << "Local Newton iteration unconverged (!) after " << iter + 1 << " iterations" << std::endl;
  }

  if(numdim_ == 2)
  {
    if (numnodescont_ == 4)
    {
      if (fabs(xi[0]) > limit || fabs(xi[1]) > limit || fabs(eta) > limit)
        projection_valid = false;
    }
    else
      dserror("only quad4 ele valid in 2D so far");
  }
  else if (numdim_ == 3)
  {
    if (numnodescont_ == 8)
    {
      if (fabs(xi[0]) > limit || fabs(xi[1]) > limit || fabs(xi[2]) > limit || fabs(eta) > limit)
        projection_valid = false;
    }
    else
      dserror("only hex8 ele valid in 3D so far");
  }
  else
    dserror("Only numdim_ = 2, 3 is valid");

  if(PROJOUTPUT)
  {
    if(projection_valid)
      std::cout << "Projection allowed" << std::endl;
    else
      std::cout << "Projection not allowed" << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 | project into 2D/3D-element                          kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>template<typename T>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::Projection(
    const T& eta,
    std::vector<T>& xi,
    bool& projection_valid)
{

  projection_valid = true;
  bool parallel = false;

  // Initialize limit for parameter values (interval [-limit, limit])
  const double limit = 1.0 + VALIDPROJTOL;

  if(numdim_ == 2)
  {
    if (numnodescont_ == 4)
    {
      xi[0] = 0.0;
      xi[1] = 0.0;
    }
    else
      dserror("only quad4 ele valid in 2D so far");
  }
  else if (numdim_ == 3)
  {
    if (numnodescont_ == 8)
    {
      xi[0] = 0.0;
      xi[1] = 0.0;
      xi[2] = 0.0;
    }
    else
      dserror("only hex8 ele valid in 3D so far");
  }
  else
    dserror("Only numdim_ = 2, 3 is valid");

  if (PROJOUTPUT)
  {
    std::cout << "Projection output:" << std::endl;
    std::cout << "Start parameters ";
    if(numdim_ == 2)
      std::cout << "xi1: " << xi[0] << ", xi2: " << xi[1] << std::endl;
    else if(numdim_ == 3)
      std::cout << "xi1: " << xi[0] << ", xi2: " << xi[1] << ", xi3: " << xi[2] << std::endl;
  }

  // Initialize function f and Jacobian J for Newton iteration
  LINALG::TMatrix<T,numdim_,1> f(true);
  LINALG::TMatrix<T,numdim_,numdim_> J(true);
  LINALG::TMatrix<T,numdim_,numdim_> Jinv(true);

  // Vectors for shape functions and their derivatives
  static LINALG::TMatrix<T, 1, numnodesart_> N1(true);         // = N1
  static LINALG::TMatrix<T, 1, numnodesart_> N1_eta(true);     // = N1,eta

  static LINALG::TMatrix<T, 1, numnodescont_> N2(true);                     // = N2
  static LINALG::TMatrix<T, numdim_, numnodescont_> N2_xi(true);             // = N2,xi1

  // Coords and derivatives of 1D and 2D/3D element
  static LINALG::TMatrix<T, numdim_, 1> r1(true);                                 // = r1
  static LINALG::TMatrix<T, numdim_, 1> r1_eta(true);                             // = r1,eta

  static LINALG::TMatrix<T, numdim_, 1> x2(true);                                 // = x2
  static LINALG::TMatrix<T, numdim_, numdim_> x2_xi(true);                   // = x2,xi

  // Initial scalar residual (L2-norm of f)
  T residual;

  // Local newton iteration
  // -----------------------------------------------------------------

  int iter;

  for (iter = 0; iter < PROJMAXITER; iter++)
  {
    // Update shape functions and their derivatives for 1D and 2D/3D element
    Get1DShapeFunctions<T>(N1, N1_eta, eta);
    Get2D3DShapeFunctions<T>(N2, N2_xi, xi);

    // Update coordinates and derivatives for 1D and 2D/3D element
    ComputeArteryCoordsAndDerivsRef<T>(r1, r1_eta, N1, N1_eta);
    Compute2D3DCoordsAndDerivsRef<T>(x2, x2_xi, N2, N2_xi);

    // Evaluate f at current xi1, xi2, alpha
    f.Clear();
    for (unsigned int i = 0; i < numdim_; i++)
      f(i) = x2(i) - r1(i);

    residual = FADUTILS::VectorNorm<numdim_>(f);

    // Reset matrices
    for (unsigned int i = 0; i < numdim_; i++)
      for (unsigned int j = 0; j < numdim_; j++)
        J(i,j) = x2_xi(i,j);

    const double jacdet = FADUTILS::CastToDouble<T,numdim_,numdim_>(J).Determinant();

    // If det_J = 0 we assume, that the artery element and the surface edge are parallel.
    // These projection is not needed due the fact that the contact interval can also be
    // identified by two contact interval borders found with the GetContactLines method
    parallel = fabs(jacdet) < COLINEARTOL;
    if (!parallel)
      J.Invert();

    // Check if the local Newton iteration has converged
    // If the start point fulfills the orthogonalty conditions (residual < CONVTOLNEWTONPROJ), we also check if
    // the artery element and the surface edge are parallel. This is done by calculating det_J before checking
    // if the local Newton iteration has converged by fulfilling the condition residual < CONVTOLNEWTONPROJ
    if (residual < CONVTOLNEWTONPROJ && !parallel)
    {
      if (PROJOUTPUT)
      {
        std::cout << "Local Newton iteration converged after " << iter << " iterations" << std::endl;
        std::cout << "Found point at ";
        if(numdim_ == 2)
          std::cout << "xi1: " << xi[0] << ", xi2: " << xi[1] << std::endl;
        else if(numdim_ == 3)
          std::cout << "xi1: " << xi[0] << ", xi2: " << xi[1] << ", xi3: " << xi[2] << std::endl;
        std::cout <<  " with residual: " << residual << std::endl;
        std::cout << "r1:\n" << r1 << ", x2:\n" << x2 << std::endl;
      }
      // Local Newton iteration converged
      break;
    }
    else if (PROJOUTPUT && iter > 0)
    {
      std::cout << "New point at xi1: ";
      if(numdim_ == 2)
        std::cout << "xi1: " << xi[0] << ", xi2: " << xi[1] << std::endl;
      else if(numdim_ == 3)
        std::cout << "xi1: " << xi[0] << ", xi2: " << xi[1] << ", xi3: " << xi[2] << std::endl;
      std::cout <<  " with residual: " << residual << std::endl;
    }

    // Singular J
    if (parallel)
    {
      // Sort out
      if (PROJOUTPUT)
      {
        std::cout << "elementscolinear: det_J = " << jacdet << std::endl;
      }
      break;
    }
    // Regular J (inversion possible)
    for(unsigned int idim = 0; idim < numdim_; idim++)
      for(unsigned int jdim = 0; jdim < numdim_; jdim++)
        xi[idim] += -J(idim, jdim) * f(jdim);

    //xi1 += -J(0, 0) * f(0) - J(0, 1) * f(1) - J(0, 2) * f(2);
    //xi2 += -J(1, 0) * f(0) - J(1, 1) * f(1) - J(1, 2) * f(2);
    //xi3 += -J(2, 0) * f(0) - J(2, 1) * f(1) - J(2, 2) * f(2);

  }
  // -----------------------------------------------------------------
  // End: Local Newton iteration

  // Local Newton iteration unconverged after PROJMAXITER
  if (residual > CONVTOLNEWTONPROJ || parallel)
  {

    for(unsigned int idim = 0; idim < numdim_; idim++)
      xi[idim] = 1e+12;

    if (PROJOUTPUT)
      std::cout << "Local Newton iteration unconverged (!) after " << iter + 1 << " iterations" << std::endl;
  }

  if(numdim_ == 2)
  {
    if (numnodescont_ == 4)
    {
      if (fabs(xi[0]) > limit || fabs(xi[1]) > limit)
        projection_valid = false;
    }
    else
      dserror("only quad4 ele valid in 2D so far");
  }
  else if (numdim_ == 3)
  {
    if (numnodescont_ == 8)
    {
      if (fabs(xi[0]) > limit || fabs(xi[1]) > limit || fabs(xi[2]) > limit)
        projection_valid = false;
    }
    else
      dserror("only hex8 ele valid in 3D so far");
  }
  else
    dserror("Only numdim_ = 2, 3 is valid");

  if(PROJOUTPUT)
  {
    if(projection_valid)
      std::cout << "Projection allowed" << std::endl;
    else
      std::cout << "Projection not allowed" << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 | get artery shape-function                           kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>template<typename T>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::Get1DShapeFunctions(
    LINALG::TMatrix<T, 1, numnodesart_>& N1,
    LINALG::TMatrix<T, 1, numnodesart_>& N1_eta,
    const T& eta)
{

  // Clear shape functions and derivatives
  N1.Clear();
  N1_eta.Clear();

  // Get discretization type
  const DRT::Element::DiscretizationType distype = element1_->Shape();

  // Get values and derivatives of shape functions
  DRT::UTILS::shape_function_1D(N1, eta, distype);
  DRT::UTILS::shape_function_1D_deriv1(N1_eta, eta, distype);

  return;
}

/*----------------------------------------------------------------------*
 | get shape-function of 2D/3D-element                 kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>template<typename T>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::Get2D3DShapeFunctions(
    LINALG::TMatrix<T, 1, numnodescont_>& N2,
    LINALG::TMatrix<T, numdim_, numnodescont_>& N2_xi,
    const std::vector<T>& xi)
{

  // Clear shape functions and derivatives
  N2.Clear();
  N2_xi.Clear();

  if(numdim_ == 2)
  {
    DRT::UTILS::shape_function_2D(N2, xi[0], xi[1], distypeCont);
    DRT::UTILS::shape_function_2D_deriv1(N2_xi, xi[0], xi[1], distypeCont);
  }
  else if(numdim_ == 3)
  {
    DRT::UTILS::shape_function_3D(N2, xi[0], xi[1], xi[2], distypeCont);
    DRT::UTILS::shape_function_3D_deriv1(N2_xi, xi[0], xi[1], xi[2], distypeCont);
  }
  else
    dserror("only numdim_ = 2,3 valid");

  return;
}

/*------------------------------------------------------------------------*
 | compute coordinates and derivatives of artery element kremheller 05/18 |
 *------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>template<typename T>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::ComputeArteryCoordsAndDerivsRef(
    LINALG::TMatrix<T, numdim_, 1>& r1,
    LINALG::TMatrix<T, numdim_, 1>& r1_eta,
    const LINALG::TMatrix<T, 1, numnodesart_>& N1,
    const LINALG::TMatrix<T, 1, numnodesart_>& N1_eta
    )
{

  r1.Clear();
  r1_eta.Clear();

  for (unsigned int j = 0; j < numnodesart_; j++)
  {
    for (unsigned int idim = 0; idim < numdim_; idim++)
    {
      r1(idim) += N1(j) * ele1posref_(numdim_*j+idim);
      r1_eta(idim) += N1_eta(j) * ele1posref_(numdim_*j+idim);
    }
  }
  return;
}

/*------------------------------------------------------------------------*
 | compute coordinates and derivatives of 2D/3D-element  kremheller 05/18 |
 *------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>template<typename T>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::Compute2D3DCoordsAndDerivsRef(
    LINALG::TMatrix<T, numdim_, 1>& x2,
    LINALG::TMatrix<T, numdim_, numdim_>& x2_xi,
    const LINALG::TMatrix<T, 1, numnodescont_>& N2,
    const LINALG::TMatrix<T, numdim_, numnodescont_>& N2_xi
    )
{

  x2.Clear();
  x2_xi.Clear();

  for (unsigned int j = 0; j < numnodescont_; j++)
  {
    for (unsigned int idim = 0; idim < numdim_; idim++)
    {
      x2(idim) += N2(j) * ele2posref_(idim, j);
      for (unsigned int jdim = 0; jdim < numdim_; jdim++)
      {
        x2_xi(idim, jdim) += N2_xi(jdim, j) * ele2posref_(idim, j);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | fill the function vector                            kremheller 05/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::FillFunctionVector(
    std::vector<DRT::UTILS::VariableExprFunction*>* my_funct_vec,
    const std::vector<int>& funct_vec,
    const std::vector<int>& scale_vec
    )
{

  for(unsigned int i = 0; i < funct_vec.size(); i++)
  {
    if(funct_vec[i] >= 0 && abs(scale_vec[i]) > 0)
    {
      my_funct_vec->at(i) = Function(funct_vec[i]);
      funct_coupl_active_ = true;
    }
    else
      my_funct_vec->at(i) = 0;
  }

  return;
}

/*----------------------------------------------------------------------*
 | initialize the functions                            kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::InitializeFunction(
   DRT::UTILS::VariableExprFunction* funct
   )
{

  // safety check
  if(funct->NumberComponents()!=1)
    dserror("expected only one component for coupling function!");

  for (int k=0;k<numscalcont_;k++)
  {
    // add scalar names
    if(not funct->IsVariable(0,scalarnames_[k]))
      funct->AddVariable(0,scalarnames_[k],0.0);
  }

  for (int k=0;k<numscalart_;k++)
  {
    // add artery-scalar names
    if(not funct->IsVariable(0,artscalarnames_[k]))
      funct->AddVariable(0,artscalarnames_[k],0.0);
  }

  for (int k=0;k<numfluidphases_;k++)
  {
    // add pressures
    if(not funct->IsVariable(0,pressurenames_[k]))
      funct->AddVariable(0,pressurenames_[k],0.0);
    // add saturations
    if(not funct->IsVariable(0,saturationnames_[k]))
      funct->AddVariable(0,saturationnames_[k],0.0);
  }

  // add porosity
  if(not funct->IsVariable(0,porosityname_))
    funct->AddVariable(0,porosityname_,0.0);

  // add additional volume fractions
  for (int k=0;k<numvolfrac_;k++)
  {
    // add volume fraction names
    if(not funct->IsVariable(0,volfracnames_[k]))
      funct->AddVariable(0,volfracnames_[k],0.0);
    // add volume fraction pressure names
    if(not funct->IsVariable(0,volfracpressurenames_[k]))
      funct->AddVariable(0,volfracpressurenames_[k],0.0);
  }

  // add artery-pressure
  if(not funct->IsVariable(0,artpressname_))
    funct->AddVariable(0,artpressname_,0.0);

  // add diameter
  if(not funct->IsVariable(0,"D"))
    funct->AddVariable(0,"D",0.0);

  return;
}

/*----------------------------------------------------------------------*
 | initialize the names used for function coupling     kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
void POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::InitializeFunctionNames()
{

  pressurenames_.resize(numfluidphases_);
  saturationnames_.resize(numfluidphases_);
  volfracnames_.resize(numvolfrac_);
  volfracpressurenames_.resize(numvolfrac_);
  scalarnames_.resize(numscalcont_);
  artscalarnames_.resize(numscalart_);

  for (int k=0;k<numscalcont_;k++)
  {
    // add scalar names
    {
      std::ostringstream temp;
      temp << k+1;
      scalarnames_[k] = "phi"+temp.str();
    }
  }

  for (int k=0;k<numscalart_;k++)
  {
    // add artery-scalar names
    {
      std::ostringstream temp;
      temp << k+1;
      artscalarnames_[k] = "phi_art"+temp.str();
    }
  }

  for (int k=0;k<numfluidphases_;k++)
  {
    // add pressure names
    {
      std::ostringstream temp;
      temp << k+1;
      pressurenames_[k] = "p"+temp.str();
    }

    // add saturation names
    {
      std::ostringstream temp;
      temp << k+1;
      saturationnames_[k] = "S"+temp.str();
    }
  }

  // add additional volume fractions
  for (int k=0;k<numvolfrac_;k++)
  {
    // add volume fraction names
    {
      std::ostringstream temp;
      temp << k+1;
      volfracnames_[k] = "VF"+temp.str();
    }
    // add volume fraction pressure names
    {
      std::ostringstream temp;
      temp << k+1;
      volfracpressurenames_[k] = "VFP"+temp.str();
    }
  }

  return;

}

/*----------------------------------------------------------------------*
 | get function                                        kremheller 07/18 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeArt,DRT::Element::DiscretizationType distypeCont>
DRT::UTILS::VariableExprFunction* POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<distypeArt,distypeCont>::Function(
    int functnum) const
{
  try
  {
    DRT::UTILS::VariableExprFunction* funct =
        dynamic_cast<DRT::UTILS::VariableExprFunction*>(&DRT::Problem::Instance()->Funct(functnum));
    return funct;
  }
  catch(std::bad_cast * exp)
  {
    dserror("Cast to VarExp Function failed! For coupling functions only 'VARFUNCTION' functions are allowed!\n"
        "Check your input file!");
    return dynamic_cast<DRT::UTILS::VariableExprFunction*>(&DRT::Problem::Instance()->Funct(functnum));
  }
}

//explicit template instantiations
template class POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<DRT::Element::line2,DRT::Element::quad4>;
template class POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<DRT::Element::line2,DRT::Element::hex8>;
