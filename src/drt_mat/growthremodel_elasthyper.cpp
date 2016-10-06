/*----------------------------------------------------------------------*/
/*!
\file growthremodel_elasthyper.cpp

\brief This file is used to manage the homogenized constraint mixture during growth and remodeling

The input line should read
MAT 1 MAT_GrowthRemodel_ElastHyper NUMMATRF 2 MATIDSRF 11 21 NUMMATEL3D 1 MATIDSEL3D 31 NUMMATEL2D 1 MATIDSEL2D 41
MATIDELPENALTY 51 ELMASSFRAC 0.23 DENS 1050 PRESTRESSTIME 10.0 PRESTRETCHELASTINCIR 1.34 PRESTRETCHELASTINAX 1.25
THICKNESS 1.4099822832e-3 MEANPRESSURE 13332.2668 RADIUS 1e-2 DAMAGE 1 GROWTHTYPE 1 LOCTIMEINT 1

\level 3

\maintainer Fabian Braeu
*/
/*----------------------------------------------------------------------*/
/* headers */
#include "growthremodel_elasthyper.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_matelast/elast_summand.H"
#include "../drt_matelast/elast_remodelfiber.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"
#include "Epetra_SerialDenseSolver.h"


//#define GROWTHPENALTY

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::GrowthRemodel_ElastHyper::GrowthRemodel_ElastHyper(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  nummat_remodelfiber_(matdata->GetInt("NUMMATRF")),
  nummat_elastiniso_(matdata->GetInt("NUMMATEL3D")),
  nummat_elastinmem_(matdata->GetInt("NUMMATEL2D")),
  matids_remodelfiber_(matdata->Get<std::vector<int> >("MATIDSRF")),
  matids_elastiniso_(matdata->Get<std::vector<int> >("MATIDSEL3D")),
  matids_elastinmem_(matdata->Get<std::vector<int> >("MATIDSEL2D")),
  matid_penalty_(matdata->GetInt("MATIDELPENALTY")),
  init_w_el_(matdata->GetDouble("ELMASSFRAC")),
  density_(matdata->GetDouble("DENS")),
  t_prestress_(matdata->GetDouble("PRESTRESSTIME")),
  lamb_prestretch_cir_(matdata->GetDouble("PRESTRETCHELASTINCIR")),
  lamb_prestretch_ax_(matdata->GetDouble("PRESTRETCHELASTINAX")),
  t_ref_(matdata->GetDouble("THICKNESS")),
  p_mean_(matdata->GetDouble("MEANPRESSURE")),
  ri_(matdata->GetDouble("RADIUS")),
  damage_(matdata->GetInt("DAMAGE")),
  growthtype_(matdata->GetInt("GROWTHTYPE")),
  loctimeint_(matdata->GetInt("LOCTIMEINT"))
{
  // check if sizes fit
  if (nummat_remodelfiber_ != (int)matids_remodelfiber_->size())
    dserror("number of remodelfiber materials %d does not fit to size of remodelfiber material vector %d", nummat_remodelfiber_, matids_remodelfiber_->size());

  if (nummat_elastiniso_ != (int)matids_elastiniso_->size())
    dserror("number of elastin materials %d does not fit to size of elastin material vector %d", nummat_elastiniso_, matids_elastiniso_->size());

  if (nummat_elastinmem_ != (int)matids_elastinmem_->size())
    dserror("number of elastin materials %d does not fit to size of elastin material vector %d", nummat_elastinmem_, matids_elastinmem_->size());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthRemodel_ElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::GrowthRemodel_ElastHyper(this));
}


MAT::GrowthRemodel_ElastHyperType MAT::GrowthRemodel_ElastHyperType::instance_;


DRT::ParObject* MAT::GrowthRemodel_ElastHyperType::Create( const std::vector<char> & data )
{
  MAT::GrowthRemodel_ElastHyper* gr_elhy = new MAT::GrowthRemodel_ElastHyper();
  gr_elhy->Unpack(data);

  return gr_elhy;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::GrowthRemodel_ElastHyper::GrowthRemodel_ElastHyper()
  : params_(NULL),
    potsumrf_(0),
    potsumeliso_(0),
    potsumelmem_(0),
    potsumelpenalty_(0),
    t_tot_(0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::GrowthRemodel_ElastHyper::GrowthRemodel_ElastHyper(MAT::PAR::GrowthRemodel_ElastHyper* params)
  : params_(params),
    potsumrf_(0),
    potsumeliso_(0),
    potsumelmem_(0),
    potsumelpenalty_(0),
    t_tot_(0)
{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;

  // RemodelFiber
  for (m=params_->matids_remodelfiber_->begin(); m!=params_->matids_remodelfiber_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::RemodelFiber> sum = Teuchos::rcp_static_cast<MAT::ELASTIC::RemodelFiber>(MAT::ELASTIC::Summand::Factory(matid));
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumrf_.push_back(sum);
  }

  // 3d Elastin matrix
  for (m=params_->matids_elastiniso_->begin(); m!=params_->matids_elastiniso_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumeliso_.push_back(sum);
  }
  for(unsigned p=0;p<potsumeliso_.size();++p)
    if(potsumeliso_[p]->MaterialType() != INPAR::MAT::mes_isoneohooke)
      dserror("3D Elastin Material: So far, prestressing is only implemented for IsoNeoHooke material but it can be easily expanded!");

  // 2d Elastin matrix
  for (m=params_->matids_elastinmem_->begin(); m!=params_->matids_elastinmem_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumelmem_.push_back(sum);
  }
  for(unsigned p=0;p<potsumelmem_.size();++p)
    if(potsumelmem_[p]->MaterialType() != INPAR::MAT::mes_coupneohooke)
      dserror("2D Elastin Material: So far, prestressing is only implemented for CoupNeoHooke material but it can be easily expanded!");

  // VolPenalty
  if(params_->matid_penalty_ != -1)
  {
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(params_->matid_penalty_);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumelpenalty_ = sum;
  }

  // initialize total simulation time
  t_tot_ = 0.0;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);


  // mass fraction of elastin
  AddtoPack(data,cur_rho_el_);
  AddtoPack(data,init_rho_el_);

  AddtoPack(data,v_);
  AddtoPack(data,gp_ax_);
  AddtoPack(data,gp_rad_);
  AddtoPack(data,AcirM_);
  AddtoPack(data,AaxM_);
  AddtoPack(data,AradM_);
  AddtoPack(data,Aradv_);
  AddtoPack(data,AplM_);
  AddtoPack(data,AgM_);
  AddtoPack(data,GM_);
  AddtoPack(data,cylcoords_);
  AddtoPack(data,Q_proj_loc_);
  AddtoPack(data,Q_proj_glob_);
  AddtoPack(data,mue_frac_);
  AddtoPack(data,first_);

  if (params_ != NULL) // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumrf_.size(); ++p)
     potsumrf_[p]->PackSummand(data);

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumeliso_.size(); ++p)
     potsumeliso_[p]->PackSummand(data);

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumelmem_.size(); ++p)
     potsumelmem_[p]->PackSummand(data);

    if(params_->matid_penalty_ != -1)
      potsumelpenalty_->PackSummand(data);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = NULL;
  potsumrf_.clear();
  potsumeliso_.clear();
  potsumelmem_.clear();

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const unsigned int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::GrowthRemodel_ElastHyper*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }
  }


  // mass fractions of elastin and ground matrix
  ExtractfromPack(position,data,cur_rho_el_);
  ExtractfromPack(position,data,init_rho_el_);

  ExtractfromPack(position,data,v_);
  ExtractfromPack(position,data,gp_ax_);
  ExtractfromPack(position,data,gp_rad_);
  ExtractfromPack(position,data,AcirM_);
  ExtractfromPack(position,data,AaxM_);
  ExtractfromPack(position,data,AradM_);
  ExtractfromPack(position,data,Aradv_);
  ExtractfromPack(position,data,AplM_);
  ExtractfromPack(position,data,AgM_);
  ExtractfromPack(position,data,GM_);
  ExtractfromPack(position,data,cylcoords_);
  ExtractfromPack(position,data,Q_proj_loc_);
  ExtractfromPack(position,data,Q_proj_glob_);
  ExtractfromPack(position,data,mue_frac_);
  ExtractfromPack(position,data,first_);

  if (params_ != NULL) // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;

    // RemodelFiber
    for (m=params_->matids_remodelfiber_->begin(); m!=params_->matids_remodelfiber_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::RemodelFiber> sum = Teuchos::rcp_static_cast<MAT::ELASTIC::RemodelFiber>(MAT::ELASTIC::Summand::Factory(matid));
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsumrf_.push_back(sum);
    }

    // 3D Elastin matrix
    for (m=params_->matids_elastiniso_->begin(); m!=params_->matids_elastiniso_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsumeliso_.push_back(sum);
    }

    // 2D Elastin matrix
    for (m=params_->matids_elastinmem_->begin(); m!=params_->matids_elastinmem_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsumelmem_.push_back(sum);
    }

    // VolPenalty
    if(params_->matid_penalty_ != -1)
    {
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(params_->matid_penalty_);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsumelpenalty_ = sum;
    }

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumrf_.size(); ++p)
     potsumrf_[p]->UnpackSummand(data,position);

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumeliso_.size(); ++p)
     potsumeliso_[p]->UnpackSummand(data,position);

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumelmem_.size(); ++p)
     potsumelmem_[p]->UnpackSummand(data,position);

    if(params_->matid_penalty_ != -1)
      potsumelpenalty_->UnpackSummand(data,position);

    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Initialize some variables
  v_.resize(numgp,1.0);
  gp_ax_.resize(numgp,0.0);
  gp_rad_.resize(numgp,0.0);
  cur_rho_el_.resize(numgp);
  init_rho_el_.resize(numgp);
  GM_.resize(numgp,LINALG::Matrix<3,3>(true));
  first_.resize(numgp,1);

  for(int gp=0;gp<numgp;++gp)
  {
    init_rho_el_[gp] = params_->density_ * params_->init_w_el_;
    cur_rho_el_[gp] = init_rho_el_[gp];
  }


  // Setup summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->Setup(numgp,params_->density_,linedef);

  // 3D elastin matrix
  for (unsigned int p=0; p<potsumeliso_.size(); ++p)
    potsumeliso_[p]->Setup(linedef);

  // 2D elastin matrix
  for (unsigned int p=0; p<potsumelmem_.size(); ++p)
    potsumelmem_[p]->Setup(linedef);

  // volpenalty
  if(params_->matid_penalty_ != -1)
    potsumelpenalty_->Setup(linedef);


  // Setup circumferential, radial and axial structural tensor
  SetupAxiCirRadStructuralTensor(linedef);

  //TODO: Here we assume that all fibers are in one plane, so that the growth
  //      direction results from the cross product of two fiber directions.
  //      Later we should implement the evolution equation for the growth tensor but
  //      for now this is enough
  //      structural tensor of growth direction
  AgM_.Update(1.0,AradM_,0.0);

  // structural tensor of the plane in which all fibers are located
  LINALG::Matrix<3,3> id(true);
  for(int i=0;i<3;++i)
    id(i,i) = 1.0;
  AplM_.Update(1.0,AradM_,0.0);
  AplM_.Update(1.0,id,-1.0);


  // setup projection in orthogonal coordinate frame (for the special case of a cylindrical geometry, order: ax,cir,rad)
  LINALG::Matrix<3,3> Q_trafo(true);
  Q_trafo(0,0) = cylcoords_[0](0);
  Q_trafo(1,0) = cylcoords_[0](1);
  Q_trafo(2,0) = cylcoords_[0](2);
  Q_trafo(0,1) = cylcoords_[1](0);
  Q_trafo(1,1) = cylcoords_[1](1);
  Q_trafo(2,1) = cylcoords_[1](2);
  Q_trafo(0,2) = cylcoords_[2](0);
  Q_trafo(1,2) = cylcoords_[2](1);
  Q_trafo(2,2) = cylcoords_[2](2);

  SetupTrafoMatrices(Q_trafo);


  // variable which is multiplied with the material parameter of the elastin sheets
  mue_frac_.resize(numgp,1.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::SetupTrafoMatrices(LINALG::Matrix<3,3> Q_trafo)
{
  // strain like transformation matrices
  Q_proj_loc_(0,0) = Q_proj_glob_(0,0) = Q_trafo(0,0)*Q_trafo(0,0);
  Q_proj_loc_(0,1) = Q_proj_glob_(0,1) = Q_trafo(1,0)*Q_trafo(1,0);
  Q_proj_loc_(0,2) = Q_proj_glob_(0,2) = Q_trafo(2,0)*Q_trafo(2,0);
  Q_proj_loc_(0,3) = 2.0*Q_trafo(0,0)*Q_trafo(1,0);
  Q_proj_glob_(0,3) = 0.5*Q_proj_loc_(0,3);
  Q_proj_loc_(0,4) = 2.0*Q_trafo(1,0)*Q_trafo(2,0);
  Q_proj_glob_(0,4) = 0.5*Q_proj_loc_(0,4);
  Q_proj_loc_(0,5) = 2.0*Q_trafo(0,0)*Q_trafo(2,0);
  Q_proj_glob_(0,5) = 0.5*Q_proj_loc_(0,5);

  Q_proj_loc_(1,0) = Q_proj_glob_(1,0) = Q_trafo(0,1)*Q_trafo(0,1);
  Q_proj_loc_(1,1) = Q_proj_glob_(1,1) = Q_trafo(1,1)*Q_trafo(1,1);
  Q_proj_loc_(1,2) = Q_proj_glob_(1,2) = Q_trafo(2,1)*Q_trafo(2,1);
  Q_proj_loc_(1,3) = 2.0*Q_trafo(0,1)*Q_trafo(1,1);
  Q_proj_glob_(1,3) = 0.5*Q_proj_loc_(1,3);
  Q_proj_loc_(1,4) = 2.0*Q_trafo(1,1)*Q_trafo(2,1);
  Q_proj_glob_(1,4) = 0.5*Q_proj_loc_(1,4);
  Q_proj_loc_(1,5) = 2.0*Q_trafo(0,1)*Q_trafo(2,1);
  Q_proj_glob_(1,5) = 0.5*Q_proj_loc_(1,5);

  Q_proj_loc_(2,0) = Q_proj_glob_(2,0) = Q_trafo(0,2)*Q_trafo(0,2);
  Q_proj_loc_(2,1) = Q_proj_glob_(2,1) = Q_trafo(1,2)*Q_trafo(1,2);
  Q_proj_loc_(2,2) = Q_proj_glob_(2,2) = Q_trafo(2,2)*Q_trafo(2,2);
  Q_proj_loc_(2,3) = 2.0*Q_trafo(0,2)*Q_trafo(1,2);
  Q_proj_glob_(2,3) = 0.5*Q_proj_loc_(2,3);
  Q_proj_loc_(2,4) = 2.0*Q_trafo(1,2)*Q_trafo(2,2);
  Q_proj_glob_(2,4) = 0.5*Q_proj_loc_(2,4);
  Q_proj_loc_(2,5) = 2.0*Q_trafo(0,2)*Q_trafo(2,2);
  Q_proj_glob_(2,5) = 0.5*Q_proj_loc_(2,5);

  Q_proj_loc_(3,0) = Q_trafo(0,0)*Q_trafo(0,1);
  Q_proj_glob_(3,0) = 2.0*Q_proj_loc_(3,0);
  Q_proj_loc_(3,1) = Q_trafo(1,0)*Q_trafo(1,1);
  Q_proj_glob_(3,1) = 2.0*Q_proj_loc_(3,1);
  Q_proj_loc_(3,2) = Q_trafo(2,0)*Q_trafo(2,1);
  Q_proj_glob_(3,2) = 2.0*Q_proj_loc_(3,2);
  Q_proj_loc_(3,3) = Q_proj_glob_(3,3) = Q_trafo(0,0)*Q_trafo(1,1) + Q_trafo(1,0)*Q_trafo(0,1);
  Q_proj_loc_(3,4) = Q_proj_glob_(3,4) = Q_trafo(1,0)*Q_trafo(2,1) + Q_trafo(2,0)*Q_trafo(1,1);
  Q_proj_loc_(3,5) = Q_proj_glob_(3,5) = Q_trafo(0,0)*Q_trafo(2,1) + Q_trafo(2,0)*Q_trafo(0,1);

  Q_proj_loc_(4,0) = Q_trafo(0,1)*Q_trafo(0,2);
  Q_proj_glob_(4,0) = 2.0*Q_proj_loc_(4,0);
  Q_proj_loc_(4,1) = Q_trafo(1,1)*Q_trafo(1,2);
  Q_proj_glob_(4,1) = 2.0*Q_proj_loc_(4,1);
  Q_proj_loc_(4,2) = Q_trafo(2,1)*Q_trafo(2,2);
  Q_proj_glob_(4,2) = 2.0*Q_proj_loc_(4,2);
  Q_proj_loc_(4,3) = Q_proj_glob_(4,3) = Q_trafo(0,1)*Q_trafo(1,2) + Q_trafo(1,1)*Q_trafo(0,2);
  Q_proj_loc_(4,4) = Q_proj_glob_(4,4) = Q_trafo(1,1)*Q_trafo(2,2) + Q_trafo(2,1)*Q_trafo(1,2);
  Q_proj_loc_(4,5) = Q_proj_glob_(4,5) = Q_trafo(0,1)*Q_trafo(2,2) + Q_trafo(2,1)*Q_trafo(0,2);

  Q_proj_loc_(5,0) = Q_trafo(0,0)*Q_trafo(0,2);
  Q_proj_glob_(5,0) = 2.0*Q_proj_loc_(5,0);
  Q_proj_loc_(5,1) = Q_trafo(1,0)*Q_trafo(1,2);
  Q_proj_glob_(5,1) = 2.0*Q_proj_loc_(5,1);
  Q_proj_loc_(5,2) = Q_trafo(2,0)*Q_trafo(2,2);
  Q_proj_glob_(5,2) = 2.0*Q_proj_loc_(5,2);
  Q_proj_loc_(5,3) = Q_proj_glob_(5,3) = Q_trafo(0,0)*Q_trafo(1,2) + Q_trafo(1,0)*Q_trafo(0,2);
  Q_proj_loc_(5,4) = Q_proj_glob_(5,4) = Q_trafo(1,0)*Q_trafo(2,2) + Q_trafo(2,0)*Q_trafo(1,2);
  Q_proj_loc_(5,5) = Q_proj_glob_(5,5) = Q_trafo(0,0)*Q_trafo(2,2) + Q_trafo(2,0)*Q_trafo(0,2);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluatePrestretch(const LINALG::Matrix<3,3>* defgrd,
                                                       const int gp,
                                                       const int eleGID)
{
  // setup prestretch of elastin (radial prestretch is determined in "EvaluatePrestretch(const int gp)")
  GM_[gp].Update(params_->lamb_prestretch_cir_,AcirM_,0.0);
  GM_[gp].Update(params_->lamb_prestretch_ax_,AaxM_,1.0);

  double R = 1.0;
  double dRdlamb_pre = 0.0;
  double lamb_pre = 1./(1.34*1.25);
  int first = 1;

  while(fabs(R) > 1.0e-8)
  {
    if(first == 0)
      lamb_pre = lamb_pre + (-R/dRdlamb_pre);

    R = 72.0*0.23*1050.0*std::pow(1.34*1.34*1.25*1.25*lamb_pre*lamb_pre,-2./3.)*(lamb_pre*lamb_pre - (1./3.)*(1.34*1.34 + 1.25*1.25 + lamb_pre*lamb_pre)) +
        720*0.23*1050*((1.34*1.25*lamb_pre)*(1.34*1.25*lamb_pre)-(1.34*1.25*lamb_pre)) +
        ((1.0-(gp_rad_[gp]-10.0e-3)/params_->t_ref_)*params_->p_mean_);

    dRdlamb_pre = 72.0*0.23*1050.0*(-(4./3.)*std::pow(1.34*1.25*lamb_pre,-7./3.)*1.34*1.25)*(lamb_pre*lamb_pre -(1./3.)*(1.34*1.34 + 1.25*1.25 + lamb_pre*lamb_pre)) +
                  72.0*0.23*1050.0*std::pow(1.34*1.34*1.25*1.25*lamb_pre*lamb_pre,-2./3.)*(2.0*lamb_pre - (1./3.)*(2.0*lamb_pre)) +
                  720*0.23*1050*(2.0*(1.34*1.25*lamb_pre)*1.34*1.25 - 1.34*1.25);

    first = 0;
  }

  // update radial prestretch of elastin
  GM_[gp].Update(lamb_pre,AradM_,1.0);


  // calculate circumferential residual stress
  double sig = 0.0;
  LINALG::Matrix<6,1> stress_aniso(true);
  LINALG::Matrix<6,6> cmat_aniso(true);
  LINALG::Matrix<3,3> id(true);
  for(int i=0;i<3;++i)
    id(i,i) = 1.0;
  LINALG::Matrix<6,1> Acir_strain(true);
  Acir_strain(0) = AcirM_(0,0);
  Acir_strain(1) = AcirM_(1,1);
  Acir_strain(2) = AcirM_(2,2);
  Acir_strain(3) = AcirM_(0,1)+AcirM_(1,0);
  Acir_strain(4) = AcirM_(1,2)+AcirM_(2,1);
  Acir_strain(5) = AcirM_(0,2)+AcirM_(2,0);

  // total circumferential Cauchy stress ("membrane stress")
  double sig_tot = (params_->p_mean_*params_->ri_)/params_->t_ref_;

  // evaluate anisotropic remodel fibers
  LINALG::Matrix<3,3> CM(true);
  CM.MultiplyTN(1.0,*defgrd,*defgrd,0.0);
  for(unsigned p=0;p<potsumrf_.size();++p)
  {
    potsumrf_[p]->EvaluateAnisotropicStressCmat(CM,id,cmat_aniso,stress_aniso,gp,eleGID);
    sig += stress_aniso.Dot(Acir_strain);
  }

  // build stress response and elasticity tensor of 3D material
  LINALG::Matrix<NUM_STRESS_3D,1> stressiso(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);

  EvaluateStressCmat(defgrd,GM_[gp],gp,stressiso,cmatiso,eleGID);

  sig += stressiso.Dot(Acir_strain);

  // residual stress
  double sig_res = sig_tot - sig;

  if((potsumelmem_.size() != 1) || (potsumelmem_[0]->MaterialType() != INPAR::MAT::mes_coupneohooke))
    dserror("So far, only ONE CoupNeoHooke material can be chosen for the membrane \"material\"");
  else
  {
    LINALG::Matrix<3,1> stress_loc(true);
    LINALG::Matrix<3,3> cmat_loc(true);

    // total inelastic deformation gradient for elastin in local system
    LINALG::Matrix<3,3> iFin_loc(true);
    LINALG::Matrix<6,1> iFinv_glob(true);
    LINALG::Matrix<6,1> iFinv_loc(true);
    MatrixtoStressLikeVoigtNotation(GM_[gp],iFinv_glob);
    iFinv_loc.MultiplyNN(1.0,Q_proj_loc_,iFinv_glob,0.0);
    StressLikeVoigtNotationtoMatrix(iFinv_loc,iFin_loc);

    // local right Cauchy-Green tensor
    LINALG::Matrix<3,3> C_loc(true);
    LINALG::Matrix<3,3> C_glob(true);
    LINALG::Matrix<6,1> Cv_glob(true);
    LINALG::Matrix<6,1> Cv_loc(true);
    C_glob.MultiplyTN(1.0,*defgrd,*defgrd,0.0);
    MatrixtoStressLikeVoigtNotation(C_glob,Cv_glob);
    Cv_loc.MultiplyNN(1.0,Q_proj_loc_,Cv_glob,0.0);
    StressLikeVoigtNotationtoMatrix(Cv_loc,C_loc);

    EvaluateStressCmatMembrane(C_loc,iFin_loc,stress_loc,cmat_loc,gp,eleGID);

    // convert stress and cmat to global coordinate system
    LINALG::Matrix<6,1> stress_glob(true);
    LINALG::Matrix<6,6> cmat_glob(true);
    ConvertStressCmatLoctoGlob(stress_loc,cmat_loc,stress_glob,cmat_glob);

    mue_frac_[gp] = sig_res/stress_glob.Dot(Acir_strain);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::SetupAxiCirRadStructuralTensor(DRT::INPUT::LineDefinition* linedef)
{
  // CIR-AXI-RAD nomenclature
  if (linedef->HaveNamed("RAD") and
      linedef->HaveNamed("AXI") and
      linedef->HaveNamed("CIR"))
  {
    // Read in of data
    // read local (cylindrical) cosy-directions at current element
    LINALG::Matrix<3,1> dir;
    cylcoords_.resize(3,LINALG::Matrix<1,3> (true));

    // Axial direction
    ReadDir(linedef,"AXI",dir);
    cylcoords_[0].UpdateT(1.0,dir,0.0);
    AaxM_.MultiplyNT(1.0,dir,dir,0.0);

    // Circumferential direction
    ReadDir(linedef,"CIR",dir);
    cylcoords_[1].UpdateT(1.0,dir,0.0);
    AcirM_.MultiplyNT(1.0,dir,dir,0.0);

    // Radial direction
    ReadDir(linedef,"RAD",dir);
    cylcoords_[2].UpdateT(1.0,dir,0.0);
    AradM_.MultiplyNT(1.0,dir,dir,0.0);

    // radial structural tensor in "stress-like" Voigt notation
    for(int i=0;i<3;++i)
      Aradv_(i) = AradM_(i,i);
    Aradv_(3) = 0.5*(AradM_(0,1)+AradM_(1,0));
    Aradv_(4) = 0.5*(AradM_(1,2)+AradM_(2,1));
    Aradv_(5) = 0.5*(AradM_(0,2)+AradM_(2,0));
  }
  // No AXI CIR RAD direction defined in .dat file
  else
    dserror("Homogenized Constrained Mixture Model can so far only be used by defining AXI-, CIR- and RAD-direction in the .dat file!");

  return;
}


/*----------------------------------------------------------------------*
 * Function which reads in the AXI CIR RAD directions
 *----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::ReadDir(
    DRT::INPUT::LineDefinition* linedef,
    std::string specifier,
    LINALG::Matrix<3,1> &dir)
{
  std::vector<double> fiber;
  linedef->ExtractDoubleVector(specifier,fiber);
  double fnorm=0.;
  //normalization
  for (int i = 0; i < 3; ++i)
  {
    fnorm += fiber[i]*fiber[i];
  }
  fnorm = sqrt(fnorm);

  // fill final normalized vector
  for (int i = 0; i < 3; ++i)
    dir(i) = fiber[i]/fnorm;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Update(const LINALG::Matrix<3,3> defgrd,
                                           const int gp,
                                           Teuchos::ParameterList& params,
                                           const int eleGID)
{
  // Update individual volume of elastin
  if(t_tot_ > params_->t_prestress_)
  {
    if(params_->damage_ == 1)
      EvaluateElastinDamage();
  }

  // loop map of associated potential summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->Update();

  // 3D elastin matrix
  for (unsigned int p=0; p<potsumeliso_.size(); ++p)
    potsumeliso_[p]->Update();

  // 2D elastin matrix
  for (unsigned int p=0; p<potsumelmem_.size(); ++p)
    potsumelmem_[p]->Update();

  // volpenalty
  if(params_->matid_penalty_ != -1)
    potsumelpenalty_->Update();


  if(t_tot_ > params_->t_prestress_)
  {
    // build inelastic growth deformation gradient
    LINALG::Matrix<3,3> FgM(true);
    LINALG::Matrix<3,3> iFgM(true);

    EvaluateGrowthDefGrad(FgM,iFgM,gp);

    // time step size
    double dt = params.get<double>("delta time");

    switch(params_->loctimeint_)
    {
    case 0:
      for (unsigned int p=0; p<potsumrf_.size(); ++p)
        potsumrf_[p]->EvaluateGrowthAndRemodelingExpl(defgrd,dt,iFgM,gp,eleGID);
      break;
    case 1: // do nothing
      break;
    default:
      dserror("LOCTIMEINT has to be either 1 (Backward Euler) or 0 (Forward Euler)");
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::UpdateMembrane(const LINALG::Matrix<3,3> defgrd_loc,
                                                   const LINALG::Matrix<3,3> Q_trafo,
                                                   const int gp,
                                                   Teuchos::ParameterList& params,
                                                   const int eleGID)
{
  // Update individual volume of elastin
  if(t_tot_ > params_->t_prestress_)
  {
    if(params_->damage_ == 1)
      EvaluateElastinDamage();
  }

  // loop map of associated potential summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->Update();

  // 2D elastin matrix
  for (unsigned int p=0; p<potsumelmem_.size(); ++p)
    potsumelmem_[p]->Update();


  // just update G&R deformation gradient here, when no membrane formulation is used
  if(t_tot_ > params_->t_prestress_)
  {
    // build inelastic growth deformation gradient
    LINALG::Matrix<3,3> FgM(true);
    LINALG::Matrix<3,3> iFgM(true);

    EvaluateGrowthDefGrad(FgM,iFgM,gp);

    // time step size
    double dt = params.get<double>("delta time");

    switch(params_->loctimeint_)
    {
    case 0:
      for (unsigned int p=0; p<potsumrf_.size(); ++p)
        potsumrf_[p]->EvaluateGrowthAndRemodelingExplMembrane(defgrd_loc,Q_trafo,v_[gp],dt,iFgM,gp,eleGID);
      break;
    case 1: // do nothing
      break;
    default:
      dserror("LOCTIMEINT has to be either 1 (Backward Euler) or 0 (Forward Euler)");
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Evaluate(const LINALG::Matrix<3,3>* defgrd,
                                             const LINALG::Matrix<6,1>* glstrain,
                                             Teuchos::ParameterList& params,
                                             LINALG::Matrix<6,1>* stress,
                                             LINALG::Matrix<6,6>* cmat,
                                             const int eleGID)
{
  // save current simulation time (used for the evaluation of elastin degradation)
  t_tot_ = params.get<double>("total time");

  // current gp
  int gp = params.get<int>("gp");

  // time step size
  double dt = params.get<double>("delta time");

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress->Clear();
  cmat->Clear();

  // right Cauchy Green tensor
  LINALG::Matrix<3,3> CM(true);
  CM.MultiplyTN(1.0,*defgrd,*defgrd,0.0);

  // evaluate remodelfiber (remodeling of collagen)
  // build stress response and elasticity tensor
  LINALG::Matrix<NUM_STRESS_3D,1> stressremodel(true) ;
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatremodel(true) ;

  // sum over all fiber families of the derivatives of the individual mass density w.r.t. right Cauchy Green tensor
  LINALG::Matrix<1,6> sum_drhodC(true);

  // build growth deformation gradient
  LINALG::Matrix<3,3> FgM(true);
  LINALG::Matrix<3,3> iFgM(true);
  LINALG::Matrix<3,3> iFin(true);
  EvaluateGrowthDefGrad(FgM,iFgM,gp);

  // prestressing
  if(t_tot_ <= params_->t_prestress_)
  {
    if(first_[gp] == 1)
    {
      // gp coordinates in reference configuration
      LINALG::Matrix<1,3> gprefecoord(true);
      gprefecoord = params.get<LINALG::Matrix<1,3> >("gprefecoord");
      gp_ax_[gp] = cylcoords_[0].Dot(gprefecoord);
      gp_rad_[gp] = cylcoords_[2].Dot(gprefecoord);

      //
      // TODO: BE CAREFULL! So far, this prestretching procedure is only valid for certain materials and a cylindrical geometry.
      //       The principle of the prestretching routine can easily be adapted to other materials or general geometries!!!
      //
      EvaluatePrestretch(defgrd,gp,eleGID);
      first_[gp] = 0;
    }

    // evaluate anisotropic remodel fibers
    for(unsigned p=0;p<potsumrf_.size();++p)
    {
      potsumrf_[p]->EvaluateAnisotropicStressCmat(CM,iFgM,cmatremodel,stressremodel,gp,eleGID);
      stress->Update(1.0,stressremodel,1.0);
      cmat->Update(1.0,cmatremodel,1.0);
    }


    // build stress response and elasticity tensor of 3D material
    LINALG::Matrix<NUM_STRESS_3D,1> stressiso(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);

    EvaluateStressCmat(defgrd,GM_[gp],gp,stressiso,cmatiso,eleGID);

    stress->Update(1.0,stressiso,1.0);
    cmat->Update(1.0,cmatiso,1.0);


    // build stress response and elasticity tensor of quasi 2D material
    LINALG::Matrix<3,1> stress_loc(true) ;
    LINALG::Matrix<3,3> cmat_loc(true) ;

    // total inelastic deformation gradient for elastin in local system
    LINALG::Matrix<3,3> iFin_loc(true);
    LINALG::Matrix<6,1> iFinv_glob(true);
    LINALG::Matrix<6,1> iFinv_loc(true);
    iFin.MultiplyNN(1.0,iFgM,GM_[gp],0.0);
    MatrixtoStressLikeVoigtNotation(iFin,iFinv_glob);
    iFinv_loc.MultiplyNN(1.0,Q_proj_loc_,iFinv_glob,0.0);
    StressLikeVoigtNotationtoMatrix(iFinv_loc,iFin_loc);

    // local right Cauchy-Green tensor
    LINALG::Matrix<3,3> C_loc(true);
    LINALG::Matrix<3,3> C_glob(true);
    LINALG::Matrix<6,1> Cv_glob(true);
    LINALG::Matrix<6,1> Cv_loc(true);
    C_glob.MultiplyTN(1.0,*defgrd,*defgrd,0.0);
    MatrixtoStressLikeVoigtNotation(C_glob,Cv_glob);
    Cv_loc.MultiplyNN(1.0,Q_proj_loc_,Cv_glob,0.0);
    StressLikeVoigtNotationtoMatrix(Cv_loc,C_loc);

    EvaluateStressCmatMembrane(C_loc,iFin_loc,stress_loc,cmat_loc,gp,eleGID);

    // convert stress and cmat to global coordinate system
    LINALG::Matrix<6,1> stress_glob(true);
    LINALG::Matrix<6,6> cmat_glob(true);
    ConvertStressCmatLoctoGlob(stress_loc,cmat_loc,stress_glob,cmat_glob);

    stress->Update(1.0,stress_glob,1.0);
    cmat->Update(1.0,cmat_glob,1.0);

    return;
  }
  // Growth and Remodeling
  else
  {
    switch(params_->loctimeint_)
    {
    case 0:
    {
      // evaluate anisotropic remodel fibers
      for(unsigned p=0;p<potsumrf_.size();++p)
      {
        potsumrf_[p]->EvaluateAnisotropicStressCmat(CM,iFgM,cmatremodel,stressremodel,gp,eleGID);
        stress->Update(1.0,stressremodel,1.0);
        cmat->Update(1.0,cmatremodel,1.0);
      }

      // current inelastic deformation gradient for elastin
      iFin.MultiplyNN(1.0,iFgM,GM_[gp],0.0);

      // build stress response and elasticity tensor
      LINALG::Matrix<NUM_STRESS_3D,1> stressiso(true) ;
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true) ;

      EvaluateStressCmat(defgrd,iFin,gp,stressiso,cmatiso,eleGID);

      stress->Update(1.0,stressiso,1.0);
      cmat->Update(1.0,cmatiso,1.0);

      // build stress response and elasticity tensor of quasi 2D material
      LINALG::Matrix<3,1> stress_loc(true) ;
      LINALG::Matrix<3,3> cmat_loc(true) ;

      // total inelastic deformation gradient for elastin in local system
      LINALG::Matrix<3,3> iFin_loc(true);
      LINALG::Matrix<6,1> iFinv_glob(true);
      LINALG::Matrix<6,1> iFinv_loc(true);
      iFin.MultiplyNN(1.0,iFgM,GM_[gp],0.0);
      MatrixtoStressLikeVoigtNotation(iFin,iFinv_glob);
      iFinv_loc.MultiplyNN(1.0,Q_proj_loc_,iFinv_glob,0.0);
      StressLikeVoigtNotationtoMatrix(iFinv_loc,iFin_loc);

      // local right Cauchy-Green tensor
      LINALG::Matrix<3,3> C_loc(true);
      LINALG::Matrix<3,3> C_glob(true);
      LINALG::Matrix<6,1> Cv_glob(true);
      LINALG::Matrix<6,1> Cv_loc(true);
      C_glob.MultiplyTN(1.0,*defgrd,*defgrd,0.0);
      MatrixtoStressLikeVoigtNotation(C_glob,Cv_glob);
      Cv_loc.MultiplyNN(1.0,Q_proj_loc_,Cv_glob,0.0);
      StressLikeVoigtNotationtoMatrix(Cv_loc,C_loc);

      EvaluateStressCmatMembrane(C_loc,iFin_loc,stress_loc,cmat_loc,gp,eleGID);

      // convert stress and cmat to global coordinate system
      LINALG::Matrix<6,1> stress_glob(true);
      LINALG::Matrix<6,6> cmat_glob(true);
      ConvertStressCmatLoctoGlob(stress_loc,cmat_loc,stress_glob,cmat_glob);

      stress->Update(1.0,stress_glob,1.0);
      cmat->Update(1.0,cmat_glob,1.0);

      break;
    }
    case 1:
    {
      // total number of fibers
      unsigned nr_f_tot = 0;

      for(unsigned p=0;p<potsumrf_.size();++p)
        nr_f_tot += potsumrf_[p]->GetNumFibers();


      // allocate some variables
      std::vector<std::vector<double> > dWdrho(nr_f_tot,std::vector<double>(nr_f_tot,0.0));
      std::vector<std::vector<double> > dWdlamb(nr_f_tot,std::vector<double>(nr_f_tot,0.0));
      std::vector<double> W(nr_f_tot,0.0);
      std::vector<std::vector<double> > dEdrho(nr_f_tot,std::vector<double>(nr_f_tot,0.0));
      std::vector<std::vector<double> > dEdlamb(nr_f_tot,std::vector<double>(nr_f_tot,0.0));
      std::vector<double> E(nr_f_tot,0.0);


      // epetra serial dense solver
      Epetra_SerialDenseSolver solver;

      // "tangential matrix"
      LINALG::SerialDenseMatrix K_T(2*nr_f_tot,2*nr_f_tot,true);

      // residual vector of assembled system of equation
      LINALG::SerialDenseMatrix R(2*nr_f_tot,1);
      for(unsigned i=0;i<2*nr_f_tot;++i)
        R(i,0) = 1.0;

      // solution vector of assembled system of equation
      LINALG::SerialDenseMatrix dsol(2*nr_f_tot,1);

      int nr_grf_proc = 0;
      int nr_loop = 0;
      int l=0;
      // identity tensor
      LINALG::Matrix<3,3> id(true);
      for(int i=0;i<3;++i)
        id(i,i) = 1.0;
      while(fabs(R.NormInf()) > 1.0e-6)
      {
        if(nr_loop != 0)
        {
          // Solve linearized system of equations
          solver.SetMatrix(K_T);
          solver.SetVectors(dsol,R);
          solver.SolveToRefinedSolution(true);
          solver.ApplyRefinement();
          solver.Solve();


          l=0;
          for(unsigned p=0;p<potsumrf_.size();++p)
            for(unsigned k=0;k<potsumrf_[p]->GetNumFibers();++k)
            {
              // update inelastic fiber stretch caused by remodeling
              potsumrf_[p]->UpdateCurRemodelStretch(dsol(l,0),k,gp);

              // update collagen density
              potsumrf_[p]->UpdateCurMassDensity(dsol(nr_f_tot+l,0),k,gp);

              l++;
            }

          EvaluateGrowthDefGrad(FgM,iFgM,gp);
        }

        // global number of fibers which were already processed
        nr_grf_proc = 0;
        for(unsigned p=0;p<potsumrf_.size();++p)
        {
          potsumrf_[p]->EvaluateDerivativesInternalNewton(defgrd,id,nr_grf_proc,params_->density_,gp,dt,v_[gp],FgM,iFgM,AgM_,AcirM_,
              AradM_,AaxM_,dWdrho,dWdlamb,W,dEdrho,dEdlamb,E,eleGID,params_->growthtype_);
          nr_grf_proc += potsumrf_[p]->GetNumFibers();
        }

        // Assembly
        // R
        for(unsigned i=0;i<nr_f_tot;++i)
        {
          R(i,0) = -E[i];
          R(nr_f_tot+i,0) = -W[i];
        }

        // K_T
        for(unsigned i=0;i<nr_f_tot;++i)
          for(unsigned j=0;j<nr_f_tot;++j)
          {
            K_T(i,j) = dEdlamb[i][j];
            K_T(nr_f_tot+i,j) = dWdlamb[i][j];
          }

        for(unsigned i=0;i<nr_f_tot;++i)
          for(unsigned j=0;j<nr_f_tot;++j)
          {
            K_T(i,nr_f_tot+j) = dEdrho[i][j];
            K_T(nr_f_tot+i,nr_f_tot+j) = dWdrho[i][j];
          }

        nr_loop++;
      }


      std::vector<LINALG::Matrix<1,6> > dWdC(nr_f_tot,LINALG::Matrix<1,6>(true));
      std::vector<LINALG::Matrix<1,6> > dEdC(nr_f_tot,LINALG::Matrix<1,6>(true));

      nr_grf_proc = 0;
      for(unsigned p=0;p<potsumrf_.size();++p)
      {
        potsumrf_[p]->EvaluateDerivativesCauchyGreen(defgrd,id,nr_grf_proc,gp,dt,FgM,iFgM,dWdC,dEdC,eleGID);
        nr_grf_proc += potsumrf_[p]->GetNumFibers();
      }

      // Assembly
      // Rcmat
      LINALG::SerialDenseMatrix Rcmat(2*nr_f_tot,6);
      for(unsigned i=0;i<nr_f_tot;++i)
        for(unsigned j=0;j<6;++j)
        {
          Rcmat(i,j) = -dEdC[i](0,j);
          Rcmat(nr_f_tot+i,j) = -dWdC[i](0,j);
        }

      // Solve
      LINALG::SerialDenseMatrix dsolcmat(2*nr_f_tot,6);
      solver.SetMatrix(K_T);
      solver.SetVectors(dsolcmat,Rcmat);
      solver.SolveToRefinedSolution(true);
      solver.ApplyRefinement();
      solver.EstimateSolutionErrors(true);
      solver.Solve();



      // split solution vector in dlambda_r/dC and drho/dC
      std::vector<LINALG::Matrix<1,6> > drhodC(nr_f_tot,LINALG::Matrix<1,6>(true));
      std::vector<LINALG::Matrix<1,6> > dlambdC(nr_f_tot,LINALG::Matrix<1,6>(true));
      sum_drhodC.Clear();
      for(unsigned i=0;i<nr_f_tot;++i)
      {
        for(unsigned j=0;j<6;++j)
        {
          dlambdC[i](0,j) = dsolcmat(i,j);
          drhodC[i](0,j) = dsolcmat(nr_f_tot+i,j);
        }
        sum_drhodC.Update(1.0,drhodC[i],1.0);
      }


      // update stress and elasticity tensor
      nr_grf_proc = 0;
      for(unsigned p=0;p<potsumrf_.size();++p)
      {
        potsumrf_[p]->AddStressCmatGrowthRemodel(defgrd,id,nr_f_tot,nr_grf_proc,gp,v_[gp],params_->density_,iFgM,AgM_,AcirM_,AradM_,AaxM_,
            drhodC,dlambdC,stressremodel,cmatremodel,eleGID,params_->growthtype_);
        stress->Update(1.0,stressremodel,1.0);
        cmat->Update(1.0,cmatremodel,1.0);
        nr_grf_proc += potsumrf_[p]->GetNumFibers();
      }


      // Evaluate elastin
      // some variables
      LINALG::Matrix<6,1> iCg(true);
      LINALG::Matrix<6,1> iCgCiCg(true);
      LINALG::Matrix<6,1> iC(true);
      LINALG::Matrix<3,1> prinv(true);
      LINALG::Matrix<3,3> iCgC(true);
      LINALG::Matrix<3,3> iFgCe(true);
      LINALG::Matrix<9,1> CiFg(true);
      LINALG::Matrix<9,1> CiFgCe(true);
      LINALG::Matrix<9,1> CiFgiCe(true);

      // first and second derivatives w.r.t. invariants weighted with the corresponding volume fraction
      LINALG::Matrix<3,1> dPIw(true);
      LINALG::Matrix<6,1> ddPIIw(true);

      iFin.MultiplyNN(1.0,iFgM,GM_[gp],0.0);
      EvaluateKinQuantElast(defgrd,iFin,gp,iCg,iCgCiCg,iC,iCgC,iFgCe,CiFg,CiFgCe,CiFgiCe,prinv);
      EvaluateInvariantDerivatives(prinv,gp,dPIw,ddPIIw,eleGID);


      // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
      LINALG::Matrix<3,1> gamma(true);
      // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
      LINALG::Matrix<8,1> delta(true);

      // compose coefficients
      CalculateGammaDelta(gamma,delta,prinv,dPIw,ddPIIw);


      // build stress response and elasticity tensor of 3D material
      LINALG::Matrix<NUM_STRESS_3D,1> stressiso(true) ;
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true) ;

      EvaluateIsotropicPrincElast(stressiso,cmatiso,iCg,iCgCiCg,iC,gamma,delta);

      stress->Update(1.0,stressiso,1.0);
      cmat->Update(1.0,cmatiso,1.0);


      // build additional terms for elasticity tensor
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatadd(true);
      EvaluateAdditionalTermsCmat(cmatadd,gamma,delta,v_[gp],gp,AgM_,sum_drhodC,iFin,iCgC,iFgCe,iCg,CiFg,CiFgCe,iCgCiCg,CiFgiCe,iC);

      cmat->Update(1.0,cmatadd,1.0);


      // build stress response and elasticity tensor of quasi 2D material
      LINALG::Matrix<3,1> stress_loc(true) ;
      LINALG::Matrix<3,3> cmat_loc(true) ;

      // total inelastic deformation gradient for elastin in local system
      LINALG::Matrix<3,3> iFin_loc(true);
      LINALG::Matrix<6,1> iFinv_glob(true);
      LINALG::Matrix<6,1> iFinv_loc(true);
      iFin.MultiplyNN(1.0,iFgM,GM_[gp],0.0);
      MatrixtoStressLikeVoigtNotation(iFin,iFinv_glob);
      iFinv_loc.MultiplyNN(1.0,Q_proj_loc_,iFinv_glob,0.0);
      StressLikeVoigtNotationtoMatrix(iFinv_loc,iFin_loc);

      // local right Cauchy-Green tensor
      LINALG::Matrix<3,3> C_loc(true);
      LINALG::Matrix<3,3> C_glob(true);
      LINALG::Matrix<6,1> Cv_glob(true);
      LINALG::Matrix<6,1> Cv_loc(true);
      C_glob.MultiplyTN(1.0,*defgrd,*defgrd,0.0);
      MatrixtoStressLikeVoigtNotation(C_glob,Cv_glob);
      Cv_loc.MultiplyNN(1.0,Q_proj_loc_,Cv_glob,0.0);
      StressLikeVoigtNotationtoMatrix(Cv_loc,C_loc);

      EvaluateStressCmatMembrane(C_loc,iFin_loc,stress_loc,cmat_loc,gp,eleGID);

      // convert stress and cmat to global coordinate system
      LINALG::Matrix<6,1> stress_glob(true);
      LINALG::Matrix<6,6> cmat_glob(true);
      ConvertStressCmatLoctoGlob(stress_loc,cmat_loc,stress_glob,cmat_glob);

      stress->Update(1.0,stress_glob,1.0);
      cmat->Update(1.0,cmat_glob,1.0);

      break;
    }
    default:
      dserror("LOCTIMEINT has to be 1 (Backward Euler Method) or 0 (Forward Euler Method)");
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateStressCmat(const LINALG::Matrix<3,3>* defgrd,
                                                       const LINALG::Matrix<3,3> iFinel,
                                                       const int gp,
                                                       LINALG::Matrix<6,1>& stressiso,
                                                       LINALG::Matrix<6,6>& cmatiso,
                                                       const int eleGID)
{
  // Evaluate elastin matrix
  // some variables
  LINALG::Matrix<6,1> iCinel(true);
  LINALG::Matrix<6,1> iCinelCiCinel(true);
  LINALG::Matrix<6,1> iC(true);
  LINALG::Matrix<3,1> prinv(true);
  LINALG::Matrix<3,3> iCinelC(true);
  LINALG::Matrix<3,3> iFinelCe(true);
  LINALG::Matrix<9,1> CiFinel(true);
  LINALG::Matrix<9,1> CiFinelCe(true);
  LINALG::Matrix<9,1> CiFineliCe(true);

  // first and second derivatives w.r.t. invariants weighted with the corresponding volume fraction
  LINALG::Matrix<3,1> dPIe(true);
  LINALG::Matrix<6,1> ddPIIe(true);
  LINALG::Matrix<3,1> dPmodI(true);
  LINALG::Matrix<6,1> ddPmodII(true);


  EvaluateKinQuantElast(defgrd,iFinel,gp,iCinel,iCinelCiCinel,iC,iCinelC,iFinelCe,CiFinel,CiFinelCe,CiFineliCe,prinv);

  LINALG::Matrix<3,1> modinv(true);
  InvariantsModified(modinv,prinv);
  LINALG::Matrix<3,1> dPgrowthI(true);
  LINALG::Matrix<6,1> ddPgrowthII(true);
  // loop map of associated potential summands
  // elastin matrix
  // derivatives of strain energy function w.r.t. principal invariants
  for(unsigned p=0;p<potsumeliso_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumeliso_[p]->AddDerivativesPrincipal(dPgrowthI,ddPgrowthII,prinv,eleGID);
    dPIe.Update(cur_rho_el_[gp],dPgrowthI,1.0);
    ddPIIe.Update(cur_rho_el_[gp],ddPgrowthII,1.0);
  }

  // derivatives of strain energy function w.r.t. modified invariants
  for(unsigned p=0;p<potsumeliso_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumeliso_[p]->AddDerivativesModified(dPgrowthI,ddPgrowthII,modinv,0);
    dPmodI.Update(cur_rho_el_[gp],dPgrowthI,1.0);
    ddPmodII.Update(cur_rho_el_[gp],ddPgrowthII,1.0);
  }

  // evaluate penalty term
  dPgrowthI.Clear();
  ddPgrowthII.Clear();
  if(params_->matid_penalty_ != -1)
    potsumelpenalty_->AddDerivativesModified(dPgrowthI,ddPgrowthII,modinv,0);
  dPmodI.Update(cur_rho_el_[gp],dPgrowthI,1.0);
  ddPmodII.Update(cur_rho_el_[gp],ddPgrowthII,1.0);

  // convert decoupled derivatives to principal derivatives
  ConvertModToPrinc(prinv,dPmodI,ddPmodII,dPIe,ddPIIe);

  // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  LINALG::Matrix<3,1> gamma(true);
  // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  LINALG::Matrix<8,1> delta(true);

  // compose coefficients
  CalculateGammaDelta(gamma,delta,prinv,dPIe,ddPIIe);

  EvaluateIsotropicPrincElast(stressiso,cmatiso,iCinel,iCinelCiCinel,iC,gamma,delta);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateKinQuantElast(
    const LINALG::Matrix<3,3>* defgrd,
    LINALG::Matrix<3,3> iFinel,
    const int gp,
    LINALG::Matrix<6,1>& iCinel,
    LINALG::Matrix<6,1>& iCinelCiCinel,
    LINALG::Matrix<6,1>& iC,
    LINALG::Matrix<3,3>& iCinelC,
    LINALG::Matrix<3,3>& iFinelCe,
    LINALG::Matrix<9,1>& CiFinel,
    LINALG::Matrix<9,1>& CiFinelCe,
    LINALG::Matrix<9,1>& CiFineliCe,
    LINALG::Matrix<3,1>& prinv)
{
  // inverse inelastic right Cauchy-Green
  LINALG::Matrix<3,3> iCinelM(true);
  iCinelM.MultiplyNT(iFinel,iFinel);
  // stress-like Voigt notation
  for (int i=0; i<3; i++) iCinel(i) = iCinelM(i,i);
  iCinel(3) = (iCinelM(0,1)+iCinelM(1,0))/2.;
  iCinel(4) = (iCinelM(2,1)+iCinelM(1,2))/2.;
  iCinel(5) = (iCinelM(0,2)+iCinelM(2,0))/2.;

  // inverse RCG
  LINALG::Matrix<3,3> iRCG(true);
  LINALG::Matrix<3,3> RCG(true);
  RCG.MultiplyTN(*defgrd,*defgrd);
  iRCG.Invert(RCG);
  // stress-like Voigt notation
  for (int i=0; i<3; i++) iC(i) = iRCG(i,i);
  iC(3) = (iRCG(0,1)+iRCG(1,0))/2.;
  iC(4) = (iRCG(2,1)+iRCG(1,2))/2.;
  iC(5) = (iRCG(0,2)+iRCG(2,0))/2.;

  // C_g^-1 * C * C_g^-1
  LINALG::Matrix<3,3> tmp(true);
  LINALG::Matrix<3,3> iCinelCiCinelM;
  tmp.Multiply(iCinelM,RCG);
  iCinelCiCinelM.Multiply(tmp,iCinelM);
  // stress-like Voigt notation
  for (int i=0; i<3; i++) iCinelCiCinel(i) = iCinelCiCinelM(i,i);
  iCinelCiCinel(3) = (iCinelCiCinelM(0,1)+iCinelCiCinelM(1,0))/2.;
  iCinelCiCinel(4) = (iCinelCiCinelM(2,1)+iCinelCiCinelM(1,2))/2.;
  iCinelCiCinel(5) = (iCinelCiCinelM(0,2)+iCinelCiCinelM(2,0))/2.;

  tmp.Multiply(*defgrd,iFinel);
  LINALG::Matrix<3,3> CeM(true);
  CeM.MultiplyTN(tmp,tmp);
  // elastic right Cauchy-Green in strain-like Voigt notation.
  LINALG::Matrix<6,1> elasticRCGv(true);
  for (int i=0; i<3; i++)
    elasticRCGv(i)=CeM(i,i);
  elasticRCGv(3) = (CeM(0,1)+CeM(1,0));
  elasticRCGv(4) = (CeM(2,1)+CeM(1,2));
  elasticRCGv(5) = (CeM(0,2)+CeM(2,0));

  // principal invariants of elastic Cauchy-Green strain
  InvariantsPrincipal(prinv,elasticRCGv);

  // C_g^-1 * C
  iCinelC.MultiplyNN(1.0,iCinelM,RCG,0.0);

  // F_g^-1 * C_e
  iFinelCe.MultiplyNN(1.0,iFinel,CeM,0.0);

  // C * F_g^-1
  LINALG::Matrix<3,3> CiFgM(true);
  CiFgM.MultiplyNN(1.0,RCG,iFinel,0.0);
  Matrix3x3to9x1(CiFgM,CiFinel);

  // C * F_g^-1 * C_e
  LINALG::Matrix<3,3> CiFgCeM(true);
  tmp.MultiplyNN(1.0,RCG,iFinel,0.0);
  CiFgCeM.MultiplyNN(1.0,tmp,CeM,0.0);
  Matrix3x3to9x1(CiFgCeM,CiFinelCe);

  // C * F_g^-1 * C_e^-1
  LINALG::Matrix<3,3> CiFgiCeM(true);
  LINALG::Matrix<3,3> iCeM(true);
  iCeM.Invert(CeM);
  tmp.MultiplyNN(1.0,RCG,iFinel,0.0);
  CiFgiCeM.MultiplyNN(1.0,tmp,iCeM,0.0);
  Matrix3x3to9x1(CiFgiCeM,CiFineliCe);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::InvariantsPrincipal(
    LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<6,1>& C_strain)
{
  // 1st invariant, trace
  prinv(0) = C_strain(0) + C_strain(1) + C_strain(2);
  // 2nd invariant
  prinv(1) = 0.5*( prinv(0)*prinv(0)
                   - C_strain(0)*C_strain(0) - C_strain(1)*C_strain(1) - C_strain(2)*C_strain(2)
                   - .5*C_strain(3)*C_strain(3) - .5*C_strain(4)*C_strain(4) - .5*C_strain(5)*C_strain(5) );
  // 3rd invariant, determinant
  prinv(2) = C_strain(0)*C_strain(1)*C_strain(2)
    + 0.25 * C_strain(3)*C_strain(4)*C_strain(5)
    - 0.25 * C_strain(1)*C_strain(5)*C_strain(5)
    - 0.25 * C_strain(2)*C_strain(3)*C_strain(3)
    - 0.25 * C_strain(0)*C_strain(4)*C_strain(4);
}


/*----------------------------------------------------------------------/
 * Reads derivatives with respect to invariants and modified invariants
 * from all materials of the elasthyper-toolbox                         */
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateInvariantDerivatives(
    const LINALG::Matrix<3,1>& prinv,
    const int& gp,
    LINALG::Matrix<3,1>& dPIw,
    LINALG::Matrix<6,1>& ddPIIw,
    int eleGID
)

{
  // derivatives of principal materials weighted with their mass fraction in the constraint mixture
  LINALG::Matrix<3,1> dPgrowthI(true);
  LINALG::Matrix<6,1> ddPgrowthII(true);

  // loop map of associated potential summands
  // elastin matrix
  // derivatives of strain energy function w.r.t. principal invariants
  for(unsigned p=0;p<potsumeliso_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumeliso_[p]->AddDerivativesPrincipal(dPgrowthI,ddPgrowthII,prinv,eleGID);
    dPIw.Update(cur_rho_el_[gp],dPgrowthI,1.0);
    ddPIIw.Update(cur_rho_el_[gp],ddPgrowthII,1.0);
  }


  // derivatives of decoupled (volumetric or isochoric) materials weighted with their mass fraction in the constraint mixture
  LINALG::Matrix<3,1> modinv(true);
  InvariantsModified(modinv,prinv);
  LINALG::Matrix<3,1> dPmodI(true);
  LINALG::Matrix<6,1> ddPmodII(true);

  for(unsigned p=0;p<potsumeliso_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumeliso_[p]->AddDerivativesModified(dPgrowthI,ddPgrowthII,modinv,eleGID);
    dPmodI.Update(cur_rho_el_[gp],dPgrowthI,1.0);
    ddPmodII.Update(cur_rho_el_[gp],ddPgrowthII,1.0);
  }

  // volpenalty
  dPgrowthI.Clear();
  ddPgrowthII.Clear();
  if(params_->matid_penalty_ != -1)
    potsumelpenalty_->AddDerivativesModified(dPgrowthI,ddPgrowthII,modinv,0);
  dPmodI.Update(cur_rho_el_[gp],dPgrowthI,1.0);
  ddPmodII.Update(cur_rho_el_[gp],ddPgrowthII,1.0);

  // convert decoupled derivatives to principal derivatives
  ConvertModToPrinc(prinv,dPmodI,ddPmodII,dPIw,ddPIIw);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::InvariantsModified(
    LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<3,1>& prinv
    )
{
  // 1st invariant, trace
  modinv(0) = prinv(0)*std::pow(prinv(2),-1./3.);
  // 2nd invariant
  modinv(1) = prinv(1)*std::pow(prinv(2),-2./3.);
  // J
  modinv(2) = std::pow(prinv(2),1./2.);

  return;
}


/*----------------------------------------------------------------------/
 * Converts derivatives with respect to modified invariants in derivatives
 * with respect to principal invariants                                 */
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::ConvertModToPrinc(
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& dPmodI,
    const LINALG::Matrix<6,1>& ddPmodII,
    LINALG::Matrix<3,1>& dPI,
    LINALG::Matrix<6,1>& ddPII
    )
{

  // Conversions to dPI
  dPI(0) += std::pow(prinv(2),-1./3.)*dPmodI(0);
  dPI(1) += std::pow(prinv(2),-2./3.)*dPmodI(1);
  dPI(2) += 0.5*std::pow(prinv(2),-0.5)*dPmodI(2) - 1./3.*prinv(0)*std::pow(prinv(2),-4./3.)*dPmodI(0)
      - 2./3.*prinv(1)*std::pow(prinv(2),-5./3.)*dPmodI(1);

  // Conversions to ddPII
  ddPII(0) += std::pow(prinv(2),-2./3.)*ddPmodII(0);
  ddPII(1) += std::pow(prinv(2),-4./3.)*ddPmodII(1);
  ddPII(2) += (1./9.)*std::pow(prinv(2),-8./3.)*prinv(0)*prinv(0)*ddPmodII(0)
      + (4./9.)*prinv(0)*prinv(1)*std::pow(prinv(2),-3.)*ddPmodII(5) - (1./3.)*std::pow(prinv(2),-11./6.)*prinv(0)*ddPmodII(4)
      + (4./9.)*std::pow(prinv(2),-7./3.)*prinv(0)*dPmodI(0) + (4./9.)*std::pow(prinv(2),-10./3.)*prinv(1)*prinv(1)*ddPmodII(1)
      - (2./3.)*std::pow(prinv(2),-13./6.)*prinv(1)*ddPmodII(3) + (10./9.)*std::pow(prinv(2),-8./3.)*prinv(1)*dPmodI(1)
      + 0.25*std::pow(prinv(2),-1.)*ddPmodII(2) - 0.25*std::pow(prinv(2),-1.5)*dPmodI(2);
  ddPII(3) += -(1./3.)*std::pow(prinv(2),-2.)*prinv(0)*ddPmodII(5) - (2./3.)*std::pow(prinv(2),-7./3.)*prinv(1)*ddPmodII(1)
          + 0.5*std::pow(prinv(2),-7./6.)*ddPmodII(3) - (2./3.)*std::pow(prinv(2),-5./3.)*dPmodI(1);
  ddPII(4) += -(1./3.)*std::pow(prinv(2),-5./3.)*prinv(0)*ddPmodII(0) - (2./3.)*std::pow(prinv(2),-2.)*prinv(1)*ddPmodII(5)
  + 0.5*std::pow(prinv(2),-5./6.)*ddPmodII(4) - (1./3.)*std::pow(prinv(2),-4./3.)*dPmodI(0);
  ddPII(5) += std::pow(prinv(2),-1.)*ddPmodII(5);

  return;
}


/*----------------------------------------------------------------------
 * Evaluates the 2nd Piola-Kirchhoff Stress and constitutive tensor
 * with use of first and second derivatives according to invariants;
 * use principal calculation for all materials, as modified material
 * derivatives are converted before
 *                                                                      */
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateIsotropicStressCmat(
    LINALG::Matrix<6,1>& stress,
    LINALG::Matrix<6,6>& cmat,
    LINALG::Matrix<6,1> scg,
    LINALG::Matrix<6,1> id,
    LINALG::Matrix<6,1> icg,
    LINALG::Matrix<6,6> id4sharp,
    LINALG::Matrix<3,1> prinv,
    LINALG::Matrix<3,1> dPI,
    LINALG::Matrix<6,1> ddPII
    )
{
  // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  LINALG::Matrix<3,1> gamma(true);
  // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  LINALG::Matrix<8,1> delta(true);

  // compose coefficients
  CalculateGammaDelta(gamma,delta,prinv,dPI,ddPII);

  // 2nd Piola Kirchhoff stress
  stress.Update(gamma(0), id, 1.0);
  stress.Update(gamma(1), scg, 1.0);
  stress.Update(gamma(2), icg, 1.0);

  // constitutive tensor
  // contribution: Id \otimes Id
  cmat.MultiplyNT(delta(0), id, id, 1.0);
  // contribution: Id \otimes C + C \otimes Id
  cmat.MultiplyNT(delta(1), id, scg, 1.0);
  cmat.MultiplyNT(delta(1), scg, id, 1.0);
  // contribution: Id \otimes Cinv + Cinv \otimes Id
  cmat.MultiplyNT(delta(2), id, icg, 1.0);
  cmat.MultiplyNT(delta(2), icg, id, 1.0);
  // contribution: C \otimes C
  cmat.MultiplyNT(delta(3), scg, scg, 1.0);
  // contribution: C \otimes Cinv + Cinv \otimes C
  cmat.MultiplyNT(delta(4), scg, icg, 1.0);
  cmat.MultiplyNT(delta(4), icg, scg, 1.0);
  // contribution: Cinv \otimes Cinv
  cmat.MultiplyNT(delta(5), icg, icg, 1.0);
  // contribution: Cinv \odot Cinv
  AddtoCmatHolzapfelProduct(cmat, icg, delta(6));
  // contribution: Id4^#
  cmat.Update(delta(7), id4sharp, 1.0);

  return ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateIsotropicPrincElast(
    LINALG::Matrix<6,1>& stressisoprinc,
    LINALG::Matrix<6,6>& cmatisoprinc,
    LINALG::Matrix<6,1> iCg,
    LINALG::Matrix<6,1> iCgCiCg,
    LINALG::Matrix<6,1> iC,
    LINALG::Matrix<3,1> gamma,
    LINALG::Matrix<8,1> delta
    )
{
  // 2nd Piola Kirchhoff stresses
  stressisoprinc.Update(gamma(0), iCg, 1.0);
  stressisoprinc.Update(gamma(1), iCgCiCg, 1.0);
  stressisoprinc.Update(gamma(2), iC, 1.0);

  // constitutive tensor
  cmatisoprinc.MultiplyNT(delta(0),iCg,iCg,1.);
  cmatisoprinc.MultiplyNT(delta(1),iCgCiCg,iCg,1.);
  cmatisoprinc.MultiplyNT(delta(1),iCg,iCgCiCg,1.);
  cmatisoprinc.MultiplyNT(delta(2),iCg,iC,1.);
  cmatisoprinc.MultiplyNT(delta(2),iC,iCg,1.);
  cmatisoprinc.MultiplyNT(delta(3),iCgCiCg,iCgCiCg,1.);
  cmatisoprinc.MultiplyNT(delta(4),iCgCiCg,iC,1.);
  cmatisoprinc.MultiplyNT(delta(4),iC,iCgCiCg,1.);
  cmatisoprinc.MultiplyNT(delta(5),iC,iC,1.);
  AddtoCmatHolzapfelProduct(cmatisoprinc,iC,delta(6));
  AddtoCmatHolzapfelProduct(cmatisoprinc,iCg,delta(7));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateAdditionalTermsCmat(
    LINALG::Matrix<6,6>& cmatadd,
    LINALG::Matrix<3,1> gamma,
    LINALG::Matrix<8,1> delta,
    double v,
    int gp,
    LINALG::Matrix<3,3> AgM,
    LINALG::Matrix<1,6> sum_drhodC,
    LINALG::Matrix<3,3> iFgM,
    LINALG::Matrix<3,3> iCgC,
    LINALG::Matrix<3,3> iFgCe,
    LINALG::Matrix<6,1> iCg,
    LINALG::Matrix<9,1> CiFg,
    LINALG::Matrix<9,1> CiFgCe,
    LINALG::Matrix<6,1> iCgCiCg,
    LINALG::Matrix<9,1> CiFgiCe,
    LINALG::Matrix<6,1> iC
)
{
  //Initialize some variables
  LINALG::Matrix<6,9> dPK2diFin(true);

  LINALG::Matrix<3,3> id(true);
  for(int i=0;i<3;++i)
    id(i,i) = 1.0;

  // derivative of second Piola Kirchhoff stress w.r.t. inverse growth deformation gradient
  AddtodPK2diFin(dPK2diFin,id,iFgM,gamma(0));
  AddtodPK2diFin(dPK2diFin,iCgC,iFgM,gamma(1));
  dPK2diFin.MultiplyNT(delta(0),iCg,CiFg,1.);
  dPK2diFin.MultiplyNT(delta(1),iCg,CiFgCe,1.);
  dPK2diFin.MultiplyNT(delta(1),iCgCiCg,CiFg,1.);
  dPK2diFin.MultiplyNT(delta(2),iCg,CiFgiCe,1.);
  dPK2diFin.MultiplyNT(delta(2),iC,CiFg,1.);
  dPK2diFin.MultiplyNT(delta(3),iCgCiCg,CiFgCe,1.);
  dPK2diFin.MultiplyNT(delta(4),iCgCiCg,CiFgiCe,1.);
  dPK2diFin.MultiplyNT(delta(4),iC,CiFgCe,1.);
  dPK2diFin.MultiplyNT(delta(5),iC,CiFgiCe,1.);
  AddtodPK2diFin(dPK2diFin,id,iFgCe,0.5*delta(7));

  // diFinelastic/diFg
  LINALG::Matrix<3,3> tmp3x3(true);
  LINALG::Matrix<9,9> diFineldiFg(true);

  tmp3x3.UpdateT(1.0,GM_[gp],0.0);
  AddtoMatrixCrossProduct(1.0,id,tmp3x3,diFineldiFg);

  // diFg/dC
  LINALG::Matrix<9,6> diFgdC(true);
  LINALG::Matrix<9,1> Ag9x1(true);
  Matrix3x3to9x1(AgM,Ag9x1);

  switch(params_->growthtype_)
  {
  case 1:
  {
    diFgdC.MultiplyNN(-(1./(v*v))*(1./params_->density_),Ag9x1,sum_drhodC,0.0);

    break;
  }
  case 0:
  {
    LINALG::Matrix<9,1> Acir9x1(true);
    LINALG::Matrix<9,1> Aax9x1(true);
    Matrix3x3to9x1(AcirM_,Acir9x1);
    Matrix3x3to9x1(AaxM_,Aax9x1);
    diFgdC.MultiplyNN(-1./3.*std::pow(v,-4./3.)*(1./params_->density_),Ag9x1,sum_drhodC,0.0);
    diFgdC.MultiplyNN(-1./3.*std::pow(v,-4./3.)*(1./params_->density_),Acir9x1,sum_drhodC,1.0);
    diFgdC.MultiplyNN(-1./3.*std::pow(v,-4./3.)*(1./params_->density_),Aax9x1,sum_drhodC,0.0);

    break;
  }
  default:
    dserror("growthtype has to be either 1: anisotropic growth or 0: isotropic growth");
    break;
  }

  // update elasticity tensor
  LINALG::Matrix<6,9> tmp6x9(true);

  tmp6x9.MultiplyNN(1.0,dPK2diFin,diFineldiFg,0.0);
  cmatadd.MultiplyNN(2.0,tmp6x9,diFgdC,0.0);

  return;
}


/*----------------------------------------------------------------------*
 * Calculate the coefficients gamma and delta from the partial          *
 * derivatives w.r.t. invariants                                        *
 *----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::CalculateGammaDelta(
        LINALG::Matrix<3,1>& gamma,
        LINALG::Matrix<8,1>& delta,
        const LINALG::Matrix<3,1>& prinv,
        const LINALG::Matrix<3,1>& dPI,
        const LINALG::Matrix<6,1>& ddPII)
{
  // according to Holzapfel-Nonlinear Solid Mechanics p. 216
  gamma(0) = 2.*(dPI(0)+prinv(0)*dPI(1));
  gamma(1) = -2.*dPI(1);
  gamma(2) = 2.*prinv(2)*dPI(2);

  // according to Holzapfel-Nonlinear Solid Mechanics p. 261
  delta(0) = 4.*(ddPII(0) +2.*prinv(0)*ddPII(5) +dPI(1) +prinv(0)*prinv(0)*ddPII(1));
  delta(1) = -4.*(ddPII(5) +prinv(0)*ddPII(1));
  delta(2) = 4.*(prinv(2)*ddPII(4) +prinv(0)*prinv(2)*ddPII(3));
  delta(3) = 4.*ddPII(1);
  delta(4) = -4.*prinv(2)*ddPII(3);
  delta(5) = 4.*(prinv(2)*dPI(2) +prinv(2)*prinv(2)*ddPII(2));
  delta(6) = -4.*prinv(2)*dPI(2);
  delta(7) = -4.*dPI(1);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateElastinDamage()
{
  double T_el = 101*365.25; // in days
  double t_dam = 40.0;

  for(unsigned gp=0;gp<cur_rho_el_.size();++gp)
    cur_rho_el_[gp] = init_rho_el_[gp]*(exp(-(t_tot_-params_->t_prestress_)/T_el) - 0.5*(1.0-exp(-(t_tot_-params_->t_prestress_)/t_dam))*exp(-0.5*(100.0*(fabs(gp_ax_[gp])-0.09))*(100.0*(fabs(gp_ax_[gp])-0.09))));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateGrowthDefGrad(LINALG::Matrix<3,3>& FgM,
                                                          LINALG::Matrix<3,3>& iFgM,
                                                          const int gp)
{
   double rho_col_sum = 0.0;
   // evaluate volume change
   for(unsigned p=0;p<potsumrf_.size();++p)
     for(unsigned k=0;k<potsumrf_[p]->GetNumFibers();++k)
       rho_col_sum += potsumrf_[p]->GetCurMassDensity(k,gp);
   v_[gp] = (rho_col_sum+cur_rho_el_[gp])/params_->density_;

   switch(params_->growthtype_)
   {
   case 1:
     // build growth and inverse growth deformation gradient for anisotropic growth
     FgM.Update(1.0,AplM_,0.0);
     FgM.Update(v_[gp],AgM_,1.0);
     iFgM.Invert(FgM);
     break;
   case 0:
     // build growth and inverse growth deformation gradient for isotropic growth
     FgM.Update(std::pow(v_[gp],1./3.),AcirM_,0.0);
     FgM.Update(std::pow(v_[gp],1./3.),AradM_,1.0);
     FgM.Update(std::pow(v_[gp],1./3.),AaxM_,1.0);
     iFgM.Invert(FgM);
     break;
   default:
     dserror("growthtype has to be either 1: anisotropic growth or 0: isotropic growth");
     break;
   }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::ConvertStressCmatLoctoGlob(LINALG::Matrix<3,1> stress_loc,
                                                               LINALG::Matrix<3,3> cmat_loc,
                                                               LINALG::Matrix<6,1>& stress_glob,
                                                               LINALG::Matrix<6,6>& cmat_glob)
{
  // convert stress to global coordinate system
  LINALG::Matrix<6,1> tmp6x1(true);
  tmp6x1(0) = stress_loc(0);
  tmp6x1(1) = stress_loc(1);
  tmp6x1(3) = stress_loc(2);

  stress_glob.MultiplyTN(1.0,Q_proj_glob_,tmp6x1,0.0);


  // convert cmat to global coordinate system
  LINALG::Matrix<6,6> tmp6x6_1(true);
  LINALG::Matrix<6,6> tmp6x6_2(true);
  tmp6x6_1(0,0) = cmat_loc(0,0);
  tmp6x6_1(0,1) = cmat_loc(0,1);
  tmp6x6_1(0,3) = cmat_loc(0,2);
  tmp6x6_1(1,0) = cmat_loc(1,0);
  tmp6x6_1(1,1) = cmat_loc(1,1);
  tmp6x6_1(1,3) = cmat_loc(1,2);
  tmp6x6_1(3,0) = cmat_loc(2,0);
  tmp6x6_1(3,1) = cmat_loc(2,1);
  tmp6x6_1(3,3) = cmat_loc(2,2);

  tmp6x6_2.MultiplyTN(1.0,Q_proj_glob_,tmp6x6_1,0.0);
  cmat_glob.MultiplyNN(1.0,tmp6x6_2,Q_proj_glob_,1.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::AddtoMatrixCrossProduct(const double& fac,
                                                            const LINALG::Matrix<3,3>& A,
                                                            const LINALG::Matrix<3,3>& B,
                                                            LINALG::Matrix<9,9>& X_ABCD)
{
  X_ABCD(0,0) = fac * A(0,0) * B(0,0);
  X_ABCD(0,3) = fac * A(0,0) * B(0,1);
  X_ABCD(0,7) = fac * A(0,0) * B(0,2);
  X_ABCD(0,4) = fac * A(0,1) * B(0,0);
  X_ABCD(0,1) = fac * A(0,1) * B(0,1);
  X_ABCD(0,5) = fac * A(0,1) * B(0,2);
  X_ABCD(0,8) = fac * A(0,2) * B(0,0);
  X_ABCD(0,6) = fac * A(0,2) * B(0,1);
  X_ABCD(0,2) = fac * A(0,2) * B(0,2);

  X_ABCD(3,0) = fac * A(0,0) * B(1,0);
  X_ABCD(3,3) = fac * A(0,0) * B(1,1);
  X_ABCD(3,7) = fac * A(0,0) * B(1,2);
  X_ABCD(3,4) = fac * A(0,1) * B(1,0);
  X_ABCD(3,1) = fac * A(0,1) * B(1,1);
  X_ABCD(3,5) = fac * A(0,1) * B(1,2);
  X_ABCD(3,8) = fac * A(0,2) * B(1,0);
  X_ABCD(3,6) = fac * A(0,2) * B(1,1);
  X_ABCD(3,2) = fac * A(0,2) * B(1,2);

  X_ABCD(7,0) = fac * A(0,0) * B(2,0);
  X_ABCD(7,3) = fac * A(0,0) * B(2,1);
  X_ABCD(7,7) = fac * A(0,0) * B(2,2);
  X_ABCD(7,4) = fac * A(0,1) * B(2,0);
  X_ABCD(7,1) = fac * A(0,1) * B(2,1);
  X_ABCD(7,5) = fac * A(0,1) * B(2,2);
  X_ABCD(7,8) = fac * A(0,2) * B(2,0);
  X_ABCD(7,6) = fac * A(0,2) * B(2,1);
  X_ABCD(7,2) = fac * A(0,2) * B(2,2);

  X_ABCD(4,0) = fac * A(1,0) * B(0,0);
  X_ABCD(4,3) = fac * A(1,0) * B(0,1);
  X_ABCD(4,7) = fac * A(1,0) * B(0,2);
  X_ABCD(4,4) = fac * A(1,1) * B(0,0);
  X_ABCD(4,1) = fac * A(1,1) * B(0,1);
  X_ABCD(4,5) = fac * A(1,1) * B(0,2);
  X_ABCD(4,8) = fac * A(1,2) * B(0,0);
  X_ABCD(4,6) = fac * A(1,2) * B(0,1);
  X_ABCD(4,2) = fac * A(1,2) * B(0,2);

  X_ABCD(1,0) = fac * A(1,0) * B(1,0);
  X_ABCD(1,3) = fac * A(1,0) * B(1,1);
  X_ABCD(1,7) = fac * A(1,0) * B(1,2);
  X_ABCD(1,4) = fac * A(1,1) * B(1,0);
  X_ABCD(1,1) = fac * A(1,1) * B(1,1);
  X_ABCD(1,5) = fac * A(1,1) * B(1,2);
  X_ABCD(1,8) = fac * A(1,2) * B(1,0);
  X_ABCD(1,6) = fac * A(1,2) * B(1,1);
  X_ABCD(1,2) = fac * A(1,2) * B(1,2);

  X_ABCD(5,0) = fac * A(1,0) * B(2,0);
  X_ABCD(5,3) = fac * A(1,0) * B(2,1);
  X_ABCD(5,7) = fac * A(1,0) * B(2,2);
  X_ABCD(5,4) = fac * A(1,1) * B(2,0);
  X_ABCD(5,1) = fac * A(1,1) * B(2,1);
  X_ABCD(5,5) = fac * A(1,1) * B(2,2);
  X_ABCD(5,8) = fac * A(1,2) * B(2,0);
  X_ABCD(5,6) = fac * A(1,2) * B(2,1);
  X_ABCD(5,2) = fac * A(1,2) * B(2,2);

  X_ABCD(8,0) = fac * A(2,0) * B(0,0);
  X_ABCD(8,3) = fac * A(2,0) * B(0,1);
  X_ABCD(8,7) = fac * A(2,0) * B(0,2);
  X_ABCD(8,4) = fac * A(2,1) * B(0,0);
  X_ABCD(8,1) = fac * A(2,1) * B(0,1);
  X_ABCD(8,5) = fac * A(2,1) * B(0,2);
  X_ABCD(8,8) = fac * A(2,2) * B(0,0);
  X_ABCD(8,6) = fac * A(2,2) * B(0,1);
  X_ABCD(8,2) = fac * A(2,2) * B(0,2);

  X_ABCD(6,0) = fac * A(2,0) * B(1,0);
  X_ABCD(6,3) = fac * A(2,0) * B(1,1);
  X_ABCD(6,7) = fac * A(2,0) * B(1,2);
  X_ABCD(6,4) = fac * A(2,1) * B(1,0);
  X_ABCD(6,1) = fac * A(2,1) * B(1,1);
  X_ABCD(6,5) = fac * A(2,1) * B(1,2);
  X_ABCD(6,8) = fac * A(2,2) * B(1,0);
  X_ABCD(6,6) = fac * A(2,2) * B(1,1);
  X_ABCD(6,2) = fac * A(2,2) * B(1,2);

  X_ABCD(2,0) = fac * A(2,0) * B(2,0);
  X_ABCD(2,3) = fac * A(2,0) * B(2,1);
  X_ABCD(2,7) = fac * A(2,0) * B(2,2);
  X_ABCD(2,4) = fac * A(2,1) * B(2,0);
  X_ABCD(2,1) = fac * A(2,1) * B(2,1);
  X_ABCD(2,5) = fac * A(2,1) * B(2,2);
  X_ABCD(2,8) = fac * A(2,2) * B(2,0);
  X_ABCD(2,6) = fac * A(2,2) * B(2,1);
  X_ABCD(2,2) = fac * A(2,2) * B(2,2);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::AddtodPK2diFin(
    LINALG::Matrix<6,9>& dPK2diFin,
    LINALG::Matrix<3,3> A,
    LINALG::Matrix<3,3> B,
    double fac)
{
  dPK2diFin(0,0) += 2 * fac * A(0,0) * B(0,0);
  dPK2diFin(0,3) += 2 * fac * A(0,0) * B(0,1);
  dPK2diFin(0,5) += 2 * fac * A(0,0) * B(0,2);
  dPK2diFin(0,6) += 2 * fac * A(0,1) * B(0,0);
  dPK2diFin(0,1) += 2 * fac * A(0,1) * B(0,1);
  dPK2diFin(0,4) += 2 * fac * A(0,1) * B(0,2);
  dPK2diFin(0,8) += 2 * fac * A(0,2) * B(0,0);
  dPK2diFin(0,7) += 2 * fac * A(0,2) * B(0,1);
  dPK2diFin(0,2) += 2 * fac * A(0,2) * B(0,2);

  dPK2diFin(1,0) += 2 * fac * A(1,0) * B(1,0);
  dPK2diFin(1,3) += 2 * fac * A(1,0) * B(1,1);
  dPK2diFin(1,5) += 2 * fac * A(1,0) * B(1,2);
  dPK2diFin(1,6) += 2 * fac * A(1,1) * B(1,0);
  dPK2diFin(1,1) += 2 * fac * A(1,1) * B(1,1);
  dPK2diFin(1,4) += 2 * fac * A(1,1) * B(1,2);
  dPK2diFin(1,8) += 2 * fac * A(1,2) * B(1,0);
  dPK2diFin(1,7) += 2 * fac * A(1,2) * B(1,1);
  dPK2diFin(1,2) += 2 * fac * A(1,2) * B(1,2);

  dPK2diFin(2,0) += 2 * fac * A(2,0) * B(2,0);
  dPK2diFin(2,3) += 2 * fac * A(2,0) * B(2,1);
  dPK2diFin(2,5) += 2 * fac * A(2,0) * B(2,2);
  dPK2diFin(2,6) += 2 * fac * A(2,1) * B(2,0);
  dPK2diFin(2,1) += 2 * fac * A(2,1) * B(2,1);
  dPK2diFin(2,4) += 2 * fac * A(2,1) * B(2,2);
  dPK2diFin(2,8) += 2 * fac * A(2,2) * B(2,0);
  dPK2diFin(2,7) += 2 * fac * A(2,2) * B(2,1);
  dPK2diFin(2,2) += 2 * fac * A(2,2) * B(2,2);

  dPK2diFin(3,0) += fac * (A(0,0) * B(1,0) + A(1,0) * B(0,0));
  dPK2diFin(3,3) += fac * (A(0,0) * B(1,1) + A(1,0) * B(0,1));
  dPK2diFin(3,5) += fac * (A(0,0) * B(1,2) + A(1,0) * B(0,2));
  dPK2diFin(3,6) += fac * (A(0,1) * B(1,0) + A(1,1) * B(0,0));
  dPK2diFin(3,1) += fac * (A(0,1) * B(1,1) + A(1,1) * B(0,1));
  dPK2diFin(3,4) += fac * (A(0,1) * B(1,2) + A(1,1) * B(0,2));
  dPK2diFin(3,8) += fac * (A(0,2) * B(1,0) + A(1,2) * B(0,0));
  dPK2diFin(3,7) += fac * (A(0,2) * B(1,1) + A(1,2) * B(0,1));
  dPK2diFin(3,2) += fac * (A(0,2) * B(1,2) + A(1,2) * B(0,2));

  dPK2diFin(4,0) += fac * (A(1,0) * B(2,0) + A(2,0) * B(1,0));
  dPK2diFin(4,3) += fac * (A(1,0) * B(2,1) + A(2,0) * B(1,1));
  dPK2diFin(4,5) += fac * (A(1,0) * B(2,2) + A(2,0) * B(1,2));
  dPK2diFin(4,6) += fac * (A(1,1) * B(2,0) + A(2,1) * B(1,0));
  dPK2diFin(4,1) += fac * (A(1,1) * B(2,1) + A(2,1) * B(1,1));
  dPK2diFin(4,4) += fac * (A(1,1) * B(2,2) + A(2,1) * B(1,2));
  dPK2diFin(4,8) += fac * (A(1,2) * B(2,0) + A(2,2) * B(1,0));
  dPK2diFin(4,7) += fac * (A(1,2) * B(2,1) + A(2,2) * B(1,1));
  dPK2diFin(4,2) += fac * (A(1,2) * B(2,2) + A(2,2) * B(1,2));

  dPK2diFin(5,0) += fac * (A(0,0) * B(2,0) + A(2,0) * B(0,0));
  dPK2diFin(5,3) += fac * (A(0,0) * B(2,1) + A(2,0) * B(0,1));
  dPK2diFin(5,5) += fac * (A(0,0) * B(2,2) + A(2,0) * B(0,2));
  dPK2diFin(5,6) += fac * (A(0,1) * B(2,0) + A(2,1) * B(0,0));
  dPK2diFin(5,1) += fac * (A(0,1) * B(2,1) + A(2,1) * B(0,1));
  dPK2diFin(5,4) += fac * (A(0,1) * B(2,2) + A(2,1) * B(0,2));
  dPK2diFin(5,8) += fac * (A(0,2) * B(2,0) + A(2,2) * B(0,0));
  dPK2diFin(5,7) += fac * (A(0,2) * B(2,1) + A(2,2) * B(0,1));
  dPK2diFin(5,2) += fac * (A(0,2) * B(2,2) + A(2,2) * B(0,2));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Matrix3x3to9x1(LINALG::Matrix<3,3> A, LINALG::Matrix<9,1>& Out)
{
  Out(0) = A(0,0);
  Out(1) = A(1,1);
  Out(2) = A(2,2);
  Out(3) = A(0,1);
  Out(4) = A(1,2);
  Out(5) = A(0,2);
  Out(6) = A(1,0);
  Out(7) = A(2,1);
  Out(8) = A(2,0);
  return;
}


/*-------------------------------------------------------------------------------*
 |  evaluate isotropic and anisotropic stress and cmat (membrane formulation)    |
 *-------------------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateStressCmatMembrane(
    LINALG::Matrix<3,3>&              C_loc,
    const LINALG::Matrix<3,3>&        iFin_loc,
    LINALG::Matrix<3,1>&              stress_loc,
    LINALG::Matrix<3,3>&              cmat_loc,
    const int                         gp,
    const int                         eleGID
    )
{
  // some quantities
  LINALG::Matrix<3,1> prinv(true);      // elastic invariants

  LINALG::Matrix<3,1> iC(true);         // inverse right cauchygreen
  LINALG::Matrix<3,1> iCrg(true);       // inverse growth right cauchygreen
  LINALG::Matrix<3,1> iCrgCiCrg(true);

  double rcg_e_33 = 0.0;                // entry 33 of elastic right cauchygreen tensor

  LINALG::Matrix<3,1> dPI(true);       // first derivatives of strain energy functin w.r.t. invariants weightes with the corresponding volume fraction
  LINALG::Matrix<6,1> ddPII(true);     // first derivatives of strain energy functin w.r.t. invariants weightes with the corresponding volume fraction

  LINALG::Matrix<3,1> gamma_mem(true);  // 2nd Piola Kirchhoff stress factors according to Fakhreddine2011 equation (11)
  LINALG::Matrix<8,1> delta_mem(true);  // constitutive tensor factors according to Fakhreddine2011 equation (15)


  /*===============================================================================*
   | determine isotropic stress and elasticity tensor for elastin matrix           |
   *===============================================================================*/
  // evaluate elastic kinematic quatities
  EvaluateKinQuantElastMembrane(iFin_loc,prinv,iC,iCrg,iCrgCiCrg,C_loc,rcg_e_33);

  // first and second derivatives w.r.t. invariants weighted with the corresponding volume fraction (assuming 2D NeoHooke material)
  EvaluateInvariantDerivativesMembrane(prinv,gp,dPI,ddPII,eleGID);

  // compose coefficients
  CalculateGammaDeltaMembrane(gamma_mem,delta_mem,prinv,dPI,ddPII,rcg_e_33);

  // isotropic elastic stress response and elasticity tensor
  EvaluateIsotropicPrincElastMembrane(stress_loc,cmat_loc,iCrg,iCrgCiCrg,iC,gamma_mem,delta_mem);

  return;
}


/*-------------------------------------------------------------------------------*
 |  evaluate kinematic elastic quantities (membrane formulation)                 |
 *-------------------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateKinQuantElastMembrane(
    const LINALG::Matrix<3,3>& iFin_loc,
    LINALG::Matrix<3,1>& prinv,
    LINALG::Matrix<3,1>& iC_loc,
    LINALG::Matrix<3,1>& iCin_loc,
    LINALG::Matrix<3,1>& iCinCiCin_loc,
    LINALG::Matrix<3,3>& C_loc,
    double& rcg_e_33
    )
{
  // temporary variable
  LINALG::Matrix<3,3> tmp(true);

  // elastic right Cauchy-Green in matrix notation
  LINALG::Matrix<3,3> CeM_loc(true);
  tmp.Clear();
  tmp.MultiplyNN(1.0,C_loc,iFin_loc,0.0);
  CeM_loc.MultiplyTN(1.0,iFin_loc,tmp,0.0);

  // impose incompressibility for component in thickness direction of elastic right Cauchy-Green
  rcg_e_33 = 1.0/(CeM_loc(0,0)*CeM_loc(1,1)-CeM_loc(0,1)*CeM_loc(1,0));
  CeM_loc(2,2) = rcg_e_33;

  // inverse remodel growth right Cauchy-Green in local matrix notation
  LINALG::Matrix<3,3> iCinM_loc(true);
  iCinM_loc.MultiplyNT(iFin_loc,iFin_loc);

  // right Cauchy-Green in matrix notation
  LINALG::Matrix<3,3> Fin_loc(true);

  Fin_loc.Invert(iFin_loc);
  tmp.MultiplyNN(1.0,CeM_loc,Fin_loc,0.0);
  C_loc.MultiplyTN(1.0,Fin_loc,tmp,0.0);

  // inverse right Cauchy-Green in matrix notation
  LINALG::Matrix<3,3> iCM_loc(true);
  iCM_loc.Invert(C_loc);

  // C_in^-1 * C * C_in^-1
  tmp.Clear();
  LINALG::Matrix<3,3> iCinCiCinM_loc(true);
  tmp.MultiplyNN(1.0,iCinM_loc,C_loc,0.0);
  iCinCiCinM_loc.MultiplyNN(1.0,tmp,iCinM_loc,0.0);

  // stress-like voigt notation in two dimensions
  iCinCiCin_loc(0) = iCinCiCinM_loc(0,0);
  iCinCiCin_loc(1) = iCinCiCinM_loc(1,1);
  iCinCiCin_loc(2) = 0.5*(iCinCiCinM_loc(0,1)+iCinCiCinM_loc(1,0));

  // inverse right Cauchy-Green in stress-like voigt notation
  iC_loc(0) = iCM_loc(0,0);
  iC_loc(1) = iCM_loc(1,1);
  iC_loc(2) = 0.5*(iCM_loc(0,1)+iCM_loc(1,0));

  // inverse remodel growth right Cauchy-Green in stress-like voigt notation
  iCin_loc(0) = iCinM_loc(0,0);
  iCin_loc(1) = iCinM_loc(1,1);
  iCin_loc(2) = 0.5*(iCinM_loc(0,1)+iCinM_loc(1,0));

  // principal (isotropic) elastic invariants
  prinv(0) = CeM_loc(0,0)+CeM_loc(1,1)+CeM_loc(2,2);
  prinv(1) = 0.5*(std::pow(prinv(0),2.0) - (std::pow(CeM_loc(0,0),2.0) + std::pow(CeM_loc(1,1),2.0) + std::pow(CeM_loc(2,2),2.0) + 2.0*std::pow(CeM_loc(0,1),2.0)));
  prinv(2) = 1.0; // incompressibility condition

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateInvariantDerivativesMembrane(
    const LINALG::Matrix<3,1>& prinv,
    const int& gp,
    LINALG::Matrix<3,1>& dPI,
    LINALG::Matrix<6,1>& ddPII,
    int eleGID
)
{
  // derivatives of principal materials weighted with their mass fraction in the constraint mixture
  LINALG::Matrix<3,1> dPgrowthI(true);
  LINALG::Matrix<6,1> ddPgrowthII(true);

  // loop map of associated potential summands
  // elastin matrix
  // derivatives of strain energy function w.r.t. principal invariants
  for(unsigned p=0;p<potsumelmem_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumelmem_[p]->AddDerivativesPrincipal(dPgrowthI,ddPgrowthII,prinv,eleGID);
    dPI.Update(cur_rho_el_[gp]*mue_frac_[gp],dPgrowthI,1.0);
    ddPII.Update(cur_rho_el_[gp]*mue_frac_[gp],ddPgrowthII,1.0);
  }

  // as we use an incompressible membrane formulation derivation w.r.t. the third invariant have to be zero
  dPI(2) = 0.0;
  ddPII(2) = 0.0;
  ddPII(3) = 0.0;
  ddPII(4) = 0.0;

  // derivatives of decoupled (volumetric or isochoric) materials weighted with their mass fraction in the constraint mixture
  LINALG::Matrix<3,1> modinv(true);
  InvariantsModified(modinv,prinv);
  LINALG::Matrix<3,1> dPmodI(true);
  LINALG::Matrix<6,1> ddPmodII(true);

  for(unsigned p=0;p<potsumelmem_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumelmem_[p]->AddDerivativesModified(dPgrowthI,ddPgrowthII,modinv,eleGID);
    dPmodI.Update(cur_rho_el_[gp]*mue_frac_[gp],dPgrowthI,1.0);
    ddPmodII.Update(cur_rho_el_[gp]*mue_frac_[gp],ddPgrowthII,1.0);
  }

  // as we use an incompressible membrane formulation derivation w.r.t. the third invariant have to be zero
  dPmodI(2) = 0.0;
  ddPmodII(2) = 0.0;
  ddPmodII(3) = 0.0;
  ddPmodII(4) = 0.0;

  // convert decoupled derivatives to principal derivatives
  ConvertModToPrinc(prinv,dPmodI,ddPmodII,dPI,ddPII);

  return;
}


/*-------------------------------------------------------------------------------------*
 |  gamma and delta from partial derivatives w.r.t. invariants (membrane formulation)  |
 *-------------------------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::CalculateGammaDeltaMembrane(
    LINALG::Matrix<3,1>& gamma_mem,
    LINALG::Matrix<8,1>& delta_mem,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& dPI,
    const LINALG::Matrix<6,1>& ddPII,
    const double& rcg33)
{
  // compose coefficients
  // according to Fakhreddine2011 equation (11)
  gamma_mem(0) = 2.0*(dPI(0) + prinv(0)*dPI(1));
  gamma_mem(1) = -2.0*dPI(1);
  gamma_mem(2) = -rcg33*gamma_mem(0) - rcg33*rcg33*gamma_mem(1);

  // according to Fakhreddine2011 equation (15)
  delta_mem(0) = 4.0*(ddPII(0) + 2.0*prinv(0)*ddPII(5) + dPI(1) + prinv(0)*prinv(0)*ddPII(1));
  delta_mem(1) = -4.0*(ddPII(5) + prinv(0)*ddPII(1));
  delta_mem(2) = -4.0*rcg33*(ddPII(0) + prinv(0)*ddPII(5) + dPI(1) + (prinv(0)-rcg33)*(ddPII(5)+prinv(0)*ddPII(1)));
  delta_mem(3) = 4.0*ddPII(1);
  delta_mem(4) = 4.0*rcg33*(ddPII(5) + (prinv(0)-rcg33)*ddPII(1));
  delta_mem(5) = -2.0*gamma_mem(2) + 4.0*rcg33*rcg33*(ddPII(0) + 2.0*(prinv(0)-rcg33)*ddPII(5) + std::pow((prinv(0)-rcg33),2.0)*ddPII(1));
  delta_mem(6) = -2.0*gamma_mem(2);
  delta_mem(7) = 2.0*gamma_mem(1);

  return;
}


/*------------------------------------------------------------------------------------------*
 |  isotropic stress and elasticity tensor for coupled configuration (membrane formulation) |
 *------------------------------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateIsotropicPrincElastMembrane(
    LINALG::Matrix<3,1>& stress_mem,
    LINALG::Matrix<3,3>& cmat_mem,
    LINALG::Matrix<3,1>& iCrg,
    LINALG::Matrix<3,1>& iCrgCiCrg,
    LINALG::Matrix<3,1>& iC,
    LINALG::Matrix<3,1>& gamma_mem,
    LINALG::Matrix<8,1>& delta_mem
    )
{
  // clear entries
  stress_mem.Clear();
  cmat_mem.Clear();

  // isotropic 2nd Piola Kirchhoff stress
  stress_mem.Update(gamma_mem(0),iCrg,1.0);
  stress_mem.Update(gamma_mem(1),iCrgCiCrg,1.0);
  stress_mem.Update(gamma_mem(2),iC,1.0);

  // isotropic constitutive tensor
  // delta 1
  cmat_mem.MultiplyNT(delta_mem(0), iCrg, iCrg, 1.0);
  // delta 2
  cmat_mem.MultiplyNT(delta_mem(1), iCrg, iCrgCiCrg, 1.0);
  cmat_mem.MultiplyNT(delta_mem(1), iCrgCiCrg, iCrg, 1.0);
  // delta 3
  cmat_mem.MultiplyNT(delta_mem(2), iCrg, iC, 1.0);
  cmat_mem.MultiplyNT(delta_mem(2), iC, iCrg, 1.0);
  // delta 4
  cmat_mem.MultiplyNT(delta_mem(3), iCrgCiCrg, iCrgCiCrg, 1.0);
  // delta 5
  cmat_mem.MultiplyNT(delta_mem(4), iCrgCiCrg, iC, 1.0);
  cmat_mem.MultiplyNT(delta_mem(4), iC, iCrgCiCrg, 1.0);
  // delta 6
  cmat_mem.MultiplyNT(delta_mem(5), iC, iC, 1.0);
  // delta 7
  cmat_mem(0,0) += delta_mem(6) * iC(0)*iC(0);
  cmat_mem(0,1) += delta_mem(6) * iC(2)*iC(2);
  cmat_mem(0,2) += delta_mem(6) * iC(0)*iC(2);
  cmat_mem(1,0) += delta_mem(6) * iC(2)*iC(2);
  cmat_mem(1,1) += delta_mem(6) * iC(1)*iC(1);
  cmat_mem(1,2) += delta_mem(6) * iC(1)*iC(2);
  cmat_mem(2,0) += delta_mem(6) * iC(0)*iC(2);
  cmat_mem(2,1) += delta_mem(6) * iC(1)*iC(2);
  cmat_mem(2,2) += delta_mem(6) * 0.5 * (iC(0)*iC(1)+iC(2)*iC(2));
  // delta 8
  cmat_mem(0,0) += delta_mem(7) * iCrg(0)*iCrg(0);
  cmat_mem(0,1) += delta_mem(7) * iCrg(2)*iCrg(2);
  cmat_mem(0,2) += delta_mem(7) * iCrg(0)*iCrg(2);
  cmat_mem(1,0) += delta_mem(7) * iCrg(2)*iCrg(2);
  cmat_mem(1,1) += delta_mem(7) * iCrg(1)*iCrg(1);
  cmat_mem(1,2) += delta_mem(7) * iCrg(1)*iCrg(2);
  cmat_mem(2,0) += delta_mem(7) * iCrg(0)*iCrg(2);
  cmat_mem(2,1) += delta_mem(7) * iCrg(1)*iCrg(2);
  cmat_mem(2,2) += delta_mem(7) * 0.5 * (iCrg(0)*iCrg(1)+iCrg(2)*iCrg(2));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::MatrixtoStressLikeVoigtNotation(LINALG::Matrix<3,3> A,
                                                                    LINALG::Matrix<6,1>& B)
{
  for(int i=0;i<3;++i)
    B(i) = A(i,i);
  B(3) = 0.5*(A(0,1)+A(1,0));
  B(4) = 0.5*(A(1,2)+A(2,1));
  B(5) = 0.5*(A(0,2)+A(2,0));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::StressLikeVoigtNotationtoMatrix(LINALG::Matrix<6,1> A,
                                                                    LINALG::Matrix<3,3>& B)
{
  for(int i=0;i<3;++i)
    B(i,i) = A(i);
  B(0,1) = B(1,0) = A(3);
  B(1,2) = B(2,1) = A(4);
  B(0,2) = B(2,0) = A(5);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateMembrane(LINALG::Matrix<3,3>& defgrd_loc,
                                                     LINALG::Matrix<3,3>& C_loc,
                                                     LINALG::Matrix<3,3>& Q_trafo,
                                                     Teuchos::ParameterList& params,
                                                     LINALG::Matrix<3,1>* stress,
                                                     LINALG::Matrix<3,3>* cmat,
                                                     const int eleGID)
{
  // check local time integration method
  if(params_->loctimeint_ != 0)
    dserror("So far, G&R with a membrane finite element formulation is only implemented with explicit local time integration method");

  // Evaluate transformation matrix for Voigt notation (glob <-> loc coordinate frame)
  SetupTrafoMatrices(Q_trafo);

  // blank resulting quantities
  stress->Clear();
  cmat->Clear();

  // current gauss point
  int gp = params.get<int>("gp");

  // get simulation time
  t_tot_ = params.get<double>("total time");

  // identity matrix
  LINALG::Matrix<3,3> id(true);
  id(0,0) = id(1,1) = id(2,2) = 1.0;

  // growth deformation gradient
  LINALG::Matrix<3,3> Fg_glob(true);
  LINALG::Matrix<3,3> iFg_glob(true);
  if(first_[gp] == 1)
  {
    // gp coordinates in reference configuration
    LINALG::Matrix<1,3> gprefecoord(true);
    gprefecoord = params.get<LINALG::Matrix<1,3> >("gprefecoord");
    gp_ax_[gp] = cylcoords_[0].Dot(gprefecoord);

    // update elastin prestretch in radial direction
    GM_[gp].Update(params_->lamb_prestretch_cir_,AcirM_,0.0);
    GM_[gp].Update(params_->lamb_prestretch_ax_,AaxM_,1.0);
    GM_[gp].Update(1./(params_->lamb_prestretch_ax_*params_->lamb_prestretch_cir_),AradM_,1.0);
    first_[gp] = 0;
  }

  // Evaluate growth deformation gradient
  EvaluateGrowthDefGrad(Fg_glob,iFg_glob,gp);

  // total inelastic deformation gradient for elastin in local system
  LINALG::Matrix<3,3> Fin_loc(true);
  LINALG::Matrix<3,3> iGM(true);
  LINALG::Matrix<3,3> Fin_glob(true);
  LINALG::Matrix<6,1> Finv_loc(true);
  LINALG::Matrix<6,1> Finv_glob(true);

  iGM.Invert(GM_[gp]);
  Fin_glob.Multiply(1.0,iGM,Fg_glob,0.0);
  MatrixtoStressLikeVoigtNotation(Fin_glob,Finv_glob);
  Finv_loc.MultiplyNN(1.0,Q_proj_loc_,Finv_glob,0.0);
  StressLikeVoigtNotationtoMatrix(Finv_loc,Fin_loc);



 /*=============================================================================*
  | determine isotropic stress and elasticity tensor                            |
  *=============================================================================*/
  // total inelastic deformation gradient for elastin in local system
  LINALG::Matrix<3,1> stress_loc(true);
  LINALG::Matrix<3,3> cmat_loc(true);
  LINALG::Matrix<3,3> iFin_loc(true);
  iFin_loc.Invert(Fin_loc);

  // local right Cauchy-Green tensor
  EvaluateStressCmatMembrane(C_loc,iFin_loc,stress_loc,cmat_loc,gp,eleGID);

  stress->Update(1.0,stress_loc,1.0);
  cmat->Update(1.0,cmat_loc,1.0);

 /*===============================================================================*
  | determine anisotropic stress and elasticity tensor                            |
  *===============================================================================*/
  // stress response and elasticity tensor for remodeling of collagen fiber
  LINALG::Matrix<6,1> stress_aniso_glob(true);
  LINALG::Matrix<6,6> cmat_aniso_glob(true);

  // updated right Cauchy Green tensor in global coordinates
  LINALG::Matrix<3,3> CM_glob(true);
  LINALG::Matrix<6,1> Cv_loc(true);
  LINALG::Matrix<6,1> Cv_glob(true);
  MatrixtoStressLikeVoigtNotation(C_loc,Cv_loc);
  Cv_glob.MultiplyTN(1.0,Q_proj_glob_,Cv_loc,0.0);
  StressLikeVoigtNotationtoMatrix(Cv_glob,CM_glob);


  // call evaluation of stress and cmat for all fiber types
  LINALG::Matrix<6,1> stress_aniso_loc(true);
  LINALG::Matrix<6,6> tmp(true);
  LINALG::Matrix<6,6> cmat_aniso_loc(true);
  LINALG::Matrix<3,1> stress_aniso_loc_red(true);
  LINALG::Matrix<3,3> cmat_aniso_loc_red(true);

  for(unsigned p=0;p<potsumrf_.size();++p)
  {
    potsumrf_[p]->EvaluateAnisotropicStressCmat(CM_glob,iFg_glob,cmat_aniso_glob,stress_aniso_glob,gp,eleGID);

    // update
    stress_aniso_loc.MultiplyNN(1.0,Q_proj_loc_,stress_aniso_glob,0.0);
    tmp.MultiplyNN(1.0,Q_proj_loc_,cmat_aniso_glob,0.0);
    cmat_aniso_loc.MultiplyNT(1.0,tmp,Q_proj_loc_,0.0);

    // only use the in plane stress and stiffness
    stress_aniso_loc_red(0) = stress_aniso_loc(0);
    stress_aniso_loc_red(1) = stress_aniso_loc(1);
    stress_aniso_loc_red(2) = stress_aniso_loc(3);

    cmat_aniso_loc_red(0,0) = cmat_aniso_loc(0,0);
    cmat_aniso_loc_red(0,1) = cmat_aniso_loc(0,1);
    cmat_aniso_loc_red(0,2) = cmat_aniso_loc(0,3);
    cmat_aniso_loc_red(1,0) = cmat_aniso_loc(1,0);
    cmat_aniso_loc_red(1,1) = cmat_aniso_loc(1,1);
    cmat_aniso_loc_red(1,2) = cmat_aniso_loc(1,3);
    cmat_aniso_loc_red(2,0) = cmat_aniso_loc(3,0);
    cmat_aniso_loc_red(2,1) = cmat_aniso_loc(3,1);
    cmat_aniso_loc_red(2,2) = cmat_aniso_loc(3,3);

    stress->Update(1.0,stress_aniso_loc_red,1.0);
    cmat->Update(1.0,cmat_aniso_loc_red,1.0);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::VisNames(std::map<std::string,int>& names)
{
  std::string result_mass_fraction_el;
  result_mass_fraction_el = "mass_fraction_el";

  names[result_mass_fraction_el] = 1;


  std::string result_v_growth;
  result_v_growth = "v_growth";

  names[result_v_growth] = 1;


  // loop map of associated potential summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->VisNames(names,p);

  // 3D elastin matrix
  for (unsigned int p=0; p<potsumeliso_.size(); ++p)
    potsumeliso_[p]->VisNames(names);

  // 2D elastin matrix
  for (unsigned int p=0; p<potsumelmem_.size(); ++p)
    potsumelmem_[p]->VisNames(names);

  // volpenalty
  if(params_->matid_penalty_ != -1)
    potsumelpenalty_->VisNames(names);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::GrowthRemodel_ElastHyper::VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  int return_val= 0;

  if(name == "mass_fraction_el")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned i=0;i<cur_rho_el_.size();++i)
    {
      data[0] += cur_rho_el_[i];
    }
    data[0] = data[0]/cur_rho_el_.size();

    return true;
  }

  if(name == "v_growth")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned i=0;i<v_.size();++i)
    {
      data[0] += v_[i];
    }
    data[0] = data[0]/v_.size();

    return true;
  }

  // loop map of associated potential summands
  // remodelfiber
  if(name.at( name.length() - 3 ) == '0')
    return_val+=potsumrf_[0]->VisData(name,data,numgp,eleID);
  if(name.at( name.length() - 3 ) == '1')
    return_val+=potsumrf_[1]->VisData(name,data,numgp,eleID);

  // 3D elastin matrix
  for (unsigned int p=0; p<potsumeliso_.size(); ++p)
    return_val+=potsumeliso_[p]->VisData(name,data,numgp,eleID);

  // 2D elastin matrix
  for (unsigned int p=0; p<potsumelmem_.size(); ++p)
    return_val+=potsumelmem_[p]->VisData(name,data,numgp,eleID);

  // volpenalty
  if(params_->matid_penalty_ != -1)
    return_val+=potsumelpenalty_->VisData(name,data,numgp,eleID);

  return (bool)return_val;
}
