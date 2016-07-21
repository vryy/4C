/*----------------------------------------------------------------------*/
/*!
\file growthremodel_elasthyper.cpp
\brief This file is used to manage the homogenized constraint mixture during growth and remodeling

The input line should read
MAT 0   MAT_Growthremodel_ElastHyper   NUMMATRF 0 MATIDSRF NUMMATEL 0 MATIDSEL NUMMATGR 0 MATIDSGR MATIDPENALTY ELMASSFRAC GRMASSFRAC DENS 0

\level 3

<pre>
\maintainer Fabian Bräu
            braeu@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

/*----------------------------------------------------------------------*/

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
  nummat_elastin_(matdata->GetInt("NUMMATEL")),
  matids_remodelfiber_(matdata->Get<std::vector<int> >("MATIDSRF")),
  matids_elastin_(matdata->Get<std::vector<int> >("MATIDSEL")),
  matid_penalty_(matdata->GetInt("MATIDPENALTY")),
  init_w_el_(matdata->Get<std::vector<double> >("ELMASSFRAC")),
  density_(matdata->GetDouble("DENS")),
  t_prestress_(matdata->GetDouble("PRESTRESSTIME")),
  lamb_prestretch_cir_(matdata->GetDouble("PRESTRETCHELASTINCIR")),
  lamb_prestretch_ax_(matdata->GetDouble("PRESTRETCHELASTINAX")),
  damage_(matdata->GetInt("DAMAGE")),
  growthtype_(matdata->GetInt("GROWTHTYPE"))
{
  // check if sizes fit
  if (nummat_remodelfiber_ != (int)matids_remodelfiber_->size())
    dserror("number of remodelfiber materials %d does not fit to size of remodelfiber material vector %d", nummat_remodelfiber_, matids_remodelfiber_->size());

  if (nummat_elastin_ != (int)matids_elastin_->size())
    dserror("number of elastin materials %d does not fit to size of elastin material vector %d", nummat_elastin_, matids_elastin_->size());
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
    potsumel_(0),
    potsumpenalty_(0),
    t_tot_(0),
    t_ref_(0),
    p_mean_(0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::GrowthRemodel_ElastHyper::GrowthRemodel_ElastHyper(MAT::PAR::GrowthRemodel_ElastHyper* params)
  : params_(params),
    potsumrf_(0),
    potsumel_(0),
    potsumpenalty_(0),
    t_tot_(0),
    t_ref_(0),
    p_mean_(0)
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

  // Elastin matrix
  for (m=params_->matids_elastin_->begin(); m!=params_->matids_elastin_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumel_.push_back(sum);
  }
  for(unsigned p=0;p<potsumel_.size();++p)
    if(potsumel_[p]->MaterialType() != INPAR::MAT::mes_isoneohooke)
      dserror("So far, prestressing is only implemented for IsoNeoHooke material but it can be easily expanded!");

  // VolPenalty
  if(params_->matid_penalty_ != -1)
  {
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(params_->matid_penalty_);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumpenalty_ = sum;
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


  // mass fractions of elastin and ground matrix
  int num_el = 0;
  num_el = cur_rho_el_.size();

  AddtoPack(data,num_el);

  for(int i=0;i<num_el;++i)
  {
    AddtoPack(data,cur_rho_el_[i]);
    AddtoPack(data,init_rho_el_[i]);
  }

  AddtoPack(data,v_);
  AddtoPack(data,gp_ax_);
  AddtoPack(data,gp_rad_);
  AddtoPack(data,t_ref_);
  AddtoPack(data,p_mean_);
  AddtoPack(data,AcirM_);
  AddtoPack(data,AaxM_);
  AddtoPack(data,AradM_);
  AddtoPack(data,Aradv_);
  AddtoPack(data,GM_);
  AddtoPack(data,cylcoords_);

  if (params_ != NULL) // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumrf_.size(); ++p)
     potsumrf_[p]->PackSummand(data);

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumel_.size(); ++p)
     potsumel_[p]->PackSummand(data);

    if(params_->matid_penalty_ != -1)
      potsumpenalty_->PackSummand(data);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = NULL;
  potsumrf_.clear();
  potsumel_.clear();

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
  int num_el = 0;

  ExtractfromPack(position,data,num_el);
  cur_rho_el_.resize(num_el);
  init_rho_el_.resize(num_el);

  for(int i=0;i<num_el;++i)
  {
    ExtractfromPack(position,data,cur_rho_el_[i]);
    ExtractfromPack(position,data,init_rho_el_[i]);
  }

  ExtractfromPack(position,data,v_);
  ExtractfromPack(position,data,gp_ax_);
  ExtractfromPack(position,data,gp_rad_);
  ExtractfromPack(position,data,t_ref_);
  ExtractfromPack(position,data,p_mean_);
  ExtractfromPack(position,data,AcirM_);
  ExtractfromPack(position,data,AaxM_);
  ExtractfromPack(position,data,AradM_);
  ExtractfromPack(position,data,Aradv_);
  ExtractfromPack(position,data,GM_);
  ExtractfromPack(position,data,cylcoords_);

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

    // Elastin matrix
    for (m=params_->matids_elastin_->begin(); m!=params_->matids_elastin_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsumel_.push_back(sum);
    }

    // VolPenalty
    if(params_->matid_penalty_ != -1)
    {
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(params_->matid_penalty_);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsumpenalty_ = sum;
    }

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumrf_.size(); ++p)
     potsumrf_[p]->UnpackSummand(data,position);

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumel_.size(); ++p)
     potsumel_[p]->UnpackSummand(data,position);

    if(params_->matid_penalty_ != -1)
      potsumpenalty_->UnpackSummand(data,position);

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
  cur_rho_el_.resize(potsumel_.size());
  init_rho_el_.resize(potsumel_.size());
  GM_.resize(numgp,LINALG::Matrix<3,3>(true));


  for(unsigned p=0;p<potsumel_.size();++p)
  {
    init_rho_el_[p] = params_->density_ * params_->init_w_el_->at(p);
    cur_rho_el_[p].resize(numgp,init_rho_el_[p]);
  }


  // Setup summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->Setup(numgp,params_->density_,linedef);

  // elastin matrix
  for (unsigned int p=0; p<potsumel_.size(); ++p)
    potsumel_[p]->Setup(linedef);

  // volpenalty
  if(params_->matid_penalty_ != -1)
    potsumpenalty_->Setup(linedef);


  // Setup circumferential, radial and axial structural tensor
  SetupAxiCirRadStructuralTensor(linedef);

  // setup prestretch of elastin
  for(int gp=0;gp<numgp;++gp)
  {
    GM_[gp].Update(params_->lamb_prestretch_cir_,AcirM_,0.0);
    GM_[gp].Update(params_->lamb_prestretch_ax_,AaxM_,1.0);
    GM_[gp].Update(1./(params_->lamb_prestretch_cir_*params_->lamb_prestretch_ax_),AradM_,1.0);
  }

  // TODO: Make prestressing independent of the geometry and the material
  // setup wall thickness in reference configuration (used in prestressing routine)
  t_ref_ = 0.001682814382266;
  // setup mean blood pressure (used in prestressing routine)
  p_mean_ = 13332.2668;

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
void MAT::GrowthRemodel_ElastHyper::Update()
{
  // Update individual volume of elastin
  if(params_->damage_ == 1)
    EvaluateElastinDamage();

  // loop map of associated potential summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->Update();

  // elastin matrix
  for (unsigned int p=0; p<potsumel_.size(); ++p)
    potsumel_[p]->Update();

  // volpenalty
  if(params_->matid_penalty_ != -1)
    potsumpenalty_->Update();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateStressCmat(const LINALG::Matrix<3,3>* defgrd,
                                                       const int gp,
                                                       LINALG::Matrix<6,1>& stressiso,
                                                       LINALG::Matrix<6,6>& cmatiso,
                                                       const int eleGID)
{
  // Evaluate elastin matrix
  // some variables
  LINALG::Matrix<6,1> iCr(true);
  LINALG::Matrix<6,1> iCrCiCr(true);
  LINALG::Matrix<6,1> iC(true);
  LINALG::Matrix<3,1> prinv(true);
  LINALG::Matrix<3,3> iCrC(true);
  LINALG::Matrix<3,3> iFrCe(true);
  LINALG::Matrix<9,1> CiFr(true);
  LINALG::Matrix<9,1> CiFrCe(true);
  LINALG::Matrix<9,1> CiFriCe(true);

  // first and second derivatives w.r.t. invariants weighted with the corresponding volume fraction
  LINALG::Matrix<3,1> dPIe(true);
  LINALG::Matrix<6,1> ddPIIe(true);
  LINALG::Matrix<3,1> dPmodI(true);
  LINALG::Matrix<6,1> ddPmodII(true);


  EvaluateKinQuantElast(defgrd,GM_[gp],gp,iCr,iCrCiCr,iC,iCrC,iFrCe,CiFr,CiFrCe,CiFriCe,prinv);

  LINALG::Matrix<3,1> modinv(true);
  InvariantsModified(modinv,prinv);
  LINALG::Matrix<3,1> dPgrowthI(true);
  LINALG::Matrix<6,1> ddPgrowthII(true);
  // loop map of associated potential summands
  // elastin matrix
  // derivatives of strain energy function w.r.t. principal invariants
  for(unsigned p=0;p<potsumel_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumel_[p]->AddDerivativesPrincipal(dPgrowthI,ddPgrowthII,prinv,eleGID);
    dPIe.Update(cur_rho_el_[p][gp],dPgrowthI,1.0);
    ddPIIe.Update(cur_rho_el_[p][gp],ddPgrowthII,1.0);
  }

  // derivatives of strain energy function w.r.t. modified invariants
  for(unsigned p=0;p<potsumel_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumel_[p]->AddDerivativesModified(dPgrowthI,ddPgrowthII,modinv,0);
    dPmodI.Update(cur_rho_el_[p][gp],dPgrowthI,1.0);
    ddPmodII.Update(cur_rho_el_[p][gp],ddPgrowthII,1.0);
  }

  if(params_->matid_penalty_ != -1)
    potsumpenalty_->AddDerivativesModified(dPmodI,ddPmodII,modinv,0);

  // convert decoupled derivatives to principal derivatives
  ConvertModToPrinc(prinv,dPmodI,ddPmodII,dPIe,ddPIIe);

  // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  LINALG::Matrix<3,1> gamma(true);
  // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  LINALG::Matrix<8,1> delta(true);

  // compose coefficients
  CalculateGammaDelta(gamma,delta,prinv,dPIe,ddPIIe);

  EvaluateIsotropicPrincElast(stressiso,cmatiso,iCr,iCrCiCr,iC,gamma,delta);

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


  // evaluate remodelfiber (remodeling of collagen)
  // build stress response and elasticity tensor
  LINALG::Matrix<NUM_STRESS_3D,1> stressremodel(true) ;
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatremodel(true) ;

  // factors which are used for the calculation of volume change v
  double fac_el = 0.0;

  for(unsigned p=0;p<potsumel_.size();++p)
    fac_el += init_rho_el_[p];

  // sum over all fiber families of the derivatives of the individual mass density w.r.t. right Cauchy Green tensor
  LINALG::Matrix<1,6> sum_drhodC(true);

  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<3,3> id(true);
  for (int i=0; i<3; i++)
    id(i,i) = 1.0;

  // growth deformation gradient
  LINALG::Matrix<3,3> FgM(true);
  LINALG::Matrix<3,3> iFgM(true);
  FgM.Update(1.0,id,0.0);
  iFgM.Update(1.0,id,0.0);


  //TODO: Here we assume that all fibers are in one plane, so that the growth
  //      direction results from the cross product of two fiber directions.
  //      Later we should implement the evolution equation for the growth tensor but
  //      for now this is enough
  // structural tensor of growth direction
  LINALG::Matrix<3,3> AgM(true);
  AgM.Update(1.0,AradM_,0.0);

  // structural tensor of the plane in which all fibers are located
  LINALG::Matrix<3,3> AplM(true);
  AplM.Update(1.0,AradM_,0.0);
  AplM.Update(1.0,id,-1.0);


  // prestressing
  if(t_tot_ <= params_->t_prestress_)
  {
    // gp coordinates in reference configuration
    LINALG::Matrix<1,3> gprefecoord(true);
    gprefecoord = params.get<LINALG::Matrix<1,3> >("gprefecoord");
    gp_ax_[gp] = cylcoords_[0].Dot(gprefecoord);
    gp_rad_[gp] = cylcoords_[2].Dot(gprefecoord);


    for(unsigned p=0;p<potsumrf_.size();++p)
    {
      potsumrf_[p]->EvaluatePrestressing(defgrd,cmatremodel,stressremodel,gp,eleGID);
      stress->Update(1.0,stressremodel,1.0);
      cmat->Update(1.0,cmatremodel,1.0);
    }

    // build stress response and elasticity tensor
    LINALG::Matrix<NUM_STRESS_3D,1> stressiso(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true) ;

    EvaluateStressCmat(defgrd,gp,stressiso,cmatiso,eleGID);

    stress->Update(1.0,stressiso,1.0);
    cmat->Update(1.0,cmatiso,1.0);

    // update radial stress because of prestressing
    if(params_->matid_penalty_ != -1)
      stress->Update(-72.0*init_rho_el_[0]*std::pow(GM_[gp].Determinant(),-4./3.)*(std::pow(GM_[gp].Dot(AradM_),2.0)-
          (1./3.)*(GM_[gp].Dot(GM_[gp])))-(1.0-(gp_rad_[gp]-10.0e-3)/t_ref_)*p_mean_,Aradv_,1.0);
    else


    return;
  }
  // Growth and Remodeling
  else
  {
    // total number of fibers
    unsigned nr_f_tot = 0;

    for(unsigned p=0;p<potsumrf_.size();++p)
      nr_f_tot += potsumrf_[p]->GetNumFibers();

    double tmp = 0.0;
    // evaluate volume change
    for(unsigned p=0;p<potsumrf_.size();++p)
      for(unsigned k=0;k<potsumrf_[p]->GetNumFibers();++k)
        tmp += potsumrf_[p]->GetCurMassDensity(k,gp);
    v_[gp] = (tmp+fac_el)/params_->density_;

    switch(params_->growthtype_)
    {
    case 1:
      // build growth and inverse growth deformation gradient for anisotropic growth
      FgM.Update(1.0,AplM,0.0);
      FgM.Update(v_[gp],AradM_,1.0);
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

        tmp = 0.0;
        // evaluate volume change
        for(unsigned p=0;p<potsumrf_.size();++p)
          for(unsigned k=0;k<potsumrf_[p]->GetNumFibers();++k)
            tmp += potsumrf_[p]->GetCurMassDensity(k,gp);
        v_[gp] = (tmp+fac_el)/params_->density_;

        switch(params_->growthtype_)
        {
        case 1:
          // build growth and inverse growth deformation gradient for anisotropic growth
          FgM.Update(1.0,AplM,0.0);
          FgM.Update(v_[gp],AradM_,1.0);
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
      }

      // global number of fibers which were already processed
      nr_grf_proc = 0;
      for(unsigned p=0;p<potsumrf_.size();++p)
      {
        potsumrf_[p]->EvaluateDerivativesInternalNewton(defgrd,id,nr_grf_proc,params_->density_,gp,dt,v_[gp],FgM,iFgM,AgM,AcirM_,
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
      potsumrf_[p]->AddStressCmatGrowthRemodel(defgrd,id,nr_f_tot,nr_grf_proc,gp,v_[gp],params_->density_,iFgM,AgM,AcirM_,AradM_,AaxM_,
          drhodC,dlambdC,stressremodel,cmatremodel,eleGID,params_->growthtype_);
      stress->Update(1.0,stressremodel,1.0);
      cmat->Update(1.0,cmatremodel,1.0);
      nr_grf_proc += potsumrf_[p]->GetNumFibers();
    }
  }


  // Evaluate elastin and ground matrix
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
  LINALG::Matrix<3,3> iFinel(true);

  iFinel.MultiplyNN(1.0,iFgM,GM_[gp],0.0);
  EvaluateKinQuantElast(defgrd,iFinel,gp,iCg,iCgCiCg,iC,iCgC,iFgCe,CiFg,CiFgCe,CiFgiCe,prinv);
  EvaluateInvariantDerivatives(prinv,gp,dPIw,ddPIIw,eleGID);


  // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  LINALG::Matrix<3,1> gamma(true);
  // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  LINALG::Matrix<8,1> delta(true);

  // compose coefficients
  CalculateGammaDelta(gamma,delta,prinv,dPIw,ddPIIw);


  // build stress response and elasticity tensor
  LINALG::Matrix<NUM_STRESS_3D,1> stressiso(true) ;
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true) ;

  EvaluateIsotropicPrincElast(stressiso,cmatiso,iCg,iCgCiCg,iC,gamma,delta);

  stress->Update(1.0,stressiso,1.0);
  cmat->Update(1.0,cmatiso,1.0);


  // build additional terms for elasticity tensor
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatadd(true);
  EvaluateAdditionalTermsCmat(cmatadd,gamma,delta,v_[gp],gp,AgM,sum_drhodC,iFinel,iCgC,iFgCe,iCg,CiFg,CiFgCe,iCgCiCg,CiFgiCe,iC);

  cmat->Update(1.0,cmatadd,1.0);

  // update radial stress because of prestressing
  stress->Update(-72.0*init_rho_el_[0]*std::pow(GM_[gp].Determinant(),-4./3.)*(std::pow(GM_[gp].Dot(AradM_),2.0)-
      (1./3.)*(GM_[gp].Dot(GM_[gp])))-(1.0-(gp_rad_[gp]-10.0e-3)/t_ref_)*p_mean_,Aradv_,1.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateKinQuant(
    const LINALG::Matrix<6,1>& glstrain,
    LINALG::Matrix<6,1>& id,
    LINALG::Matrix<6,1>& C_voigt,
    LINALG::Matrix<6,1>& C_strain,
    LINALG::Matrix<6,1>& iC_voigt,
    LINALG::Matrix<6,6>& id4sharp,
    LINALG::Matrix<3,1>& prinv)

{
  // build Cartesian identity 2-tensor I_{AB}
  for (int i=0; i<3; i++) id(i) = 1.0;

  // right Cauchy-Green Tensor  C_{AB} = 2 * E_{AB} + I_{AB}
  // REMARK: strain-like 6-Voigt vector
  C_strain.Update(2.0,glstrain,1.0);
  C_strain.Update(1.0, id, 1.0);

  // 'contra-variant' right Cauchy-Green Tensor C^{AB}
  // REMARK: stress-like 6-Voigt vector of right CG
  C_voigt.Update(1.0,C_strain,1.0);
  for (int i=3; i<6; i++) C_voigt(i) *= 0.5;

  // principal invariants of right Cauchy-Green strain
  InvariantsPrincipal(prinv,C_strain);

  // invert right Cauchy-Green tensor
  // REMARK: stress-like 6-Voigt vector
  {
    iC_voigt(0) = ( C_strain(1)*C_strain(2) - 0.25*C_strain(4)*C_strain(4) ) / prinv(2);
    iC_voigt(1) = ( C_strain(0)*C_strain(2) - 0.25*C_strain(5)*C_strain(5) ) / prinv(2);
    iC_voigt(2) = ( C_strain(0)*C_strain(1) - 0.25*C_strain(3)*C_strain(3) ) / prinv(2);
    iC_voigt(3) = ( 0.25*C_strain(5)*C_strain(4) - 0.5*C_strain(3)*C_strain(2) ) / prinv(2);
    iC_voigt(4) = ( 0.25*C_strain(3)*C_strain(5) - 0.5*C_strain(0)*C_strain(4) ) / prinv(2);
    iC_voigt(5) = ( 0.25*C_strain(3)*C_strain(4) - 0.5*C_strain(5)*C_strain(1) ) / prinv(2);
  }

  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  for (int i=0; i<3; i++) id4sharp(i,i) = 1.0;
  for (int i=3; i<6; i++) id4sharp(i,i) = 0.5;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateKinQuantElast(
    const LINALG::Matrix<3,3>* defgrd,
    LINALG::Matrix<3,3> iFgM,
    const int gp,
    LINALG::Matrix<6,1>& iCg,
    LINALG::Matrix<6,1>& iCgCiCg,
    LINALG::Matrix<6,1>& iC,
    LINALG::Matrix<3,3>& iCgC,
    LINALG::Matrix<3,3>& iFgCe,
    LINALG::Matrix<9,1>& CiFg,
    LINALG::Matrix<9,1>& CiFgCe,
    LINALG::Matrix<9,1>& CiFgiCe,
    LINALG::Matrix<3,1>& prinv)
{
  // inverse growth right Cauchy-Green
  LINALG::Matrix<3,3> iCgM(true);
  iCgM.MultiplyNT(iFgM,iFgM);
  // stress-like Voigt notation
  for (int i=0; i<3; i++) iCg(i) = iCgM(i,i);
  iCg(3) = (iCgM(0,1)+iCgM(1,0))/2.;
  iCg(4) = (iCgM(2,1)+iCgM(1,2))/2.;
  iCg(5) = (iCgM(0,2)+iCgM(2,0))/2.;

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
  LINALG::Matrix<3,3> iCgCiCgM;
  tmp.Multiply(iCgM,RCG);
  iCgCiCgM.Multiply(tmp,iCgM);
  // stress-like Voigt notation
  for (int i=0; i<3; i++) iCgCiCg(i) = iCgCiCgM(i,i);
  iCgCiCg(3) = (iCgCiCgM(0,1)+iCgCiCgM(1,0))/2.;
  iCgCiCg(4) = (iCgCiCgM(2,1)+iCgCiCgM(1,2))/2.;
  iCgCiCg(5) = (iCgCiCgM(0,2)+iCgCiCgM(2,0))/2.;

  tmp.Multiply(*defgrd,iFgM);
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
  iCgC.MultiplyNN(1.0,iCgM,RCG,0.0);

  // F_g^-1 * C_e
  iFgCe.MultiplyNN(1.0,iFgM,CeM,0.0);

  // C * F_g^-1
  LINALG::Matrix<3,3> CiFgM(true);
  CiFgM.MultiplyNN(1.0,RCG,iFgM,0.0);
  Matrix3x3to9x1(CiFgM,CiFg);

  // C * F_g^-1 * C_e
  LINALG::Matrix<3,3> CiFgCeM(true);
  tmp.MultiplyNN(1.0,RCG,iFgM,0.0);
  CiFgCeM.MultiplyNN(1.0,tmp,CeM,0.0);
  Matrix3x3to9x1(CiFgCeM,CiFgCe);

  // C * F_g^-1 * C_e^-1
  LINALG::Matrix<3,3> CiFgiCeM(true);
  LINALG::Matrix<3,3> iCeM(true);
  iCeM.Invert(CeM);
  tmp.MultiplyNN(1.0,RCG,iFgM,0.0);
  CiFgiCeM.MultiplyNN(1.0,tmp,iCeM,0.0);
  Matrix3x3to9x1(CiFgiCeM,CiFgiCe);

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
  for(unsigned p=0;p<potsumel_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumel_[p]->AddDerivativesPrincipal(dPgrowthI,ddPgrowthII,prinv,eleGID);
    dPIw.Update(cur_rho_el_[p][gp],dPgrowthI,1.0);
    ddPIIw.Update(cur_rho_el_[p][gp],ddPgrowthII,1.0);
  }


  // derivatives of decoupled (volumetric or isochoric) materials weighted with their mass fraction in the constraint mixture
  LINALG::Matrix<3,1> modinv(true);
  InvariantsModified(modinv,prinv);
  LINALG::Matrix<3,1> dPmodI(true);
  LINALG::Matrix<6,1> ddPmodII(true);

  for(unsigned p=0;p<potsumel_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumel_[p]->AddDerivativesModified(dPgrowthI,ddPgrowthII,modinv,eleGID);
    dPmodI.Update(cur_rho_el_[p][gp],dPgrowthI,1.0);
    ddPmodII.Update(cur_rho_el_[p][gp],ddPgrowthII,1.0);
  }

  // volpenalty
  if(params_->matid_penalty_ != -1)
    potsumpenalty_->AddDerivativesModified(dPmodI,ddPmodII,modinv,eleGID);

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
  LINALG::Matrix<6,9> dPK2diFinel(true);

  LINALG::Matrix<3,3> id(true);
  for(int i=0;i<3;++i)
    id(i,i) = 1.0;

  // derivative of second Piola Kirchhoff stress w.r.t. inverse growth deformation gradient
  AddtodPK2diFinel(dPK2diFinel,id,iFgM,gamma(0));
  AddtodPK2diFinel(dPK2diFinel,iCgC,iFgM,gamma(1));
  dPK2diFinel.MultiplyNT(delta(0),iCg,CiFg,1.);
  dPK2diFinel.MultiplyNT(delta(1),iCg,CiFgCe,1.);
  dPK2diFinel.MultiplyNT(delta(1),iCgCiCg,CiFg,1.);
  dPK2diFinel.MultiplyNT(delta(2),iCg,CiFgiCe,1.);
  dPK2diFinel.MultiplyNT(delta(2),iC,CiFg,1.);
  dPK2diFinel.MultiplyNT(delta(3),iCgCiCg,CiFgCe,1.);
  dPK2diFinel.MultiplyNT(delta(4),iCgCiCg,CiFgiCe,1.);
  dPK2diFinel.MultiplyNT(delta(4),iC,CiFgCe,1.);
  dPK2diFinel.MultiplyNT(delta(5),iC,CiFgiCe,1.);
  AddtodPK2diFinel(dPK2diFinel,id,iFgCe,0.5*delta(7));

  // diFinelastic/diFg
  LINALG::Matrix<3,3> tmp3x3(true);
  LINALG::Matrix<9,9> diFineldiFg(true);

  tmp3x3.UpdateT(1.0,GM_[gp],0.0);
  AddtoMatrixCrossProduct(1.0,id,tmp3x3,diFineldiFg);

  // diFg/dC
  LINALG::Matrix<9,6> diFgdC(true);
  LINALG::Matrix<9,1> Ag9x1(true);
  Matrix3x3to9x1(AgM,Ag9x1);

  diFgdC.MultiplyNN(-(1./(v*v))*(1./params_->density_),Ag9x1,sum_drhodC,0.0);

  // update elasticity tensor
  LINALG::Matrix<6,9> tmp6x9(true);

  tmp6x9.MultiplyNN(1.0,dPK2diFinel,diFineldiFg,0.0);
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
  double time_end = 40.0;

  if(t_tot_ > params_->t_prestress_)
  {
    for(unsigned p=0;p<potsumel_.size();++p)
      for(unsigned gp=0;gp<cur_rho_el_[p].size();++gp)
       cur_rho_el_[p][gp] = init_rho_el_[p] - 0.7*init_rho_el_[p]*(1.0-exp(-(t_tot_-params_->t_prestress_)/time_end))*exp(-0.5*(100.0*gp_ax_[gp])*(100.0*gp_ax_[gp]));
  }

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
void MAT::GrowthRemodel_ElastHyper::AddtodPK2diFinel(
    LINALG::Matrix<6,9>& dPK2diFinel,
    LINALG::Matrix<3,3> A,
    LINALG::Matrix<3,3> B,
    double fac)
{
  dPK2diFinel(0,0) += 2 * fac * A(0,0) * B(0,0);
  dPK2diFinel(0,3) += 2 * fac * A(0,0) * B(0,1);
  dPK2diFinel(0,5) += 2 * fac * A(0,0) * B(0,2);
  dPK2diFinel(0,6) += 2 * fac * A(0,1) * B(0,0);
  dPK2diFinel(0,1) += 2 * fac * A(0,1) * B(0,1);
  dPK2diFinel(0,4) += 2 * fac * A(0,1) * B(0,2);
  dPK2diFinel(0,8) += 2 * fac * A(0,2) * B(0,0);
  dPK2diFinel(0,7) += 2 * fac * A(0,2) * B(0,1);
  dPK2diFinel(0,2) += 2 * fac * A(0,2) * B(0,2);

  dPK2diFinel(1,0) += 2 * fac * A(1,0) * B(1,0);
  dPK2diFinel(1,3) += 2 * fac * A(1,0) * B(1,1);
  dPK2diFinel(1,5) += 2 * fac * A(1,0) * B(1,2);
  dPK2diFinel(1,6) += 2 * fac * A(1,1) * B(1,0);
  dPK2diFinel(1,1) += 2 * fac * A(1,1) * B(1,1);
  dPK2diFinel(1,4) += 2 * fac * A(1,1) * B(1,2);
  dPK2diFinel(1,8) += 2 * fac * A(1,2) * B(1,0);
  dPK2diFinel(1,7) += 2 * fac * A(1,2) * B(1,1);
  dPK2diFinel(1,2) += 2 * fac * A(1,2) * B(1,2);

  dPK2diFinel(2,0) += 2 * fac * A(2,0) * B(2,0);
  dPK2diFinel(2,3) += 2 * fac * A(2,0) * B(2,1);
  dPK2diFinel(2,5) += 2 * fac * A(2,0) * B(2,2);
  dPK2diFinel(2,6) += 2 * fac * A(2,1) * B(2,0);
  dPK2diFinel(2,1) += 2 * fac * A(2,1) * B(2,1);
  dPK2diFinel(2,4) += 2 * fac * A(2,1) * B(2,2);
  dPK2diFinel(2,8) += 2 * fac * A(2,2) * B(2,0);
  dPK2diFinel(2,7) += 2 * fac * A(2,2) * B(2,1);
  dPK2diFinel(2,2) += 2 * fac * A(2,2) * B(2,2);

  dPK2diFinel(3,0) += fac * (A(0,0) * B(1,0) + A(1,0) * B(0,0));
  dPK2diFinel(3,3) += fac * (A(0,0) * B(1,1) + A(1,0) * B(0,1));
  dPK2diFinel(3,5) += fac * (A(0,0) * B(1,2) + A(1,0) * B(0,2));
  dPK2diFinel(3,6) += fac * (A(0,1) * B(1,0) + A(1,1) * B(0,0));
  dPK2diFinel(3,1) += fac * (A(0,1) * B(1,1) + A(1,1) * B(0,1));
  dPK2diFinel(3,4) += fac * (A(0,1) * B(1,2) + A(1,1) * B(0,2));
  dPK2diFinel(3,8) += fac * (A(0,2) * B(1,0) + A(1,2) * B(0,0));
  dPK2diFinel(3,7) += fac * (A(0,2) * B(1,1) + A(1,2) * B(0,1));
  dPK2diFinel(3,2) += fac * (A(0,2) * B(1,2) + A(1,2) * B(0,2));

  dPK2diFinel(4,0) += fac * (A(1,0) * B(2,0) + A(2,0) * B(1,0));
  dPK2diFinel(4,3) += fac * (A(1,0) * B(2,1) + A(2,0) * B(1,1));
  dPK2diFinel(4,5) += fac * (A(1,0) * B(2,2) + A(2,0) * B(1,2));
  dPK2diFinel(4,6) += fac * (A(1,1) * B(2,0) + A(2,1) * B(1,0));
  dPK2diFinel(4,1) += fac * (A(1,1) * B(2,1) + A(2,1) * B(1,1));
  dPK2diFinel(4,4) += fac * (A(1,1) * B(2,2) + A(2,1) * B(1,2));
  dPK2diFinel(4,8) += fac * (A(1,2) * B(2,0) + A(2,2) * B(1,0));
  dPK2diFinel(4,7) += fac * (A(1,2) * B(2,1) + A(2,2) * B(1,1));
  dPK2diFinel(4,2) += fac * (A(1,2) * B(2,2) + A(2,2) * B(1,2));

  dPK2diFinel(5,0) += fac * (A(0,0) * B(2,0) + A(2,0) * B(0,0));
  dPK2diFinel(5,3) += fac * (A(0,0) * B(2,1) + A(2,0) * B(0,1));
  dPK2diFinel(5,5) += fac * (A(0,0) * B(2,2) + A(2,0) * B(0,2));
  dPK2diFinel(5,6) += fac * (A(0,1) * B(2,0) + A(2,1) * B(0,0));
  dPK2diFinel(5,1) += fac * (A(0,1) * B(2,1) + A(2,1) * B(0,1));
  dPK2diFinel(5,4) += fac * (A(0,1) * B(2,2) + A(2,1) * B(0,2));
  dPK2diFinel(5,8) += fac * (A(0,2) * B(2,0) + A(2,2) * B(0,0));
  dPK2diFinel(5,7) += fac * (A(0,2) * B(2,1) + A(2,2) * B(0,1));
  dPK2diFinel(5,2) += fac * (A(0,2) * B(2,2) + A(2,2) * B(0,2));

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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::VisNames(std::map<std::string,int>& names)
{
  std::string mass_fraction_el = "mass_fraction_el";
  std::string result_mass_fraction_el;

  for (unsigned int p=0; p<potsumel_.size(); ++p)
  {
    std::stringstream sstm;
    sstm << mass_fraction_el <<"_" << p;
    result_mass_fraction_el = sstm.str();

    names[result_mass_fraction_el] = 1;
  }


  std::string v_growth = "v_growth";
  std::string result_v_growth;

  std::stringstream sstm;
  sstm << v_growth;
  result_v_growth = sstm.str();

  names[result_v_growth] = 1;


  // loop map of associated potential summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->VisNames(names,p);

  // elastin matrix
  for (unsigned int p=0; p<potsumel_.size(); ++p)
    potsumel_[p]->VisNames(names);

  // volpenalty
  if(params_->matid_penalty_ != -1)
    potsumpenalty_->VisNames(names);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::GrowthRemodel_ElastHyper::VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  int return_val= 0;

  if(name == "mass_fraction_el_0")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned i=0;i<cur_rho_el_[0].size();++i)
    {
      data[0] += cur_rho_el_[0][i];
    }
    data[0] = data[0]/cur_rho_el_[0].size();

    return true;
  }
  if(name == "mass_fraction_el_1")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned i=0;i<cur_rho_el_[1].size();++i)
    {
      data[0] += cur_rho_el_[1][i];
    }
    data[0] = data[0]/cur_rho_el_[1].size();

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

  // elastin matrix
  for (unsigned int p=0; p<potsumel_.size(); ++p)
    return_val+=potsumel_[p]->VisData(name,data,numgp,eleID);

  // volpenalty
  if(params_->matid_penalty_ != -1)
    return_val+=potsumpenalty_->VisData(name,data,numgp,eleID);

  return (bool)return_val;
}
