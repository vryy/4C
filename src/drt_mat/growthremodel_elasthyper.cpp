/*----------------------------------------------------------------------*/
/*!
\file growthremodel_elasthyper.cpp
\brief
This file is used to manage the constraint mixture in growth and remodeling

The input line should read
MAT 0   MAT_Growthremodel_ElastHyper   NUMMATRF 0 MATIDSRF NUMMATEL 0 MATIDSEL NUMMATGR 0 MATIDSGR MATIDPENALTY ELMASSFRAC GRMASSFRAC DENS 0

<pre>
Maintainer: Fabian Bräu
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
#include "../drt_matelast/elast_volgrowthpenalty.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::GrowthRemodel_ElastHyper::GrowthRemodel_ElastHyper(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  nummat_remodelfiber_(matdata->GetInt("NUMMATRF")),
  nummat_elastin_(matdata->GetInt("NUMMATEL")),
  nummat_ground_(matdata->GetInt("NUMMATGR")),
  matids_remodelfiber_(matdata->Get<std::vector<int> >("MATIDSRF")),
  matids_elastin_(matdata->Get<std::vector<int> >("MATIDSEL")),
  matids_ground_(matdata->Get<std::vector<int> >("MATIDSGR")),
  matid_penalty_(matdata->GetInt("MATIDPENALTY")),
  cur_w_elastin_(matdata->Get<std::vector<double> >("ELMASSFRAC")),
  cur_w_ground_(matdata->Get<std::vector<double> >("GRMASSFRAC")),
  density_(matdata->GetDouble("DENS"))
{
  // check if sizes fit
  if (nummat_remodelfiber_ != (int)matids_remodelfiber_->size())
    dserror("number of remodelfiber materials %d does not fit to size of remodelfiber material vector %d", nummat_remodelfiber_, matids_remodelfiber_->size());

  if (nummat_elastin_ != (int)matids_elastin_->size())
    dserror("number of elastin materials %d does not fit to size of elastin material vector %d", nummat_elastin_, matids_elastin_->size());

  if (nummat_ground_ != (int)matids_ground_->size())
    dserror("number of ground materials %d does not fit to size of ground material vector %d", nummat_ground_, matids_ground_->size());

  if(matid_penalty_ == -1)
    dserror("growth and remodeling is only working with a volumetric growth penalty material (in the current implementation)");
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


/*----------------------------------------------------------------------*
 |  initialise static arrays                                 seitz 05/14|
 *----------------------------------------------------------------------*/
const int MAT::GrowthRemodel_ElastHyper::VOIGT3X3_[3][3]       = {{0,3,5},{3,1,4},{5,4,2}};
const int MAT::GrowthRemodel_ElastHyper::VOIGT3X3NONSYM_[3][3] = {{0,3,5},{6,1,4},{8,7,2}};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::GrowthRemodel_ElastHyper::GrowthRemodel_ElastHyper()
  : params_(NULL),
    potsumrf_(0),
    potsumel_(0),
    potsumgr_(0),
    potsumpenalty_(0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::GrowthRemodel_ElastHyper::GrowthRemodel_ElastHyper(MAT::PAR::GrowthRemodel_ElastHyper* params)
  : params_(params),
    potsumrf_(0),
    potsumel_(0),
    potsumgr_(0),
    potsumpenalty_(0)
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

  // Ground matrix
  for (m=params_->matids_ground_->begin(); m!=params_->matids_ground_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumgr_.push_back(sum);
  }

  // VolGrowthPenalty
  Teuchos::RCP<MAT::ELASTIC::VolGrowthPenalty> sum = Teuchos::rcp_static_cast<MAT::ELASTIC::VolGrowthPenalty>(MAT::ELASTIC::Summand::Factory(params_->matid_penalty_));
  if (sum == Teuchos::null) dserror("Failed to allocate");
  potsumpenalty_ = sum;
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
  int num_gr = 0;
  num_el = current_w_elastin_.size();
  num_gr = current_w_ground_.size();

  AddtoPack(data,num_el);
  AddtoPack(data,num_gr);

  for(int i=0;i<num_el;++i)
    AddtoPack(data,current_w_elastin_[i]);

  for(int i=0;i<num_gr;++i)
    AddtoPack(data,current_w_ground_[i]);


  if (params_ != NULL) // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumrf_.size(); ++p)
     potsumrf_[p]->PackSummand(data);

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumel_.size(); ++p)
     potsumel_[p]->PackSummand(data);

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumgr_.size(); ++p)
     potsumgr_[p]->PackSummand(data);

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
  potsumgr_.clear();

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
  int num_gr = 0;

  ExtractfromPack(position,data,num_el);
  current_w_elastin_.resize(num_el);

  ExtractfromPack(position,data,num_gr);
  current_w_ground_.resize(num_gr);

  for(int i=0;i<num_el;++i)
    ExtractfromPack(position,data,current_w_elastin_[i]);

  for(int i=0;i<num_gr;++i)
    ExtractfromPack(position,data,current_w_ground_[i]);


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

    // Ground matrix
    for (m=params_->matids_ground_->begin(); m!=params_->matids_ground_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsumgr_.push_back(sum);
    }

    // VolGrowthPenalty
    Teuchos::RCP<MAT::ELASTIC::VolGrowthPenalty> sum = Teuchos::rcp_static_cast<MAT::ELASTIC::VolGrowthPenalty>(MAT::ELASTIC::Summand::Factory(params_->matid_penalty_));
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumpenalty_ = sum;

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumrf_.size(); ++p)
     potsumrf_[p]->UnpackSummand(data,position);

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumel_.size(); ++p)
     potsumel_[p]->UnpackSummand(data,position);

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumgr_.size(); ++p)
     potsumgr_[p]->UnpackSummand(data,position);

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
  current_w_elastin_.resize(potsumel_.size());
  current_w_ground_.resize(potsumgr_.size());

  for(unsigned p=0;p<potsumel_.size();++p)
    current_w_elastin_[p].resize(numgp,1.0);

  for(unsigned p=0;p<potsumgr_.size();++p)
    current_w_ground_[p].resize(numgp,1.0);


  // Setup summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->Setup(numgp,linedef);

  // elastin matrix
  for (unsigned int p=0; p<potsumel_.size(); ++p)
    potsumel_[p]->Setup(linedef);

  // ground matrix
  for (unsigned int p=0; p<potsumgr_.size(); ++p)
    potsumgr_[p]->Setup(linedef);

  // volgrowthpenalty
  potsumpenalty_->Setup(linedef);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::Update()
{
  // loop map of associated potential summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->Update();

  // elastin matrix
  for (unsigned int p=0; p<potsumel_.size(); ++p)
    potsumel_[p]->Update();

  // ground matrix
  for (unsigned int p=0; p<potsumgr_.size(); ++p)
    potsumgr_[p]->Update();

  // volgrowthpenalty
  potsumpenalty_->Update();

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
  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress->Clear();
  cmat->Clear();

  // Evaluate elasin and ground matrix and volumetric growth penalty
  LINALG::Matrix<6,1> id2(true) ;
  LINALG::Matrix<6,1> C_strain(true) ;
  LINALG::Matrix<6,1> C_voigt(true) ;
  LINALG::Matrix<6,1> iC_voigt(true) ;
  LINALG::Matrix<6,6> id4(true) ;
  LINALG::Matrix<6,6> id4sharp(true) ;

  LINALG::Matrix<3,1> prinv(true);
  LINALG::Matrix<3,1> dPI(true);
  LINALG::Matrix<6,1> ddPII(true);

  EvaluateKinQuant(*glstrain,id2,C_voigt,C_strain,iC_voigt,id4,id4sharp,prinv);
  EvaluateInvariantDerivatives(prinv,dPI,ddPII,params,eleGID);

  // build stress response and elasticity tensor
  LINALG::Matrix<NUM_STRESS_3D,1> stressiso(true) ;
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true) ;

  EvaluateIsotropicStressCmat(stressiso,cmatiso,C_voigt,id2,iC_voigt,id4sharp,prinv,dPI,ddPII);

  stress->Update(1.0,stressiso,1.0);
  cmat->Update(1.0,cmatiso,1.0);


  // evaluate remodelfiber (remodeling of collagen)
  // build stress response and elasticity tensor
  LINALG::Matrix<NUM_STRESS_3D,1> stressremodel(true) ;
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatremodel(true) ;

  for(unsigned p=0;p<potsumrf_.size();++p)
  {
    potsumrf_[p]->AddStressCmatRemodel(defgrd,params,cmatremodel,stressremodel,eleGID);
    stress->Update(1.0,stressremodel,1.0);
    cmat->Update(1.0,cmatremodel,1.0);
  }

  return ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::GrowthRemodel_ElastHyper::EvaluateKinQuant(
    const LINALG::Matrix<6,1>& glstrain,
    LINALG::Matrix<6,1>& id2,
    LINALG::Matrix<6,1>& C_voigt,
    LINALG::Matrix<6,1>& C_strain,
    LINALG::Matrix<6,1>& iC_voigt,
    LINALG::Matrix<6,6>& id4,
    LINALG::Matrix<6,6>& id4sharp,
    LINALG::Matrix<3,1>& prinv)

{
  // build Cartesian identity 2-tensor I_{AB}
  for (int i=0; i<3; i++) id2(i) = 1.0;

  // right Cauchy-Green Tensor  C_{AB} = 2 * E_{AB} + I_{AB}
  // REMARK: strain-like 6-Voigt vector
  C_strain.Update(2.0,glstrain,1.0);
  C_strain.Update(1.0, id2, 1.0);

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

  // set Cartesian identity 4-tensor in 6x6-matrix notation (stress-like)
  // this is a 'mixed co- and contra-variant' identity 4-tensor, ie I^{AB}_{CD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are strain-like 6-Voigt
  for (int i=0; i<6; i++) id4(i,i) = 1.0;


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
    LINALG::Matrix<3,1>& dPI,
    LINALG::Matrix<6,1>& ddPII,
    Teuchos::ParameterList& params,
    int eleGID
)

{
  // current Gauß-Point
  int gp = params.get<int>("gp");

  // derivatives of principal materials weighted with their mass fraction in the constraint mixture
  LINALG::Matrix<3,1> dPgrowthI(true);
  LINALG::Matrix<6,1> ddPgrowthII(true);

  // loop map of associated potential summands
  // elastin matrix
  for(unsigned p=0;p<potsumel_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumel_[p]->AddDerivativesPrincipal(dPgrowthI,ddPgrowthII,prinv,eleGID);
    dPI.Update(current_w_elastin_[p][gp],dPgrowthI,1.0);
    ddPII.Update(current_w_elastin_[p][gp],ddPgrowthII,1.0);
  }
  //ground matrix
  for(unsigned p=0;p<potsumgr_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumgr_[p]->AddDerivativesPrincipal(dPgrowthI,ddPgrowthII,prinv,eleGID);
    dPI.Update(current_w_ground_[p][gp],dPgrowthI,1.0);
    ddPII.Update(current_w_ground_[p][gp],ddPgrowthII,1.0);
  }



  // derivatives of decoupled (volumetric or isochoric) materials weighted with their mass fraction in the constraint mixture
  LINALG::Matrix<3,1> modinv(true);
  InvariantsModified(modinv,prinv);
  LINALG::Matrix<3,1> dPmodI(true);
  LINALG::Matrix<6,1> ddPmodII(true);

  // loop map of associated potential summands
  // elastin matrix
  for(unsigned p=0;p<potsumel_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumel_[p]->AddDerivativesModified(dPgrowthI,ddPgrowthII,modinv,eleGID);
    dPmodI.Update(current_w_elastin_[p][gp],dPgrowthI,1.0);
    ddPmodII.Update(current_w_elastin_[p][gp],ddPgrowthII,1.0);
  }
  // ground matrix
  for(unsigned p=0;p<potsumgr_.size();++p)
  {
    dPgrowthI.Clear();
    ddPgrowthII.Clear();

    potsumgr_[p]->AddDerivativesModified(dPgrowthI,ddPgrowthII,prinv,eleGID);
    dPmodI.Update(current_w_ground_[p][gp],dPgrowthI,1.0);
    ddPmodII.Update(current_w_ground_[p][gp],ddPgrowthII,1.0);
  }
  // volgrowthpenalty
  potsumpenalty_->AddDerivativesModified(dPmodI,ddPmodII,modinv,eleGID);

  // convert decoupled derivatives to principal derivatives
  ConvertModToPrinc(prinv,dPmodI,ddPmodII,dPI,ddPII);


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
    LINALG::Matrix<6,1> id2,
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
  stress.Update(gamma(0), id2, 1.0);
  stress.Update(gamma(1), scg, 1.0);
  stress.Update(gamma(2), icg, 1.0);

  // constitutive tensor
  // contribution: Id \otimes Id
  cmat.MultiplyNT(delta(0), id2, id2, 1.0);
  // contribution: Id \otimes C + C \otimes Id
  cmat.MultiplyNT(delta(1), id2, scg, 1.0);
  cmat.MultiplyNT(delta(1), scg, id2, 1.0);
  // contribution: Id \otimes Cinv + Cinv \otimes Id
  cmat.MultiplyNT(delta(2), id2, icg, 1.0);
  cmat.MultiplyNT(delta(2), icg, id2, 1.0);
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
void MAT::GrowthRemodel_ElastHyper::VisNames(std::map<std::string,int>& names)
{
  // loop map of associated potential summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->VisNames(names);

  // elastin matrix
  for (unsigned int p=0; p<potsumel_.size(); ++p)
    potsumel_[p]->VisNames(names);

  // ground matrix
  for (unsigned int p=0; p<potsumgr_.size(); ++p)
    potsumgr_[p]->VisNames(names);

  // volgrowthpenalty
  potsumpenalty_->VisNames(names);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::GrowthRemodel_ElastHyper::VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  int return_val= 0;

  // loop map of associated potential summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    return_val+=potsumrf_[p]->VisData(name,data,numgp,eleID);

  // elastin matrix
  for (unsigned int p=0; p<potsumel_.size(); ++p)
    return_val+=potsumel_[p]->VisData(name,data,numgp,eleID);

  // ground matrix
  for (unsigned int p=0; p<potsumgr_.size(); ++p)
    return_val+=potsumgr_[p]->VisData(name,data,numgp,eleID);

  // volgrowthpenalty
  return_val+=potsumpenalty_->VisData(name,data,numgp,eleID);

  return (bool)return_val;
}
