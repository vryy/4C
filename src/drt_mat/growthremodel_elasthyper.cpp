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
  // Setup summands
  // remodelfiber
  for (unsigned int p=0; p<potsumrf_.size(); ++p)
    potsumrf_[p]->Setup(linedef);

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
  dserror("Not implemented so far!");

  return ;
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
