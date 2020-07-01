/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the hyperelastic constituent

\level 3


*/
/*----------------------------------------------------------------------*/

#include "mixture_constituent_elasthyperbase.H"
#include "../drt_mat/material_service.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/multiplicative_split_defgrad_elasthyper_service.H"
#include "../drt_mat/mixture_elasthyper.H"
#include "mixture_prestress_strategy.H"

// Constructor for the parameter class
MIXTURE::PAR::MixtureConstituent_ElastHyperBase::MixtureConstituent_ElastHyperBase(
    const Teuchos::RCP<MAT::PAR::Material>& matdata, const double ref_mass_fraction)
    : MixtureConstituent(matdata, ref_mass_fraction),
      matid_prestress_strategy_(matdata->GetInt("PRESTRESS_STRATEGY")),
      nummat_(matdata->GetInt("NUMMAT")),
      matids_(matdata->Get<std::vector<int>>("MATIDS"))
{
  // check, if size of summands fits to the number of summands
  if (nummat_ != (int)matids_->size())
  {
    dserror(
        "number of summands %d does not fit to the size of the summands vector"
        " %d",
        nummat_, matids_->size());
  }


  // Create Prestress strategy
  if (GetPrestressingMatId() > 0)
  {
    prestressStrategy_ =
        MIXTURE::PAR::PrestressStrategy::Factory(GetPrestressingMatId())->CreatePrestressStrategy();
  }
}

// Constructor of the constituent holding the material parameters
MIXTURE::MixtureConstituent_ElastHyperBase::MixtureConstituent_ElastHyperBase(
    MIXTURE::PAR::MixtureConstituent_ElastHyperBase* params)
    : MixtureConstituent(params),
      summandProperties_(),
      params_(params),
      potsum_(0),
      cosyAnisotropyExtension_()
{
  // Create summands
  for (const auto& matid : *params_->matids_)
  {
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to read elastic summand.");
    potsum_.push_back(sum);
  }
}

// Pack the constituent
void MIXTURE::MixtureConstituent_ElastHyperBase::PackConstituent(DRT::PackBuffer& data) const
{
  MixtureConstituent::PackConstituent(data);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  DRT::ParObject::AddtoPack(data, matid);
  summandProperties_.Pack(data);

  DRT::ParObject::AddtoPack(data, prestretch_);

  cosyAnisotropyExtension_.PackAnisotropy(data);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (const auto& p : potsum_) p->PackSummand(data);
  }
}

// Unpack the constituent
void MIXTURE::MixtureConstituent_ElastHyperBase::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureConstituent::UnpackConstituent(position, data);

  // make sure we have a pristine material
  params_ = nullptr;
  potsum_.clear();

  // matid and recover params_
  int matid;
  DRT::ParObject::ExtractfromPack(position, data, matid);

  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const unsigned int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        params_ = dynamic_cast<MIXTURE::PAR::MixtureConstituent_ElastHyperBase*>(mat);
      }
      else
      {
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
      }
    }
  }

  summandProperties_.Unpack(position, data);

  DRT::ParObject::ExtractfromPack(position, data, prestretch_);

  cosyAnisotropyExtension_.UnpackAnisotropy(data, position);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = params_->matids_->begin(); m != params_->matids_->end(); ++m)
    {
      const int summatid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(summatid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsum_.push_back(sum);
    }


    /*if (params_->GetPrestressingMatId() > 0)
    {
      prestressStrategy_ = MIXTURE::PAR::PrestressStrategy::Factory(params_->GetPrestressingMatId())
                               ->CreatePrestressStrategy();
    }*/

    // loop map of associated potential summands
    for (auto& summand : potsum_) summand->UnpackSummand(data, position);
  }
}

void MIXTURE::MixtureConstituent_ElastHyperBase::RegisterAnisotropyExtensions(
    MAT::Anisotropy& anisotropy)
{
  // Setup summands
  for (const auto& summand : potsum_) summand->RegisterAnisotropyExtensions(anisotropy);

  anisotropy.RegisterAnisotropyExtension(cosyAnisotropyExtension_);
}

// Reads the element from the input file
void MIXTURE::MixtureConstituent_ElastHyperBase::ReadElement(
    int numgp, DRT::INPUT::LineDefinition* linedef)
{
  MixtureConstituent::ReadElement(numgp, linedef);

  // Setup summands
  for (const auto& summand : potsum_) summand->Setup(numgp, linedef);

  // find out which formulations are used
  MAT::ElastHyperProperties(potsum_, summandProperties_);

  if (summandProperties_.viscoGeneral)
  {
    dserror("Never use viscoelastic materials in Elasthyper-Toolbox.");
  }
}

// Returns the reference mass fraction of the constituent
double MIXTURE::MixtureConstituent_ElastHyperBase::CurrentRefDensity(int gp) const
{
  return params_->RefMassFraction() * InitialRefDensity();
}

// Updates all summands
void MIXTURE::MixtureConstituent_ElastHyperBase::Update(LINALG::Matrix<3, 3> const& defgrd,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  MixtureConstituent::Update(defgrd, params, gp, eleGID);

  // loop map of associated potential summands
  for (auto& summand : potsum_) summand->Update();
}

// Add names for each summand for the quantities for post processing
void MIXTURE::MixtureConstituent_ElastHyperBase::VisNames(std::map<std::string, int>& names)
{
  MixtureConstituent::VisNames(names);

  // loop map of associated potential summands
  for (auto& summand : potsum_) summand->VisNames(names);
}

// Add values for each summand of the quantities for post processing
bool MIXTURE::MixtureConstituent_ElastHyperBase::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleGID)
{
  // loop map of associated potential summands
  for (auto& summand : potsum_)
  {
    bool vis = summand->VisData(name, data, numgp, eleGID);
    if (vis) return true;
  }
  return MixtureConstituent::VisData(name, data, numgp, eleGID);
}

void MIXTURE::MixtureConstituent_ElastHyperBase::Setup(
    Teuchos::ParameterList& params, const int eleGID)
{
  MixtureConstituent::Setup(params, eleGID);
  if (params_->GetPrestressingMatId() > 0)
  {
    prestretch_.resize(NumGP());
  }
}

void MIXTURE::MixtureConstituent_ElastHyperBase::PreEvaluate(
    MixtureRule& mixtureRule, Teuchos::ParameterList& params, int gp, int eleGID)
{
  // do nothing in the default case
  if (params_->GetPrestressingMatId() > 0)
  {
    LINALG::Matrix<3, 3> prestretch;
    params_->PrestressStrategy()->EvaluatePrestress(
        cosyAnisotropyExtension_.GetCylinderCoordinateSystem(gp), *this, prestretch_[gp], params,
        gp, eleGID);
  }
}
