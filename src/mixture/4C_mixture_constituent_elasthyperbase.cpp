/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the hyperelastic constituent

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_mixture_constituent_elasthyperbase.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_mixture.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper_service.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_mixture_prestress_strategy.hpp"

#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

// Constructor for the parameter class
MIXTURE::PAR::MixtureConstituentElastHyperBase::MixtureConstituentElastHyperBase(
    const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : MixtureConstituent(matdata),
      matid_prestress_strategy_(matdata->Get<int>("PRESTRESS_STRATEGY")),
      nummat_(matdata->Get<int>("NUMMAT")),
      matids_(matdata->Get<std::vector<int>>("MATIDS"))
{
  // check, if size of summands fits to the number of summands
  if (nummat_ != (int)matids_.size())
  {
    FOUR_C_THROW(
        "number of summands %d does not fit to the size of the summands vector"
        " %d",
        nummat_, matids_.size());
  }
}

// Constructor of the constituent holding the material parameters
MIXTURE::MixtureConstituentElastHyperBase::MixtureConstituentElastHyperBase(
    MIXTURE::PAR::MixtureConstituentElastHyperBase* params, int id)
    : MixtureConstituent(params, id),
      summand_properties_(),
      params_(params),
      potsum_(0),
      cosy_anisotropy_extension_()
{
  // Create summands
  for (const auto& matid : params_->matids_)
  {
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) FOUR_C_THROW("Failed to read elastic summand.");
    potsum_.push_back(sum);
  }

  // Create Prestress strategy
  if (params->get_prestressing_mat_id() > 0)
  {
    prestress_strategy_ =
        MIXTURE::PAR::PrestressStrategy::Factory(params->get_prestressing_mat_id())
            ->create_prestress_strategy();
  }
}

// Pack the constituent
void MIXTURE::MixtureConstituentElastHyperBase::PackConstituent(CORE::COMM::PackBuffer& data) const
{
  MixtureConstituent::PackConstituent(data);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  CORE::COMM::ParObject::AddtoPack(data, matid);
  summand_properties_.Pack(data);

  CORE::COMM::ParObject::AddtoPack(data, prestretch_);

  cosy_anisotropy_extension_.PackAnisotropy(data);

  if (prestress_strategy_ != nullptr) prestress_strategy_->Pack(data);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (const auto& p : potsum_) p->PackSummand(data);
  }
}

// Unpack the constituent
void MIXTURE::MixtureConstituentElastHyperBase::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureConstituent::UnpackConstituent(position, data);

  // make sure we have a pristine material
  params_ = nullptr;
  potsum_.clear();

  // matid and recover params_
  int matid;
  CORE::COMM::ParObject::ExtractfromPack(position, data, matid);

  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const unsigned int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        params_ = dynamic_cast<MIXTURE::PAR::MixtureConstituentElastHyperBase*>(mat);
      }
      else
      {
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
      }
    }
  }

  summand_properties_.Unpack(position, data);

  CORE::COMM::ParObject::ExtractfromPack(position, data, prestretch_);

  cosy_anisotropy_extension_.UnpackAnisotropy(data, position);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    if (params_->get_prestressing_mat_id() > 0)
    {
      prestress_strategy_ =
          MIXTURE::PAR::PrestressStrategy::Factory(params_->get_prestressing_mat_id())
              ->create_prestress_strategy();

      prestress_strategy_->Unpack(position, data);
    }

    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = params_->matids_.begin(); m != params_->matids_.end(); ++m)
    {
      const int summatid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(summatid);
      if (sum == Teuchos::null) FOUR_C_THROW("Failed to allocate");
      potsum_.push_back(sum);
    }

    // loop map of associated potential summands
    for (auto& summand : potsum_) summand->UnpackSummand(data, position);
  }
}

void MIXTURE::MixtureConstituentElastHyperBase::register_anisotropy_extensions(
    MAT::Anisotropy& anisotropy)
{
  // Setup summands
  for (const auto& summand : potsum_) summand->register_anisotropy_extensions(anisotropy);

  anisotropy.register_anisotropy_extension(cosy_anisotropy_extension_);
}

// Reads the element from the input file
void MIXTURE::MixtureConstituentElastHyperBase::ReadElement(
    int numgp, INPUT::LineDefinition* linedef)
{
  MixtureConstituent::ReadElement(numgp, linedef);

  // Setup summands
  for (const auto& summand : potsum_) summand->Setup(numgp, linedef);

  // find out which formulations are used
  MAT::ElastHyperProperties(potsum_, summand_properties_);

  if (summand_properties_.viscoGeneral)
  {
    FOUR_C_THROW("Never use viscoelastic materials in Elasthyper-Toolbox.");
  }
}

// Updates all summands
void MIXTURE::MixtureConstituentElastHyperBase::Update(CORE::LINALG::Matrix<3, 3> const& defgrd,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  MixtureConstituent::Update(defgrd, params, gp, eleGID);

  // loop map of associated potential summands
  for (auto& summand : potsum_) summand->Update();

  // do nothing in the default case
  if (params_->get_prestressing_mat_id() > 0)
  {
    prestress_strategy_->Update(cosy_anisotropy_extension_.get_coordinate_system_provider(gp),
        *this, defgrd, prestretch_[gp], params, gp, eleGID);
  }
}

void MIXTURE::MixtureConstituentElastHyperBase::Setup(
    Teuchos::ParameterList& params, const int eleGID)
{
  MixtureConstituent::Setup(params, eleGID);
  if (params_->get_prestressing_mat_id() > 0)
  {
    prestretch_.resize(num_gp());

    prestress_strategy_->Setup(*this, params, num_gp(), eleGID);
  }
}

void MIXTURE::MixtureConstituentElastHyperBase::pre_evaluate(
    MixtureRule& mixtureRule, Teuchos::ParameterList& params, int gp, int eleGID)
{
  // do nothing in the default case
  if (params_->get_prestressing_mat_id() > 0)
  {
    prestress_strategy_->EvaluatePrestress(mixtureRule,
        cosy_anisotropy_extension_.get_coordinate_system_provider(gp), *this, prestretch_[gp],
        params, gp, eleGID);
  }
}

void MIXTURE::MixtureConstituentElastHyperBase::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  if (prestress_strategy_ != nullptr)
  {
    names_and_size["mixture_constituent_" + std::to_string(Id()) + "_elasthyper_prestretch"] = 9;
  }
}

bool MIXTURE::MixtureConstituentElastHyperBase::EvaluateOutputData(
    const std::string& name, CORE::LINALG::SerialDenseMatrix& data) const
{
  if (prestress_strategy_ != nullptr &&
      name == "mixture_constituent_" + std::to_string(Id()) + "_elasthyper_prestretch")
  {
    for (int gp = 0; gp < num_gp(); ++gp)
    {
      static CORE::LINALG::Matrix<9, 1> tmp(false);
      tmp.Clear();
      CORE::LINALG::VOIGT::Matrix3x3to9x1(prestretch_[gp], tmp);

      for (int i = 0; i < 9; ++i)
      {
        data(gp, i) = tmp(i, 0);
      }
    }
    return true;
  }
  return false;
}
FOUR_C_NAMESPACE_CLOSE
