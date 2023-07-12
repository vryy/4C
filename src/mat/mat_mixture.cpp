/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the mixture material holding a general mixturerule and mixture constituents

\level 3

*/
/*----------------------------------------------------------------------*/
#include <memory>

#include "mat_mixture.H"
#include "mat_service.H"
#include "lib_globalproblem.H"
#include "mat_par_bundle.H"

// constructor of the parameters
MAT::PAR::Mixture::Mixture(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), constituents_(0)
{
  const int num_constituents = matdata->GetInt("NUMCONST");
  const auto* constituent_matids = matdata->Get<std::vector<int>>("MATIDSCONST");

  // check, if size of constituents fits to the number of constituents
  if (num_constituents != (int)constituent_matids->size())
  {
    dserror(
        "number of constituents %d does not fit to the size of the constituents material vector"
        " %d",
        num_constituents, constituent_matids->size());
  }

  // Create constituents
  for (int i = 0; i < num_constituents; ++i)
  {
    // Create constituent material
    constituents_.emplace_back(MIXTURE::PAR::MixtureConstituent::Factory((*constituent_matids)[i]));
  }

  // Create mixture rule
  mixture_rule_ = MIXTURE::PAR::MixtureRule::Factory(matdata->GetInt("MATIDMIXTURERULE"));
}

// Create a material instance from parameters
Teuchos::RCP<MAT::Material> MAT::PAR::Mixture::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Mixture(this));
}

MAT::MixtureType MAT::MixtureType::instance_;

// Create a material instance from packed data
DRT::ParObject* MAT::MixtureType::Create(const std::vector<char>& data)
{
  auto* mix_elhy = new MAT::Mixture();
  mix_elhy->Unpack(data);

  return mix_elhy;
}

// constructor
MAT::Mixture::Mixture()
    : params_(nullptr),
      constituents_(std::make_shared<std::vector<std::unique_ptr<MIXTURE::MixtureConstituent>>>(0)),
      setup_(false),
      anisotropy_()
{
}

// constructor
MAT::Mixture::Mixture(MAT::PAR::Mixture* params)
    : params_(params),
      constituents_(std::make_shared<std::vector<std::unique_ptr<MIXTURE::MixtureConstituent>>>(0)),
      setup_(false),
      anisotropy_()
{
  // create instances of constituents
  int id = 0;
  for (auto const& constituent : params_->constituents_)
  {
    constituents_->emplace_back(constituent->CreateConstituent(id));
    constituents_->back()->RegisterAnisotropyExtensions(anisotropy_);

    ++id;
  }

  // create instance of mixture rule
  mixture_rule_ = params->mixture_rule_->CreateRule();
  mixture_rule_->SetConstituents(constituents_);
  mixture_rule_->RegisterAnisotropyExtensions(anisotropy_);
}

// Pack data
void MAT::Mixture::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // Pack material id
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // pack setup flag
  AddtoPack(data, static_cast<int>(setup_));

  // Pack isPreEvaluated flag
  std::vector<int> isPreEvaluatedInt;
  isPreEvaluatedInt.resize(isPreEvaluated_.size());
  for (unsigned i = 0; i < isPreEvaluated_.size(); ++i)
  {
    isPreEvaluatedInt[i] = static_cast<int>(isPreEvaluated_[i]);
  }
  AddtoPack(data, isPreEvaluatedInt);

  anisotropy_.PackAnisotropy(data);

  // pack all constituents
  // constituents are not accessible during post processing
  if (params_ != nullptr)
  {
    for (const auto& constituent : *constituents_)
    {
      constituent->PackConstituent(data);
    }

    // pack mixturerule
    mixture_rule_->PackMixtureRule(data);
  }
}

// Unpack data
void MAT::Mixture::Unpack(const std::vector<char>& data)
{
  params_ = nullptr;
  constituents_->clear();
  setup_ = false;

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        params_ = dynamic_cast<MAT::PAR::Mixture*>(mat);
      }
      else
      {
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
      }
    }

    // Extract setup flag
    setup_ = (bool)ExtractInt(position, data);


    // Extract is isPreEvaluated
    std::vector<int> isPreEvaluatedInt(0);
    DRT::ParObject::ExtractfromPack(position, data, isPreEvaluatedInt);
    isPreEvaluated_.resize(isPreEvaluatedInt.size());
    for (unsigned i = 0; i < isPreEvaluatedInt.size(); ++i)
    {
      isPreEvaluated_[i] = static_cast<bool>(isPreEvaluatedInt[i]);
    }

    anisotropy_.UnpackAnisotropy(data, position);

    // extract constituents
    // constituents are not accessible during post processing
    if (params_ != nullptr)
    {
      // create instances of constituents
      int id = 0;
      for (auto const& constituent : params_->constituents_)
      {
        constituents_->emplace_back(constituent->CreateConstituent(id));

        ++id;
      }

      // create instance of mixture rule
      mixture_rule_ = params_->mixture_rule_->CreateRule();

      // make sure the referenced materials in material list have quick access parameters
      for (const auto& constituent : *constituents_)
      {
        constituent->UnpackConstituent(position, data);
        constituent->RegisterAnisotropyExtensions(anisotropy_);
      }

      // unpack mixturerule
      mixture_rule_->UnpackMixtureRule(position, data);
      mixture_rule_->SetConstituents(constituents_);
      mixture_rule_->RegisterAnisotropyExtensions(anisotropy_);

      // position checking is not available in post processing mode
      if (position != data.size())
      {
        dserror("Mismatch in size of data to unpack (%d <-> %d)", data.size(), position);
      }
    }
  }
}

// Read element and create arrays for the quantities at the Gauss points
void MAT::Mixture::Setup(const int numgp, DRT::INPUT::LineDefinition* linedef)
{
  So3Material::Setup(numgp, linedef);

  // resize preevaluation flag
  isPreEvaluated_.resize(numgp, false);

  // Setup anisotropy
  anisotropy_.SetNumberOfGaussPoints(numgp);
  anisotropy_.ReadAnisotropyFromElement(linedef);

  // Let all constituents read the line definition
  for (const auto& constituent : *constituents_)
  {
    constituent->ReadElement(numgp, linedef);
  }

  mixture_rule_->ReadElement(numgp, linedef);
}

// Post setup routine -> Call Setup of constituents and mixture rule
void MAT::Mixture::PostSetup(Teuchos::ParameterList& params, const int eleGID)
{
  So3Material::PostSetup(params, eleGID);
  anisotropy_.ReadAnisotropyFromParameterList(params);
  if (constituents_ != nullptr)
  {
    for (const auto& constituent : *constituents_)
    {
      constituent->Setup(params, eleGID);
    }
  }

  if (mixture_rule_ != nullptr)
  {
    mixture_rule_->Setup(params, eleGID);
  }

  setup_ = true;
}

// This method is called between two timesteps
void MAT::Mixture::Update(CORE::LINALG::Matrix<3, 3> const& defgrd, const int gp,
    Teuchos::ParameterList& params, const int eleGID)
{
  // Update all constituents
  for (const auto& constituent : *constituents_)
  {
    constituent->Update(defgrd, params, gp, eleGID);
  }

  mixture_rule_->Update(defgrd, params, gp, eleGID);
}

// This method is called between two timesteps during prestress
void MAT::Mixture::UpdatePrestress(CORE::LINALG::Matrix<3, 3> const& defgrd, const int gp,
    Teuchos::ParameterList& params, const int eleGID)
{
  // Update all constituents
  for (const auto& constituent : *constituents_)
  {
    constituent->UpdatePrestress(defgrd, params, gp, eleGID);
  }

  mixture_rule_->UpdatePrestress(defgrd, params, gp, eleGID);
}

// Evaluates the material
void MAT::Mixture::Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
    const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  // check, whether the PostSetup method was already called
  if (!setup_) dserror("The material's PostSetup() method has not been called yet.");

  if (!isPreEvaluated_[gp])
  {
    for (const auto& constituent : *constituents_)
    {
      isPreEvaluated_[gp] = true;
      constituent->PreEvaluate(*mixture_rule_, params, gp, eleGID);
    }

    mixture_rule_->PreEvaluate(params, gp, eleGID);
  }

  // Evaluate mixturerule
  mixture_rule_->Evaluate(*defgrd, *glstrain, params, *stress, *cmat, gp, eleGID);
}

void MAT::Mixture::RegisterVtkOutputDataNames(
    std::unordered_map<std::string, int>& names_and_size) const
{
  mixture_rule_->RegisterVtkOutputDataNames(names_and_size);
  for (const auto& constituent : *constituents_)
  {
    constituent->RegisterVtkOutputDataNames(names_and_size);
  }
}

bool MAT::Mixture::EvaluateVtkOutputData(
    const std::string& name, Epetra_SerialDenseMatrix& data) const
{
  bool out = mixture_rule_->EvaluateVtkOutputData(name, data);
  for (const auto& constituent : *constituents_)
  {
    out = out || constituent->EvaluateVtkOutputData(name, data);
  }

  return out;
}