/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the mixture material holding a general mixturerule and mixture constituents

\level 3

*/
/*----------------------------------------------------------------------*/

#include "mixture_elasthyper.H"
#include "../drt_mat/material_service.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

// constructor of the parameters
MAT::PAR::Mixture_ElastHyper::Mixture_ElastHyper(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      mass_fractions_(matdata->Get<std::vector<double>>("MASSFRAC")),
      density_(matdata->GetDouble("DENS")),
      constituents_(0)
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

  // check, if the size of the mass fractions fits to the number of constituents
  if (num_constituents != (int)mass_fractions_->size())
  {
    dserror("Number of constituents %d does not fit to the size of the mass fraction vector %d",
        num_constituents, mass_fractions_->size());
  }

  // Create constituents
  for (int i = 0; i < num_constituents; ++i)
  {
    // Create constituent material
    MIXTURE::PAR::MixtureConstituent* mix_const =
        MIXTURE::PAR::MixtureConstituent::Factory((*constituent_matids)[i], (*mass_fractions_)[i]);

    constituents_.emplace_back(mix_const);
  }

  // Create mixture rule
  mixture_rule_ = MIXTURE::PAR::MixtureRule::Factory(matdata->GetInt("MATIDMIXTURERULE"));

  // check, whether the mass frac sums up to 1
  double sum = 0.0;
  for (double massfrac : *mass_fractions_) sum += massfrac;

  if (std::abs(1.0 - sum) > 1e-8) dserror("Mass fractions don't sum up to 1, which is unphysical.");
}

// Create a material instance from parameters
Teuchos::RCP<MAT::Material> MAT::PAR::Mixture_ElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Mixture_ElastHyper(this));
}

MAT::Mixture_ElastHyperType MAT::Mixture_ElastHyperType::instance_;

// Create a material instance from packed data
DRT::ParObject* MAT::Mixture_ElastHyperType::Create(const std::vector<char>& data)
{
  auto* mix_elhy = new MAT::Mixture_ElastHyper();
  mix_elhy->Unpack(data);

  return mix_elhy;
}

// constructor
MAT::Mixture_ElastHyper::Mixture_ElastHyper()
    : params_(nullptr),
      constituents_(Teuchos::rcp(new std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>(0))),
      setup_(false),
      anisotropy_()
{
}

// constructor
MAT::Mixture_ElastHyper::Mixture_ElastHyper(MAT::PAR::Mixture_ElastHyper* params)
    : params_(params),
      constituents_(Teuchos::rcp(new std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>(0))),
      setup_(false),
      anisotropy_()
{
  // create instances of constituents
  for (auto const& constituent : params_->constituents_)
  {
    Teuchos::RCP<MIXTURE::MixtureConstituent> c = constituent->CreateConstituent();
    constituents_->emplace_back(Teuchos::rcp_static_cast<MIXTURE::MixtureConstituent>(c));
    c->SetInitialReferenceDensity(params_->density_);
    c->RegisterAnisotropyExtensions(anisotropy_);
  }

  // create instance of mixture rule
  mixture_rule_ =
      Teuchos::rcp_static_cast<MIXTURE::MixtureRule>(params->mixture_rule_->CreateRule());
  mixture_rule_->SetConstituents(constituents_);
  mixture_rule_->SetMaterialMassDensity(params_->density_);
  mixture_rule_->RegisterAnisotropyExtensions(anisotropy_);
}

// Pack data
void MAT::Mixture_ElastHyper::Pack(DRT::PackBuffer& data) const
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
  AddtoPack(data, static_cast<const int>(setup_));

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
    for (auto const& constituent : *constituents_)
    {
      constituent->PackConstituent(data);
    }

    // pack mixturerule
    mixture_rule_->PackMixtureRule(data);
  }
}

// Unpack data
void MAT::Mixture_ElastHyper::Unpack(const std::vector<char>& data)
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
        params_ = dynamic_cast<MAT::PAR::Mixture_ElastHyper*>(mat);
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
      for (auto const& constituent : params_->constituents_)
      {
        Teuchos::RCP<MIXTURE::MixtureConstituent> c = constituent->CreateConstituent();
        c->SetInitialReferenceDensity(params_->density_);
        constituents_->emplace_back(c);
      }

      // create instance of mixture rule
      mixture_rule_ =
          Teuchos::rcp_static_cast<MIXTURE::MixtureRule>(params_->mixture_rule_->CreateRule());

      // make sure the referenced materials in material list have quick access parameters
      for (auto const& constituent : *constituents_)
      {
        constituent->UnpackConstituent(position, data);
        constituent->RegisterAnisotropyExtensions(anisotropy_);
      }

      // unpack mixturerule
      mixture_rule_->UnpackMixtureRule(position, data);
      mixture_rule_->SetConstituents(constituents_);
      mixture_rule_->SetMaterialMassDensity(params_->density_);
      mixture_rule_->RegisterAnisotropyExtensions(anisotropy_);

      // position checking is not available in post processing mode
      if (position != data.size())
      {
        dserror("Mismatch in size of data to unpack (%d <-> %d)", data.size(), position);
      }
    }
  }
}

// Read element and create arrays for the quantities at the Gauß points
void MAT::Mixture_ElastHyper::Setup(const int numgp, DRT::INPUT::LineDefinition* linedef)
{
  So3Material::Setup(numgp, linedef);

  // resize preevaluation flag
  isPreEvaluated_.resize(numgp, false);

  // Setup anisotropy
  anisotropy_.SetNumberOfGaussPoints(numgp);
  anisotropy_.ReadAnisotropyFromElement(linedef);

  // Let all constituents read the line definition
  for (auto const& constituent : *constituents_)
  {
    constituent->ReadElement(numgp, linedef);
  }

  mixture_rule_->ReadElement(numgp, linedef);
}

// Post setup routine -> Call Setup of constituents and mixture rule
void MAT::Mixture_ElastHyper::PostSetup(Teuchos::ParameterList& params, const int eleGID)
{
  So3Material::PostSetup(params, eleGID);
  anisotropy_.ReadAnisotropyFromParameterList(params);

  for (auto const& constituent : *constituents_)
  {
    constituent->Setup(params, eleGID);
  }

  mixture_rule_->Setup(params, eleGID);

  setup_ = true;
}

// This method is called between two timesteps
void MAT::Mixture_ElastHyper::Update(LINALG::Matrix<3, 3> const& defgrd, const int gp,
    Teuchos::ParameterList& params, const int eleGID)
{
  // Update all constituents
  for (auto const& constituent : *constituents_)
  {
    constituent->Update(defgrd, params, gp, eleGID);
  }

  mixture_rule_->Update(defgrd, params, gp, eleGID);
}

// Evaluates the material
void MAT::Mixture_ElastHyper::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int gp, const int eleGID)
{
  // check, whether the PostSetup method was already called
  if (!setup_) dserror("The material's PostSetup() method has not been called yet.");

  if (!isPreEvaluated_[gp])
  {
    for (auto const& constituent : *constituents_)
    {
      isPreEvaluated_[gp] = true;
      constituent->PreEvaluate(*mixture_rule_, params, gp, eleGID);
    }

    mixture_rule_->PreEvaluate(params, gp, eleGID);
  }

  // Evaluate mixturerule
  mixture_rule_->Evaluate(*defgrd, *glstrain, params, *stress, *cmat, gp, eleGID);
}

// Returns the names of the quantities written during post-processing
void MAT::Mixture_ElastHyper::VisNames(std::map<std::string, int>& names)
{
  mixture_rule_->VisNames(names);
  for (auto const& constituent : *constituents_)
  {
    constituent->VisNames(names);
  }
}

// Returns the names of the quantities written during post-processing
bool MAT::Mixture_ElastHyper::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleGID)
{
  bool vis = mixture_rule_->VisData(name, data, numgp, eleGID);
  if (vis) return true;
  for (auto const& constituent : *constituents_)
  {
    vis = constituent->VisData(name, data, numgp, eleGID);
    if (vis) return true;
  }
  return false;
}
