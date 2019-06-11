/*----------------------------------------------------------------------*/
/*!
\brief Implementation of the mixture material holding a general mixture law and mixture constituents

\level 3

\maintainer Amadeus Gebauer
*/
/*----------------------------------------------------------------------*/

#include "mixture_elasthyper.H"
#include "../drt_mat/material_service.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/// constructor of the parameters
MAT::PAR::Mixture_ElastHyper::Mixture_ElastHyper(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      mass_fractions_(matdata->Get<std::vector<double>>("MASSFRAC")),
      density_(matdata->GetDouble("DENS")),
      constituents_(0)
{
  const int num_constituents = matdata->GetInt("NUMCONST");
  const auto* constituent_matids = matdata->Get<std::vector<int>>("MATIDSCONST");

  // check, if size of constituents fits to the number of constituents
  if (num_constituents != (int)constituent_matids->size())
    dserror(
        "number of constituents %d does not fit to the size of the constituents material vector"
        " %d",
        num_constituents, constituent_matids->size());

  // check, if the size of the mass fractions fits to the number of constituents
  if (num_constituents != (int)mass_fractions_->size())
    dserror("Number of constituents %d does not fit to the size of the mass fraction vector %d",
        num_constituents, mass_fractions_->size());

  // Create constituents
  for (int i = 0; i < num_constituents; ++i)
  {
    // Create constituent material
    MIXTURE::PAR::MixtureConstituent* mix_const =
        MIXTURE::PAR::MixtureConstituent::Factory((*constituent_matids)[i], (*mass_fractions_)[i]);

    constituents_.emplace_back(mix_const);
  }

  mixture_rule_ = MIXTURE::PAR::MixtureRule::Factory(matdata->GetInt("MATIDMIXTURELAW"));


  // check, whether the mass frac sums up to 1
  double sum = 0.0;
  for (double massfrac : *mass_fractions_) sum += massfrac;

  if (std::abs(1.0 - sum) > 1e-8) dserror("Mass fractions don't sum up to 1, which is unphysical.");
}

/// Create a material instance from parameters
Teuchos::RCP<MAT::Material> MAT::PAR::Mixture_ElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Mixture_ElastHyper(this));
}

/// Create a material instance from packed data
DRT::ParObject* MAT::Mixture_ElastHyperType::Create(const std::vector<char>& data)
{
  auto* mix_elhy = new MAT::Mixture_ElastHyper();
  mix_elhy->Unpack(data);

  return mix_elhy;
}

MAT::Mixture_ElastHyperType MAT::Mixture_ElastHyperType::instance_;

/// constructor
MAT::Mixture_ElastHyper::Mixture_ElastHyper()
    : params_(nullptr),
      constituents_(Teuchos::rcp(new std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>(0))),
      setup_(0),
      constituent_init_(false)
{
}

/// constructor
MAT::Mixture_ElastHyper::Mixture_ElastHyper(MAT::PAR::Mixture_ElastHyper* params)
    : params_(params),
      constituents_(Teuchos::rcp(new std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>(0))),
      setup_(0),
      constituent_init_(false)
{
  // create instances of constituents
  for (auto const& constituent : params_->constituents_)
  {
    constituents_->emplace_back(
        Teuchos::rcp_static_cast<MIXTURE::MixtureConstituent>(constituent->CreateConstituent()));
  }

  // create instance of mixture rule
  mixture_rule_ =
      Teuchos::rcp_static_cast<MIXTURE::MixtureRule>(params->mixture_rule_->CreateRule());
}

/// Pack data
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
  AddtoPack(data, setup_);

  // pack init flag
  AddtoPack(data, constituent_init_);

  // pack all constituents
  for (auto const& constituent : *constituents_)
  {
    constituent->PackConstituent(data);
  }

  // pack mixture law
  mixture_rule_->PackMixtureLaw(data);
}

// Unpack data
void MAT::Mixture_ElastHyper::Unpack(const std::vector<char>& data)
{
  params_ = nullptr;
  constituents_->clear();
  setup_.clear();

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
      const unsigned int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<MAT::PAR::Mixture_ElastHyper*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

    // Extract setup flag
    ExtractfromPack(position, data, setup_);

    // Extract init flag
    ExtractfromPack(position, data, constituent_init_);

    // extract constituents
    // constituents are not accessible during post processing
    if (params_ != nullptr)
    {
      // create instances of constituents
      for (auto const& constituent : params_->constituents_)
      {
        constituents_->emplace_back(Teuchos::rcp_static_cast<MIXTURE::MixtureConstituent>(
            constituent->CreateConstituent()));
      }

      // create instance of mixture rule
      mixture_rule_ =
          Teuchos::rcp_static_cast<MIXTURE::MixtureRule>(params_->mixture_rule_->CreateRule());

      // make sure the referenced materials in material list have quick access parameters
      for (auto const& constituent : *constituents_)
      {
        constituent->UnpackConstituent(position, data);
      }

      // unpack mixture law
      mixture_rule_->UnpackMixtureLaw(position, data);
    }
  }
}

///! Read element and create arrays for the quantities at the Gauß points
void MAT::Mixture_ElastHyper::Setup(const int numgp, DRT::INPUT::LineDefinition* linedef)
{
  So3Material::Setup(numgp, linedef);

  setup_.resize(numgp, 0);

  // Let all constituents read the line definition
  for (auto const& constituent : *constituents_)
  {
    constituent->ReadElement(numgp, linedef);
  }

  mixture_rule_->ReadElement(numgp, linedef);
}

/// This method is called between two timesteps
void MAT::Mixture_ElastHyper::Update(LINALG::Matrix<3, 3> const& defgrd, const int gp,
    Teuchos::ParameterList& params, const int eleGID)
{
  // Update all constituents
  for (auto const& constituent : *constituents_)
  {
    constituent->Update(defgrd, params, gp);
  }

  mixture_rule_->Update(defgrd, params, gp);
}

/// Evaluates the material
void MAT::Mixture_ElastHyper::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  // read Gauß point
  int gp = params.get<int>("gp");

  if (!constituent_init_)
  {
    // constituents are not initialized yet, do it
    for (auto const& constituent : *constituents_)
    {
      constituent->Init(params);
    }
    mixture_rule_->Init(params, constituents_);
    constituent_init_ = true;
  }

  // check if material is already set up
  if (!setup_[gp])
  {
    // constituents are not set up yet, do it
    for (auto const& constituent : *constituents_)
    {
      constituent->Setup(gp, params);
    }
    mixture_rule_->Setup(gp, params);

    // mark constituents as setup
    setup_[gp] = 1;
  }

  // Evaluate mixture law
  mixture_rule_->Evaluate(defgrd, glstrain, params, stress, cmat, gp, eleGID);
}

/// Returns the names of the quantities written during post-processing
void MAT::Mixture_ElastHyper::VisNames(std::map<std::string, int>& names)
{
  mixture_rule_->VisNames(names);
  for (auto const& constituent : *constituents_)
  {
    constituent->VisNames(names);
  }
}

/// Returns the names of the quantities written during post-processing
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
