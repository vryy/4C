/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the hyperelastic constituent

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#include "mixture_constituent_elasthyper.H"
#include "../drt_mat/material_service.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/elasthyper.H"

// Constructor for the parameter class
MIXTURE::PAR::MixtureConstituent_ElastHyper::MixtureConstituent_ElastHyper(
    const Teuchos::RCP<MAT::PAR::Material>& matdata, const double ref_mass_fraction)
    : MixtureConstituent(matdata, ref_mass_fraction),
      nummat_(matdata->GetInt("NUMMAT")),
      matids_(matdata->Get<std::vector<int>>("MATIDS"))
{
  // check, if size of summands fits to the number of summands
  if (nummat_ != (int)matids_->size())
    dserror(
        "number of summands %d does not fit to the size of the summands vector"
        " %d",
        nummat_, matids_->size());
}

// Create an instance of MIXTURE::MixtureConstituent_ElastHyper from the parameters
Teuchos::RCP<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituent_ElastHyper::CreateConstituent()
{
  return Teuchos::rcp(new MIXTURE::MixtureConstituent_ElastHyper(this));
}

// Empty constructor of the constituent
MIXTURE::MixtureConstituent_ElastHyper::MixtureConstituent_ElastHyper()
    : MixtureConstituent(), summandProperties_(), params_(nullptr), potsum_(0)
{
  // empty constructor
}

// Constructor of the constituent holding the material parameters
MIXTURE::MixtureConstituent_ElastHyper::MixtureConstituent_ElastHyper(
    MIXTURE::PAR::MixtureConstituent_ElastHyper* params)
    : MixtureConstituent(), summandProperties_(), params_(params), potsum_(0)
{
  std::vector<int>::const_iterator m;

  // Create constituents
  for (m = params_->matids_->begin(); m != params_->matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to read elastic summand.");
    potsum_.push_back(sum);
  }
}

// Pack the constituent
void MIXTURE::MixtureConstituent_ElastHyper::PackConstituent(DRT::PackBuffer& data) const
{
  MixtureConstituent::PackConstituent(data);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  DRT::ParObject::AddtoPack(data, matid);
  summandProperties_.Pack(data);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (const auto& p : potsum_) p->PackSummand(data);
  }
}

// Unpack the constituent
void MIXTURE::MixtureConstituent_ElastHyper::UnpackConstituent(
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
        params_ = dynamic_cast<MIXTURE::PAR::MixtureConstituent_ElastHyper*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  summandProperties_.Unpack(position, data);

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

    // loop map of associated potential summands
    for (auto& summand : potsum_) summand->UnpackSummand(data, position);
  }
}

// Returns the material type
INPAR::MAT::MaterialType MIXTURE::MixtureConstituent_ElastHyper::MaterialType() const
{
  return INPAR::MAT::mix_elasthyper;
}

// Reads the element from the input file
void MIXTURE::MixtureConstituent_ElastHyper::ReadElement(
    int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Call base class
  MixtureConstituent::ReadElement(numgp, linedef);

  // Setup summands
  for (const auto& summand : potsum_) summand->Setup(numgp, linedef);

  // find out which formulations are used
  MAT::ElastHyperProperties(potsum_, summandProperties_);

  if (summandProperties_.viscoGeneral)
    dserror("Never use viscoelastic materials in Elasthyper-Toolbox.");
}

// Evaluates the stress of the constituent
void MIXTURE::MixtureConstituent_ElastHyper::Evaluate(const LINALG::Matrix<3, 3>& F,
    const LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, const int gp, const int eleGID)
{
  // 2nd Piola-Kirchhoff stress tensor in stress-like Voigt notation of the constituent
  static LINALG::Matrix<6, 1> Sc_stress(false);
  Sc_stress.Clear();
  // Constitutive tensor of constituent
  static LINALG::Matrix<6, 6> ccmat(false);
  ccmat.Clear();
  // Evaluate stresses using ElastHyper service functions
  MAT::ElastHyperEvaluate(
      F, E_strain, params, Sc_stress, ccmat, eleGID, potsum_, summandProperties_, false);

  S_stress.Update(CurrentRefDensity(), Sc_stress, 1.0);
  cmat.Update(CurrentRefDensity(), ccmat, 1.0);
}

// Returns the reference mass fraction of the constituent
double MIXTURE::MixtureConstituent_ElastHyper::CurrentRefDensity() const
{
  return params_->RefMassFraction() * InitialRefDensity();
}

// Updates all summands
void MIXTURE::MixtureConstituent_ElastHyper::Update(LINALG::Matrix<3, 3> const& defgrd,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  MixtureConstituent::Update(defgrd, params, gp, eleGID);

  // loop map of associated potential summands
  for (auto& summand : potsum_) summand->Update();
}

// Add names for each summand for the quantities for post processing
void MIXTURE::MixtureConstituent_ElastHyper::VisNames(std::map<std::string, int>& names)
{
  MixtureConstituent::VisNames(names);

  // loop map of associated potential summands
  for (auto& summand : potsum_) summand->VisNames(names);
}

// Add values for each summand of the quantities for post processing
bool MIXTURE::MixtureConstituent_ElastHyper::VisData(
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
