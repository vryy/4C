/*----------------------------------------------------------------------*/
/*!

\brief Implementation of the base mixture rule

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#include "mixture_rule.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"

// Constructor of the material parameters
MIXTURE::PAR::MixtureRule::MixtureRule(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata)
{
}

// Create rule from the material parameters
Teuchos::RCP<MIXTURE::MixtureRule> MIXTURE::PAR::MixtureRule::CreateRule()
{
  return Teuchos::rcp(new MIXTURE::MixtureRule(this));
}

// Mixture rule factory generates the mixturerule parameters for a specific material id
MIXTURE::PAR::MixtureRule* MIXTURE::PAR::MixtureRule::Factory(int matid)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
    dserror("List of materials cannot be access in the global problem instance.");

  // yet another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matid);

  switch (curmat->Type())
  {
    case INPAR::MAT::mix_rule_base:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MIXTURE::PAR::MixtureRule(curmat));
      auto* params = dynamic_cast<MIXTURE::PAR::MixtureRule*>(curmat->Parameter());
      return params;
    }
    default:
      dserror("The referenced material with id %d is not registered as a Mixture Law!", matid);
  }
  return 0;
}

// Create a mixtureRule from packed data
DRT::ParObject* MIXTURE::MixtureRuleType::Create(const std::vector<char>& data)
{
  auto* mix_elhy = new MIXTURE::MixtureRule();
  mix_elhy->Unpack(data);

  return mix_elhy;
}

MIXTURE::MixtureRuleType MIXTURE::MixtureRuleType::instance_;

// Empty constructor
MIXTURE::MixtureRule::MixtureRule()
    : constituents_(Teuchos::null), numgp_(0), is_init_(false), is_setup_(0)
{
}

// Constructor with parameters
MIXTURE::MixtureRule::MixtureRule(MIXTURE::PAR::MixtureRule* params)
    : constituents_(Teuchos::null), numgp_(0), is_init_(false), is_setup_(0)
{
}

// Pack -> Do not pack this class directly, it should only be packed with MAT::Mixture_ElastHyper
void MIXTURE::MixtureRule::Pack(DRT::PackBuffer& data) const
{
  dserror(
      "This class should not be packed independently. It should only be used in the context "
      "of a Mixture_ElastHyper material");
}

// Unpack -> Do not unpack this class directly, it should only be unpacked with
// MAT::Mixture_ElastHyper
void MIXTURE::MixtureRule::Unpack(const std::vector<char>& data)
{
  dserror(
      "This class should not be unpacked independently. It should only be used in the context "
      "of a Mixture_ElastHyper material");
}

// Pack the mixture law
void MIXTURE::MixtureRule::PackMixtureLaw(DRT::PackBuffer& data) const
{
  // Add flag whether it is initialized
  AddtoPack(data, is_init_);

  // Add flags whether it is setup
  AddtoPack(data, is_setup_);
}

// Unpack the mixture rule
void MIXTURE::MixtureRule::UnpackMixtureLaw(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // Read initialized flag
  ExtractfromPack(position, data, is_init_);

  // Read is setup flag
  ExtractfromPack(position, data, is_setup_);
}

// reads the element definition and set up all quantities
void MIXTURE::MixtureRule::ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Init must only be called once
  if (!is_setup_.empty()) dserror("ReadElement() is called multiple times. Just once allowed.");
  numgp_ = numgp;

  is_setup_.resize(numgp, 0);
}

// Initialize the material rule with the constituents - Called once per element
void MIXTURE::MixtureRule::Init(Teuchos::ParameterList& params,
    const Teuchos::RCP<std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>>& constituents)
{
  if (is_setup_.empty()) dserror("ReadElement() must be called before Init()");
  if (is_init_) dserror("Init() is called more than once. Just once allowed.");
  constituents_ = constituents;

  is_init_ = true;
}

// Setup the mixture rule - Called once per gp
void MIXTURE::MixtureRule::Setup(int gp, Teuchos::ParameterList& params)
{
  // Setup must be called after Init()
  if (!is_init_) dserror("Init() must be called before Setup()!");

  // Setup must only be called once
  if (is_setup_[gp]) dserror("Setup() is called multiple times. Just once per GP allowed.");
  is_setup_[gp] = 1;
}

// Evaluates the stresses of the mixture
void MIXTURE::MixtureRule::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int gp, const int eleGID)
{
  // define temporary matrices
  static LINALG::Matrix<6, 1> cstress;
  static LINALG::Matrix<6, 6> ccmat;

  // This is the simplest mixture law
  // Just iterate over all constituents and add all stress/cmat contributions
  for (auto const& constituent : *Constituents())
  {
    cstress.Clear();
    ccmat.Clear();
    constituent->Evaluate(defgrd, glstrain, params, &cstress, &ccmat, gp, eleGID);

    // Add stress contribution to global stress
    // In this basic mixture law, the mass fractions do not change
    stress->Update(constituent->RefMassFraction(), cstress, 1.0);
    cmat->Update(constituent->RefMassFraction(), ccmat, 1.0);
  }
}
